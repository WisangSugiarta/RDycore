#ifndef SWE_OPERATORS_IMPL_H
#define SWE_OPERATORS_IMPL_H

#include "swe_ceed.h"

#ifndef Square
#define Square(x) ((x) * (x))
#endif
#ifndef SafeDiv
#define SafeDiv(a, b, c, tiny) ((c) > (tiny) ? (a) / (b) : 0.0)
#endif

// we disable compiler warnings for implicitly-declared math functions known to
// the JIT compiler
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-function-declaration"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wimplicit-function-declaration"

// Q-function context with data attached
typedef struct SWEContext_ *SWEContext;
struct SWEContext_ {
  CeedScalar dtime;
  CeedScalar tiny_h;
  CeedScalar h_anuga_regular;
  CeedScalar gravity;
  CeedScalar xq2018_threshold;
};

struct SWEState_ {
  CeedScalar h, hu, hv;
};
typedef struct SWEState_ SWEState;

// supported Riemann solver types
#include "swe_roe_ceed_impl.h"
typedef enum {
  RIEMANN_FLUX_ROE,
} RiemannFluxType;

// The following Q functions use C99 VLA features for shaping multidimensional
// arrays, which don't have the same drawbacks as VLA allocations.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wvla"

// SWE interior flux operator Q-function
CEED_QFUNCTION_HELPER int SWEFlux(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[], RiemannFluxType flux_type) {
  const CeedScalar(*geom)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // sn, cn, weight_L, weight_R
  const CeedScalar(*q_L)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];
  const CeedScalar(*q_R)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[2];
  CeedScalar(*cell_L)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*cell_R)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[1];
  CeedScalar(*accum_flux)[CEED_Q_VLA]  = (CeedScalar(*)[CEED_Q_VLA])out[2];
  CeedScalar(*courant_num)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[3];
  const SWEContext context             = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    SWEState   qL = {q_L[0][i], q_L[1][i], q_L[2][i]};
    SWEState   qR = {q_R[0][i], q_R[1][i], q_R[2][i]};
    CeedScalar flux[3], amax;
    if (qL.h > tiny_h || qR.h > tiny_h) {
      switch (flux_type) {
        case RIEMANN_FLUX_ROE:
          SWERiemannFlux_Roe(gravity, tiny_h, h_anuga, qL, qR, geom[0][i], geom[1][i], flux, &amax);
          break;
      }
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = flux[j] * geom[2][i];
        cell_R[j][i]     = flux[j] * geom[3][i];
        accum_flux[j][i] = flux[j];
      }
      courant_num[0][i] = -amax * geom[2][i] * dt;
      courant_num[1][i] = amax * geom[3][i] * dt;
    } else {
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = 0.0;
        cell_R[j][i]     = 0.0;
        accum_flux[j][i] = 0.0;
      }
      courant_num[0][i] = 0.0;
      courant_num[1][i] = 0.0;
    }
  }
  return 0;
}

CEED_QFUNCTION(SWEFlux_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEFlux(ctx, Q, in, out, RIEMANN_FLUX_ROE);
}

CEED_QFUNCTION_HELPER CeedScalar MinModLimiter_CEED(CeedScalar a, CeedScalar b) {
  if (a * b <= 0.0) return 0.0;  // Different signs or zero
  return (fabs(a) < fabs(b)) ? a : b;
}
CEED_QFUNCTION_HELPER int SWEFluxReconstructed(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[], RiemannFluxType flux_type) {
  // Inputs (must match the operator field order)
  const CeedScalar(*geom)[CEED_Q_VLA]           = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // sn, cn, weight_L, weight_R
  const CeedScalar(*grad_geom)[CEED_Q_VLA]      = (const CeedScalar(*)[CEED_Q_VLA])in[1];  // dx_L, dy_L, dx_R, dy_R, dist_sq, edge_len
  const CeedScalar(*q_L)[CEED_Q_VLA]            = (const CeedScalar(*)[CEED_Q_VLA])in[2];  // left cell values
  const CeedScalar(*q_R)[CEED_Q_VLA]            = (const CeedScalar(*)[CEED_Q_VLA])in[3];  // right cell values
  const CeedScalar(*q_L_neighbors)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[4];  // left cell neighbors (2 neighbors * 3 components)
  const CeedScalar(*q_R_neighbors)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[5];  // right cell neighbors (2 neighbors * 3 components)

  // Outputs
  CeedScalar(*cell_L)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*cell_R)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[1];
  CeedScalar(*accum_flux)[CEED_Q_VLA]  = (CeedScalar(*)[CEED_Q_VLA])out[2];
  CeedScalar(*courant_num)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[3];
  
  const SWEContext context = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    // Extract geometric factors for gradient computation
    const CeedScalar dxL[2]  = {grad_geom[0][i], grad_geom[1][i]};  // displacement from left centroid to midpoint
    const CeedScalar dxR[2]  = {grad_geom[2][i], grad_geom[3][i]};  // displacement from right centroid to midpoint
    const CeedScalar dist_sq = grad_geom[4][i];                     // distance squared between cell centers
    const CeedScalar edge_len = grad_geom[5][i];                    // edge length

    // Compute geometric weights for gradients
    const CeedScalar dx_cell = sqrt(dist_sq);
    const CeedScalar grad_weight_x = (dxR[0] - dxL[0]) / dist_sq;  // normalized direction in x
    const CeedScalar grad_weight_y = (dxR[1] - dxL[1]) / dist_sq;  // normalized direction in y

    // Get current cell values
    SWEState qL_center = {q_L[0][i], q_L[1][i], q_L[2][i]};
    SWEState qR_center = {q_R[0][i], q_R[1][i], q_R[2][i]};

    // Get neighbor values for slope limiting
    SWEState qL_neighbor = {q_L_neighbors[0][i], q_L_neighbors[1][i], q_L_neighbors[2][i]};
    SWEState qR_neighbor = {q_R_neighbors[3][i], q_R_neighbors[4][i], q_R_neighbors[5][i]};

    // Reconstruct each component with slope limiting
    SWEState qL_reconstructed, qR_reconstructed;
    
    // Component 0: h (water height)
    {
      const CeedScalar dq = qR_center.h - qL_center.h;
      
      // Unlimited gradient
      const CeedScalar grad_unlimited[2] = {dq * grad_weight_x, dq * grad_weight_y};
      
      // Slope limiting using neighbor information
      CeedScalar grad_limited[2];
      
      // Compute neighbor gradients for limiting
      const CeedScalar dq_L_neighbor = qL_center.h - qL_neighbor.h;
      const CeedScalar dq_R_neighbor = qR_neighbor.h - qR_center.h;
      
      // Apply MinMod limiter in each direction
      const CeedScalar grad_L_x = dq_L_neighbor * grad_weight_x;
      const CeedScalar grad_R_x = dq_R_neighbor * grad_weight_x;
      grad_limited[0] = MinModLimiter_CEED(grad_unlimited[0], 
                                          MinModLimiter_CEED(2.0 * grad_L_x, 2.0 * grad_R_x));
      
      const CeedScalar grad_L_y = dq_L_neighbor * grad_weight_y;
      const CeedScalar grad_R_y = dq_R_neighbor * grad_weight_y;
      grad_limited[1] = MinModLimiter_CEED(grad_unlimited[1], 
                                          MinModLimiter_CEED(2.0 * grad_L_y, 2.0 * grad_R_y));
      
      // Reconstruct values at edge midpoint
      qL_reconstructed.h = fmax(0.0, qL_center.h + grad_limited[0] * dxL[0] + grad_limited[1] * dxL[1]);
      qR_reconstructed.h = fmax(0.0, qR_center.h + grad_limited[0] * dxR[0] + grad_limited[1] * dxR[1]);
    }
    
    // Component 1: hu (x-momentum)  
    {
      const CeedScalar dq = qR_center.hu - qL_center.hu;
      const CeedScalar grad_unlimited[2] = {dq * grad_weight_x, dq * grad_weight_y};
      
      const CeedScalar dq_L_neighbor = qL_center.hu - qL_neighbor.hu;
      const CeedScalar dq_R_neighbor = qR_neighbor.hu - qR_center.hu;
      
      const CeedScalar grad_L_x = dq_L_neighbor * grad_weight_x;
      const CeedScalar grad_R_x = dq_R_neighbor * grad_weight_x;
      const CeedScalar grad_limited_x = MinModLimiter_CEED(grad_unlimited[0], 
                                                          MinModLimiter_CEED(2.0 * grad_L_x, 2.0 * grad_R_x));
      
      const CeedScalar grad_L_y = dq_L_neighbor * grad_weight_y;
      const CeedScalar grad_R_y = dq_R_neighbor * grad_weight_y;
      const CeedScalar grad_limited_y = MinModLimiter_CEED(grad_unlimited[1], 
                                                          MinModLimiter_CEED(2.0 * grad_L_y, 2.0 * grad_R_y));
      
      qL_reconstructed.hu = qL_center.hu + grad_limited_x * dxL[0] + grad_limited_y * dxL[1];
      qR_reconstructed.hu = qR_center.hu + grad_limited_x * dxR[0] + grad_limited_y * dxR[1];
    }
    
    // Component 2: hv (y-momentum)
    {
      const CeedScalar dq = qR_center.hv - qL_center.hv;
      const CeedScalar grad_unlimited[2] = {dq * grad_weight_x, dq * grad_weight_y};
      
      const CeedScalar dq_L_neighbor = qL_center.hv - qL_neighbor.hv;
      const CeedScalar dq_R_neighbor = qR_neighbor.hv - qR_center.hv;
      
      const CeedScalar grad_L_x = dq_L_neighbor * grad_weight_x;
      const CeedScalar grad_R_x = dq_R_neighbor * grad_weight_x;
      const CeedScalar grad_limited_x = MinModLimiter_CEED(grad_unlimited[0], 
                                                          MinModLimiter_CEED(2.0 * grad_L_x, 2.0 * grad_R_x));
      
      const CeedScalar grad_L_y = dq_L_neighbor * grad_weight_y;
      const CeedScalar grad_R_y = dq_R_neighbor * grad_weight_y;
      const CeedScalar grad_limited_y = MinModLimiter_CEED(grad_unlimited[1], 
                                                          MinModLimiter_CEED(2.0 * grad_L_y, 2.0 * grad_R_y));
      
      qL_reconstructed.hv = qL_center.hv + grad_limited_x * dxL[0] + grad_limited_y * dxL[1];
      qR_reconstructed.hv = qR_center.hv + grad_limited_x * dxR[0] + grad_limited_y * dxR[1];
    }

    // Compute Riemann flux using reconstructed states
    if (qL_reconstructed.h > tiny_h || qR_reconstructed.h > tiny_h) {
      CeedScalar flux[3], amax;
      switch (flux_type) {
        case RIEMANN_FLUX_ROE:
          SWERiemannFlux_Roe(gravity, tiny_h, h_anuga, qL_reconstructed, qR_reconstructed, geom[0][i], geom[1][i], flux, &amax);
          break;
      }
      
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = flux[j] * geom[2][i];
        cell_R[j][i]     = flux[j] * geom[3][i];
        accum_flux[j][i] = flux[j];
      }
      courant_num[0][i] = -amax * geom[2][i] * dt;
      courant_num[1][i] = amax * geom[3][i] * dt;
    } else {
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = 0.0;
        cell_R[j][i]     = 0.0;
        accum_flux[j][i] = 0.0;
      }
      courant_num[0][i] = 0.0;
      courant_num[1][i] = 0.0;
    }
  }
  return 0;
}

CEED_QFUNCTION(SWEFluxReconstructed_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEFluxSlopeLimited(ctx, Q, in, out, RIEMANN_FLUX_ROE);
}


// SWE boundary flux operator Q-function (Dirichlet condition)
CEED_QFUNCTION_HELPER int SWEBoundaryFlux_Dirichlet(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[],
                                                    RiemannFluxType flux_type) {
  const CeedScalar(*geom)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // sn, cn, weight_L
  const CeedScalar(*q_L)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];
  const CeedScalar(*q_R)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[2];  // Dirichlet boundary values
  CeedScalar(*cell_L)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*accum_flux)[CEED_Q_VLA]  = (CeedScalar(*)[CEED_Q_VLA])out[1];
  CeedScalar(*courant_num)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[2];
  const SWEContext context             = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    SWEState qL = {q_L[0][i], q_L[1][i], q_L[2][i]};
    SWEState qR = {q_R[0][i], q_R[1][i], q_R[2][i]};
    if (qL.h > tiny_h) {
      CeedScalar flux[3], amax;
      switch (flux_type) {
        case RIEMANN_FLUX_ROE:
          SWERiemannFlux_Roe(gravity, tiny_h, h_anuga, qL, qR, geom[0][i], geom[1][i], flux, &amax);
          break;
      }
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = flux[j] * geom[2][i];
        accum_flux[j][i] = flux[j];
      }
      courant_num[0][i] = -amax * geom[2][i] * dt;
    } else {
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = 0.0;
        accum_flux[j][i] = 0.0;
      }
      courant_num[0][i] = 0.0;
    }
  }
  return 0;
}

CEED_QFUNCTION(SWEBoundaryFlux_Dirichlet_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEBoundaryFlux_Dirichlet(ctx, Q, in, out, RIEMANN_FLUX_ROE);
}

// SWE boundary flux operator Q-function (reflecting condition)
CEED_QFUNCTION_HELPER int SWEBoundaryFlux_Reflecting(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[],
                                                     RiemannFluxType flux_type) {
  const CeedScalar(*geom)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // sn, cn, weight_L
  const CeedScalar(*q_L)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];
  CeedScalar(*cell_L)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*courant_num)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[2];
  const SWEContext context             = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    CeedScalar sn = geom[0][i], cn = geom[1][i];
    SWEState   qL = {q_L[0][i], q_L[1][i], q_L[2][i]};
    if (qL.h > tiny_h) {
      CeedScalar dum1 = sn * sn - cn * cn;
      CeedScalar dum2 = 2.0 * sn * cn;
      SWEState   qR   = {qL.h, qL.hu * dum1 - qL.hv * dum2, -qL.hu * dum2 - qL.hv * dum1};
      CeedScalar flux[3], amax;
      switch (flux_type) {
        case RIEMANN_FLUX_ROE:
          SWERiemannFlux_Roe(gravity, tiny_h, h_anuga, qL, qR, sn, cn, flux, &amax);
          break;
      }
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i] = flux[j] * geom[2][i];
      }
      courant_num[0][i] = -amax * geom[2][i] * dt;
    } else {
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i] = 0.0;
      }
      courant_num[0][i] = 0.0;
    }
  }
  return 0;
}

CEED_QFUNCTION(SWEBoundaryFlux_Reflecting_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEBoundaryFlux_Reflecting(ctx, Q, in, out, RIEMANN_FLUX_ROE);
}

// SWE boundary flux operator Q-function (outflow condition)
CEED_QFUNCTION_HELPER int SWEBoundaryFlux_Outflow(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[],
                                                  RiemannFluxType flux_type) {
  const CeedScalar(*geom)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // sn, cn, weight_L
  const CeedScalar(*q_L)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];
  CeedScalar(*cell_L)[CEED_Q_VLA]      = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*accum_flux)[CEED_Q_VLA]  = (CeedScalar(*)[CEED_Q_VLA])out[1];
  CeedScalar(*courant_num)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[2];
  const SWEContext context             = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    CeedScalar sn = geom[0][i], cn = geom[1][i];
    SWEState   qL    = {q_L[0][i], q_L[1][i], q_L[2][i]};
    CeedScalar q     = fabs(qL.hu * cn + qL.hv * sn);
    CeedScalar hR    = pow(q * q / gravity, 1.0 / 3.0);
    CeedScalar speed = sqrt(gravity * hR);
    SWEState   qR    = {hR, hR * speed * cn, hR * speed * sn};
    if (qL.h > tiny_h || qR.h > tiny_h) {
      CeedScalar flux[3], amax;
      switch (flux_type) {
        case RIEMANN_FLUX_ROE:
          SWERiemannFlux_Roe(gravity, tiny_h, h_anuga, qL, qR, sn, cn, flux, &amax);
          break;
      }
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = flux[j] * geom[2][i];
        accum_flux[j][i] = flux[j];
      }
      courant_num[0][i] = -amax * geom[2][i] * dt;
    } else {
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i]     = 0.0;
        accum_flux[j][i] = 0.0;
      }
      courant_num[0][i] = 0.0;
    }
  }
  return 0;
}

CEED_QFUNCTION(SWEBoundaryFlux_Outflow_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEBoundaryFlux_Outflow(ctx, Q, in, out, RIEMANN_FLUX_ROE);
}

// SWE regional source operator Q-function
CEED_QFUNCTION(SWESourceTermSemiImplicit)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  const CeedScalar(*geom)[CEED_Q_VLA]      = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // dz/dx, dz/dy
  const CeedScalar(*ext_src)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];  // external source (e.g. rain rate)
  const CeedScalar(*mat_props)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2];  // material properties
  const CeedScalar(*riemannf)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[3];  // riemann flux
  const CeedScalar(*q)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[4];
  CeedScalar(*cell)[CEED_Q_VLA]            = (CeedScalar(*)[CEED_Q_VLA])out[0];
  const SWEContext context                 = (SWEContext)ctx;

  const CeedScalar dt      = context->dtime;
  const CeedScalar tiny_h  = context->tiny_h;
  const CeedScalar h_anuga = context->h_anuga_regular;
  const CeedScalar gravity = context->gravity;

  for (CeedInt i = 0; i < Q; i++) {
    SWEState         state = {q[0][i], q[1][i], q[2][i]};
    const CeedScalar h     = state.h;
    const CeedScalar hu    = state.hu;
    const CeedScalar hv    = state.hv;
    const CeedScalar denom = Square(h) + Square(h_anuga);

    const CeedScalar u = SafeDiv(state.hu * h, denom, h, tiny_h);
    const CeedScalar v = SafeDiv(state.hv * h, denom, h, tiny_h);

    const CeedScalar dz_dx = geom[0][i];
    const CeedScalar dz_dy = geom[1][i];

    const CeedScalar bedx = dz_dx * gravity * h;
    const CeedScalar bedy = dz_dy * gravity * h;

    const CeedScalar Fsum_x = riemannf[1][i];
    const CeedScalar Fsum_y = riemannf[2][i];

    CeedScalar tbx = 0.0, tby = 0.0;
    if (h > tiny_h) {
      const CeedScalar mannings_n = mat_props[MATERIAL_PROPERTY_MANNINGS][i];
      const CeedScalar Cd         = gravity * Square(mannings_n) * pow(h, -1.0 / 3.0);

      const CeedScalar velocity = sqrt(Square(u) + Square(v));

      const CeedScalar tb = Cd * velocity / h;

      const CeedScalar factor = tb / (1.0 + dt * tb);

      tbx = (hu + dt * Fsum_x - dt * bedx) * factor;
      tby = (hv + dt * Fsum_y - dt * bedy) * factor;
    }

    cell[0][i] = riemannf[0][i] + ext_src[0][i];
    cell[1][i] = riemannf[1][i] - bedx - tbx + ext_src[1][i];
    cell[2][i] = riemannf[2][i] - bedy - tby + ext_src[2][i];
  }
  return 0;
}

/// @brief Adds contribution of the source-term using implicit time integration approach of:
///        Xia, Xilin, and Qiuhua Liang. "A new efficient implicit scheme for discretising the stiff
///        friction terms in the shallow water equations." Advances in water resources 117 (2018): 87-97.
///        https://www.sciencedirect.com/science/article/pii/S0309170818302124?ref=cra_js_challenge&fr=RR-1
CEED_QFUNCTION(SWESourceTermImplicitXQ2018)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  const CeedScalar(*geom)[CEED_Q_VLA]      = (const CeedScalar(*)[CEED_Q_VLA])in[0];  // dz/dx, dz/dy
  const CeedScalar(*ext_src)[CEED_Q_VLA]   = (const CeedScalar(*)[CEED_Q_VLA])in[1];  // external source (e.g. rain rate)
  const CeedScalar(*mat_props)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2];  // material properties
  const CeedScalar(*riemannf)[CEED_Q_VLA]  = (const CeedScalar(*)[CEED_Q_VLA])in[3];  // riemann flux
  const CeedScalar(*q)[CEED_Q_VLA]         = (const CeedScalar(*)[CEED_Q_VLA])in[4];
  CeedScalar(*cell)[CEED_Q_VLA]            = (CeedScalar(*)[CEED_Q_VLA])out[0];
  const SWEContext context                 = (SWEContext)ctx;

  const CeedScalar dt               = context->dtime;
  const CeedScalar tiny_h           = context->tiny_h;
  const CeedScalar gravity          = context->gravity;
  const CeedScalar xq2018_threshold = context->xq2018_threshold;

  for (CeedInt i = 0; i < Q; i++) {
    SWEState         state = {q[0][i], q[1][i], q[2][i]};
    const CeedScalar h     = state.h;
    const CeedScalar hu    = state.hu;
    const CeedScalar hv    = state.hv;

    const CeedScalar dz_dx = geom[0][i];
    const CeedScalar dz_dy = geom[1][i];

    const CeedScalar bedx = dz_dx * gravity * h;
    const CeedScalar bedy = dz_dy * gravity * h;

    const CeedScalar Fsum_x = riemannf[1][i];
    const CeedScalar Fsum_y = riemannf[2][i];

    CeedScalar tbx = 0.0, tby = 0.0;
    if (h > tiny_h) {
      // defined in the text below equation 22 of XQ2018
      const CeedScalar Ax = Fsum_x - bedx;
      const CeedScalar Ay = Fsum_y - bedy;

      // equation 27 of XQ2018
      const CeedScalar mx = hu + Ax * dt;
      const CeedScalar my = hv + Ay * dt;

      const CeedScalar mannings_n = mat_props[MATERIAL_PROPERTY_MANNINGS][i];
      const CeedScalar lambda     = gravity * Square(mannings_n) * pow(h, -4.0 / 3.0) * pow(Square(mx / h) + Square(my / h), 0.5);

      CeedScalar qx_nplus1 = 0.0, qy_nplus1 = 0.0;

      // equation 36 and 37 of XQ2018
      if (dt * lambda < xq2018_threshold) {
        qx_nplus1 = mx;
        qy_nplus1 = my;
      } else {
        qx_nplus1 = (mx - mx * pow(1.0 + 4.0 * dt * lambda, 0.5)) / (-2.0 * dt * lambda);
        qy_nplus1 = (my - my * pow(1.0 + 4.0 * dt * lambda, 0.5)) / (-2.0 * dt * lambda);
      }

      const CeedScalar q_magnitude = pow(Square(qx_nplus1) + Square(qy_nplus1), 0.5);

      // equation 21 and 22 of XQ2018
      tbx = gravity * Square(mannings_n) * pow(h, -7.0 / 3.0) * qx_nplus1 * q_magnitude;
      tby = gravity * Square(mannings_n) * pow(h, -7.0 / 3.0) * qy_nplus1 * q_magnitude;
    }

    cell[0][i] = riemannf[0][i] + ext_src[0][i];
    cell[1][i] = riemannf[1][i] - bedx - tbx + ext_src[1][i];
    cell[2][i] = riemannf[2][i] - bedy - tby + ext_src[2][i];
  }
  return 0;
}

#pragma GCC diagnostic   pop
#pragma clang diagnostic pop

#endif
