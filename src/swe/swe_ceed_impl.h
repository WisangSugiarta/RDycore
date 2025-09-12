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

CEED_QFUNCTION_HELPER int SWEFluxReconstructionKernel(void *ctx, CeedInt Q, 
                                                     const CeedScalar *const in[], 
                                                     CeedScalar *const out[]) {
  // Inputs
  const CeedScalar(*edge_data)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[0];
  const CeedScalar(*cell_data)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[1];
  const CeedScalar(*neighbor_coords)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[2];
  const CeedScalar(*neighbor_values)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[3];
  const CeedScalar(*q_left)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[4];
  const CeedScalar(*q_right)[CEED_Q_VLA] = (const CeedScalar(*)[CEED_Q_VLA])in[5];
  
  // Outputs
  CeedScalar(*cell_L)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[0];
  CeedScalar(*cell_R)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[1];
  CeedScalar(*flux)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[2];
  CeedScalar(*courant)[CEED_Q_VLA] = (CeedScalar(*)[CEED_Q_VLA])out[3];

  const SWEContext context = (SWEContext)ctx;
  const CeedScalar gravity = context->gravity;
  const CeedScalar tiny_h = context->tiny_h;
  const CeedScalar dt = context->dtime;
  const CeedScalar h_anuga = context->h_anuga_regular;

  for (CeedInt i = 0; i < Q; i++) {
    // Extract edge geometry
    CeedScalar sn = edge_data[0][i];
    CeedScalar cn = edge_data[1][i];
    CeedScalar edge_length = edge_data[2][i];
    CeedScalar edge_mid_x = edge_data[3][i];
    CeedScalar edge_mid_y = edge_data[4][i];
    
    // Extract cell geometry
    CeedScalar xl = cell_data[0][i], yl = cell_data[1][i], Al = cell_data[2][i];
    CeedScalar xr = cell_data[3][i], yr = cell_data[4][i], Ar = cell_data[5][i];
    
    // Get left and right cell solutions
    CeedScalar qL[3] = {q_left[0][i], q_left[1][i], q_left[2][i]};
    CeedScalar qR[3] = {q_right[0][i], q_right[1][i], q_right[2][i]};
    
    // Initialize reconstructed values with cell center values
    CeedScalar qL_recon[3] = {qL[0], qL[1], qL[2]};
    CeedScalar qR_recon[3] = {qR[0], qR[1], qR[2]};
    
    // **LEFT CELL RECONSTRUCTION**
    if (qL[0] > tiny_h) {
      CeedScalar A[2][2] = {{0}}, b[3][2] = {{0}};
      CeedInt valid_neighbors = 0;
      
      // Loop through left cell neighbors (first 4 neighbors, coordinates in first 8 components)
      for (CeedInt n = 0; n < 4; n++) {
        CeedScalar nx = neighbor_coords[2*n][i];     // x coordinate of neighbor n
        CeedScalar ny = neighbor_coords[2*n+1][i];   // y coordinate of neighbor n
        
        // Check if this is a valid neighbor (zeros indicate no neighbor)
        if (nx == 0.0 && ny == 0.0 && n > 0) break;  // First neighbor might be at origin
        
        // Get neighbor solution values
        CeedScalar qN[3] = {
          neighbor_values[3*n][i],     // h of neighbor n
          neighbor_values[3*n+1][i],   // hu of neighbor n
          neighbor_values[3*n+2][i]    // hv of neighbor n
        };
        
        // Distance from left cell center to neighbor center
        CeedScalar dx = nx - xl;
        CeedScalar dy = ny - yl;
        
        // Only use if neighbor is actually different from center
        if (fabs(dx) > 1e-12 || fabs(dy) > 1e-12) {
          // Add to least squares system
          A[0][0] += dx * dx;
          A[0][1] += dx * dy;
          A[1][0] += dy * dx;
          A[1][1] += dy * dy;
          
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar dq = qN[comp] - qL[comp];
            b[comp][0] += dq * dx;
            b[comp][1] += dq * dy;
          }
          valid_neighbors++;
        }
      }
      
      // Solve least squares system and reconstruct (reduced requirement to 1 neighbor for now)
      if (valid_neighbors >= 1) {
        CeedScalar det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (fabs(det) > 1e-12) {
          CeedScalar inv_A[2][2] = {
            { A[1][1] / det, -A[0][1] / det},
            {-A[1][0] / det,  A[0][0] / det}
          };
          
          // Compute gradients and reconstruct at edge midpoint
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar grad_x = inv_A[0][0] * b[comp][0] + inv_A[0][1] * b[comp][1];
            CeedScalar grad_y = inv_A[1][0] * b[comp][0] + inv_A[1][1] * b[comp][1];
            
            CeedScalar dx_edge = edge_mid_x - xl;
            CeedScalar dy_edge = edge_mid_y - yl;
            
            // Reconstruct value at edge midpoint
            qL_recon[comp] = qL[comp] + grad_x * dx_edge + grad_y * dy_edge;
          }
          
          // Apply MinMod slope limiting
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar delta_recon = qL_recon[comp] - qL[comp];
            CeedScalar limited_delta = delta_recon;  // Start with unlimited
            
            // Check against all valid neighbors for MinMod limiting
            for (CeedInt n = 0; n < valid_neighbors; n++) {
              CeedScalar qN = neighbor_values[3*n + comp][i];
              CeedScalar delta_neighbor = qN - qL[comp];
              
              // MinMod limiting: if signs differ, limit to zero
              if (delta_recon * delta_neighbor <= 0.0) {
                limited_delta = 0.0;
                break;
              } else {
                // Same signs - take minimum magnitude
                if (fabs(delta_neighbor) < fabs(limited_delta)) {
                  limited_delta = delta_neighbor;
                }
              }
            }
            qL_recon[comp] = qL[comp] + limited_delta;
          }
          
          // Ensure positive depth
          if (qL_recon[0] < tiny_h) {
            qL_recon[0] = tiny_h;
            qL_recon[1] = 0.0;
            qL_recon[2] = 0.0;
          }
        }
      }
    }
    
    // **RIGHT CELL RECONSTRUCTION**
    if (qR[0] > tiny_h) {
      CeedScalar A[2][2] = {{0}}, b[3][2] = {{0}};
      CeedInt valid_neighbors = 0;
      
      // Loop through right cell neighbors (last 4 neighbors, offset by 8 in coords, 12 in values)
      for (CeedInt n = 0; n < 4; n++) {
        CeedScalar nx = neighbor_coords[8 + 2*n][i];     // x coordinate of neighbor n
        CeedScalar ny = neighbor_coords[8 + 2*n+1][i];   // y coordinate of neighbor n
        
        if (nx == 0.0 && ny == 0.0 && n > 0) break;
        
        CeedScalar qN[3] = {
          neighbor_values[12 + 3*n][i],     // h of neighbor n
          neighbor_values[12 + 3*n+1][i],   // hu of neighbor n
          neighbor_values[12 + 3*n+2][i]    // hv of neighbor n
        };
        
        CeedScalar dx = nx - xr;
        CeedScalar dy = ny - yr;
        
        if (fabs(dx) > 1e-12 || fabs(dy) > 1e-12) {
          A[0][0] += dx * dx;
          A[0][1] += dx * dy;
          A[1][0] += dy * dx;
          A[1][1] += dy * dy;
          
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar dq = qN[comp] - qR[comp];
            b[comp][0] += dq * dx;
            b[comp][1] += dq * dy;
          }
          valid_neighbors++;
        }
      }
      
      if (valid_neighbors >= 1) {
        CeedScalar det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
        if (fabs(det) > 1e-12) {
          CeedScalar inv_A[2][2] = {
            { A[1][1] / det, -A[0][1] / det},
            {-A[1][0] / det,  A[0][0] / det}
          };
          
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar grad_x = inv_A[0][0] * b[comp][0] + inv_A[0][1] * b[comp][1];
            CeedScalar grad_y = inv_A[1][0] * b[comp][0] + inv_A[1][1] * b[comp][1];
            
            CeedScalar dx_edge = edge_mid_x - xr;
            CeedScalar dy_edge = edge_mid_y - yr;
            qR_recon[comp] = qR[comp] + grad_x * dx_edge + grad_y * dy_edge;
          }
          
          // Apply MinMod limiting for right cell
          for (CeedInt comp = 0; comp < 3; comp++) {
            CeedScalar delta_recon = qR_recon[comp] - qR[comp];
            CeedScalar limited_delta = delta_recon;
            
            for (CeedInt n = 0; n < valid_neighbors; n++) {
              CeedScalar qN = neighbor_values[12 + 3*n + comp][i];
              CeedScalar delta_neighbor = qN - qR[comp];
              
              if (delta_recon * delta_neighbor <= 0.0) {
                limited_delta = 0.0;
                break;
              } else {
                if (fabs(delta_neighbor) < fabs(limited_delta)) {
                  limited_delta = delta_neighbor;
                }
              }
            }
            qR_recon[comp] = qR[comp] + limited_delta;
          }
          
          if (qR_recon[0] < tiny_h) {
            qR_recon[0] = tiny_h;
            qR_recon[1] = 0.0;
            qR_recon[2] = 0.0;
          }
        }
      }
    }
    
    // Compute flux using reconstructed values
    CeedScalar f[3], amax;
    if (qL_recon[0] > tiny_h || qR_recon[0] > tiny_h) {
      SWERiemannFlux_Roe(gravity, tiny_h, h_anuga,
                        (SWEState){qL_recon[0], qL_recon[1], qL_recon[2]},
                        (SWEState){qR_recon[0], qR_recon[1], qR_recon[2]},
                        sn, cn, f, &amax);
      
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i] = f[j] * (-edge_length / Al);
        cell_R[j][i] = f[j] * (edge_length / Ar);
        flux[j][i] = f[j];
      }
      courant[0][i] = -amax * edge_length / Al * dt;
      courant[1][i] = amax * edge_length / Ar * dt;
    } else {
      // Both cells are dry
      for (CeedInt j = 0; j < 3; j++) {
        cell_L[j][i] = 0.0;
        cell_R[j][i] = 0.0;
        flux[j][i] = 0.0;
      }
      courant[0][i] = 0.0;
      courant[1][i] = 0.0;
    }
  }
  
  return 0;
}

// Entry point for Roe flux with full reconstruction
CEED_QFUNCTION(SWEFluxReconstruction_Roe)(void *ctx, CeedInt Q, const CeedScalar *const in[], CeedScalar *const out[]) {
  return SWEFluxReconstructionKernel(ctx, Q, in, out);
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
    if (qL.h > tiny_h || qR.h > tiny_h) {
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
