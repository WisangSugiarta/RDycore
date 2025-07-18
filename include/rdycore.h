// RDycore -- a compound flooding dynamics simulator

#ifndef RDYCORE_H
#define RDYCORE_H

#include <petsc.h>
#include <private/config.h>

//-------------------
// RDycore Interface
//-------------------

typedef struct _p_RDy *RDy;

// Version information
PETSC_EXTERN PetscErrorCode RDyGetVersion(int *, int *, int *, PetscBool *);
PETSC_EXTERN PetscErrorCode RDyGetBuildConfiguration(const char **);

// Process initialization
PETSC_EXTERN PetscErrorCode RDyInit(int, char *[], const char *);
PETSC_EXTERN PetscErrorCode RDyInitFortran(void);
PETSC_EXTERN PetscErrorCode RDyOnFinalize(void (*)(void));
PETSC_EXTERN PetscInt       RDyFinalize(void);
PETSC_EXTERN PetscBool      RDyInitialized(void);

// RDycore online configuration
PETSC_EXTERN PetscErrorCode RDySetLogFile(RDy, const char *);

// RDycore setup/breakdown
PETSC_EXTERN PetscErrorCode RDyCreate(MPI_Comm, const char *, RDy *);
PETSC_EXTERN PetscErrorCode RDySetup(RDy);
PETSC_EXTERN PetscErrorCode RDyDestroy(RDy *);

// RDycore support for Method of Manufactured Solutions (MMS)
PETSC_EXTERN PetscErrorCode RDyMMSSetup(RDy);
PETSC_EXTERN PetscErrorCode RDyMMSComputeSolution(RDy, PetscReal, Vec);
PETSC_EXTERN PetscErrorCode RDyMMSEnforceBoundaryConditions(RDy, PetscReal);
PETSC_EXTERN PetscErrorCode RDyMMSComputeSourceTerms(RDy, PetscReal);
PETSC_EXTERN PetscErrorCode RDyMMSUpdateMaterialProperties(RDy);
PETSC_EXTERN PetscErrorCode RDyMMSComputeErrorNorms(RDy, PetscReal, PetscReal *, PetscReal *, PetscReal *, PetscInt *, PetscReal *);
PETSC_EXTERN PetscErrorCode RDyMMSEstimateConvergenceRates(RDy, PetscReal *, PetscReal *, PetscReal *);
PETSC_EXTERN PetscErrorCode RDyMMSRun(RDy);

// RDycore support for AMR
PETSC_EXTERN PetscErrorCode RDyRefine(RDy);

// time integration
PETSC_EXTERN PetscErrorCode DestroyOutputViewer(RDy);
PETSC_EXTERN PetscErrorCode RDyAdvance(RDy);

// Accessing data
PETSC_EXTERN PetscBool RDyFinished(RDy);
PETSC_EXTERN PetscBool RDyRestarted(RDy);

// time units
typedef enum { RDY_TIME_UNSET = 0, RDY_TIME_SECONDS, RDY_TIME_MINUTES, RDY_TIME_HOURS, RDY_TIME_DAYS, RDY_TIME_MONTHS, RDY_TIME_YEARS } RDyTimeUnit;

PETSC_EXTERN PetscErrorCode RDyGetTimeUnit(RDy, RDyTimeUnit *);
PETSC_EXTERN PetscErrorCode RDyGetTime(RDy, RDyTimeUnit, PetscReal *);
PETSC_EXTERN PetscErrorCode RDyGetTimeStep(RDy, RDyTimeUnit, PetscReal *);
PETSC_EXTERN PetscErrorCode RDyConvertTime(RDyTimeUnit, PetscReal, RDyTimeUnit, PetscReal *);
PETSC_EXTERN PetscErrorCode RDyGetStep(RDy, PetscInt *);
PETSC_EXTERN PetscErrorCode RDyGetCouplingInterval(RDy, RDyTimeUnit, PetscReal *);
PETSC_EXTERN PetscErrorCode RDySetCouplingInterval(RDy, RDyTimeUnit, PetscReal);

PETSC_EXTERN PetscErrorCode RDyGetNumGlobalCells(RDy, PetscInt *);
PETSC_EXTERN PetscErrorCode RDyGetNumLocalCells(RDy, PetscInt *);
PETSC_EXTERN PetscErrorCode RDyGetNumBoundaryConditions(RDy, PetscInt *);
PETSC_EXTERN PetscErrorCode RDyGetNumBoundaryEdges(RDy, const PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryConditionFlowType(RDy, const PetscInt, PetscInt *);

PETSC_EXTERN PetscErrorCode RDyGetLocalCellHeights(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellXMomenta(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellYMomenta(RDy rdy, const PetscInt size, PetscReal *values);

PETSC_EXTERN PetscErrorCode RDyGetLocalCellXCentroids(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellYCentroids(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellZCentroids(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellAreas(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetLocalCellNaturalIDs(RDy rdy, const PetscInt size, PetscInt *values);

PETSC_EXTERN PetscErrorCode RDyGetBoundaryEdgeXCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryEdgeYCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryEdgeZCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryCellXCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryCellYCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDyGetBoundaryCellZCentroids(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscReal *values);

PETSC_EXTERN PetscErrorCode RDyGetBoundaryID(RDy rdy, const PetscInt boundary_index, PetscInt *value);

PETSC_EXTERN PetscErrorCode RDyGetBoundaryCellNaturalIDs(RDy rdy, const PetscInt boundary_index, const PetscInt size, PetscInt *values);

PETSC_EXTERN PetscErrorCode RDySetFlowDirichletBoundaryValues(RDy rdy, const PetscInt boundary_index, const PetscInt num_edges, const PetscInt ndof,
                                                              PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetSedimentDirichletBoundaryValues(RDy rdy, const PetscInt boundary_index, const PetscInt num_edges,
                                                                  const PetscInt ndof, PetscReal *values);

PETSC_EXTERN PetscErrorCode RDySetRegionalWaterSource(RDy rdy, const PetscInt region_idx, PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetRegionalXMomentumSource(RDy rdy, const PetscInt region_idx, PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetRegionalYMomentumSource(RDy rdy, const PetscInt region_idx, PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetRegionalSedimentSource(RDy rdy, const PetscInt region_idx, PetscInt sediment_class, PetscInt size,
                                                         PetscReal *values);

PETSC_EXTERN PetscErrorCode RDySetHomogeneousRegionalWaterSource(RDy rdy, const PetscInt region_idx, PetscReal value);
PETSC_EXTERN PetscErrorCode RDySetHomogeneousRegionalXMomentumSource(RDy rdy, const PetscInt region_idx, PetscReal value);
PETSC_EXTERN PetscErrorCode RDySetHomogeneousRegionalYMomentumSource(RDy rdy, const PetscInt region_idx, PetscReal value);

PETSC_EXTERN PetscErrorCode RDySetDomainWaterSource(RDy rdy, PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetDomainXMomentumSource(RDy rdy, PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetDomainYMomentumSource(RDy rdy, PetscInt size, PetscReal *values);

PETSC_EXTERN PetscErrorCode RDySetRegionalManningsN(RDy rdy, const PetscInt region_idx, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetDomainManningsN(RDy rdy, const PetscInt size, PetscReal *values);
PETSC_EXTERN PetscErrorCode RDySetInitialConditions(RDy rdy, Vec ic);

PETSC_EXTERN PetscErrorCode RDyCreatePrognosticVec(RDy rdy, Vec *prog_vec);
PETSC_EXTERN PetscErrorCode RDyReadOneDOFLocalVecFromBinaryFile(RDy rdy, const char *, Vec *local_vec);
PETSC_EXTERN PetscErrorCode RDyReadOneDOFGlobalVecFromBinaryFile(RDy rdy, const char *, Vec *local_vec);
PETSC_EXTERN PetscErrorCode RDyWriteOneDOFGlobalVecToBinaryFile(RDy rdy, const char *, Vec *global);
PETSC_EXTERN PetscErrorCode RDyCreateOneDOFGlobalVec(RDy rdy, Vec *global);

// "kinds" of initial/boundary/source conditions applied to regions/boundaries
typedef enum {
  CONDITION_DIRICHLET = 0,     // Dirichlet condition (value is specified)
  CONDITION_NEUMANN,           // Neumann condition (derivative is specified)
  CONDITION_REFLECTING,        // Reflecting condition
  CONDITION_CRITICAL_OUTFLOW,  // Critical flow
  CONDITION_RUNOFF             // Runoff that is a source to the dh/dt equation
} RDyConditionType;

#endif
