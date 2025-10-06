# 1 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/rdycore.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/wsugiarta/Documents/github/RDycore/build_test//"
# 1 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/rdycore.F90"
! This file defines the Fortran interface for RDycore.

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h" 1
!
!  Manually maintained part of the base include file for Fortran use of PETSc.
!  Note: This file should contain only define statements
!



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petscconf.h" 1



# 8 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h" 2




# 1 "/Users/wsugiarta/Documents/github/petsc/include/petscversion.h" 1















































# 12 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h" 2

!
! The real*8,complex*16 notatiton is used so that the
! PETSc double/complex variables are not affected by
! compiler options like -r4,-r8, that are sometimes invoked
! by the user. NAG compiler does not like integer*4,real*8






























!
! Fortran does not support unsigned, though ISO_C_BINDING
! supports INTEGER(KIND=C_SIZE_T). We don't use that here
! only to avoid importing the module.





!


!


!
# 74 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h"

# 96 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h"
!
!     Macro for templating between real and complex
!
# 119 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h"













!
!     Macros for error checking
!
# 155 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h"

# 165 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscsysbase.h"





# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmacros.h" 1







# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmath.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 1
# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmath.h" 2







# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscerror.h" 1












# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscviewer.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdraw.h" 1


















# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscviewer.h" 2











# 9 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscoptions.h" 1














# 10 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclog.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctime.h" 1








# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclog.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbt.h" 1






# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscstring.h" 1








# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbt.h" 2




# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclog.h" 2

# 16 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclog.h"







# 11 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2


# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdevice.h" 1





# 15 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdevice.h"





# 13 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h" 2

# 37 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h"



# 56 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsys.h"

# 3 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/rdycore.F90" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscvecbase.h" 1
!
!  Part of the base include file for Fortran use of PETSc Vec
!  Note: This file should contain only define statements

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!










# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 2


# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 1
# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h" 2

# 18 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h"






# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscisbase.h" 1
!
!  Part of the base include file for Fortran use of PETSc IS
!  Note: This file should contain only define statements

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!





# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 2



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsection.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 1
# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsection.h" 2








# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 2















# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 2



















# 4 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/rdycore.F90" 2


module rdycore

  use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_loc
  use petscsys

  implicit none

  public :: RDyDouble, RDy, RDyGetVersion, RDyGetBuildConfiguration, RDyInit, RDyFinalize, &
            RDyInitialized, RDyCreate, RDySetup, RDyAdvance, RDyDestroy, &
            RDyMMSSetup, RDyMMSComputeSolution, RDyMMSEnforceBoundaryConditions, &
            RDyMMSComputeSourceTerms, RDyMMSUpdateMaterialProperties, &
            RDyMMSComputeErrorNorms, RDyMMSEstimateConvergenceRates, RDyMMSRun, &
            RDyGetNumGlobalCells, RDyGetNumLocalCells, RDyGetNumBoundaryConditions, &
            RDyGetNumBoundaryEdges, RDyGetBoundaryConditionFlowType, &
            RDySetFlowDirichletBoundaryValues, &
            RDyGetTimeUnit, RDyGetTime, RDyGetTimeStep, RDyConvertTime, &
            RDyGetStep, RDyGetCouplingInterval, RDySetCouplingInterval, &
            RDyGetLocalCellHeights, RDyGetLocalCellXMomenta, RDyGetLocalCellYMomenta, &
            RDyGetLocalCellXCentroids, RDyGetLocalCellYCentroids, RDyGetLocalCellZCentroids, &
            RDyGetLocalCellAreas, RDyGetLocalCellNaturalIDs, &
            RDyGetBoundaryEdgeXCentroids, RDyGetBoundaryEdgeYCentroids, RDyGetBoundaryEdgeZCentroids, &
            RDyGetBoundaryCellNaturalIDs, &
            RDySetDomainWaterSource, RDySetRegionalWaterSource, RDySetDomainXMomentumSource, &
            RDySetDomainYMomentumSource, RDySetDomainManningsN, RDySetInitialConditions, &
            RDyCreatePrognosticVec, RDyReadOneDOFLocalVecFromBinaryFile, RDyReadOneDOFGlobalVecFromBinaryFile, &
            RDyCreateOneDOFGlobalVec, RDyWriteOneDOFGlobalVecToBinaryFile

  ! RDycore uses double-precision floating point numbers
  integer, parameter :: RDyDouble = selected_real_kind(12)

  ! Supported time units (must be synchronized with RDyTimeUnit in rdycore.h)
  integer, parameter :: RDyTimeUnit  = c_int
  integer, parameter :: RDY_TIME_UNSET   = 0
  integer, parameter :: RDY_TIME_SECONDS = 1
  integer, parameter :: RDY_TIME_MINUTES = 2
  integer, parameter :: RDY_TIME_HOURS   = 3
  integer, parameter :: RDY_TIME_DAYS    = 4
  integer, parameter :: RDY_TIME_MONTHS  = 5
  integer, parameter :: RDY_TIME_YEARS   = 6

  type :: RDy
    ! C pointer to RDy type
    type(c_ptr) :: c_rdy
  end type RDy

  interface
    integer(c_int) function rdygetversion_(major, minor, patch, release) bind(c, name="RDyGetVersion")
      use iso_c_binding, only: c_int
      integer(kind=selected_int_kind(5)),  intent(out) :: major, minor, patch
      logical(kind=4), intent(out) :: release
    end function

    integer(c_int) function rdygetbuildconfiguration_(c_build_config) bind(c, name="RDyGetBuildConfigurationF90")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value :: c_build_config
    end function

    integer(c_int) function rdyinitfortran_() bind(c, name="RDyInitFortran")
      use iso_c_binding, only: c_int
    end function

    logical(c_bool) function rdyinitialized_() bind(c, name="RDyInitialized")
      use iso_c_binding, only: c_bool
    end function

    integer(c_int) function rdyfinalize_() bind(c, name="RDyFinalize")
      use iso_c_binding, only: c_int
    end function

    integer(c_int) function rdycreate_(comm, filename, rdy) bind(c, name="RDyCreateF90")
      use iso_c_binding, only: c_int, c_ptr
      integer,            intent(in)  :: comm
      type(c_ptr), value, intent(in)  :: filename
      type(c_ptr),        intent(out) :: rdy
    end function

    integer(c_int) function rdysetlogfile_(rdy, filename) bind(c, name="RDySetLogFile")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      type(c_ptr), value, intent(in) :: filename
    end function

    integer(c_int) function rdysetup_(rdy) bind(c, name="RDySetup")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    integer(c_int) function rdymmssetup_(rdy) bind(c, name="RDyMMSSetup")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    integer(c_int) function rdymmscomputesolution_(rdy, time, soln) bind(c, name="RDyMMSComputeSolution")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr),    value, intent(in)   :: rdy
      real(c_double), value, intent(in)   :: time
      integer(kind=selected_int_kind(10)), value, intent(in) :: soln
    end function

    integer(c_int) function rdymmsenforceboundaryconditions_(rdy, time) bind(c, name="RDyMMSEnforceBoundaryConditions")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in)   :: rdy
      real(c_double), value, intent(in)   :: time
    end function

    integer(c_int) function rdymmscomputesourceterms_(rdy, time) bind(c, name="RDyMMSComputeSourceTerms")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in)   :: rdy
      real(c_double), value, intent(in)   :: time
    end function

    integer(c_int) function rdymmsupdatematerialproperties_(rdy) bind(c, name="RDyMMSUpdateMaterialProperties")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    integer(c_int) function rdymmscomputeerrornorms_(rdy, time, l1_norms, l2_norms, linf_norms, &
                                                     num_global_cells, global_area) bind(c, name="RDyMMSComputeErrorNorms")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in)  :: rdy
      real(c_double), value, intent(in)  :: time
      type(c_ptr),    value, intent(in)  :: l1_norms, l2_norms, linf_norms
      integer(kind=selected_int_kind(5)),              intent(out) :: num_global_cells
      real(c_double),        intent(out) :: global_area
    end function

    integer(c_int) function rdymmsestimateconvergencerates_(rdy, num_refinements, l1_rates, l2_rates, linf_rates) &
                            bind(c, name="RDyMMSEstimateConvergenceRates")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),    value, intent(in)  :: num_refinements
      type(c_ptr), value, intent(in)  :: l1_rates, l2_rates, linf_rates
    end function

    integer(c_int) function rdymmsrun_(rdy) bind(c, name="RDyMMSRun")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr), value, intent(in)  :: rdy
    end function

    integer(c_int) function rdygetnumglobalcells_(rdy, num_cells_global) bind(c, name="RDyGetNumGlobalCells")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),           intent(out) :: num_cells_global
    end function

    integer(c_int) function rdygetnumlocalcells_(rdy, num_cells) bind(c, name="RDyGetNumLocalCells")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),           intent(out) :: num_cells
    end function

    integer(c_int) function rdygetnumboundaryconditions_(rdy, num_bnd_conds) bind(c, name="RDyGetNumBoundaryConditions")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),           intent(out) :: num_bnd_conds
    end function

    integer(c_int) function rdygetnumboundaryedges_(rdy, boundary_id, num_edges) bind(c, name="RDyGetNumBoundaryEdges")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5)),    value, intent(in) :: boundary_id
      integer(kind=selected_int_kind(5)),    intent(out)       :: num_edges
    end function

    integer(c_int) function RDySetFlowDirichletBoundaryValues_(rdy, boundary_id, num_edges, ndof, bc_values) bind(c, name="RDySetFlowDirichletBoundaryValues")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),    value, intent(in)  :: boundary_id
      integer(kind=selected_int_kind(5)),    value, intent(in)  :: num_edges
      integer(kind=selected_int_kind(5)),    value, intent(in)  :: ndof
      type(c_ptr), value, intent(in)  :: bc_values
    end function

    integer(c_int) function rdygetboundaryconditionflowtype_(rdy, boundary_id, bnd_cond_type) bind(c, name="RDyGetBoundaryConditionFlowType")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr),    value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)), value, intent(in)  :: boundary_id
      integer(kind=selected_int_kind(5)),        intent(out) :: bnd_cond_type
    end function

    integer(c_int) function rdygettimeunit_(rdy, time_unit) bind(c, name="RDyGetTimeUnit")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in) :: rdy
      integer(c_int), intent(out)       :: time_unit
    end function

    integer(c_int) function rdygettime_(rdy, time_unit, time) bind(c, name="RDyGetTime")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in)  :: rdy
      integer(c_int), value, intent(in)  :: time_unit
      real(c_double),        intent(out) :: time
    end function

    integer(c_int) function rdygettimestep_(rdy, time_unit, dt) bind(c, name="RDyGetTimeStep")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in)  :: rdy
      integer(c_int), value, intent(in)  :: time_unit
      real(c_double),        intent(out) :: dt
    end function

    integer(c_int) function rdyconverttime_(unit_from, t_from, unit_to, t_to) bind(c, name="RDyConvertTime")
      use iso_c_binding, only: c_int, c_ptr, c_double
      integer(c_int), value, intent(in)  :: unit_from
      real(c_double), value, intent(in)  :: t_from
      integer(c_int), value, intent(in)  :: unit_to
      real(c_double),        intent(out) :: t_to
    end function

    integer(c_int) function rdygetstep_(rdy, step) bind(c, name="RDyGetStep")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(5)),           intent(out) :: step
    end function

    integer(c_int) function rdygetcouplinginterval_(rdy, time_unit, interval) bind(c, name="RDyGetCouplingInterval")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in) :: rdy
      integer(c_int), value, intent(in) :: time_unit
      real(c_double), intent(out)       :: interval
    end function

    integer(c_int) function rdysetcouplinginterval_(rdy, time_unit, interval) bind(c, name="RDySetCouplingInterval")
      use iso_c_binding, only: c_int, c_ptr, c_double
      type(c_ptr),    value, intent(in) :: rdy
      integer(c_int), value, intent(in) :: time_unit
      real(c_double), value, intent(in) :: interval
    end function

    integer(c_int) function rdygetlocalcellheights_(rdy, size, h) bind(c, name="RDyGetLocalCellHeights")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: h
    end function

    integer(c_int) function rdygetlocaclellxmomenta_(rdy, size, hu) bind(c, name="RDyGetLocalCellXMomenta")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5)), value,    intent(in) :: size
      type(c_ptr), value, intent(in) :: hu
    end function

    integer(c_int) function rdygetlocaclellymomenta_(rdy, size, hv) bind(c, name="RDyGetLocalCellYMomenta")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5)), value,    intent(in) :: size
      type(c_ptr), value, intent(in) :: hv
    end function

    integer(c_int) function rdygetlocalcellxcentroids_(rdy, size, values) bind(c, name="RDyGetLocalCellXCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetlocalcellycentroids_(rdy, size, values) bind(c, name="RDyGetLocalCellYCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetlocalcellzcentroids_(rdy, size, values) bind(c, name="RDyGetLocalCellZCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetlocalcellareas_(rdy, size, values) bind(c, name="RDyGetLocalCellAreas")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetlocalcellnaturalids_(rdy, size, values) bind(c, name="RDyGetLocalCellNaturalIDs")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundaryedgexcentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryEdgeXCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundaryedgeycentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryEdgeYCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundaryedgezcentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryEdgeZCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundarycellxcentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryCellXCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundarycellycentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryCellYCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundarycellzcentroids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryCellZCentroids")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdygetboundarycellnaturalids_(rdy, boundary_index, size, values) bind(c, name="RDyGetBoundaryCellNaturalIDs")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: boundary_index
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: values
    end function

    integer(c_int) function rdysetdomainwatersource_(rdy, size, watsrc) bind(c, name="RDySetDomainWaterSource")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: watsrc
    end function

    integer(c_int) function rdysetregionalwatersource_(rdy, region_idx, size, watsrc) bind(c, name="RDySetRegionalWaterSource")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: region_idx
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: watsrc
    end function

    integer(c_int) function rdysetdomainxmomentumsource_(rdy, size, xmomsrc) bind(c, name="RDySetDomainXMomentumSource")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: xmomsrc
    end function

    integer(c_int) function rdysetdomainymomentumsource_(rdy, size, ymomsrc) bind(c, name="RDySetDomainYMomentumSource")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: ymomsrc
    end function

    integer(c_int) function rdysetdomainmanningsn_(rdy, size, n) bind(c, name="RDySetDomainManningsN")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
      integer(kind=selected_int_kind(5))   , value, intent(in) :: size
      type(c_ptr), value, intent(in) :: n
    end function

    integer(c_int) function rdysetinitialconditions_(rdy, ic) bind(c, name="RDySetInitialConditions")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr),      value, intent(in) :: rdy
      integer(kind=selected_int_kind(10)), value, intent(in) :: ic
    end function

    integer(c_int) function rdycreateprognosticvec_(rdy, prog_vec) bind(c, name="RDyCreatePrognosticVec")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(10)),   intent(out) :: prog_vec
    end function

    integer(c_int) function rdycreateonedofglobalvec_(rdy, global_vec) bind(c, name="RDyCreateOneDOFGlobalVec")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr), value, intent(in)  :: rdy
      integer(kind=selected_int_kind(10)),   intent(out) :: global_vec
    end function

    integer(c_int) function rdyreadonedoflocalvecfrombinaryfile_(rdy, filename, local_vec) bind(c, name="RDyReadOneDOFLocalVecFromBinaryFile")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr), value, intent(in)  :: rdy
      type(c_ptr), value, intent(in)  :: filename
      integer(kind=selected_int_kind(10)),   intent(out) :: local_vec
    end function

    integer(c_int) function rdyreadonedofglobalvecfrombinaryfile_(rdy, filename, global_vec) bind(c, name="RDyReadOneDOFGlobalVecFromBinaryFile")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr), value, intent(in)  :: rdy
      type(c_ptr), value, intent(in)  :: filename
      integer(kind=selected_int_kind(10)),   intent(out) :: global_vec
    end function

    integer(c_int) function rdywriteonedofglobalvectobinaryfile_(rdy, filename, global_vec) bind(c, name="RDyWriteOneDOFGlobalVecToBinaryFile")
      use iso_c_binding, only: c_int, c_ptr
      use petscvec
      type(c_ptr), value, intent(in)  :: rdy
      type(c_ptr), value, intent(in)  :: filename
      integer(kind=selected_int_kind(10)),   intent(out) :: global_vec
    end function

    integer(c_int) function rdyadvance_(rdy) bind(c, name="RDyAdvance")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    logical(c_bool) function rdyfinished_(rdy) bind(c, name="RDyFinished")
      use iso_c_binding, only: c_bool, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    logical(c_bool) function rdyrestarted_(rdy) bind(c, name="RDyRestarted")
      use iso_c_binding, only: c_bool, c_ptr
      type(c_ptr), value, intent(in) :: rdy
    end function

    integer(c_int) function rdydestroy_(rdy) bind(c, name="RDyDestroy")
      use iso_c_binding, only: c_int, c_ptr
      type(c_ptr), intent(inout) :: rdy
    end function

  end interface

contains

  subroutine RDyGetVersion(major, minor, patch, release, ierr)
    use petscsys
    integer(kind=selected_int_kind(5)),  intent(out) :: major, minor, patch
    logical(kind=4), intent(out) :: release
    integer,   intent(out) :: ierr
  end subroutine

  subroutine RDyGetBuildConfiguration(build_config, ierr)
    use petscsys
    use iso_c_binding
    character(len=1024), target, intent(out) :: build_config
    integer,                     intent(out) :: ierr

    ierr = rdygetbuildconfiguration_(c_loc(build_config))
  end subroutine

  subroutine RDyInit(ierr)
    use petscsys
    implicit none
    integer, intent(out) :: ierr
    call PetscInitialize(ierr)
    if ((ierr == 0) .and. .not. RDyInitialized()) ierr = rdyinitfortran_()
  end subroutine

  subroutine RDyFinalize(ierr)
    integer, intent(out) :: ierr
    ierr = rdyfinalize_()
  end subroutine

  logical function RDyInitialized()
    RDyInitialized = rdyinitialized_()
  end function

  subroutine RDyCreate(comm, filename, rdy_, ierr)
    use iso_c_binding, only: c_null_char
    character(len=1024), intent(in) :: filename
    integer,   intent(in)  :: comm
    type(RDy), intent(out) :: rdy_
    integer,   intent(out) :: ierr

    integer                      :: n
    character(len=1024), pointer :: config_file

    n = len_trim(filename)
    allocate(config_file)
    config_file(1:n) = filename(1:n)
    config_file(n+1:n+1) = c_null_char
    ierr = rdycreate_(comm, c_loc(config_file), rdy_%c_rdy)
    deallocate(config_file)
  end subroutine

  subroutine RDySetLogFile(rdy_, filename, ierr)
    use iso_c_binding, only: c_null_char
    character(len=1024), intent(in)  :: filename
    type(RDy),           intent(out) :: rdy_
    integer,             intent(out) :: ierr

    integer                      :: n
    character(len=1024), pointer :: log_file

    n = min(len_trim(filename), 1024)
    allocate(log_file)
    log_file(1:n) = filename(1:n)
    log_file(n+1:n+1) = c_null_char
    ierr = rdysetlogfile_(rdy_%c_rdy, c_loc(log_file))
    deallocate(log_file)
  end subroutine

  subroutine RDySetup(rdy_, ierr)
    type(RDy), intent(inout) :: rdy_
    integer,   intent(out)   :: ierr
    ierr = rdysetup_(rdy_%c_rdy);if (ierr .ne. 0) then;call PetscErrorF(ierr,523,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/rdycore.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
  end subroutine

  subroutine RDyMMSSetup(rdy_, ierr)
    type(RDy), intent(inout) :: rdy_
    integer,   intent(out)   :: ierr
    ierr = rdymmssetup_(rdy_%c_rdy)
  end subroutine

  subroutine RDyMMSComputeSolution(rdy_, time, solution, ierr)
    use petscvec
    type(RDy),       intent(inout) :: rdy_
    real(RDyDouble), intent(in)    :: time
    type(tVec),      intent(in)    :: solution  ! type(tVec)
    integer,         intent(out)   :: ierr
    ierr = rdymmscomputesolution_(rdy_%c_rdy, time, solution%v)
  end subroutine

  subroutine RDyMMSEnforceBoundaryConditions(rdy_, time, ierr)
    type(RDy),       intent(inout) :: rdy_
    real(RDyDouble), intent(in)    :: time
    integer,         intent(out)   :: ierr
    ierr = rdymmsenforceboundaryconditions_(rdy_%c_rdy, time)
  end subroutine

  subroutine RDyMMSComputeSourceTerms(rdy_, time, ierr)
    type(RDy),       intent(inout) :: rdy_
    real(RDyDouble), intent(in)    :: time
    integer,         intent(out)   :: ierr
    ierr = rdymmscomputesourceterms_(rdy_%c_rdy, time)
  end subroutine

  subroutine RDyMMSUpdateMaterialProperties(rdy_, ierr)
    type(RDy), intent(inout) :: rdy_
    integer,   intent(out)   :: ierr
    ierr = rdymmsupdatematerialproperties_(rdy_%c_rdy)
  end subroutine

  subroutine RDyMMSComputeErrorNorms(rdy_, time, l1_norms, l2_norms, linf_norms, &
                                     num_global_cells, global_area, ierr)
    type(RDy),                intent(inout) :: rdy_
    real(RDyDouble),          intent(in)    :: time
    real(RDyDouble), pointer, intent(in)    :: l1_norms(:), l2_norms(:), linf_norms(:)
    integer(kind=selected_int_kind(5)),                 intent(out)   :: num_global_cells
    real(RDyDouble),          intent(out)   :: global_area
    integer,                  intent(out)   :: ierr
    ierr = rdymmscomputeerrornorms_(rdy_%c_rdy, time, c_loc(l1_norms), c_loc(l2_norms), c_loc(linf_norms), &
                                    num_global_cells, global_area)
  end subroutine

  subroutine RDyMMSEstimateConvergenceRates(rdy_, num_refinements, &
                                            l1_conv_rates, l2_conv_rates, linf_conv_rates, ierr)
    type(RDy),                intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),                 intent(in)    :: num_refinements
    real(RDyDouble), pointer, intent(in)    :: l1_conv_rates(:), l2_conv_rates(:), linf_conv_rates(:)
    integer,                  intent(out)   :: ierr
    ierr = rdymmsestimateconvergencerates_(rdy_%c_rdy, num_refinements, &
                                           c_loc(l1_conv_rates), c_loc(l2_conv_rates), c_loc(linf_conv_rates))
  end subroutine

  subroutine RDyMMSRun(rdy_, ierr)
    type(RDy),                intent(inout) :: rdy_
    integer,                  intent(out)   :: ierr
    ierr = rdymmsrun_(rdy_%c_rdy)
  end subroutine

  subroutine RDyGetNumGlobalCells(rdy_, num_cells_global, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(out)   :: num_cells_global
    integer,   intent(out)   :: ierr
    ierr = rdygetnumglobalcells_(rdy_%c_rdy, num_cells_global)
  end subroutine

  subroutine RDyGetNumLocalCells(rdy_, num_cells, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(out)   :: num_cells
    integer,   intent(out)   :: ierr
    ierr = rdygetnumlocalcells_(rdy_%c_rdy, num_cells)
  end subroutine

  subroutine RDyGetNumBoundaryConditions(rdy_, num_bnd_conds, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(out)   :: num_bnd_conds
    integer,   intent(out)   :: ierr
    ierr = rdygetnumboundaryconditions_(rdy_%c_rdy, num_bnd_conds)
  end subroutine

  subroutine RDyGetNumBoundaryEdges(rdy_, boundary_id, num_edges, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(in)    :: boundary_id
    integer(kind=selected_int_kind(5)),  intent(out)   :: num_edges
    integer,   intent(out)   :: ierr
    ierr = rdygetnumboundaryedges_(rdy_%c_rdy, boundary_id-1, num_edges)
  end subroutine

  subroutine RDySetFlowDirichletBoundaryValues(rdy_, boundary_id, num_edges, ndof, bc_values, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: boundary_id
    integer(kind=selected_int_kind(5)),        intent(in)          :: num_edges
    integer(kind=selected_int_kind(5)),        intent(in)          :: ndof
    real(RDyDouble), pointer, intent(in) :: bc_values(:)
    integer,         intent(out)         :: ierr
    ierr = RDySetFlowDirichletBoundaryValues_(rdy_%c_rdy, boundary_id-1, num_edges, ndof, c_loc(bc_values))
  end subroutine

  subroutine RDyGetBoundaryConditionFlowType(rdy_, boundary_id, bnd_cond_type, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(in)    :: boundary_id
    integer(kind=selected_int_kind(5)),  intent(out)   :: bnd_cond_type
    integer,   intent(out)   :: ierr
    ierr = rdygetboundaryconditionflowtype_(rdy_%c_rdy, boundary_id-1, bnd_cond_type)
  end subroutine

  subroutine RDyGetTimeUnit(rdy_, time_unit, ierr)
    type(RDy),            intent(inout) :: rdy_
    integer(RDyTimeUnit), intent(out)   :: time_unit
    integer,              intent(out)   :: ierr
    ierr = rdygettimeunit_(rdy_%c_rdy, time_unit)
  end subroutine

  subroutine RDyGetTime(rdy_, time_unit, time, ierr)
    type(RDy),            intent(inout) :: rdy_
    integer(RDyTimeUnit), intent(in)    :: time_unit
    real(RDyDouble),      intent(out)   :: time
    integer,              intent(out)   :: ierr
    ierr = rdygettime_(rdy_%c_rdy, time_unit, time)
  end subroutine

  subroutine RDyGetTimeStep(rdy_, time_unit, timestep, ierr)
    type(RDy),            intent(inout) :: rdy_
    integer(RDyTimeUnit), intent(in)    :: time_unit
    real(RDyDouble),      intent(out)   :: timestep
    integer,              intent(out)   :: ierr
    ierr = rdygettimestep_(rdy_%c_rdy, time_unit, timestep)
  end subroutine

  subroutine RDyConvertTime(unit_from, t_from, unit_to, t_to, ierr)
    integer(RDyTimeUnit), intent(in)    :: unit_from
    real(RDyDouble),      intent(in)    :: t_from
    integer(RDyTimeUnit), intent(in)    :: unit_to
    real(RDyDouble),      intent(out)   :: t_to
    integer,              intent(out)   :: ierr
    ierr = rdyconverttime_(unit_from, t_from, unit_to, t_to)
  end subroutine

  subroutine RDyGetCouplingInterval(rdy_, time_unit, interval, ierr)
    type(RDy),            intent(inout) :: rdy_
    integer(RDyTimeUnit), intent(in)    :: time_unit
    real(RDyDouble),      intent(out)   :: interval
    integer,              intent(out)   :: ierr
    ierr = rdygetcouplinginterval_(rdy_%c_rdy, time_unit, interval)
  end subroutine

  subroutine RDySetCouplingInterval(rdy_, time_unit, interval, ierr)
    type(RDy),            intent(inout) :: rdy_
    integer(RDyTimeUnit), intent(in)    :: time_unit
    real(RDyDouble),      intent(in)    :: interval
    integer,              intent(out)   :: ierr
    ierr = rdysetcouplinginterval_(rdy_%c_rdy, time_unit, interval)
  end subroutine

  subroutine RDyGetStep(rdy_, step, ierr)
    type(RDy), intent(inout) :: rdy_
    integer(kind=selected_int_kind(5)),  intent(out)   :: step
    integer,   intent(out)   :: ierr
    ierr = rdygetstep_(rdy_%c_rdy, step)
  end subroutine

  subroutine RDyGetLocalCellHeights(rdy_, size, h, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: h(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellheights_(rdy_%c_rdy, size, c_loc(h))
  end subroutine

  subroutine RDyGetLocalCellXMomenta(rdy_, size, hu, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: hu(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocaclellxmomenta_(rdy_%c_rdy, size, c_loc(hu))
  end subroutine

  subroutine RDyGetLocalCellYMomenta(rdy_, size, hv, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: hv(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocaclellymomenta_(rdy_%c_rdy, size, c_loc(hv))
  end subroutine

  subroutine RDyGetLocalCellXCentroids(rdy_, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellxcentroids_(rdy_%c_rdy, size, c_loc(values))
  end subroutine

  subroutine RDyGetLocalCellYCentroids(rdy_, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellycentroids_(rdy_%c_rdy, size, c_loc(values))
  end subroutine

  subroutine RDyGetLocalCellZCentroids(rdy_, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellzcentroids_(rdy_%c_rdy, size, c_loc(values))
  end subroutine

  subroutine RDyGetLocalCellAreas(rdy_, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellareas_(rdy_%c_rdy, size, c_loc(values))
  end subroutine

  subroutine RDyGetLocalCellNaturalIDs(rdy_, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    integer(kind=selected_int_kind(5)),        pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetlocalcellnaturalids_(rdy_%c_rdy, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryEdgeXCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundaryedgexcentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryEdgeYCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundaryedgeycentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryEdgeZCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundaryedgezcentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryCellXCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundarycellxcentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryCellYCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundarycellycentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryCellZCentroids(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    real(RDyDouble), pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundarycellzcentroids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDyGetBoundaryCellNaturalIDs(rdy_, boundary_index, size, values, ierr)
    type(RDy),       intent(inout)          :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)             :: size
    integer(kind=selected_int_kind(5)),        intent(in)             :: boundary_index
    integer(kind=selected_int_kind(5)),        pointer, intent(inout) :: values(:)
    integer,         intent(out)            :: ierr
    ierr = rdygetboundarycellnaturalids_(rdy_%c_rdy, boundary_index - 1, size, c_loc(values))
  end subroutine

  subroutine RDySetDomainWaterSource(rdy_, size, watsrc, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: size
    real(RDyDouble), pointer, intent(in) :: watsrc(:)
    integer,         intent(out)         :: ierr
    ierr = rdysetdomainwatersource_(rdy_%c_rdy, size, c_loc(watsrc))
  end subroutine

  subroutine RDySetRegionalWaterSource(rdy_, region_idx, size, watsrc, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: region_idx
    integer(kind=selected_int_kind(5)),        intent(in)          :: size
    real(RDyDouble), pointer, intent(in) :: watsrc(:)
    integer,         intent(out)         :: ierr
    ierr = rdysetregionalwatersource_(rdy_%c_rdy, region_idx, size, c_loc(watsrc))
  end subroutine

  subroutine RDySetDomainXMomentumSource(rdy_, size, xmomsrc, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: size
    real(RDyDouble), pointer, intent(in) :: xmomsrc(:)
    integer,         intent(out)         :: ierr
    ierr = rdysetdomainxmomentumsource_(rdy_%c_rdy, size, c_loc(xmomsrc))
  end subroutine

  subroutine RDySetDomainYMomentumSource(rdy_, size, ymomsrc, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: size
    real(RDyDouble), pointer, intent(in) :: ymomsrc(:)
    integer,         intent(out)         :: ierr
    ierr = rdysetdomainymomentumsource_(rdy_%c_rdy, size, c_loc(ymomsrc))
  end subroutine

  subroutine RDySetDomainManningsN(rdy_, size, n, ierr)
    type(RDy),       intent(inout)       :: rdy_
    integer(kind=selected_int_kind(5)),        intent(in)          :: size
    real(RDyDouble), pointer, intent(in) :: n(:)
    integer,         intent(out)         :: ierr
    ierr = rdysetdomainmanningsn_(rdy_%c_rdy, size, c_loc(n))
  end subroutine

  subroutine RDySetInitialConditions(rdy_, ic, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    type(tVec), intent(in)    :: ic  ! type(tVec)
    integer,    intent(out)   :: ierr
    ierr = rdysetinitialconditions_(rdy_%c_rdy, ic%v)
  end subroutine

  subroutine RDyCreatePrognosticVec(rdy_, prog_vec, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    type(tVec), intent(inout) :: prog_vec  ! type(tVec)
    integer,    intent(out)   :: ierr
    ierr = rdycreateprognosticvec_(rdy_%c_rdy, prog_vec%v)
  end subroutine

  subroutine RDyCreateOneDOFGlobalVec(rdy_, global_vec, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    type(tVec), intent(inout) :: global_vec  ! type(tVec)
    integer,    intent(out)   :: ierr
    ierr = rdycreateonedofglobalvec_(rdy_%c_rdy, global_vec%v)
  end subroutine

  subroutine RDyReadOneDOFLocalVecFromBinaryFile(rdy_, filename, local_vec, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    character(len=1024), intent(in) :: filename
    type(tVec), intent(inout) :: local_vec  ! type(tVec)
    integer,    intent(out)   :: ierr

    integer                      :: n
    character(len=1024), pointer :: binary_file

    n = len_trim(filename)
    allocate(binary_file)
    binary_file(1:n) = filename(1:n)
    binary_file(n+1:n+1) = c_null_char
    ierr = rdyreadonedoflocalvecfrombinaryfile_(rdy_%c_rdy, c_loc(binary_file), local_vec%v)
    deallocate(binary_file)
  end subroutine

  subroutine RDyReadOneDOFGlobalVecFromBinaryFile(rdy_, filename, global_vec, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    character(len=1024), intent(in) :: filename
    type(tVec), intent(inout) :: global_vec  ! type(tVec)
    integer,    intent(out)   :: ierr

    integer                      :: n
    character(len=1024), pointer :: binary_file

    n = len_trim(filename)
    allocate(binary_file)
    binary_file(1:n) = filename(1:n)
    binary_file(n+1:n+1) = c_null_char
    ierr = rdyreadonedofglobalvecfrombinaryfile_(rdy_%c_rdy, c_loc(binary_file), global_vec%v)
    deallocate(binary_file)
  end subroutine

  subroutine RDyWriteOneDOFGlobalVecToBinaryFile(rdy_, filename, global_vec, ierr)
    use petscvec
    type(RDy),  intent(inout) :: rdy_
    character(len=1024), intent(in) :: filename
    type(tVec), intent(inout) :: global_vec  ! type(tVec)
    integer,    intent(out)   :: ierr

    integer                      :: n
    character(len=1024), pointer :: binary_file

    n = len_trim(filename)
    allocate(binary_file)
    binary_file(1:n) = filename(1:n)
    binary_file(n+1:n+1) = c_null_char
    ierr = rdywriteonedofglobalvectobinaryfile_(rdy_%c_rdy, c_loc(binary_file), global_vec%v)
    deallocate(binary_file)
  end subroutine

  subroutine RDyAdvance(rdy_, ierr)
    type(RDy), intent(inout) :: rdy_
    integer,   intent(out)   :: ierr
    ierr = rdyadvance_(rdy_%c_rdy)
  end subroutine

  function RDyFinished(rdy_)
    type(RDy), intent(inout) :: rdy_
    logical :: RDyFinished
    RDyFinished = rdyfinished_(rdy_%c_rdy)
  end function

  function RDyRestarted(rdy_)
    type(RDy), intent(inout) :: rdy_
    logical :: RDyRestarted
    RDyRestarted = rdyrestarted_(rdy_%c_rdy)
  end function

  subroutine RDyDestroy(rdy_, ierr)
    type(RDy), intent(inout) :: rdy_
    integer,   intent(out)   :: ierr
    ierr = rdydestroy_(rdy_%c_rdy)
  end subroutine

end module rdycore
