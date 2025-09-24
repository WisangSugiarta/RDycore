# 1 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/wsugiarta/Documents/github/RDycore/build_test//"
# 1 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90"
! This module was adapted from https://github.com/RDycore/E3SM/blob/bishtgautam/mosart-rdycore/78a5ea1613-2023-08-07/components/mosart/src/rdycore/rdycoreMod.F90

module rdycoreMod


# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbag.h" 1




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

# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbag.h" 2






# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmatlab.h" 1








# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbm.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/include/petsc/finclude/petscisbase.h" 1
!
!  Part of the base include file for Fortran use of PETSc IS
!  Note: This file should contain only define statements

! No spaces for #defines as some compilers (PGI) also adds
! those additional spaces during preprocessing - bad for fixed format
!





# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 2


# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h" 1





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
# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 1
# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscvec.h" 2



















# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h" 2

# 18 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsf.h"






# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsection.h" 1













# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscis.h" 2















# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscbm.h" 2






# 11 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmda.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmat.h" 1





# 38 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscmat.h"

















# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmlabel.h" 1










# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfe.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h" 1
# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfe.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdt.h" 1

















# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfe.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscds.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfe.h" 1
# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscds.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfv.h" 1






# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscspace.h" 1











# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfv.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdualspace.h" 1







# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscds.h" 1
# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdualspace.h" 2









# 8 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfv.h" 2











# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscds.h" 2










# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscfe.h" 2











# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h" 2



# 23 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdm.h"









# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmda.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscpf.h" 1










# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmda.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscao.h" 1










# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmda.h" 2










# 14 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmcomposite.h" 1








# 15 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmpatch.h" 1








# 16 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmplex.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscpartitioner.h" 1











# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmplex.h" 2







# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmfield.h" 1












# 13 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmplex.h" 2














# 17 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmplextransform.h" 1










# 18 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmredundant.h" 1








# 19 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmshell.h" 1








# 20 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmsliced.h" 1








# 21 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmswarm.h" 1


















# 22 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmstag.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmproduct.h" 1








# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmstag.h" 2







# 23 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmforest.h" 1










# 24 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmnetwork.h" 1












# 25 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmadaptor.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsnes.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscksp.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscpc.h" 1






# 29 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscpc.h"







# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscksp.h" 2


# 19 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscksp.h"







# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsnes.h" 2




# 24 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscsnes.h"








# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmadaptor.h" 2






# 26 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscdmlabelephemeral.h" 1









# 27 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclandau.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscts.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscconvest.h" 1









# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscts.h" 2

# 15 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscts.h"

# 31 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscts.h"

# 45 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscts.h"

# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsclandau.h" 2










# 28 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2



# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsccharacteristic.h" 1











# 31 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2


# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctao.h" 1





# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctaolinesearch.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctao.h" 1
# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctaolinesearch.h" 2







# 6 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctao.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctao_deprecations.h" 1







# 7 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsctao.h" 2














# 33 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2

# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscml.h" 1




# 1 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscregressor.h" 1











# 5 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petscml.h" 2




# 34 "/Users/wsugiarta/Documents/github/petsc/arch-darwin-c-debug/include/petsc/finclude/petsc.h" 2




# 6 "/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90" 2

  use petsc
  use rdycore

  implicit none

  private

  type(RDy)          :: rdy_


  integer(kind=selected_int_kind(5))           :: num_cells_owned
  real(kind=selected_real_kind(10)), pointer :: rain_data(:)
  real(kind=selected_real_kind(10))          :: dtime ! in units expressed within config file
  integer(kind=selected_int_kind(5))           :: nstep

  public :: rdycore_usage
  public :: rdycore_init
  public :: rdycore_run
  public :: rdycore_final

  public :: dtime, nstep

contains

  !-----------------------------------------------------------------------
  subroutine rdycore_usage()
    print *, "test_coupling: usage:"
    print *, "test_coupling <input.yaml>"
    print *, ""
  end subroutine

  !-----------------------------------------------------------------------
  subroutine rdycore_init(config_file)
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !
    implicit none
    !
    character(len=1024), intent(in) :: config_file
    !
    !
    ! !LOCAL VARIABLES:
    character(len=1024) :: log_file = 'rof_modelio.log'
    type(tPetscViewer)         :: viewer
    integer(kind=selected_int_kind(5))            :: num_cells_global
    integer(kind=selected_int_kind(5))      :: ierr

    ! initialize subsystems
    call RDyInit(ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,58,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

    ! create rdycore and set it up with the given file
    call RDyCreate(PETSC_COMM_WORLD, config_file, rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,61,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
    call RDySetLogFile(rdy_, log_file, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,62,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
    call RDySetup(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,63,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

    ! allocate memory for grid-level rain data
    call RDyGetNumLocalCells(rdy_, num_cells_owned, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,66,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
    allocate(rain_data(num_cells_owned))

    call RDyGetNumLocalCells(rdy_, num_cells_global, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,69,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

  end subroutine rdycore_init

  !-----------------------------------------------------------------------
  subroutine rdycore_run()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    real(kind=selected_real_kind(10)), pointer :: rain_p(:)
    integer(kind=selected_int_kind(5))             :: t
    integer(RDyTimeUnit) :: time_unit
    real(kind=selected_real_kind(10))            :: time_dn, time_up, cur_time, cur_rain
    integer(kind=selected_int_kind(5))       :: ierr

    cur_time = (nstep-1)*dtime

    ! Set spatially homogeneous rainfall for all grid cells
    ! (the domain is region 0)
    rain_data(:) = 0.d0
    call RDySetDomainWaterSource(rdy_, num_cells_owned, rain_data, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,93,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

    ! Set the coupling time step
    call RDyGetTimeUnit(rdy_, time_unit, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,96,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
    call RDySetCouplingInterval(rdy_, time_unit, dtime, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,97,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

    ! Run the simulation to completion.
    call RDyAdvance(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,100,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

  end subroutine rdycore_run

  !-----------------------------------------------------------------------
  subroutine rdycore_final()
    !
    ! !DESCRIPTION:
    ! Destroy RDy object
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer(kind=selected_int_kind(5)) :: ierr

    ! deallocate memory for rain data
    deallocate(rain_data)

    ! destroy RDy object
    call RDyDestroy(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,121,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif;

    ! finalize
    call RDyFinalize(ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,124,"/Users/wsugiarta/Documents/github/RDycore/src/f90-mod/tests/test_coupling.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif;

  end subroutine rdycore_final

end module rdycoreMod

program test_coupling

  use rdycore
  use rdycoreMod
  use petsc

  implicit none

  character(len=1024) :: config_file

  if (command_argument_count() < 1) then
    call rdycore_usage()
  else
    ! fetch config file name
    call get_command_argument(1, config_file)
  end if

  nstep = 0
  dtime = 1 ! hours for ex2b.yaml

  call rdycore_init(config_file)
  call rdycore_run()
  call rdycore_final()

end program
