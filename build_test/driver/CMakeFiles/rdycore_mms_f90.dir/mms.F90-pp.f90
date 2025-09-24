# 1 "/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/wsugiarta/Documents/github/RDycore/build_test//"
# 1 "/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90"
module mms_driver

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




# 3 "/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90" 2
  use petsc
  implicit none

contains
  subroutine usage()
    print *, "rdycore_mms_f90: usage:"
    print *, "rdycore_mms_f90 <input.yaml>"
    print *, ""
  end subroutine

end module mms_driver

program mms_f90



# 1 "/Users/wsugiarta/Documents/github/RDycore/include/finclude/rdycore.h" 1




# 1 "/Users/wsugiarta/Documents/github/RDycore/build_test/include/private/config.h" 1




















# 5 "/Users/wsugiarta/Documents/github/RDycore/include/finclude/rdycore.h" 2







# 19 "/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90" 2

  use rdycore
  use mms_driver
  use petsc

  implicit none

  character(len=1024)  :: config_file
  type(RDy)            :: rdy_
  integer(kind=selected_int_kind(5))       :: ierr

  if (command_argument_count() < 1) then
    call usage()
  else
    call get_command_argument(1, config_file)

    ! initialize subsystems
    call RDyInit(ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,36,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

    if (trim(config_file) /= trim('-help')) then
      ! create rdycore and set it up with the given file
      call RDyCreate(PETSC_COMM_WORLD, config_file, rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,40,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
      call RDyMMSSetup(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,41,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif

      ! run the problem according to the given configuration
      call RDyMMSRun(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,44,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif;
      call RDyDestroy(rdy_, ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,45,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
    endif
    ! shut off
    call RDyFinalize(ierr); if (ierr .ne. 0) then;call PetscErrorF(ierr,48,"/Users/wsugiarta/Documents/github/RDycore/driver/mms.F90");call MPIU_Abort(PETSC_COMM_SELF,ierr);endif
  endif

end program
