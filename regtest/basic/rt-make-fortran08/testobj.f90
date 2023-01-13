function create_a_new_instance() result(pp)
  use plumed_f08_module
  IMPLICIT NONE
  type(plumed) :: pp
  call plumed_create(pp)
end function create_a_new_instance

SUBROUTINE TEST3A()
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED) :: pl1
  TYPE(PLUMED) :: pl2
  call pl1%cmd("init")
  pl2=pl1
  ! this is ok, both references are destroyed and the object is deleted
END SUBROUTINE TEST3A

SUBROUTINE TEST3B()
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED), allocatable :: all(:)
  allocate(all(1))
  call all(1)%cmd("init")
  ! deallocate not needed in fortran
  ! this is ok, the reference is destroyed and the object is deleted
END SUBROUTINE TEST3B

SUBROUTINE TEST3C()
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED), pointer :: all(:)
  allocate(all(1))
  call all(1)%cmd("init")
  deallocate(all)
  ! deallocate is needed here
  ! this is ok, the reference is destroyed and the object is deleted
END SUBROUTINE TEST3C

SUBROUTINE TEST3D()
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED) :: create_a_new_instance
  TYPE(PLUMED) :: p
  ! here, the temporary object returned by create_a_new_instance() is not destroyed correctly
  ! therefore, there is a leak
  p=create_a_new_instance()
  call p%cmd("init")
  ! this is the recommended workaround: it will decrease the number of references until there is just one
  ! it should work with all compilers (those with the bug, and those without)
  call p%decref(1)
END SUBROUTINE TEST3D

SUBROUTINE TEST3E()
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE(PLUMED), ALLOCATABLE :: p(:)
  INTEGER :: i
  allocate(p(10))
  call p%cmd("init")
  call p%finalize()
end SUBROUTINE TEST3E

SUBROUTINE TEST3()
  IMPLICIT NONE
  open(10,file="log")
  CALL TEST3A()
  write(10,*) "3A"
  CALL TEST3B()
  write(10,*) "3B"
  CALL TEST3C()
  write(10,*) "3C"
  CALL TEST3D()
  write(10,*) "3D"
  CALL TEST3E()
  write(10,*) "3E"
  close(10)
END SUBROUTINE

MODULE TEST_DERIVED
  USE PLUMED_F08_MODULE
  IMPLICIT NONE
  TYPE, EXTENDS(PLUMED) :: PLUMEDX
    INTEGER :: i=77
  END TYPE PLUMEDX
END MODULE TEST_DERIVED

SUBROUTINE TEST4()
  USE TEST_DERIVED
  TYPE(PLUMEDX) :: px
  open(10,file="log",position='append')
  call px%cmd("init")
  write(10,*) "4",px%i
  close(10)
END SUBROUTINE TEST4

SUBROUTINE TEST5()
  USE PLUMED_F08_MODULE
  TYPE(PLUMED), ALLOCATABLE :: b(:)
  TYPE(PLUMED), ALLOCATABLE :: c(:)
  ALLOCATE(b(3))
  ALLOCATE(c(3))
  open(10,file="log",position='append')
  write(10,*)"B",b%use_count()
  write(10,*)"C",c%use_count()
  c=b
  write(10,*)"B",b%use_count()
  write(10,*)"C",c%use_count()
  deallocate(b)
  write(10,*)"C",c%use_count()
  close(10)
END SUBROUTINE TEST5

SUBROUTINE TEST6()
  USE PLUMED_F08_MODULE
  USE ISO_C_BINDING
  IMPLICIT NONE
  TYPE(PLUMED) :: pippo
  TYPE(PLUMED_ERROR), target  :: error
  TYPE(PLUMED_ERROR), target  :: error2
  TYPE(PLUMED_ERROR), pointer :: error_nested
  open(10,file="error_codes")
  error%code=-1 ! check if this is overwritten
  call pippo%cmd("init",error=error)
  write(10,*) "should be zero",error%code
  error%code=-1 ! check if this is overwritten
  call pippo%cmd("initxx",error=error)
  write(10,*) "should be nonzero",error%code
 
  call plumed_create(pippo) ! reset instance
  error%code=-1 ! check if this is overwritten
  call pippo%cmd("init",error=error)
  write(10,*) "should be zero",error%code
  error%code=-1 ! check if this is overwritten
  call pippo%cmd("initxx",error=error)
  write(10,*) "should be nonzero",error%code
 
  ! test nested exceptions
 
  ! first, disabled
   call pippo%cmd_val("throw","test_nested2",error=error)
  write(10,*) error%code,allocated(error%nested)
  write(10,'(A)') error%what
 
  ! then, enabled
  call pippo%cmd_val("setNestedExceptions",1)
  call pippo%cmd_val("throw","test_nested2",error=error)
  write(10,*) "code:",error%code,"nested:",allocated(error%nested)
  write(10,'(A)') error%what
 
 
  ! this is a deep copy, not needed in general
  ! however, the default copy is buggy in gfortran
  ! and here we test our implementation
  error2=error
  error=error2
 
  ! here we write the linked list in forward order
  error_nested => error
  do while(allocated(error_nested%nested))
    error_nested => error_nested%nested
    write(10,*) "code:",error_nested%code,"nested:",allocated(error_nested%nested)
    write(10,'(A)') error_nested%what
  end do

  call plumed_create(pippo) ! reset instance
  call pippo%cmd_val("setNatoms",999)
  call pippo%cmd("init")
  close(10)
END SUBROUTINE TEST6

