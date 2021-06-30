function create_a_new_instance() result(pp)
  use plumed_module_f08
  IMPLICIT NONE
  type(plumed) :: pp
  call pp%create()
end function create_a_new_instance

SUBROUTINE TEST3A()
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED) :: pl1
  TYPE(PLUMED) :: pl2
  call pl1%cmd("init")
  pl2=pl1
  ! this is ok, both references are destroyed and the object is deleted
END SUBROUTINE TEST3A

SUBROUTINE TEST3B()
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED), allocatable :: all(:)
  allocate(all(1))
  call all(1)%cmd("init")
  ! deallocate not needed in fortran
  ! this is ok, the reference is destroyed and the object is deleted
END SUBROUTINE TEST3B

SUBROUTINE TEST3C()
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED), pointer :: all(:)
  allocate(all(1))
  call all(1)%cmd("init")
  deallocate(all)
  ! deallocate is needed here
  ! this is ok, the reference is destroyed and the object is deleted
END SUBROUTINE TEST3C

SUBROUTINE TEST3D()
  USE PLUMED_MODULE_F08
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
  USE PLUMED_MODULE_F08
  IMPLICIT NONE
  TYPE(PLUMED), ALLOCATABLE :: p(:)
  INTEGER :: i
  allocate(p(10))
  call p%cmd("init")
  call p%finalize()
end SUBROUTINE TEST3E

SUBROUTINE TEST3()
  IMPLICIT NONE
  CALL TEST3A()
  CALL TEST3B()
  CALL TEST3C()
  CALL TEST3D()
  CALL TEST3E()
END SUBROUTINE
