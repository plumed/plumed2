! vim:ft=fortran




module plumed_module
  use iso_c_binding
  implicit none

  ! names are private by default
  private

  ! only these names are public
  public :: plumed_f_create
  public :: plumed_f_create_dlopen
  public :: plumed_f_create_reference
  public :: plumed_f_create_invalid
  public :: plumed_f_cmd
  public :: plumed_f_finalize
  public :: plumed_f_installed
  public :: plumed_f_valid
  public :: plumed_f_use_count
  public :: plumed_f_global
  public :: plumed_f_ginitialized
  public :: plumed_f_gcreate
  public :: plumed_f_gcmd
  public :: plumed_f_gfinalize
  public :: plumed_f_gvalid

  ! this type maps to the struct plumed defined in src/wrapper/Plumed.h
  type, bind(C) :: cplumed
    type(c_ptr) :: ptr
  end type cplumed

  ! this type maps to the struct plumed_safeptr defined in src/wrapper/Plumed.h
  type, bind(C) :: cplumed_safeptr
    type(c_ptr)            :: ptr
    integer(kind=c_size_t) :: nelem
    type(c_ptr)            :: shape
    integer(kind=c_size_t) :: flags
    type(c_ptr)            :: opt
  end type cplumed_safeptr

  integer(kind=c_size_t), parameter :: flags_ptr = 67108864 ! 0x2000000*2

  ! this function is used to translate 32-char to c identifiers, only used internally
  interface
    function plumed_f2c(c) bind(C)
      import
      character(kind=c_char), intent(in) :: c(32)
      type(cplumed) :: plumed_f2c
    end function plumed_f2c
  end interface

  ! this subroutine provides a typesafe interface to plumed
  interface
    subroutine plumed_cmd_safe(p,key,safeptr) bind(C)
      import
      type(cplumed),                  value :: p
      character(kind=c_char), intent(in)    :: key(*)
      type(cplumed_safeptr),          value :: safeptr
    end subroutine plumed_cmd_safe
  end interface

  ! now there are interfaces to the classic Fortran functions

  interface
    subroutine plumed_f_create(c) bind(C)
      import
      character(kind=c_char), intent(out) :: c(32)
    end subroutine plumed_f_create
  end interface

  interface
    subroutine plumed_f_create_dlopen(path,c) bind(C)
      import
      character(kind=c_char), intent(in)  :: path(*)
      character(kind=c_char), intent(out) :: c(32)
    end subroutine plumed_f_create_dlopen
  end interface

  interface
    subroutine plumed_f_create_reference(r,c) bind(C)
      import
      character(kind=c_char), intent(in)  :: r(32)
      character(kind=c_char), intent(out) :: c(32)
    end subroutine plumed_f_create_reference
  end interface

  interface
    subroutine plumed_f_create_invalid(c) bind(C)
      import
      character(kind=c_char), intent(out) :: c(32)
    end subroutine plumed_f_create_invalid
  end interface

  interface
    subroutine plumed_f_finalize(c) bind(C)
      import
      character(kind=c_char), intent(in) :: c(32)
    end subroutine plumed_f_finalize
  end interface

  interface
    subroutine plumed_f_installed(i) bind(C)
      import
      integer(kind=c_int), intent(out) :: i
    end subroutine plumed_f_installed
  end interface

  interface
    subroutine plumed_f_valid(c,i) bind(C)
      import
      character(kind=c_char), intent(in) :: c(32)
      integer(kind=c_int),    intent(out) :: i
    end subroutine plumed_f_valid
  end interface

  interface
    subroutine plumed_f_use_count(c,i) bind(C)
      import
      character(kind=c_char), intent(in)  :: c(32)
      integer(kind=c_int),    intent(out) :: i
    end subroutine plumed_f_use_count
  end interface

  interface
    subroutine plumed_f_global(c) bind(C)
      import
      character(kind=c_char), intent(out) :: c(32)
    end subroutine plumed_f_global
  end interface 

  interface
    subroutine plumed_f_ginitialized(i) bind(C)
      import
      integer(kind=c_int), intent(out) :: i
    end subroutine plumed_f_ginitialized
  end interface

  interface
    subroutine plumed_f_gcreate() bind(C)
    end subroutine plumed_f_gcreate
  end interface
  
  interface
    subroutine plumed_f_gfinalize() bind(C)
    end subroutine plumed_f_gfinalize
  end interface

  interface
    subroutine plumed_f_gvalid(i) bind(C)
      import
      integer(kind=c_int),    intent(out) :: i
    end subroutine plumed_f_gvalid
  end interface

  ! here the interfaces to C functions to construct the plumed_safeptr object
  interface
    function plumed_f_safeptr_ptr(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_ptr
    end function plumed_f_safeptr_ptr
  end interface

  interface
    function plumed_f_safeptr_char(val,nelem,pass_shape,flags,opt) bind(C)
      import
      character(kind=c_char)                :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_char
    end function plumed_f_safeptr_char
  end interface

  interface
    function plumed_f_safeptr_int_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_int)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_int_scalar
    end function plumed_f_safeptr_int_scalar
  end interface
  interface
    function plumed_f_safeptr_int(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_int)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_int
    end function plumed_f_safeptr_int
  end interface
  interface
    function plumed_f_safeptr_short_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_short)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_short_scalar
    end function plumed_f_safeptr_short_scalar
  end interface
  interface
    function plumed_f_safeptr_short(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_short)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_short
    end function plumed_f_safeptr_short
  end interface
  interface
    function plumed_f_safeptr_long_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_long)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_long_scalar
    end function plumed_f_safeptr_long_scalar
  end interface
  interface
    function plumed_f_safeptr_long(val,nelem,pass_shape,flags,opt) bind(C)
      import
      integer(kind=c_long)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_long
    end function plumed_f_safeptr_long
  end interface
  interface
    function plumed_f_safeptr_float_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_float)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_float_scalar
    end function plumed_f_safeptr_float_scalar
  end interface
  interface
    function plumed_f_safeptr_float(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_float)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_float
    end function plumed_f_safeptr_float
  end interface
  interface
    function plumed_f_safeptr_double_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_double)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_double_scalar
    end function plumed_f_safeptr_double_scalar
  end interface
  interface
    function plumed_f_safeptr_double(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_double)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_double
    end function plumed_f_safeptr_double
  end interface

  ! here are the interfaces used for overloading
  interface plumed_f_cmd
    module procedure plumed_f_cmd_char
    module procedure plumed_f_cmd_integer_0_0
    module procedure plumed_f_cmd_integer_0_1
    module procedure plumed_f_cmd_integer_0_2
    module procedure plumed_f_cmd_integer_0_3
    module procedure plumed_f_cmd_integer_0_4
    module procedure plumed_f_cmd_integer_1_0
    module procedure plumed_f_cmd_integer_1_1
    module procedure plumed_f_cmd_integer_1_2
    module procedure plumed_f_cmd_integer_1_3
    module procedure plumed_f_cmd_integer_1_4
    module procedure plumed_f_cmd_integer_2_0
    module procedure plumed_f_cmd_integer_2_1
    module procedure plumed_f_cmd_integer_2_2
    module procedure plumed_f_cmd_integer_2_3
    module procedure plumed_f_cmd_integer_2_4
    module procedure plumed_f_cmd_real_0_0
    module procedure plumed_f_cmd_real_0_1
    module procedure plumed_f_cmd_real_0_2
    module procedure plumed_f_cmd_real_0_3
    module procedure plumed_f_cmd_real_0_4
    module procedure plumed_f_cmd_real_1_0
    module procedure plumed_f_cmd_real_1_1
    module procedure plumed_f_cmd_real_1_2
    module procedure plumed_f_cmd_real_1_3
    module procedure plumed_f_cmd_real_1_4
  end interface plumed_f_cmd

  interface plumed_f_gcmd
    module procedure plumed_f_gcmd_char
    module procedure plumed_f_gcmd_integer_0_0
    module procedure plumed_f_gcmd_integer_0_1
    module procedure plumed_f_gcmd_integer_0_2
    module procedure plumed_f_gcmd_integer_0_3
    module procedure plumed_f_gcmd_integer_0_4
    module procedure plumed_f_gcmd_integer_1_0
    module procedure plumed_f_gcmd_integer_1_1
    module procedure plumed_f_gcmd_integer_1_2
    module procedure plumed_f_gcmd_integer_1_3
    module procedure plumed_f_gcmd_integer_1_4
    module procedure plumed_f_gcmd_integer_2_0
    module procedure plumed_f_gcmd_integer_2_1
    module procedure plumed_f_gcmd_integer_2_2
    module procedure plumed_f_gcmd_integer_2_3
    module procedure plumed_f_gcmd_integer_2_4
    module procedure plumed_f_gcmd_real_0_0
    module procedure plumed_f_gcmd_real_0_1
    module procedure plumed_f_gcmd_real_0_2
    module procedure plumed_f_gcmd_real_0_3
    module procedure plumed_f_gcmd_real_0_4
    module procedure plumed_f_gcmd_real_1_0
    module procedure plumed_f_gcmd_real_1_1
    module procedure plumed_f_gcmd_real_1_2
    module procedure plumed_f_gcmd_real_1_3
    module procedure plumed_f_gcmd_real_1_4
  end interface plumed_f_gcmd

  contains

     ! we then define all the functions needed for overloading

     subroutine plumed_f_cmd_char(p,key,val)
       character(kind=c_char,len=32), intent(in)    :: p
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*)                 :: val
       integer(kind=c_size_t) :: pass_shape(2)
       integer(kind=c_size_t) :: nelem
       pass_shape=(/len(val),0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_char(val,nelem,pass_shape,flags_ptr,c_null_ptr))
     end subroutine plumed_f_cmd_char

     subroutine plumed_f_gcmd_char(key,val)
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*)                 :: val
       character(kind=c_char,len=32) :: global
       call plumed_f_global(global)
       call plumed_f_cmd(global,key,val)
     end subroutine plumed_f_gcmd_char

    subroutine plumed_f_cmd_integer_0_0(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_int_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_0_0
    subroutine plumed_f_gcmd_integer_0_0(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_int)                            :: val
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_0_0
    subroutine plumed_f_cmd_integer_0_1(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_0_1
    subroutine plumed_f_gcmd_integer_0_1(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_int)                            :: val(:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_0_1
    subroutine plumed_f_cmd_integer_0_2(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_0_2
    subroutine plumed_f_gcmd_integer_0_2(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_int)                            :: val(:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_0_2
    subroutine plumed_f_cmd_integer_0_3(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_0_3
    subroutine plumed_f_gcmd_integer_0_3(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_int)                            :: val(:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_0_3
    subroutine plumed_f_cmd_integer_0_4(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_0_4
    subroutine plumed_f_gcmd_integer_0_4(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_int)                            :: val(:,:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_0_4
    subroutine plumed_f_cmd_integer_1_0(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_short_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_1_0
    subroutine plumed_f_gcmd_integer_1_0(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_short)                            :: val
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_1_0
    subroutine plumed_f_cmd_integer_1_1(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_1_1
    subroutine plumed_f_gcmd_integer_1_1(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_short)                            :: val(:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_1_1
    subroutine plumed_f_cmd_integer_1_2(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_1_2
    subroutine plumed_f_gcmd_integer_1_2(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_short)                            :: val(:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_1_2
    subroutine plumed_f_cmd_integer_1_3(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_1_3
    subroutine plumed_f_gcmd_integer_1_3(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_short)                            :: val(:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_1_3
    subroutine plumed_f_cmd_integer_1_4(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_1_4
    subroutine plumed_f_gcmd_integer_1_4(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_short)                            :: val(:,:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_1_4
    subroutine plumed_f_cmd_integer_2_0(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_long_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_2_0
    subroutine plumed_f_gcmd_integer_2_0(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_long)                            :: val
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_2_0
    subroutine plumed_f_cmd_integer_2_1(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_2_1
    subroutine plumed_f_gcmd_integer_2_1(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_long)                            :: val(:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_2_1
    subroutine plumed_f_cmd_integer_2_2(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_2_2
    subroutine plumed_f_gcmd_integer_2_2(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_long)                            :: val(:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_2_2
    subroutine plumed_f_cmd_integer_2_3(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_2_3
    subroutine plumed_f_gcmd_integer_2_3(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_long)                            :: val(:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_2_3
    subroutine plumed_f_cmd_integer_2_4(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_integer_2_4
    subroutine plumed_f_gcmd_integer_2_4(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(kind=c_long)                            :: val(:,:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_integer_2_4
    subroutine plumed_f_cmd_real_0_0(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_float_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_0_0
    subroutine plumed_f_gcmd_real_0_0(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_float)                            :: val
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_0_0
    subroutine plumed_f_cmd_real_0_1(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_0_1
    subroutine plumed_f_gcmd_real_0_1(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_float)                            :: val(:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_0_1
    subroutine plumed_f_cmd_real_0_2(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_0_2
    subroutine plumed_f_gcmd_real_0_2(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_float)                            :: val(:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_0_2
    subroutine plumed_f_cmd_real_0_3(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_0_3
    subroutine plumed_f_gcmd_real_0_3(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_float)                            :: val(:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_0_3
    subroutine plumed_f_cmd_real_0_4(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_0_4
    subroutine plumed_f_gcmd_real_0_4(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_float)                            :: val(:,:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_0_4
    subroutine plumed_f_cmd_real_1_0(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_double_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_1_0
    subroutine plumed_f_gcmd_real_1_0(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_double)                            :: val
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_1_0
    subroutine plumed_f_cmd_real_1_1(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_1_1
    subroutine plumed_f_gcmd_real_1_1(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_double)                            :: val(:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_1_1
    subroutine plumed_f_cmd_real_1_2(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_1_2
    subroutine plumed_f_gcmd_real_1_2(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_double)                            :: val(:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_1_2
    subroutine plumed_f_cmd_real_1_3(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_1_3
    subroutine plumed_f_gcmd_real_1_3(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_double)                            :: val(:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_1_3
    subroutine plumed_f_cmd_real_1_4(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:,:)
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       call plumed_cmd_safe(plumed_f2c(p),key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr))
    end subroutine plumed_f_cmd_real_1_4
    subroutine plumed_f_gcmd_real_1_4(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      real(kind=c_double)                            :: val(:,:,:,:)
      character(kind=c_char,len=32) :: global
      call plumed_f_global(global)
      call plumed_f_cmd(global,key,val)
    end subroutine plumed_f_gcmd_real_1_4

end module plumed_module

