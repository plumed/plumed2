! vim:ft=fortran




module plumed_f08_module
  use iso_c_binding
  implicit none

  ! names are private by default
  private

  ! only these names are public
  public :: plumed
  public :: plumed_create
  public :: plumed_installed
  public :: plumed_error

  ! this type maps to the struct plumed defined in src/wrapper/Plumed.h
  type, bind(C) :: cplumed
    type(c_ptr) :: ptr
  end type cplumed

  ! this type maps to the struct plumed_nothrow_handler defined in src/wrapper/Plumed.h
  type, bind(C) :: cplumed_nothrow_handler
    type(c_ptr)    :: ptr
    type(c_funptr) :: handler
  end type cplumed_nothrow_handler

  ! this type maps to the struct plumed_safeptr defined in src/wrapper/Plumed.h
  type, bind(C) :: cplumed_safeptr
    type(c_ptr)            :: ptr
    integer(kind=c_size_t) :: nelem
    type(c_ptr)            :: shape
    integer(kind=c_size_t) :: flags
    type(c_ptr)            :: opt
  end type cplumed_safeptr

  ! this subroutine provides a typesafe, not throwing interfafe to plumed
  interface
    subroutine plumed_cmd_safe_nothrow(p,key,safeptr,nothrow) bind(C)
      import
      type(cplumed),                  value :: p
      character(kind=c_char), intent(in)    :: key(*)
      type(cplumed_safeptr),          value :: safeptr
      type(cplumed_nothrow_handler),  value :: nothrow
    end subroutine plumed_cmd_safe_nothrow
  end interface

  integer(kind=c_size_t), parameter :: flags_ptr = 67108864 ! 0x2000000*2
  integer(kind=c_size_t), parameter :: flags_const_ptr = 100663296 ! 0x2000000*3
  integer(kind=c_size_t), parameter :: flags_nocopy = 268435456 ! 0x10000000

  ! this type is used to manipulate a plumed instance
  ! since it has a destructor, if contains a logical label to check if it has been initialized
  ! it also contains a number of methods
  ! notice that constructor (plumed_create) is NOT a member
  type plumed
    type(cplumed), private :: handle
    logical,       private :: initialized = .false.
  contains
    private

    generic, public :: cmd => pl_cmd
    procedure       ::        pl_cmd

    generic, public :: cmd_val => pl_cmd_val_integer_0_0
    procedure       ::             pl_cmd_val_integer_0_0
    generic, public :: cmd_val => pl_cmd_val_integer_0_1
    procedure       ::             pl_cmd_val_integer_0_1
    generic, public :: cmd_val => pl_cmd_val_integer_0_2
    procedure       ::             pl_cmd_val_integer_0_2
    generic, public :: cmd_val => pl_cmd_val_integer_0_3
    procedure       ::             pl_cmd_val_integer_0_3
    generic, public :: cmd_val => pl_cmd_val_integer_0_4
    procedure       ::             pl_cmd_val_integer_0_4
    generic, public :: cmd_val => pl_cmd_val_integer_1_0
    procedure       ::             pl_cmd_val_integer_1_0
    generic, public :: cmd_val => pl_cmd_val_integer_1_1
    procedure       ::             pl_cmd_val_integer_1_1
    generic, public :: cmd_val => pl_cmd_val_integer_1_2
    procedure       ::             pl_cmd_val_integer_1_2
    generic, public :: cmd_val => pl_cmd_val_integer_1_3
    procedure       ::             pl_cmd_val_integer_1_3
    generic, public :: cmd_val => pl_cmd_val_integer_1_4
    procedure       ::             pl_cmd_val_integer_1_4
    generic, public :: cmd_val => pl_cmd_val_integer_2_0
    procedure       ::             pl_cmd_val_integer_2_0
    generic, public :: cmd_val => pl_cmd_val_integer_2_1
    procedure       ::             pl_cmd_val_integer_2_1
    generic, public :: cmd_val => pl_cmd_val_integer_2_2
    procedure       ::             pl_cmd_val_integer_2_2
    generic, public :: cmd_val => pl_cmd_val_integer_2_3
    procedure       ::             pl_cmd_val_integer_2_3
    generic, public :: cmd_val => pl_cmd_val_integer_2_4
    procedure       ::             pl_cmd_val_integer_2_4
    generic, public :: cmd_val => pl_cmd_val_real_0_0
    procedure       ::             pl_cmd_val_real_0_0
    generic, public :: cmd_val => pl_cmd_val_real_0_1
    procedure       ::             pl_cmd_val_real_0_1
    generic, public :: cmd_val => pl_cmd_val_real_0_2
    procedure       ::             pl_cmd_val_real_0_2
    generic, public :: cmd_val => pl_cmd_val_real_0_3
    procedure       ::             pl_cmd_val_real_0_3
    generic, public :: cmd_val => pl_cmd_val_real_0_4
    procedure       ::             pl_cmd_val_real_0_4
    generic, public :: cmd_val => pl_cmd_val_real_1_0
    procedure       ::             pl_cmd_val_real_1_0
    generic, public :: cmd_val => pl_cmd_val_real_1_1
    procedure       ::             pl_cmd_val_real_1_1
    generic, public :: cmd_val => pl_cmd_val_real_1_2
    procedure       ::             pl_cmd_val_real_1_2
    generic, public :: cmd_val => pl_cmd_val_real_1_3
    procedure       ::             pl_cmd_val_real_1_3
    generic, public :: cmd_val => pl_cmd_val_real_1_4
    procedure       ::             pl_cmd_val_real_1_4
    generic, public :: cmd_ref => pl_cmd_ref_integer_0_0
    procedure       ::             pl_cmd_ref_integer_0_0
    generic, public :: cmd_ref => pl_cmd_ref_integer_0_1
    procedure       ::             pl_cmd_ref_integer_0_1
    generic, public :: cmd_ref => pl_cmd_ref_integer_0_2
    procedure       ::             pl_cmd_ref_integer_0_2
    generic, public :: cmd_ref => pl_cmd_ref_integer_0_3
    procedure       ::             pl_cmd_ref_integer_0_3
    generic, public :: cmd_ref => pl_cmd_ref_integer_0_4
    procedure       ::             pl_cmd_ref_integer_0_4
    generic, public :: cmd_ref => pl_cmd_ref_integer_1_0
    procedure       ::             pl_cmd_ref_integer_1_0
    generic, public :: cmd_ref => pl_cmd_ref_integer_1_1
    procedure       ::             pl_cmd_ref_integer_1_1
    generic, public :: cmd_ref => pl_cmd_ref_integer_1_2
    procedure       ::             pl_cmd_ref_integer_1_2
    generic, public :: cmd_ref => pl_cmd_ref_integer_1_3
    procedure       ::             pl_cmd_ref_integer_1_3
    generic, public :: cmd_ref => pl_cmd_ref_integer_1_4
    procedure       ::             pl_cmd_ref_integer_1_4
    generic, public :: cmd_ref => pl_cmd_ref_integer_2_0
    procedure       ::             pl_cmd_ref_integer_2_0
    generic, public :: cmd_ref => pl_cmd_ref_integer_2_1
    procedure       ::             pl_cmd_ref_integer_2_1
    generic, public :: cmd_ref => pl_cmd_ref_integer_2_2
    procedure       ::             pl_cmd_ref_integer_2_2
    generic, public :: cmd_ref => pl_cmd_ref_integer_2_3
    procedure       ::             pl_cmd_ref_integer_2_3
    generic, public :: cmd_ref => pl_cmd_ref_integer_2_4
    procedure       ::             pl_cmd_ref_integer_2_4
    generic, public :: cmd_ref => pl_cmd_ref_real_0_0
    procedure       ::             pl_cmd_ref_real_0_0
    generic, public :: cmd_ref => pl_cmd_ref_real_0_1
    procedure       ::             pl_cmd_ref_real_0_1
    generic, public :: cmd_ref => pl_cmd_ref_real_0_2
    procedure       ::             pl_cmd_ref_real_0_2
    generic, public :: cmd_ref => pl_cmd_ref_real_0_3
    procedure       ::             pl_cmd_ref_real_0_3
    generic, public :: cmd_ref => pl_cmd_ref_real_0_4
    procedure       ::             pl_cmd_ref_real_0_4
    generic, public :: cmd_ref => pl_cmd_ref_real_1_0
    procedure       ::             pl_cmd_ref_real_1_0
    generic, public :: cmd_ref => pl_cmd_ref_real_1_1
    procedure       ::             pl_cmd_ref_real_1_1
    generic, public :: cmd_ref => pl_cmd_ref_real_1_2
    procedure       ::             pl_cmd_ref_real_1_2
    generic, public :: cmd_ref => pl_cmd_ref_real_1_3
    procedure       ::             pl_cmd_ref_real_1_3
    generic, public :: cmd_ref => pl_cmd_ref_real_1_4
    procedure       ::             pl_cmd_ref_real_1_4
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_0_0
    procedure       ::             pl_cmd_ptr_integer_0_0
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_0_1
    procedure       ::             pl_cmd_ptr_integer_0_1
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_0_2
    procedure       ::             pl_cmd_ptr_integer_0_2
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_0_3
    procedure       ::             pl_cmd_ptr_integer_0_3
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_0_4
    procedure       ::             pl_cmd_ptr_integer_0_4
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_1_0
    procedure       ::             pl_cmd_ptr_integer_1_0
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_1_1
    procedure       ::             pl_cmd_ptr_integer_1_1
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_1_2
    procedure       ::             pl_cmd_ptr_integer_1_2
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_1_3
    procedure       ::             pl_cmd_ptr_integer_1_3
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_1_4
    procedure       ::             pl_cmd_ptr_integer_1_4
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_2_0
    procedure       ::             pl_cmd_ptr_integer_2_0
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_2_1
    procedure       ::             pl_cmd_ptr_integer_2_1
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_2_2
    procedure       ::             pl_cmd_ptr_integer_2_2
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_2_3
    procedure       ::             pl_cmd_ptr_integer_2_3
    generic, public :: cmd_ptr => pl_cmd_ptr_integer_2_4
    procedure       ::             pl_cmd_ptr_integer_2_4
    generic, public :: cmd_ptr => pl_cmd_ptr_real_0_0
    procedure       ::             pl_cmd_ptr_real_0_0
    generic, public :: cmd_ptr => pl_cmd_ptr_real_0_1
    procedure       ::             pl_cmd_ptr_real_0_1
    generic, public :: cmd_ptr => pl_cmd_ptr_real_0_2
    procedure       ::             pl_cmd_ptr_real_0_2
    generic, public :: cmd_ptr => pl_cmd_ptr_real_0_3
    procedure       ::             pl_cmd_ptr_real_0_3
    generic, public :: cmd_ptr => pl_cmd_ptr_real_0_4
    procedure       ::             pl_cmd_ptr_real_0_4
    generic, public :: cmd_ptr => pl_cmd_ptr_real_1_0
    procedure       ::             pl_cmd_ptr_real_1_0
    generic, public :: cmd_ptr => pl_cmd_ptr_real_1_1
    procedure       ::             pl_cmd_ptr_real_1_1
    generic, public :: cmd_ptr => pl_cmd_ptr_real_1_2
    procedure       ::             pl_cmd_ptr_real_1_2
    generic, public :: cmd_ptr => pl_cmd_ptr_real_1_3
    procedure       ::             pl_cmd_ptr_real_1_3
    generic, public :: cmd_ptr => pl_cmd_ptr_real_1_4
    procedure       ::             pl_cmd_ptr_real_1_4
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_0_0
    procedure       ::             pl_cmd_const_ptr_integer_0_0
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_0_1
    procedure       ::             pl_cmd_const_ptr_integer_0_1
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_0_2
    procedure       ::             pl_cmd_const_ptr_integer_0_2
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_0_3
    procedure       ::             pl_cmd_const_ptr_integer_0_3
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_0_4
    procedure       ::             pl_cmd_const_ptr_integer_0_4
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_1_0
    procedure       ::             pl_cmd_const_ptr_integer_1_0
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_1_1
    procedure       ::             pl_cmd_const_ptr_integer_1_1
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_1_2
    procedure       ::             pl_cmd_const_ptr_integer_1_2
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_1_3
    procedure       ::             pl_cmd_const_ptr_integer_1_3
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_1_4
    procedure       ::             pl_cmd_const_ptr_integer_1_4
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_2_0
    procedure       ::             pl_cmd_const_ptr_integer_2_0
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_2_1
    procedure       ::             pl_cmd_const_ptr_integer_2_1
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_2_2
    procedure       ::             pl_cmd_const_ptr_integer_2_2
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_2_3
    procedure       ::             pl_cmd_const_ptr_integer_2_3
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_integer_2_4
    procedure       ::             pl_cmd_const_ptr_integer_2_4
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_0_0
    procedure       ::             pl_cmd_const_ptr_real_0_0
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_0_1
    procedure       ::             pl_cmd_const_ptr_real_0_1
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_0_2
    procedure       ::             pl_cmd_const_ptr_real_0_2
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_0_3
    procedure       ::             pl_cmd_const_ptr_real_0_3
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_0_4
    procedure       ::             pl_cmd_const_ptr_real_0_4
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_1_0
    procedure       ::             pl_cmd_const_ptr_real_1_0
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_1_1
    procedure       ::             pl_cmd_const_ptr_real_1_1
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_1_2
    procedure       ::             pl_cmd_const_ptr_real_1_2
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_1_3
    procedure       ::             pl_cmd_const_ptr_real_1_3
    generic, public :: cmd_const_ptr => pl_cmd_const_ptr_real_1_4
    procedure       ::             pl_cmd_const_ptr_real_1_4
    generic, public :: cmd_val => pl_cmd_val_char
    procedure       ::            pl_cmd_val_char
    generic, public :: cmd_ptr => pl_cmd_ptr_c
    procedure       ::            pl_cmd_ptr_c
    generic, public :: cmd_const_ptr => pl_cmd_ptr_c
    procedure       ::                  pl_cmd_const_ptr_c

    procedure, public :: finalize => pl_finalize
    procedure, public :: incref => pl_incref
    procedure, public :: decref => pl_decref
    generic,   public :: assignment(=) => pl_assign
    final     :: pl_destructor
    procedure, public :: valid => pl_valid
    procedure, public :: use_count => pl_use_count
    procedure :: pl_assign
  end type plumed

  ! this type holds the information associated to a thrown exception
  type :: plumed_error
    integer                         :: code=0
    character(len = :), allocatable :: what
    ! nested error, if present
    type(plumed_error), allocatable :: nested
  contains
    private
    generic,   public :: assignment(=) => pl_error_assign
    procedure :: pl_error_assign
  end type plumed_error

  ! now there are interfaces to some of the classic C functions, only used internally

  interface
    function cplumed_create() bind(C,name="plumed_create")
      import
      type(cplumed) :: cplumed_create
    end function cplumed_create
  end interface

  interface
    function cplumed_create_dlopen(path) bind(C,name="plumed_create_dlopen")
      import
      character(kind=c_char), intent(in)  :: path(*)
      type(cplumed) :: cplumed_create_dlopen
    end function cplumed_create_dlopen
  end interface

  interface
    function cplumed_create_reference(p) bind(C,name="plumed_create_reference")
      import
      type(cplumed), value ::  p
      type(cplumed)        :: cplumed_create_reference
    end function cplumed_create_reference
  end interface

  interface
    subroutine cplumed_finalize(p) bind(C,name="plumed_finalize")
      import
      type(cplumed), value :: p
    end subroutine cplumed_finalize
  end interface

  interface
    function cplumed_installed() bind(C,name="plumed_installed")
      import
      integer(kind=c_int) :: cplumed_installed
    end function cplumed_installed
  end interface

  interface
    function cplumed_valid(p) bind(C,name="plumed_valid")
      import
      type(cplumed), value :: p
      integer(kind=c_int) :: cplumed_valid
    end function cplumed_valid
  end interface

  interface
    function cplumed_use_count(p) bind(C,name="plumed_use_count")
      import
      type(cplumed), value :: p
      integer(kind=c_int) :: cplumed_use_count
    end function cplumed_use_count
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
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_char
    end function plumed_f_safeptr_char
  end interface

  interface
    function plumed_f_safeptr_int(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_int
    end function plumed_f_safeptr_int
  end interface
  interface
    function plumed_f_safeptr_short(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_short
    end function plumed_f_safeptr_short
  end interface
  interface
    function plumed_f_safeptr_long(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_long
    end function plumed_f_safeptr_long
  end interface
  interface
    function plumed_f_safeptr_float(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_float
    end function plumed_f_safeptr_float
  end interface
  interface
    function plumed_f_safeptr_double(val,nelem,pass_shape,flags,opt) bind(C)
      import
      type(c_ptr),                    value :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_double
    end function plumed_f_safeptr_double
  end interface

  contains

     ! this is a callback function.
     ! notice that it ends up in global namespace (no protection for being a module function!)
     ! be careful with name thus
     subroutine plumed_f_f08_eh(error_ptr,code,what_ptr,opt_ptr) bind(C)
       type(c_ptr),         value :: error_ptr
       integer(kind=c_int), value :: code
       type(c_ptr),         value :: what_ptr
       type(c_ptr),         value :: opt_ptr
       type(plumed_error), pointer :: error
       character(len=1, kind=C_CHAR), pointer :: p_chars(:)
       type(c_ptr),        pointer :: opt(:)
       character(len=1, kind=C_CHAR), pointer :: opt_key
       type(c_ptr),        pointer :: error_nested
       integer :: i,j
       call c_f_pointer(error_ptr,error)
       error%code=code
       if (.not. C_associated(what_ptr)) then
         error%what=""
       else
         call C_F_pointer(what_ptr, p_chars, [huge(0)])
         do i = 1, huge(0)
           if (p_chars(i) == C_NULL_CHAR) exit
         enddo
         allocate(character(i-1) :: error%what)
         do j = 1,i-1
           error%what(j:j)=p_chars(j)
         enddo
       endif
       if (C_associated(opt_ptr)) then
         call C_F_pointer(opt_ptr,opt, [huge(0)])
         do i = 1, huge(0),2
           if (.not. c_associated(opt(i))) exit
           if (c_associated(opt(i+1))) then
             call C_F_pointer(opt(i),opt_key)
             if (opt_key == "n") then
               call C_F_pointer(opt(i+1),error_nested)
               allocate(error%nested)
               error_nested=c_loc(error%nested)
               exit ! make sure only the first "n" pointer is used
             endif
           endif
         enddo
       endif
     end subroutine plumed_f_f08_eh

     ! we then define all the functions needed for overloading

     subroutine plumed_f_cmd_ptr(p,key,val,const,nocopy,error)
       type(cplumed),                 intent(in)    :: p
       character(kind=c_char,len=*),  intent(in)    :: key
       type(c_ptr),                   intent(in)    :: val
       logical,                       intent(in)    :: const
       logical,                       intent(in)    :: nocopy
       type(plumed_error), optional,  target, intent(out)   :: error
       integer(kind=c_size_t) :: pass_shape(1)
       integer(kind=c_size_t) :: flags
       type(cplumed_nothrow_handler) :: nothrow
       integer(kind=c_size_t) :: nelem
       pass_shape=[0]
       nelem=0
       flags=flags_ptr
       if (const) flags=flags_const_ptr
       if (nocopy) flags=flags+flags_nocopy
       if(present(error)) then
         nothrow%ptr = c_loc(error)
         nothrow%handler = c_funloc(plumed_f_f08_eh)
       else
         nothrow%ptr = c_null_ptr
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_ptr(val,nelem,pass_shape,flags,c_null_ptr),nothrow)
     end subroutine plumed_f_cmd_ptr

     subroutine plumed_f_cmd_char(p,key, val, const, nocopy, error)
       type(cplumed),                 intent(in)    :: p
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*),  target, intent(in) :: val
       logical,                       intent(in)    :: const
       logical,                       intent(in)    :: nocopy
       type(plumed_error), optional, target, intent(out)   :: error
       integer(kind=c_size_t) :: pass_shape(2)
       integer(kind=c_size_t) :: flags
       type(cplumed_nothrow_handler) :: nothrow
       integer(kind=c_size_t) :: nelem
       pass_shape=[len(val) ,0]
       nelem=0
       flags=flags_ptr
       if (const) flags=flags_const_ptr
       if (nocopy) flags=flags + flags_nocopy
       if(present(error)) then
         nothrow%ptr=c_loc(error)
         nothrow%handler = c_funloc(plumed_f_f08_eh)
       else
         nothrow%ptr = c_null_ptr
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_char(c_loc(val),nelem,pass_shape,flags,c_null_ptr),nothrow)
     end subroutine plumed_f_cmd_char

    subroutine plumed_f_cmd_integer_0(p, key, valptr, valshape, const, nocopy, error)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      type(c_ptr),                   intent(in)    :: valptr
      integer,                       intent(in)    :: valshape(:)
      logical,                       intent(in)    :: const
      logical,                       intent(in)    :: nocopy
      type(plumed_error), optional, target, intent(out)   :: error

      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(max(size(valshape) + 1, 2))

      if (size(valshape) == 0) then
        pass_shape(:) = [1, 0]
      else
        pass_shape(:) = [valshape(size(valshape):1:-1), 0]
      endif
      nelem=0
      flags = flags_ptr
      if (const) flags = flags_const_ptr
      if(nocopy) flags = flags + flags_nocopy
      if(present(error)) then
        nothrow%ptr = c_loc(error)
        nothrow%handler = c_funloc(plumed_f_f08_eh)
      else
        nothrow%ptr = c_null_ptr
        nothrow%handler=c_null_funptr
      endif
      call plumed_cmd_safe_nothrow(p, key, &
         plumed_f_safeptr_int(valptr, nelem, pass_shape, flags, c_null_ptr), nothrow)

    end subroutine plumed_f_cmd_integer_0
    subroutine plumed_f_cmd_integer_1(p, key, valptr, valshape, const, nocopy, error)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      type(c_ptr),                   intent(in)    :: valptr
      integer,                       intent(in)    :: valshape(:)
      logical,                       intent(in)    :: const
      logical,                       intent(in)    :: nocopy
      type(plumed_error), optional, target, intent(out)   :: error

      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(max(size(valshape) + 1, 2))

      if (size(valshape) == 0) then
        pass_shape(:) = [1, 0]
      else
        pass_shape(:) = [valshape(size(valshape):1:-1), 0]
      endif
      nelem=0
      flags = flags_ptr
      if (const) flags = flags_const_ptr
      if(nocopy) flags = flags + flags_nocopy
      if(present(error)) then
        nothrow%ptr = c_loc(error)
        nothrow%handler = c_funloc(plumed_f_f08_eh)
      else
        nothrow%ptr = c_null_ptr
        nothrow%handler=c_null_funptr
      endif
      call plumed_cmd_safe_nothrow(p, key, &
         plumed_f_safeptr_short(valptr, nelem, pass_shape, flags, c_null_ptr), nothrow)

    end subroutine plumed_f_cmd_integer_1
    subroutine plumed_f_cmd_integer_2(p, key, valptr, valshape, const, nocopy, error)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      type(c_ptr),                   intent(in)    :: valptr
      integer,                       intent(in)    :: valshape(:)
      logical,                       intent(in)    :: const
      logical,                       intent(in)    :: nocopy
      type(plumed_error), optional, target, intent(out)   :: error

      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(max(size(valshape) + 1, 2))

      if (size(valshape) == 0) then
        pass_shape(:) = [1, 0]
      else
        pass_shape(:) = [valshape(size(valshape):1:-1), 0]
      endif
      nelem=0
      flags = flags_ptr
      if (const) flags = flags_const_ptr
      if(nocopy) flags = flags + flags_nocopy
      if(present(error)) then
        nothrow%ptr = c_loc(error)
        nothrow%handler = c_funloc(plumed_f_f08_eh)
      else
        nothrow%ptr = c_null_ptr
        nothrow%handler=c_null_funptr
      endif
      call plumed_cmd_safe_nothrow(p, key, &
         plumed_f_safeptr_long(valptr, nelem, pass_shape, flags, c_null_ptr), nothrow)

    end subroutine plumed_f_cmd_integer_2
    subroutine plumed_f_cmd_real_0(p, key, valptr, valshape, const, nocopy, error)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      type(c_ptr),                   intent(in)    :: valptr
      integer,                       intent(in)    :: valshape(:)
      logical,                       intent(in)    :: const
      logical,                       intent(in)    :: nocopy
      type(plumed_error), optional, target, intent(out)   :: error

      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(max(size(valshape) + 1, 2))

      if (size(valshape) == 0) then
        pass_shape(:) = [1, 0]
      else
        pass_shape(:) = [valshape(size(valshape):1:-1), 0]
      endif
      nelem=0
      flags = flags_ptr
      if (const) flags = flags_const_ptr
      if(nocopy) flags = flags + flags_nocopy
      if(present(error)) then
        nothrow%ptr = c_loc(error)
        nothrow%handler = c_funloc(plumed_f_f08_eh)
      else
        nothrow%ptr = c_null_ptr
        nothrow%handler=c_null_funptr
      endif
      call plumed_cmd_safe_nothrow(p, key, &
         plumed_f_safeptr_float(valptr, nelem, pass_shape, flags, c_null_ptr), nothrow)

    end subroutine plumed_f_cmd_real_0
    subroutine plumed_f_cmd_real_1(p, key, valptr, valshape, const, nocopy, error)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      type(c_ptr),                   intent(in)    :: valptr
      integer,                       intent(in)    :: valshape(:)
      logical,                       intent(in)    :: const
      logical,                       intent(in)    :: nocopy
      type(plumed_error), optional, target, intent(out)   :: error

      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(max(size(valshape) + 1, 2))

      if (size(valshape) == 0) then
        pass_shape(:) = [1, 0]
      else
        pass_shape(:) = [valshape(size(valshape):1:-1), 0]
      endif
      nelem=0
      flags = flags_ptr
      if (const) flags = flags_const_ptr
      if(nocopy) flags = flags + flags_nocopy
      if(present(error)) then
        nothrow%ptr = c_loc(error)
        nothrow%handler = c_funloc(plumed_f_f08_eh)
      else
        nothrow%ptr = c_null_ptr
        nothrow%handler=c_null_funptr
      endif
      call plumed_cmd_safe_nothrow(p, key, &
         plumed_f_safeptr_double(valptr, nelem, pass_shape, flags, c_null_ptr), nothrow)

    end subroutine plumed_f_cmd_real_1

     ! this is a soft wrapper to a C function
     function plumed_installed() result(res)
       logical             :: res
       res=cplumed_installed()>0
     end function plumed_installed

     ! this is the constructor
     impure elemental subroutine plumed_create(this,kernel)
       type(plumed),    intent(out)          :: this
       character(len=*), intent(in), optional :: kernel
       if(present(kernel)) then
         this%handle=cplumed_create_dlopen(kernel // c_null_char)
       else
         this%handle=cplumed_create()
       endif
       this%initialized=.true.
     end subroutine plumed_create

     ! then we define all member functions

     impure elemental subroutine pl_finalize(this)
       class(plumed), intent(inout) :: this
       if(this%initialized) then
         call cplumed_finalize(this%handle)
         this%initialized=.false.
       endif
     end subroutine pl_finalize

     impure elemental subroutine pl_incref(this)
       class(plumed), intent(inout) :: this
       type(cplumed) :: that
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       that=cplumed_create_reference(this%handle)
     end subroutine pl_incref

     impure elemental subroutine pl_decref(this,to)
       class(plumed),     intent(inout) :: this
       integer, optional, intent(in)    :: to
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if(present(to)) then
         do while(this%use_count()>to)
           call cplumed_finalize(this%handle)
         end do
       else
         call cplumed_finalize(this%handle)
       endif
     end subroutine pl_decref

     ! "impure elemental" needed for the destructor to work on arrays
     impure elemental subroutine pl_destructor(this)
       type(plumed), intent(inout) :: this
       call this%finalize()
     end subroutine pl_destructor

     impure elemental function pl_valid(this) result(valid)
       class(plumed), intent(inout) :: this
       logical :: valid
       integer(c_int) :: i
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valid=cplumed_valid(this%handle)>0
     end function pl_valid

     impure elemental function pl_use_count(this) result(use_count)
       class(plumed), intent(inout) :: this
       integer(c_int) :: use_count
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       use_count=cplumed_use_count(this%handle)
     end function pl_use_count

     impure elemental subroutine pl_assign(this,that)
       class(plumed),intent(out) :: this
       class(plumed),intent(in)  :: that
       if(that%initialized) then
         this%handle=cplumed_create_reference(that%handle)
         this%initialized=.true.
       endif
     end subroutine pl_assign

     impure elemental subroutine pl_cmd(this,key,error)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       type(plumed_error), optional,  intent(out)   :: error
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_ptr(this%handle, key // c_null_char, c_null_ptr, const=.false.,&
           & nocopy=.false., error=error)
     end subroutine pl_cmd

     subroutine pl_cmd_val_char(this,key,val,error)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*),  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_char(this%handle,key // c_null_char,val // c_null_char, &
           & const=.true., nocopy=.true., error=error)
     end subroutine pl_cmd_val_char

     subroutine pl_cmd_ptr_c(this,key,val,error)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       type(c_ptr),                        value    :: val
       type(plumed_error), optional,  intent(out)   :: error
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_ptr(this%handle,key // c_null_char, val, const=.false., nocopy=.false.,&
           & error=error)
     end subroutine pl_cmd_ptr_c

     subroutine pl_cmd_const_ptr_c(this,key,val,error)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       type(c_ptr),                        value    :: val
       type(plumed_error), optional,  intent(out)   :: error
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_ptr(this%handle,key // c_null_char,val, const=.true., nocopy=.false.,&
           & error=error)
     end subroutine pl_cmd_const_ptr_c


    subroutine pl_cmd_val_integer_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),  target, intent(in) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_0_0

    subroutine pl_cmd_ref_integer_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),  target, intent(inout) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
        call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_0_0

    subroutine pl_cmd_ptr_integer_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_0_0

    subroutine pl_cmd_const_ptr_integer_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_0_0


    subroutine pl_cmd_val_integer_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(in) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_0_1

    subroutine pl_cmd_ref_integer_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(inout) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_0_1

    subroutine pl_cmd_ptr_integer_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_0_1

    subroutine pl_cmd_const_ptr_integer_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_0_1


    subroutine pl_cmd_val_integer_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(in) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_0_2

    subroutine pl_cmd_ref_integer_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(inout) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_0_2

    subroutine pl_cmd_ptr_integer_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_0_2

    subroutine pl_cmd_const_ptr_integer_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_0_2


    subroutine pl_cmd_val_integer_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(in) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_0_3

    subroutine pl_cmd_ref_integer_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(inout) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_0_3

    subroutine pl_cmd_ptr_integer_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_0_3

    subroutine pl_cmd_const_ptr_integer_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_0_3


    subroutine pl_cmd_val_integer_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(in) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_0_4

    subroutine pl_cmd_ref_integer_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int), contiguous, target, intent(inout) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_0_4

    subroutine pl_cmd_ptr_integer_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_0_4

    subroutine pl_cmd_const_ptr_integer_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_0_4


    subroutine pl_cmd_val_integer_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),  target, intent(in) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_1_0

    subroutine pl_cmd_ref_integer_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),  target, intent(inout) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
        call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_1_0

    subroutine pl_cmd_ptr_integer_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_1_0

    subroutine pl_cmd_const_ptr_integer_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_1_0


    subroutine pl_cmd_val_integer_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(in) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_1_1

    subroutine pl_cmd_ref_integer_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(inout) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_1_1

    subroutine pl_cmd_ptr_integer_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_1_1

    subroutine pl_cmd_const_ptr_integer_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_1_1


    subroutine pl_cmd_val_integer_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(in) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_1_2

    subroutine pl_cmd_ref_integer_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(inout) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_1_2

    subroutine pl_cmd_ptr_integer_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_1_2

    subroutine pl_cmd_const_ptr_integer_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_1_2


    subroutine pl_cmd_val_integer_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(in) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_1_3

    subroutine pl_cmd_ref_integer_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(inout) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_1_3

    subroutine pl_cmd_ptr_integer_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_1_3

    subroutine pl_cmd_const_ptr_integer_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_1_3


    subroutine pl_cmd_val_integer_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(in) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_1_4

    subroutine pl_cmd_ref_integer_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short), contiguous, target, intent(inout) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_1_4

    subroutine pl_cmd_ptr_integer_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_1_4

    subroutine pl_cmd_const_ptr_integer_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_1_4


    subroutine pl_cmd_val_integer_2_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),  target, intent(in) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_2_0

    subroutine pl_cmd_ref_integer_2_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),  target, intent(inout) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
        call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_2_0

    subroutine pl_cmd_ptr_integer_2_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_2_0

    subroutine pl_cmd_const_ptr_integer_2_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_2_0


    subroutine pl_cmd_val_integer_2_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(in) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_2_1

    subroutine pl_cmd_ref_integer_2_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(inout) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_2_1

    subroutine pl_cmd_ptr_integer_2_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_2_1

    subroutine pl_cmd_const_ptr_integer_2_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_2_1


    subroutine pl_cmd_val_integer_2_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(in) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_2_2

    subroutine pl_cmd_ref_integer_2_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(inout) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_2_2

    subroutine pl_cmd_ptr_integer_2_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_2_2

    subroutine pl_cmd_const_ptr_integer_2_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_2_2


    subroutine pl_cmd_val_integer_2_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(in) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_2_3

    subroutine pl_cmd_ref_integer_2_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(inout) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_2_3

    subroutine pl_cmd_ptr_integer_2_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_2_3

    subroutine pl_cmd_const_ptr_integer_2_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_2_3


    subroutine pl_cmd_val_integer_2_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(in) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_integer_2_4

    subroutine pl_cmd_ref_integer_2_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long), contiguous, target, intent(inout) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_integer_2_4

    subroutine pl_cmd_ptr_integer_2_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_integer_2_4

    subroutine pl_cmd_const_ptr_integer_2_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_integer_2(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_integer_2_4


    subroutine pl_cmd_val_real_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),  target, intent(in) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_0_0

    subroutine pl_cmd_ref_real_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),  target, intent(inout) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
        call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_0_0

    subroutine pl_cmd_ptr_real_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_0_0

    subroutine pl_cmd_const_ptr_real_0_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_0_0


    subroutine pl_cmd_val_real_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(in) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_0_1

    subroutine pl_cmd_ref_real_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(inout) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_0_1

    subroutine pl_cmd_ptr_real_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_0_1

    subroutine pl_cmd_const_ptr_real_0_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_0_1


    subroutine pl_cmd_val_real_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(in) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_0_2

    subroutine pl_cmd_ref_real_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(inout) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_0_2

    subroutine pl_cmd_ptr_real_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_0_2

    subroutine pl_cmd_const_ptr_real_0_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_0_2


    subroutine pl_cmd_val_real_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(in) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_0_3

    subroutine pl_cmd_ref_real_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(inout) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_0_3

    subroutine pl_cmd_ptr_real_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_0_3

    subroutine pl_cmd_const_ptr_real_0_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_0_3


    subroutine pl_cmd_val_real_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(in) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_0_4

    subroutine pl_cmd_ref_real_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float), contiguous, target, intent(inout) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_0_4

    subroutine pl_cmd_ptr_real_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_0_4

    subroutine pl_cmd_const_ptr_real_0_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_0(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_0_4


    subroutine pl_cmd_val_real_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),  target, intent(in) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_1_0

    subroutine pl_cmd_ref_real_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),  target, intent(inout) :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
        call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_1_0

    subroutine pl_cmd_ptr_real_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_1_0

    subroutine pl_cmd_const_ptr_real_1_0(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(0)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_1_0


    subroutine pl_cmd_val_real_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(in) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_1_1

    subroutine pl_cmd_ref_real_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(inout) :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_1_1

    subroutine pl_cmd_ptr_real_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_1_1

    subroutine pl_cmd_const_ptr_real_1_1(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(1)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_1_1


    subroutine pl_cmd_val_real_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(in) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_1_2

    subroutine pl_cmd_ref_real_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(inout) :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_1_2

    subroutine pl_cmd_ptr_real_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_1_2

    subroutine pl_cmd_const_ptr_real_1_2(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(2)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_1_2


    subroutine pl_cmd_val_real_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(in) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_1_3

    subroutine pl_cmd_ref_real_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(inout) :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_1_3

    subroutine pl_cmd_ptr_real_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_1_3

    subroutine pl_cmd_const_ptr_real_1_3(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(3)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_1_3


    subroutine pl_cmd_val_real_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(in) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
      valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.true., error=error)
    end subroutine pl_cmd_val_real_1_4

    subroutine pl_cmd_ref_real_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double), contiguous, target, intent(inout) :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       valshape = shape(val)
        call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
            & const=.false., nocopy=.true., error=error)
    end subroutine pl_cmd_ref_real_1_4

    subroutine pl_cmd_ptr_real_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained pointer to non-contiguous data"
           return
         else
           error stop "Obtained pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.false., nocopy=.false., error=error)
    end subroutine pl_cmd_ptr_real_1_4

    subroutine pl_cmd_const_ptr_real_1_4(this,key,val,error)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double),   pointer,  intent(in)    :: val(:,:,:,:)
       type(plumed_error), optional,  intent(out)   :: error
       integer :: valshape(4)
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       if (.not. is_contiguous(val)) then
         if (present(error)) then
           error%code=20300
           error%what="Obtained const pointer to non-contiguous data"
           return
         else
           error stop "Obtained const pointer to non-contiguous data"
         end if
       end if
       valshape = shape(val)
       call plumed_f_cmd_real_1(this%handle, key // c_null_char, c_loc(val), valshape,&
           & const=.true., nocopy=.false., error=error)
    end subroutine pl_cmd_const_ptr_real_1_4


     ! copy, including deep copy of nested errors
     ! this should be fortran default copy but
     ! fails on gfortran, so it is reimplemented here
     ! see https://godbolt.org/z/9dM6vj5bo
     impure elemental subroutine pl_error_assign(this,that)
       class(plumed_error),target,intent(out) :: this
       class(plumed_error),target,intent(in)  :: that
       type(plumed_error), pointer :: that_ptr, this_ptr
       that_ptr => that
       this_ptr => this
       do
         this_ptr%code = that_ptr%code
         this_ptr%what = that_ptr%what
         if(.not. allocated(that_ptr%nested)) exit
         allocate(this_ptr%nested)
         this_ptr => this_ptr%nested
         that_ptr => that_ptr%nested
       end do
     end subroutine pl_error_assign

end module plumed_f08_module
