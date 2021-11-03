! vim:ft=fortran



module plumed_module_f08
  use iso_c_binding
  implicit none

  ! names are private by default
  private

  ! only these names are public
  public :: plumed
  public :: plumed_create
  public :: plumed_installed
  public :: plumed_error

  ! used to enforce keyword-only arguments
  type dummy_type
  end type dummy_type

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
    generic, public :: cmd => &
    pl_cmd_integer_0_0, &
    pl_cmd_integer_0_1, &
    pl_cmd_integer_0_2, &
    pl_cmd_integer_0_3, &
    pl_cmd_integer_0_4, &
    pl_cmd_integer_1_0, &
    pl_cmd_integer_1_1, &
    pl_cmd_integer_1_2, &
    pl_cmd_integer_1_3, &
    pl_cmd_integer_1_4, &
    pl_cmd_integer_2_0, &
    pl_cmd_integer_2_1, &
    pl_cmd_integer_2_2, &
    pl_cmd_integer_2_3, &
    pl_cmd_integer_2_4, &
    pl_cmd_real_0_0, &
    pl_cmd_real_0_1, &
    pl_cmd_real_0_2, &
    pl_cmd_real_0_3, &
    pl_cmd_real_0_4, &
    pl_cmd_real_1_0, &
    pl_cmd_real_1_1, &
    pl_cmd_real_1_2, &
    pl_cmd_real_1_3, &
    pl_cmd_real_1_4, &
    pl_cmd_real_2_0, &
    pl_cmd_real_2_1, &
    pl_cmd_real_2_2, &
    pl_cmd_real_2_3, &
    pl_cmd_real_2_4, &
    pl_cmd, &
    pl_cmd_char, &
    pl_cmd_ptr

    procedure :: pl_cmd_integer_0_0
    procedure :: pl_cmd_integer_0_1
    procedure :: pl_cmd_integer_0_2
    procedure :: pl_cmd_integer_0_3
    procedure :: pl_cmd_integer_0_4
    procedure :: pl_cmd_integer_1_0
    procedure :: pl_cmd_integer_1_1
    procedure :: pl_cmd_integer_1_2
    procedure :: pl_cmd_integer_1_3
    procedure :: pl_cmd_integer_1_4
    procedure :: pl_cmd_integer_2_0
    procedure :: pl_cmd_integer_2_1
    procedure :: pl_cmd_integer_2_2
    procedure :: pl_cmd_integer_2_3
    procedure :: pl_cmd_integer_2_4
    procedure :: pl_cmd_real_0_0
    procedure :: pl_cmd_real_0_1
    procedure :: pl_cmd_real_0_2
    procedure :: pl_cmd_real_0_3
    procedure :: pl_cmd_real_0_4
    procedure :: pl_cmd_real_1_0
    procedure :: pl_cmd_real_1_1
    procedure :: pl_cmd_real_1_2
    procedure :: pl_cmd_real_1_3
    procedure :: pl_cmd_real_1_4
    procedure :: pl_cmd_real_2_0
    procedure :: pl_cmd_real_2_1
    procedure :: pl_cmd_real_2_2
    procedure :: pl_cmd_real_2_3
    procedure :: pl_cmd_real_2_4
    procedure :: pl_cmd
    procedure :: pl_cmd_char
    procedure :: pl_cmd_ptr

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
  interface
    function plumed_f_safeptr_long_double_scalar(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_long_double)                     :: val
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_long_double_scalar
    end function plumed_f_safeptr_long_double_scalar
  end interface
  interface
    function plumed_f_safeptr_long_double(val,nelem,pass_shape,flags,opt) bind(C)
      import
      real(kind=c_long_double)                     :: val(*)
      integer(kind=c_size_t),         value :: nelem
      integer(kind=c_size_t)                :: pass_shape(*)
      integer(kind=c_size_t),         value :: flags
      type(c_ptr),                    value :: opt
      type(cplumed_safeptr)                 :: plumed_f_safeptr_long_double
    end function plumed_f_safeptr_long_double
  end interface

  ! here are the interfaces used for overloading
  interface plumed_f_cmd
    module procedure plumed_f_cmd_ptr
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
    module procedure plumed_f_cmd_real_2_0
    module procedure plumed_f_cmd_real_2_1
    module procedure plumed_f_cmd_real_2_2
    module procedure plumed_f_cmd_real_2_3
    module procedure plumed_f_cmd_real_2_4
  end interface plumed_f_cmd

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
     end subroutine plumed_f_f08_eh

     ! we then define all the functions needed for overloading

     subroutine plumed_f_cmd_ptr(p,key,val,dummy,error,const,nocopy)
       type(cplumed),                 intent(in)    :: p
       character(kind=c_char,len=*),  intent(in)    :: key
       type(c_ptr),                     value       :: val
       type(dummy_type),   optional                 :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       type(plumed_error), target :: myerror
       integer(kind=c_size_t) :: pass_shape(1)
       integer(kind=c_size_t) :: flags
       type(cplumed_nothrow_handler) :: nothrow
       integer(kind=c_size_t) :: nelem
       pass_shape=(/0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_ptr(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
     end subroutine plumed_f_cmd_ptr

     subroutine plumed_f_cmd_char(p,key,val,dummy,error,const,nocopy)
       type(cplumed),                 intent(in)    :: p
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*)                 :: val
       type(dummy_type),   optional                 :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       type(plumed_error), target :: myerror
       integer(kind=c_size_t) :: pass_shape(2)
       integer(kind=c_size_t) :: flags
       type(cplumed_nothrow_handler) :: nothrow
       integer(kind=c_size_t) :: nelem
       pass_shape=(/len(val),0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_char(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
     end subroutine plumed_f_cmd_char

    subroutine plumed_f_cmd_integer_0_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_int_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_0_0
    subroutine plumed_f_cmd_integer_0_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_0_1
    subroutine plumed_f_cmd_integer_0_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_0_2
    subroutine plumed_f_cmd_integer_0_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_0_3
    subroutine plumed_f_cmd_integer_0_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_int(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_0_4
    subroutine plumed_f_cmd_integer_1_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_short_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_1_0
    subroutine plumed_f_cmd_integer_1_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_1_1
    subroutine plumed_f_cmd_integer_1_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_1_2
    subroutine plumed_f_cmd_integer_1_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_1_3
    subroutine plumed_f_cmd_integer_1_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_short(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_1_4
    subroutine plumed_f_cmd_integer_2_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_2_0
    subroutine plumed_f_cmd_integer_2_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_2_1
    subroutine plumed_f_cmd_integer_2_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_2_2
    subroutine plumed_f_cmd_integer_2_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_2_3
    subroutine plumed_f_cmd_integer_2_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_integer_2_4
    subroutine plumed_f_cmd_real_0_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_float_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_0_0
    subroutine plumed_f_cmd_real_0_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_0_1
    subroutine plumed_f_cmd_real_0_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_0_2
    subroutine plumed_f_cmd_real_0_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_0_3
    subroutine plumed_f_cmd_real_0_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_float(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_0_4
    subroutine plumed_f_cmd_real_1_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_double_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_1_0
    subroutine plumed_f_cmd_real_1_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_1_1
    subroutine plumed_f_cmd_real_1_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_1_2
    subroutine plumed_f_cmd_real_1_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_1_3
    subroutine plumed_f_cmd_real_1_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_1_4
    subroutine plumed_f_cmd_real_2_0(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape=(/1,0/)
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_double_scalar(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_2_0
    subroutine plumed_f_cmd_real_2_1(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(2)
      pass_shape(1)=size(val,1)
      pass_shape(2)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_2_1
    subroutine plumed_f_cmd_real_2_2(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(3)
      pass_shape(1)=size(val,2)
      pass_shape(2)=size(val,1)
      pass_shape(3)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_2_2
    subroutine plumed_f_cmd_real_2_3(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(4)
      pass_shape(1)=size(val,3)
      pass_shape(2)=size(val,2)
      pass_shape(3)=size(val,1)
      pass_shape(4)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_2_3
    subroutine plumed_f_cmd_real_2_4(p,key,val,dummy,error,const,nocopy)
      type(cplumed),                 intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:,:,:)
       type(dummy_type),   optional                 :: dummy
      type(plumed_error), optional,  intent(out)   :: error
      logical,            optional,  intent(in)    :: const
      logical,            optional,  intent(in)    :: nocopy
      type(plumed_error), target :: myerror
      integer(kind=c_size_t) :: flags
      type(cplumed_nothrow_handler) :: nothrow
      integer(kind=c_size_t) :: nelem
      integer(kind=c_size_t) :: pass_shape(5)
      pass_shape(1)=size(val,4)
      pass_shape(2)=size(val,3)
      pass_shape(3)=size(val,2)
      pass_shape(4)=size(val,1)
      pass_shape(5)=0
       nelem=product(pass_shape)
       flags=flags_ptr
       if(present(const)) then
         if(const) then
           flags=flags_const_ptr
         endif
       endif
       if(present(nocopy)) then
         if(nocopy) then
           flags=flags+flags_nocopy
         endif
       endif
       nothrow%ptr=c_loc(myerror)
       if(present(error)) then
         nothrow%handler=c_funloc(plumed_f_f08_eh)
       else
         nothrow%handler=c_null_funptr
       endif
       call plumed_cmd_safe_nothrow(p,key, &
         plumed_f_safeptr_long_double(val,nelem,pass_shape,flags_ptr,c_null_ptr),nothrow)
       if(present(error)) then
         error=myerror
       endif
    end subroutine plumed_f_cmd_real_2_4

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

     impure elemental subroutine pl_cmd(this,key,dummy,error,const,nocopy)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       type(dummy_type),   optional,  intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,c_null_ptr,error=error,const=const,nocopy=nocopy)
     end subroutine pl_cmd

     subroutine pl_cmd_char(this,key,val,dummy,error,const,nocopy)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       character(kind=c_char,len=*)                 :: val
       type(dummy_type),   optional,  intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val // c_null_char,error=error,const=const,nocopy=nocopy)
     end subroutine pl_cmd_char

     subroutine pl_cmd_ptr(this,key,val,dummy,error,const,nocopy)
       class(plumed),                 intent(inout) :: this ! inout to allow for initialization
       character(kind=c_char,len=*),  intent(in)    :: key
       type(c_ptr),                        value    :: val
       type(dummy_type),   optional,  intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
     end subroutine pl_cmd_ptr

    subroutine pl_cmd_integer_0_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_0_0
    subroutine pl_cmd_integer_0_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_0_1
    subroutine pl_cmd_integer_0_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_0_2
    subroutine pl_cmd_integer_0_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_0_3
    subroutine pl_cmd_integer_0_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_int)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_0_4
    subroutine pl_cmd_integer_1_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_1_0
    subroutine pl_cmd_integer_1_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_1_1
    subroutine pl_cmd_integer_1_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_1_2
    subroutine pl_cmd_integer_1_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_1_3
    subroutine pl_cmd_integer_1_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_short)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_1_4
    subroutine pl_cmd_integer_2_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_2_0
    subroutine pl_cmd_integer_2_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_2_1
    subroutine pl_cmd_integer_2_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_2_2
    subroutine pl_cmd_integer_2_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_2_3
    subroutine pl_cmd_integer_2_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      integer(KIND=c_long)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_integer_2_4
    subroutine pl_cmd_real_0_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_0_0
    subroutine pl_cmd_real_0_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_0_1
    subroutine pl_cmd_real_0_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_0_2
    subroutine pl_cmd_real_0_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_0_3
    subroutine pl_cmd_real_0_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_float)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_0_4
    subroutine pl_cmd_real_1_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_1_0
    subroutine pl_cmd_real_1_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_1_1
    subroutine pl_cmd_real_1_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_1_2
    subroutine pl_cmd_real_1_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_1_3
    subroutine pl_cmd_real_1_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_double)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_1_4
    subroutine pl_cmd_real_2_0(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_2_0
    subroutine pl_cmd_real_2_1(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_2_1
    subroutine pl_cmd_real_2_2(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_2_2
    subroutine pl_cmd_real_2_3(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_2_3
    subroutine pl_cmd_real_2_4(this,key,val,dummy,error,const,nocopy)
      class(plumed),                 intent(inout) :: this ! inout to allow for initialization
      character(kind=c_char,len=*),  intent(in)    :: key
      real(KIND=c_long_double)                            :: val(:,:,:,:)
       type(dummy_type),   optional, intent(inout) :: dummy
       type(plumed_error), optional,  intent(out)   :: error
       logical,            optional,  intent(in)    :: const
       logical,            optional,  intent(in)    :: nocopy
       if(.not.this%initialized) then
         call plumed_create(this)
       endif
       call plumed_f_cmd(this%handle,key // c_null_char,val,error=error,const=const,nocopy=nocopy)
    end subroutine pl_cmd_real_2_4

end module plumed_module_f08

