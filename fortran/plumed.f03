! vim:ft=fortran
module plumed_module
  use iso_c_binding
  implicit none

  private
  public :: plumed, plumed_eos, plumed_error
  public :: plumed_f03_create, plumed_f03_create_dlopen, plumed_f03_create_reference, plumed_f03_create_invalid
  public :: plumed_f03_cmd
  public :: plumed_f03_finalize
  public :: plumed_f03_installed
  public :: plumed_f03_valid, plumed_f03_use_count
  public :: plumed_f03_global
  public :: plumed_f03_ginitialized
  public :: plumed_f03_gcreate
  public :: plumed_f03_gcmd
  public :: plumed_f03_gfinalize
  public :: plumed_f03_gvalid

  type :: plumed
    character(32) :: ptr
  end type plumed
  character, parameter :: plumed_eos=achar(0)

  type :: plumed_error
    integer                         :: code
    character(len = :), allocatable :: what
  end type plumed_error


  interface plumed_f03_cmd
    module procedure plumed_f03_cmd_null
    module procedure plumed_f03_cmd_char
    module procedure plumed_f03_cmd_integer_0_0
    module procedure plumed_f03_cmd_integer_0_1
    module procedure plumed_f03_cmd_integer_0_2
    module procedure plumed_f03_cmd_integer_0_3
    module procedure plumed_f03_cmd_integer_0_4
    module procedure plumed_f03_cmd_integer_1_0
    module procedure plumed_f03_cmd_integer_1_1
    module procedure plumed_f03_cmd_integer_1_2
    module procedure plumed_f03_cmd_integer_1_3
    module procedure plumed_f03_cmd_integer_1_4
    module procedure plumed_f03_cmd_integer_2_0
    module procedure plumed_f03_cmd_integer_2_1
    module procedure plumed_f03_cmd_integer_2_2
    module procedure plumed_f03_cmd_integer_2_3
    module procedure plumed_f03_cmd_integer_2_4
    module procedure plumed_f03_cmd_real_0_0
    module procedure plumed_f03_cmd_real_0_1
    module procedure plumed_f03_cmd_real_0_2
    module procedure plumed_f03_cmd_real_0_3
    module procedure plumed_f03_cmd_real_0_4
    module procedure plumed_f03_cmd_real_1_0
    module procedure plumed_f03_cmd_real_1_1
    module procedure plumed_f03_cmd_real_1_2
    module procedure plumed_f03_cmd_real_1_3
    module procedure plumed_f03_cmd_real_1_4
    module procedure plumed_f03_cmd_real_2_0
    module procedure plumed_f03_cmd_real_2_1
    module procedure plumed_f03_cmd_real_2_2
    module procedure plumed_f03_cmd_real_2_3
    module procedure plumed_f03_cmd_real_2_4
  end interface plumed_f03_cmd
  interface plumed_f03_gcmd
    module procedure plumed_f03_gcmd_null
    module procedure plumed_f03_gcmd_char
    module procedure plumed_f03_gcmd_integer_0_0
    module procedure plumed_f03_gcmd_integer_0_1
    module procedure plumed_f03_gcmd_integer_0_2
    module procedure plumed_f03_gcmd_integer_0_3
    module procedure plumed_f03_gcmd_integer_0_4
    module procedure plumed_f03_gcmd_integer_1_0
    module procedure plumed_f03_gcmd_integer_1_1
    module procedure plumed_f03_gcmd_integer_1_2
    module procedure plumed_f03_gcmd_integer_1_3
    module procedure plumed_f03_gcmd_integer_1_4
    module procedure plumed_f03_gcmd_integer_2_0
    module procedure plumed_f03_gcmd_integer_2_1
    module procedure plumed_f03_gcmd_integer_2_2
    module procedure plumed_f03_gcmd_integer_2_3
    module procedure plumed_f03_gcmd_integer_2_4
    module procedure plumed_f03_gcmd_real_0_0
    module procedure plumed_f03_gcmd_real_0_1
    module procedure plumed_f03_gcmd_real_0_2
    module procedure plumed_f03_gcmd_real_0_3
    module procedure plumed_f03_gcmd_real_0_4
    module procedure plumed_f03_gcmd_real_1_0
    module procedure plumed_f03_gcmd_real_1_1
    module procedure plumed_f03_gcmd_real_1_2
    module procedure plumed_f03_gcmd_real_1_3
    module procedure plumed_f03_gcmd_real_1_4
    module procedure plumed_f03_gcmd_real_2_0
    module procedure plumed_f03_gcmd_real_2_1
    module procedure plumed_f03_gcmd_real_2_2
    module procedure plumed_f03_gcmd_real_2_3
    module procedure plumed_f03_gcmd_real_2_4
  end interface plumed_f03_gcmd

  contains

    subroutine error_handler(error,code,what,length)
      type(plumed_error), intent(out) :: error
      integer,            intent(in)  :: code
      character,          intent(in)  :: what(*)
      integer,            intent(in)  :: length
      integer :: i
      error%code=code
      allocate(character(length) :: error%what)
      do i=1,length
        error%what(i:i)=what(i)
      end do
    end subroutine error_handler

    function plumed_f03_create() result(p)
      type(plumed) :: p
      call plumed_f_create(p%ptr)
    end function plumed_f03_create

    function plumed_f03_create_dlopen(path) result(p)
      character(len=*), intent(in) :: path
      type(plumed)                 :: p
      call plumed_f_create_dlopen(trim(path)//plumed_eos,p%ptr)
    end function plumed_f03_create_dlopen

    function plumed_f03_create_reference(ref) result(p)
      type(plumed), intent(in):: ref
      type(plumed):: p
      call plumed_f_create_reference(ref%ptr,p%ptr)
    end function plumed_f03_create_reference

    function plumed_f03_create_invalid() result(p)
      type(plumed):: p
      call plumed_f_create_invalid(p%ptr)
    end function plumed_f03_create_invalid

    subroutine plumed_f03_cmd_null(p,key,error)
      type(plumed)                  :: p
      character(len=*),  intent(in) :: key
      type(plumed_error), optional  :: error
      if(present(error)) then
        error%code=0
        error%what=""
        call plumed_f_cmd_safe_nothrow_null(p%ptr,trim(key)//plumed_eos,error,error_handler)
      else
        call plumed_f_cmd_safe_null(p%ptr,trim(key)//plumed_eos)
      endif
    end subroutine plumed_f03_cmd_null

    subroutine plumed_f03_finalize(p)
      type(plumed), intent(in) :: p
      call plumed_f_finalize(p%ptr)
    end subroutine plumed_f03_finalize

    function plumed_f03_installed() result(installed)
      logical :: installed
      integer :: i
      i=0
      call plumed_f_installed(i)
      installed = i>0
    end function plumed_f03_installed

    function plumed_f03_valid(p) result(valid)
      type(plumed) :: p
      logical      :: valid
      integer      :: i
      i=0
      call plumed_f_valid(p%ptr,i)
      valid = i>0
    end function plumed_f03_valid

    function plumed_f03_use_count(p) result(count)
      type(plumed) :: p
      integer      :: count
      count=0
      call plumed_f_use_count(p%ptr,count)
    end function plumed_f03_use_count

    function plumed_f03_global() result(p)
      type(plumed) :: p
      call plumed_f_global(p%ptr)
    end function plumed_f03_global

    function plumed_f03_ginitialized() result(initialized)
      logical :: initialized
      integer :: i
      call plumed_f_ginitialized(i)
    end function plumed_f03_ginitialized

    subroutine plumed_f03_gcreate()
      call plumed_f_gcreate()
    end subroutine plumed_f03_gcreate

    subroutine plumed_f03_gcmd_null(key,error)
      character(len=*),  intent(in)    :: key
      type(plumed_error), optional  :: error
      call plumed_f03_cmd_null(plumed_f03_global(),key,error)
    end subroutine plumed_f03_gcmd_null

    subroutine plumed_f03_gfinalize()
      call plumed_f_gfinalize()
    end subroutine plumed_f03_gfinalize

    function plumed_f03_gvalid() result(valid)
      logical :: valid
      integer :: i
      i=0
      call plumed_f_gvalid(i)
      valid=i>0
    end function plumed_f03_gvalid

    subroutine plumed_f03_cmd_char(p,key,val,error)
      type(plumed)                     :: p
      character(len=*),  intent(in)    :: key
      character(len=*)                 :: val
      type(plumed_error), optional  :: error
      integer(kind=C_SIZE_T) :: pass_shape(1)
      integer(kind=C_SIZE_T) :: pass_nelem
      pass_shape=(/0/)
      pass_nelem=0
      if(present(error)) then
        error%code=0
        error%what=""
        call plumed_f_cmd_safe_nothrow_char(p%ptr,trim(key)//plumed_eos,val//plumed_eos,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_char(p%ptr,trim(key)//plumed_eos,val//plumed_eos,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_char

    subroutine plumed_f03_gcmd_char(key,val,error)
      character(len=*),  intent(in)    :: key
      character(len=*)                 :: val
      type(plumed_error), optional  :: error
      call plumed_f03_cmd_char(plumed_f03_global(),key,val,error)
    end subroutine plumed_f03_gcmd_char


    subroutine plumed_f03_cmd_integer_0_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_int(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_0_0

    subroutine plumed_f03_gcmd_integer_0_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_0_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_0_0

    subroutine plumed_f03_cmd_integer_0_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_int(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_0_1

    subroutine plumed_f03_gcmd_integer_0_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_0_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_0_1

    subroutine plumed_f03_cmd_integer_0_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_int(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_0_2

    subroutine plumed_f03_gcmd_integer_0_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_0_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_0_2

    subroutine plumed_f03_cmd_integer_0_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_int(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_0_3

    subroutine plumed_f03_gcmd_integer_0_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_0_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_0_3

    subroutine plumed_f03_cmd_integer_0_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_int(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_0_4

    subroutine plumed_f03_gcmd_integer_0_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_0_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_0_4

    subroutine plumed_f03_cmd_integer_1_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_short(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_1_0

    subroutine plumed_f03_gcmd_integer_1_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_1_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_1_0

    subroutine plumed_f03_cmd_integer_1_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_short(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_1_1

    subroutine plumed_f03_gcmd_integer_1_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_1_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_1_1

    subroutine plumed_f03_cmd_integer_1_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_short(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_1_2

    subroutine plumed_f03_gcmd_integer_1_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_1_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_1_2

    subroutine plumed_f03_cmd_integer_1_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_short(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_1_3

    subroutine plumed_f03_gcmd_integer_1_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_1_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_1_3

    subroutine plumed_f03_cmd_integer_1_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_short(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_1_4

    subroutine plumed_f03_gcmd_integer_1_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_1_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_1_4

    subroutine plumed_f03_cmd_integer_2_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_2_0

    subroutine plumed_f03_gcmd_integer_2_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_2_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_2_0

    subroutine plumed_f03_cmd_integer_2_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_2_1

    subroutine plumed_f03_gcmd_integer_2_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_2_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_2_1

    subroutine plumed_f03_cmd_integer_2_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_2_2

    subroutine plumed_f03_gcmd_integer_2_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_2_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_2_2

    subroutine plumed_f03_cmd_integer_2_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_2_3

    subroutine plumed_f03_gcmd_integer_2_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_2_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_2_3

    subroutine plumed_f03_cmd_integer_2_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_integer_2_4

    subroutine plumed_f03_gcmd_integer_2_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_integer_2_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_integer_2_4

    subroutine plumed_f03_cmd_real_0_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_float(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_0_0

    subroutine plumed_f03_gcmd_real_0_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_0_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_0_0

    subroutine plumed_f03_cmd_real_0_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_float(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_0_1

    subroutine plumed_f03_gcmd_real_0_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_0_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_0_1

    subroutine plumed_f03_cmd_real_0_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_float(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_0_2

    subroutine plumed_f03_gcmd_real_0_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_0_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_0_2

    subroutine plumed_f03_cmd_real_0_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_float(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_0_3

    subroutine plumed_f03_gcmd_real_0_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_0_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_0_3

    subroutine plumed_f03_cmd_real_0_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_float(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_0_4

    subroutine plumed_f03_gcmd_real_0_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_0_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_0_4

    subroutine plumed_f03_cmd_real_1_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_1_0

    subroutine plumed_f03_gcmd_real_1_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_1_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_1_0

    subroutine plumed_f03_cmd_real_1_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_1_1

    subroutine plumed_f03_gcmd_real_1_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_1_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_1_1

    subroutine plumed_f03_cmd_real_1_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_1_2

    subroutine plumed_f03_gcmd_real_1_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_1_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_1_2

    subroutine plumed_f03_cmd_real_1_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_1_3

    subroutine plumed_f03_gcmd_real_1_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_1_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_1_3

    subroutine plumed_f03_cmd_real_1_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_1_4

    subroutine plumed_f03_gcmd_real_1_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_1_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_1_4

    subroutine plumed_f03_cmd_real_2_0(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_2_0

    subroutine plumed_f03_gcmd_real_2_0(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_2_0(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_2_0

    subroutine plumed_f03_cmd_real_2_1(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_2_1

    subroutine plumed_f03_gcmd_real_2_1(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_2_1(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_2_1

    subroutine plumed_f03_cmd_real_2_2(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_2_2

    subroutine plumed_f03_gcmd_real_2_2(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_2_2(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_2_2

    subroutine plumed_f03_cmd_real_2_3(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_2_3

    subroutine plumed_f03_gcmd_real_2_3(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_2_3(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_2_3

    subroutine plumed_f03_cmd_real_2_4(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_real_2_4

    subroutine plumed_f03_gcmd_real_2_4(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:,:)
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_real_2_4(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_real_2_4

end module plumed_module
