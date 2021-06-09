! vim:ft=fortran
module plumed_module
  use iso_c_binding
  implicit none

  private
  public :: plumed, plumed_eos
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

    subroutine plumed_f03_cmd_null(p,key)
      type(plumed)                  :: p
      character(len=*),  intent(in) :: key
      call plumed_f_cmd_safe_null(p%ptr,trim(key)//plumed_eos)
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

    subroutine plumed_f03_gcmd_null(key)
      character(len=*),  intent(in)    :: key
      call plumed_f_gcmd_safe_null(trim(key)//plumed_eos)
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

    subroutine plumed_f03_cmd_char(p,key,val)
      type(plumed)                     :: p
      character(len=*),  intent(in)    :: key
      character(len=*)                 :: val
      integer(kind=C_SIZE_T) :: pass_shape(1)
      integer(kind=C_SIZE_T) :: pass_nelem
      pass_shape=(/0/)
      pass_nelem=0
      call plumed_f_cmd_safe_char(p%ptr,trim(key)//plumed_eos,val//plumed_eos,pass_shape(1))
    end subroutine plumed_f03_cmd_char

    subroutine plumed_f03_gcmd_char(key,val)
      character(len=*),  intent(in)    :: key
      character(len=*)                 :: val
      integer(kind=C_SIZE_T) :: pass_shape(1)
      integer(kind=C_SIZE_T) :: pass_nelem
      pass_shape=(/0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_char(trim(key)//plumed_eos,val//plumed_eos,pass_shape(1))
    end subroutine plumed_f03_gcmd_char


    subroutine plumed_f03_cmd_integer_0_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_integer_0_0

    subroutine plumed_f03_gcmd_integer_0_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_int(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_0_0

    subroutine plumed_f03_cmd_integer_0_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_0_1

    subroutine plumed_f03_gcmd_integer_0_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_int(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_0_1

    subroutine plumed_f03_cmd_integer_0_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_0_2

    subroutine plumed_f03_gcmd_integer_0_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_int(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_0_2

    subroutine plumed_f03_cmd_integer_0_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_0_3

    subroutine plumed_f03_gcmd_integer_0_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_int(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_0_3

    subroutine plumed_f03_cmd_integer_0_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_int(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_0_4

    subroutine plumed_f03_gcmd_integer_0_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_int)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_int(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_0_4

    subroutine plumed_f03_cmd_integer_1_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_integer_1_0

    subroutine plumed_f03_gcmd_integer_1_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_short(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_1_0

    subroutine plumed_f03_cmd_integer_1_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_1_1

    subroutine plumed_f03_gcmd_integer_1_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_short(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_1_1

    subroutine plumed_f03_cmd_integer_1_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_1_2

    subroutine plumed_f03_gcmd_integer_1_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_short(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_1_2

    subroutine plumed_f03_cmd_integer_1_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_1_3

    subroutine plumed_f03_gcmd_integer_1_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_short(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_1_3

    subroutine plumed_f03_cmd_integer_1_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_short(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_1_4

    subroutine plumed_f03_gcmd_integer_1_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_short)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_short(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_1_4

    subroutine plumed_f03_cmd_integer_2_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_integer_2_0

    subroutine plumed_f03_gcmd_integer_2_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_long(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_2_0

    subroutine plumed_f03_cmd_integer_2_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_2_1

    subroutine plumed_f03_gcmd_integer_2_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_2_1

    subroutine plumed_f03_cmd_integer_2_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_2_2

    subroutine plumed_f03_gcmd_integer_2_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_2_2

    subroutine plumed_f03_cmd_integer_2_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_2_3

    subroutine plumed_f03_gcmd_integer_2_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_2_3

    subroutine plumed_f03_cmd_integer_2_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_integer_2_4

    subroutine plumed_f03_gcmd_integer_2_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      integer(KIND=c_long)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_integer_2_4

    subroutine plumed_f03_cmd_real_0_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_real_0_0

    subroutine plumed_f03_gcmd_real_0_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_float(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_real_0_0

    subroutine plumed_f03_cmd_real_0_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_0_1

    subroutine plumed_f03_gcmd_real_0_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_float(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_0_1

    subroutine plumed_f03_cmd_real_0_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_0_2

    subroutine plumed_f03_gcmd_real_0_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_float(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_0_2

    subroutine plumed_f03_cmd_real_0_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_0_3

    subroutine plumed_f03_gcmd_real_0_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_float(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_0_3

    subroutine plumed_f03_cmd_real_0_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_float(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_0_4

    subroutine plumed_f03_gcmd_real_0_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_float)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_float(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_0_4

    subroutine plumed_f03_cmd_real_1_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_real_1_0

    subroutine plumed_f03_gcmd_real_1_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_double(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_real_1_0

    subroutine plumed_f03_cmd_real_1_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_1_1

    subroutine plumed_f03_gcmd_real_1_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_double(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_1_1

    subroutine plumed_f03_cmd_real_1_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_1_2

    subroutine plumed_f03_gcmd_real_1_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_double(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_1_2

    subroutine plumed_f03_cmd_real_1_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_1_3

    subroutine plumed_f03_gcmd_real_1_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_double(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_1_3

    subroutine plumed_f03_cmd_real_1_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_1_4

    subroutine plumed_f03_gcmd_real_1_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_double)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_double(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_1_4

    subroutine plumed_f03_cmd_real_2_0(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_cmd_real_2_0

    subroutine plumed_f03_gcmd_real_2_0(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=(/1,0/)
      pass_nelem=0
      call plumed_f_gcmd_safe_long_double(trim(key)//plumed_eos,val,pass_shape(1))
    end subroutine plumed_f03_gcmd_real_2_0

    subroutine plumed_f03_cmd_real_2_1(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_2_1

    subroutine plumed_f03_gcmd_real_2_1(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:)
      integer(kind=C_SIZE_T) :: pass_shape(2) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long_double(trim(key)//plumed_eos,val(1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_2_1

    subroutine plumed_f03_cmd_real_2_2(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_2_2

    subroutine plumed_f03_gcmd_real_2_2(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:)
      integer(kind=C_SIZE_T) :: pass_shape(3) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long_double(trim(key)//plumed_eos,val(1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_2_2

    subroutine plumed_f03_cmd_real_2_3(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_2_3

    subroutine plumed_f03_gcmd_real_2_3(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(4) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long_double(trim(key)//plumed_eos,val(1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_2_3

    subroutine plumed_f03_cmd_real_2_4(P,KEY,VAL)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_cmd_safe_long_double(P%ptr,trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_cmd_real_2_4

    subroutine plumed_f03_gcmd_real_2_4(KEY,VAL)
      character(len=*),  intent(in)    :: KEY
      real(KIND=c_long_double)                      :: VAL(:,:,:,:)
      integer(kind=C_SIZE_T) :: pass_shape(5) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      pass_shape=[ shape(val), (/0/)]
      pass_nelem=0
      call plumed_f_gcmd_safe_long_double(trim(key)//plumed_eos,val(1,1,1,1),pass_shape(1))
    end subroutine plumed_f03_gcmd_real_2_4

end module plumed_module
