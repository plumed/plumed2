! vim:ft=fortran
module plumed_module
  use iso_c_binding
  implicit none

  private
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

  interface
    subroutine plumed_f_create(c)
      import
      character(kind=c_char,len=32), intent(out) :: c
    end subroutine plumed_f_create
  end interface

  interface
    subroutine plumed_f_create_dlopen(path,c)
      import
      character(kind=c_char),        intent(in)  :: path(*)
      character(kind=c_char,len=32), intent(out) :: c
    end subroutine plumed_f_create_dlopen
  end interface

  interface
    subroutine plumed_f_create_reference(r,c)
      import
      character(kind=c_char,len=32), intent(in)  :: r
      character(kind=c_char,len=32), intent(out) :: c
    end subroutine plumed_f_create_reference
  end interface

  interface
    subroutine plumed_f_create_invalid(c)
      import
      character(kind=c_char,len=32), intent(out) :: c
    end subroutine plumed_f_create_invalid
  end interface

  interface
    subroutine plumed_f_finalize(c)
      import
      character(kind=c_char,len=32), intent(in) :: c
    end subroutine plumed_f_finalize
  end interface

  interface
    subroutine plumed_f_installed(i)
      import
      integer(kind=c_int), intent(out) :: i
    end subroutine plumed_f_installed
  end interface

  interface
    subroutine plumed_f_valid(c,i)
      import
      character(kind=c_char,len=32), intent(in) :: c
      integer(kind=c_int),           intent(out) :: i
    end subroutine plumed_f_valid
  end interface

  interface
    subroutine plumed_f_use_count(c,i)
      import
      character(kind=c_char,len=32), intent(in) :: c
      integer(kind=c_int),           intent(out) :: i
    end subroutine plumed_f_use_count
  end interface

  interface
    subroutine plumed_f_global(c)
      import
      character(kind=c_char,len=32), intent(out) :: c
    end subroutine plumed_f_global
  end interface 

  interface
    subroutine plumed_f_ginitialized(i)
      import
      integer(kind=c_int), intent(out) :: i
    end subroutine plumed_f_ginitialized
  end interface

  interface
    subroutine plumed_f_gcreate()
    end subroutine plumed_f_gcreate
  end interface
  
  interface
    subroutine plumed_f_gfinalize()
    end subroutine plumed_f_gfinalize
  end interface

  interface
    subroutine plumed_f_gvalid(i)
      import
      integer(kind=c_int),    intent(out) :: i
    end subroutine plumed_f_gvalid
  end interface


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
    module procedure plumed_f_cmd_real_2_0
    module procedure plumed_f_cmd_real_2_1
    module procedure plumed_f_cmd_real_2_2
    module procedure plumed_f_cmd_real_2_3
    module procedure plumed_f_cmd_real_2_4
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
    module procedure plumed_f_gcmd_real_2_0
    module procedure plumed_f_gcmd_real_2_1
    module procedure plumed_f_gcmd_real_2_2
    module procedure plumed_f_gcmd_real_2_3
    module procedure plumed_f_gcmd_real_2_4
  end interface plumed_f_gcmd

  contains

     subroutine plumed_f_cmd_char(p,key,val)
       character(len=32), intent(in)    :: p
       character(len=*),  intent(in)    :: key
       character(len=*)                 :: val
       integer(kind=c_size_t) :: pass_shape(1)
       pass_shape=(/0/)
       call plumed_f_cmd_safe_char(p,key,val,pass_shape)
     end subroutine plumed_f_cmd_char

     subroutine plumed_f_gcmd_char(key,val)
       character(len=*),  intent(in)    :: key
       character(len=*)                 :: val
       integer(kind=c_size_t) :: pass_shape(1)
       pass_shape=(/0/)
       call plumed_f_gcmd_safe_char(key,val,pass_shape)
     end subroutine plumed_f_gcmd_char


    subroutine plumed_f_cmd_integer_0_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_int)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_int(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_0_0

    subroutine plumed_f_gcmd_integer_0_0(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_int)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_int(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_0_0

    subroutine plumed_f_cmd_integer_0_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_int)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_int(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_0_1

    subroutine plumed_f_gcmd_integer_0_1(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_int)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_int(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_0_1

    subroutine plumed_f_cmd_integer_0_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_int)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_int(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_0_2

    subroutine plumed_f_gcmd_integer_0_2(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_int)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_int(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_0_2

    subroutine plumed_f_cmd_integer_0_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_int)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_int(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_0_3

    subroutine plumed_f_gcmd_integer_0_3(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_int)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_int(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_0_3

    subroutine plumed_f_cmd_integer_0_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_int)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_int(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_0_4

    subroutine plumed_f_gcmd_integer_0_4(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_int)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_int(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_0_4

    subroutine plumed_f_cmd_integer_1_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_short)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_short(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_1_0

    subroutine plumed_f_gcmd_integer_1_0(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_short)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_short(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_1_0

    subroutine plumed_f_cmd_integer_1_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_short)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_short(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_1_1

    subroutine plumed_f_gcmd_integer_1_1(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_short)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_short(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_1_1

    subroutine plumed_f_cmd_integer_1_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_short)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_short(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_1_2

    subroutine plumed_f_gcmd_integer_1_2(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_short)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_short(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_1_2

    subroutine plumed_f_cmd_integer_1_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_short)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_short(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_1_3

    subroutine plumed_f_gcmd_integer_1_3(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_short)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_short(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_1_3

    subroutine plumed_f_cmd_integer_1_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_short)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_short(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_1_4

    subroutine plumed_f_gcmd_integer_1_4(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_short)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_short(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_1_4

    subroutine plumed_f_cmd_integer_2_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_long)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_long(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_2_0

    subroutine plumed_f_gcmd_integer_2_0(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_long)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_long(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_2_0

    subroutine plumed_f_cmd_integer_2_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_long)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_2_1

    subroutine plumed_f_gcmd_integer_2_1(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_long)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_2_1

    subroutine plumed_f_cmd_integer_2_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_long)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_2_2

    subroutine plumed_f_gcmd_integer_2_2(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_long)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_2_2

    subroutine plumed_f_cmd_integer_2_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_long)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_2_3

    subroutine plumed_f_gcmd_integer_2_3(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_long)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_2_3

    subroutine plumed_f_cmd_integer_2_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      integer(KIND=c_long)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_integer_2_4

    subroutine plumed_f_gcmd_integer_2_4(key,val)
      character(len=*),  intent(in)    :: key
      integer(kind=c_long)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long(key,val,pass_shape)
    end subroutine plumed_f_gcmd_integer_2_4

    subroutine plumed_f_cmd_real_0_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_float)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_float(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_0_0

    subroutine plumed_f_gcmd_real_0_0(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_float)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_float(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_0_0

    subroutine plumed_f_cmd_real_0_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_float)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_float(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_0_1

    subroutine plumed_f_gcmd_real_0_1(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_float)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_float(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_0_1

    subroutine plumed_f_cmd_real_0_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_float)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_float(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_0_2

    subroutine plumed_f_gcmd_real_0_2(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_float)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_float(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_0_2

    subroutine plumed_f_cmd_real_0_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_float)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_float(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_0_3

    subroutine plumed_f_gcmd_real_0_3(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_float)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_float(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_0_3

    subroutine plumed_f_cmd_real_0_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_float)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_float(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_0_4

    subroutine plumed_f_gcmd_real_0_4(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_float)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_float(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_0_4

    subroutine plumed_f_cmd_real_1_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_double)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_1_0

    subroutine plumed_f_gcmd_real_1_0(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_double)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_1_0

    subroutine plumed_f_cmd_real_1_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_double)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_1_1

    subroutine plumed_f_gcmd_real_1_1(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_double)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_1_1

    subroutine plumed_f_cmd_real_1_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_double)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_1_2

    subroutine plumed_f_gcmd_real_1_2(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_double)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_1_2

    subroutine plumed_f_cmd_real_1_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_double)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_1_3

    subroutine plumed_f_gcmd_real_1_3(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_double)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_1_3

    subroutine plumed_f_cmd_real_1_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_double)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_1_4

    subroutine plumed_f_gcmd_real_1_4(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_double)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_1_4

    subroutine plumed_f_cmd_real_2_0(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_long_double)                      :: val
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_cmd_safe_long_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_2_0

    subroutine plumed_f_gcmd_real_2_0(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_long_double)                      :: val
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=(/1,0/)
      call plumed_f_gcmd_safe_long_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_2_0

    subroutine plumed_f_cmd_real_2_1(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_long_double)                      :: val(:)
      integer(kind=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_2_1

    subroutine plumed_f_gcmd_real_2_1(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_long_double)                      :: val(:)
      integer(kIND=c_size_t) :: pass_shape(2) ! fix this
      pass_shape=[(/size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_2_1

    subroutine plumed_f_cmd_real_2_2(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_long_double)                      :: val(:,:)
      integer(kind=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_2_2

    subroutine plumed_f_gcmd_real_2_2(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_long_double)                      :: val(:,:)
      integer(kIND=c_size_t) :: pass_shape(3) ! fix this
      pass_shape=[(/size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_2_2

    subroutine plumed_f_cmd_real_2_3(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_long_double)                      :: val(:,:,:)
      integer(kind=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_2_3

    subroutine plumed_f_gcmd_real_2_3(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_long_double)                      :: val(:,:,:)
      integer(kIND=c_size_t) :: pass_shape(4) ! fix this
      pass_shape=[(/size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_2_3

    subroutine plumed_f_cmd_real_2_4(p,key,val)
      character(len=32), intent(in)    :: p
      character(len=*),  intent(in)    :: key
      real(KIND=c_long_double)                      :: val(:,:,:,:)
      integer(kind=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_cmd_safe_long_double(p,key,val,pass_shape)
    end subroutine plumed_f_cmd_real_2_4

    subroutine plumed_f_gcmd_real_2_4(key,val)
      character(len=*),  intent(in)    :: key
      real(kind=c_long_double)                      :: val(:,:,:,:)
      integer(kIND=c_size_t) :: pass_shape(5) ! fix this
      pass_shape=[(/size(val,4),size(val,3),size(val,2),size(val,1)/),(/0/)]
      call plumed_f_gcmd_safe_long_double(key,val,pass_shape)
    end subroutine plumed_f_gcmd_real_2_4

end module plumed_module

