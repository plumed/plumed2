# usage: python3 gen.py < plumed.f90.template > plumed.f90
#types=["char","int","real"]
maxrank=4
kinds={
  "integer":["c_int","c_short","c_long"],
  "real":["c_float","c_double","c_long_double"]
}

types=kinds.keys()

shapestr=[""]
for i in range(maxrank):
  s="(:"
  s+=",:"*i
  s+=")"
  shapestr.append(s)

shapestr1=[""]
for i in range(maxrank):
  s="(1"
  s+=",1"*i
  s+=")"
#  shapestr1.append(s)
  shapestr1.append("")

interface="""
  interface plumed_f_cmd
    module procedure plumed_f_cmd_char
"""

for t in types:
  ik=0
  for k in kinds[t]:
    for r in range(maxrank+1):
      interface+="    module procedure plumed_f_cmd_{}_{}_{}\n".format(t,ik,r)
    ik+=1
interface+="  end interface plumed_f_cmd"

interface+="""
  interface plumed_f_gcmd
    module procedure plumed_f_gcmd_char
"""
for t in types:
  ik=0
  for k in kinds[t]:
    for r in range(maxrank+1):
      interface+="    module procedure plumed_f_gcmd_{}_{}_{}\n".format(t,ik,r)
    ik+=1
interface+="  end interface plumed_f_gcmd"

implementation=""

for t in types:
  ik=0
  for k in kinds[t]:
    cname=k.replace("c_","")
    interface+="""
  interface
    subroutine plumed_f_cmd_safe_{}_scalar(p,key,val,pass_shape) bind(C)
      import
      character(kind=c_char), intent(in)    :: p(32)
      character(kind=c_char), intent(in)    :: key(*)
      {}(kind={})                           :: val
      integer(kind=c_size_t) :: pass_shape(*)
    end subroutine plumed_f_cmd_safe_{}_scalar
  end interface
  interface
    subroutine plumed_f_cmd_safe_{}(p,key,val,pass_shape) bind(C)
      import
      character(kind=c_char), intent(in)    :: p(32)
      character(kind=c_char), intent(in)    :: key(*)
      {}(kind={})                           :: val(*)
      integer(kind=c_size_t) :: pass_shape(*)
    end subroutine plumed_f_cmd_safe_{}
  end interface
  interface
    subroutine plumed_f_gcmd_safe_{}_scalar(key,val,pass_shape) bind(C)
      import
      character(kind=c_char), intent(in)    :: key(*)
      {}(kind={})                           :: val
      integer(kind=c_size_t) :: pass_shape(*)
    end subroutine plumed_f_gcmd_safe_{}_scalar
  end interface
  interface
    subroutine plumed_f_gcmd_safe_{}(key,val,pass_shape) bind(C)
      import
      character(kind=c_char), intent(in)    :: key(*)
      {}(kind={})                           :: val(*)
      integer(kind=c_size_t) :: pass_shape(*)
    end subroutine plumed_f_gcmd_safe_{}
  end interface
""".format(cname,t,k,cname,cname,t,k,cname,cname,t,k,cname,cname,t,k,cname)
    for r in range(maxrank+1):
      compute_shape="pass_shape=(/1,0/)"
      cname=k.replace("c_","")
      if r==0:
         cname+="_scalar"
      if r>0:
        compute_shape="pass_shape=[(/"
        for i in range(r):
           if i>0:
               compute_shape+=","
           compute_shape+="size(val,{})".format(r-i)
        compute_shape+="/),(/0/)]"
      implementation+="""
    subroutine plumed_f_cmd_{}_{}_{}(p,key,val)
      character(kind=c_char,len=32), intent(in)    :: p
      character(kind=c_char,len=*),  intent(in)    :: key
      {}(KIND={})                      :: val{}
      integer(kind=c_size_t) :: pass_shape({})
      {}
      call plumed_f_cmd_safe_{}(p,key,val{},pass_shape)
    end subroutine plumed_f_cmd_{}_{}_{}
""".format(t,ik,r,t,k,shapestr[r],max(r+1,2),compute_shape,cname,shapestr1[r],t,ik,r)
      implementation+="""
    subroutine plumed_f_gcmd_{}_{}_{}(key,val)
      character(kind=c_char,len=*),  intent(in)    :: key
      {}(kind={})                      :: val{}
      integer(kind=c_size_t) :: pass_shape({})
      {}
      call plumed_f_gcmd_safe_{}(key,val{},pass_shape)
    end subroutine plumed_f_gcmd_{}_{}_{}
""".format(t,ik,r,t,k,shapestr[r],max(r+1,2),compute_shape,cname,shapestr1[r],t,ik,r)
    ik+=1

with open("plumed.f90.template") as f:
  for l in f:
    if l=="__INTERFACE__\n":
      print(interface)
    elif l=="__IMPLEMENTATION__\n":
      print(implementation)
    else:
      print(l,end="")

