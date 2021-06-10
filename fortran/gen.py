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
  shapestr1.append(s)

interface="""
  interface plumed_f03_cmd
    module procedure plumed_f03_cmd_null
    module procedure plumed_f03_cmd_char
"""

for t in types:
  ik=0
  for k in kinds[t]:
    for r in range(maxrank+1):
      interface+="    module procedure plumed_f03_cmd_{}_{}_{}\n".format(t,ik,r)
    ik+=1
interface+="  end interface plumed_f03_cmd"

interface+="""
  interface plumed_f03_gcmd
    module procedure plumed_f03_gcmd_null
    module procedure plumed_f03_gcmd_char
"""
for t in types:
  ik=0
  for k in kinds[t]:
    for r in range(maxrank+1):
      interface+="    module procedure plumed_f03_gcmd_{}_{}_{}\n".format(t,ik,r)
    ik+=1
interface+="  end interface plumed_f03_gcmd"

implementation=""

for t in types:
  ik=0
  for k in kinds[t]:
    for r in range(maxrank+1):
      compute_shape="pass_shape=(/1,0/)"
      cname=k.replace("c_","")
      if r>0:
        compute_shape="pass_shape=[ shape(val), (/0/)]"
      implementation+="""
    subroutine plumed_f03_cmd_{}_{}_{}(P,KEY,VAL,ERROR)
      type(plumed)                     :: P
      character(len=*),  intent(in)    :: KEY
      {}(KIND={})                      :: VAL{}
      type(plumed_error), optional     :: ERROR
      integer(kind=C_SIZE_T) :: pass_shape({}) ! fix this
      integer(kind=C_SIZE_T) :: pass_nelem
      integer(kind=C_LONG)   :: flags
      {}
      pass_nelem=0
      if(present(ERROR)) then
        ERROR%code=0
        ERROR%what=""
        call plumed_f_cmd_safe_nothrow_{}(P%ptr,trim(key)//plumed_eos,val{},pass_shape(1),error,error_handler)
      else
        call plumed_f_cmd_safe_{}(P%ptr,trim(key)//plumed_eos,val{},pass_shape(1))
      endif
    end subroutine plumed_f03_cmd_{}_{}_{}
""".format(t,ik,r,t,k,shapestr[r],max(r+1,2),compute_shape,cname,shapestr1[r],cname,shapestr1[r],t,ik,r)
      implementation+="""
    subroutine plumed_f03_gcmd_{}_{}_{}(KEY,VAL,ERROR)
      character(len=*),  intent(in)    :: KEY
      {}(KIND={})                      :: VAL{}
      type(plumed_error), optional     :: ERROR
      call plumed_f03_cmd_{}_{}_{}(plumed_f03_global(),KEY,VAL,ERROR)
    end subroutine plumed_f03_gcmd_{}_{}_{}
""".format(t,ik,r,t,k,shapestr[r],t,ik,r,t,ik,r)
    ik+=1

with open("plumed.f03.template") as f:
  for l in f:
    if l=="__INTERFACE__\n":
      print(interface)
    elif l=="__IMPLEMENTATION__\n":
      print(implementation)
    else:
      print(l,end="")

