
import unittest
import plumed

ib=plumed.InputBuilder()

try:
  import MDAnalysis
  _HAS_MDANALYSIS=True
except:
  _HAS_MDANALYSIS=False

class Test(unittest.TestCase):

 def check(self,s1,s2):
  self.assertEqual(s1,s2)

 def checkfiles(self,f1,f2):
  import filecmp
  import difflib
  import sys
  if not filecmp.cmp(f1,f2):
     s1=[]
     with open(f1,"r") as file1:
       for l in file1:
         s1.append(l)
     s2=[]
     with open(f2,"r") as file2:
       for l in file2:
         s2.append(l)
     for line in difflib.context_diff(s1, s2, fromfile=f1, tofile=f2):
         sys.stdout.write(line)
     self.assertTrue(False)

 def test1(self):
   self.check(ib.TORSION("phi",ATOMS="5,7,9,15") , 'phi: TORSION ATOMS=5,7,9,15\n')

 def test1b(self):
   self.check(ib.TORSION__("phi",ATOMS="5,7,9,15") , 'phi: TORSION ATOMS=5,7,9,15\n')

 def test2(self):
   self.check(ib.TORSION("phi",ATOMS="5 7 9 15") , 'phi: TORSION ATOMS={5 7 9 15}\n')

 def test3(self):
   self.check(ib.COORDINATION(GROUPA="1-10",GROUPB="11-20",SWITCH="RATIONAL NN=6 R_0=1") , 'COORDINATION GROUPA=1-10 GROUPB=11-20 SWITCH={RATIONAL NN=6 R_0=1}\n')

 def test4(self):
   self.check(ib.COORDINATION(GROUPA="""
  1 2 3 4 5 6 7 8 9 10
  """,GROUPB="""
  11 12 13 14 15 16 17 18 19 20
  """,SWITCH="""
  RATIONAL NN=6 R_0=1
  """) , 'COORDINATION GROUPA={   1 2 3 4 5 6 7 8 9 10   } GROUPB={   11 12 13 14 15 16 17 18 19 20   } SWITCH={   RATIONAL NN=6 R_0=1   }\n')

 def test5(self):
   self.check(ib.DISTANCE("d",ATOMS="11 21",NOPBC=True) , 'd: DISTANCE ATOMS={11 21} NOPBC\n')

 def test6(self):
   self.check(ib.DISTANCE("d",ATOMS="11 21",NOPBC=False),'d: DISTANCE ATOMS={11 21}\n')

 def test7(self):
   self.check(ib.METAD(ARG="phi psi",PACE=500,HEIGHT=1.2,SIGMA="0.35 0.35",FILE="HILLS"),'METAD ARG={phi psi} FILE=HILLS HEIGHT=1.2 PACE=500 SIGMA={0.35 0.35}\n')

 def test8(self):
   self.check(ib.GROUP("g",ATOMS=range(1,101)),'g: GROUP ATOMS={1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100}\n')

 def test9(self):
  self.check(ib.METAD(ARG=("phi","psi"),PACE=500,HEIGHT=1.2,SIGMA=(0.35,"pi/8"),FILE="HILLS"),'METAD ARG={phi psi} FILE=HILLS HEIGHT=1.2 PACE=500 SIGMA={0.35 pi/8}\n')

 def test10(self):
  self.check(ib.MOVINGRESTRAINT(ARG="d1",KAPPA0=0,KAPPA1=10.0,AT0=20,AT1=20,STEP0=1,STEP1=10000),'MOVINGRESTRAINT ARG=d1 AT0=20 AT1=20 KAPPA0=0 KAPPA1=10.0 STEP0=1 STEP1=10000\n')

 def test11(self):
  self.check(
    ib.MOVINGRESTRAINT(ARG="d1",KAPPA=ib.numbered([0,10.0]),AT=ib.numbered([20,20]),STEP=ib.numbered([1,10000]))
  ,
    'MOVINGRESTRAINT ARG=d1 AT0=20 AT1=20 KAPPA0=0 KAPPA1=10.0 STEP0=1 STEP1=10000\n'
  )

 def test12(self):
  self.check(
    ib.MOVINGRESTRAINT(ARG="d1",KAPPA=ib.numbered([100]),AT=ib.numbered({0:0.0,2:10.0}),STEP=ib.numbered((0,10,20,30)))
  ,
    'MOVINGRESTRAINT ARG=d1 AT0=0.0 AT2=10.0 KAPPA0=100 STEP0=0 STEP1=10 STEP2=20 STEP3=30\n'
  )

 def test13(self):
  self.check(
     ib.MOVINGRESTRAINT(ARG="d1,d2",
                     KAPPA=ib.numbered([11]),
                     AT=ib.numbered(((0.0,1.0),(2.0,3.0))),
                     STEP=ib.numbered((0,100)))
  ,
    'MOVINGRESTRAINT ARG=d1,d2 AT0={0.0 1.0} AT1={2.0 3.0} KAPPA0=11 STEP0=0 STEP1=100\n'
  )

 def test14(self):
  self.check(
     ib.MOVINGRESTRAINT(ARG="d1,d2",
                     KAPPA=ib.numbered([100]),AT=ib.numbered(([0.0,"pi"],[2.0,"pi"])),
                     STEP=ib.numbered((0,100)))
  ,
    'MOVINGRESTRAINT ARG=d1,d2 AT0={0.0 pi} AT1={2.0 pi} KAPPA0=100 STEP0=0 STEP1=100\n'
  )

 def test15(self):
  self.check(
    ib.RESTRAINT(ARG="d1",KAPPA=10,AT=ib.replicas((0.0,1.0,2.0,3.0)))
  ,
    'RESTRAINT ARG=d1 AT=@replicas:{0.0 1.0 2.0 3.0} KAPPA=10\n'
  )

 def test16(self):
  try:
    import numpy
  except:
    print("This test requires numpy module installed.")
  self.check(
     ib.RESTRAINT(ARG="d1",KAPPA=10,AT=ib.replicas(numpy.linspace(3.0,5.0,17)))
  ,
     'RESTRAINT ARG=d1 AT=@replicas:{3.0 3.125 3.25 3.375 3.5 3.625 3.75 3.875 4.0 4.125 4.25 4.375 4.5 4.625 4.75 4.875 5.0} KAPPA=10\n'
  )

 def test17(self):
  self.check(
    ib.RESTRAINT(ARG="d1,d2",KAPPA=(10,10),AT=ib.replicas(([0.0,1.0],[10.0,11.0])))
  ,
     'RESTRAINT ARG=d1,d2 AT=@replicas:{{0.0 1.0} {10.0 11.0}} KAPPA={10 10}\n'
  )

 def test18(self):
  ib1=plumed.InputBuilder(comma_separator=True)
  self.check(
    ib1.GROUP("g1",ATOMS=[1,2,3,4,5,6,7,8,9,10])
  ,
    'g1: GROUP ATOMS=1,2,3,4,5,6,7,8,9,10\n'
  )

 def test19(self):
  self.check(ib.at.phi(6),'@phi-6')

 def test20(self):
  self.check(ib.at.phi(range(1,11),"A"),'@phi-A1 @phi-A2 @phi-A3 @phi-A4 @phi-A5 @phi-A6 @phi-A7 @phi-A8 @phi-A9 @phi-A10')

 def test20(self):
  self.check(ib.at.phi(range(1,11),["A","B"]),'@phi-A1 @phi-A2 @phi-A3 @phi-A4 @phi-A5 @phi-A6 @phi-A7 @phi-A8 @phi-A9 @phi-A10 @phi-B1 @phi-B2 @phi-B3 @phi-B4 @phi-B5 @phi-B6 @phi-B7 @phi-B8 @phi-B9 @phi-B10')

 def test21(self):
  self.check(ib.at.phi(4,[1,2]),'@phi-1_4 @phi-2_4')

 def test22(self):
  self.check(ib.at("OW",range(20,40)),'@OW-20 @OW-21 @OW-22 @OW-23 @OW-24 @OW-25 @OW-26 @OW-27 @OW-28 @OW-29 @OW-30 @OW-31 @OW-32 @OW-33 @OW-34 @OW-35 @OW-36 @OW-37 @OW-38 @OW-39')

 def test23(self):
  self.check(ib.at.mdatoms,'@mdatoms')

 def test24(self):
  self.check(ib.RESTRAINT(ARG="d1,d2",verbatim="AT={10 20} KAPPA={5 6}"),'RESTRAINT ARG=d1,d2 AT={10 20} KAPPA={5 6}\n')

 def test25(self):
  self.check(ib.verbatim("# here is a comment"),'# here is a comment\n')

 if _HAS_MDANALYSIS:
  def test_mdanalysis(self):
    u=MDAnalysis.Universe("test/ref.pdb")
    self.check(
       ib.GROUP(ATOMS=u.select_atoms("name C2 C4 C6")),
       'GROUP ATOMS={16 20 25 50 54 59 84 88 93 118 122 127 147 151 156 178 182 187 209 213 218 240 244 249}\n'
    )
    self.check(
       ib.GROUP(ATOMS=u.select_atoms("name C2","name C4","name C6")),
       'GROUP ATOMS={20 54 88 122 156 187 218 249 25 59 93 127 151 182 213 244 16 50 84 118 147 178 209 240}\n'
    )

if __name__ == "__main__":
    unittest.main()

