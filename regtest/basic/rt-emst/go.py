import random
def print_pdb(coords,file):
    i=0
    for coord in coords:
      i+=1
      print("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(
        "ATOM",i,"type","","","",1,"",coord[0],coord[1],coord[2],1.0,1.0,"",""),file=f)

def print_xyz(box,coords,file):
    print(len(coords),file=f)
    print(box[0],box[1],box[2],file=f)
    for coord in coords:
      print("type ",coord[0],coord[1],coord[2],file=f)

coords=[]


box=[11,11,11]

for i in range(-20,21):
  if i%2==0: j=i
  else: j=-i
  coords.append([j,0.0,0.0])
for i in range(-20,21):
  if i%2==0: j=i
  else: j=-i
  coords.append([0.0,j,0.0])
for i in range(-20,21):
  if i%2==0: j=i
  else: j=-i
  coords.append([0.0,0.0,j])

with open("ref.pdb","w") as f:
  print_pdb(coords,f)


with open("traj.xyz","w") as f:
  print_xyz(box,coords,f)
  for i in range(len(coords)):
    coords[i][0]+=box[0]*(i%3)
  print_xyz(box,coords,f)


