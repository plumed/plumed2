#!/usr/bin/env python


###################################
## 1 # OPTIONS AND DOCUMENTATION ##  -> @DOC <-
###################################
from __future__ import division
import sys,logging
logging.basicConfig(format='%(levelname)-7s    %(message)s',level=9)

    
# This is a simple and versatily option class that allows easy
# definition and parsing of options. 
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self): 
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]
    
options = [
    ("-f",        Option(str,1,None, "PDB Input file.\n\t\tDifferent chains (or chain breaks) should be separated by TER or named with different chain-ID.\n\t\tOnly one model is expected. ENDMDL is accepted only at the end of the file.")),
    ("-dat",      Option(str,1,None,"Experimental Data Input file.\n\t\tEach line should contain momentum transfer (in inverse Ang) and non-zero intensity, separated by blanks or comma.\n\t\tIf momentum transfer is in inverse nm, use the option \"-unit nm\"")),
    ("-unit",     Option(str,1,"Ang","Unit of the momentum transfer [Ang/nm]. Default: Ang")),
    ("-nq",       Option(str,1,15,"Number of scattering points to evaluate. Default:15")),
    ("-warn",     Option(str,1,"NO","Accept warnings and do not exit.")),
    ("-h",        Option(bool,                     0,    False, "Display this help."))
  ]


def option_parser(args,options):

    # Check whether there is a request for help
    if '-h' in args or '--help' in args:
        help()


    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])
    options['Arguments']           = args[:]
    while args:
        ar = args.pop(0)
        options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])
    
    return options 


def help():
    """Print help text and list of options and end the program."""
    import sys
    print("\n")
    for item in options:
        if type(item) == str:
            print(item)
    for item in options:
        if type(item) != str:
            print("%10s  %s" % (item[0], item[1].description))
    print
    sys.exit()


#################################################
## 2 # HELPER FUNCTIONS, CLASSES AND SHORTCUTS ##  -> @FUNC <-
#################################################

import math

# Split a string                                                              
def spl(x):                                                                   
    return x.split()                                                          

# Split each argument in a list                                               
def nsplit(*x):                                                               
    return [i.split() for i in x]                                             

# Make a dictionary from two lists                                            
def hash(x,y):                                                                
    return dict(zip(x,y))                                                     

                                            
def norm2(a):
    return sum([i*i for i in a])


def norm(a):
    return math.sqrt(norm2(a))


def distance2(a,b):
    return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


##########################
## 3 # FG -> CG MAPPING ##  -> @MAP <-
##########################


dnares3 = " DA DC DG DT" 
dnares1 = " dA dC dG dT"
rnares3 = "  A  C  G  U"
rnares1 = " rA rC rG rU" # 

# Amino acid nucleic acid codes:                                                                                 
# The naming (AA and '3') is not strictly correct when adding DNA/RNA, but we keep it like this for consistincy.
AA3     = spl("TRP TYR PHE HIS HIH ARG LYS CYS ASP GLU ILE LEU MET ASN PRO HYP GLN SER THR VAL ALA GLY"+dnares3+rnares3) #@#
AA1     = spl("  W   Y   F   H   H   R   K   C   D   E   I   L   M   N   P   O   Q   S   T   V   A   G"+dnares1+rnares1) #@#


# Dictionaries for conversion from one letter code to three letter code v.v.                         
AA123, AA321 = hash(AA1,AA3),hash(AA3,AA1)                                                           


# Residue classes:
protein = AA3[:-8]   # remove eight to get rid of DNA/RNA here.
water   = spl("HOH SOL TIP")
lipids  = spl("DPP DHP DLP DMP DSP POP DOP DAP DUP DPP DHP DLP DMP DSP PPC DSM DSD DSS")
nucleic = spl("DAD DCY DGU DTH ADE CYT GUA THY URA DA DC DG DT A C G U")


residueTypes = dict(
    [(i,"Protein") for i in protein ]+
    [(i,"Water")   for i in water   ]+
    [(i,"Lipid")   for i in lipids  ]+
    [(i,"Nucleic") for i in nucleic ]
    )

class CoarseGrained:
    # Class for mapping an atomistic residue list to a coarsegrained one
    # Should get an __init__ function taking a residuelist, atomlist, Pymol selection or ChemPy model
    # The result should be stored in a list-type attribute
    # The class should have pdbstr and grostr methods

    # Standard mapping groups
    bb        = "N CA C O H H1 H2 H3 O1 O2 OC1 OC2"                                                                    #@#  
    dna_bb = "P OP1 OP2 O5' O3'","C5' O4' C4'","C3' C2' C1'"
    rna_bb = "P OP1 OP2 O5' O3'","C5' O4' C4'","C3' C2' O2' C1'"


    # This is the mapping dictionary
    # For each residue it returns a list, each element of which
    # lists the atom names to be mapped to the corresponding bead.
    # The order should be the standard order of the coarse grained
    # beads for the residue. Only atom names matching with those 
    # present in the list of atoms for the residue will be used
    # to determine the bead position. This adds flexibility to the
    # approach, as a single definition can be used for different 
    # states of a residue (e.g., GLU/GLUH).
    # For convenience, the list can be specified as a set of strings,
    # converted into a list of lists by 'nsplit' defined above.
    mapping = {
        "ALA":  nsplit(bb + " CB"),
        "CYS":  nsplit(bb,"CB SG"),
        "ASP":  nsplit(bb,"CB CG OD1 OD2"),
        "GLU":  nsplit(bb,"CB CG CD OE1 OE2"),
        "PHE":  nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ HZ"),
        "GLY":  nsplit(bb),
        "HIS":  nsplit(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),
        "HIH":  nsplit(bb,"CB CG","CD2 HD2 NE2 HE2","ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  nsplit(bb,"CB CG1 CG2 CD CD1"),
        "LYS":  nsplit(bb,"CB CG CD","CE NZ HZ1 HZ2 HZ3"),
        "LEU":  nsplit(bb,"CB CG CD1 CD2"),
        "MET":  nsplit(bb,"CB CG SD CE"),
        "ASN":  nsplit(bb,"CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  nsplit(bb,"CB CG CD"),
        "HYP":  nsplit(bb,"CB CG CD OD"),
        "GLN":  nsplit(bb,"CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"), # 1HE2 2HE2 1HE1 1HE2
        "ARG":  nsplit(bb,"CB CG CD","NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),    #1HH1 1HH2 2HH1 2HH2
        "SER":  nsplit(bb,"CB OG HG"),
        "THR":  nsplit(bb,"CB OG1 HG1 CG2"),
        "VAL":  nsplit(bb,"CB CG1 CG2"),
        "TRP":  nsplit(bb,"CB CG CD2","CD1 HD1 NE1 HE1 CE2","CE3 HE3 CZ3 HZ3","CZ2 HZ2 CH2 HH2"),
        "TYR":  nsplit(bb,"CB CG CD1 HD1","CD2 HD2 CE2 HE2","CE1 HE1 CZ OH HH"),
        "DA": nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' C1'",
                          "N9 C4",
                          "C2 N3",
                          "C6 N6 N1",
                          "C8 N7 C5"),
        "DG": nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' C1'",
                          "N9 C4",
                          "C2 N2 N3",
                          "C6 O6 N1",
                          "C8 N7 C5"),
        "DC": nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' C1'",
                          "N1 C6",
                          "N3 C2 O2",
                          "C5 C4 N4"),
        "DT": nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' C1'",
                          "N1 C6",
                          "N3 C2 O2",
                          "C5 C4 O4 C7 C5M"),
         "A":  nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' O2' C1'",
                          "N9 C4",
                          "C2 N3",
                          "C6 N6 N1",
                          "C8 N7 C5"),
        "G":  nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' O2' C1'",
                          "N9 C4",
                          "C2 N2 N3",
                          "C6 O6 N1",
                          "C8 N7 C5"),
        "C":  nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' O2' C1'",
                          "N1 C6",
                          "N3 C2 O2",
                          "C5 C4 N4"),
        "U":  nsplit("P OP1 OP2 O5' O3' O1P O2P",
                          "C5' O4' C4'",
                          "C3' C2' O2' C1'",
                          "N1 C6",
                          "N3 C2 O2",
                          "C5 C4 O4 C7 C5M"),
        }



 

    # Generic names for side chain beads
    residue_bead_names = spl("BB SC1 SC2 SC3 SC4")
    # Generic names for DNA/RNA beads
    residue_bead_names_dna = spl("BB1 BB2 BB3 SC1 SC2 SC3 SC4")
    residue_bead_names_rna = spl("BB1 BB2 BB3 SC1 SC2 SC3 SC4")

    # This dictionary contains the bead names for all residues,
    # following the order in 'mapping'
    # Add default bead names for all amino acids
    names = {}
    names.update([(i,("BB","SC1","SC2","SC3","SC4")) for i in AA3])
    # Add the default bead names for all DNA nucleic acids
    names.update([(i,("BB1","BB2","BB3","SC1","SC2","SC3","SC4")) for i in nucleic])

    # Crude mass for weighted average. No consideration of united atoms.
    # This will probably give only minor deviations, while also giving less headache
    mass = {'H': 1,'C': 12,'N': 14,'O': 16,'S': 32,'P': 31,'M': 0}
    
# Determine average position for a set of weights and coordinates
# This is a rather specific function that requires a list of items
# [(m,(x,y,z),id),..] and returns the weighted average of the 
# coordinates and the list of ids mapped to this bead
def aver(b):
    mwx,ids,atom,aids = zip(*[((m*x,m*y,m*z),i,at,aid) for m,(x,y,z),i,at,aid in b])              # Weighted coordinates     
    tm  = sum(zip(*b)[0])                                                 # Sum of weights           
    return [sum(i)/tm for i in zip(*mwx)],ids,atom,aids                            # Centre of mass           

                         
# Return the CG beads for an atomistic residue, using the mapping specified above
# The residue 'r' is simply a list of atoms, and each atom is a list:
# [ name, resname, resid, chain, x, y, z ]
def mapCG(r):
    p = CoarseGrained.mapping[r[0][1]]     # Mapping for this residue 
    # Get the atom_name, mass, coordinates (x,y,z), atom id for all atoms in the residue
    a = [(i[0],CoarseGrained.mass.get(i[0][0],0),i[4:7],i[7]) for i in r]               
    # Store weight, coordinate and index for atoms that match a bead
    q = [[(m,coord,a.index((atom,m,coord,aid)),atom,aid) for atom,m,coord,aid in a if atom in i] for i in p]
    # Return 
    # pos: bead positions, index, atom name, atom id  for atoms in beads 
    return zip(*[aver(i) for i in q if len(i)>0]) 




#######################
## 4 # STRUCTURE I/O ##  -> @IO <-
#######################
import math,sys

#----+---------+
## A | PDB I/O |
#----+---------+

# Reformatting of lines in structure file                                     
pdbAtomLine = "ATOM  %5s %4s%4s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n"        

def pdbAtom(a):
    ##01234567890123456789012345678901234567890123456789012345678901234567890123456789
    ##ATOM   2155 HH11 ARG C 203     116.140  48.800   6.280  1.00  0.00
    if a.startswith("TER"):
        return 0
    # NOTE: The 27th field of an ATOM line in the PDB definition can contain an
    #       insertion code. We shift that 20 bits and add it to the residue number
    #       to ensure that residue numbers will be unique.
    ## ===> atom name,       res name,        res id,                        chain,
    return (a[12:16].strip(),a[17:20].strip(),int(a[22:26])+(ord(a[26])<<20),a[21],
    ##            x,              y,              z       
    float(a[30:38]),float(a[38:46]),float(a[46:54]),int(a[6:12]))

# Function for splitting a PDB file in chains, based
# on chain identifiers and TER statements
def pdbChains(pdbAtomList):
    chain = []
    for atom in pdbAtomList:
        if not atom: # Was a "TER" statement
            if chain:
                yield chain
            else:
                logging.info("Skipping empty chain definition")
            chain = [] 
            continue
        if not chain or chain[-1][3] == atom[3]:
            chain.append(atom)
        else:
            yield chain
            chain = [atom]
    if chain:
        yield chain


# Simple PDB iterator
def pdbFrameIterator(streamIterator):  
    title, atoms = [], []
    c = 0 
    for i in streamIterator:
        if i.startswith("ENDMDL"):
            yield "".join(title), atoms, c
            title, atoms = [], []         
        elif i.startswith("TITLE"):
            title.append(i)
        elif i.startswith("ATOM") or i.startswith("HETATM"):
            atoms.append(pdbAtom(i))
            c += 1
        elif i.startswith("TER"):
            atoms.append("")

    if atoms:
        yield "".join(title), atoms, c



#----+---------+
## B | GRO I/O |
#----+---------+

groline = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"                                    

def groAtom(a):
    # In PDB files, there might by an insertion code. To handle this, we internally add
    # constant to all resids. To be consistent, we have to do the same for gro files.
    # 32 equal ord(' '), eg an empty insertion code
    constant = 32<<20
    #012345678901234567890123456789012345678901234567890
    #    1PRN      N    1   4.168  11.132   5.291
    ## ===> atom name,        res name,          res id,    chain,
    return (a[10:15].strip(), a[5:10].strip(),   int(a[:5])+constant, " ", 
    ##               x,                 y,                 z,  atomid      
    10*float(a[20:28]),10*float(a[28:36]),10*float(a[36:44]), int(a[15:20]) )

# Simple GRO iterator
def groFrameIterator(streamIterator):
    while True:
        try:
            title = streamIterator.next()
        except StopIteration:
            break
        natoms = streamIterator.next().strip()
        if not natoms:
            break
        natoms = int(natoms)
        atoms  = [groAtom(streamIterator.next())  for i in range(natoms)] 
        yield title, atoms, natoms

#----+-------------+
## C | GENERAL I/O |
#----+-------------+

# *NOTE*: This should probably be a CheckableStream class that
# reads in lines until either of a set of specified conditions
# is met, then setting the type and from thereon functioning as
# a normal stream.
def streamTag(stream):
    # Tag the stream with the type of structure file
    # If necessary, open the stream, taking care of 
    # opening using gzip for gzipped files

    # First check whether we have have an open stream or a file
    # If it's a file, check whether it's zipped and open it
    if type(stream) == str:
        if stream.endswith("gz"):
            logging.info('Read input structure from zipped file.')
            s = gzip.open(stream)
        else:
            logging.info('Read input structure from file.')
            s = open(stream)
    else:
        logging.info('Read input structure from command-line')
        s = stream

    # Read a few lines, but save them
    x = [s.readline(), s.readline()]
    if x[-1].strip().isdigit():
        # Must be a GRO file
        logging.info("Input structure is a GRO file. Chains will be labeled consecutively.")
        yield "GRO"
    else:
        # Must be a PDB file then
        # Could wind further to see if we encounter an "ATOM" record
        logging.info("Input structure is a PDB file.")
        yield "PDB"
    
    # Hand over the lines that were stored
    for i in x:
        yield i

    # Now give the rest of the lines from the stream
    for i in s:
        yield i





# Converts an integer to a base36 string

digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
digits_lower = digits_upper.lower()
digits_upper_values = dict([pair for pair in zip(digits_upper, range(36))])
digits_lower_values = dict([pair for pair in zip(digits_lower, range(36))])

def encode_pure(digits, value):
  "encodes value using the given digits"
  assert value >= 0
  if (value == 0): return digits[0]
  n = len(digits)
  result = []
  while (value != 0):
    rest = value // n
    result.append(digits[value - rest * n])
    value = rest
  result.reverse()
  return "".join(result)

def hy36encode(width, value):
  "encodes value as base-10/upper-case base-36/lower-case base-36 hybrid"
  i = value
  if (i >= 1-10**(width-1)):
    if (i < 10**width):
      return ("%%%dd" % width) % i
    i -= 10**width
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_upper, i)
    i -= 26*36**(width-1)
    if (i < 26*36**(width-1)):
      i += 10*36**(width-1)
      return encode_pure(digits_lower, i)
  raise ValueError("value out of range.")



# if atm-number >= 100000, starts again from 0
def renumberPDB(filename,firstAtm):
    filepdb = open(filename, 'r')
    count=firstAtm
    for line in filepdb:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            yield line[:6] + str(count).rjust(5) + line[11:]
            count += 1
            if count == 100000:
                count = 0
        else:
            yield line

# if atm-number >= 100000, use hy36 encode
# delete ENDMDL
def renumberPDBh36(filename,firstAtm):
    filepdb = open(filename, 'r')
    count=firstAtm
    for line in filepdb:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if count >= 10000:
                newc=hy36encode(5,count)
            else:
                newc=count
            yield line[:6] + str(newc).rjust(5) + line[11:]
            count += 1
        else:
            yield line



#----+-----------------+
## D | STRUCTURE STUFF |
#----+-----------------+


# This list allows to retrieve atoms based on the name or the index
# If standard, dictionary type indexing is used, only exact matches are
# returned. Alternatively, partial matching can be achieved by setting
# a second 'True' argument. 
class Residue(list):
    def __getitem__(self,tag): 
        if type(tag) == int:
            # Call the parent class __getitem__
            return list.__getitem__(self,tag)
        if type(tag) == str:
            for i in self:
                if i[0] == tag:
                    return i
            return
        if tag[1]:
            return [i for i in self if tag[0] in i[0]] # Return partial matches
        else:
            return [i for i in self if i[0] == tag[0]] # Return exact matches only


def residues(atomList):
    residue = [atomList[0]]
    for atom in atomList[1:]:
        if (atom[1] == residue[-1][1] and # Residue name check
            atom[2] == residue[-1][2] and # Residue id check
            atom[3] == residue[-1][3]):   # Chain id check
            residue.append(atom)
        else:
            yield Residue(residue)
            residue = [atom]
    yield Residue(residue)


def residueDistance2(r1,r2):
    return min([distance2(i,j) for i in r1 for j in r2])


# CUT-OFF_CHANGE
# Increased the cut-off from 2.5 to 3.0 for DNA. Might have to be adjusted for other DNA structures.
# Not guaranteed to work for proteins, recommended to use the regular martinize for them. 
def breaks(residuelist,selection=("N","CA","C","P","C2'","C3'","O3'","C4'","C5'","O5'", "OP1", "OP2"),cutoff=3.0):
    # Extract backbone atoms coordinates
    bb = [[atom[4:] for atom in residue if atom[0] in selection] for residue in residuelist]
    # Needed to remove waters residues from mixed residues.
    bb = [res for res in bb if res != []]

    # We cannot rely on some standard order for the backbone atoms.
    # Therefore breaks are inferred from the minimal distance between
    # backbone atoms from adjacent residues.
    return [ i+1 for i in range(len(bb)-1) if residueDistance2(bb[i],bb[i+1]) > cutoff]


###################################################


## !! NOTE !! ##
## XXX The chain class needs to be simplified by extracting things to separate functions/classes
class Chain:
    # Attributes defining a chain
    # When copying a chain, or slicing, the attributes in this list have to
    # be handled accordingly.
    _attributes = ("residues","sequence","seq")

    def __init__(self,options,residuelist=[],name=None):
        self.residues   = residuelist
        self._atoms     = [atom[:3] for residue in residuelist for atom in residue]
        self.sequence   = [residue[0][1] for residue in residuelist]
        # *NOTE*: Check for unknown residues and remove them if requested
        #         before proceeding.
        self.seq        = "".join([AA321.get(i,"X") for i in self.sequence])
        self.mapping    = []
        self.options    = options

        # Unknown residues
        self.unknowns   = "X" in self.seq

        # Determine the type of chain
        self._type      = ""
        self.type()

        # Determine number of atoms
        self.natoms     = len(self._atoms) 

        # BREAKS: List of indices of residues where a new fragment starts
        # Only when polymeric (protein, DNA, RNA, ...)
        # For now, let's remove it for the Nucleic acids...
        self.breaks     = self.type() in ("Protein","Mixed") and breaks(self.residues) or []

        # Chain identifier; try to read from residue definition if no name is given
        self.id         = name or residuelist and residuelist[0][0][3] or ""

        # Container for coarse grained beads
        self._cg        = None

        
    def __len__(self):
        # Return the number of residues
        # DNA/RNA contain non-CAP d/r to indicate type. We remove those first.
        return len(''.join(i for i in self.seq if i.isupper()))

    def __add__(self,other):
        newchain = Chain(name=self.id+"+"+other.id)
        # Combine the chain items that can be simply added
        for attr in self._attributes:
            setattr(newchain, attr, getattr(self,attr) + getattr(other,attr))
        # Set chain items, shifting the residue numbers
        shift  = len(self)
        newchain.breaks     = self.breaks + [shift] + [i+shift for i in other.breaks]
        newchain.natoms     = len(newchain.atoms())
        # Return the merged chain
        return newchain

    def __eq__(self,other):
        return (self.seq        == other.seq    and 
                self.breaks     == other.breaks )
    
    def __ne__(self, other):
        return not self.__eq__(other)


    # Extract a residue by number or the list of residues of a given type
    # This facilitates selecting residues for links, like chain["CYS"]
    def __getitem__(self,other):
        if type(other) == str:
            if not other in self.sequence:
                return []
            return [i for i in self.residues if i[0][1] == other]
        elif type(other) == tuple:
            # This functionality is set up for links
            # between coarse grained beads. So these are
            # checked first,
            for i in self.cg():
                if other == i[:4]:
                    return i
            for i in self.atoms():
                if other[:3] == i[:3]:
                    return i
            return []
        elif type(other) == slice:
            # This implements the old __getslice__ method 
            i, j = other.start, other.stop
            newchain = Chain(self.options, name=self.id)
            # Extract the slices from all lists
            for attr in self._attributes:
                setattr(newchain, attr, getattr(self, attr)[i:j])
            # Breaks that fall within the start and end of this chain need to be passed on.
            # Residue numbering is increased by 20 bits!!
            ch_sta, ch_end      = newchain.residues[0][0][2], newchain.residues[-1][0][2]
            newchain.breaks     = [crack for crack in self.breaks if ch_sta < (crack << 20) < ch_end]
            newchain.natoms     = len(newchain.atoms())
            # newchain.type()
            # Return the chain slice
            return newchain
        return self.sequence[other]

    def _contains(self,atomlist,atom):
        atnm,resn,resi,chn = atom
        
        # If the chain does not match, bail out
        if chn != self.id:
            return False

        # Check if the whole tuple is in
        if atnm and resn and resi:
            return (atnm,resn,resi) in self.atoms()

        # Fetch atoms with matching residue id
        match = (not resi) and atomlist or [j for j in atomlist if j[2] == resi]
        if not match:
            return False

        # Select atoms with matching residue name
        match = (not resn) and match or [j for j in match if j[1] == resn]
        if not match:
            return False

        # Check whether the atom is given and listed
        if not atnm or [j for j in match if j[0] == atnm]:
            return True

        # It just is not in the list!
        return False

    def __contains__(self,other):
        return self._contains(self.atoms(),other) or self._contains(self.cg(),other)

    def __hash__(self):
        return id(self)

    def atoms(self):
        if not self._atoms:
            self._atoms = [atom[:3] for residue in self.residues for atom in residue]
        return self._atoms

    # Split a chain based on residue types; each subchain can have only one type
    def split(self):
        chains = []
        chainStart = 0
        for i in range(len(self.sequence)-1):
            if residueTypes.get(self.sequence[i],"Unknown") != residueTypes.get(self.sequence[i+1],"Unknown"):
                # Use the __getslice__ method to take a part of the chain.
                chains.append(self[chainStart:i+1])
                chainStart = i+1
        if chains:
            logging.debug('Splitting chain %s in %s chains'%(self.id,len(chains)+1))
        return chains + [self[chainStart:]]

    def getname(self,basename=None):
        name = []
        if basename:                      name.append(basename)
        if self.type() and not basename:  name.append(self.type())
        if type(self.id) == int:
            name.append(chr(64+self.id))
        elif self.id.strip():               
            name.append(str(self.id))
        return "_".join(name)

    def type(self,other=None):
        if other:
            self._type = other
        elif not self._type and len(self):
            # Determine the type of chain
            self._type     = set([residueTypes.get(i,"Unknown") for i in set(self.sequence)])
            self._type     = list(self._type)[0]
        return self._type


    # XXX The following (at least the greater part of it) should be made a separate function, put under "MAPPING"
    def cg(self,force=False,com=False,dna=False):
        # Generate the coarse grained structure
        # Set the b-factor field to something that reflects the secondary structure
        
        # If the coarse grained structure is set already, just return, 
        # unless regeneration is forced.
        if self._cg and not force:
            return self._cg
        self._cg = []
        atid     = 1
        bb       = [1]
        fail     = False
        previous = ''

        for residue,resname in zip(self.residues,self.sequence):
            # For DNA we need to get the O3' to the following residue when calculating COM
            # The force and com options ensure that this part does not affect itp generation or anything else
            if com:
                # Just an initialization, this should complain if it isn't updated in the loop
                store = 0
                for ind, i in enumerate(residue):
                    if i[0] == "O3'":
                        if previous != '':
                            residue[ind] = previous
                            previous = i
                        else:
                            store = ind
                            previous = i
                # We couldn't remove the O3' from the 5' end residue during the loop so we do it now
                if store > 0:
                    del residue[store]

            # Check if residues names has changed, for example because user has set residues interactively.
            residue = [(atom[0],resname)+atom[2:] for atom in residue]
            if residue[0][1] in ("SOL","HOH","TIP"):
                continue
            if not residue[0][1] in CoarseGrained.mapping.keys():
                logging.warning("Skipped unknown residue %s\n"%residue[0][1])
                continue

            # Get the mapping for this residue
            # CG.map returns bead coordinates and mapped atoms
            # This will fail if there are (too many) atoms missing, which is
            # only problematic if a mapped structure is written; the topology
            # is inferred from the sequence. So this is the best place to raise 
            # an error
            try:
                beads, ids, atm, aid = mapCG(residue)
                beads      = zip(CoarseGrained.names[residue[0][1]],beads,ids,atm,aid)

            except ValueError:
                logging.error("Too many atoms missing from residue %s %d(ch:%s):",residue[0][1],residue[0][2]>>20,residue[0][3])
                logging.error(repr([ i[0] for i in residue ]))
                fail = True

            for name,(x,y,z),ids,atm,aid in beads:                    
                # Add the bead with coordinates and secondary structure id to the list
                self._cg.append((name,residue[0][1][:3],residue[0][2],residue[0][3],x,y,z,atm,aid))

           
            # Increment the atom id; This pertains to the atoms that are included in the output.
            atid += len(residue)

            # Keep track of the numbers for CONECTing
            bb.append(bb[-1]+len(beads))

        if fail:
            logging.error("Unable to generate coarse grained structure due to missing atoms.")
            sys.exit(1)

        return self._cg

 

#############
## 8 # MAIN #  -> @MAIN <-
#############
import sys,math

def main(options):
    # Check whether to read from a gro/pdb file or from stdin
    # We use an iterator to wrap around the stream to allow
    # inferring the file type, without consuming lines already
    inStream = streamTag(options["-f"] and options["-f"].value or sys.stdin)
    # The streamTag iterator first yields the file type, which 
    # is used to specify the function for reading frames
    fileType = inStream.next()
    if fileType == "GRO":
        # We do not want gro files since we need to create a full PDB model AA/CG to be used as plumed input
        logging.error("Input file must be a PDB!")
        sys.exit(1)
        #frameIterator = groFrameIterator
        # if file is gro I expect atom numbers to be ordered starting from 1
    else:
        frameIterator = pdbFrameIterator
        # reorder atom numbers starting from 1
        inStream=renumberPDB(options["-f"].value,1)

    ## ITERATE OVER FRAMES IN STRUCTURE FILE ##
    # Now iterate over the frames in the stream
    # This should become a StructureFile class with a nice .next method
    model     = 1 # we expect only 1 model here
    
    for title,atoms,NATM in frameIterator(inStream):
        if model > 1:
            logging.error("Only one model (and no ENDMDL) should be present in the PDB file!")
            sys.exit(1)

        if fileType == "PDB":
            # The PDB file can have chains, in which case we list and process them specifically
            # TER statements are also interpreted as chain separators
            # A chain may have breaks in which case the breaking residues are flagged
            chains = [ Chain(options,[i for i in residues(chain)]) for chain in pdbChains(atoms) ]    
        else:
            # We do not want gro files since we need to create a full PDB model AA/CG to be used as plumed input
            logging.error("Input file must be a PDB!")
            sys.exit(1)
            # DELETE COMMENT IF YOU WANT TO WORK WITH GRO
            # The GRO file does not define chains. Here breaks in the backbone are
            # interpreted as chain separators. 
            # residuelist = [residue for residue in residues(atoms)]
            # The breaks are indices to residues
            # broken = breaks(residuelist)
            # Reorder, such that each chain is specified with (i,j,k)
            # where i and j are the start and end of the chain, and 
            # k is a chain identifier
            # chains = zip([0]+broken,broken+[len(residuelist)],range(len(broken)+1))
            # chains = [ Chain(options,residuelist[i:j],name=chr(65+k)) for i,j,k in chains ]
    
        # Check the chain identifiers
        if model == 1 and len(chains) != len(set([i.id for i in chains])):
            # Ending down here means that non-consecutive blocks of atoms in the 
            # PDB file have the same chain ID. The warning pertains to PDB files only, 
            # since chains from GRO files get a unique chain identifier assigned.
            logging.warning("Several chains have identical chain identifiers in the PDB file.")
   
        # Check if chains are of mixed type. If so, split them.
        # Note that in some cases HETATM residues are part of a 
        # chain. This will get problematic. But we cannot cover
        # all, probably.
        demixedChains = []
        for chain in chains:
            demixedChains.extend(chain.split())
        chains = demixedChains

        n = 1
        NATOMSTOT=0
        logging.info("Found %d chains:"%len(chains))
        for chain in chains:
            logging.info("  %2d:   %s (%s), %d atoms in %d residues."%(n,chain.id,chain._type,chain.natoms,len(chain)))
            NATOMSTOT += chain.natoms
            n += 1
    
        # Check all chains
        keep = []
        for chain in chains:
            if chain.type() == "Water":
                logging.info("Removing %d water molecules (chain %s)."%(len(chain),chain.id))
            elif chain.type() in ("Protein","Nucleic"):
                keep.append(chain)
            else:
                logging.info("Removing HETATM chain %s consisting of %d residues."%(chain.id,len(chain)))
        chains = keep


        if model == 1:
            order = []
            order.extend([j for j in range(len(chains)) if not j in order])

        # Get the total length of the sequence
        seqlength = sum([len(chain) for chain in chains])


        if NATOMSTOT == NATM:
            logging.info('Total size of the system: %s residues and %d atoms'% (seqlength,NATOMSTOT))
        else:
            logging.info('Total size of the system: %s residues'% (seqlength))
            if options["-warn"].value != "NO":
                logging.warning('Not sure about the total number of atoms in input file. Here using %d atoms. Double check the output files.' %NATOMSTOT )
            else:
                logging.warning('Not sure about the total number of atoms in input file. Exiting...')
                logging.warning('If you want to skip this warning and use %d atoms you can run the script with the option \"-warn 1\"' %NATOMSTOT)
                sys.exit(1)

        # Write the coarse grained structure if requested
        cgOutPDB = open("aacg_template.pdb","w")
        cgMAP = open("plumed_beads.dat","w")


        for i in list((renumberPDBh36(options["-f"].value,1))):
            if i[:6] != "ENDMDL":
                cgOutPDB.write(i)
        
        cgOutPDB.write("MODEL %8d\n"%model)
        cgOutPDB.write(title)
           
        atid = 1;
        write_start = 0
        stringbeads = 'martini: GROUP ATOMS='
        for i in order:
            ci = chains[i]
            coarseGrained = ci.cg(com=True)
            if coarseGrained:
                # For DNA we need to remove the first bead on the 5' end and shift the atids. 
                if ci.type() == 'Nucleic':
                    write_start = 1
                else:
                    write_start = 0

                # CG pdb and index

                lastr=coarseGrained[-1][2]
                insc = lastr>>20
                lastr -= insc<<20
                firstr=coarseGrained[write_start][2]
                insc = firstr>>20
                firstr -= insc<<20

                for name,resn,resi,chain,x,y,z,atm,aid in coarseGrained[write_start:]:
                    insc  = resi>>20
                    resi -= insc<<20
                    
                    if atid+NATOMSTOT >= 100000:
                        newid=hy36encode(5,atid+NATOMSTOT)
                    else:
                        newid=atid+NATOMSTOT

                    #write pdb
                    if  resi == firstr and name=="BB2":
                        cgOutPDB.write(pdbAtomLine%(newid,"TE5",resn[:3],"",resi,chr(insc),x,y,z,1,0))
                    elif resi == lastr and name=="BB3":
                        cgOutPDB.write(pdbAtomLine%(newid,"TE3",resn[:3],"",resi,chr(insc),x,y,z,1,0))
                    else:
                        cgOutPDB.write(pdbAtomLine%(newid,name,resn[:3],"",resi,chr(insc),x,y,z,1,0))

                    #cgOutIND.write
                    linemap = 'bead' + str(atid) + ': CENTER ATOMS='
                    for k in aid:
                        linemap+='%d,'%k
                    linemap = linemap[:-1] + ' WEIGHTS='
                    for k in atm:
                        linemap+='%s,'%CoarseGrained.mass.get(k[0][0],0)

                    linemap = linemap[:-1] + '\n'
                    cgMAP.write(linemap)
                        
                    stringbeads = stringbeads + 'bead' + str(atid)  + ','
                    atid += 1
                cgOutPDB.write("TER\n") 
                     
            else:
                logging.warning("No mapping for coarse graining chain %s (%s); chain is skipped."%(ci.id,ci.type()))
        cgMAP.write("\n")
        cgMAP.write(stringbeads[:-1])
        cgMAP.write("\n")
        cgOutPDB.write("ENDMDL\n")
        cgOutPDB.close(); cgMAP.close();
       
        model += 1

    if options["-dat"].value:
        
        datfile=open(options["-dat"].value,"r")
        q=[]; expdata=[];

        plumedmain = open("plumed.dat","w")
        plumedmain.write("MOLINFO STRUCTURE=aacg_template.pdb\n\n# BEADS DEFINITION\nINCLUDE FILE=plumed_beads.dat\n\n")
        plumedmain.write("# SAXS\nSAXS ...\n\n\tLABEL=saxsdata\n\tATOMS=martini\n\tMARTINI\n\n\t# You can use SCALEINT keyword to set appropriate scaling factor.\n\t# SCALEINT is expected to correspond to the intensity in q=0\n\t# SCALEINT=\n\n")

        for line in datfile:
            i=line.replace(';',' ').replace(',',' ').split()
            if is_number(i[0]) and is_number(i[1]):
                if options["-unit"].value == "nm":
                    q.append(float(i[0])/10.)
                elif options["-unit"].value == "Ang":
                    q.append(float(i[0]))
                else:
                    logging.error("Option unit not understood. Use -unit nm if your input file is in inverse nm. Use -unit Ang (default) if it is in inverse Angstrom")
                    sys.exit(1)
                expdata.append(float(i[1]))

        nqp=int(options["-nq"].value)
        if len(q)<nqp:
            for i in range(len(q)):
                plumedmain.write("\tQVALUE%d=%3.10f\tEXPINT%d=%3.10f\n"%(i+1,q[i],i+1,expdata[i]))
        else:
            i=0
            for k in range(0,nqp*int(len(q)/nqp),int(len(q)/nqp)):
                plumedmain.write("\tQVALUE%d=%3.10f\tEXPINT%d=%3.10f\n"%(i+1,q[k],i+1,expdata[k]))
                i+=1
        

        plumedmain.write("\n\t# METAINFERENCE\n\t# Uncomment the following keywords and adjust parameters to activate METAINFERENCE\n\t# DOSCORE NOENSEMBLE SIGMA_MEAN0=0\n\t# REGRES_ZERO=500\n\t# SIGMA0=5 SIGMA_MIN=0.001 SIGMA_MAX=5.00\n\t# NOISETYPE=MGAUSS\n\n... SAXS\n")
        plumedmain.write("\n# METAINFERENCE\n# Uncomment the following keyword to activate METAINF\n# saxsbias: BIASVALUE ARG=(saxsdata\.score) STRIDE=10\n# STATISTICS\nstatcg: STATS ARG=(saxsdata\.q_.*) PARARG=(saxsdata\.exp_.*)\n")
        plumedmain.write("\n# PRINT\n\n# Uncomment the following line to print METAINFERENCE output\n# PRINT ARG=(saxsdata\.score),(saxsdata\.biasDer),(saxsdata\.weight),(saxsdata\.scale),(saxsdata\.offset),(saxsdata\.acceptSigma),(saxsdata\.sigma.*) STRIDE=500 FILE=BAYES.SAXS\n\n# change stride if you are using METAINFERENCE\nPRINT ARG=(saxsdata\.q_.*) STRIDE=1 FILE=SAXSINT\nPRINT ARG=statcg.corr STRIDE=1 FILE=ST.SAXSCG\n")



if __name__ == '__main__':

    import sys,time

    start = time.time()
    stloc = time.localtime(start)
    logging.info("Start at time: {}:{}:{} of {}/{}/{}.".format(stloc[3],stloc[4],stloc[5],stloc[2],stloc[1],stloc[0]))

    args = sys.argv[1:]
    if not args:
        logging.error("NO INPUT! Try to use the option -h for help.")
        sys.exit(1)

    options = option_parser(args,options)

    if options["-f"].value is None:
        logging.error("No input PDB file. Try to use the option -h for help.")
        sys.exit(1)

    if options["-dat"].value is None:
        logging.error("No input data file. Try to use the option -h for help.")
        sys.exit(1)

    main(options)

    stop = time.time()
    stoploc = time.localtime(stop)
    logging.info("\n\nEnded at time: {}:{}:{} of {}/{}/{}.".format(stoploc[3],stoploc[4],stoploc[5],stoploc[2],stoploc[1],stoploc[0]))
    logging.info("Total time needed: {} sec\n".format(stop-start))
