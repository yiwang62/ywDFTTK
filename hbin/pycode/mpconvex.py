"""The purpose of this python code is for the convennience to download energy and structure from Materials Project and then calculate convex hull
   This code follow the linux command line style
   command line option:
   -within elist : elist = a list of elements such as "Al Mo Ni Ta Mo Nb"
   -containall elist : download all structues containing all "elist" within the "within elist" defined by option of "-within list"
   -containany elist : download all structues containing any elements in "elist" within the "within elist" defined by option of "-within list"
   -ehull arg : arg is a float number in eV/atom to control the output of downloaded structures
   -mpid list : list = a list of the materials project id such as "mp-1234 mp-1235 mp-1236"
   -pdir dir : dir = a folder that the downloaded VASP files will be saved

   output:
   A folder name "PDIRxx" where xx represent an automated number contain the VASP setting (named "xxx.VASP"=POSCAR "xxx.INCAR" "xxx.KPOINTS")
   total energy in eV/atom, formation energy, energy above convex hull if not stable, convex hull phases refered

   if no contain provided : the phase composition for the given within
"""

""" To run this code, you need to first install pymatgen follow the instruction from http://pymatgen.org/
   You may need first install anaconda https://www.anaconda.com/
"""
""" Make sure that you have the Materials API key.
To use the MAPI, you need to first have an API key. Please keep your API key safe and you should never share your key with anyone. If your key is compromised, you can regenerate the key again by visiting your dashboard page. you need first to login into your account with gmail or facebook,  them click the API tab followed by "dashborad" located under the "API keyes" heading. The click "Generate API Key". Finally copy it. 
"""
MY_API_KEY = 'YourMPkey'


import os
import shutil
import sys
import json
import numpy as np
from scipy.optimize import linprog
from numpy.linalg import solve
from fractions import Fraction
#import pymatgen as mg
from pymatgen import MPRester, Composition
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

def myjsonout(data,fp,indent="",comma=""):
	#print (data)
	mj = ''
	if (isinstance(data,dict)):
		fp.write('{}\n'.format('{'))
		#fp.write('{}'.format('{'))
			#sys.stdout.write('\n{}{}\n'.format(indent, '{'))
		nkey = 0
		for key in data:
			nkey += 1
			if nkey!=len(data):
				comma1 = ","
			else:
				comma1 = ""
			val = data[key]
			jval = json.dumps(val)
			jkey = json.dumps(key)
			#print (val)
			if (isinstance(val,dict)):
				fp.write('{}{}: '.format(indent+"    ",jkey))
				myjsonout(val,fp,indent+"    ",comma1)
			elif (isinstance(val,tuple)):
				#print (val)
				out = list(val)
				#print(out)
				fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, out, comma1))
			elif (isinstance(val,str)):
				if (indent == ""):
					#fp.write('{}{}: {}{}\n\n'.format(indent + "    ", jkey, jval, comma1))
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
			else:
				if (indent==""):
					#fp.write('{}{}: {}{}\n\n'.format(indent + "    ", jkey, jval, comma1))
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))
				else:
					fp.write('{}{}: {}{}\n'.format(indent + "    ", jkey, jval, comma1))

				#print(val)
				"""
				if (nkey!=len(data)):
					sys.stdout.write('{}{}: {},\n'.format(indent+"    ", key, val))
				else:
					sys.stdout.write('{}{}: {}\n'.format(indent+"    ", key, val))
				"""
		if comma==',':
			fp.write('{}{}{}\n\n'.format(indent,'}', comma))
		else:
			fp.write('{}{}{}\n'.format(indent, '}', comma))

"""convert decimal number q into fractional"""
def frac(q):
  a = Fraction(q).limit_denominator()
  if (a.denominator >1024):
    return('{:+.5f}*'.format(q))
  else:
    return('{:+.0f}/{:.0f}*'.format(a.numerator, a.denominator))

"""check if value is a float number"""
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

"""add 1.0 as summation constrain for the list of number in line
line - a list of composition
return:
a list lead lead by 1.0 followed by a list of composition
"""
def bMaker(line):
  b = [1.0]
  b.extend(line)
  b = np.array(list(map(float,b)))
  return (b)

"""combine phase fraction with phase names
X - phase composition
Phases - list of phase name
return:
a string with composition followed by phase name (multiple)
"""
def getPhase(X,Phases):
  decompose = " "
  for x in X:
    if (x[1] == 1.0):
      decompose = Phases[x[0]]
      break
    else:
      decompose = decompose + frac(x[1])+Phases[x[0]]
  return(decompose.replace(" +",""))

"""combine elemental symbol and phase composition into a string
Elements - 
compistion - 
return:
"""
def getPhCom(Elements,composition):
  decompose = " "
  for e,c in enumerate(composition):
    if (c == 1.0):
      decompose = Elements[e]
      break
    else:
      decompose = decompose + frac(c)+Elements[e]
  return(decompose.replace(" +",""))

"""convert a chemical formula into element list and composition
formula - chemical formula
return:
element list and composition list
"""
def formula2composition(formula):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'')
  newc = ""
  """Follow the convention, elemental symbol must start from capital letter"""
  for c in formula:
    if c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
      newc = newc + '|'
    newc = newc + c
  els = newc.split('|')
  els = [k for k in els if k != '']

  """now get the composition for each element"""
  ele = []
  com = []
  for el in els:
    newel = ""
    newcc = ""
    for c in el:
      if c.isalpha():
        newel = newel + c
      else:
        newcc = newcc + c

    if (newel not in periodictable):
      print('"',newel,'" is not an element! your formula is wrong!')
      sys.exit(1)
    ele.append(newel)

    if (len(newcc)!=0):
      if (isfloat(newcc)):
        com.append(float(newcc))
      else:
        print('"',newcc,'" is not a float number! your formula is wrong!')
        sys.exit(1)
    else:
      com.append(1.0)
  com = np.array(list(map(float,com)))
  com = com/sum(com)

  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=float)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]
      
  return elist,clist


MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0}


"""convert a chemical formula into element list and composition
formula - chemical formula
return:
element list and composition list
"""
def wCom2composition(formula):
  elist, clist = formula2composition(formula)
  for i,el in enumerate(elist):
    clist[i] = clist[i]/MM_of_Elements[el]
  clist = clist/sum(clist)
  return elist,clist


"""make the reduced formula
els - a list of elements
natype - a list of number of elements
return:
"""
def reduced_formula(els, nat):
  _els=els.split(' ')
  _els = [k for k in _els if k != '']
  _nat = np.array(list(map(int,nat)))
  #print(_els, _nat)
  els = sorted(set(_els))
  nat = np.zeros(len(els),dtype=int)
  for i,el in enumerate(_els):
    ix = els.index(el)
    nat[ix] += _nat[i]
  
  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  form = ""
  for j,el in enumerate(els):
    ix = nat[j]//i
    form = form+el
    if ix!=1:
      form = form+str(ix)
  return form

"""extract incar, kpoints, poscar, energy from bulk download data
dstack - data stack contain all information, such as phasename, energy, elements, etc
mstack - list of mp id
return:
"""
def bulk_extractMPdata(mstack,dstack):
  for rec in mstack:
    structure = rec['structure']
    poscar = structure.to(fmt="poscar")
    split_l = str(poscar).split('\n')
    nat=split_l[6].split(' ')
    nat = [int(k) for k in nat if k != '']
    natom=sum(nat)
    nCom = np.array(list(map(float,nat)))/natom

    energy = rec['final_energy']
    energy = energy/natom

    mpid = rec['task_id']
    spacegroup = rec["spacegroup"]
    nsites = rec["nsites"]
    p_formula = rec["pretty_formula"]
    pretty_formula = reduced_formula(split_l[5], nat)
    punitcell = makePCF(split_l[5], split_l[6])
    #inout = rec["input"]
    #incar = inout["incar"]
    #potcar = inout["potcar"]
    #potcar = inout["potcar_spec"]
    #potcar = "potcar"
    try:
      kpoints = inout["kpoints"]
    except:
      kpoints = "kpoints in error by materials project currently doe to Ag"
    incar = rec["input.incar"]
    potcar = rec["input.potcar_spec"]

    pname = spacegroup['crystal_system'] + "|" + str(spacegroup['number']) + "|" \
    + spacegroup['symbol'].replace('/','_') + "#" +  punitcell+'#'+mpid
    #pname = pretty_formula+'#'+mpid
    #pname = pretty_formula+'_'+mpid

    line = {"elements":split_l[5], "natype":split_l[6], "phasename":pname, "pretty_formula":p_formula,"energy":energy,
      "composition":nCom, "formation_energy":"N/A", "hull_energy":"N/A", "hull": "N/A",
      "INCAR":incar, "KPOINTS":kpoints, "POSCAR":poscar, "POTCAR":potcar, "NSITES":nsites}
    dstack.append(line)
    #sys.stdout.write("{} {} {} {} E0={:.6f} eV/atom, downloaded\n".format(mpid, split_l[5], split_l[6], pname, energy))
  sys.stdout.write("\n")

"""
make the primitive cell formula
"""
def makePCF(_els, _nat):
  els = _els.split(' ')
  els = [k for k in els if k != '']
  nat = _nat.split(' ')
  nat = [k for k in nat if k != '']
  elist = sorted(set(els))
  clist = np.zeros(len(elist), dtype=int)
  for j,el in enumerate(els):
    ix = elist.index(el)
    clist[ix] += int(nat[j])

  punitcell = ""
  for i,el in enumerate(elist):
    punitcell += el+str(clist[i])
  return (punitcell)


def ExcelPhaseName():
  df = pd.read_excel('2019-02-08-List-Phases.xlsx', sheet_name='model_MP')
  #df = pd.read_excel('2019-02-22-List-Phases.xlsx', sheet_name='model_MP')
  columns = df.columns.values
  for cc in columns:
    if "Phases".upper() == cc.strip().upper() : phases = cc
    elif "Ratio".upper() == cc.strip().upper() : ratio = cc
    elif "Spacegroup".upper() == cc.strip().upper() : spacegroup = cc
  phases = df[phases]
  ratio = df[ratio]
  spacegroup = df[spacegroup]
  return phases,ratio,spacegroup

def natratio(nat):
  nat = np.array(list(map(int,nat)))
  Nd = min(nat)
  for i in range(Nd,0,-1):
    n = 0
    for x in nat:
      n += x%i
    if n==0: break
  nat = nat//i
  return nat

def GetPhaseName(phasename, ratio, PhaseNameList, SublatticeRatio, SpaceGroup):
  ratio = natratio(ratio)
  pname = phasename.replace('|','+').replace('#','+')
  spacegroup = int(phasename.split('+')[1])
  for i,pp in enumerate(PhaseNameList):
    try:
      sps = int(pp.strip().split('_')[-1])
    except:
      try:
        #print ("************",SpaceGroup[i])
        sps = int(SpaceGroup[i])
        #print (sps)
      except:
        continue
    if spacegroup!=sps: continue
    sratio = SublatticeRatio[i].strip().replace('[','').replace(']','').replace(' ','').split(',')
    sratio = np.array(list(map(int,sratio)))
    #print(sorted(ratio),sorted(sratio))
    if sorted(ratio)!=sorted(sratio): continue
    return pp
  return pname

def outESPEI(posdir,_dstack,ehull):
  jdir = posdir+"json/"
  if not os.path.exists(jdir):
    os.mkdir(jdir)
  PhaseNameList, SublatticeRatio, SpaceGroup = ExcelPhaseName()
  outE = []
  for ff in _dstack:
    ff0 = ff["phasename"].replace('|','#').split('#')
    ff["phasename"] = '+'.join(ff0)
    if ff["hull_energy"] > ehull: continue
    rec = {}
    rec["components"] = [ss.upper() for ss in ff["elements"].split(" ")]
    solver = {}
    solver["mode"] = "manual"
    solver["sublattice_site_ratios"] = [int(n) for n in ff["natype"].split(" ")]
    solver["sublattice_configurations"] = [[ss.upper() for ss in ff["elements"].split(" ")]]
    rec["solver"] = solver 
    conditions = {}
    conditions["P"] = 0
    conditions["T"] = 298.15
    rec["conditions"] = conditions
    rec["output"] = "HM_FORM"
    rec["values"] = [[[float('{:.6f}'.format(ff["formation_energy"]*96.485))]]]
    #rec["values"] = [[[ff["formation_energy"]*96.485]]]
    rec["reference"] = "Materials Project "+ff["phasename"].split("#")[-1]
    structure = ff["phasename"].split('+')
    rec["structure"] = {"crystal system":structure[0], "spacegroup":int(structure[1]), "Pearson symmetry":structure[2], "primitive unit cell":structure[3], "Materials Project id":structure[4]}
    rec["fomula"] = ff["pretty_formula"]
    rec["hull_energy"] = ff["hull_energy"]
    rec["POSCAR"] = ff["POSCAR"]
    rec["INCAR"] = ff["INCAR"]
    rec["KPOINTS"] = ff["KPOINTS"]
    rec["POTCAR"] = ff["POTCAR"]
    #rec["phases"] = [ff["phasename"]]
    rec["phases"] = [GetPhaseName(ff["phasename"], solver["sublattice_site_ratios"], PhaseNameList, SublatticeRatio,SpaceGroup)]
    spname = rec["phases"][0]
    ff0 = ff["phasename"].split('+')
    if not spname.endswith(ff0[-1]):
      spname = ff0[-2]+'+'+ff0[2].replace('_','')+'+'+spname+'+'+ff0[1]+'+'+ff0[-1]
    else:
      spname = ff0[-2]+'+'+ff0[2].replace('_','')+'+'+ff0[0]+'+'+ff0[1]+'+'+ff0[-1]

    ff["phases"] = spname
    outE.append(rec)

    #with open(jdir + rec["phases"][0]+'+'+rec["structure"]["primitive unit cell"]+".json", 'w') as fp:
    with open(jdir + spname+".json", 'w') as fp:
      myjsonout(rec, fp, indent="", comma="")
      if not True:
        myjsonout(rec, sys.stdout, indent="", comma="")

"""
    line = {"elements":split_l[5], "natype":split_l[6], "phasename":pname, "pretty_formula":p_formula,"energy":energy,
      "composition":nCom, "formation_energy":"N/A", "hull_energy":"N/A", "hull": "N/A",
      "INCAR":incar, "KPOINTS":kpoints, "POSCAR":poscar, "POTCAR":potcar, "NSITES":nsites}
"""


"""output POSCAR, INCAR, KPOINTS and stability of each download structure
posdir - dir for files to be outputted into
dstack - data stack contain all information, such as phasename, energy, elements, etc
calhull - bool to instruct if calculatiing convex hull
ehull - how high above the hull the structures will be included in the output
return:
""" 

def outMPdata(posdir, _dstack, calhull, ehull):
  outESPEI(posdir,_dstack,ehull)
  print ("Compound         Formula           SG Symmetry Tenergy(eV/atom) Fenergy : Primitive POSCAR")
  dstack = []
  for ff in _dstack:
    if ff["hull_energy"] < 1.e-6:
      dstack.append(ff)
  for ff in _dstack:
    if not ff["hull_energy"] < 1.e-6:
      dstack.append(ff)

  for ff in dstack:
    pname = ff["phasename"]
    energy = ff["energy"]
    pformula = ff["pretty_formula"]
    out = (pname.split('#')[0]).split('+')
    
    punitcell =  makePCF(ff["elements"], ff["natype"])
    if (not calhull):
      sys.stdout.write("{} {} {:.6f}\n".format(punitcell, pname, energy))
    else:
      eh = ff["hull_energy"]
      ef = ff["formation_energy"]
      
      if eh < ehull:
        if eh < 1.e-6:
          sys.stdout.write("{:<16s} {:<16s} {:>3s} {:10s} {:11.6f} {:10.6f} : ".format(punitcell, pformula, out[1],out[2], energy, ef))
        else:
          sys.stdout.write("{:>16s} {:>16s} {:>3s} {:10s} {:11.6f} {:10.6f} above hull {:.6f} : "
            .format(punitcell, pformula, out[1],out[2], energy, ef, eh))
          #sys.stdout.write("{:<10s} {:<3s} {:10s} {:11.6f} {:10.6f} {:.6f} above hull {}\n"
          #  .format(punitcell,  out[1],out[2], energy, ef, eh, ff["hull"]))
      else:
        continue

    spname = posdir + ff["phases"]
    print(spname+".VASP")
    with open(spname+".VASP",'w') as f0:
      f0.write("{}\n".format(ff["POSCAR"]))
    kustoutfiles(spname+".VASP")

    with open(spname+".KPOINTS",'w') as f0:
      f0.write("{}\n".format(ff["KPOINTS"]))

    with open(spname+".POTCAR",'w') as f0:
      f0.write("{}\n".format(ff["POTCAR"]))

    with open(spname+".INCAR",'w') as f0:
      incar = ff["INCAR"]
      for key, value in incar.items() :
        f0.write("{}= ".format(key))
        if type(value) is list:
          for v in value:
            f0.write("{} ".format(v))
        elif type(value) is bool:
          f0.write(".{}.".format(value))
        else:
          f0.write("{}".format(value))
        f0.write("\n")

"""exclude elements from the periodictable for the purpose of $nin operation in Mogodb"""
def excludeelements(within):
  els = []
  for el in periodictable:
    if el not in within:
      els.append(el)
  #bugs in Materials project
  return els

"""download records for  pure elements
within - list of elements to be downloaded
excludeelement - list of elements to be excluded from within
dstack - stack for data to be saved excludeelementing all information, such as phasename, energy, elements, etc
return:
dstack
"""
def bulk_get_elements(excludeelement, elements):
  dstack = []
  els = []
  for el in elements:
    if el not in excludeelement:
      els.append(el)
      #result = mpr.query(criteria={"elements":{"$eq":[el]}}, properties=["task_id","structure", "input", "final_energy"])
  result = mpr.query(criteria={"elements":{"$in":els}, "nelements":1}, properties=properties)
  dstack.extend(result)
  return dstack

"""bulk download records for structure with non of their elements contained not in nlist and contained in contain
excludeelements - list of elements for $nin condition
contain - list of elements to be used by the containoperator
within - list of elements for all compounds within it to be downloaded
containoperator - $all or $in
dstack - stack for data to be saved containing all information, such as phasename, energy, elements, etc
return:
"""
def bulk_getPOSj(contain,excludeelements,dstack, within, containoperator):
  data=bulk_get_elements(contain, within)
  result = mpr.query(criteria={"elements":{"$nin":excludeelements, containoperator:contain}}, properties=properties)
  data.extend(result)
  ndata = len(data)
  bulk_extractMPdata(data,dstack)
  print(ndata, " Compounds found for the system\n")

"""bulk download records for structure with non of their elements contained in nlist
excludeelements - list of elements for $nin condition
dstack - stack for data to be saved containing all information, such as phasename, energy, elements, etc
return:
"""
def bulk_getPOSl(excludeelements,dstack):
  #print(excludeelements)
  data = mpr.query(criteria={"elements":{"$nin":excludeelements}}, properties=properties)
  bulk_extractMPdata(data,dstack)
  ndata = len(data)
  print(ndata, " Compounds found for the system\n")

"""bulk download records with mp id
mstack - list of mp id
dstack - stack for data to be saved containing all information, such as phasename, energy, elements, etc
return:
"""
def bulk_getMPdata(mstack,dstack):
  data = mpr.query(criteria={"task_id":{"$in":mstack}}, properties=properties)
  bulk_extractMPdata(data,dstack)
  ndata = len(data)
  print(ndata, " Compounds found for the system\n")

#data = mpr.query(criteria={"elements":{"$in":["Au", "Ir","Os", "Pd", "Pt", "Rh", "Ru"], "$in":["Ag"]}, "nelements":5}, properties=["task_id"])
#criteria='{"elements":{"$in":["Li", "Na", "K"], "$all": ["O"]}, "nelements":2}'
#tr -cd '\11\12\15\40-\176' <

"""make the convex hull and find the stability of all phases
dstack - data stack contain all information, such as phasename, energy, elements, etc
return:
"""
def convexhull(_dstack,fastcode=False):
  Elements = []
  Phases = []
  Gstack = []
  dstack = []
  for ff in _dstack:
    els = ff["elements"].split(' ')
    els = [k for k in els if k != '']
    els = sorted(sorted(set(els)))
    for el in els:
      if el not in Elements:
        Elements.append(el)
    if len(els)==1:
      if els[0] in suspend:
        continue
    dstack.append(ff)
    Phases.append(ff["phasename"])
    Gstack.append(ff["energy"])

  Elements = sorted(set(Elements))
  print('#Elemental components:', Elements)

  amatrix = []
  Amatrix = []
  refG = np.zeros(len(Elements), dtype=float)

  """Make the A matrix for simplex algebra"""
  for ff in dstack:
    tComponents = ff["elements"].split(' ')
    tComponents = [k for k in tComponents if k != '']
    tnComponents = np.array(list(map(float,ff["composition"])))
    Components = sorted(set(tComponents))
    nComponents = np.zeros(len(Components))
    for i0,el in enumerate(tComponents):
      ix = Components.index(el)
      nComponents[ix] = nComponents[ix] +  tnComponents[i0]

    aaa = []
    Amatrix.append(1.0)
    for i,element in enumerate(Elements):
      if element in Components:
        j = Components.index(element)
        aaa.append(nComponents[j])
        Amatrix.append(nComponents[j])
      else:
        aaa.append(0.0)
        Amatrix.append(0.0)

    amatrix.append(aaa)
    for m,el in enumerate(Elements):
      if (float(aaa[m]) == 1.0):
        refG[m] = min(refG[m], float(ff["energy"]))

  row= len(Phases)
  col= len(Elements)
  
  Gstack = np.array(list(map(float,Gstack)))
  Amatrix = np.array(list(map(float,Amatrix))).reshape(row,col+1)
  Amatrix=Amatrix.T

  if fastcode:
    fastGstack, fastPhases, fastAmatrix = fast_preparing(amatrix, Elements, Gstack, Phases, Amatrix)

  ndata = 0
  for j,line in enumerate(amatrix):
    """using simplex to find if is a stable phase""" 
    comp = np.array(list(map(float,line)))
    if fastcode:
      """optimizing code by remove redunt conditions"""
      res = fast_linprog(fastGstack, fastAmatrix, comp)
      x = [[i,q] for i,q in enumerate(res.x) if (q>1.e-8)]
      hull = getPhase(x,fastPhases)
    else:
      """use normal minimization"""
      b = bMaker(comp)
      res = linprog(Gstack, Amatrix, b, bounds=(0.0, 1.0),method='Simplex')
      x = [[i,q] for i,q in enumerate(res.x) if (q>1.e-8)]
      hull = getPhase(x,Phases)

    """deltaF is the formation energy"""
    deltaF = (float(Gstack[j]) - sum(comp*refG))
    dstack[j]["formation_energy"] = deltaF
    x2 = sum(res.x*res.x)
    if (abs(x2-1.0) < 1.e-4) and abs(res.fun-Gstack[j])<1.e-6:
      """stable phase"""
      dstack[j]["hull"] = ""
      dstack[j]["hull_energy"] = 0.0
      ndata += 1
    else:
      """unstable phase, give the phase separation"""
      dstack[j]["hull"] = hull
      """unstable phase, energy above the hull"""
      dstack[j]["hull_energy"] = (Gstack[j]-res.fun)
    print(j,"/",len(amatrix), " Phases: ", Phases[j], "Total energy=",Gstack[j])
  print("\n", ndata, " Compounds made the convex hull for the system\n")
  return dstack


def datacheck(_dstack):
  pots = {}
  for rec in _dstack:
    lines = rec["POTCAR"]
    for line in lines :
      elp = [k for k in line['titel'].split(' ') if k != '']
      pot = elp[0]+'/'+elp[1]
      if pot in pots.keys():
        pots[pot] += 1
      else:
        pots[pot] = 1
  print ("\nPOTCAR staticstics :", pots)
  pp = {}
  pot = {}
  for key in pots.keys():
    ppp = key.split('/')[0]
    pott = key.split('/')[1]
    pp[ppp] = 1
    el = pott.split('_')[0]
    if el in pot.keys():
      pot[el] = pot[key].append(pott)
    else:   
      pot[el] = [pott]
  if len(pp) > 1:
    print("************WARNING!", pp, "pseudopotentila Messed up")
  for el in pot.keys():
    if len(pot[el]) > 1:
      print("************WARNING!", pot[el], "pseudopotentila Messed up")
  print('\n')

    
"""attempt to speed up the speed. I am not sure if I will continue on it"""
def fast_preparing(amatrix, Elements, Gstack, Phases, Amatrix):
  comPhases = [getPhCom(Elements,l) for l in amatrix]
  g = []
  ccg = []
  nx = []
  for j,gg in enumerate(Gstack):
    cc = comPhases[j]
    if cc in ccg:
      ix = ccg.index(cc)
      if gg<g[ix]:
        g[ix] = gg
        nx[ix] = j
    else:
      ccg.append(cc)
      g.append(gg)
      nx.append(j)

  nx = np.array(list(map(int,nx)))
  fastGstack = Gstack[nx]
  fastPhases = [Phases[k] for k in nx]
  fastAmatrix = Amatrix[:,nx]
  return fastGstack, fastPhases, fastAmatrix

def fast_linprog(Gstack, Amatrix, line):
  ny = []
  b = bMaker(line)
  zzz = np.zeros(len(b))
  for i,l in enumerate(b):
    if l==0.0:
      zzz[i] = 1.0
    else:
      ny.append(i)
  ny = np.array(list(map(int,ny)))

  nr = []
  for j in range(len(Gstack)):
    aaa = Amatrix[:,j]
    if sum(zzz*aaa)==0.0:
      nr.append(j)
  nr = np.array(list(map(int,nr)))

  nb = b[ny]
  nrg = Gstack[nr]
  nra = Amatrix[:,nr]
  nra = nra[ny,:]

  res = linprog(nrg, nra, nb, bounds=(0.0, 1.0),method='Simplex')
  class res:
    x = np.zeros(len(Gstack))
    x[nr] = res.x
    fun = res.fun
  return res 


"""find the phase combination for a given composition
dstack - data stack
formula - formula
elist - element list
clist - composition list
return:
"""
def calcphasedecompostion(dstack,formula,elist,clist):
  Elements = []
  Phases = []
  Gstack = []
  for ff in dstack:
    Phases.append(ff["phasename"])
    Gstack.append(ff["energy"])
    els = ff["elements"].split(' ')
    els = [k for k in els if k != '']
    for el in els:
      if el not in Elements:
        Elements.append(el)

  Elements = sorted(set(Elements))

  blist = np.zeros(len(Elements))
  for i,el in enumerate(elist):
    if el in Elements:
      ix = Elements.index(el)
      blist[ix] = clist[i]

  Amatrix = []

  for ff in dstack:
    """accumulate the compostion"""
    tComponents = ff["elements"].split(' ')
    tComponents = [k for k in tComponents if k != '']
    tnComponents = np.array(list(map(float,ff["composition"])))
    Components = sorted(set(tComponents))
    nComponents = np.zeros(len(Components))
    for i0,el in enumerate(tComponents):
      ix = Components.index(el)
      nComponents[ix] = nComponents[ix] +  tnComponents[i0]

    """accumulate the simplex matrix"""
    Amatrix.append(1.0)
    for i,element in enumerate(Elements):
      if element in Components:
        j = Components.index(element)
        Amatrix.append(nComponents[j])
      else:
        Amatrix.append(0.0)

  row= len(Phases) #"""number of phases"""
  col= len(Elements) #"""number of elements"""
  
  Amatrix = np.array(list(map(float,Amatrix))).reshape(row,col+1)
  Amatrix=Amatrix.T

  b = bMaker(blist)
  """call the simplex algebra"""
  res = linprog(Gstack, Amatrix, b, bounds=(0.0, 1.0),method='Simplex')
  """X is the phase compostion from simplex"""
  X = [[i,q] for i,q in enumerate(res.x) if (q>1.e-8)]
  #print("\n", formula, " is made of ",getPhase(X,Phases),"\n")

  pcom = ""
  print("\nBy atomic percentages, the phases are:\n")
  for x in X:
    phases = Phases[x[0]].split("#")
    print('{:6.2f}  {:12s}  {:12s}  {}'.format(x[1]*100, phases[1], phases[2], phases[0]))
    if pcom == "":
      pcom = '{:6.2f}*{}'.format(x[1]*100, phases[1])
    else:
      pcom += ' +' + '{:6.2f}*{}'.format(x[1]*100, phases[1])
      
  print("\n", formula, " is made of ", pcom, "\n")

  print("\n")

def setoutdir(pdir):
  if pdir=="":
    pdir = "MP"
    i = 0
    posdir = pdir + str(i)
    #print ( os.path.isdir(posdir), os.path.isfile(posdir), os.path.exists(posdir) )
    while ( os.path.exists(posdir) ):
      i += 1
      posdir = pdir + str(i)
    os.mkdir(posdir)
  elif (os.path.isdir(pdir)):
    posdir = pdir
  elif os.path.isfile(pdir):
    print("file with the same name exist! please give another name!")
    sys.exit(1)
  elif not os.path.exists(pdir):
    posdir = pdir
    os.mkdir(posdir)
  return posdir+"/"


"""unit conversion"""
eVtoGPa = 160.21766208
eVtoJ = 96486.9

MM_of_Elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'ZERO': 0}

"""property list to be downloaded"""
#properties = ["task_id","structure", "input", "nsites", "potcar", "spacegroup", "pretty_formula", "final_energy"]
properties = ["prototype", "task_id", "structure", "input.incar", "input.potcar_spec", "nsites", "potcar", "spacegroup", "pretty_formula", "final_energy"]
"""for code developing test"""
#within = ["Ag", "Au", "Ir","Os", "Pd", "Pt", "Rh", "Ru"]
within = ['Pt', 'Ir', 'Rh', 'Ni', 'Zr', 'Hf', 'Si', 'Cr']
contain = []
"""handle command line option"""

fastcode = True #""" key to cotrol the search method"""
#periodictable = [key for key,value in MM_of_Elements.items()] #""" list of all elements from the periodic table"""
periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""
ehull = 1.e-6
input_mpid = False #""" key to cotrol the mp id input"""
calhull = False #""" key to cotrol the if calculating the convex hull"""
pdir = "MP0"
suspend = []

mpid = []
input_within = False #""" key to cotrol the within input"""
input_wCom = False #""" key to cotrol the within input"""
input_contain = False #""" key to cotrol the contain input"""
formula_within = "" #"""chemical formula"""
formula_wCom = "" #"""chemical formula"""
formula_contain = "" #"""chemical formula"""
containoperator = "$all" #contain any elements within the list of contain
KUST0 = "SSS/"

count = 1
while (count < len(sys.argv)):
  if (sys.argv[count] == "-within"):
    """for code developing test"""
    within = ['Fe', 'Ni', 'Al', 'C', 'Co', 'Cr', 'Cu', 'Hf', 'Mn', 'Mo', 'Nb', 'Re', 'Si', 'Ta', 'Ti', 'W']
    input_within = True
    input_mpid = False
    input_contain = False
    input_wCom = False
  elif (sys.argv[count] == "-wCom"):
    input_wCom = True
    input_within = False
    input_mpid = False
    input_contain = False
  elif (sys.argv[count] == "-containall"):
    containoperator = "$all"
    input_contain = True
    input_within = False
    input_mpid = False
    input_wCom = False
  elif (sys.argv[count] == "-containany"):
    containoperator = "$in"
    input_contain = True
    input_within = False
    input_mpid = False
    input_wCom = False
  elif (sys.argv[count] == "-mpid"):
    input_mpid = True
    input_within = False
    input_contain = False
    input_wCom = False
  elif (sys.argv[count] == "-fastcode"):
    fastcode = not fastcode
  elif (sys.argv[count] == "-kust"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    KUST0 = sys.argv[count].replace('/','') + '/'
  elif (sys.argv[count] == "-suspend"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    suspend.append(sys.argv[count])
  elif (sys.argv[count] == "-ehull"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    ehull = float(sys.argv[count])
  elif (sys.argv[count] == "-pdir"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    pdir = str(sys.argv[count])
  elif (input_within):
    formula_within = formula_within+sys.argv[count]
  elif (input_wCom):
    formula_wCom = formula_wCom+sys.argv[count]
  elif (input_contain):
    formula_contain = formula_contain+sys.argv[count]
  elif (input_mpid):
    if (count >= len(sys.argv)):
      break
    id = sys.argv[count]
    mpid.append(str(id))
  else:
    print ("************* UNKOWN option",'"',sys.argv[count],'"')
    sys.exit(1)
  count = count + 1

def kustdirlist():
  if not os.path.exists(KUST0):
    os.mkdir(KUST0)
  list = os.listdir(KUST0)
  list = [i for i in list if i.endswith('.VASP')]
  return list

def kustoutfiles(spname):
  ff0 = spname.split('/')[1].split('+')
  ff0 = ff0[0]+'+'+ff0[1]
  for file in kustfiles:
    if file.startswith(ff0) : return

  ff0 = spname.split('/')[0]+'/'
  shutil.copy(spname, spname.replace(ff0,KUST0))

kustfiles = kustdirlist()

mpr = MPRester(MY_API_KEY)
posdir = setoutdir(pdir)

dstack = []
if formula_wCom!="":
  within,compos_within = wCom2composition(formula_wCom)
if formula_within!="":
  within,compos_within = formula2composition(formula_within)
if formula_contain!="":
  contain,compos_contain = formula2composition(formula_contain)

if len(mpid)!=0:
  #if mp id's are given by -formula command line option"""
  bulk_getMPdata(mpid, dstack)
  outMPdata(posdir, dstack, False, ehull)
elif len(within)!=0 and len(contain)!=0:
  #if both -with and -contain elements are given by -contain command line option
  bulk_getPOSj(contain, excludeelements(within), dstack, within, containoperator)
  hstack = convexhull(dstack,fastcode)
  outMPdata(posdir, hstack, True, ehull)
  datacheck(hstack)
elif len(within)!=0:
  #if compostion is given by -formula command line option
  bulk_getPOSl(excludeelements(within), dstack)
  hstack = convexhull(dstack,fastcode)
  outMPdata(posdir, hstack, True, ehull)
  datacheck(hstack)
  try: 
    compos_within
    calcphasedecompostion(hstack,formula_within,within,compos_within)
  except:
    pass

