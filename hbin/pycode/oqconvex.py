import os
import shutil
import string
import numpy as np
from scipy.optimize import linprog
import sys
import urllib.parse
import requests
from urllib.request import urlopen,Request,urlretrieve

printable = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&\'()*+,-./:;<=>?@[\\]^_`{|}~ \n"
cell = 'Primitive'

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


import MyPhaseClass
def outESPEI(posdir,_dstack,ehull):
  jdir = posdir+"json/"
  if not os.path.exists(jdir):
    os.mkdir(jdir)
  PhaseNameList, SublatticeRatio, SpaceGroup = MyPhaseClass.ShunliExcelPhaseName()
  outE = []
  for ff in _dstack:
    ff0 = ff["phasename"].replace('|','#').split('#')
    ff["phasename"] = '+'.join(ff0)
    if ff["hull_energy"] > ehull: continue
    rec = {}
    rec["components"] = [ss.upper() for ss in ff["elements"].split(" ")]
    solver = {}
    solver["mode"] = "manual"
    #solver["sublattice_site_ratios"] = [int(n) for n in ff["composition"].split(" ")]
    solver["sublattice_site_ratios"] = ff["composition"]
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

"""check if value is a float number"""
def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

"""
reorder formula
"""
def punitcell_formula(formula,N):
  elist,clist = formula2composition(formula)
  punitcell = ""
  for i,el in enumerate(elist):
    punitcell += el+str(int(clist[i]*float(N)))
  return (punitcell)

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


def setoutdir(pdir):
  if pdir=="":
    pdir = "OQ"
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

def downloadsetting(setting,pos,formula):
  req = urllib.request.Request(setting)
  xmlstr = ''
  for ii,bline in enumerate(urlopen(req).readlines()):
    line = bline.decode("utf-8").strip()
    if line=='<div id="sidebar" class="box">':
      xmlstr += line
    elif xmlstr!='':
      xmlstr += line
      if line=='</div>': break
    
  lines = xmlstr.replace('<br />','\n').strip().replace('  <h2>','#').replace('<p>','').replace(' </h2>','').replace('</p>','\n')
  lines = lines.replace('</code><h2>','\n#').replace('<code>','\n').replace('</code></div>','\n')
  lines = [l.strip() for l in lines.split('\n')]
  outname = posdir+formula+'+oq-'+pos+".PAR"
  f = open(outname, 'w')

  for line in lines:
    if line.find("xml")!=-1 : continue
    if line.find("<div")!=-1 : continue
    if line.find("</div")!=-1 : continue
    if line.find("<code")!=-1 : continue
    if line.find("</code")!=-1 : continue
    f.write(line)
    f.write('\n')
  f.close()

def downloadstructure(entry,formula):
  req = urllib.request.Request('http://oqmd.org/materials/entry/'+entry)
  blines = urlopen(req).readlines()
  outname = ""
  for bb,bline in enumerate(blines):
    line = bline.decode("utf-8").strip()
    if line.startswith('<a href="/materials/export/primitive/poscar/'):
      pos = line.replace('<a href="/materials/export/primitive/poscar/',"").split('"')[0]
      if cell=='Primitive':
        urlpos = 'http://oqmd.org/materials/export/primitive/poscar/'+pos
        outname = posdir+formula+'+oq-'+pos+".VASP"
        with open(outname, 'w') as f:
          req = Request(urlpos)
          lines = urlopen(req).readlines()
          for ii,line in enumerate(lines):
            #f.write(line.decode('utf-16').encode('utf-8'))
            line = line.decode("utf-8")
            for c in line:
              if c in printable: f.write(c)
            """
            """
        kustoutfiles(outname)
      else:
        urlpos = 'http://oqmd.org/materials/export/conventional/poscar/'+pos
        urllib.request.urlretrieve(urlpos, posdir+pos+".VASP")
        sys.stdout.write(" : "+pos+".VASP")
    elif line.startswith('<tr class=clickableRow href="/analysis/calculation/'):
      static = blines[bb+1].decode("utf-8").strip()
      if static=='<td>static</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        setting = line.replace('<tr class=clickableRow href="','http://oqmd.org/').split('"')[0]
        print(' {:>8s} : {}'.format(energy,outname))
        downloadsetting(setting,pos,formula)
        return energy
      elif static=='<td>standard</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        setting = line.replace('<tr class=clickableRow href="','http://oqmd.org/').split('"')[0]
        print(' {:>8s} : standard {}'.format(energy,outname))
        downloadsetting(setting,pos,formula)
        return energy
      elif static=='<td>static_lda</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        setting = line.replace('<tr class=clickableRow href="','http://oqmd.org/').split('"')[0]
        print(' {:>8s} : static_lda {}'.format(energy,outname))
        downloadsetting(setting,pos,formula)
        return energy
  print(' WARNING no setting FOUND :', outname)

def downloadtotalenergy(entry,formula,materialsid):
  req = urllib.request.Request('http://oqmd.org/materials/entry/'+entry)
  blines = urlopen(req).readlines()
  for bb,bline in enumerate(blines):
    line = bline.decode("utf-8").strip()
    if line.startswith('<a href="/materials/export/primitive/poscar/'):
      pos = line.replace('<a href="/materials/export/primitive/poscar/',"").split('"')[0]
      if cell=='Primitive':
        urlpos = 'http://oqmd.org/materials/export/primitive/poscar/'+pos
        outname = posdir+formula+'_oq-'+pos+".VASP"
    elif line.startswith('<tr class=clickableRow href="/analysis/calculation/'):
      static = blines[bb+1].decode("utf-8").strip()
      if static=='<td>static</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        if os.path.exists(outname):
          print(' {:>8s} : {}'.format(energy,outname))
        else:
          print(' {:>8s} : {}'.format(energy,'static'))
        return
      elif static=='<td>static_lda</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        print(' {:>8s} : {}'.format(energy,'static_lda'))
        return
      elif static=='<td>standard</td>':
        energy = blines[bb+2].decode("utf-8").strip().replace('<td>','').replace('</td>','')
        if os.path.exists(outname):
          print(' {:>8s} : {}'.format(energy,outname))
        else:
          print(' {:>8s} : {}'.format(energy,'standard'))
        return
  print(' : WARNING no data FOUND!')



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
    nComponents /= nComponents.sum()

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

  pcom = ""
  print("\nBy atomic percentages, the phases are:\n")
  for x in X:
    #phases = Phases[x[0]].split("#")
    phases = Phases[x[0]]
    print('{:6.2f}  {:12s}'.format(x[1]*100, phases))
    if pcom == "":
      pcom = '{:6.2f}*{}'.format(x[1]*100, phases)
    else:
      pcom += ' +' + '{:6.2f}*{}'.format(x[1]*100, phases)

  print("\n", formula, " is made of ", pcom, "\n")


periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""

within = ['Fe', 'Ni']
within = ['Pt', 'Ir', 'Rh', 'Ni', 'Zr', 'Hf', 'Si', 'Cr']
fastcode = True #""" key to cotrol the search method"""
ehull = 1.e-6
calhull = False #""" key to cotrol the if calculating the convex hull"""
suspend = []

input_within = False #""" key to cotrol the within input"""
input_wCom = False #""" key to cotrol the within input"""
input_contain = False #""" key to cotrol the contain input"""
formula_within = "" #"""chemical formula"""
formula_wCom = "" #"""chemical formula"""
formula_contain = "" #"""chemical formula"""
containoperator = "$all" #contain any elements within the list of contain
KUST0 = "SSS/"

downloadenergy = False
deltaH = True

alldata = False

count = 1
while (count < len(sys.argv)):
  if (sys.argv[count] == "-within"):
    """for code developing test"""
    within = ['Fe', 'Ni']
    within = ['Pt', 'Ir', 'Rh', 'Ni', 'Zr', 'Hf', 'Si', 'Cr']
    within = ['Fe', 'Ni', 'Al', 'C', 'Co', 'Cr', 'Cu', 'Hf', 'Mn', 'Mo', 'Nb', 'Re', 'Si', 'Ta', 'Ti', 'W']
    input_within = True
    input_contain = False
    input_wCom = False
  elif (sys.argv[count] == "-kust"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    KUST0 = sys.argv[count].replace('/','') + '/'
  elif (sys.argv[count] == "-energy"):
    downloadenergy = True
  elif (sys.argv[count] == "-alldata"):
    alldata = not alldata
  elif (sys.argv[count] == "-deltaH"):
    deltaH = not deltaH
  elif (sys.argv[count] == "-wCom"):
    input_wCom = True
    input_within = False
    input_contain = False
  elif (sys.argv[count] == "-containall"):
    containoperator = "$all"
    input_contain = True
    input_within = False
    input_wCom = False
  elif (sys.argv[count] == "-containany"):
    containoperator = "$in"
    input_contain = True
    input_within = False
    input_wCom = False
  elif (sys.argv[count] == "-fastcode"):
    fastcode = not fastcode
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


posdir = setoutdir("OQ0")
if formula_within=="":
  for el in within:
    formula_within += el
if formula_wCom!="":
  within,compos_within = wCom2composition(formula_wCom)
if formula_within!="":
  within,compos_within = formula2composition(formula_within)
if formula_contain!="":
  contain,compos_contain = formula2composition(formula_contain)

system='-'.join(within)

Stable = "Stable phases"
Allphases = "Compounds contained in this region of phase space"

req = Request('http://oqmd.org/materials/composition/'+system)
blines = urlopen(req).readlines()
dld = False
itot = 0
istable = 0

hstack = []
for bidx,bline in enumerate(blines):
  line = bline.decode("utf-8").strip() 
  if line.find(Stable)!=-1:
    print("\n downloading",Stable,"for", system,"\n")
    print("Compound  Symmetry Fenergy(eV/atom)  Prototype  Cell_Size  Tenergy :",cell, "POSCAR\n")
    dld = True != alldata
    irange = 5
    istart = 0
    continue
  if line.find(Allphases)!=-1:
    print("\n ",Allphases,"\n")
    #dld = False
    dld = True == alldata
    istart = 1
    irange = 6
    continue

  if line.startswith('<tr class=clickableRow href="/materials/entry/'):
    itot += 1
    entry = line.replace('<tr class=clickableRow href="/materials/entry/',"").replace('">','')
    out = []
    for i in range(irange):
      bb = blines[bidx+1+i].decode("utf-8").strip().replace("<td>","").replace("</td>","").replace("<sub>","").replace("</sub>","")
      out.append(bb)
    materialid = out[0]
    out = out[istart:]
    formula = out[0]
    punitcell= punitcell_formula(out[0],out[4])
    els,com = formula2composition(formula)
    sys.stdout.write('{:<16s} {:<10s} {:>7s} {:>16s} {:>4s}'.format(punitcell,out[1],out[2],out[3],out[4]))
      
    if dld:
      #energy = downloadstructure(entry,formula)
      energy = downloadstructure(entry,punitcell+'+'+out[1].replace('/','').replace('_','')+'+'+out[3].replace('/',''))
      istable += 1
      if deltaH:
        drec = {"hull_energy":0,"formula":formula, "entry":entry, "elements":' '.join(els), "phasename":punitcell, "composition":com, "energy":float(out[2])-1.e-7}
      else:
        drec = {"hull_energy":0,"formula":formula, "entry":entry, "elements":' '.join(els), "phasename":punitcell, "composition":com, "energy":float(energy)}
      hstack.append(drec)
    elif downloadenergy:
      downloadtotalenergy(entry,formula,materialid)
    else:
      #print(" : "+bb+"#oq-"+entry)
      print(" : id=", materialid)

print ("\n", istable, "out of", itot, "that are stable downloaded POSCAR", '\n')

#outESPEI(posdir,hstack,ehull)

if formula_within!="":
  calcphasedecompostion(hstack,formula_within,within,compos_within)
