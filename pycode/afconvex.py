import  numpy as np
from scipy.optimize import linprog
from numpy.linalg import solve
from fractions import Fraction
import traceback
import signal
import sys
import os
import shutil
import time
import json
import urllib
import requests
from urllib.request import urlopen,Request,urlretrieve
from aflow import *

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

vaspPot = {"H" :  "H" ,
"He" :  "He" ,
"Li" :  "Li_sv" ,
"Be" :  "Be" ,
"B" :  "B" ,
"C" :  "C" ,
"N" :  "N" ,
"O" :  "O" ,
"F" :  "F" ,
"Ne" :  "Ne" ,
"Na" :  "Na_pv" ,
"Mg" :  "Mg" ,
"Al" :  "Al" ,
"Si" :  "Si" ,
"P" :  "P" ,
"S" :  "S" ,
"Cl" :  "Cl" ,
"Ar" :  "Ar" ,
"K" :  "K_sv" ,
"Ca" :  "Ca_sv" ,
"Sc" :  "Sc_sv" ,
"Ti" :  "Ti_sv" ,
"V" :  "V_sv" ,
"Cr" :  "Cr_pv" ,
"Mn" :  "Mn_pv" ,
"Fe" :  "Fe" ,
"Co" :  "Co" ,
"Ni" :  "Ni" ,
"Cu" :  "Cu" ,
"Zn" :  "Zn" ,
"Ga" :  "Ga_d" ,
"Ge" :  "Ge_d" ,
"As" :  "As" ,
"Se" :  "Se" ,
"Br" :  "Br" ,
"Kr" :  "Kr" ,
"Rb" :  "Rb_sv" ,
"Sr" :  "Sr_sv" ,
"Y" :  "Y_sv" ,
"Zr" :  "Zr_sv" ,
"Nb" :  "Nb_sv" ,
"Mo" :  "Mo_pv" ,
"Tc" :  "Tc_pv" ,
"Ru" :  "Ru_pv" ,
"Rh" :  "Rh_pv" ,
"Pd" :  "Pd" ,
"Ag" :  "Ag" ,
"Cd" :  "Cd" ,
"In" :  "In_d" ,
"Sn" :  "Sn_d" ,
"Sb" :  "Sb" ,
"Te" :  "Te" ,
"I" :  "I" ,
"Xe" :  "Xe" ,
"Cs" :  "Cs_sv" ,
"Ba" :  "Ba_sv" ,
"La" :  "La" ,
"Ce" :  "Ce" ,
"Pr" :  "Pr_3" ,
"Nd" :  "Nd_3" ,
"Pm" :  "Pm_3" ,
"Sm" :  "Sm_3" ,
"Eu" :  "Eu_2" ,
"Gd" :  "Gd_3" ,
"Tb" :  "Tb_3" ,
"Dy" :  "Dy_3" ,
"Ho" :  "Ho_3" ,
"Er" :  "Er_3" ,
"Tm" :  "Tm_3" ,
"Yb" :  "Yb_2" ,
"Lu" :  "Lu_3" ,
"Hf" :  "Hf_pv" ,
"Ta" :  "Ta_pv" ,
"W" :  "W_pv" ,
"Re" :  "Re" ,
"Os" :  "Os" ,
"Ir" :  "Ir" ,
"Pt" :  "Pt" ,
"Au" :  "Au" ,
"Hg" :  "Hg" ,
"Tl" :  "Tl_d" ,
"Pb" :  "Pb_d" ,
"Bi" :  "Bi_d" ,
"Po" :  "Po_d" ,
"At" :  "At_d" ,
"Rn" :  "Rn" ,
"Fr" :  "Fr_sv" ,
"Ra" :  "Ra_sv" ,
"Ac" :  "Ac" ,
"Th" :  "Th" ,
"Pa" :  "Pa" ,
"U" :  "U" ,
"Np" :  "Np" ,
"Pu" :  "Pu" ,
"Am" :  "Am" ,
"Cm" :  "Cm" }

def aflowPotSet():
  aflowPot = {}
  for i,rec in enumerate(data):
    if not rec.get("code").startswith("vasp") : continue
    if rec.get("dft_type")!="PAW_PBE" : continue
    elp = rec.get("species_pp").split(",")

    for pot in elp:
      if pot in aflowPot.keys():
        aflowPot[pot] += 1
      else:
        aflowPot[pot] = 1
  for key in vaspPot:
    kk = -1
    for k1 in aflowPot:
      el = k1.split('_')[0]
      if el==key:
        if aflowPot[k1]>kk:
          pot = k1
          kk = aflowPot[k1]
    if kk > 0: vaspPot[key] = pot

def signal_handler(sig, frame):
        print('You pressed Ctrl+C!')
        sys.exit(0)

def aflow_missing(entry):
  if os.path.exists("aflow_missing.json") :
    with open("aflow_missing.json", encoding='utf-8') as json_file:
      missing_data = json.load(json_file)
      missing_data.append(str(entry))
  else:
      missing_data = [str(entry)]
  with open('aflow_missing.json', 'w') as outfile:
    json.dump(missing_data, outfile)
  sys.exit(0)

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
def formula2nat(formula):
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
        com.append(int(newcc))
      else:
        print('"',newcc,'" is not a float number! your formula is wrong!')
        sys.exit(1)
    else:
      com.append(1)
  com = np.array(list(map(int,com)))

  #sorted the sequence and merge the duplicate
  elist = sorted(set(ele))
  clist = np.zeros(len(elist), dtype=int)
  for j,el in enumerate(ele):
    ix = elist.index(el)
    clist[ix] += com[j]

  return elist,clist

def prety_formula(_els,_nat):
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

def downloads(filters):
  #http://aflowlib.duke.edu/search/API/?species((Hf,Pt),Rh),$nspecies(3),energy_atom,$paging(0)
  #{"compound":"Hf1Pt1Rh2","auid":"aflow:bb46bd634476e7fb","aurl":"aflowlib.duke.edu:AFLOWDATA/LIB3_RAW/Hf_pvPtRh_pv/TBCC015.CAB","species":"Hf,Pt,Rh"}
  #   ).orderby(K.energy_atom)
  sys.stdout.write('\nScreening {} ... '.format(filters))
  result = search(batch_size=9999
     ).filter(filters
     ).filter((K.code % 'vasp')
     ).select(K.energy_atom, K.species, K.composition, K.compound, K.code, K.dft_type, K.volume_atom, K.Pearson_symbol_relax, K.prototype, K.sg2, K.spacegroup_relax, K.species_pp)
  return processresult(result)

def processresult(result):
  global allrec, nentries, totaltime, startime, missing_data
  allcompds = {}
  try:
    sys.stdout.write('{}+{} entries returned.\n'.format(nentries,len(result)))
    nentries += len(result)
  except:
    sys.stdout.write('{}+{} entries returned.\n'.format(nentries,0))
    return allcompds

  for entry in result:
    if str(entry) in missing_data: 
      print(str(entry), "not FOUND!")
      continue
    try:
      if not all(elem in runPot for elem in entry.species_pp) : continue
      composition = entry.composition
      compound = prety_formula(entry.species,composition)
      if compound in allcompds.keys():
        energy_atom = entry.energy_atom
        if float(energy_atom) >= float(allcompds.get(compound)[0]) : continue
      """
      ret = requests.head(str(entry))
      if ret.status_code >= 400 : continue #error in AFLOW, link not exist
      """

      _dft_type = entry.dft_type[0].strip()
      if _dft_type!=dft_type : continue
      # prec = [entry.compound,entry.code,entry.dft_type,entry.enthalpy_formation_atom,entry.natoms,entry.nspecies,entry.Pearson_symbol_orig,entry.prototype,entry.sg,entry.spacegroup_orig,entry.species,entry.species_pp,entry.stoichiometry,entry.energy_atom,entry.volume_atom]
      arec = [entry.energy_atom, entry, entry.species, composition, entry.compound, entry.code, _dft_type, entry.volume_atom, entry.Pearson_symbol_relax, entry.prototype, entry.sg2, entry.spacegroup_relax, entry.auid]
      prec = {compound:arec}
      #print(prec)
      allcompds.update(prec)
    except:
      traceback.print_exc()
      pass
  sys.stdout.write ('{} + {} Secs. '.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime
  if len(allcompds) != 0:
    sys.stdout.write('{}+{} compounds accumulated\n'.format(len(allrec),len(allcompds)))
  return allcompds
            

def outelement(els):
  if containoperator=="$all":
    if all(elem in els for elem in contain) : return True
    else : return False
  elif containoperator=="$in":
    for el in els:
      if el in contain: return True
    return False
  else:
    return True


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
  print('\n\n#Elemental components:', Elements)

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
    nComponents /= nComponents.sum()

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
  hstack = []
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
    #print(sum(res.x*res.x))
    #print(dstack[j])
    #if (abs(sum(res.x*res.x)-1.0) < 1.e-3) and abs(res.fun-Gstack[j])<ehullthr:
    if abs(res.fun-Gstack[j])<ehullthr:
      """stable phase"""
      dstack[j]["hull"] = ""
      hstack.append(dstack[j])
      #print(dstack[j])
      ndata += 1
    else:
      """unstable phase, give the phase separation"""
      dstack[j]["hull"] = hull
    """energy above the hull"""
    dstack[j]["hull_energy"] = (Gstack[j]-res.fun)
    #print(j,"/",len(amatrix), " Phases: ", Phases[j])
    if allstr :
      sys.stdout.write("{:5d}/{} : ".format(j,len(amatrix)))
      ff = dstack[j]
      try:
        sys.stdout.write("{:<16s} {:3d} {:12s} {:7s} {:<16s} E0={:11.6f} above hull {:.6f}\n".format(ff["phasename"], ff["spacegroup"], ff["sg2"].split(' ')[0], ff["Pearson_symbol"], ff["prototype"],ff["energy"], ff["hull_energy"]))
      except:
        try:
          webget(ff["entry"],[],[],ff)
          sys.stdout.write(" E0={:11.6f} above hull {:.6f}\n".format(ff["energy"], ff["hull_energy"]))
        except:
          sys.stdout.write("{:<16s} {:3d} {:12s} {:7s} {:<16s} E0={:11.6f} {:.6f}\n".format(ff["phasename"], spacegroup, sg2, Pearson_symbol, prototype,ff["energy"], ff["hull_energy"]))

  print("\n", ndata, " Compounds made the convex hull for the system\n")
  return hstack,dstack

def webget(entry,contcars=[],kpoints=[],ff={}):
  #req = urllib.request.Request(entry)
  for i in range(16) :
    try:
      #url = urlopen(req,timeout = 1)
      page = urllib.request.urlopen(entry,timeout=16)
      records = page.readlines()
      break
    except:
      if i==15 :
        #print("CANNOT open", entry)
    #raise ValueError('A very specific bad thing happened.')
        sys.stdout.write("{:<16s}".format(ff["phasename"]))
        return entry
  
  con = False
  kpo = False
  poscarfile = ""
  
  for ii,bline in enumerate(records) :
    aline = bline.decode("utf-8").strip()
    line = aline.split("</a>=")
    if line[0].endswith("Pearson_symbol_relax") : Pearson_symbol = line[1]
    elif line[0].endswith("prototype") : prototype = line[1]
    elif line[0].endswith("sg2") : sg2 = line[1].split(" ")[0]
    elif line[0].endswith("spacegroup_relax") : spacegroup = int(line[1])
    elif line[0].endswith("energy_cutoff") : energy_cutoff = line[1]
    elif line[0].endswith("species_pp") : species_pp = line[1]
  spname = ff["phasename"]+'+'+sg2.replace('{','').replace('}','').replace('/','').replace('_','')+'+'+prototype+'+'+str(spacegroup)+'+'+ff["afid"]
     
  for ii,bline in enumerate(records) :
    aline = bline.decode("utf-8").strip()
    line = aline.split("</a>=")
      
    if len(kpoints)!=0 and not kpo:
      for kpoint in kpoints:
        try:
          aline.index(kpoint+'</a>')
          kpo = True
        except:
          pass
        if kpo:
          urllib.request.urlretrieve(entry+'/'+kpoint, posdir + spname + ".KPOINTS")
          break

    if len(contcars)!=0 and not con:
      for contcar in contcars:
        try:
          aline.index(contcar+'</a>')
          con = True
        except:
          pass
        if con:
          poscarfile = posdir + spname + ".VASP"
          urllib.request.urlretrieve(entry+'/'+contcar, poscarfile)
          alines = open(poscarfile, 'r').readlines()
          alines.insert(5,"   "+ff["elements"]+'\n')
          N = 8+sum(ff["composition"])
          open(poscarfile, 'w').writelines(alines[:N])
          kustoutfiles(poscarfile)
          break

  open(posdir + spname + ".POTCAR", 'w').write(str(species_pp)+'\n')
  open(posdir + spname + ".ENCUT", 'w').write(str(energy_cutoff)+'\n')
  if ff["hull_energy"]<=ehullthr:
    sys.stdout.write("{:<16s} {:<16s} {:3d} {:12s} {:7s} {:<22s}".format(ff["phasename"], ff["pretty_formula"], spacegroup, sg2, Pearson_symbol, prototype))
  else:
    sys.stdout.write("{:>16s} {:>16s} {:3d} {:12s} {:7s} {:<22s}".format(ff["phasename"], ff["pretty_formula"], spacegroup, sg2, Pearson_symbol, prototype))
  return poscarfile
  

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

  prettycom = ""
  print("\nBy atomic percentages, the phases are:\n")
  for x in X:
    #phases = Phases[x[0]].split("#")
    phases = Phases[x[0]]
    print('{:6.2f}  {:12s}'.format(x[1]*100, phases))
    if prettycom == "":
      prettycom = '{:6.2f}*{}'.format(x[1]*100, phases)
    else:
      prettycom += ' +' + '{:6.2f}*{}'.format(x[1]*100, phases)
      
  print("\n", formula, " is made of ", prettycom, "\n")

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


"""output POSCAR, INCAR, KPOINTS and stability of each download structure
posdir - dir for files to be outputted into
dstack - data stack contain all information, such as phasename, energy, elements, etc
calhull - bool to instruct if calculatiing convex hull
ehull - how high above the hull the structures will be included in the output
return:
"""
def outMPdata(posdir, _dstack, calhull, ehull):
  KPOINTS = ['KPOINTS.relax2', 'KPOINTS.relax1', 'KPOINTS.relax', 'KPOINTS.static']
  CONTCARS = ['CONTCAR.relax2', 'CONTCAR.relax1', 'CONTCAR.relax', 'CONTCAR.relax.vasp', 'POSCAR.bands', 'POSCAR.orig']
  print ("Compound         Formula           SG Symmetry     Pearson Prototype         Tenergy(eV/atom)    Fenergy : Primitive POSCAR")
  dstack = []
  for ff in _dstack:
    if ff["hull_energy"] < 1.e-6:
      dstack.append(ff)
  for ff in _dstack:
    if not ff["hull_energy"] < 1.e-6:
      dstack.append(ff)

  for ff in dstack:
    if not outelement(ff["elements"]) : continue
    pname = ff["phasename"]
    energy = ff["energy"]
    entry = ff["entry"]

    #sys.stdout.write ('{} Secs. '.format(round(time.time()-startime,3)))
    if (not calhull):
      sys.stdout.write("{} {} {} {} {:.6f}\n".format(ff["pretty_formula"], ff["elements"], ff["composition"], pname, energy))
    else:
      eh = ff["hull_energy"]
      ef = ff["formation_energy"]
      if eh < ehull:
        if eh < ehullthr:
          if preload:
            try:
              poscarfile = webget(ff["entry"],CONTCARS,KPOINTS,ff)
              sys.stdout.write(" {:11.6f} {:10.6f} : ".format(energy,ef))
            except:
              traceback.print_exc()
              sys.stdout.write("{:<16s} {:11.6f} {:10.6f} : ".format(pname, energy, ef))
            if eh > 0.0000005: sys.stdout.write("{:.6f} above hull {}, ".format(eh,ff["hull"]))
            sys.stdout.write("{}\n".format(poscarfile))
            continue
          else:
            sys.stdout.write("{:<16s} {:3d} {:12s} {:7s} {:<22s} {:11.6f} {:10.6f} : ".format(pname, ff["spacegroup"], ff["sg2"].split(' ')[0], ff["Pearson_symbol"], ff["prototype"],energy, ef))
            if eh > 0.0000005: sys.stdout.write("{:.6f} above hull {}, ".format(eh,ff["hull"]))
        else:
          poscarfile = ""
          if preload:
            try:
              poscarfile = webget(ff["entry"],CONTCARS,KPOINTS,ff)
              sys.stdout.write(" {:11.6f} {:10.6f} {:.6f} above hull : ".format(energy,ef,eh))
            except:
              sys.stdout.write("{:<16s} {:11.6f} {:10.6f} {:.6f} above hull : ".format(pname, energy, ef, eh))
            if eh > 0.0000005: sys.stdout.write("{:.6f} above hull {}, ".format(eh,ff["hull"]))
            sys.stdout.write("{}\n".format(poscarfile))
            continue
          else:
            sys.stdout.write("{:<16s} {:3d} {:12s} {:7s} {:<22s} {:11.6f} {:10.6f} {:.6f} above hull : ".format(pname, ff["spacegroup"], ff["sg2"].split(' ')[0], ff["Pearson_symbol"], ff["prototype"],energy, ef, eh))
            if eh > 0.0000005: sys.stdout.write("{:.6f} above hull {}, ".format(eh,ff["hull"]))
      else:
        continue

    spname = ff["phasename"]+'+'+ff["sg2"].split(' ')[0].replace('{','').replace('}','').replace('/','').replace('_','')+'+'+ff["prototype"]+'+'+str(ff["spacegroup"])+'+'+ff["afid"]
    for cc in KPOINTS:
      try:
        entry.files[cc](posdir + spname + ".KPOINTS")
      except:
        """ cannot find CONTCAR cc from entry, goto web download """
        kpoints = str(entry)+'/'+cc
        ret = requests.head(kpoints)
        if ret.status_code >= 400 : continue
        urllib.request.urlretrieve(kpoints, posdir + spname + ".KPOINTS")
      break

    found = False
    for cc in CONTCARS:
      try:
        lines = entry.files[cc]()
        lines = lines.split('\n')
        lines = [k for k in lines if k != '']
      except:
        """ cannot find CONTCAR cc from entry, goto web download """
        contcar = str(entry)+'/'+cc
        ret = requests.head(contcar)
        if ret.status_code >= 400 : continue
        req = urllib.request.Request(contcar)
        blines = urlopen(req).readlines()
        lines = []
        for bb,bline in enumerate(blines):
          lines.append(bline.decode("utf-8").strip())

      with open(posdir + spname + ".VASP", 'w') as f0:
        for l in range(0,5):
          f0.write("{}\n".format(lines[l]))
        f0.write("   {}\n".format(ff["elements"]))
        natom = sum(ff["composition"])
        for l in range(5,7+natom):
          f0.write("{}\n".format(lines[l]))
        print(posdir + spname + ".VASP")
        found = True
      kustoutfiles(posdir + spname + ".VASP")
      break

    if not found :
      print("WARNING! CONTCAR not FOUND", str(entry))
      #aflow_missing(entry)

    """
    with open(pname+".POTCAR",'w') as f0:
      f0.write("{}\n".format(ff["POTCAR"]))

    with open(pname+".INCAR",'w') as f0:
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
    """
  print()

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

startime = time.time()
totaltime = 0

"""unit conversion"""
eVtoGPa = 160.21766208
eVtoJ = 96486.9
preload = False
within = ['Pt', 'Ni', 'Cr']
within = ['Co', 'Yb', 'Mn', 'Sb']
within = ['Pt', 'Ir', 'Rh', 'Ni', 'Zr', 'Hf', 'Si', 'Cr']
contain = []
"""handle command line option"""

fastcode = True #""" key to cotrol the search method"""
#periodictable = [key for key,value in MM_of_Elements.items()] #""" list of all elements from the periodic table"""
periodictable = MM_of_Elements.keys() #""" list of all elements from the periodic table"""
ehullthr = 1.e-4
ehull = ehullthr
calhull = False #""" key to cotrol the if calculating the convex hull"""
suspend = []

input_within = False #""" key to cotrol the within input"""
input_wCom = False #""" key to cotrol the within input"""
input_contain = False #""" key to cotrol the contain input"""
formula_within = "" #"""chemical formula"""
formula_wCom = "" #"""chemical formula"""
formula_contain = "" #"""chemical formula"""
containoperator = "$all" #contain any elements within the list of contain
dft_type = "PAW_PBE"
quaternary = False
allpot = False
vasppot = False
aflowpot = True
deltaH = False
allstr = False
KUST0 = "SSS/"

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
  elif (sys.argv[count] == "-preload"):
    preload = not preload
  elif (sys.argv[count] == "-kust"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    KUST0 = sys.argv[count].replace('/','') + '/'
  elif (sys.argv[count] == "-allstr"):
    allstr = not allstr
  elif (sys.argv[count] == "-allpot"):
    allpot = not allpot
  elif (sys.argv[count] == "-vasppot"):
    vasppot = not vasppot
  elif (sys.argv[count] == "-aflowpot"):
    aflowpot = not aflowpot
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
  elif (sys.argv[count] == "-ehullthr"):
    count = count + 1
    if (count >= len(sys.argv)):
      break
    ehullthr = float(sys.argv[count])
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

signal.signal(signal.SIGINT, signal_handler)
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

posdir = setoutdir("AF0")
if formula_within=="":
  for el in within:
    formula_within += el
if formula_wCom!="":
  within,compos_within = wCom2composition(formula_wCom)
if formula_within!="":
  within,compos_within = formula2composition(formula_within)
if formula_contain!="":
  contain,compos_contain = formula2composition(formula_contain)

nentries = 0
allrec = {}

missing_data = []
if os.path.exists("aflow_missing.json") :
  with open("aflow_missing.json", encoding='utf-8') as json_file:
    missing_data = json.load(json_file)

if aflowpot:
  """
  with open("download_code.json", encoding='utf-8') as json_file:
    code_data = json.load(json_file)
  with open("download_dft_type.json", encoding='utf-8') as json_file:
    dft_data = json.load(json_file)
  with open("download_pp.json", encoding='utf-8') as json_file:
    pp_data = json.load(json_file)
  with open("download_efa.json", encoding='utf-8') as json_file:
    efa_data = json.load(json_file)
  for i,rec in enumerate(code_data):
    if rec.get("aurl")!=dft_data[i].get("aurl") :
      print ("ERROR!",rec.get("aurl"),"not equal to", dft_data[i].get("aurl"))
      sys.exit()
    code_data[i]["dft_type"] = dft_data[i].get("dft_type")
  for i,rec in enumerate(pp_data):
    if rec.get("aurl")!=dft_data[i].get("aurl") :
      print ("ERROR!",rec.get("aurl"),"not equal to", dft_data[i].get("aurl"))
      sys.exit()
    code_data[i]["species_pp"] = rec.get("species_pp")
  for i,rec in enumerate(efa_data):
    if rec.get("aurl")!=dft_data[i].get("aurl") :
      print ("ERROR!",rec.get("aurl"),"not equal to", dft_data[i].get("aurl"))
      sys.exit()
    code_data[i]["enthalpy_formation_atom"] = rec.get("enthalpy_formation_atom")
  with open("download_spacegroup_relax.json", encoding='utf-8') as json_file:
    efa_data = json.load(json_file)
  for i,rec in enumerate(efa_data):
    if rec.get("aurl")!=dft_data[i].get("aurl") :
      print ("ERROR!",rec.get("aurl"),"not equal to", dft_data[i].get("aurl"))
      sys.exit()
    code_data[i]["spacegroup_relax"] = rec.get("spacegroup_relax")
  with open('aflow.json', 'w') as outfile:
    json.dump(code_data, outfile)
  sys.exit()
  """

  with open("aflow.json", encoding='utf-8') as json_file:
    data = json.load(json_file)
  if not vasppot: aflowPotSet()

runPot = vaspPot.values()

if preload:
  entries = []
  otherPot = []
  for i,rec in enumerate(data):
    if not rec.get("code").startswith("vasp") : continue
    if rec.get("dft_type")!="PAW_PBE" : continue
    spc = rec.get("species").replace("'","").split(",")
    if all(elem in within for elem in spc) :
      if allpot: entries.append(rec)
      else:
        elp = rec.get("species_pp").split(",")
        if all(elem in runPot for elem in elp) :
          entries.append(rec)
        else:
          for el in elp:
            if el not in runPot:
              if el not in otherPot:
                otherPot.append(el)
  if not allpot: print("\n other Potential found", otherPot, "\n")
  nentries = len(entries)
 
  for entry in entries:
    try:
      els,composition = formula2nat(entry.get("compound"))
      prettycompound = prety_formula(els,composition)
      if deltaH:
        energy_atom = entry.get("enthalpy_formation_atom")
        if energy_atom is None: continue
        if not isfloat(energy_atom): continue
        """
        if (len(els)==1) and (float(energy_atom)<0.0) : 
          print(entry)
          continue
        """
      else:
        energy_atom = entry.get("energy_atom")

      #print(energy_atom)
      link = "http://"+entry.get("aurl").replace("duke.edu:","duke.edu/")
      if link in missing_data: continue
      """
      ret = requests.head(link)
      if ret.status_code >= 400 : continue #error in AFLOW, link not exist
      """

      afid = entry["auid"].replace(':','-')
      if ehull <= ehullthr:
        if prettycompound in allrec.keys():
          if float(energy_atom) >= float(allrec.get(prettycompound)[0]) : continue
        arec = [float(energy_atom), link, els, composition, entry.get("compound"),prettycompound,afid]
        allrec.update({prettycompound:arec})
      else:
        compound_spacegroup = entry.get("compound")+"_"+entry.get("spacegroup_relax")
        if compound_spacegroup in allrec.keys():
          if float(energy_atom) >= float(allrec.get(compound_spacegroup)[0]) : continue
        arec = [float(energy_atom), link, els, composition, entry.get("compound"),prettycompound,afid]
        allrec.update({compound_spacegroup:arec})
      """
      else:
        arec = [float(energy_atom), link, els, composition, entry.get("compound"),prettycompound]
        allrec.update({link:arec})
      """
    except:
      traceback.print_exc()
      pass

  if deltaH: 
    for kk in allrec.keys():
      rr = allrec.get(kk)
      if len(rr[2])==1 :
        rr[0] = -0.000000001
        allrec[kk] = rr
  """
  """

  sys.stdout.write ('{} + {} Secs. '.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  if len(allrec) != 0:
    sys.stdout.write('{} compounds found from {} AFLOW entries\n'.format(len(allrec),len(data)))

  dstack = []
  for kk in allrec.keys():
    line = allrec.get(kk)
    drec = {"pretty_formula":line[5], "entry":line[1], "elements":' '.join(line[2]), "phasename":allrec.get(kk)[4], "composition":line[3], "energy":line[0], "afid":afid}
    dstack.append(drec)
  hstack,dstack = convexhull(dstack,fastcode)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  outMPdata(posdir, dstack, True, ehull)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  if formula_within!="":
    calcphasedecompostion(hstack,formula_within,within,compos_within)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime
  print(nentries," entries processed. Total time cost=",round(time.time()-startime,3),"\n")

else:
  for i1,el1 in enumerate(within):
    allrec.update(downloads((K.nspecies == 1) & (K.species==el1)))
    for i2,el2 in enumerate(within):
      if i2 <= i1: continue
      allrec.update(downloads((K.nspecies == 2) & ((K.species==el1) & (K.species==el2))))
      for i3,el3 in enumerate(within):
        if i3 <= i2: continue
        allrec.update(downloads((K.nspecies == 3) & (((K.species==el1) & (K.species==el2)) & (K.species==el3))))
        if not quaternary : continue
        for i4,el4 in enumerate(within):
          if i4 <= i3: continue
          allrec.update(downloads((K.nspecies == 4) & ((((K.species==el1) & (K.species==el2)) & (K.species==el3)) & (K.species==el4))))

  dstack = []
  for kk in allrec.keys():
    line = allrec.get(kk)
    #print (line)
    #drec = {"pretty_formula":kk, "entry":line[1], "elements":' '.join(line[2]), "phasename":allrec.get(kk)[4], "composition":line[3], "energy":line[0], "Pearson_symbol":line[8], "prototype":line[9],"sg2":line[10][0],"spacegroup":line[11]}
    drec = {"pretty_formula":kk, "entry":line[1], "elements":' '.join(line[2]), "phasename":allrec.get(kk)[4], "composition":line[3], "energy":line[0], "Pearson_symbol":line[8], "prototype":line[9],"sg2":line[10][0],"spacegroup":line[11], "afid":line[12].replace(":","-")}
    dstack.append(drec)

  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  hstack,dstack = convexhull(dstack,fastcode)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  outMPdata(posdir, dstack, True, ehull)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime

  if formula_within!="":
    calcphasedecompostion(hstack,formula_within,within,compos_within)
  sys.stdout.write ('{} + {} Secs.\n\n'.format(round(totaltime,3), round(time.time()-startime -totaltime,3)))
  totaltime = time.time()-startime
  print(nentries," entries processed. Total time cost=",round(time.time()-startime,3),"\n")
#http://aflowlib.duke.edu/search/API/?species(!Pu),energy_atom,dft_type,$paging(0)
#http://aflowlib.duke.edu/search/API/?species(!Pu),energy_atom,code,$paging(0)
#http://aflowlib.duke.edu/AFLOWDATA/LIB3_RAW/CoMg_pvSi/TFCC010.ABChttp://aflowlib.duke.edu/AFLOWDATA/LIB3_RAW/CoMg_pvSi/TFCC010.ABC#http://aflowlib.duke.edu/search/API/?species(!Pu),enthalpy_formation_atom,$paging(0)
