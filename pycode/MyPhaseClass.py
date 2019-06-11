import  numpy as np
import sys
import os
import json
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import shutil


class MyPhaseClass:
  def __init__(self):
    self.sphases,self.sratio,self.sspacegroup = ShunliExcelPhaseName()
    self.aphases,self.aratio,self.aspacegroupnumber, self.aspacegroupsymbol, \
      self.astrucdesign, self.aprototype = AFLOWExcelPhaseName()
    
  def KustPhaseNameByFileName(self, filename):
    xx = filename.split('+')
    pname = self.KustPhaseName(xx[0], xx[1])
    if pname==None:
      return 't-'+filename
    else:
      return pname+'+'+xx[-1]

  def KustPhaseName(self, formula, spacegroup, siteoc=[]):
    if len(siteoc)==0 : siteoc = formula2siteoc(formula)
    try:
      spacegroup = int(spacegroup)
    except:
      for i,s in enumerate(self.aspacegroupsymbol):
        if spacegroup==s:
          spacegroup = self.aspacegroupnumber[i]

  #print(siteoc)
    for i,s in enumerate(self.sspacegroup):
      if spacegroup == s:
        if siteoc == self.sratio[i]:
          return self.sphases[i]+'+'+formula

    for i,s in enumerate(self.aspacegroupnumber):
      if spacegroup == s:
        if siteoc == self.aratio[i]:
          if self.astrucdesign[i]!="":
            try:
              return 'af-'+self.astrucdesign[i].upper()+'_'+str(s)+'+'+formula
              #return self.astrucdesign[i].upper()+'_'+str(s)+'+'+formula
            except:
              return 'p-'+self.aprototype[i].upper()+'_'+str(s)+'+'+formula
              #return 'af-'+self.aphases[i]+'+'+formula
              #return 'af-'+self.aphases[i]+'+'+formula
          else:
            return 'af-'+self.aphases[i]+'+'+formula



def ShunliExcelPhaseName():
  #
  df = pd.read_excel('2019-02-08-List-Phases.xlsx', sheet_name='model_MP')
  #df = pd.read_excel('2019-02-22-List-Phases.xlsx', sheet_name='model_MP')
  columns = df.columns.values
  for cc in columns:
    if "Phases".upper() == cc.strip().upper() : phases = cc
    elif "Ratio".upper() == cc.strip().upper() : ratio = cc
    elif "Spacegroup".upper() == cc.strip().upper() : spacegroup = cc
  phases = df[phases]
  sratio = []
  for x in df[ratio]:
    s = []
    for xx in x.replace('[','').replace(']','').replace(' ','').split(','):
      s.append(int(xx))
    sratio.append(sorted(s))
  ss = []
  for i,x in enumerate(df[spacegroup]):
    try:
      ss.append(int(x))
    except:
      ss.append(-1)
  spacegroup = ss
  for i,x in enumerate(phases):
    try:
      xx = x.split('_')
      if len(xx) > 1: spacegroup[i] = int(xx[-1])
    except:
      pass
  return phases,sratio,spacegroup


def AFLOWExcelPhaseName():
  df = pd.read_excel('2019-02-08-List-Phases.xlsx', sheet_name='AFLOW')
  #df = pd.read_excel('2019-02-22-List-Phases.xlsx', sheet_name='model_MP')
  columns = df.columns.values
  for cc in columns:
    if "AFLOW Prototype".upper() == cc.strip().upper() : phases = cc
    elif "Prototype".upper() == cc.strip().upper() : prototype = cc
    elif "Space Group Number".upper() == cc.strip().upper() : spacegroup = cc
    elif "Space Group Symbol".upper() == cc.strip().upper() : spacegroupS = cc
    elif "Strukturbericht Designation".upper() == cc.strip().upper() : strucdesign = cc
  phases = df[phases]
  prototype = df[prototype] 
  ratio = [formula2siteoc(x) for x in prototype]
  spacegroupnumber = df[spacegroup]
  #spacegroupnumber = [int(x) for x in df[spacegroup]]
  spacegroupsymbol = [x.replace('/','') for x in df[spacegroupS]]
  return phases,ratio,spacegroupnumber, spacegroupsymbol,df[strucdesign],prototype


def formula2siteoc(formula):
  formula = formula.replace(" ",'').replace("-",'').replace(",",'').replace("_",'')
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
    ele.append(newel)

    if (len(newcc)!=0):
      try:
        com.append(int(newcc))
      except:
        print('"',newcc,'" is not a float number! your formula is wrong!')
        sys.exit(1)
    else:
      com.append(1)

  #sorted the sequence and merge the duplicate
  nat = sorted(com)
  Nd = min(nat)
  for i in range(Nd,0,-1):
    out = True
    for j in range(len(nat)):
      if ((nat[j]//i)*i!=nat[j]):
        out = False
        break
    if out:
      break
  for j,n in enumerate(nat):
    nat[j] //= i
  return nat

"""
  test call

KUST0 = "KUST0_all_compd_2019_2_25/"
KUST1 = "kustdir/"
myphaseclass= MyPhaseClass()
if not os.path.exists(KUST1):
  os.mkdir(KUST1)

for x in os.listdir(KUST0):
  xx = x.split('+')
  pname = myphaseclass.KustPhaseNameByFileName(x)
  if pname == None: pname=x
  print (pname)
  shutil.copy2(KUST0+x,KUST1+pname)
"""
