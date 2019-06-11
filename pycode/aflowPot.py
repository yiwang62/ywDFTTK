import  numpy as np
import sys
import os
import json

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

All = vaspPot.keys()
with open("aflow.json", encoding='utf-8') as json_file:
    data = json.load(json_file)
aflowPot = {}
for i,rec in enumerate(data):
    if not rec.get("code").startswith("vasp") : continue
    if rec.get("dft_type")!="PAW_PBE" : continue
    elp = rec.get("species_pp").split(",")

    for pot in elp:
      if pot in aflowPot.keys():
        aflowPot[pot] += 1
      else:
        #aflowPot.update({pot:1})
        aflowPot[pot] = 1
print (json.dumps(vaspPot, indent=4, sort_keys=True))
print (json.dumps(aflowPot, indent=4, sort_keys=True))
for key in vaspPot:
    kk = -1
    for k1 in aflowPot:
      el = k1.split('_')[0]
      if el==key:
        if aflowPot[k1]>kk:
          pot = k1
          kk = aflowPot[k1]
    if kk > 0: vaspPot[key] = pot
    print(kk,pot)

print (json.dumps(vaspPot, indent=4, sort_keys=True))
#print (aflowPot)
