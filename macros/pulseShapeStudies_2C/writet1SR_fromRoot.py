#! /usr/bin/env python3
#macro to print the relevant paramters of the t1, t2, SR and tF fit 
import os
import shutil
import glob
import math
import array
import sys
import time

import ROOT
import tdrstyle

from typing import NamedTuple

from SiPM import *
from ROOT import TMath


inFileSR = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')


print ("---------------- SR ------------------------")
print (inFileSR)
for thRef in [15, 20]:

    print ("")
    print ("-------- thRef = %d"%thRef)
    fSR = inFileSR.Get('f_SRglob_vs_gainNpe_linlog_th%02d'%thRef)
    for ipar in range(0,4): 
        if ipar == 0 or ipar == 2:  print ("p_%d = %E   +/-   %E"%(ipar, fSR.GetParameter(ipar), fSR.GetParError(ipar)))
        else: print ("p_%d = %f   +/-   %f"%(ipar, fSR.GetParameter(ipar), fSR.GetParError(ipar)))

    print ("")

print ("")
print ("")
print ("---------------- t1 ------------------------")

print (inFileSR)
for thRef in [15, 20]:

    print ("")
    print ("-------- thRef = %d"%thRef)
    ft1 = inFileSR.Get('f_t1glob_vs_gainNpe_th%02d'%thRef)
    for ipar in range(0,3): 
        if ipar == 1: print ("p_%d = %E   +/-   %E"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
        else: print ("p_%d = %f   +/-   %f"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
    print ("")


###### t2
print ("---------------- t2 ------------------------")

inFilet2 = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_t2_LYSO818_pulseShapeStudy_vsNpe.root')
print (inFilet2)
for delta in [4,8]:
    print ("Scenario delta = %d"%delta)
    for thRef in [15, 20]:
       print ("")
       print ("--------thRef = %d"%(thRef+delta))
       ft2 = inFilet2.Get('f_t2glob_vs_gainNpe_th%02d'%(thRef+delta))
       for ipar in range(0,3): 
           if ipar == 1: print ("p_%d = %E   +/-   %E"%(ipar, ft2.GetParameter(ipar), ft2.GetParError(ipar)))
           else: print ("p_%d = %f   +/-   %f"%(ipar, ft2.GetParameter(ipar), ft2.GetParError(ipar)))
       print ("")

print ("")
print ("")
print ("---------------- tF ------------------------")


inFiletF = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_tFmint1_LYSO818_pulseShapeStudy_vsNpe.root')
print (inFiletF)
for thRef in [15, 20]:

    print ("")
    print ("-------- thRef = %d"%thRef)
    ft1 = inFiletF.Get('f_tF-t1glob_vs_gainNpe_th%02d'%thRef)
    for ipar in range(0,7): 
        if ipar == 0 or ipar == 2 or ipar == 3 or ipar == 4 or ipar == 5: print ("p_%d = %E   +/-   %E"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
        else: print ("p_%d = %f   +/-   %f"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
    print ("")


