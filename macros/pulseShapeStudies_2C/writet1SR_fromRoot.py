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


print (inFileSR)
print ("---------------- SR ------------------------")
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
print ("---------------- t1 and t2 ------------------------")
for thRef in [15, 20, 19, 24, 23, 28]:

    print ("")
    print ("-------- thRef = %d"%thRef)
    ft1 = inFileSR.Get('f_t1glob_vs_gainNpe_th%02d'%thRef)
    for ipar in range(0,3): 
        if ipar == 1: print ("p_%d = %E   +/-   %E"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
        else: print ("p_%d = %f   +/-   %f"%(ipar, ft1.GetParameter(ipar), ft1.GetParError(ipar)))
    print ("")

