#! /usr/bin/env python3
#macro to fit t2-t1 ad tF-t1 so that they are all defined such as t2 > t1 and tF > t1 in all Npes range
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
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.07,'X')
ROOT.gStyle.SetTitleSize(0.07,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning



colors ={
        15: ROOT.kAzure,
        19: ROOT.kOrange,
        20: ROOT.kGreen +1,
        23: ROOT.kPink,
        24: ROOT.kRed+1,
        28: ROOT.kBlue,
        }
outdir = "/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/laserMeasurements_LYSO818_HPK25/t2Fits_p0FromT1/"
inFileSR = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')

outfilename = '%s/plots_t2_LYSO818_pulseShapeStudy_vsNpe.root'%(outdir)
outfile = ROOT.TFile(outfilename, 'RECREATE')



c = ROOT.TCanvas()
print ("")
print ("")
print ("---------------- t1 and t2 ------------------------")
c.Clear()

deltas = [4,8] # delta from the first th



for delta in deltas:
    print ("---------------------- deltaTh =  %d"%delta)
    for thRef in [15, 20]:
        print ("")
        print ("-------- th1 = %d"%thRef)
    
    
        #getting functions for later initialization
        ft1 = inFileSR.Get('f_t1glob_vs_gainNpe_th%02d'%(thRef))
        ft2 = inFileSR.Get('f_t1glob_vs_gainNpe_th%02d'%(thRef+delta))
        #g_tFglob_vs_gainNpe_th35
       
        #getting t1 and t2 graphs from file (for t2 choose ith1+4 and ith1+8)
        gt1 = inFileSR.Get('g_t1glob_vs_gainNpe_th%02d'%(thRef))
        gt2 = inFileSR.Get('g_t1glob_vs_gainNpe_th%02d'%(thRef+delta))
        
    
        # drawing and fitting tF
    
        c.Clear()
        hPad = ROOT.gPad.DrawFrame(0,0,12E09, 1.5)
        c.SetGridx()
        c.SetGridy()
    
        #dummy = ROOT.TF1 ("dummy_th%d"%(thRef+delta), "[0]*x^[1] ", 0, 12E09)
        fnew = ROOT.TF1 ("fnew_th%d"%(thRef+delta), "[1]*x^[2]+[0] ", 0, 12E09)
        fnew.SetParameters(1, ft2.GetParameter(1), ft2.GetParameter(2))
        fnew.FixParameter(0, ft1.GetParameter(0))
    
        gt2.Fit(fnew, "N", "", 0.E09, 11E09)
         
        c.Clear()
        hPad = ROOT.gPad.DrawFrame(0,0,12E09, 15)
        c.SetGridx()
        c.SetGridy()
    
    
        gt2.GetYaxis().SetTitle("t2 [ns]")
        gt2.Draw("psame")
    
        fnew.SetLineStyle(2)
        fnew.Draw("same")
    
    
        c.SaveAs("%s/t2_th2_%02d.png"%(outdir,thRef+delta))
    
        outfile.cd()
        gt2.Write('g_t2glob_vs_gainNpe_th%02d'%(thRef+delta))
        fnew.Write('f_t2glob_vs_gainNpe_th%02d'%(thRef+delta))

              
