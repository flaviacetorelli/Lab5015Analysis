#! /usr/bin/env python3
import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import json

from collections import OrderedDict

import ROOT
import CMS_lumi, tdrstyle                                                                                                                                               

from SiPM import *

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



colors = {
    'L' : ROOT.kRed + 1,
    'R' : ROOT.kBlue + 2,
    'L-R' : ROOT.kBlack,
        }

markers = {
    0.60 : 20, 
    0.80 : 25, 
    1.00 : 22, 
    1.25 : 23, 
    1.50 : 24, 
    2.00 : 21, 
    3.00 : 27, 
    3.50 : 26 
        }
colors = {
    7 : ROOT.kPink +1 , 
    11: ROOT.kOrange, 
    15: ROOT.kRed,  
    20: ROOT.kGreen -4, 
    25: ROOT.kCyan, 
        }
dac_to_uA = {
        1 : 0.313
        }
ithMode = 1
label = 'HPK_2E14_C25_LYSO100056'
sipmType = 'HPK-PIT-C25-ES2'
LO_at3p5 = 2390
temp = 'T-35C'
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/compareLaser/'%(label, temp)
g_Npe_vs_Vov_ref = ROOT.TGraph()

thBestFromTB = { #for bar07
        0.645: 7, 
        0.76:  11,
        0.875: 11,
        0.965: 15,
        }
dcrs = { 
        0.645: 10, #5
        0.76:  13,
        0.875: 20,
        0.965: 25,
        }



#thFixed = [7, 11, 15, 20, 25 ]
bars = [0,1,2,3,4,5,6,7,8,9,10,11,13,14,15]
Vovs = [ 0.645, 0.76, 0.875, 0.965]

def tot(sipmType, vov, Npe, fSR, dcr):
    noise_single = math.sqrt( pow(420/f_SR.Eval(Npe*Gain(sipmType, vov)*0.95),2) + 16.7*16.7 ) #added gain loss 5%
    stoch = 25.7 * pow(7000/Npe,0.5)
    dcr_noise = 34 * (6000 / Npe) * pow(dcr/30, 0.41)
    return math.sqrt( pow(noise_single/math.sqrt(2),2) + pow(stoch,2) + pow(dcr_noise,2)) 

def totNEW(sipmType, vov, Npe, fSR, dcr):
    noise_single = math.sqrt( pow(420/f_SR.Eval(Npe*Gain(sipmType, vov)*0.95),2) + 16.7*16.7 ) #added gain loss 5%
    stoch = 30 * pow(7000/Npe,0.7)
    dcr_noise = 34 * (6000 / Npe) * pow(dcr/30, 0.41)
    return math.sqrt( pow(noise_single/math.sqrt(2),2) + pow(stoch,2) + pow(dcr_noise,2)) 



inFileTB = ROOT.TFile.Open('summaryPlots_HPK_2E14_LYSO100056_angle52_T-35C.root')
inFileLaser = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')
print (inFileLaser)

tRes_Laser = ROOT.TGraphErrors()
tRes_exp = ROOT.TGraphErrors()
tRes_exp2 = ROOT.TGraphErrors()

### scaling for angle factor since in TB Sept offest of 3deg
angleFact = math.cos(49 * math.pi / 180) / math.cos(52*math.pi / 180)


tRes_vs_Vov = inFileTB.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_bar07_enBin01')
tRes_vs_Vov.Scale(1./angleFact)


for vov in Vovs:

    dcr = dcrs[vov]
    thRef = thBestFromTB[vov]
    print (vov, thRef)
    
    f_SR = inFileLaser.Get('f_SRglob_vs_gainNpe_linlog_th%02d'%thRef)

    Npe = LO_at3p5* (PDE(sipmType, vov) * 0.85 /PDE(sipmType, 3.50))* (4.2) * 1.25 #PDE loss of 15% for 2E14
    
    tRes_exp.SetPoint (tRes_exp.GetN(), vov, tot(sipmType, vov, Npe, f_SR, dcr))
    tRes_exp.SetPointError (tRes_exp.GetN()-1, 0., 0.5 * (tot(sipmType, vov, Npe*1.1, f_SR, dcr) - tot(sipmType, vov, Npe*0.9, f_SR, dcr)))

    tRes_exp2.SetPoint (tRes_exp2.GetN(), vov, totNEW(sipmType, vov, Npe, f_SR, dcr))
    tRes_exp2.SetPointError (tRes_exp2.GetN()-1, 0., 0.5 * (totNEW(sipmType, vov, Npe*1.1, f_SR, dcr) - totNEW(sipmType, vov, Npe*0.9, f_SR, dcr)))



c = ROOT.TCanvas("c_tRes_vs_Vov", "c_tRes_vs_Vov")
c.cd()
c.SetGridy()
hPad = ROOT.gPad.DrawFrame(0,0.,2.0, 100.)
hPad.SetTitle(";  V_{ov}; time resolution [ps]")

leg = ROOT.TLegend(0.7, 0.65, 0.92 , 0.92)


tRes_vs_Vov.SetLineColor(ROOT.kBlack)
tRes_vs_Vov.SetMarkerColor(ROOT.kBlack)
tRes_vs_Vov.Draw("PL")

#### draw tRes expected with stoch alpha = 0.5
#tRes_exp.SetLineWidth(2)
#tRes_exp.SetLineColor(ROOT.kPink+1)
#tRes_exp.SetFillColor(ROOT.kPink+1)
#tRes_exp.SetFillColorAlpha(ROOT.kPink+1,0.5)
#tRes_exp.SetFillStyle(3004)
#tRes_exp.Draw('E3lsame')    
#leg.AddEntry(tRes_exp , "alpha =  0.5", "L" )
#leg.AddEntry(tRes_vs_Vov, "TB points", "PL" )
#leg.Draw("same")
#
#
##c.SaveAs(outdir+c.GetName()+'_thBest_alpha05_lowerPDE.png')
#c.SaveAs(outdir+c.GetName()+'_thBest_alpha05.png')
#
##### draw tRes expected with stoch alpha = 0.7
tRes_exp2.SetLineWidth(2)
tRes_exp2.SetLineColor(ROOT.kCyan)
tRes_exp2.SetFillColor(ROOT.kCyan)
tRes_exp2.SetFillColorAlpha(ROOT.kCyan+1,0.5)
tRes_exp2.SetFillStyle(3004)
tRes_exp2.Draw('E3lsame')    

leg.AddEntry(tRes_exp2 , "alpha =  0.7, sigma = 30 ps ", "L" )
#leg.AddEntry(tRes_exp2 , "alpha =  0.7", "L" )
leg.AddEntry(tRes_vs_Vov, "TB points", "PL" )
leg.Draw("same")


#c.SaveAs(outdir+c.GetName()+'_thBest_alpha07.png')
c.SaveAs(outdir+c.GetName()+'_thBest_alpha07_p130.png')
    


