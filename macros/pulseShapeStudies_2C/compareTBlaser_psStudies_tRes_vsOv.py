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

#def noise_vs_Vov(x, par):
#    xx = x[0]
#    Npe = LO_at3p5* PDE(sipmType, xx) /PDE(sipmType, 3.50))* 4.2 * 1.25
#    noise_single = math.sqrt( pow(par[3]/f_SR_vs_Npe_ref.Eval(Gain((sipmType, xx))),2) + 16.7*16.7 )
#    return noise_single / math.sqrt(2)
#
#def stoch_vs_Vov(x, par):
#        xx = x[0]
#        Npe = LO_at3p5* PDE(sipmType, xx) /PDE(sipmType, 3.50))* 4.2 * 1.25
#        return par[1] * pow(par[0]/Npe,par[2]) 
#
#def tot_vs_Vov(x, par):
#    xx = x[0]
#    return math.sqrt( pow(noise_vs_Vov(x,par),2) + pow(stoch_vs_Vov(x,par),2) )



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
label = 'HPK_nonIrr_C25_LYSO818'
sipmType = 'HPK-PIT-C25-ES2'
LO_at3p5 = 2390
temp = 'T5C'
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/compareLaser/'%(label, temp)
g_Npe_vs_Vov_ref = ROOT.TGraph()
thBestFromTB = {
        0.60: 5,
        0.80: 11,
        1.00: 11,
        1.25: 15,
        1.50: 15,
        2.00: 20,
        3.00: 20,
        3.50: 25,
        }
thFixed = [7, 11, 15, 20, 25 ]
bars = [0,1,2,3,4,5,6,7,8,9,10,11,13,14,15]
Vovs = [ 0.60, 0.80, 1.00, 1.25, 1.50, 2.00]

def tot(sipmType, vov, Npe, fSR):
    noise_single = math.sqrt( pow(457/f_SR.Eval(Npe*Gain(sipmType, vov)),2) + 16.7*16.7 )
    stoch = 25.7 * pow(7000/Npe,0.5)
    return math.sqrt( pow(noise_single/math.sqrt(2),2) + pow(stoch,2) ) 

def totNEW(sipmType, vov, Npe, fSR):
    noise_single = math.sqrt( pow(457/f_SR.Eval(Npe*Gain(sipmType, vov)),2) + 16.7*16.7 )
    stoch = (35.8 - 1.1E-05 *Gain(sipmType, vov) )* pow(7000/Npe, 0.77)
    return math.sqrt( pow(noise_single/math.sqrt(2),2) + pow(stoch,2) ) 



inFileTB = ROOT.TFile.Open('/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/HPK_nonIrr_C25_LYSO818_T5C/summaryPlots_HPK_nonIrr_C25_LYSO818_T5C.root')
inFileLaser = ROOT.TFile.Open('plots_SRt1_LYSO818_pulseShapeStudy_vsNpe.root')
inFileLaserTRes = ROOT.TFile.Open('plots_tRes_LYSO818_pulseShapeStudy_vsNpe.root')
#tRes_vs_Vov = inFileTB.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin01_average')
print (inFileLaser)

tRes_Laser = ROOT.TGraphErrors()
tRes_exp = ROOT.TGraphErrors()
tRes_expUp = ROOT.TGraphErrors()
tRes_expDown = ROOT.TGraphErrors()

tRes_exp2 = ROOT.TGraphErrors()
tRes_expUp2 = ROOT.TGraphErrors()
tRes_expDown2 = ROOT.TGraphErrors()

angleFact = math.cos(49 * math.pi / 180) / math.cos(52*math.pi / 180)
thRef = -999

if thRef > 0:

    tRes_vs_Vov = inFileTB.Get('g_deltaT_totRatioCorr_vs_vov_bar07_th%02d_enBin01'%thRef)
    tRes_vs_Vov.Scale(1./angleFact)
    for vov in Vovs:
        print (vov, thBestFromTB[vov])
        
        f_SR = inFileLaser.Get('f_SRglob_vs_gainNpe_th%02d'%thRef)
        #print ('f_SRglob_vs_gainNpe_th%02d'%thRef)
        Npe = LO_at3p5* (PDE(sipmType, vov) /PDE(sipmType, 3.50))* (4.2) * 1.25
        gLaserRes = inFileLaserTRes.Get('tRes_vs_Npe_Vov%0.2f_th%02d'%(vov, thRef)) 
        if gLaserRes != None: 
        
            #gLaserRes.Print()
            for point in range(0, gLaserRes.GetN()):
                #print (Npe, gLaserRes.GetPointX(point))
                if gLaserRes.GetPointX(point) >= Npe - Npe*0.20 : 
                    print (Npe, gLaserRes.GetPointX(point), gLaserRes.GetPointY(point))
                    tRes_Laser.SetPoint(tRes_Laser.GetN(), vov, gLaserRes.GetPointY(point))    
                    break
        tRes_exp.SetPoint (tRes_exp.GetN(), vov, tot(sipmType, vov, Npe, f_SR))
        tRes_expUp.SetPoint (tRes_expUp.GetN(), vov, tot(sipmType, vov, Npe*1.1, f_SR))
        tRes_expDown.SetPoint (tRes_expDown.GetN(), vov, tot(sipmType, vov, Npe*0.9, f_SR))
    
    c = ROOT.TCanvas("c_tRes_vs_Vov", "c_tRes_vs_Vov")
    c.cd()
    c.SetGridy()
    hPad = ROOT.gPad.DrawFrame(0,0.,4.0, 100.)
    hPad.SetTitle(";  V_{ov}; time resolution [ps]")
    
    leg = ROOT.TLegend(0.8, 0.5, 0.95 , 0.95)
    
    
    tRes_vs_Vov.SetLineColor(ROOT.kBlack)
    tRes_vs_Vov.SetMarkerColor(ROOT.kBlack)
    tRes_vs_Vov.Draw("PL")
    
    tRes_Laser.SetLineColor(ROOT.kRed)
    tRes_Laser.SetMarkerColor(ROOT.kRed)
    tRes_Laser.SetMarkerStyle(21)
    tRes_Laser.SetLineWidth(2)
    tRes_Laser.Draw("PLsame")
    
    
    
    tRes_exp.SetLineColor(ROOT.kPink+1)
    tRes_exp.SetLineWidth(2)
    tRes_exp.Draw("Lsame")
    
    tRes_expUp.SetLineColor(ROOT.kPink+1)
    tRes_expUp.SetLineStyle(7)
    tRes_expUp.SetLineWidth(2)
    tRes_expUp.Draw("Lsame")
    
    tRes_expDown.SetLineColor(ROOT.kPink+1)
    tRes_expDown.SetLineStyle(7)
    tRes_expDown.SetLineWidth(2)
    tRes_expDown.Draw("Lsame")
    
    
    
    #leg.AddEntry(g, "th = %02d"%th, "PL" )
    #leg.Draw("same")
    c.SaveAs(outdir+c.GetName()+'thRef%02d.png'%thRef)
        

else:

    tRes_vs_Vov = inFileTB.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_bar07_enBin01')
    tRes_vs_Vov.Print()
    tRes_vs_Vov.Scale(1./angleFact)
    tRes_vs_Vov.Print()

    for vov in Vovs:
        print (vov, thBestFromTB[vov])
        
        f_SR = inFileLaser.Get('f_SRglob_vs_gainNpe_linlog_th%02d'%thBestFromTB[vov])
        #f_SR = inFileLaser.Get('g_SR_vs_gainNpe_Vov%.2f_th%02d'%(vov, thBestFromTB[vov]))
        Npe = LO_at3p5* (PDE(sipmType, vov) /PDE(sipmType, 3.50))* (4.2) * 1.25
        gLaserRes = inFileLaserTRes.Get('tRes_vs_Npe_Vov%0.2f_th%02d'%(vov, thBestFromTB[vov])) 
        if gLaserRes != None: tRes_Laser.SetPoint(tRes_Laser.GetN(), vov, gLaserRes.GetPointY(gLaserRes.GetN()-1))    
        tRes_exp.SetPoint (tRes_exp.GetN(), vov, tot(sipmType, vov, Npe, f_SR))
        tRes_expUp.SetPoint (tRes_expUp.GetN(), vov, tot(sipmType, vov, Npe*1.1, f_SR))
        tRes_expDown.SetPoint (tRes_expDown.GetN(), vov, tot(sipmType, vov, Npe*0.9, f_SR))

        tRes_exp2.SetPoint (tRes_exp2.GetN(), vov, totNEW(sipmType, vov, Npe, f_SR))
        tRes_expUp2.SetPoint (tRes_expUp2.GetN(), vov, totNEW(sipmType, vov, Npe*1.1, f_SR))
        tRes_expDown2.SetPoint (tRes_expDown2.GetN(), vov, totNEW(sipmType, vov, Npe*0.9, f_SR))
 

    c = ROOT.TCanvas("c_tRes_vs_Vov", "c_tRes_vs_Vov")
    c.cd()
    c.SetGridy()
    hPad = ROOT.gPad.DrawFrame(0,0.,4.0, 100.)
    hPad.SetTitle(";  V_{ov}; time resolution [ps]")
    
    leg = ROOT.TLegend(0.8, 0.5, 0.95 , 0.95)
    
    
    tRes_vs_Vov.Draw("PL")
    
    tRes_Laser.SetLineColor(ROOT.kRed)
    tRes_Laser.SetLineWidth(2)
    #tRes_Laser.Draw("Lsame")
    
    
    
    tRes_exp.SetLineColor(ROOT.kPink+1)
    tRes_exp.SetLineWidth(2)
    tRes_exp.Draw("Lsame")
    
    tRes_expUp.SetLineColor(ROOT.kPink+1)
    tRes_expUp.SetLineStyle(7)
    tRes_expUp.SetLineWidth(2)
    tRes_expUp.Draw("Lsame")
    
    tRes_expDown.SetLineColor(ROOT.kPink+1)
    tRes_expDown.SetLineStyle(7)
    tRes_expDown.SetLineWidth(2)
    tRes_expDown.Draw("Lsame")
    
    tRes_exp2.SetLineColor(ROOT.kCyan)
    tRes_exp2.SetLineWidth(2)
    tRes_exp2.Draw("Lsame")
    
    tRes_expUp2.SetLineColor(ROOT.kCyan)
    tRes_expUp2.SetLineStyle(7)
    tRes_expUp2.SetLineWidth(2)
    tRes_expUp2.Draw("Lsame")
    
    tRes_expDown2.SetLineColor(ROOT.kCyan)
    tRes_expDown2.SetLineStyle(7)
    tRes_expDown2.SetLineWidth(2)
    tRes_expDown2.Draw("Lsame")
    
    
    #leg.AddEntry(g, "th = %02d"%th, "PL" )
    #leg.Draw("same")
    c.SaveAs(outdir+c.GetName()+'_bestTh.png')
        

