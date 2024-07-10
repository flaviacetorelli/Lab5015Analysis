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
ROOT.gStyle.SetOptFit(0)
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
    3.50 : 26 
        }

dac_to_uA = {
        1 : 0.313
        }
ithMode = 1
label = 'HPK_nonIrr_C25_LYSO818'
sipmType = 'HPK-PIT-C25-ES2'
LO_at3p5 = 2390
temp = 'T5C'
inputdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/'%(label, temp) 
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/compareLaser/'%(label, temp)
#outdir = './'
thFixed = 11
bars = [0,1,2,3,4,5,6,7,8,9,10,11,13,14,15]
Vovs = [0.60, 0.80, 1.00, 1.25, 1.50, 2.00]
#bar = 3
#for bar in bars:
#    inFileTB = ROOT.TFile.Open(inputdir+'c_SR_vs_gainNpe_allVovs_bar%02d.root'%bar)
#    
#    gs = {}
#    for i,Vov in enumerate(Vovs):
#        g = inFileTB.Get('c_SR_vs_gainNpe_allVovs_bar%02d'%bar).FindObject('g_pulseShape%s_bar%02d_Vov%.2f'%('L-R', bar,Vov))
#        g.SetMarkerColor(ROOT.kPink+10)
#        g.SetMarkerStyle(markers[Vov])
#        g.SetMarkerSize(1.25)
#        g.Print() 
#        gs[Vov] = g
#        
#        #g.Draw("same PL")       
#    
#    c = ROOT.TCanvas('c_SR_vs_gainNpe_allVovs_TBLaser','c_SR_vs_gainNpe_allVovs_TBLaser')
#    inFileLaser = ROOT.TFile.Open('c_SR_vs_gainNpe_allVovs.root')
#    mgLaser = inFileLaser.Get('c_SR_vs_gainNpe_allVovs_LYSO818_pulseShapeStudy_vsNpe').FindObject('g_SRglob_vs_gainNpe')
#    
#    c.SetGrid()
#    c.cd()
#    hPad = ROOT.gPad.DrawFrame(0,0.,14000E06, 70.)
#    hPad.SetTitle("; Gain x N_{pe}; Slew Rate [#mu A / ns]")
#     
#    
#    mgLaser.Draw("apl")
#    
#    for Vov in Vovs:
#        gs[Vov].Draw("PL same")
#    
#    
#    c.Update()
#    c.SaveAs(outdir+c.GetName()+'_bar%02d.png'%bar)
        
inFileTB = ROOT.TFile.Open(inputdir+'c_SR_vs_gainNpe_allVovs_barsAve.root')

gs = {}
for i,Vov in enumerate(Vovs):
    g = inFileTB.Get('c_SR_vs_gainNpe_allVovs_barsAve').FindObject('g_pulseShapeL-R_barsAve_Vov%.2f'%Vov)
    g.SetMarkerColor(ROOT.kPink+10)
    g.SetMarkerStyle(markers[Vov])
    g.SetMarkerSize(1.25)
    g.Print() 
    gs[Vov] = g
    
    #g.Draw("same PL")       

c = ROOT.TCanvas('c_SR_vs_gainNpe_allVovs_TBLaser','c_SR_vs_gainNpe_allVovs_TBLaser')

inFileLaser = ROOT.TFile.Open('plots_SRt1_LYSO818_pulseShapeStudy_vsNpe.root')
mgLaser = inFileLaser.Get('g_SRglob_vs_gainNpe_th%02d'%thFixed)
fLaser = inFileLaser.Get('f_SRglob_vs_gainNpe_th%02d'%thFixed)
#mgLaser.FindObject('f_SRglob_vs_gainNpe_th%02d'%thFixed).Delete()
c.SetGrid()
c.cd()
hPad = ROOT.gPad.DrawFrame(0,0.,14000E06, 70.)
hPad.SetTitle("; Gain x N_{pe}; Slew Rate [#mu A / ns]")
 

mgLaser.Draw("apl")
fLaser.SetLineColor(ROOT.kBlack)
fLaser.SetLineStyle(2)

fLaser.Draw("same")
for Vov in Vovs:
    gs[Vov].Draw("PL same")


c.SaveAs(outdir+c.GetName()+'_barsAve_th%02d_v1.png'%thFixed)
        

