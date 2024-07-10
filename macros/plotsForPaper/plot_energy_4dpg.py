#! /usr/bin/env python
# very harcoded macro to produce plots of energy VS impact position (impact position known asking for coincidence with REF module)
# produce linear fit VS x and store in a histograms to extract average slope and spread
# can be used as input for DPG simulation

import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import json

import ROOT
import CMS_lumi, tdrstyle
from utils import *
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
ROOT.gStyle.SetLegendTextSize(0.045)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

#### SELECT MAY or SEPTEMBER

### SEPT TB
#label = 'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' 
label = 'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C' 
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/energyUniformity_4DPG/May24/%s/'%label
inputdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/uniformityStudy_4TBpaper_May24/%s/'%label
offset = 5 
barConversionFact = 0.3 / math.cos(49*math.pi/180) # [cm]
gnames = { 'L'   : 'energyLNorm_vs_refbar_Vov1.50_th15_bar', 
           'R'   : 'energyRNorm_vs_refbar_Vov1.50_th15_bar', 
           'L-R' : 'energyL-RNorm_vs_refbar_Vov1.50_th15_bar', 
          }




### MAY TB
#label = 'HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C' 
#outdir = '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/energyUniformity_4DPG/May24/%s/'%label
#inputdir = '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization//uniformityStudy_TBpaper_Feb24/%s/'%label
#offset = 6 
#barConversionFact = 0.3 / math.cos(52*math.pi/180) # [cm]
#gnames = { 'L'   : 'energyLNorm_vs_refbar_Vov1.00_th05_bar', 
#           'R'   : 'energyRNorm_vs_refbar_Vov1.00_th05_bar', 
#           'L-R' : 'energyL-RNorm_vs_refbar_Vov1.00_th05_bar', 
#          }


fnames = {}
gnames = {}

irradiation = ''
inFile = ROOT.TFile.Open( inputdir+'uniformityCheck_%s.root'%label , "READ")
SiPM = 'HPK, 25 #mum'
irradiation = 'non irradiated'
hslope = {}
for key in ['L', 'R', 'L-R']:
    hslope [key] = ROOT.TH1F("hslope","hslope",100, -0.25, 0.25)

plotAttrs = { 'L' :   [23, ROOT.kBlue,   'Left'],
              'R' :   [20, ROOT.kRed,    'Right'],
              'L-R' : [21, ROOT.kBlack,  'Average'],
}


print (barConversionFact)

g = {}
g_mm = {}
f = {}
goodbars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
# plots, vs bar... later vs mm
fitLin = {}

print (inputdir)

for bar in goodbars:
    c = ROOT.TCanvas('c_energyUniformity_bar%02d'%bar,'c_energyUniformity_bar%02d'%bar ,  600, 500)
    leg = ROOT.TLegend(0.69, 0.72, 0.89, 0.89)
    
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045) 
    
    
    for key,gname in gnames.items():
        g[key] = inFile.Get(gname + str(bar))
        #print (gname + str(bar))
        g[key].SetName(gname + str(bar))
        g_mm[key] = ROOT.TGraphErrors() # to make axis in cm
        for i in range(0,g[key].GetN()+1):
            g_mm[key].SetPoint(g_mm[key].GetN(),(g[key].GetPointX(i)-offset)*barConversionFact, g[key].GetPointY(i))
            g_mm[key].SetPointError(g_mm[key].GetN()-1, 0, g[key].GetErrorY(i) )
    
   
    
    
    
    hPad1 = ROOT.gPad.DrawFrame(-3.2,0.4,3.2, 1.6)
    hPad1.SetTitle("; distance from bar center [cm]; energy [a.u.]")
    
    c.SetGridy()
    hPad1.Draw()
    ROOT.gPad.SetTicks(1)
    tl = ROOT.TLatex()
    tl.SetNDC()
    tl.SetTextFont(42)
    tl.SetTextSize(0.045)
    
    
    i = 0
    for key,gname in gnames.items():
        g_mm[key].SetMarkerStyle(plotAttrs[key][0])
        g_mm[key].SetMarkerColor(plotAttrs[key][1])
        g_mm[key].SetMarkerSize(1.15)
        if (g_mm[key].GetMarkerStyle() == 22): g_mm[key].SetMarkerSize(1.25)
        g_mm[key].SetLineColor(plotAttrs[key][1])
        g_mm[key].SetLineWidth(1)
        g_mm[key].Draw("p same")
        #if key == "L-R": continue
        fitLin["%02d%s"%(bar,key)] = ROOT.TF1("pol1_%02d%s"%(bar,key), "pol1", -3, 3) 
        fitLin["%02d%s"%(bar,key)].SetLineColor((plotAttrs[key][1]))
        g_mm[key].Fit(fitLin["%02d%s"%(bar,key)], "QRS")
        tl.DrawLatex(0.20,0.85-i,"%s fit %.4f x [a.u. / cm] + %.4f"%(key, fitLin["%02d%s"%(bar,key)].GetParameter(1), fitLin["%02d%s"%(bar,key)].GetParameter(0)))
        i = i +0.05

        #if key == 'L-R': continue
        hslope[key].Fill(fitLin["%02d%s"%(bar,key)].GetParameter(1))

    cms_logo = draw_logo()
    cms_logo.Draw()
    c.SaveAs(outdir+'%s.png'%c.GetName())
    c.SaveAs(outdir+'%s.pdf'%c.GetName())
    c.SaveAs(outdir+'%s.C'%c.GetName())

for key,h in hslope.items():
    c.Clear()
    c.cd()
    h.Draw("")
    gaus = ROOT.TF1("gaus","gaus",h.GetMean()-3*h.GetRMS(), h.GetMean()+3*h.GetRMS())
    print (h.GetMean(), h.GetRMS())
    h.Fit(gaus, "R")
    h.GetXaxis().SetRangeUser(h.GetMean()-5*h.GetRMS(),  h.GetMean()+5*h.GetRMS())
    print (gaus.GetParameter(2)/gaus.GetParameter(1))
    c.SaveAs(outdir+"/slope%s.png"%key)
