#! /usr/bin/env python
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


outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/uniformity/paper1_Oct24//'

comparison = 'tRes'
#comparison = 'tRes_nonIrrVov3p5'
#comparison = 'energy'
fnames = {}
gnames = {}
labels = {}

irradiation = ''
inputdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/uniformityStudy_4TBpaper_May24/'


if (comparison == 'tRes_nonIrrVov3p5'): 

    SiPM = 'HPK, 25 #mum'
    fnames = { 
               818 : inputdir+'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C/uniformityCheck_HPK_nonIrr_C25_LYSO818_Vov3.50_T5C.root',
               100056 : inputdir+'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C/uniformityCheck_HPK_2E14_C25_LYSO100056_Vov1.50_T-35C.root',
             }

    gnames = { 818 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
               100056 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
              }

    labels = { 818 : 'HPK 25 μm T2 non irr',
               100056 : 'HPK 25 μm T2 2E+14',
              }
    
    plotAttrs = { 
                  818 : [20, ROOT.kRed+1, 'non irradiated, V_{OV} = 3.50 V'],
                  100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}, V_{OV} = 0.96 V'],
                }
 

if (comparison == 'tRes'): 

    SiPM = 'HPK, 25 #mum'
    fnames = { 
               818 : inputdir+'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C/uniformityCheck_HPK_nonIrr_C25_LYSO818_Vov1.00_T5C.root',
               100056 : inputdir+'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C/uniformityCheck_HPK_2E14_C25_LYSO100056_Vov1.50_T-35C.root',
             }

    gnames = { 818 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
               100056 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
              }

    labels = { 818 : 'HPK 25 μm T2 non irr',
               100056 : 'HPK 25 μm T2 2E+14',
              }
    
    plotAttrs = { 
                  818 : [20, ROOT.kGreen+2, 'non irradiated, V_{OV} = 1.00 V'],
                  100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}, V_{OV} = 0.96 V'],
                }
 
if (comparison == 'energy'):  

    SiPM = 'HPK, 25 #mum'
    irradiation = 'non irradiated'
    
    fnames = { 'L'   : inputdir+'uniformityCheck_HPK_nonIrr_C25_LYSO100056_Vov1.50_T-35C.root',
               'R'   : inputdir+'uniformityCheck_HPK_nonIrr_C25_LYSO100056_Vov1.50_T-35C.root',
               'L-R' : inputdir+'uniformityCheck_HPK_nonIrr_C25_LYSO100056_Vov1.50_T-35C.root',
               }

    gnames = { 'L'   : 'energyLNorm_vs_refbar_Vov1.50_th15_bar2', 
               'R'   : 'energyRNorm_vs_refbar_Vov1.50_th15_bar2', 
               'L-R' : 'energyL-RNorm_vs_refbar_Vov1.50_th15_bar2', 
              }

    labels = { 'L' : 'HPK 25 μm T1 2E14',
               'R' : 'HPK 25 μm T1 2E14',
               'L-R' : 'HPK 25 μm T1 2E14',
              }
    
    plotAttrs = { 'L' :   [23, ROOT.kBlue,   'Left'],
                  'R' :   [20, ROOT.kRed,    'Right'],
                  'L-R' : [21, ROOT.kBlack,  'Average'],
                }
 
g = {}
g_mm = {}
f = {}

# plots, vs bar... later vs mm
        

c = ROOT.TCanvas('c_%s_barUniformity'%comparison, 'c_%s_barUniformity'%comparison,  600, 500)
if comparison == 'tRes':
    leg = ROOT.TLegend(0.20, 0.74, 0.50, 0.87)
    hPad = ROOT.gPad.DrawFrame(-0.5,0.,15.5, 120)
    hPad.SetTitle("; reference bar; time resolution [ps]")
elif comparison == 'energy':
    leg = ROOT.TLegend(0.69, 0.72, 0.89, 0.89)
    hPad = ROOT.gPad.DrawFrame(-0.5,0.6,15.5, 1.4)
    hPad.SetTitle("; reference bar; energy [a.u.]")

leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045) 



c.SetGridy()
hPad.Draw()
ROOT.gPad.SetTicks(1)

#since the TB sept data have 3 deg offset
angleEff = 49 
angle = 52 
enScale = math.cos(angleEff*math.pi/180) / math.cos(angle*math.pi/180)

# converting from ref bar to mm
offset = 0 # put to zero, since we don't know where center is actually... counting from the ref module beginning
barConversionFact = 3.12 / math.cos(angleEff*math.pi/180) # [mm]

for key,gname in gnames.items():
    f[key] = ROOT.TFile.Open(fnames[key])
    g[key] = f[key].Get(gname)
    g[key].SetName(gname + "_" + str(key))

    g_mm[key] = ROOT.TGraphErrors() # to conver x-axis in mm
    for i in range(0,g[key].GetN()+1):
        g_mm[key].SetPoint(g_mm[key].GetN(),(g[key].GetPointX(i)-offset)*barConversionFact, g[key].GetPointY(i)/enScale) #accounting for angle offset
        g_mm[key].SetPointError(g_mm[key].GetN()-1, 0, g[key].GetErrorY(i)/enScale )

    ### vs bar plots
    g[key].Scale(1./enScale) # accounting for angle offset
    g[key].SetMarkerStyle(plotAttrs[key][0])
    g[key].SetMarkerColor(plotAttrs[key][1])
    g[key].SetMarkerSize(1.15)

    if (g[key].GetMarkerStyle() == 22): g[key].SetMarkerSize(1.25)
    g[key].SetLineColor(plotAttrs[key][1])
    g[key].SetLineWidth(1)
    g[key].Draw("p same")
    leg.AddEntry(g[key], '%s'%plotAttrs[key][2],'P')


leg.Draw()

tl2 = ROOT.TLatex()
tl2.SetNDC()
tl2.SetTextFont(42)
tl2.SetTextSize(0.045)
if 'tRes' in comparison: tl2.DrawLatex(0.20,0.20,SiPM)
else: tl2.DrawLatex(0.20,0.85,SiPM)

tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextFont(42)
tl.SetTextSize(0.045)
tl.DrawLatex(0.20,0.79,irradiation)

cms_logo = draw_logo()
cms_logo.Draw()

c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())



####### vs mm plots
c1 = ROOT.TCanvas('c_%s_barUniformity_mm'%comparison, 'c_%s_barUniformity_mm'%comparison,  600, 500)
if comparison == 'tRes':
    hPad1 = ROOT.gPad.DrawFrame(10.,0.,60, 120)
    hPad1.SetTitle("; x [mm]; time resolution [ps]")
elif comparison == 'energy':
    hPad1 = ROOT.gPad.DrawFrame(10.,0.6,60., 1.4)
    hPad1.SetTitle("; x [mm]; energy [a.u.]")


c1.SetGridy()
hPad1.Draw()
ROOT.gPad.SetTicks(1)


 
for key,gname in gnames.items():
    g_mm[key].SetMarkerStyle(plotAttrs[key][0])
    g_mm[key].SetMarkerColor(plotAttrs[key][1])
    g_mm[key].SetMarkerSize(1.15)
    if (g_mm[key].GetMarkerStyle() == 22): g_mm[key].SetMarkerSize(1.25)
    g_mm[key].SetLineColor(plotAttrs[key][1])
    g_mm[key].SetLineWidth(1)
    g_mm[key].Draw("p same")
 

####

leg.Draw("same")
tl2 = ROOT.TLatex()
tl2.SetNDC()
tl2.SetTextFont(42)
tl2.SetTextSize(0.045)
if comparison == 'tRes': tl2.DrawLatex(0.20,0.20,SiPM)
elif comparison == 'energy': tl2.DrawLatex(0.20,0.85,SiPM)

tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextFont(42)
tl.SetTextSize(0.045)
tl.DrawLatex(0.20,0.79,irradiation)

cms_logo = draw_logo()
cms_logo.Draw()


c1.SaveAs(outdir+'%s.png'%c1.GetName())
c1.SaveAs(outdir+'%s.pdf'%c1.GetName())
c1.SaveAs(outdir+'%s.C'%c1.GetName())





