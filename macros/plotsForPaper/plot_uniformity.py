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



outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/uniformity/test/'

comparison = 'tRes'
#comparison = 'energy'
fnames = {}
gnames = {}
labels = {}

irradiation = ''

if (comparison == 'tRes'): 

    SiPM = 'HPK, 25 #mum'
    fnames = { 813 : '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Feb24/HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C/uniformityCheck_HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C.root',
               815 : '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Feb24/HPK_2E14_C25_LYSO815_Vov1.50_T-30C/uniformityCheck_HPK_2E14_C25_LYSO815_Vov1.50_T-30C.root',
             }

    gnames = { 813 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
               815 : 'g_tRes_average_deltaT_totRatioCorr_bestTh_vs_bar', 
              }

    labels = { 813 : 'HPK 25 μm T2 non irr',
               815 : 'HPK 25 μm T2 2E+14',
              }
    
    plotAttrs = { 813 : [23, ROOT.kOrange+1, 'non irradiated, V_{OV} = 1.00 V'],
                  815 : [20, ROOT.kGreen+2,  '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}, V_{OV} = 0.93 V'],
                }
 
if (comparison == 'energy'):  

    SiPM = 'HPK, 25 #mum'
    irradiation = 'non irradiated'
    
    fnames = { 'L'   : '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Feb24/HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C/uniformityCheck_HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C.root',
               'R'   : '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Feb24/HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C/uniformityCheck_HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C.root',
               'L-R' : '/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Feb24/HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C/uniformityCheck_HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C.root',
               }

    gnames = { 'L'   : 'energyLNorm_vs_refbar_Vov1.00_th05_bar2', 
               'R'   : 'energyRNorm_vs_refbar_Vov1.00_th05_bar2', 
               'L-R' : 'energyL-RNorm_vs_refbar_Vov1.00_th05_bar2', 
              }

    labels = { 'L' : 'HPK 25 μm T2 non irr',
               'R' : 'HPK 25 μm T2 2E+14',
               'L-R' : 'HPK 25 μm T2 2E+14',
              }
    
    plotAttrs = { 'L' :   [23, ROOT.kBlue,   'Left'],
                  'R' :   [20, ROOT.kRed,    'Right'],
                  'L-R' : [21, ROOT.kBlack,  'Average'],
                }
 
g = {}
f = {}

# plot
        

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

for key,gname in gnames.items():
    f[key] = ROOT.TFile.Open(fnames[key])
    g[key] = f[key].Get(gname)
    g[key].SetMarkerStyle(plotAttrs[key][0])
    g[key].SetMarkerColor(plotAttrs[key][1])
    g[key].SetMarkerSize(1.15)
    g[key].SetLineColor(plotAttrs[key][1])
    g[key].SetLineWidth(1)
    g[key].Draw("p same")
    leg.AddEntry(g[key], '%s'%plotAttrs[key][2],'P')
leg.Draw()

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

c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())



