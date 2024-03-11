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
ROOT.gStyle.SetOptFit(1)
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

parser = argparse.ArgumentParser(description='Plots IV or DCR curves')
parser.add_argument("--comparison",   required=True, type=str, help="Type of comparison: cellsize, irradiation, vendor")
parser.add_argument("--IVorDCR",       required=True, type=str, help="Choose IVEff_ch or DCR")
parser.add_argument("--outFolder", required=True, type=str, help="out folder: choose an existing ones")
args = parser.parse_args()
#usage: python3 plot_IVcurve.py --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/test/ --comparison irradiation --IVorDCR IVEff_ch
#outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/test/'
outdir = args.outFolder
comparison = args.comparison
IVorDCR = args.IVorDCR

fnames = {}
gnames = {}
labels = {}

irradiation = ''
SiPM = 'HPK'

if (comparison == 'cellsize'):  # USE FILES FROM MIA
    yminleg = 0.6
    fnames = { 30 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               25 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               20 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               15 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root'}

    gnames = { 30 : 'g_%s_conf54.01_ASIC2_ALDO'%IVorDCR,
               25 : 'g_%s_conf52.00_ASIC2_ALDO'%IVorDCR,
               20 : 'g_%s_conf55.01_ASIC2_ALDO'%IVorDCR,
               15 : 'g_%s_conf6.00_ASIC0_ALDO'%IVorDCR,
              }

    labels = { 30 : 'HPK 30 μm T2 2E+14',
               25 : 'HPK 25 μm T2 2E+14',
               20 : 'HPK 20 μm T2 2E+14',
               15 : 'HPK 15 μm T2 2E+14',
              }

    plotAttrs = { 30 : [23, ROOT.kOrange+1, '30 #mum'],
                  25 : [20, ROOT.kGreen+2,  '25 #mum'],
                  20 : [21, ROOT.kBlue,     '20 #mum'],
                  15 : [22, ROOT.kRed,      '15 #mum']
                }
    irradiation = '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'
    ypadIV = 3000 
    ypadDCR = 80
    xpad = 2.5

elif comparison == 'irradiation': # 214 --> 2E14

    yminleg =  0.69
    fnames = { 214 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root',
               114 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root',
               113 : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root'
             }

    gnames = { 214 : 'g_%s_conf51.04_ASIC2_ALDO'%IVorDCR, # T = -35 C # to be scaled for 3 C
               114 : 'g_%s_conf56.00_ASIC2_ALDO'%IVorDCR, # T = -32 C
               113 : 'g_%s_conf28.01_ASIC2_ALDO'%IVorDCR, # T = -32 C
              }

    labels = { 214 : 'HPK 25 μm T1 2E+14',
               114 : 'HPK 25 μm T1 1E+14',
               113 : 'HPK 25 um T1 1E+13',
              }
    
    plotAttrs = { 214 : [23, ROOT.kOrange+1, '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'],
                 114 : [20, ROOT.kGreen+2,  '1 #times 10^{14} 1 MeV n_{eq}/cm^{2}'],
                 113 : [21, ROOT.kBlue,     '1 #times 10^{13} 1 MeV n_{eq}/cm^{2}'],
                }
    ypadIV = 4000 
    ypadDCR = 50
    xpad = 4.5

elif comparison == 'vendor': # 214 --> 2E14

    yminleg =  0.69
    fnames = { 'fbk' : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_vendors.root',
               'hpk' : '/afs/cern.ch/work/f/fcetorel/private/work2/TBMay2023/Lab5015Analysis/plots/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_vendors.root',
             }

    gnames = { 'fbk' : 'g_%s_conf73.00_ASIC2_ALDO'%IVorDCR, 
               'hpk' : 'g_%s_conf52.00_ASIC2_ALDO'%IVorDCR, 
              }

    labels = { 'fbk' : 'FBK 25 μm T2 2E+14',
               'hpk' : 'HPK 25 μm T2 1E+14',
              }
    
    plotAttrs = { 'fbk' : [23, ROOT.kOrange+1, 'FBK'],
                  'hpk' : [20, ROOT.kGreen+2,  'HPK'],
                }
 
    irradiation = '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'
    SiPM = ''
    ypadIV = 3000 
    ypadDCR = 50
    xpad = 2.5

             


g = {}
f = {}

# plot
        
for ALDO in ['A','B']:

    leg = ROOT.TLegend(0.20, yminleg, 0.50, 0.89) #aligned on the left
    if comparison == 'irradiation': leg = ROOT.TLegend(0.50, yminleg, 0.89, 0.89)  #aligned on the right
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045) 
    
    c = ROOT.TCanvas('c_%s_%s_ALDO%s'%(IVorDCR,comparison, ALDO),'c_%s_%s_ALDO%s'%(IVorDCR, comparison, ALDO), 600, 500)
    if IVorDCR == 'IVEff_ch': 
        ypad = ypadIV
        hPad = ROOT.gPad.DrawFrame(0.,0.,xpad,ypad)
        hPad.SetTitle(";V_{OV} [V]; I [#muA]")
    else: 
        ypad = ypadDCR
        hPad = ROOT.gPad.DrawFrame(0.,0.,xpad,ypad)
        hPad.SetTitle(";V_{OV} [V]; DCR [GHz]")
    hPad.Draw()
    ROOT.gPad.SetTicks(1)

    for key,gname in gnames.items():
        f[key] = ROOT.TFile.Open(fnames[key])
        g[key] = f[key].Get(gname +ALDO)
        g[key].SetMarkerStyle(plotAttrs[key][0])
        g[key].SetMarkerColor(plotAttrs[key][1])
        g[key].SetMarkerSize(1.15)
        g[key].SetLineColor(plotAttrs[key][1])
        g[key].SetLineWidth(1)
        g[key].Draw("pl same")
        leg.AddEntry(g[key], '%s'%plotAttrs[key][2],'PL')
    leg.Draw()
    
    tl2 = ROOT.TLatex()
    tl2.SetNDC()
    tl2.SetTextFont(42)
    tl2.SetTextSize(0.045)
    tl2.DrawLatex(0.20,0.20,SiPM)
    
    tl = ROOT.TLatex()
    tl.SetNDC()
    tl.SetTextFont(42)
    tl.SetTextSize(0.045)
    tl.DrawLatex(0.58,0.20,irradiation)
    
    cms_logo = draw_logo()
    cms_logo.Draw()
    
    c.SaveAs(outdir+'%s.png'%c.GetName())
    c.SaveAs(outdir+'%s.pdf'%c.GetName())
    c.SaveAs(outdir+'%s.C'%c.GetName())



