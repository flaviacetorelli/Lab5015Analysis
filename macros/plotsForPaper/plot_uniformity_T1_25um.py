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
ROOT.gStyle.SetLabelSize(0.052,'X') #0.055 before
ROOT.gStyle.SetLabelSize(0.052,'Y')
ROOT.gStyle.SetTitleSize(0.066,'X') #0.07 before
ROOT.gStyle.SetTitleSize(0.067,'Y')
ROOT.gStyle.SetTitleOffset(0.95,'X')
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.05) 
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadTopMargin(0.13)
ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning



outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/uniformity/paper1_Feb25/'

#comparison = 'tRes'
comparison = 'tRes_nonIrrVov3p5'
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
    
    #plotAttrs = { 
    #              818 : [20, ROOT.kGreen+2, 'non irradiated, V_{OV} = 3.50 V'],
    #              100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}, V_{OV} = 0.96 V'],
    #            }
    plotAttrs = { 
                  818 : [20, ROOT.kGreen+2, 'non-irradiated'],
                  100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} n_{eq}/cm^{2}'],
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
    
    #plotAttrs = { 
    #              818 : [20, ROOT.kBlue, 'non irradiated, V_{OV} = 1.00 V'],
    #              100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}, V_{OV} = 0.96 V'],
    #            }
    plotAttrs = { 
                  818 : [20, ROOT.kBlue, 'non-irradiated'],
                  100056 : [22, ROOT.kOrange+1,  '2 #times 10^{14} n_{eq}/cm^{2}'],
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
if 'tRes' in comparison:
    leg = ROOT.TLegend(0.20, 0.66, 0.50, 0.84)
    hPad = ROOT.gPad.DrawFrame(-0.5,0.,15.5, 120)
    hPad.SetTitle("; reference bar; time resolution [ps]")
else:
    leg = ROOT.TLegend(0.69, 0.72, 0.89, 0.89)
    hPad = ROOT.gPad.DrawFrame(-0.5,0.6,15.5, 1.4)
    hPad.SetTitle("; reference bar; energy [a.u.]")

leg.SetBorderSize(0)
leg.SetFillStyle(0)



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
    for i in range(g[key].GetN()):
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
    leg.AddEntry(g[key], '%s'%plotAttrs[key][2],'PL')


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

#cms_logo = draw_logo()
#cms_logo.Draw()

c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())



####### vs mm plots
c1 = ROOT.TCanvas('c_%s_barUniformity_mm'%comparison, 'c_%s_barUniformity_mm'%comparison,  600, 500)
if 'tRes' in comparison:
    hPad1 = ROOT.gPad.DrawFrame(-0.5*barConversionFact,0.,15.5*barConversionFact, 120)
    hPad1.SetTitle("; x [mm]; time resolution [ps]")
else:
    hPad1 = ROOT.gPad.DrawFrame(10.,0.6,60., 1.4)
    hPad1.SetTitle("; x [mm]; energy [a.u.]")


c1.SetGridy()
hPad1.Draw()
#ROOT.gPad.SetTicks(1)
ROOT.gPad.SetTicky(1)
ROOT.gPad.SetTickx(0)


 
for key,gname in gnames.items():
    g_mm[key].SetMarkerStyle(plotAttrs[key][0])
    g_mm[key].SetMarkerColor(plotAttrs[key][1])
    g_mm[key].SetMarkerSize(1.15)
    if (g_mm[key].GetMarkerStyle() == 22): g_mm[key].SetMarkerSize(1.25)
    g_mm[key].SetLineColor(plotAttrs[key][1])
    g_mm[key].SetLineWidth(1)
    g_mm[key].Draw("p same")
 

#### new ax with bar info
f1 = ROOT.TF1("f1","x",-0.5 ,15.5);
xaxis2 = ROOT.TGaxis(-0.5*barConversionFact, 120 , 15.5*barConversionFact, 120,"f1",512,"-")
xaxis2.SetTitle("reference module bar")
xaxis2.Draw("same")
xaxis2.SetLabelSize(0.052)
xaxis2.SetTitleSize(0.06)
xaxis2.SetTitleOffset(1.05)
xaxis2.SetTitleFont(42)
xaxis2.SetLabelFont(42)


leg.Draw("same")
#tl2 = ROOT.TLatex()
#tl2.SetNDC()
#tl2.SetTextFont(42)
#tl2.SetTextSize(0.045)
#if 'tRes' in comparison: tl2.DrawLatex(0.20,0.18,SiPM)
#else: tl2.DrawLatex(0.20,0.85,SiPM)

tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextFont(42)
tl.SetTextSize(0.045)
tl.DrawLatex(0.20,0.79,irradiation)

#cms_logo = draw_logo()
#cms_logo.Draw()


c1.SaveAs(outdir+'%s.png'%c1.GetName())
c1.SaveAs(outdir+'%s.pdf'%c1.GetName())
c1.SaveAs(outdir+'%s.C'%c1.GetName())

#### with a pol0 fit to check consistency of tRefs numbers 
pol0_noirr = ROOT.TF1("noirr", "pol0", -0.5*barConversionFact, 15.5*barConversionFact)
pol0_irr = ROOT.TF1("irr", "pol0", -0.5*barConversionFact, 15.5*barConversionFact)

g_mm[818].Fit(pol0_noirr, "R")
pol0_noirr.SetLineColor(g[818].GetLineColor())
t1 = ROOT.TLatex()
t1.SetNDC()
t1.SetTextFont(42)
t1.SetTextSize(0.045)
t1.DrawLatex(0.20,0.35,"p_0 non irr: %.2f #pm %.2f"%(pol0_noirr.GetParameter(0), pol0_noirr.GetParError(0)))
pol0_noirr.Draw("same")

g_mm[100056].Fit(pol0_irr, "R")
pol0_irr.SetLineColor(g[100056].GetLineColor())
t2 = ROOT.TLatex()
t2.SetNDC()
t2.SetTextFont(42)
t2.SetTextSize(0.045)
t2.DrawLatex(0.20,0.55,"p_0 2E14: %.2f #pm %.2f"%(pol0_irr.GetParameter(0), pol0_irr.GetParError(0)))
pol0_irr.Draw("same")


c1.SaveAs(outdir+'%s_fit.png'%c1.GetName())
c1.SaveAs(outdir+'%s_fit.pdf'%c1.GetName())
c1.SaveAs(outdir+'%s_fit.C'%c1.GetName())






