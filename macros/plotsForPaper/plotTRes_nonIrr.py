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



outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/fbk/minorRevJINST_Nov24/'

#sipmProd = 'HPK'
sipmProd = 'FBK'

enScale = math.cos(49.*math.pi/180.)/math.cos(52.*math.pi/180.) # for 3 deg angle offset in Sep2023 TB
print(enScale)

fnames = {}
gnames = {}
labels = {}

cells = [15, 20, 25, 30]

if (sipmProd == 'HPK'):  # USE FILES FROM SIMONA
    fnames = { 30 : '',
               25 : '',
               20 : '',
               15 : ''}

    gnames = { 30 : 'g_data_vs_Vov_average_HPK_',
               25 : 'g_data_vs_Vov_average_HPK_',
               20 : 'g_data_vs_Vov_average_HPK_',
               15 : 'g_data_vs_Vov_average_HPK_',
              }

    labels = { 30 : '',
               25 : '',
               20 : '',
               15 : '',
              }

else: # FILES FROM JIN WANG (20,25,30) and MARTINA (15 un June22 TB)
    fnames = { 30 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/fbk_240319/tres_vs_Vov_all_result.root',
               25 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/fbk_240319/tres_vs_Vov_all_result.root',
               20 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/fbk_240319/tres_vs_Vov_all_result.root',
               15 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/RootFiles/plots_timeResolution_HPK_FBK_nonIrr_TBJun22_TOFHIR2X.root'}

    gnames = { 30 : 'module28_LYSO_200046_FBK_C30_Rq2 _T2_nonirr',
               25 : 'module27_LYSO200075_FBK_C25_Rq2_T2_nonirr',
               20 : 'module29_LYSO200074_FBK_C20_Rq2_T2_nonirr',
               15 : 'g_data_vs_Vov_average_FBK_nonIrr_LYSO800_T10C',
              }


    labels = { 30 : '',
               25 : '',
               20 : '',
               15 : 'FBK_nonIrr_LYSO800_T10C',
              }
    
    


              
plotAttrs = { 30 : [23, ROOT.kOrange+1, '30 #mum'],
              25 : [20, ROOT.kGreen+2,  '25 #mum'],
              20 : [21, ROOT.kBlue,     '20 #mum'],
              15 : [22, ROOT.kRed,      '15 #mum']}



g = {}
g_scaled = {}
g_Noise = {}
g_Stoch = {}
g_SR = {}
f = {}

for cell in cells:
    f[cell] = ROOT.TFile.Open(fnames[cell])
    g[cell] = f[cell].Get(gnames[cell])
    g_scaled[cell] = ROOT.TGraphErrors()
    g_scaled[cell].SetName(gnames[cell].replace('g_data','g_data_scaled'))

# scale (2C) to take into account angle offset in 2023 Sep TB 
for cell in [20, 25, 30]:
    if (cell not in cells): continue
    if (sipmProd == 'HPK'):
        gNoise[cell] = f[cell].Get('g_Noise_vs_Vov_average_%s'%labels[cell])
        gStoch[cell] = f[cell].Get('g_Stoch_vs_Vov_average_%s'%labels[cell])
        gSR[cell]   = f[cell].Get('g_SR_vs_Vov_average_%s'%labels[cell])
        for i in range(0, g[cell].GetN()):
            vov = g[cell].GetX()[i]
            sr = gSR[cell].Eval(vov)
            s_noise =  sigma_noise(sr*enScale, '2C')
            s_stoch = gStoch[cell].Eval(vov)/math.sqrt(enScale)
            s_dcr = 0.
            s_tot = math.sqrt(s_noise*s_noise + s_stoch*s_stoch + s_dcr*s_dcr)
            g_scaled[cell].SetPoint(i, vov, s_tot) # correct for angle offset 
            g_scaled[cell].SetPointError(i, 0, g[cell].GetErrorY(i)/enScale) # correct for angle offset
    else:
        # for FBK just scale the total resolution with enScale as we don't have the different contributions separately
        for i in range(0, g[cell].GetN()):
            vov = g[cell].GetX()[i]
            s_tot = g[cell].GetY()[i]/enScale
            g_scaled[cell].SetPoint(i, vov, s_tot) # correct for angle offset
            g_scaled[cell].SetPointError(i, 0, g[cell].GetErrorY(i)/enScale) # correct for angle offset         

        
# plot        
leg = ROOT.TLegend(0.70, 0.60, 0.89, 0.89)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.045) 

c = ROOT.TCanvas('c_timeResolution_%s_nonIrr_vs_Vov'%sipmProd,'c_timeResolution_%s_nonIrr_vs_Vov'%sipmProd, 600, 500)
xmin = 0.
xmax = 4.0
if (sipmProd == 'FBK'): xmax = 7.5
hPad = ROOT.gPad.DrawFrame(xmin, 0., xmax, 160.)
hPad.SetTitle(";V_{OV} [V];time resolution [ps]")
hPad.Draw()
ROOT.gPad.SetTicks(1)
for cell in cells:
    g_scaled[cell].SetMarkerSize(1)
    if (plotAttrs[cell][0] == 22 or plotAttrs[cell][0] == 23): g_scaled[cell].SetMarkerSize(1.15)
    g_scaled[cell].SetMarkerStyle(plotAttrs[cell][0])
    g_scaled[cell].SetMarkerColor(plotAttrs[cell][1])
    g_scaled[cell].SetLineColor(plotAttrs[cell][1])
    leg.AddEntry(g_scaled[cell], '%s'%plotAttrs[cell][2],'PL')
    g[cell].SetMarkerStyle(plotAttrs[cell][0])
    g[cell].SetMarkerColor(plotAttrs[cell][1])
    g[cell].SetMarkerSize(1.15)
    g[cell].SetLineColor(plotAttrs[cell][1])
    g[cell].SetLineWidth(1)
    if (cell != 15):
        g_scaled[cell].Draw('plsame')
    else:
        g[cell].Draw('plsame')
leg.Draw()

tl2 = ROOT.TLatex()
tl2.SetNDC()
tl2.SetTextFont(42)
tl2.SetTextSize(0.045)
tl2.DrawLatex(0.20,0.86,'%s'%sipmProd)

tl = ROOT.TLatex()
tl.SetNDC()
tl.SetTextFont(42)
tl.SetTextSize(0.045)
tl.DrawLatex(0.20,0.80,'non-irradiated')

# Draw dashed lines at y=30 and y=60
line1 = ROOT.TLine(xmin, 30, xmax, 30) 
line1.SetLineStyle(2)              
line1.SetLineColor(ROOT.kGray + 1) 
line1.Draw("same")

line2 = ROOT.TLine(xmin, 60, xmax, 60)
line2.SetLineStyle(2)              
line2.SetLineColor(ROOT.kGray + 1) 
line2.Draw("same")

#cms_logo = draw_logo()
#cms_logo.Draw()

c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())



