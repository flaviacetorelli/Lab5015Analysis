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
#usage: python3 plot_IVcurve.py --outFolder /eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/btlpaper_170424/ --comparison irradiation --IVorDCR IVEff_ch
#outdir = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/IVcurve/test/'
outdir = args.outFolder
comparison = args.comparison
IVorDCR = args.IVorDCR

fnames = {}
gnames = {}
labels = {}

irradiation = ''
SiPM = ''

if (comparison == 'cellsize'):  # USE FILES FROM MIA
    xminleg = 0.76
    yminleg = 0.6
    fnames = { 30 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               25 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               20 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root',
               15 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_size.root'}

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
    SiPM = '#splitline{HPK }{T = -35 #circC}'
    irradiation = '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'
    ypadIV = 3000 
    ypadDCR = 60
    xpad = 2.5

elif comparison == 'temperature': 
    xminleg =  0.75
    yminleg =  0.69

    fnames = {  40 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/logIVEff_TOFHIR2C_Sep23.root',
                35 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/logIVEff_TOFHIR2C_Sep23.root',
                30 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/logIVEff_TOFHIR2C_Sep23.root'
              }


    gnames = {
               40: 'g_%s_HPK_2E14_LYSO815_T-40_ALDO'%IVorDCR,
               35: 'g_%s_HPK_2E14_LYSO815_T-35_ALDO'%IVorDCR,
               30: 'g_%s_HPK_2E14_LYSO815_T-30_ALDO'%IVorDCR
             
             }

    labels = { 40 : 'HPK 25 μm T2 2E+14',
               35 : 'HPK 25 μm T2 2E+14',
               30 : 'HPK 25 um T2 2E+14',
              }
    
    plotAttrs = { 40 : [23, ROOT.kTeal-3, 'T = -40 #circC'],
                  35 : [20, ROOT.kGreen+2, 'T = -35 #circC'],
                  30 : [21, ROOT.kSpring-7, 'T = -30 #circC'],
                }
    ypadIV = 4000 
    ypadDCR = 60
    xpad = 1.4
    SiPM = 'HPK, 25 #mum'
    irradiation = '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'

elif comparison == 'irradiation': # 214 --> 2E14
    xminleg =  0.52
    yminleg =  0.69

    scaleFact = {
              214 : 1,
              114 : pow(1.9, -3. / 10.), # to scale as 1.9 ^ (deltaT / 10)
              113 : pow(1.9, -3. / 10.)
              }
    fnames = { 214 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root',
               114 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root',
               113 : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T1_irradiation.root'
             }

    gnames = { 214 : 'g_%s_conf51.04_ASIC2_ALDO'%IVorDCR, # T = -35 C
               114 : 'g_%s_conf56.00_ASIC2_ALDO'%IVorDCR, # T = -32 C # to be scaled for 3 C
               113 : 'g_%s_conf28.01_ASIC2_ALDO'%IVorDCR, # T = -32 C # to be scaled for 3 C
              }

    labels = { 214 : 'HPK 25 μm T1 2E+14',
               114 : 'HPK 25 μm T1 1E+14',
               113 : 'HPK 25 um T1 1E+13',
              }
    
    plotAttrs = { 214 : [23, ROOT.kGreen+4, '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'],
                 114 : [20, ROOT.kGreen+3,  '1 #times 10^{14} 1 MeV n_{eq}/cm^{2}'],
                 113 : [21, ROOT.kGreen-6,     '1 #times 10^{13} 1 MeV n_{eq}/cm^{2}'],
                }
    ypadIV = 4000 
    ypadDCR = 60
    xpad = 4.5
    SiPM = '#splitline{HPK, 25 #mum}{T = -35 #circC}'

elif comparison == 'vendor': # 214 --> 2E14
    xminleg =  0.66
    yminleg =  0.69
    fnames = { 'fbk' : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_vendors.root',
               'hpk' : '/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/ANALYSIS/TOFHIR2C/IVcurve_240310/logIVEff_TOFHIR2C_T2_vendors.root',
             }

    gnames = { 'fbk' : 'g_%s_conf73.00_ASIC2_ALDO'%IVorDCR, 
               'hpk' : 'g_%s_conf52.00_ASIC2_ALDO'%IVorDCR, 
              }

    labels = { 'fbk' : 'FBK 25 μm T2 2E+14',
               'hpk' : 'HPK 25 μm T2 2E+14',
              }
    
    plotAttrs = { 'fbk' : [23, ROOT.kAzure+7, 'FBK 25 #mum'],
                  'hpk' : [20, ROOT.kGreen+2,  'HPK 25 #mum'],
                }
 
    irradiation = '2 #times 10^{14} 1 MeV n_{eq}/cm^{2}'
    SiPM = 'T = -35 #circC'
    ypadIV = 3000 
    ypadDCR = 60
    xpad = 1.4

             


g = {}
f = {}
for ALDO in ['A', 'B']:
    for key,gname in gnames.items():
        f[key] = ROOT.TFile.Open(fnames[key])
        g[gname+ALDO] = f[key].Get(gname +ALDO)


# averaging
# Using eval function since aldo A and aldo B do not have exact same number of points along x.
for key, gname in gnames.items(): 
    gA = g[gname+'A'] 
    gB = g[gname+'B']
    g[gname+'Ave'] = ROOT.TGraph()
    endrange = gA.GetPointX(gA.GetMaxSize()-1) #choose the common range between aldo A and B to have meaningful evaluate
    if  gB.GetPointX(gB.GetMaxSize()-1) < gA.GetPointX(gA.GetMaxSize()-1) : endrange = gB.GetPointX(gB.GetMaxSize()-1)

    #print (gname,'gA max X = ',  gA.GetPointX(gA.GetMaxSize()-1) , 'gB max X = ', gB.GetPointX(gB.GetMaxSize()-1))
    #print ('endrange  ', endrange)
    npoints = 50 #choose number of points to eval your average
    x0 = 0
    deltax = (endrange - x0) / npoints
    for i in range(0, npoints): #omitting npoints+1 to be completely sure to be in the common range between the two 
        x = i*deltax + x0

        if (comparison == 'irradiation'  and key == 214 and (x > 0.49 and x < 0.555) ):
            print ('removing bad point', gname, x , (gA.Eval(x) +   gB.Eval(x)) / 2)
            continue #omit some no good points
        if (comparison == 'temperature'  and key == 30 and (x > 0.42 and x < 0.465) ): 
            print ('removing bad point', gname, x , (gA.Eval(x) +   gB.Eval(x)) / 2)
            continue #omit some no good points
        g[gname+'Ave'].SetPoint( g[gname+'Ave'].GetN(), x, (gA.Eval(x) +   gB.Eval(x)) / 2 )
        
        #print (gname, x ,  gA.Eval(x) ,   gB.Eval(x) ,  (gA.Eval(x) +   gB.Eval(x)) / 2)
        #print (gname, x , (gA.Eval(x) +   gB.Eval(x)) / 2)


# plot
for ALDO in ['A','B', 'Ave']:

    #leg = ROOT.TLegend(0.20, yminleg, 0.50, 0.89) #aligned on the left
    #if comparison == 'irradiation': leg = ROOT.TLegend(0.50, yminleg, 0.89, 0.89)  #aligned on the right
    leg = ROOT.TLegend(xminleg, yminleg, 0.86, 0.89)  #aligned on the right
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
        #print(g[gname+ALDO].GetN())
        g[gname+ALDO].SetMarkerStyle(plotAttrs[key][0])
        g[gname+ALDO].SetMarkerColor(plotAttrs[key][1])
        g[gname+ALDO].SetMarkerSize(1.15)
        g[gname+ALDO].SetLineColor(plotAttrs[key][1])
        g[gname+ALDO].SetLineWidth(1)

        if comparison == 'irradiation': # scale the 1E14, 1E13 to -35 C
             g[gname+ALDO].Scale(scaleFact[key], "y")
        if (comparison == 'irradiation' and ALDO == 'Ave' and key == 114):
            g[gname+'A'].Draw("pl same") #here ALDO B has too small range, so plotting ALDO A instead of average

        else: g[gname+ALDO].Draw("pl same")
        leg.AddEntry(g[gname+ALDO], '%s'%plotAttrs[key][2],'PL')
    leg.Draw()
    
    tl2 = ROOT.TLatex()
    tl2.SetNDC()
    tl2.SetTextFont(42)
    tl2.SetTextSize(0.045)
    if comparison == 'irradiation' or comparison == 'cellsize': tl2.DrawLatex(0.20,0.82,SiPM)
    else: tl2.DrawLatex(0.20,0.85,SiPM)
    #else: tl2.DrawLatex(0.20,0.20,SiPM)

    tl = ROOT.TLatex()
    tl.SetNDC()
    tl.SetTextFont(42)
    tl.SetTextSize(0.045)
    tl.DrawLatex(0.58,0.20,irradiation)
    
    #cms_logo = draw_logo()
    #cms_logo.Draw()
    
    c.SaveAs(outdir+'%s.png'%c.GetName())
    c.SaveAs(outdir+'%s.pdf'%c.GetName())
    c.SaveAs(outdir+'%s.C'%c.GetName())

#compare A/B/Ave to check average curve is meaningfull
for key,gname in gnames.items():
    
    #leg = ROOT.TLegend(0.20, yminleg, 0.50, 0.89) #aligned on the left
    #if comparison == 'irradiation': leg = ROOT.TLegend(0.50, yminleg, 0.89, 0.89)  #aligned on the right
    leg = ROOT.TLegend(xminleg, yminleg, 0.89, 0.88)  #aligned on the right
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.045) 
    
    c = ROOT.TCanvas('c_%s'%(gname),'c_%s'%(gname), 600, 500)
    if IVorDCR == 'IVEff_ch': 
        ypad = ypadIV
        hPad = ROOT.gPad.DrawFrame(0.,0.,xpad -0.8 ,ypad)
        hPad.SetTitle(";V_{OV} [V]; I [#muA]")
    else: 
        ypad = ypadDCR
        hPad = ROOT.gPad.DrawFrame(0.,0.,xpad-0.8 ,ypad)
        hPad.SetTitle(";V_{OV} [V]; DCR [GHz]")
    hPad.Draw()
    ROOT.gPad.SetTicks(1)

    for i,ALDO in enumerate(['A', 'B', 'Ave']):
        #print(g[gname+ALDO].GetN())
        g[gname+ALDO].SetMarkerStyle(plotAttrs[key][0]+i)
        g[gname+ALDO].SetMarkerColor(plotAttrs[key][1]+3*i)
        g[gname+ALDO].SetMarkerSize(1.15)
        g[gname+ALDO].SetLineColor(plotAttrs[key][1]+3*i)
        g[gname+ALDO].SetLineWidth(1)
        g[gname+ALDO].Draw("pl same")
        leg.AddEntry(g[gname+ALDO], plotAttrs[key][2]+' '+ALDO,'PL')
    leg.Draw()
    
    #if comparison == 'irradiation': tl2.DrawLatex(0.20,0.85,SiPM)
    #else: tl2.DrawLatex(0.20,0.20,SiPM)
    tl2.DrawLatex(0.20,0.85,SiPM)
    
    tl.DrawLatex(0.58,0.20,irradiation)
    
    #cms_logo = draw_logo()
    #cms_logo.Draw()
    
    c.SaveAs(outdir+'%s.png'%c.GetName())
    c.SaveAs(outdir+'%s.pdf'%c.GetName())
    c.SaveAs(outdir+'%s.C'%c.GetName())


