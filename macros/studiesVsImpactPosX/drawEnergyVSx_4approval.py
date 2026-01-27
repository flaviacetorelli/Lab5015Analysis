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
import tdrstyle
from plots_header import *

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X') #0.052
ROOT.gStyle.SetLabelSize(0.052,'Y')
ROOT.gStyle.SetTitleSize(0.07,'X') #0.067
ROOT.gStyle.SetTitleSize(0.07,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X') # 0.95
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.05) 
#ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadTopMargin(0.07) #0.13
#ROOT.gStyle.SetPadRightMargin(0.05)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning




outdir = '/eos/user/f/fcetorel/www/MTD/DPG/digiUpdate_approvalPlots_Dec25/'
inputfile = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/DPGstudies_4DigitizationModel/energy/Sept25_uncUpdate/HPK_nonIrr_C25_LYSO818_Vov1.00_T5C/c_energy_vs_x_bar11_ref15.root'

labels = {
        "L" : "left SiPM", 
        "R" : "right SiPM", 
        "L-R" : "average", 

        }

g = {}
g_mm = {}
f = {}

# plots, vs bar... later vs mm
        

c = ROOT.TCanvas('c_energy', 'c_energy',  600, 500)
hPad = ROOT.gPad.DrawFrame(0,0.7,6., 1.3)
hPad.SetTitle("; Hodoscope [cm]; Energy [a.u.]")
hPad.Draw("")


c.SetGridy()
hPad.Draw()
ROOT.gPad.SetTicks(1)
c1 = ROOT.TCanvas()
f = ROOT.TFile.Open(inputfile)
c1 = f.Get('c_energy_vs_x_bar11_ref15')
goff = {}

leg = ROOT.TLegend(0.25, 0.75, 0.5, 0.9)
#leg = ROOT.TLegend(0.69, 0.72, 0.89, 0.89)
leg.SetBorderSize(0)
leg.SetFillStyle(0)


for l in ["L", "R", "L-R"]:
    offset = 0
    g = c1.GetListOfPrimitives().FindObject('g_bar11%s'%l)
    #g.Print()
    if l == 'L': offset =  g.GetPointY(4) - 0.99884
    elif l == 'R': offset =  g.GetPointY(4) - 0.99884
    print (offset)
    goff[l] = ROOT.TGraphErrors()
    goff[l].SetMarkerStyle(g.GetMarkerStyle())
    goff[l].SetMarkerSize(g.GetMarkerSize())
    goff[l].SetMarkerColor(g.GetMarkerColor())


    for i in range(0, g.GetN()):
       goff[l].SetPoint(goff[l].GetN(), g.GetPointX(i), g.GetPointY(i)- offset)
       goff[l].SetPointError(goff[l].GetN()-1, g.GetErrorX(i), g.GetErrorY(i))
    print (l)        
    goff[l].Print()

    c.cd()
    goff[l].Draw("P same")
    pol1 = ROOT.TF1("pol%s"%l, "pol1", goff[l].GetPointX(0)-0.1, goff[l].GetPointX(8)+0.1)
    pol1.SetLineColor(goff[l].GetMarkerColor())
    pol1.SetLineStyle(2)
    goff[l].Fit(pol1, "R")
    pol1.Draw("same")
    leg.AddEntry(goff[l], labels[l], "pl")
leg.Draw("same")

cms_logo = draw_logo()
cms_logo.Draw()

c.SaveAs("%s/%s.png"%(outdir, c.GetName()))


