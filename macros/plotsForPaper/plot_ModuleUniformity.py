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


colors = {
        'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : ROOT.kBlue,  
        'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : ROOT.kGreen + 2, 
        'HPK_2E14_LYSO100056_T-35C' : ROOT.kOrange + 1, 

        }

#since the TB sept data have 3 deg offset
angleEff = 49
angle = 52
enScale = math.cos(angleEff*math.pi/180) / math.cos(angle*math.pi/180)
barConversionFact = 3.12

from VovsEff import *
# Import file with VovEff and DCR
with open('/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/VovsEff_TOFHIR2C.json', 'r') as f:
   data = json.load(f)


inputdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/'
outdir   = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/moduleUniformity/paper1_Feb25/'

#modules = ['HPK_nonIrr_C25_LYSO818_Vov1.00_T5C', 'HPK_2E14_LYSO100056_T-35C']
#modules = {'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C','HPK_nonIrr_C25_LYSO818_Vov3.50_T5C', 'HPK_2E14_LYSO100056_T-35C'}
modules = ['HPK_nonIrr_C25_LYSO818_Vov3.50_T5C', 'HPK_2E14_LYSO100056_T-35C']

fnames = { 
           'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : inputdir+'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C_refbar7/summaryPlots_HPK_nonIrr_C25_LYSO818_Vov1.00_T5C_refbar7.root',
           'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : inputdir+'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C_refbar7/summaryPlots_HPK_nonIrr_C25_LYSO818_Vov3.50_T5C_refbar7.root',
           'HPK_2E14_LYSO100056_T-35C' : inputdir+'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C_refbar7/summaryPlots_HPK_2E14_C25_LYSO100056_Vov1.50_T-35C_refbar7.root',
       }

labels = {
          'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : 'non irradiated',
          'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : 'non-irradiated',
          'HPK_2E14_LYSO100056_T-35C' : '2 x 10^{14} n_{eq}/cm^{2}',
     }


Vovs = { 
         'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : [1.00],
         'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : [3.50],
         'HPK_2E14_LYSO100056_T-35C' : [1.50],
}

bestVovs = { 
             'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : 1.00,
             'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : 3.50,
             'HPK_2E14_LYSO100056_T-35C' : 1.50,
          }


if '3.5' in modules[0] or '3.5' in modules[1]:
    c = ROOT.TCanvas('c_timeResolution_vs_bar_nonIrrVov3p5','c_timeResolution_vs_bar_nonIrrVov3p5', 600, 500)
else:
    c = ROOT.TCanvas('c_timeResolution_vs_bar','c_timeResolution_vs_bar', 600, 500)
#hPad = ROOT.TH2F('hPad','', 100, -0.5, 15.5, 100, 0, 120)
padXmin = -0.5*barConversionFact 
padXmax = 15.5*barConversionFact
hPad = ROOT.TH2F('hPad','', 100, padXmin, padXmax, 100, 0, 120)
hPad.SetTitle("; y [mm]; time resolution [ps]")
hPad.Draw()
c.SetGridy()
c.SetTicky(1)
c.SetTickx(0)


leg = ROOT.TLegend(0.20, 0.66, 0.50, 0.84)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

# second axis with bar number
f1 = ROOT.TF1("f1","x", padXmin/barConversionFact , padXmax/barConversionFact);
xaxis2 = ROOT.TGaxis(padXmin, 120 , padXmax, 120,"f1",512,"-")
xaxis2.SetTitle("DUT bar")
xaxis2.Draw("same")
xaxis2.SetLabelSize(0.052)
xaxis2.SetTitleSize(0.06)
xaxis2.SetTitleOffset(1.05)
xaxis2.SetTitleFont(42)
xaxis2.SetLabelFont(42)


h = {}
g_mm = {}

for mod in modules:

   #h [mod] = ROOT.TH1F()
   f = ROOT.TFile.Open(fnames[mod])
   g_mm[mod]= ROOT.TGraphErrors()   

   for iv,vov in enumerate(Vovs[mod]):
         
      g = f.Get('g_deltaT_totRatioCorr_bestTh_vs_bar_Vov%.2f_enBin01'%vov) 
      print('g_deltaT_totRatioCorr_bestTh_vs_bar_Vov%.2f_enBin01'%vov)
      print(g.GetN())

      for i in range(g.GetN()): #conversion in mm
          g_mm[mod].SetPoint(g_mm[mod].GetN(), g.GetPointX(i)*barConversionFact, g.GetPointY(i)/enScale) #accounting for angle offset
          g_mm[mod].SetPointError(g_mm[mod].GetN()-1, 0, g.GetErrorY(i)/enScale )


      #histogram for spread
      if (vov == bestVovs[mod]): 
          hdummy = ROOT.TH1F('h_%s'%mod,'h_%s'%mod, 60, -0.8,0.8 )
          #print (mod)
          #pol0= ROOT.TF1(mod, "pol0", -0.5*barConversionFact, 15.5*barConversionFact)
          #g_mm[mod].Fit(pol0, "R")
          for i in range(g_mm[mod].GetN()):
             #x = ( g_mm[mod].GetPointY(i) - pol0.GetParameter(0) )/pol0.GetParameter(0)
             x = ( g_mm[mod].GetPointY(i) - g_mm[mod].GetMean(2) ) / g_mm[mod].GetMean(2)
             hdummy.Fill(x)
             #print (i, " ", g_mm[mod].GetPointY(i),  " ", g_mm[mod].GetMean(2) , " ",  x)
             #print (i, " ", g_mm[mod].GetPointY(i),  " ", pol0.GetParameter(0) , " ",  x)
          hdummy.SetLineColor(g_mm[mod].GetLineColor())
          hdummy.SetFillColorAlpha(g_mm[mod].GetLineColor(),0.2)

          #Finally the spread resolution plot
          c2 = ROOT.TCanvas('c_spread_%s'%mod,'c_spread_%s'%mod, 600, 500)
          c2.cd()
          c2.Update()
          hdummy.Draw("")

          hdummy.GetXaxis().SetTitle("(#sigma_{t} - <#sigma_{t}>) / <#sigma_{t}>  ")
          tx = ROOT.TLatex()
          tx.SetNDC()
          tx.SetTextFont(42)
          tx.SetTextSize(0.045)
          tx.DrawLatex(0.7, 0.8, "RMS : %.3f "%(hdummy.GetRMS()))

          
          c2.SaveAs(outdir+'%s.png'%c2.GetName())
          c2.SaveAs(outdir+'%s.pdf'%c2.GetName())
          c2.SaveAs(outdir+'%s.C'%c2.GetName())


 
      c.cd()
      ROOT.gStyle.SetOptStat(0)
      ROOT.gStyle.SetOptFit(0)
 
      #g.SetMarkerStyle(20+iv)
      g_mm[mod].SetMarkerStyle(20+iv)
      g_mm[mod].SetMarkerColor(colors[mod])
      g_mm[mod].SetLineColor(colors[mod])
      if '2E14' in mod: g_mm[mod].SetMarkerStyle(22)
      g_mm[mod].SetMarkerSize(1)
      if (g_mm[mod].GetMarkerStyle() == 22): g_mm[mod].SetMarkerSize(1.25)
      g_mm[mod].Draw('psame')

      ovEff = vov
      if ('2E14' in mod or '1E14' in mod or '1E13' in mod):
         ovEff = getVovEffDCR(data, mod, ('%.02f'%vov))[0]
      #leg.AddEntry(g_mm[mod], '%s, V_{OV} = %.2f V'%(labels[mod],ovEff), 'PL') 
      leg.AddEntry(g_mm[mod], '%s'%(labels[mod]), 'PL') 

   leg.Draw()
   
   
   #latex = ROOT.TLatex(0.20,0.18,'HPK, 25 #mum')
   #latex.SetNDC()
   #latex.SetTextSize(0.045)
   #latex.SetTextFont(42)
   #latex.Draw()

   #cms_logo = draw_logo()
   #cms_logo.Draw()

leg.Draw()
c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())

c.cd()
#### with a pol0 fit to check consistency of tRefs numbers
for i, mod in enumerate(modules):

    pol0= ROOT.TF1(mod, "pol0", -0.5*barConversionFact, 15.5*barConversionFact)

    g_mm[mod].Fit(pol0, "R")
    pol0.SetLineColor(g_mm[mod].GetLineColor())
    t1 = ROOT.TLatex()
    t1.SetNDC()
    t1.SetTextFont(42)
    t1.SetTextSize(0.045)
    t1.DrawLatex(0.17, 0.32 + i * 0.2,"%s: %.2f #pm %.2f"%(mod, pol0.GetParameter(0), pol0.GetParError(0)))
    pol0.Draw("same")


c.SaveAs(outdir+'%s_fit.png'%c.GetName())
c.SaveAs(outdir+'%s_fit.pdf'%c.GetName())
c.SaveAs(outdir+'%s_fit.C'%c.GetName())



#hPad.Delete()

   
