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
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.052,'X') #0.055 before
ROOT.gStyle.SetLabelSize(0.052,'Y')
ROOT.gStyle.SetTitleSize(0.06,'X') #0.07 before
ROOT.gStyle.SetTitleSize(0.06,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.045)
ROOT.gStyle.SetPadBottomMargin(0.13)
ROOT.gStyle.SetPadTopMargin(0.13)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


colors = {
        'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C' : ROOT.kBlue,  
        'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C' : ROOT.kGreen + 2 , 
        'HPK_2E14_LYSO100056_T-35C' : ROOT.kOrange + 1 , 

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
outdir   = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/moduleUniformity/paper1_Nov24/'

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
          'HPK_2E14_LYSO100056_T-35C' : '2 x 10^{14} 1 MeV n_{eq}/cm^{2}',
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

h = {}
h_all = ROOT.TH1F('h_all','h_all', 40, -0.8,0.8)
h_irr = ROOT.TH1F('h_irr','h_irr', 40, -0.8,0.8)

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
#c.SetTicks()
c.SetTicky(1)
c.SetTickx(0)


leg = ROOT.TLegend(0.20, 0.73, 0.50, 0.85)
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


g_mm = {}
for mod in modules:
   f = ROOT.TFile.Open(fnames[mod])
   g_mm[mod]= ROOT.TGraphErrors()   
#   c = ROOT.TCanvas('c_timeResolution_vs_bar_%s'%mod,'c_timeResolution_vs_bar_%s'%mod, 600, 500)
#   hPad = ROOT.TH2F('hPad','', 100, -0.5, 15.5, 100, 0, 120)
#   hPad.SetTitle("; bar; time resolution [ps]")
#   hPad.Draw()
#   c.SetGridy()
#   c.SetTicks()
#
#   leg = ROOT.TLegend(0.65, 0.66, 0.95, 0.92)
#   leg.SetBorderSize(0)
#   leg.SetFillStyle(0)
#   if (len(Vovs[mod])>4):
#      leg.SetNColumns(2);
#      leg.SetColumnSeparation(0.2);

   for iv,vov in enumerate(Vovs[mod]):
         
      g = f.Get('g_deltaT_totRatioCorr_bestTh_vs_bar_Vov%.2f_enBin01'%vov) 
      print('g_deltaT_totRatioCorr_bestTh_vs_bar_Vov%.2f_enBin01'%vov)
      print(g.GetN())
      #g.Scale(1./enScale)

      for i in range(0,g.GetN()+1): #conversion in mm
          g_mm[mod].SetPoint(g_mm[mod].GetN(), g.GetPointX(i)*barConversionFact, g.GetPointY(i)/enScale) #accounting for angle offset
          g_mm[mod].SetPointError(g_mm[mod].GetN()-1, 0, g.GetErrorY(i)/enScale )


      
      #if (vov == bestVovs[mod]): 
      #   h[mod] = ROOT. TH1F('h_%s'%mod,'h_%s'%mod, 40, -0.8,0.8 )
      #   h[mod].SetLineColor(g.GetLineColor())
      #   h[mod].SetFillColorAlpha(g.GetLineColor(),0.2)
      #   for i in range(0,g.GetN()):
      #      x = (g.GetPointY(i) - g.GetMean(2) )/g.GetMean(2)
      #      h[mod].Fill(x)
      #      h_all.Fill(x)
      #      if ('nonIrr' not in mod):
      #         h_irr.Fill(x)

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
      leg.AddEntry(g_mm[mod], '%s, V_{OV} = %.2f V'%(labels[mod],ovEff), 'PL') 

   leg.Draw()
   
   
   latex = ROOT.TLatex(0.20,0.18,'HPK, 25 #mum')
   latex.SetNDC()
   latex.SetTextSize(0.045)
   latex.SetTextFont(42)
   latex.Draw()

   #cms_logo = draw_logo()
   #cms_logo.Draw()

leg.Draw()
c.SaveAs(outdir+'%s.png'%c.GetName())
c.SaveAs(outdir+'%s.pdf'%c.GetName())
c.SaveAs(outdir+'%s.C'%c.GetName())
hPad.Delete()

  
