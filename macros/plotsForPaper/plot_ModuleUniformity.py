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

from VovsEff import *
# Import file with VovEff and DCR
with open('/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2023/VovsEff_TOFHIR2C.json', 'r') as f:
   data = json.load(f)


inputdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/'
outdir   = '/eos/user/f/fcetorel/www/MTD/plot4BTLpaper/moduleUniformity/May24/'

#modules = {'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C', 'HPK_2E14_LYSO100056_T-35C'}
#modules = {'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C','HPK_nonIrr_C25_LYSO818_Vov3.50_T5C', 'HPK_2E14_LYSO100056_T-35C'}
modules = {'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C', 'HPK_2E14_LYSO100056_T-35C'}

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

c = ROOT.TCanvas('c_timeResolution_vs_bar','c_timeResolution_vs_bar', 600, 500)
hPad = ROOT.TH2F('hPad','', 100, -0.5, 15.5, 100, 0, 120)
hPad.SetTitle("; bar; time resolution [ps]")
hPad.Draw()
c.SetGridy()
c.SetTicks()

leg = ROOT.TLegend(0.20, 0.74, 0.50, 0.87)
leg.SetBorderSize(0)
leg.SetFillStyle(0)



for mod in modules:
   f = ROOT.TFile.Open(fnames[mod])
   
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
      
      if (vov == bestVovs[mod]): 
         h[mod] = ROOT. TH1F('h_%s'%mod,'h_%s'%mod, 40, -0.8,0.8 )
         h[mod].SetLineColor(g.GetLineColor())
         h[mod].SetFillColorAlpha(g.GetLineColor(),0.2)
         for i in range(0,g.GetN()):
            x = (g.GetPointY(i) - g.GetMean(2) )/g.GetMean(2)
            h[mod].Fill(x)
            h_all.Fill(x)
            if ('nonIrr' not in mod):
               h_irr.Fill(x)

      #g.SetMarkerStyle(20+iv)
      g.SetMarkerStyle(20+iv)
      if '2E14' in mod: g.SetMarkerStyle(22)
      g.SetMarkerSize(1)
      if (g.GetMarkerStyle() == 22): g.SetMarkerSize(1.25)
      g.Draw('psame')

      ovEff = vov
      if ('2E14' in mod or '1E14' in mod or '1E13' in mod):
         ovEff = getVovEffDCR(data, mod, ('%.02f'%vov))[0]
      leg.AddEntry(g, '%s, V_{OV} = %.2f V'%(labels[mod],ovEff), 'PL') 

   leg.Draw()
   
   
   latex = ROOT.TLatex(0.20,0.20,'HPK, 25 #mum')
   latex.SetNDC()
   latex.SetTextSize(0.050)
   latex.SetTextFont(42)
   latex.Draw()

   cms_logo = draw_logo()
   cms_logo.Draw()

#leg.Draw()
#c.SaveAs(outdir+'%s.png'%c.GetName())
#c.SaveAs(outdir+'%s.pdf'%c.GetName())
#c.SaveAs(outdir+'%s.C'%c.GetName())
#hPad.Delete()

leg.Draw()
c.SaveAs(outdir+'%s_v1.png'%c.GetName())
c.SaveAs(outdir+'%s_v1.pdf'%c.GetName())
c.SaveAs(outdir+'%s_v1.C'%c.GetName())
hPad.Delete()
   
