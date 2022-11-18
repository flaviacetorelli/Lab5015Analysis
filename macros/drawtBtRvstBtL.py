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

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.06,'X')
ROOT.gStyle.SetTitleSize(0.06,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.0,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning



energyCut5301 = {
"00" :"50",
"01" :"50",
"02" :"20",
"03" :"20",
"04" :"50",
"05" :"20",
"06" :"10",
"07" :"10",
"08" :"100",
"09" :"90",
"10" :"20",
"11" :"20",
"12" :"50" ,
"13" :"100" ,
"14" :"50",
"15" :"50",
}

energyCut5309 = {
"00" :"50",
"01" :"50",
"02" :"20",
"03" :"20",
"04" :"50",
"05" :"20",
"06" :"10",
"07" :"20",
"08" :"100",
"09" :"90",
"10" :"20",
"11" :"20",
"12" :"20" ,
"13" :"100" ,
"14" :"20",
"15" :"20",
}

energyCut5297 = {
"00" :"0",
"01" :"100",
"02" :"10",
"03" :"100",
"04" :"100",
"05" :"80",
"06" :"80",
"07" :"100",
"08" :"100",
"09" :"100",
"10" :"50",
"11" :"50",
"12" :"50" ,
"13" :"100" ,
"14" :"100",
"15" :"0",
}

#FIXME choose run
run = 5297
#run = 5196

#FIXME choose fit ranges
fact = 2. #for the ranges of fit, how many sigma do u want?

if run == 5301 or run == 5309: bars = ['00', '01', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14','15']
elif run == 5196 or run == 5195: bars = ['00', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14']
elif run == 5297: bars =  [ '01', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14']

if run == 5301: energyCut = energyCut5301 
elif run == 5309: energyCut = energyCut5309
elif run == 5297: energyCut = energyCut5297


th = "09"
vov = {
5301 : '1.80',
5309 : '1.60',
5297 : '1.50',
5196 : '1.50',
5195 : '3.50'
}
mod = {
5301 : 'HPK_2E14',
5309 : 'HPK_2E14',
5297 : 'HPK_nonirr',
5196 : 'HPK_nonirr',
5195 : 'HPK_nonirr'
}

#file0 =ROOT.TFile.Open("/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/moduleCharacterization_step1_run"+str(int(run))+".root"
#outdir = "/eos/user/f/fcetorel/www/MTD/TBjune22/run"+str(int(run))+"/correlationPlot/tBt/nooffset_prova/"
file0 =ROOT.TFile.Open("/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/moduleCharacterization_step1_TbToffFit_run"+str(int(run))+".root")
outdir = "/eos/user/f/fcetorel/www/MTD/TBjune22/run"+str(int(run))+"/correlationPlot/tBt/th"+th+"/offset_enCut_final/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system("cp /eos/user/f/fcetorel/www/index.php "+outdir)
#outdir = "./plotti/"
#print file0
#moduleCharacterization_step1_TbToffFit_run5301.root
gRes = ROOT.TGraphErrors()
gResDiff = ROOT.TGraphErrors()
gResDiffCut = ROOT.TGraphErrors()

ibar = 0
for bar in bars:
  ibar = int(bar)
  data = file0.Get("data_bar"+bar+"L-R_Vov"+vov[run]+"_th"+th)
  #print "data_bar"+bar+"L-R_Vov"+vov[run]+"_th05"    
  
  h_tbtL = ROOT.TH1F("h_tbtL_bar"+bar,"h_tbtL_bar"+bar,100,0,1.5)
  h_tbtDiff_cutt1fine = ROOT.TH1F("h_tbtDiff_cutt1fine_bar"+bar,"h_tbtDiff_cutt1fine_bar"+bar,100,-1,1)
  h_tbtR_vs_tbtL = ROOT.TH2F("h_tbtR_vs_tbtL_bar"+bar,"h_tbtR_vs_tbtL_bar"+bar,100,0,1.5,100,0,1.5)
  h_tbtDiff_vs_phase = ROOT.TH2F("h_tbtDiff_vs_phase_bar"+bar,"h_tbtR_vs_tbtL_bar"+bar,100,0,1000, 100,-1,1)
  h_tbtR_vs_tbtL.GetXaxis().SetTitle("TbT L [ns]")
  h_tbtR_vs_tbtL.GetYaxis().SetTitle("TbT R [ns]")
  h_tbtDiff_vs_phase.GetXaxis().SetTitle("t1fine Mean ")
  h_tbtDiff_vs_phase.GetYaxis().SetTitle("TbT Diff [ns]")
  
  h_tbtL.GetXaxis().SetTitle("TbT L [ns]")
  
                   
  h_dist = ROOT.TH1F("h_dist_bar"+bar,"h_dist_bar"+bar,100,-1,1)
  h_diff = ROOT.TH1F("h_diff_bar"+bar,"h_diff_bar"+bar,100,-1,1)
                      
  
  c1 = ROOT.TCanvas()
  # tBtL 
  c1.Clear()
  data.Draw("totL >> h_tbtL_bar"+bar,"(energyL + energyR )/2 > "+energyCut[bar]+" && (energyL + energyR )/2 < 950","GOFF")
   
  h_tbtL.Draw("")
  center = h_tbtL.GetMean()  # using this to set the range of linear fit
  rms = h_tbtL.GetRMS()  # using this to set the range of linear fit
  c1.Print(outdir + "c_bar"+bar+"_tbtL_"+mod[run]+".png")
 



  # tBt diff  vs phase plot
  c1.Clear()
  data.Draw("totR-totL:(t1fineR+t1fineL)/2 >> h_tbtDiff_vs_phase_bar"+bar,"(energyL + energyR )/2 > "+energyCut[bar]+" && (energyL + energyR )/2 < 950","GOFF")
  c1.SetGridx()
  c1.SetGridy()
  c1.SetLogz()
  h_tbtDiff_vs_phase.GetZaxis().SetRangeUser(0.1,100)
  h_tbtDiff_vs_phase.Draw("colz")
  h_tbtDiff_vs_phase_pfx = h_tbtDiff_vs_phase.ProfileX()
  h_tbtDiff_vs_phase_pfx.Draw("same")
  
  c1.Print(outdir + "c_bar"+bar+"_tbtDiff_vs_phase_"+mod[run]+".png")
  
  # tot R vs tot L
  data.Draw("totR:totL >> h_tbtR_vs_tbtL_bar"+bar,"(energyL + energyR )/2 > "+energyCut[bar]+" && (energyL + energyR )/2 < 950","GOFF")
  totR = data.GetVal(0)
  totL = data.GetVal(1)
  nent = h_tbtR_vs_tbtL.GetEntries()
  
  c1.Clear()
  c1.SetGridx()
  c1.SetGridy()
  c1.SetLogz()
  h_tbtR_vs_tbtL.GetZaxis().SetRangeUser(0.1,100)
  h_tbtR_vs_tbtL.Draw("colz")
  h_tbtR_vs_tbtL_pfx = h_tbtR_vs_tbtL.ProfileX()
  h_tbtR_vs_tbtL_pfx.Draw("same")
  
  lin = ROOT.TF1("lin","pol1", center-fact*rms,center+fact*rms)
  lin.SetLineColor(2)
  h_tbtR_vs_tbtL_pfx.Fit(lin, "QR")
  
  q = lin.GetParameter(0)
  m = lin.GetParameter(1)
  
  lin. Draw("same")
                                                      
  c1.Print(outdir + "c_bar"+bar+"_tbtR_vs_tbtL_"+mod[run]+".png")
  
  
  for en in range(0,int(nent)):
      #if (totL[en] < center +0.2 and totR[en] < center + 0.2 and totL[en] > center -0.2 and totR[en] > center -0.2):
      if (totL[en] < center + fact*rms and totR[en] < center + fact*rms and totL[en] > center - fact*rms and totR[en] > center - fact*rms):
          dist = totR[en] - (m * totL[en] + q )  / ROOT.TMath.Sqrt(1 + m*m) 
          h_dist.Fill (dist)
          h_diff.Fill( totL[en] - totR[en])
  
  # reso from tBtR vs tBtL 
  ch = ROOT.TCanvas()
  ch.cd()
   
  t = ROOT.TLatex()
  t.SetNDC()
  t.SetTextFont(32)
  t.SetTextColor(1)
  t.SetTextSize(0.05)
  t.SetTextAlign(12)
  
  ch.Clear()
  h_dist.GetXaxis().SetTitle("tBt dist [ns]")
  h_dist.Draw("histo")
  g = ROOT.TF1 ("g", "gaus", h_dist.GetMean() - 2 * h_dist.GetRMS(),h_dist.GetMean() + 2 * h_dist.GetRMS())
  h_dist.Fit(g,"RQ")
  g.Draw("same")
  sigma = g.GetParameter(2)/ROOT.TMath.Sqrt(2)
  errsigma = g.GetParError(2)/ROOT.TMath.Sqrt(2)
  t.DrawLatex(0.2, 0.5, '#sigma_{t} = %.3f #pm %.3f'%(sigma ,errsigma))
  

  gRes.SetPoint(gRes.GetN(), ibar, sigma)  
  gRes.SetPointError(gRes.GetN()-1, 0, errsigma)  
  ch.Print(outdir + "c_bar"+bar+"_Res_"+mod[run]+".png")
  
 
  
  #reso from tDiff
  ch.Clear()
  h_diff.GetXaxis().SetTitle("tBtL - tBtR [ns]")
  h_diff.Draw("histo")
  g = ROOT.TF1 ("g2", "gaus", h_diff.GetMean() - 2 * h_diff.GetRMS(),h_diff.GetMean() + 2 * h_diff.GetRMS())
  h_diff.Fit(g,"RQ")
  g.Draw("same")
  sigma = g.GetParameter(2)/ROOT.TMath.Sqrt(2)
  errsigma = g.GetParError(2)/ROOT.TMath.Sqrt(2)
  t.DrawLatex(0.2, 0.5, '#sigma_{t} = %.3f #pm %.3f'%(sigma ,errsigma))
  gResDiff.SetPoint(gResDiff.GetN(), ibar, sigma)  
  gResDiff.SetPointError(gResDiff.GetN()-1, 0, errsigma)  
 

  ch.Print(outdir + "c_bar"+bar+"_tDiffRes_"+mod[run]+".png")

  #reso from tDiff cut on t1fine
  
  data.Draw("totR-totL >> h_tbtDiff_cutt1fine_bar"+bar,"(energyL + energyR )/2 > "+energyCut[bar]+" && (energyL + energyR )/2 < 950 && (t1fineL + t1fineR )/2 > 550 && (t1fineL + t1fineR )/2 < 600","GOFF")
 
  ch.Clear()
  h_tbtDiff_cutt1fine.GetXaxis().SetTitle("tBtL - tBtR [ns]")
  h_tbtDiff_cutt1fine.Draw("histo")
  g = ROOT.TF1 ("g3", "gaus",  h_tbtDiff_cutt1fine.GetMean() - 2 *  h_tbtDiff_cutt1fine.GetRMS(), h_tbtDiff_cutt1fine.GetMean() + 2 *  h_tbtDiff_cutt1fine.GetRMS())
  h_tbtDiff_cutt1fine.Fit(g,"RQ")
  g.Draw("same")
  sigma = g.GetParameter(2)/ROOT.TMath.Sqrt(2)
  errsigma = g.GetParError(2)/ROOT.TMath.Sqrt(2)
  t.DrawLatex(0.2, 0.5, '#sigma_{t} = %.3f #pm %.3f'%(sigma ,errsigma))
  gResDiffCut.SetPoint(gResDiffCut.GetN(), ibar, sigma)  
  gResDiffCut.SetPointError(gResDiffCut.GetN()-1, 0, errsigma)  
 
  ch.Print(outdir + "c_bar"+bar+"_tDiffRes_cutt1fine550-600_"+mod[run]+".png")


  
  
cfinal = ROOT.TCanvas()
haxis = ROOT.TH1F("", "", 18, -0.9, 15.9)
haxis.GetXaxis().SetTitle("bar")
haxis.GetYaxis().SetTitle("#sigma_{t} [ns]")

cfinal.cd()
cfinal.SetGrid()
gRes.SetMarkerStyle(21)  
gResDiff.SetMarkerStyle(21) 
gRes.SetMarkerColor(1)  
gResDiff.SetMarkerColor(2)
gResDiffCut.SetMarkerColor(3)
leg = ROOT.TLegend(0.2, 0.2, 0.4, 0.3)
leg.AddEntry(gRes, "t dist  ", "epl")
leg.AddEntry(gResDiff, "tdiff", "epl")

haxis.Draw("")
if run == 5301 : haxis.GetYaxis().SetRangeUser(0.02,0.09)
else: haxis.GetYaxis().SetRangeUser(0.0,0.07)
gResDiff.Draw("P same")
gRes.Draw("P same") 

leg.Draw("same")
cfinal.Print(outdir + "c_finalResoComp_vs_bar.png")



cfinal.Clear()
cfinal.SetGrid()
leg.AddEntry(gResDiffCut, "tdiff 550 < tfine < 600", "epl")

haxis.Draw("")
if run == 5301 : haxis.GetYaxis().SetRangeUser(0.02,0.11)
else: haxis.GetYaxis().SetRangeUser(0.0,0.07)
gResDiff.Draw("P same")
gResDiffCut.Draw("P same")
gRes.Draw("P same") 

leg.Draw("same")
cfinal.Print(outdir + "c_finalResoCompAll_vs_bar.png")



 
