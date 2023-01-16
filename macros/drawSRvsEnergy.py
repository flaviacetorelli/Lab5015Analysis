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


def findMinEnergy(h_energy):

  fitLandau = ROOT.TF1("fitLandau","landau", 0, 1000)
  h_energy.GetXaxis().SetRangeUser(50,800)
  maxbin = h_energy.GetMaximumBin()
  peak = h_energy.GetBinCenter(maxbin)
  fitLandau.SetRange(peak*0.8, peak*1.2)
  fitLandau.SetParameter(1,peak)
  fitLandau.SetParameter(2,0.1*peak)
  h_energy.Fit("fitLandau","QR")
  energyMin = fitLandau.GetParameter(1)*0.75
  return energyMin

def dofinalPlot(gR, gL, label, outdir):

  ROOT.gStyle.SetOptFit(0)
  c = ROOT.TCanvas()
  haxis = ROOT.TH1F("", "", 18, -0.9, 15.9)
  haxis.GetXaxis().SetTitle("bar")
  haxis.GetYaxis().SetTitle("#sigma(Ratio)")
  
  c.cd()
  c.SetGrid()
  
  leg = ROOT.TLegend(0.2, 0.7, 0.4, 0.9)
  leg.AddEntry(gL, "Ratio L ", "P")
  leg.AddEntry(gR, "Ratio R ", "P")
  
  
  haxis.Draw("")
  haxis.GetYaxis().SetRangeUser(0.0,0.5)
  gR.Draw("P same")
  gL.Draw("P same")
  
  leg.Draw("same")


  myfitL = ROOT.TF1("myfitL","pol0", 0, 100)
  myfitL.SetLineColor(gL.GetLineColor())
  myfitL.SetParameter(0, gL.GetMean(2))
  gL.Fit("myfitL")
  latexL = ROOT.TLatex(0.50,0.85, "<#sigma> = %.3f #pm %.3f" % ( myfitL.GetParameter(0), myfitL.GetParError(0)))
  latexL.SetNDC()
  latexL.SetTextFont(42)
  latexL.SetTextSize(0.04)
  latexL.Draw("same")
  myfitR = ROOT.TF1("myfitR","pol0", 0, 100)
  myfitR.SetLineColor(gR.GetLineColor())
  myfitR.SetLineStyle(2)
  myfitR.SetParameter(0, gR.GetMean(2))
  gR.Fit("myfitR")
  latexR = ROOT.TLatex(0.50,0.72, "<#sigma> = %.3f #pm %.3f" % ( myfitR.GetParameter(0), myfitR.GetParError(0)))
  latexR.SetNDC()
  latexR.SetTextFont(42)
  latexR.SetTextSize(0.04)
  latexR.Draw("same")


  c.Print(outdir + "c_"+label+"_vs_bar.png")

#parse arguments
parser = parser = argparse.ArgumentParser()
parser.add_argument("--run",       action="store",      type=int,                         help="run number")  
parser.add_argument("--th",        action="store",      type=str,                         help="first thr value in ADC")  
parser.add_argument("--rms",       action="store",      type=float,   default=2.,         help="rms to perform fit")
parser.add_argument("--label",     action="store",      type=str,                         help="label to define output folder")
args = parser.parse_args()

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
#run = 5297
#run = 5301
#run = 5309
#run = 5196

#FIXME choose fit ranges
#fact = 2. #for the ranges of fit, how many sigma do u want?

bars = {
5301 : ['00', '01', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14','15'],
5309 : ['00', '01', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14','15'],
5195 : ['00', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14'],
5196 : ['00', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14'],
5297 : [ '01', '02', '03', '04', '05','06', '07', '08','09','10','11','12','13','14']
}


#energy cut to reject some noise events
energyCut = {
5301 : energyCut5301, 
5309 : energyCut5309, 
5297 : energyCut5297,
}

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


run = args.run
th = args.th
fact = args.rms

print "Analyzing run: " , run , " at ith: ", th
#file0 =ROOT.TFile.Open("/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/moduleCharacterization_step1_run"+str(int(run))+".root"
#outdir = "/eos/user/f/fcetorel/www/MTD/TBjune22/run"+str(int(run))+"/correlationPlot/tBt/nooffset_prova/"
file0 = ROOT.TFile.Open("/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/moduleCharacterization_step1_TbToffFit_run"+str(int(run))+".root")
#print file0
#outdir = "/eos/user/f/fcetorel/www/MTD/TBjune22/run"+str(int(run))+"/correlationPlot/SRvsEn/offset_nocut/"
outdir = "/eos/user/f/fcetorel/www/MTD/TBjune22/run"+str(int(run))+"/correlationPlot/SRvsEn/th"+th+"/offset_enCut_final/%s/"%(args.label)
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system("cp /eos/user/f/fcetorel/www/index.php "+outdir)

gLinL = ROOT.TGraphErrors()
gLinR = ROOT.TGraphErrors()
gOverL = ROOT.TGraphErrors()
gOverR = ROOT.TGraphErrors()

label = ['L','R']
iR = 0
iL = 0
ibar = 0


for bar in bars[run]:
  ibar = int(bar)
  data = file0.Get("data_bar"+bar+"L-R_Vov"+vov[run]+"_th"+th)
  #print "data_bar"+bar+"L-R_Vov"+vov[run]+"_th05"    

  h_SRratio_vs_Eratio = ROOT.TH2F("h_SRratio_vs_Eratio_bar"+bar,"h_SRratio_vs_Eratio_bar"+bar,100,0,3,100,0,3)
  h_SRratio_vs_Eratio.GetXaxis().SetTitle("energyL/energyR")
  h_SRratio_vs_Eratio.GetYaxis().SetTitle("SR L / SR R")

  # SRratio vs energyRatio
  data.Draw("totR/totL:energyL/energyR >> h_SRratio_vs_Eratio_bar"+bar,"(energyL + energyR )/2 > "+energyCut[run][bar]+" && (energyL + energyR )/2 < 950","GOFF")
  SR = data.GetVal(0)
  energy = data.GetVal(1)
  nent = h_SRratio_vs_Eratio.GetEntries()

  c1 =ROOT.TCanvas()
  c1.Clear()
  c1.SetGridx()
  c1.SetGridy()
  c1.SetLogz()
  h_SRratio_vs_Eratio.GetZaxis().SetRangeUser(0.1,100)
  h_SRratio_vs_Eratio.Draw("colz")
  h_SRratio_vs_Eratio_pfx = h_SRratio_vs_Eratio.ProfileX()
  h_SRratio_vs_Eratio_pfx.Draw("same")
  
                                                     
  c1.Print(outdir + "c_bar"+bar+"_SRratio_vs_Eratio_"+mod[run]+".png")
 


  for lab in label: 
    h_En = ROOT.TH1F("h_En_bar"+bar+lab,"h_En_bar"+bar+lab,100,0,900)
    h_En.GetXaxis().SetTitle("Energy ")                 
    
    h_ratio_lin = ROOT.TH1F("h_ratio_lin_bar"+bar+lab,"h_ratio_lin_bar"+bar+lab,100,0,2)
    h_ratio_over = ROOT.TH1F("h_ratio_over_bar"+bar+lab,"h_ratio_over_bar"+bar+lab,100,0,2)
    
    # 1/ energy
    c1.Clear()
    #data.Draw("energy"+lab+" >> h_En_bar"+bar+lab,"(energyL + energyR )/2 > "+energyCut[run][bar]+" && (energyL + energyR )/2 < 950","GOFF")
    data.Draw("energy"+lab+" >> h_En_bar"+bar+lab,"(energyL + energyR )/2 < 950","GOFF")
 
    print (energyCut[run][bar] , findMinEnergy(h_En))                    
    minEnergy = str(findMinEnergy(h_En))
 
    h_En.Draw("")
    center = h_En.GetMean()  # using this to set the range of linear fit and to set x axis of histo
    rms = h_En.GetRMS()  # using this to set the range of linear fit
    c1.Print(outdir + "c_bar"+bar+lab+"_Energy_"+mod[run]+".png")
 

    h_SR_vs_overenergy = ROOT.TH2F("h_SR_vs_overenergy_bar"+bar+lab,"h_SR_vs_overenergy_bar"+bar+lab,100,0,1/center + 3 * (rms/(center*center)),100,0,30)
    h_SR_vs_overenergy.GetXaxis().SetTitle("1 / energy [a.u.]")
    h_SR_vs_overenergy.GetYaxis().SetTitle("SR [#muA / ns]")
  
    h_SR_vs_energy = ROOT.TH2F("h_SR_vs_energy_bar"+bar+lab,"h_SR_vs_overenergy_bar"+bar,100,0,center + 3 *rms,100,0,30)
    h_SR_vs_energy.GetXaxis().SetTitle("energy [a.u.]")
    h_SR_vs_energy.GetYaxis().SetTitle("SR  [#muA / ns]")
    
    h_SR_vs_t1fine = ROOT.TH2F("h_SR_vs_t1fine_bar"+bar+lab,"h_SR_vs_t1fine_bar"+bar,50,0,1000,100,0,30)
    h_SR_vs_t1fine.GetXaxis().SetTitle("t1fine [a.u.]")
    h_SR_vs_t1fine.GetYaxis().SetTitle("SR  [#muA / ns]")
 

 
    # SR vs energy
    #data.Draw("8*0.313/tot"+lab+" : energy"+lab+" >> h_SR_vs_energy_bar"+bar+lab,"(energyL + energyR )/2 > "+energyCut[run][bar]+" && (energyL + energyR )/2 < 950","GOFF")
    data.Draw("8*0.313/tot"+lab+" : energy"+lab+" >> h_SR_vs_energy_bar"+bar+lab,"(energyL + energyR )/2 > "+minEnergy+" && (energyL + energyR )/2 < 950","GOFF")

    #energy = data.GetVal(1)

    c1.Clear()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogz()
    h_SR_vs_energy.GetZaxis().SetRangeUser(0.1,100)
    h_SR_vs_energy.Draw("colz")
    h_SR_vs_energy_pfx = h_SR_vs_energy.ProfileX()
    h_SR_vs_energy_pfx.Draw("same")
    
    lin = ROOT.TF1("lin","pol1", center-(fact)*rms,center+fact*rms)
    lin.SetLineColor(2)
    h_SR_vs_energy_pfx.Fit(lin, "QR")
    
    q = lin.GetParameter(0)
    m = lin.GetParameter(1)
    
    lin. Draw("same")
                                                        
    c1.Print(outdir + "c_bar"+bar+lab+"_SR_vs_energy_"+mod[run]+".png")

    # SR vs PHASE
    #data.Draw("8*0.313/tot"+lab+" : 1/energy"+lab+" >> h_SR_vs_overenergy_bar"+bar+lab,"(energyL + energyR )/2 > "+energyCut[run][bar]+" && (energyL + energyR )/2 < 950","GOFF")
    data.Draw("8*0.313/tot"+lab+" : t1fine"+lab+" >> h_SR_vs_t1fine_bar"+bar+lab,"(energyL + energyR )/2 > "+minEnergy+" && (energyL + energyR )/2 < 950","GOFF")
    #SR = data.GetVal(0)
    #t1fine = data.GetVal(1)
    #nent = h_SR_vs_t1fine.GetEntries()
    
    c1.Clear()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogz()
    h_SR_vs_t1fine.GetZaxis().SetRangeUser(0.1,100)
    h_SR_vs_t1fine.Draw("colz")
    h_SR_vs_t1fine_pfx = h_SR_vs_t1fine.ProfileX()
    h_SR_vs_t1fine_pfx.Draw("same")
                                                       
    c1.Print(outdir + "c_bar"+bar+lab+"_SR_vs_t1fine_"+mod[run]+".png")
  
  
    # SR vs 1 over energy
    #data.Draw("8*0.313/tot"+lab+" : 1/energy"+lab+" >> h_SR_vs_overenergy_bar"+bar+lab,"(energyL + energyR )/2 > "+energyCut[run][bar]+" && (energyL + energyR )/2 < 950","GOFF")
    data.Draw("8*0.313/tot"+lab+" : 1/energy"+lab+" >> h_SR_vs_overenergy_bar"+bar+lab,"(energyL + energyR )/2 > "+minEnergy+" && (energyL + energyR )/2 < 950","GOFF")
    SR = data.GetVal(0)
    overenergy = data.GetVal(1)
    nent = h_SR_vs_overenergy.GetEntries()
    
    c1.Clear()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetLogz()
    h_SR_vs_overenergy.GetZaxis().SetRangeUser(0.1,100)
    h_SR_vs_overenergy.Draw("colz")
    h_SR_vs_overenergy_pfx = h_SR_vs_overenergy.ProfileX()
    h_SR_vs_overenergy_pfx.Draw("same")
    
    #lin = ROOT.TF1("lin","pol1", center-fact*rms,center+fact*rms)
    #print center, "   " , rms
    #print 1/center, "   " , rms/ (center*center)
    overx = ROOT.TF1("overx","[0]/x + [1]", 1./center-(fact-1.)*(rms/(center*center)),1./center+fact*(rms/(center*center)))
    overx.SetLineColor(2)
    h_SR_vs_overenergy_pfx.Fit(overx, "QR")
    
    a = overx.GetParameter(0)
    b = overx.GetParameter(1)
    
    overx. Draw("same")
                                                        
    c1.Print(outdir + "c_bar"+bar+lab+"_SR_vs_overenergy_"+mod[run]+".png")
    
    #print nent
    
    for en in range(0,int(nent)):
        #if (totL[en] < center +0.2 and totR[en] < center + 0.2 and totL[en] > center -0.2 and totR[en] > center -0.2):
        #if (overenergy[en] > 1./center - fact*rms and overenergy[en] < 1./center + fact*rms ):
          #if (m * energy[en] + q > 0 ): ratio = SR[en] / (m * energy[en] + q )  
        if (a / overenergy[en] + b  != 0 ): ratio = SR[en] / (a / overenergy[en] + b )  
        else: ratio = -9999
        h_ratio_over.Fill (ratio)


        #if (1./overenergy[en] > center - fact*rms and 1./overenergy[en] < center + fact*rms ):
        if (m * 1/overenergy[en] + q != 0 ): ratio = SR[en] / (m * 1/overenergy[en] + q )  
        else: ratio = -9999
        h_ratio_lin.Fill (ratio)
        #print SR[en]
    
    ch = ROOT.TCanvas()
    ch.cd()
     
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(32)
    t.SetTextColor(1)
    t.SetTextSize(0.05)
    t.SetTextAlign(12)
    
    ch.Clear()
    h_ratio_over.GetXaxis().SetTitle(" ratio ")
    h_ratio_over.Draw("histo")
    g = ROOT.TF1 ("g", "gaus", h_ratio_over.GetMean() - 2 * h_ratio_over.GetRMS(),h_ratio_over.GetMean() + 2 * h_ratio_over.GetRMS())
    h_ratio_over.Fit(g,"RQ")
    g.Draw("same")
    mean1 = g.GetParameter(1)
    sigma1 = g.GetParameter(2)
    errsigma1 = g.GetParError(2)
    t.DrawLatex(0.2, 0.5, '#sigma = %.3f #pm %.3f'%(sigma1 ,errsigma1))
    
  
    ch.Print(outdir + "c_bar"+bar+lab+"_Res_"+mod[run]+".png")
    
   
    ch.Clear()
     
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(32)
    t.SetTextColor(1)
    t.SetTextSize(0.05)
    t.SetTextAlign(12)
    
    h_ratio_lin.GetXaxis().SetTitle(" ratio ")
    h_ratio_lin.Draw("histo")
    g = ROOT.TF1 ("g1", "gaus", h_ratio_lin.GetMean() - 2 * h_ratio_lin.GetRMS(),h_ratio_lin.GetMean() + 2 * h_ratio_lin.GetRMS())
    h_ratio_lin.Fit(g,"RQ")
    g.Draw("same")
    mean2 = g.GetParameter(1)
    sigma2 = g.GetParameter(2)
    errsigma2 = g.GetParError(2)
    t.DrawLatex(0.2, 0.5, '#sigma = %.3f #pm %.3f'%(sigma2 ,errsigma2))
    #print "I'm channel " , bar, lab, " that's my ch2 reduced: ", g.GetChisquare() / g.GetNDF() 

    # selections to avoid channels with terrible fits
    if mean1 < 0.8 or  mean2 < 0.8 : continue
    if g.GetChisquare() / g.GetNDF() > 2 : continue
     
    if lab == 'L':
      gOverL.SetPoint(gOverL.GetN(), ibar, sigma1)  
      gOverL.SetPointError(gOverL.GetN()-1, 0, errsigma1)  
      gLinL.SetPoint(gLinL.GetN(), ibar, sigma2)  
      gLinL.SetPointError(gLinL.GetN()-1, 0, errsigma2)  
    elif lab == 'R':
      gOverR.SetPoint(gOverR.GetN(), ibar, sigma1)  
      gOverR.SetPointError(gOverR.GetN()-1, 0, errsigma1)  
      gLinR.SetPoint(gLinR.GetN(), ibar, sigma2)  
      gLinR.SetPointError(gLinR.GetN()-1, 0, errsigma2)  
 
    ch.Print(outdir + "c_bar"+bar+lab+"_ResLin_"+mod[run]+".png")
 
 
  
cfinal = ROOT.TCanvas()
haxis = ROOT.TH1F("", "", 18, -0.9, 15.9)
haxis.GetXaxis().SetTitle("bar")
haxis.GetYaxis().SetTitle("#sigma ratio")

cfinal.cd()
cfinal.SetGrid()
gLinL.SetMarkerStyle(21)  
gOverL.SetMarkerStyle(21) 
gLinL.SetMarkerColor(1)  
gOverL.SetMarkerColor(2)

gLinR.SetMarkerStyle(24)  
gOverR.SetMarkerStyle(24) 
gLinR.SetMarkerColor(1)  
gOverR.SetMarkerColor(2)



leg = ROOT.TLegend(0.2, 0.7, 0.4, 0.9)
leg.AddEntry(gLinL, "Ratio L SR / m*E+q  ", "epl")
leg.AddEntry(gOverL, "Ratio L SR / a / invE + b", "epl")
leg.AddEntry(gLinR, "Ratio R SR / m*E+q  ", "epl")
leg.AddEntry(gOverR, "Ratio R SR / a / invE + b", "epl")


haxis.Draw("")
haxis.GetYaxis().SetRangeUser(0.0,0.5)
#if run == 5301 : haxis.GetYaxis().SetRangeUser(0.02,0.09)
#else: haxis.GetYaxis().SetRangeUser(0.0,0.07)
gLinR.Draw("P same")
gOverR.Draw("P same") 
gLinL.Draw("P same")
gOverL.Draw("P same") 


leg.Draw("same")
cfinal.Print(outdir + "c_Comp_vs_bar.png")



dofinalPlot(gLinR, gLinL, "linear", outdir)
dofinalPlot(gOverR, gOverL, "over", outdir)
