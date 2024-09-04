#! /usr/bin/env python3
## macro to fit the stochastic term from TB data and obtain alpha and sigma_stoch
## the stochastich term is obtained as the difference in quadrature between total resolution and electronics term
## use the SR vs gainXNpe from laser measruments to estimate the noise term at a certain ov
## very harcoded, have a look at all parameters and outdir and name of output file before running.
import os
import shutil
import glob
import math
import array
import sys
import time

import ROOT
import tdrstyle

from typing import NamedTuple

from SiPM import *
from ROOT import TMath
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)

def sigma_noise(sr):
    noise =  (1/math.sqrt(2)) * math.sqrt( pow ( 420/sr, 2 ) + pow (16.7, 2) ) # 420 as in TB paper
    return noise
e = 1.602e-19

sipm_type = 'HPK-PIT-C25-ES2'
ithMode = 1
label = 'HPK_nonIrr_C25_LYSO818'
#LO_at3p5 = 2390  # 10% less LO for 818
LO_at3p5 = 2390 * 0.9 # 10% less LO for 818
temp = 'T5C'

markers = {
0.80 : 20,
1.00 : 25,
1.50 : 22,
2.00 : 24,
3.00 : 21

}

colors = {}
colors['aldoA'] = ROOT.kRed
colors['aldoB'] = ROOT.kBlue

colors['ch1'] = ROOT.kRed
colors['ave'] = ROOT.kBlack
colors['ch2'] = ROOT.kBlue


inFileSR = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')

inFileTB = ROOT.TFile.Open('/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/HPK_nonIrr_C25_LYSO818_T5C/summaryPlots_HPK_nonIrr_C25_LYSO818_T5C.root')

fit = 'linlog' # u can choose between lin or log for the sr vs gain npe parametrization
plotDir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/compareLaser/'%(label, temp)

thBestFromTB = { # these are taken looking at VS th plots for bar 07 (the one im using now)
        0.60: 5, #5
        0.80: 7,
        1.00: 11,
        1.25: 11,
        1.50: 15,
        2.00: 20,
        3.00: 20, # estimated since no TB data
        3.50: 25,
        }


def stoch_vs_Vov(x, par):
    sipm_type = 'HPK-PIT-C25-ES2'
    LO_at3p5 = 2390*0.9
    #LO_at3p5 = 2390
    xx = x[0]
    Npe = LO_at3p5* (PDE(sipm_type, xx) /PDE(sipm_type, 3.50))* (4.2) * 1.25 
    return par[1] * pow(par[0]/Npe,par[2])




print (inFileTB)

alpha_vs_gain = ROOT.TGraphErrors()
p1_vs_gain = ROOT.TGraphErrors()
g_noise_vs_Vov = ROOT.TGraphErrors()
g_stoch_vs_Vov = ROOT.TGraphErrors()



angleFact = math.cos(49 * math.pi / 180) / math.cos(52*math.pi / 180)
tRes_vs_Vov = inFileTB.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_bar07_enBin01')

#tRes_vs_Vov.Print()
tRes_vs_Vov.Scale(1./angleFact)

#tRes_vs_Vov.Print()


for point in range(0,tRes_vs_Vov.GetN()):
    Vov = tRes_vs_Vov.GetPointX(point)
    gain = Gain(sipm_type, Vov)
    thRef = thBestFromTB[Vov]

    print (Vov, thRef)

    fSR = inFileSR.Get('f_SRglob_vs_gainNpe_linlog_th%02d'%thRef)

    Npe = LO_at3p5* (PDE(sipm_type, Vov) /PDE(sipm_type, 3.50))* (4.2) * 1.25 
    sr = fSR.Eval(gain*Npe)
    errSr = 0.1 # 10 % of error on slew rate
    
    tRes = tRes_vs_Vov.GetPointY(point)
    tRes_err = tRes_vs_Vov.GetErrorY(point)
    
    noise = sigma_noise(sr)
    noise_err = 0.5*(sigma_noise(sr*(1-errSr))-sigma_noise(sr*(1+errSr)))
    
    if tRes > 200 or tRes < 0: continue
    
    stoch = math.sqrt(pow(tRes, 2)- pow(noise, 2))
    stoch_err = 1./stoch*math.sqrt( pow(tRes_err*tRes,2)+pow( sigma_noise(sr)*noise_err ,2) )
    
    print (Vov, Npe, tRes, sr)
    
    g_noise_vs_Vov.SetPoint(g_noise_vs_Vov.GetN(), Vov, noise)
    g_noise_vs_Vov.SetPointError(g_noise_vs_Vov.GetN()-1, 0., noise_err)
    g_stoch_vs_Vov.SetPoint(g_stoch_vs_Vov.GetN(), Vov, stoch)
    g_stoch_vs_Vov.SetPointError(g_stoch_vs_Vov.GetN()-1, 0., stoch_err)

leg = ROOT.TLegend(0.5, 0.4, 0.92, 0.8)
leg.SetTextSize(0.035)

c = ROOT.TCanvas("c_tRes_vs_Vov_StochFit_thBest", "c_tRes_vs_Vov_StochFit_thBest" )
hPad = ROOT.gPad.DrawFrame(0,0.,4.0, 100.)
hPad.SetTitle(";  V_{ov}; time resolution [ps]")


hPad.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

tRes_vs_Vov.SetMarkerColor(ROOT.kBlack)
tRes_vs_Vov.SetLineColor(ROOT.kBlack)

tRes_vs_Vov.Draw("PL same")
leg.AddEntry(tRes_vs_Vov, "TB data", "PL")

g_noise_vs_Vov.SetLineWidth(2)
g_noise_vs_Vov.SetLineColor(ROOT.kBlue)
g_noise_vs_Vov.SetFillColor(ROOT.kBlue)
g_noise_vs_Vov.SetFillColorAlpha(ROOT.kBlue+1,0.5)
g_noise_vs_Vov.SetFillStyle(3004)
g_noise_vs_Vov.Draw('E3lsame')
leg.AddEntry(g_noise_vs_Vov, "noise, sr(gainXNpe) from laser", "L")

g_stoch_vs_Vov.SetMarkerColor(ROOT.kGreen+1)
g_stoch_vs_Vov.SetLineColor(ROOT.kGreen+1)
g_stoch_vs_Vov.Draw("PL same")
leg.AddEntry(g_stoch_vs_Vov, "stoch = sqrt(data^2-noise^2)", "PL")

f_stoch_vs_Vov = ROOT.TF1( "f_stoch_vs_Vov", stoch_vs_Vov ,0.58 , 3.6 , 3 )
f_stoch_vs_Vov.SetParameters(7000, 25.7, 0.5)
f_stoch_vs_Vov.FixParameter(0, 7000)
#f_stoch_vs_Vov.FixParameter(1, 25.7)
#g_stoch_vs_Vov.Fit(f_stoch_vs_Vov, "R" )
g_stoch_vs_Vov.Fit(f_stoch_vs_Vov, "", "", 0.89 , 3.6  )
f_stoch_vs_Vov.SetLineColor(ROOT.kGreen+1)
f_stoch_vs_Vov.Draw("same")




tExp = ROOT.TGraphErrors()
for point in range(0,tRes_vs_Vov.GetN()):
    s_noise = g_noise_vs_Vov.Eval(tRes_vs_Vov.GetPointX(point))
    err_s_noise = g_noise_vs_Vov.GetErrorY(point)

    s_stoch = f_stoch_vs_Vov.Eval( tRes_vs_Vov.GetPointX(point))
    err_s_stoch = 0.

    s_tot = math.sqrt( pow(s_stoch, 2) + pow(s_noise, 2) )
    err_s_tot = 1/s_tot * math.sqrt(math.pow(s_noise*err_s_noise,2) + math.pow(s_stoch*err_s_stoch,2) )

    tExp.SetPoint(tExp.GetN(), tRes_vs_Vov.GetPointX(point), s_tot)
    tExp.SetPointError(tExp.GetN()-1, 0. , err_s_tot)


tExp.SetLineColor(ROOT.kRed+1)
tExp.SetLineWidth(2)
tExp.SetLineColor(ROOT.kRed+1)
tExp.SetFillColor(ROOT.kRed+1)
tExp.SetFillColorAlpha(ROOT.kRed+1,0.5)
tExp.SetFillStyle(3004)
tExp.Draw('E3lsame')


leg.AddEntry(tExp, "noise + stoch (from fit)", "L")
leg.Draw("same")
#c.Print("%s/%s.png"%(plotDir, c.GetName()))
c.Print("%s/%s_LOlow.png"%(plotDir, c.GetName()))

