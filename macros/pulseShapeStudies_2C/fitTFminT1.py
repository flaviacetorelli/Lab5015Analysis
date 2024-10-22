#! /usr/bin/env python3
#macro to fit t2-t1 ad tF-t1 so that they are all defined such as t2 > t1 and tF > t1 in all Npes range
#parameters of the fit were tuned by hand...
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
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.07,'X')
ROOT.gStyle.SetTitleSize(0.07,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadRightMargin(0.09)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
def myfunc(x, par):
        xx = x[0]
        if( xx <= par[0] ):
                return par[4]*pow(xx,3) + par[3]*pow(xx,2) + par[2]*xx + par[1]
        else:
                return par[5] * xx +  par[6]
                #return par[5]*xx + par[4]*pow(par[0],3) + par[3]*pow(par[0],2) + par[2]*par[1] - par[5]*par[0] + par[1]


colors ={
        15: ROOT.kAzure,
        19: ROOT.kOrange,
        20: ROOT.kGreen +1,
        23: ROOT.kPink,
        24: ROOT.kRed+1,
        28: ROOT.kBlue,
        }
outdir = "/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/laserMeasurements_LYSO818_HPK25/tF-t1Fits/"
inFileSR = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')
c = ROOT.TCanvas()
print ("")
print ("")
print ("---------------- t1 and t2 ------------------------")
c.Clear()

# tdiff = t2-t1 ... where t1 and t2 are at th: ith1 = thRef, ith2 = thRef + delta
#multigraph to store all ovs
#tF
g_tFdiff = {}
tFdiff = {}
fFdiff = {}


#g_t1_vs_gainNpe_Vov3.00_th35
vovs = [0.8, 1.00, 1.50, 2.00, 3.00]
#vovs = [0.8]



outfilenameF = '%s/plots_tFmint1_LYSO818_pulseShapeStudy_vsNpe.root'%outdir
outfileF = ROOT.TFile(outfilenameF, 'RECREATE')



for thRef in [15, 20]:
    print ("")
    print ("-------- thRef = %d"%thRef)

    g_tFdiff [thRef] = ROOT.TMultiGraph()

    legF = ROOT.TLegend(0.56, 0.75, 0.88, 0.92)


    for vov in vovs:
   
        tFdiff[(thRef, vov)] = ROOT.TGraphErrors()

        #getting t1 and t2 graphs from file (for t2 choose ith1+4 and ith1+8)
        gt1 = inFileSR.Get('g_t1_vs_gainNpe_Vov%.2f_th%02d'%(vov,thRef))
        gtF = inFileSR.Get('g_tF_vs_gainNpe_Vov%.2f_th%02d'%(vov,thRef))
    
        print ("----------- This is Vov: %.2f   ------------"%vov)
            
        #gt1.Print()
        #gtF.Print()

        tFdiff[(thRef, vov)].SetMarkerStyle(gt1.GetMarkerStyle())

        for i in range(gt1.GetN()):
            if gt1.GetPointX(i) == gtF.GetPointX(i): 
                tFdiff[(thRef, vov)].SetPoint(tFdiff[(thRef, vov)].GetN(), gtF.GetPointX(i),  gtF.GetPointY(i) - gt1.GetPointY(i) )
                tFdiff[(thRef, vov)].SetPointError(tFdiff[(thRef, vov)].GetN()-1, gtF.GetErrorX(i),  math.sqrt(pow(gtF.GetErrorY(i),2) +pow(gt1.GetErrorY(i),2) ))

            else: print("Skipping point for tF(th1)-t1 graph %.2f entry %d"%( vov, i))



        g_tFdiff[thRef].Add(tFdiff[(thRef, vov)], "p")
        outfileF.cd()
        legF.AddEntry(tFdiff[(thRef, vov)], "V_{ov} = %.2f V"%vov, "pe")
        tFdiff[(thRef, vov)].Write('g_tF-t1_vs_gainNpe_Vov%.2f_th%02d'%(vov,thRef))
#Th 15
#****************************************
#Minimizer is Minuit2 / Migrad
#Chi2                      =      18.1689
#NDf                       =           49
#Edm                       =  3.79427e-06
#NCalls                    =          227
#p0                        =  1.47998e+09   +/-   1.41421
#p1                        =      9.58676   +/-   0.451553
#p2                        = -2.51351e-10   +/-   2.44708e-10
#p3                        =  5.98501e-18   +/-   1.78231e-18
#p4                        =  -3.2352e-27   +/-   1.0819e-27
#p5                        = -7.88332e-10   +/-   2.45698e-11
#Minimizer is Minuit2 / Migrad
#Chi2                      =      65.4516
#NDf                       =           49
#Edm                       =  8.48888e-17
#NCalls                    =          670
#p0                        =        2e+09   +/-   0.00548244      (limited)
#p1                        =      7.93403   +/-   2.67603e-10
#p2                        = -3.78578e-10   +/-   3.36894e-19
#p3                        =  5.42505e-18   +/-   8.79743e-29
#p4                        = -2.27325e-27   +/-   9.27707e-37
#p5                        = -7.32799e-10   +/-   2.34746e-18


    p0 = 2.9E09
    # drawing and fitting tF
    fFdiff[thRef]= ROOT.TF1("fFdiff_th%d"%(thRef),myfunc,0., 12E09, 7)
    # pol1 and pol3 used to find the ad hoc params
    #pol3= ROOT.TF1 ("poli3", "pol3", 0.E09, 3E09)
    #pol1= ROOT.TF1 ("poli1", "pol1", 0, 11.E09)
    #pol1.SetParameters(15, -8.E-10)

    #if thRef == 20:
    #   pol3.SetParameter(0, 7.93403) 
    #   pol3.SetParameter(1, -3.78578e-10) 
    #   pol3.FixParameter(2, 5.42505e-18) 
    #   pol3.SetParameter(3, -2.27325e-27) 


    #if thRef == 15:
    #   pol3.SetParameter(0, 9.58676) 
    #   pol3.SetParameter(1, -2.51351e-10) 
    #   pol3.SetParameter(2, 5.98501e-18) 
    #   pol3.SetParameter(3, -3.2352e-27) 
    #g_tFdiff[thRef].Fit(pol3, "N", "", 0.5E09, 3.5 )
    #g_tFdiff[thRef].Fit(pol3, "N", "", 0.5E09, 3.5 )

    #g_tFdiff[thRef].Fit(pol1, "N", "", 1.E09, 11E09)
    #fFdiff[thRef].SetParameters(1.5E09, pol3.GetParameter(0),  pol3.GetParameter(1), pol3.GetParameter(2), pol3.GetParameter(3), pol1.GetParameter(1), pol1.GetParameter(0))


    if thRef == 20:
       fFdiff[thRef].SetParameter(0, 1.4776E09)  
       fFdiff[thRef].SetParameter(1, 7.93403) 
       fFdiff[thRef].SetParameter(2, -3.78578e-10) 
       fFdiff[thRef].SetParameter(3, 5.42505e-18) 
       fFdiff[thRef].SetParameter(4, -2.27325e-27) 
       fFdiff[thRef].SetParameter(5, -7.32799e-10)
       fFdiff[thRef].SetParameter(6, 12.933)

       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.9E09, 5E09)
       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.9E09, 5E09)
       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.9E09, 5E09)

    if thRef == 15:
       fFdiff[thRef].SetParameter(0, 1.3e+09)  
       fFdiff[thRef].SetParameter(1, 9.58676) 
       fFdiff[thRef].SetParameter(2, -2.51351e-10) 
       fFdiff[thRef].SetParameter(3, 5.98501e-18) 
       fFdiff[thRef].SetParameter(4, -3.2352e-27) 
       fFdiff[thRef].SetParameter(5,-7.61576e-10)
       fFdiff[thRef].SetParameter(6,13.21)
       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.6E09, 12E09)
       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.6E09, 12E09)
       g_tFdiff[thRef].Fit(fFdiff[thRef], "N", "", 0.6E09, 12E09)
    g_tFdiff[thRef].SetTitle("; gain x N_{pe} ; t_{F}-t_{1} [ns] ")
     
    fFdiff[thRef].SetLineStyle(2)
    fFdiff[thRef].SetLineColor(4)
    g_tFdiff[thRef].Draw("psame")
    fFdiff[thRef].Draw("same")
    legF.SetNColumns(2)
    legF.Draw("same")
    c.SaveAs("%s/tF-t1_th1_%02d.png"%(outdir,thRef))
    outfileF.cd()

    g_tFdiff[thRef].Write('g_tF-t1glob_vs_gainNpe_th%02d'%(thRef))
    fFdiff[thRef].Write('f_tF-t1glob_vs_gainNpe_th%02d'%(thRef))


