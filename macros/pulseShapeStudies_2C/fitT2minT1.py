#! /usr/bin/env python3
#macro to fit t2-t1 ad tF-t1 so that they are all defined such as t2 > t1 and tF > t1 in all Npes range
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
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning



colors ={
        15: ROOT.kAzure,
        19: ROOT.kOrange,
        20: ROOT.kGreen +1,
        23: ROOT.kPink,
        24: ROOT.kRed+1,
        28: ROOT.kBlue,
        }
outdir = "/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/laserMeasurements_LYSO818_HPK25/t2-t1Fits/"
inFileSR = ROOT.TFile.Open('/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_SRt1_LYSO818_pulseShapeStudy_vsNpe_Th05to35.root')
c = ROOT.TCanvas()
print ("")
print ("")
print ("---------------- t1 and t2 ------------------------")
c.Clear()

# tdiff = t2-t1 ... where t1 and t2 are at th: ith1 = thRef, ith2 = thRef + delta
#multigraph to store all ovs
#t2
g_tdiff = {}
tdiff = {}
fdiff = {}


#g_t1_vs_gainNpe_Vov3.00_th35
vovs = [0.8, 1.00, 1.50, 2.00, 3.00]
#vovs = [0.8]
delta = 8 # delta from the first th


#outfilename = '/afs/cern.ch/work/f/fcetorel/public/BTLdigitizationModel_4DPG/plots_t2mint1_LYSO818_pulseShapeStudy_vsNpe_deltaTh%d.root'%delta
outfilename = '%s/plots_t2mint1_LYSO818_pulseShapeStudy_vsNpe_deltaTh%d.root'%(outdir,delta)
outfile = ROOT.TFile(outfilename, 'RECREATE')



for thRef in [15, 20]:
    print ("")
    print ("-------- thRef = %d"%thRef)

    g_tdiff [thRef] = ROOT.TMultiGraph()

    leg = ROOT.TLegend(0.6, 0.75, 0.92, 0.92)

    #getting functions for later initialization
    ft2 = inFileSR.Get('f_t1glob_vs_gainNpe_th%02d'%(thRef+delta))

    for vov in vovs:
   
        tdiff[(thRef, vov)] = ROOT.TGraphErrors()

        #getting t1 and t2 graphs from file (for t2 choose ith1+4 and ith1+8)
        gt1 = inFileSR.Get('g_t1_vs_gainNpe_Vov%.2f_th%02d'%(vov,thRef))
        gt2 = inFileSR.Get('g_t1_vs_gainNpe_Vov%.2f_th%02d'%(vov,thRef+delta))
    
        print ("----------- This is Vov: %.2f   ------------"%vov)
            
        #gt1.Print()
        #gt2.Print()

        tdiff[(thRef, vov)].SetMarkerStyle(gt1.GetMarkerStyle())

        #and now loopig on N entries of each to get t2-t1
        if thRef == 20 and delta == 8 and vov == 0.80: #since in this scenario missing two points at low Npe for t2
            for i in range(gt2.GetN()):
                if gt1.GetPointX(i+2) == gt2.GetPointX(i): 
                    tdiff[(thRef, vov)].SetPoint(tdiff[(thRef, vov)].GetN(), gt2.GetPointX(i),  gt2.GetPointY(i) - gt1.GetPointY(i+2) )
                    tdiff[(thRef, vov)].SetPointError(tdiff[(thRef, vov)].GetN()-1, gt2.GetErrorX(i),  math.sqrt(pow(gt2.GetErrorY(i),2) +pow(gt1.GetErrorY(i+2),2) ))
                else: print("Skipping point for t2(th1+%d)-t1 graph %.2f entry %d"%(delta, vov, i))

        else:
            for i in range(gt1.GetN()):
                if gt1.GetPointX(i) == gt2.GetPointX(i): 
                    tdiff[(thRef, vov)].SetPoint(tdiff[(thRef, vov)].GetN(), gt2.GetPointX(i),  gt2.GetPointY(i) - gt1.GetPointY(i) )
                    tdiff[(thRef, vov)].SetPointError(tdiff[(thRef, vov)].GetN()-1, gt2.GetErrorX(i),  math.sqrt(pow(gt2.GetErrorY(i),2) +pow(gt1.GetErrorY(i),2) ))
                else: print("Skipping point for t2(th1+%d)-t1 graph %.2f entry %d"%(delta, vov, i))


        g_tdiff[thRef].Add(tdiff[(thRef, vov)], "p")

        outfile.cd()
        leg.AddEntry(tdiff[(thRef, vov)], "V_{ov} = %.2f V"%vov, "pe")
        tdiff[(thRef, vov)].Write('g_t2-t1_vs_gainNpe_Vov%.2f_th%02d_delta%d'%(vov,thRef,delta))


    #print ("and now t2-t1")
    #g_tdiff[thRef].Print()
 
  
    # drawing and fitting tF

    c.Clear()
    hPad = ROOT.gPad.DrawFrame(0,0,12E09, 1.5)
    c.SetGridx()
    c.SetGridy()



    #dummy = ROOT.TF1 ("dummy_th%d"%(thRef+delta), "[0]*x^[1] ", 0, 12E09)
    fdiff[thRef]= ROOT.TF1 ("fdiff_th%d"%(thRef+delta), "[0]*x^[1] ", 0, 12E09)
    #fdiff[thRef]= ROOT.TF1 ("fdiff_th%d"%(thRef+delta), "[0]*x^([1]+[2]*x^[3]) ", 0, 12E09)
    #fdiff[thRef]= ROOT.TF1 ("fdiff_th%d"%(thRef+delta), "TMath.sqrt(TMath.pow([0]*x^[2], 2 ) + TMath.pow([0]*x^[2], 2 ))" , 0, 12E09)
    #fdiff[thRef]= ROOT.TF1 ("fdiff_th%d"%(thRef+delta), "(([0]*x^[1])^2  + ([2]*x^[3])^2)^0.5" , 0, 12E09)
    #fdiff[thRef]= ROOT.TF1 ("fdiff_th%d"%(thRef+delta), "cheb4", 0, 12E09)
    #dummy.SetParameters(1, ft2.GetParameter(2))

    fdiff[thRef].SetParameters(10E6, ft2.GetParameter(2))
    #g_tdiff[thRef].Fit(dummy, "N", "", 2E09, 9E09)
    g_tdiff[thRef].SetTitle("; gain x N_{pe} ; t2-t1 [ns] ")
     
    #fdiff[thRef].SetParameters(dummy.GetParameter(0),  dummy.GetParameter(1),  dummy.GetParameter(0),  dummy.GetParameter(1)/2.)
    #fdiff[thRef].SetParameters(dummy.GetParameter(0),  dummy.GetParameter(1),  1., dummy.GetParameter(1) )
    g_tdiff[thRef].Fit(fdiff[thRef], "N", "", 0.E09, 11E09)
    g_tdiff[thRef].Fit(fdiff[thRef], "N", "", 0.6E09, 10E09)
    #fdiff[thRef].SetLineColor(colors[thRef+delta])

    c.Clear()
    hPad = ROOT.gPad.DrawFrame(0,0,12E09, 15)
    c.SetGridx()
    c.SetGridy()


    fdiff[thRef].SetLineStyle(2)
    g_tdiff[thRef].Draw("psame")
    fdiff[thRef].Draw("same")

    leg.SetNColumns(2)
    leg.Draw("same")

    c.SaveAs("%s/t2-t1_th2_%02d.png"%(outdir,thRef+delta))
    outfile.cd()

    g_tdiff[thRef].Write('g_t2-t1glob_vs_gainNpe_th%02d_delta%d'%(thRef,delta))
    fdiff[thRef].Write('f_t2-t1glob_vs_gainNpe_th%02d_delta%d'%(thRef,delta))

