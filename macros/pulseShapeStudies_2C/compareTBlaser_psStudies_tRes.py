#! /usr/bin/env python3
# This is a macro to compare laser results and TB results with the noise + stoch term from parametrization of tres
# you need files with SR vs gainNpe stored and the ones with the tRes (from laser)
# you need also the file with TB (nb here on tb point is applied the angle scaling factor)

import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import json

from collections import OrderedDict

import ROOT
import CMS_lumi, tdrstyle                                                                                                                                               

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
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning   

def noise_vs_Npe(x, par):
    xx = x[0]*Gain(sipmType,vov)
    #print (data_struct.Vov, Gain(data_struct.sipmType,data_struct.Vov) ,  x[0])
    noise_single = math.sqrt( pow(par[3]/f_SR_vs_gainNpe_ref.Eval(xx),2) + 16.7*16.7 )
    return noise_single / math.sqrt(2)

def stoch_vs_Npe(x, par):
     xx = x[0]
     return par[1] * pow(par[0]/xx,par[2])

def tot_vs_Npe(x, par):
    xx = x[0]
    return math.sqrt( pow(noise_vs_Npe(x,par),2) + pow(stoch_vs_Npe(x,par),2) )

def noise_vs_Npe_log(x, par):
    xx = x[0]*Gain(sipmType,vov)
    #print (data_struct.Vov, Gain(data_struct.sipmType,data_struct.Vov) ,  x[0])
    noise_single = math.sqrt( pow(par[3]/f_SR_vs_gainNpe_log_ref.Eval(xx),2) + 16.7*16.7 )
    return noise_single / math.sqrt(2)
def tot_vs_Npe_log(x, par):
    xx = x[0]
    return math.sqrt( pow(noise_vs_Npe_log(x,par),2) + pow(stoch_vs_Npe(x,par),2) )





colors = {
    'L' : ROOT.kRed + 1,
    'R' : ROOT.kBlue + 2,
    'L-R' : ROOT.kBlack,
        }

markers = {
    0.60 : 20, 
    0.80 : 25, 
    1.00 : 22, 
    1.25 : 23, 
    1.50 : 24, 
    2.00 : 21, 
    3.00 : 27, 
    3.50 : 26 
        }
colors = {
    7 : ROOT.kPink +1 , 
    11: ROOT.kOrange, 
    15: ROOT.kRed,  
    20: ROOT.kGreen -4, 
    25: ROOT.kCyan, 
        }
dac_to_uA = {
        1 : 0.313
        }
ithMode = 1
label = 'HPK_nonIrr_C25_LYSO818'
sipmType = 'HPK-PIT-C25-ES2'
LO_at3p5 = 2390
temp = 'T5C'
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/compareLaser/'%(label, temp)

thBestFromTB = {
        0.60: 5,
        0.80: 7,
        1.00: 11,
        1.25: 15,
        1.50: 15,
        2.00: 20,
        3.00: 20,
        3.50: 25,
        }
#thFixed = [7, 11, 15, 20, 25 ]
thFixed = [11, 15 ]
bars = [0,1,2,3,4,5,6,7,8,9,10,11,13,14,15]
vovs = [ 0.80, 1.00, 1.50, 2.00, 3.00]
       
angleFact = math.cos(49 * math.pi / 180) / math.cos(52*math.pi / 180)
inFileTB = ROOT.TFile.Open('/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/ModuleCharacterization/HPK_nonIrr_C25_LYSO818_T5C/summaryPlots_HPK_nonIrr_C25_LYSO818_T5C.root')
inFileLaser = ROOT.TFile.Open('plots_tRes_LYSO818_pulseShapeStudy_vsNpe.root')


#for i,vov in enumerate(vovs):
#    c = ROOT.TCanvas("c_tRes_vs_Npe_vov%0.2f"%vov, "c_tRes_vs_Npe_vov%0.2f"%vov)
#    c.cd()
#    c.SetGridy()
#    hPad = ROOT.gPad.DrawFrame(0,0.,11E03, 140.)
#    hPad.SetTitle(";  N_{pe}; time resolution [ps]")
# 
#    Npe = LO_at3p5* (PDE(sipmType, vov) /PDE(sipmType, 3.50))* (4.2) * 1.25
#    leg = ROOT.TLegend(0.8, 0.5, 0.95 , 0.95)
#
#    for th in thFixed:
#         g = inFileLaser.Get('tRes_vs_Npe_vov%0.2f_th%02d'%(vov, th))
#         print ('tRes_vs_Npe_vov%0.2f_th%02d'%(vov, th))
#         g.SetMarkerColor(colors [th])
#         g.SetLineColor(colors [th])
#         g.SetMarkerStyle(markers [vov])
#         g.Draw("same PL")       
#         leg.AddEntry(g, "th = %02d"%th, "PL" )
#    leg.Draw("same")
#    c.SaveAs(outdir+c.GetName()+'.png')
        
for th in thFixed:
    inFile = ROOT.TFile.Open('plots_SRt1_LYSO818_pulseShapeStudy_vsNpe.root')
    f_SR_vs_gainNpe_ref = inFile.Get('f_SRglob_vs_gainNpe_lin_th%02d'%th)
    f_SR_vs_gainNpe_log_ref = inFile.Get('f_SRglob_vs_gainNpe_log_th%02d'%th)

   
    tRes_vs_vov = inFileTB.Get('g_deltaT_totRatioCorr_vs_vov_bar07_th%02d_enBin01'%th)

    for i,vov in enumerate(vovs):

        tRes_TB = ROOT.TGraphErrors()
        Npe = LO_at3p5* (PDE(sipmType, vov) /PDE(sipmType, 3.50))* (4.2) * 1.25

        for point in range(tRes_vs_vov.GetN()):
            if vov == tRes_vs_vov.GetPointX(point): 
                tRes_TB.SetPoint(tRes_TB.GetN(), Npe, tRes_vs_vov.GetPointY(point)/angleFact) 
                tRes_TB.SetPointError(tRes_TB.GetN()-1, Npe*0.1, tRes_vs_vov.GetErrorY(point)/angleFact) 

        c = ROOT.TCanvas("c_tRes_vs_Npe_vov%0.2f_th%02d"%(vov, th), "c_tRes_vs_Npe_vov%0.2f_th%02d"%(vov,th))
        c.cd()
        c.SetGridy()
        hPad = ROOT.gPad.DrawFrame(0,0.,11E03, 120.)
        hPad.SetTitle(";  N_{pe}; time resolution [ps]")
        leg1 = ROOT.TLegend(0.5, 0.7, 0.95, 0.95)     

        tRes_TB.SetMarkerColor(ROOT.kBlack)       
        tRes_TB.Draw("same P")
        leg1.AddEntry(tRes_TB,"TB Sept23" ,"P")

        g = inFileLaser.Get('tRes_vs_Npe_Vov%0.2f_th%02d'%(vov, th))
        g.Print()
        g.SetMarkerColor(ROOT.kRed)
        g.SetLineColor(ROOT.kRed)
        g.SetMarkerStyle(21)
        g.Draw("same P")       
        leg1.AddEntry(g,"Laser points" ,"P")

        NpeMin = g.GetPointX(0)-100 
        NpeMax = g.GetPointX(g.GetN()-1)+100

        fitFunc_tot = ROOT.TF1("fitFunc_totEGainNPE_lin",tot_vs_Npe,NpeMin, NpeMax,4) #here, SR vs gainNPE
        fitFunc_tot.FixParameter(0,7000)
        fitFunc_tot.FixParameter(1,25.7)
        fitFunc_tot.FixParameter(2,0.5)
        fitFunc_tot.FixParameter(3,457.)
        fitFunc_tot.SetLineStyle(7)
        fitFunc_tot.SetLineWidth(2)
        fitFunc_tot.SetLineColor(ROOT.kPink+6)
        fitFunc_tot.Draw("same")

        leg1.AddEntry(fitFunc_tot,"exp. tRes, SR vs gainNpe linear" ,"L")
        fitFunc_tot1 = ROOT.TF1("fitFunc_totEGainNPE_log",tot_vs_Npe_log,NpeMin,NpeMax,4) #here, SR vs gainNPE
        fitFunc_tot1.FixParameter(0,7000)
        fitFunc_tot1.FixParameter(1,25.7)
        fitFunc_tot1.FixParameter(2,0.5)
        fitFunc_tot1.FixParameter(3,457.)
        fitFunc_tot1.SetLineStyle(7)
        fitFunc_tot1.SetLineWidth(2)
        fitFunc_tot1.SetLineColor(ROOT.kViolet+2)
        fitFunc_tot1.Draw("same")
        fitFunc_tot.Draw("same")
        leg1.AddEntry(fitFunc_tot1,"exp. tRes, SR vs gainNpe log" ,"L")
        leg1.Draw("same")


 
        c.SaveAs(outdir+c.GetName()+'.png')


        c = ROOT.TCanvas("c_tRes_vs_Npe_vov%0.2f_Lin_th%02d"%(vov, th), "c_tRes_vs_Npe_vov%0.2f_Lin_th%02d"%(vov,th))
        c.cd()
        c.SetGridy()
        hPad = ROOT.gPad.DrawFrame(0,0.,11E03, 120.)
        hPad.SetTitle(";  N_{pe}; time resolution [ps]")
        leg1 = ROOT.TLegend(0.5, 0.7, 0.95, 0.95)     

        tRes_TB.Draw("same P")
        leg1.AddEntry(tRes_TB,"TB Sept23" ,"P")

        g.Draw("same P")       
        leg1.AddEntry(g,"Laser points" ,"P")

        fitFunc_tot.Draw("same")
        leg1.AddEntry(fitFunc_tot,"exp. tRes, SR vs gainNpe linear" ,"L")

        noise = ROOT.TF1("noise",noise_vs_Npe,NpeMin,NpeMax,4) #here, SR vs gainNPE
        noise.FixParameter(3,457.)
        noise.SetLineStyle(7)
        noise.SetLineWidth(2)
        noise.SetLineColor(ROOT.kBlue)
        noise.Draw("same")
        leg1.AddEntry(noise,"noise" ,"L")



        stoch = ROOT.TF1("stoch",stoch_vs_Npe,NpeMin,NpeMax,4) #here, SR vs gainNPE
        stoch.FixParameter(0,7000)
        stoch.FixParameter(1,25.7)
        stoch.FixParameter(2,0.5)
        stoch.SetLineStyle(7)
        stoch.SetLineWidth(2)
        stoch.SetLineColor(ROOT.kGreen + 1)
        stoch.Draw("same")
        leg1.AddEntry(stoch,"stochastic" ,"L")

        leg1.Draw("same")

        c.SaveAs(outdir+c.GetName()+'.png')
 

        c = ROOT.TCanvas("c_tRes_vs_Npe_vov%0.2f_Log_th%02d"%(vov, th), "c_tRes_vs_Npe_vov%0.2f_Lin_th%02d"%(vov,th))
        c.cd()
        c.SetGridy()
        hPad = ROOT.gPad.DrawFrame(0,0.,11E03, 120.)
        hPad.SetTitle(";  N_{pe}; time resolution [ps]")
        leg1 = ROOT.TLegend(0.5, 0.7, 0.95, 0.95)     

        tRes_TB.Draw("same P")
        leg1.AddEntry(tRes_TB,"TB Sept23" ,"P")

        g.Draw("same P")       
        leg1.AddEntry(g,"Laser points" ,"P")

        fitFunc_tot1.Draw("same")
        leg1.AddEntry(fitFunc_tot1,"exp. tRes, SR vs gainNpe log" ,"L")

        noise_log = ROOT.TF1("noise_log", noise_vs_Npe_log,NpeMin,NpeMax,4) #here, SR vs gainNPE
        noise_log.FixParameter(3,457.)
        noise_log.SetLineStyle(7)
        noise_log.SetLineWidth(2)
        noise_log.SetLineColor(ROOT.kBlue)
        noise_log.Draw("same")
        leg1.AddEntry(noise,"noise, SR vs gainNpe log" ,"L")

        stoch.Draw("same")
        leg1.AddEntry(stoch,"stochastic" ,"L")

        leg1.Draw("same")


        c.SaveAs(outdir+c.GetName()+'.png')
 
        c = ROOT.TCanvas("c_tRes_vs_Npe_vov%0.2f_th%02d_v2"%(vov, th), "c_tRes_vs_Npe_vov%0.2f_Lin_th%02d_v2"%(vov,th))
        c.cd()
        c.SetGridy()
        hPad = ROOT.gPad.DrawFrame(0,0.,11E03, 120.)
        hPad.SetTitle(";  N_{pe}; time resolution [ps]")
        leg1 = ROOT.TLegend(0.5, 0.65, 0.95, 0.95)     

        tRes_TB.Draw("same P")
        leg1.AddEntry(tRes_TB,"TB Sept23" ,"P")

        g.Draw("same P")       
        leg1.AddEntry(g,"Laser points" ,"P")

        fitFunc_tot.Draw("same")
        leg1.AddEntry(fitFunc_tot,"exp. tRes, SR vs gainNpe lin" ,"L")



        fitFunc_tot1.Draw("same")
        leg1.AddEntry(fitFunc_tot,"exp. tRes, SR vs gainNpe log" ,"L")

        noise_log.SetLineColor(ROOT.kBlue + 3)
        noise_log.Draw("same")
        leg1.AddEntry(noise,"noise, SR vs gainNpe log" ,"L")

        noise.Draw("same")
        leg1.AddEntry(noise,"noise, SR vs gainNpe log" ,"L")

        stoch.Draw("same")
        leg1.AddEntry(stoch,"stochastic" ,"L")

        leg1.Draw("same")

        c.SaveAs(outdir+c.GetName()+'.png')
 
 
