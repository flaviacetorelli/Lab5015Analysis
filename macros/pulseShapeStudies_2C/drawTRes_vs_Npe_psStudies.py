#! /usr/bin/env python3
# Macro to calculate SR from TB data (similar to what is done for laser data)
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

colors = {
    'L' : ROOT.kRed + 1,
    'R' : ROOT.kBlue + 2,
    'L-R' : ROOT.kBlack,
        }

markers = {
    0.60 : 20, 
    0.80 : 21, 
    1.00 : 22, 
    1.25 : 23, 
    1.50 : 24, 
    2.00 : 25, 
    3.50 : 26 
        }

dac_to_uA = {
        1 : 0.313
        }
ithMode = 1
label = 'HPK_nonIrr_C25_LYSO818'
sipmType = 'HPK-PIT-C25-ES2'
LO_at3p5 = 2390
temp = 'T5C'
inputdir = '/afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/plots/TOFHIR2C/' 
outdir = '/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/pulseShapes/vsNPE/%s_%s/'%(label, temp)
thFixed = 11
bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
Vovs = [0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 3.50]
goodbars = {
        0.60: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
        0.80: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
        1.00: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
        1.25: [0,1,2,4,5,6,7,8,9,10,11,12,13,14,15], # check 2 as well
        1.50: [0,1,2,4,5,9,10,11,12,13,14,15],
        2.00: [14,15],
        3.50: []
    }
g_SR_vs_gainNpe = {}
g_t1_vs_gainNpe = {}
g_SR_vs_bar = {}
g_SRAve_vs_gainNpe = {}

angleFact = math.cos(49 * math.pi / 180) / math.cos(52 * math.pi / 180)
#print (angleFact)
for Vov in Vovs:

        g_SRAve_vs_gainNpe[Vov] = ROOT.TGraphErrors()
        g_SRAve_vs_gainNpe[Vov].SetName('g_pulseShapeL-R_barsAve_Vov%.2f'%Vov)
        g_SR_vs_bar[Vov] = ROOT.TGraphErrors()

        Npe = (LO_at3p5* PDE(sipmType, Vov) /PDE(sipmType, 3.50))* 4.2 * 1.25  * (1./ angleFact)
        gain = Gain(sipmType,Vov)
        print ('Gain = %f and Npe = %f'%(gain,Npe))

        for bar in bars:
                # get SR for optimal or fixed th threshold
                #-----------------------------
                thRef = thFixed #FIXME choose the type of THR (fixed or best)
                inFile = None
                inFile = ROOT.TFile.Open(inputdir+'pulseShape_%s_Vov%.2f_%s.root'%(label, Vov, temp))
                print (inFile)
                print ('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,Vov))
                SR = 0.
                t1 = 0.
                channels = ['L', 'R']
                nPoints = 5

                graph1 = inFile.Get('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,Vov))
                graph2 = inFile.Get('g_pulseShapeR_bar%02d_Vov%.2f'%(bar,Vov))
                if graph1 == None or graph2 == None:
                        continue
                for ch in ['L', 'R', 'L-R']:
                        g_SR_vs_gainNpe[(Vov, bar, ch)] = ROOT.TGraphErrors()
                        g_SR_vs_gainNpe[(Vov, bar, ch)].SetMarkerStyle(markers[Vov])
                        g_SR_vs_gainNpe[(Vov, bar, ch)].SetMarkerColor(colors[ch])
                        g_SR_vs_gainNpe[(Vov, bar, ch)].SetName('g_pulseShape%s_bar%02d_Vov%.2f'%(ch, bar,Vov))

                for ch in channels:
                        graph = inFile.Get('g_pulseShape%s_bar%02d_Vov%.2f'%(ch, bar,Vov))

                        index_cen = 0
                        for point in range(0,int(graph.GetN()/2)):
                                if graph.GetPointY(point) >= 0.999 * ( thRef * dac_to_uA[ithMode] ):
                                        index_cen = point
                                        break
                        index_min = max(0,index_cen-int(nPoints/2))
                        index_max = min(index_cen+int(nPoints/2),int(min(graph1.GetN()/2+1,graph2.GetN()/2+1)))

                        fitSR = ROOT.TF1('fitSR', 'pol1', -10., 30.)
                        fitSR.SetLineWidth(1)
                        fitSR.SetRange( graph.GetPointX(index_min)-0.001, graph.GetPointX(index_max)+0.001 )
                        fitSR.SetParameters(0,((graph.GetPointY(index_max)+0.001) - (graph.GetPointY(index_min)-0.001)) /((graph.GetPointX(index_max)+0.001) -  (graph.GetPointX(index_min)-0.001)) )
                        graph.Fit(fitSR,'QRS')

                        ctemp = ROOT.TCanvas('ctemp_Vov%.2f_bar%02d%s'%(Vov,bar, ch),'ctemp_Vov%.2f_bar%02d%s'%(Vov,bar, ch))
                        #graph.GetXaxis().SetRangeUser(graph.GetX()[1]-0.5,graph.GetX()[1]+2)
                        graph.Draw('ap')
                        line1 = ROOT.TLine(graph.GetX()[1]-0.5, thRef * dac_to_uA[ithMode], graph.GetX()[1]+2,  thRef * dac_to_uA[ithMode]);
                        line1.SetLineWidth(1)
                        line1.SetLineStyle(7)
                        line1.Draw("same")
                        ctemp.Print('%s/%s.png'%(outdir, ctemp.GetName()))

                        # storing single channels and Average
                        g_SR_vs_gainNpe[(Vov, bar, ch)].SetPoint(g_SR_vs_gainNpe[(Vov,bar, ch)].GetN(),gain*Npe,fitSR.GetParameter(1))
                        g_SR_vs_gainNpe[(Vov, bar, ch)].SetPointError(g_SR_vs_gainNpe[(Vov,bar, ch)].GetN()-1,0.1*gain*Npe,0.10*fitSR.GetParameter(1))

                        #g_t1_vs_gainNpe[(Vov, bar, ch)].SetPoint(g_t1_vs_gainNpe[(Vov,bar, ch)].GetN(),Gain*Npe,t1)
                        #g_t1_vs_gainNpe[(Vov, bar, ch)].SetPointError(g_t1_vs_gainNpe[(Vov,bar, ch)].GetN()-1,0.03*Gain*Npe,0)

                        print('index_cen = %d   index_min = %d   index_max = %d  t1 = %.1f  amp = %.1f   SR = %.1f'%(index_cen,index_min,index_max,graph.GetX()[index_cen],graph.GetY()[index_cen],fitSR.GetParameter(1)))

                        SR += fitSR.GetParameter(1)
                        t1 += graph.GetX()[index_cen]
                SR /= 2.
                t1 /= 2.
                #print(g_SR_vs_gainNpe[(data_struct.label,Vov, 'ave')].GetN(),Gain(data_struct.sipmType,data_struct.Vov)*Npe,SR)
                if bar in goodbars[Vov]: 
                     g_SR_vs_bar[Vov].SetPoint(g_SR_vs_bar[Vov].GetN(), bar, SR)
                     g_SR_vs_bar[Vov].SetPointError(g_SR_vs_bar[Vov].GetN()-1, 0., SR*0.1)
                g_SR_vs_gainNpe[(Vov, bar, 'L-R')].SetPoint(g_SR_vs_gainNpe[((Vov, bar, 'L-R'))].GetN(),gain*Npe,SR)
                g_SR_vs_gainNpe[(Vov, bar, 'L-R')].SetPointError(g_SR_vs_gainNpe[((Vov, bar, 'L-R'))].GetN()-1,gain*Npe*0.1,SR*0.10)

        # SR vs bar to obtain average over bars in a module
        cvsBar = ROOT.TCanvas('c_SR_vs_bar_Vov%.2f'%Vov, 'c_SR_vs_bar_Vov%.2f'%Vov)
        hPadBar = ROOT.gPad.DrawFrame(-0.5,0.,15.5, 40.)
        hPadBar.SetTitle("; bar; Slew Rate [#mu A / ns]")
        cvsBar.SetGridy()
        pol0 = ROOT.TF1("pol0", "pol0", -0.5, 15.5)
        g_SR_vs_bar[Vov].Fit(pol0, "QRS")
        g_SR_vs_bar[Vov].Draw("pl same")
        
        cvsBar.SaveAs(outdir+cvsBar.GetName()+'.png')
        g_SRAve_vs_gainNpe[Vov].SetPoint(g_SRAve_vs_gainNpe[Vov].GetN(),gain*Npe,pol0.GetParameter(0))
        g_SRAve_vs_gainNpe[Vov].SetPointError(g_SRAve_vs_gainNpe[Vov].GetN()-1,gain*Npe*0.1,pol0.GetParError(0))



######### Now drawing for each bar several OVs
for bar in bars:
        c =  ROOT.TCanvas('c_SR_vs_gainNpe_allVovs_bar%02d'%(bar),'c_SR_vs_gainNpe_allVovs_bar%02d'%(bar),1200,700)
        hPad = ROOT.gPad.DrawFrame(0,0.,14000E06, 40.)
        hPad.SetTitle("; Gain x N_{pe}; Slew Rate [#mu A / ns]")
        c.SetGrid()
        for Vov in Vovs:
                g_SR_vs_gainNpe[(Vov, bar, 'L')].Draw("pl same")
                g_SR_vs_gainNpe[(Vov, bar, 'R')].Draw("same pl")
                g_SR_vs_gainNpe[(Vov, bar, 'L-R')].Draw("same pl")
        
        c.SaveAs(outdir+c.GetName()+'.png')
        c.SaveAs(outdir+c.GetName()+'.pdf')
        c.SaveAs(outdir+c.GetName()+'.root')
        

c =  ROOT.TCanvas('c_SR_vs_gainNpe_allVovs_barsAve','c_SR_vs_gainNpe_allVovs_barsAve',1200,700)
hPad = ROOT.gPad.DrawFrame(0,0.,14000E06, 40.)
hPad.SetTitle("; Gain x N_{pe}; Slew Rate [#mu A / ns]")
c.SetGrid()
for Vov in Vovs:
        g_SRAve_vs_gainNpe[Vov].Draw("pl same")

c.SaveAs(outdir+c.GetName()+'.png')
c.SaveAs(outdir+c.GetName()+'.pdf')
c.SaveAs(outdir+c.GetName()+'.root')





