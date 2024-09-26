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

from collections import OrderedDict

import ROOT
import CMS_lumi, tdrstyle                                                                                                                                               
                                                                                                                                                                        
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

def getDiff(h1_deltaT, saveFit, outdir, gname):
   c = ROOT.TCanvas()
   c.Clear()
   c.cd()
   tDiff = [-1,-1]

   h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())
                    
   fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
   fitFunc.SetLineColor(ROOT.kGreen+3)
   #fitFunc.SetLineWidth(2)
   fitFunc.SetLineStyle(2)
   fitFunc.SetParameters(h1_deltaT.GetMaximum(),h1_deltaT.GetMean(), h1_deltaT.GetRMS())
   
   fitXMin = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 200
   fitXMax = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 200.
   fitFunc.SetRange(fitXMin, fitXMax)
   h1_deltaT.Fit('fitFunc','QNRL','', fitXMin, fitXMax)
   h1_deltaT.Fit('fitFunc','QNRL')
   fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2))
   h1_deltaT.Fit('fitFunc','QRSL+')

   tDiff = [ fitFunc.GetParameter(1),fitFunc.GetParError(1)]
   if saveFit:
       h1_deltaT.Draw("same histo")
       c.Update()
       c.SaveAs("%s/fits/%s.png"%(outdir, gname))
   return tDiff

def removeOffset1D(h, offset, hnew):
    for ibin in range(1, h.GetNbinsX()):
        cont = h.GetBinContent(ibin)
        err  = h.GetBinError(ibin)
        if cont == 0 : continue
        hnew.SetBinContent(ibin, cont-offset)
        hnew.SetBinError(ibin, err)

# --- colors
#cols = { 
#         0 : ROOT.kMagenta+1, 
#         1 : ROOT.kOrange+8,  
#         2 : ROOT.kOrange+4,  
#         3 : ROOT.kOrange+2,  
#         4 : ROOT.kRed+1,  
#         5 : ROOT.kPink+10, 
#         6 : ROOT.kOrange,  
#         7 : ROOT.kMagenta+3, 
#         8 : ROOT.kViolet+10, 
#         9 : ROOT.kBlue+2, 
#         10 : ROOT.kAzure+10, 
#         11 :ROOT.kCyan, 
#         12 :ROOT.kTeal+10, 
#         13 :ROOT.kGreen-3, 
#         14 :ROOT.kGreen+2 , 
#         15 :ROOT.kGreen+4, 
#}

cols = { 
         0 : ROOT.kBlue+2, 
         1 : ROOT.kAzure+7,  
         2 : ROOT.kCyan-3,  
         3 : ROOT.kTeal,  
         4 : ROOT.kGreen,  
         5 : ROOT.kOrange+1, 
         6 : ROOT.kOrange+7,  
         7 : ROOT.kRed+2, 
         8 : ROOT.kRed, 
         9 : ROOT.kOrange-3, 
         10 : ROOT.kOrange, 
         11 :ROOT.kGreen +1, 
         12 :ROOT.kTeal+10, 
         13 :ROOT.kCyan, 
         14 :ROOT.kAzure+10 , 
         15 :ROOT.kBlue, 
}




parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("-l",  "--label",   required=True, type=str, help="label in the form: HPK_2E14_C25_LYSO815_Vov1.50_T-30C, HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C")
parser.add_argument("-i",  "--inputFolder",  required=True, type=str, help="input folder")
parser.add_argument("-o",  "--outFolder",   required=True, type=str, help="out folder")
#parser.add_argument( "--saveFit",   action='store_true', help="Saving fit of the tDiff")
parser.add_argument( "--debug",   action='store_true', help="Debugging mode")
args = parser.parse_args()


irradiation = ''
bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
if args.label == 'HPK_2E14_C25_LYSO815_Vov1.50_T-30C':
  goodbars = [1, 2, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 1.50
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  angle = 52 # MAY Module were at 52°
  offsetX = 0 

elif args.label == 'HPK_2E14_C25_LYSO815_Vov1.50_T-35C':
  goodbars = [1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 1.50
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 12
  cellsize = '25 #mum'
  angle = 52 # MAY Module were at 52°
  offsetX = 0 
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'

elif args.label == 'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C':
  goodBars = {
           #5: [0,1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14, 15], not particularly good
           #7: [0,1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14, 15], only few bars have good MIP peak
           11: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           15: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           20: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           25: [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15]
           }
  vov = 1.50
  refbarmin = 4
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  angle = 49 # Sept Module were at 49°
  offsetX = 0


elif args.label == 'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C':
  goodBars = {
           15: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           20: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           25: [0,1, 2, 3, 4, 5,  7, 8, 9,  11, 12, 13, 15],
           }
  vov = 3.50
  refbarmin = 5
  refbarmax = 11
  cellsize = '25 #mum'
  irradiation = 'non irradiated'
  angle = 49 # Sept Module were at 49°
  offsetX = 0 

elif args.label == 'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C':
  goodBars = {
           #5: [0,1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           #7: [0,1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           11: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           #15: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           #20: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           #25: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
           }
  vov = 1.00
  refbarmin = 4
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = 'non irradiated'
  angle = 49 # Sept Module were at 49°
  offsetX = 0 


barConversionFact = 0.3122 / math.cos(angle*math.pi/180) # [cm]


label = args.label
outdir = '%s/%s/corrections/'%(args.outFolder, label) 
#'/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/DPGstudies_4DigitizationModel/tDiff/Sept24/'%label
inputdir = args.inputFolder #'/afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/'

#--- prepare outdir
if (os.path.isdir(outdir) == False): 
    os.system('mkdir %s'%outdir)

print (outdir)
g = {}

#outfile   = ROOT.TFile.Open(outdir+'/uniformityCheck_%s.root'%label,'recreate')
graphname = 'deltaT_vs_totRatio'
print ('Doing ', graphname)
#c = {}
#gg = OrderedDict()
th = 11
bar = 15
totCorr = OrderedDict()
c = ROOT.TCanvas('c_totCorr_bar%02d_th%02d'%(bar, th),'c_totCorr_bar%02d_th%02d'%(bar, th), 500, 500)
hPad = ROOT.gPad.DrawFrame(0.6, -100, 1.3, 100)
hPad.SetTitle('; ToT ratio; #Deltat [ps]')
ROOT.gPad.SetTicks(1)
c.cd()


for refbar in range(refbarmin, refbarmax+1):
    filename = '%s/moduleCharacterization_step2_%s_refbar%d.root'%(inputdir, label, refbar)
    f1 = ROOT.TFile.Open(filename)
   
    # remove offset in y
    p1TotCorr_0 = f1.Get('p1_deltaT_vs_totRatio_bar%02dL-R_Vov%.2f_th%02d_energyBin01'%(bar,vov,th))
    
    nbinsx = p1TotCorr_0.GetNbinsX()
    xmin = p1TotCorr_0.GetBinLowEdge(1)+p1TotCorr_0.GetBinWidth(1)
    xmax = p1TotCorr_0.GetBinLowEdge(nbinsx)+p1TotCorr_0.GetBinWidth(nbinsx)
    c = ROOT.TCanvas('c_totCorr_bar%02d_th%02d_refbar%02d'%(bar, th, refbar),'c_totCorr_bar%02d_th%02d_refbar%02d'%(bar, th, refbar), 500, 500)
    hPad = ROOT.gPad.DrawFrame(xmin, -100, xmax, 100)
    hPad.SetTitle('; ToT ratio; #Deltat [ps]')
    ROOT.gPad.SetTicks(1)
    c.cd()

    p1TotCorr = ROOT.TH1F('p1TotCorr_refbar%02d'%refbar, '', nbinsx, xmin, xmax)
    removeOffset1D(p1TotCorr_0, p1TotCorr_0.GetMean(2), p1TotCorr) 
    totCorr[refbar] = p1TotCorr
    p1TotCorr.Draw("pl same")

    c.SaveAs(outdir+'%s.png'%c.GetName())
#print (totCorr)
#for refbar, p1totCorr in totCorr.items():
#    # ToT ratio corr
#    p1totCorr.Draw('same')
#    fun = ROOT.TF1('fun','pol3')
#    fun.SetLineColor(2)
#    #p1TotCorr.Fit(fun,'QRS', '', p1TotCorr_0.GetMean()-3*p1TotCorr_0.GetRMS(), p1TotCorr_0.GetMean()+3*p1TotCorr_0.GetRMS())
#    #tl.DrawLatex(0.65,0.82,'#splitline{V_{OV} = %.2f V}{HPK, 25 #mum}'%(vov))
#    #cms_logo = draw_logo()
#    #cms_logo.Draw()

#for refth, goodbars in goodBars.items():
#
#    print ("Now ref th is %02d, and these are the good DUT bars: "%refth, goodbars)
#    leg = ROOT.TLegend(0.65,0.7,0.92,0.92)
#    leg.SetBorderSize(0)
#    leg.SetFillStyle(0)
#    leg.SetNColumns(2)
#
#    #for refbar in range(4,6): #loop on coincidence bars--> REF bars
#    for bar in range(0,2): #loop on coincidence bars--> REF bars
#        #filename = '%s/moduleCharacterization_step2_%s_refbar%d.root'%(inputdir, label, refbar)
#        #f = ROOT.TFile.Open(filename)
#        #if f == None: continue
#        #if args.debug: print (filename)
#        #for bar in  goodbars:
#        for refbar in range(4,6): #loop on coincidence bars--> REF bars
#            filename = '%s/moduleCharacterization_step2_%s_refbar%d.root'%(inputdir, label, refbar)
#            f = ROOT.TFile.Open(filename)
#            if f == None: continue
#            if args.debug: print (filename)
#            gnameComplete = 'p1_%s_bar%02dL-R_Vov%.2f_th%02d_energyBin01'%( graphname, bar, vov, refth)
#            #totRatio = f.Get('h1_totRatio_bar%02dL-R_Vov%.2f_th%02d_energyBin01'%(bar, vov, refth))
#            gdummy = f.Get('h1_totRatio_bar%02dL-R_Vov%.2f_th%02d_energyBin01'%(bar, vov, refth))
#            #if args.debug: print (gnameComplete)
#            #print (bar, refbar)
#            #gdummy = f.Get(gnameComplete)
#            #gdummy.SetName("%s_REF%02d"%(gnameComplete, refbar))
#            gg ['refbar%02d_bar%02d'%(refbar,bar)] = ROOT.TProfile(f.Get(gnameComplete))
#            #if g['refbar%02d_bar%02d'%(refbar,bar)] == None: 
#            #    print ("%s not found"%gnameComplete)
#            #    continue
#
#        print (gg) 
#    for bar in goodbars:
#        c =  ROOT.TCanvas('c_%s_vs_x_bar%02d'%(graphname, bar),'c_%s_vs_x_bar%02d'%(graphname, bar),600,500)
#        c.cd()
#        c.Clear()
#        hdummy.Draw()
#        #print (gs) 
#        for refbar in range(4,6):
#            gr = gg ['refbar%02d_bar%02d'%(refbar,bar)] 
#    
#            if gr == None: 
#                print("not found")
#                continue
#            gr.Print() 
#            gr.SetMarkerColor(cols[refbar])
#            gr.SetMarkerSize(1)
#            gr.SetLineStyle(1)
#            gr.SetLineWidth(1)
#            gr.Draw('p same')
#        c.SaveAs('%s/%s_refTh%02d.png'%(outdir,c.GetName(), refth))
#
