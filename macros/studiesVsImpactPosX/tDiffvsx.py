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
parser.add_argument("-g",  "--gname",   required=True, type=str, help="Choose corrected deltaT: deltaT_totRatioPhaseCorr or deltaT_energyRatioPhaseCorr")
parser.add_argument("-l",  "--label",   required=True, type=str, help="label in the form: HPK_2E14_C25_LYSO815_Vov1.50_T-30C, HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C")
parser.add_argument("-i",  "--inputFolder",  required=True, type=str, help="input folder")
parser.add_argument("-o",  "--outFolder",   required=True, type=str, help="out folder")
parser.add_argument( "--saveFit",   action='store_true', help="Saving fit of the tDiff")
parser.add_argument( "--debug",   action='store_true', help="Debugging mode")
args = parser.parse_args()


irradiation = ''
bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
if args.label == 'HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C':
  goodBars = {
           5: [0,1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15], # 7
           7: [0,1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15], # 7
           11: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           15: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           20: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           25: [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15] #1 not there, 14
           }
  vov = 1.00
  refbarmin = 2
  refbarmax = 11
  cellsize = '25 #mum'
  irradiation = 'non irradiated'
  angle = 52 # May module at 52
  offsetX = 0


elif args.label == 'HPK_2E14_C25_LYSO815_Vov1.50_T-35C':
  goodBars = { #bar3 no visible and 15L as well
           5: [0,1, 2,  4, 5,  6, 10, 11, 12, 13, 14], #7,8,9 R not fitted
           7: [0,1, 2,  4, 5,  6, 10, 11, 12, 13, 14], #7,8,9 R not fitted
           11: [0,1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14], #8L-R
           15: [0,1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
           20: [0,1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
           25: [0, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], # 1,14 not there
           }
  vov = 1.50
  refbarmin = 2
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  angle = 52 # MAY Module were at 52°
  #offsetX = 6 #central bar
  offsetX = 0 #central bar



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
           5: [0,1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           7: [0,1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           11: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           15: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           20: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           25: [0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
           }
  vov = 1.00
  refbarmin = 4
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = 'non irradiated'
  angle = 49 # Sept Module were at 49°
  offsetX = 0 


barConversionFact = 0.312 / math.cos(angle*math.pi/180) # [cm]


label = args.label
outdir = '%s/%s'%(args.outFolder, label) 
#'/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/DPGstudies_4DigitizationModel/tDiff/Sept24/'%label
inputdir = args.inputFolder #'/afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/'

#--- prepare outdir
if (os.path.isdir(outdir) == False): 
    os.system('mkdir %s'%outdir)

print (outdir)

# Canvas things
hdummy = ROOT.TH2F('hdummy','',100,-30,30,1000,-1200,1500)
hdummy.GetXaxis().SetTitle('x [cm]')
hdummy.GetYaxis().SetTitle('#DeltaT [ps]')
hdummy.GetXaxis().SetRangeUser(0, 6 )
hdummy.GetYaxis().SetRangeUser(-1100, 1400 )

latex = ROOT.TLatex(0.65,0.84,'%s'%(irradiation))
latex.SetNDC()
latex.SetTextSize(0.038)
latex.SetTextFont(42)

c = ROOT.TCanvas("","",600,500)
#outfile   = ROOT.TFile.Open(outdir+'/uniformityCheck_%s.root'%label,'recreate')
gVsTh = ROOT.TGraphErrors()
graphname = args.gname
print ('Doing ', graphname)

for refth, goodbars in goodBars.items():
    print ("Now ref th is %02d, and these are the good DUT bars: "%refth, goodbars)

    g_tDiff_vs_x = OrderedDict()
    for bar in goodbars:
        g_tDiff_vs_x [bar] = ROOT.TGraphErrors()
    #print (g_tDiff_vs_x)

    for refbar in range(refbarmin,refbarmax+1): #loop on coincidence bars--> REF bars
        filename = '%s/moduleCharacterization_step2_%s_refbar%d.root'%(inputdir, label, refbar)
        f = ROOT.TFile.Open(filename)
        if f == None: continue
        if args.debug: print (filename)
    
        for bar, gDiff in  g_tDiff_vs_x.items():
        #for bar in goodbars:
            #print ("DOING BAR %02d"%bar)
            #g_tDiff_vs_x[bar].Print()
            gnameComplete = 'h1_%s_bar%02dL-R_Vov%.2f_th%02d_energyBin01'%( graphname, bar, vov, refth)
            if args.debug: print (gnameComplete)
            h1_deltaT = f.Get(gnameComplete)
            if h1_deltaT == None: 
                print ("%s not found"%gnameComplete)
                continue
            tDiff = getDiff( h1_deltaT, args.saveFit, outdir, gnameComplete)
            gDiff.SetPoint(gDiff.GetN(), (refbar-offsetX)*barConversionFact, tDiff[0])
            gDiff.SetPointError(gDiff.GetN()-1, 0, tDiff[1])
            
    c1 =  ROOT.TCanvas('c_%s_vs_x_bar'%graphname,'c_%s_vs_x_bar'%graphname,600,500)
    c2 =  ROOT.TCanvas('c_%s_vs_x_all'%graphname,'c_%s_vs_x_all'%graphname,600,500)
    c2.cd()
    hdummy.Draw()
    
    leg = ROOT.TLegend(0.65,0.7,0.92,0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    
    # define 1 graph for each bar of the DUT
    if '813' or '815' in args.label: 
        hSummary = ROOT.TH1F("h_refth%02d"%refth, "h_refth%02d"%refth, 50, -400, 0 )
    else:
        hSummary = ROOT.TH1F("h_refth%02d"%refth, "h_refth%02d"%refth, 100, -400, 400 )
    hSummary.GetXaxis().SetRangeUser(-250, -100)
   
    # draw Plots and do fit 
    for bar,g in g_tDiff_vs_x.items():
    #for bar in [0,8,15]:
        #g =  g_tDiff_vs_x [bar]
        #print ("DOING BAR %02d"%bar)
        #g.Print()
        c1.Clear() 
        c1.SetGridy()
        c1.SetGridx()
        c1.cd()
    
        hdummy.Draw()
        lin = ROOT.TF1("pol1", "pol1", -10, 10)
        lin.SetRange(g.GetPointX(0)-0.1,g.GetPointX(g.GetN()-1)+0.1 )
        g.SetMarkerStyle(20)
        g.SetMarkerColor(ROOT.kBlue+1)
        g.SetMarkerSize(1)
        g.SetLineStyle(2)
        g.SetLineWidth(1)
        g.Draw('p same')
        
        g.Fit(lin, "QNR")
    
        lin.SetLineStyle(2)
        lin.SetLineColor(ROOT.kBlue+1)
        lin.Draw("same")
    
        c1.SaveAs('%s/%s%02d_ref%02d.png'%(outdir,c1.GetName(), bar, refth ))
        #c1.SaveAs(outdir+c1.GetName()+'.pdf')
    
        # now all toghether
        c2.SetGridy()
        c2.SetGridx()
        c2.cd()
    
        leg.AddEntry(g, 'bar %2d'%bar, 'PL')
        g.SetMarkerColor(cols[bar])
        g.SetMarkerSize(1)
        g.SetLineStyle(1)
        g.SetLineWidth(1)
        g.Draw('p same')
    
        lin.SetLineColor(cols[bar])
        lin.Draw("same")
    
        leg.Draw("same")
     
        hSummary.Fill(lin.GetParameter(1))
    
    c2.SaveAs('%s/%s_refTh%02d.png'%(outdir,c2.GetName(), refth))
    c =  ROOT.TCanvas('LCSlope_summary_%s_refTh%02d'%(graphname, refth),'LCSlope_summary_%s_refTh%02d'%(graphname, refth),600,500)
    c.cd()
    hSummary.Draw("histo")
    hSummary.SetTitle("; Slope tDiff VS x [ps/cm] ; ")

    text = ROOT.TLatex(0.62, 0.6, "#splitline{#splitline{Histo}{#mu = %.0f #pm %.0f}}{RMS = %.0f #pm %.0f}"%(hSummary.GetMean(), hSummary.GetMeanError(), hSummary.GetRMS(), hSummary.GetRMSError()))
    text.SetNDC()
    text.SetTextSize(0.040)


    gausF = ROOT.TF1("gaus", "gaus", -500, 500)
    gausF.SetRange(hSummary.GetMean()-3*hSummary.GetRMS(), hSummary.GetMean()+3*hSummary.GetRMS())
    hSummary.Fit(gausF, "SR")
    gausF.SetLineColor(2)
    gausF.SetLineStyle(2)
    gausF.Draw("same")
    text.Draw("same")
    c.SaveAs('%s/%s.png'%(outdir,c.GetName()))
    gVsTh.SetPoint(gVsTh.GetN(), refth, abs(gausF.GetParameter(1)))
    gVsTh.SetPointError(gVsTh.GetN()-1, 0, gausF.GetParError(1))
    if refth == 15:
    
        print ("Mean   --   RMS   -- entries --  errMean")
        print ("%.3f  --   %.3f   -- %.0f  --   %.3f"%(hSummary.GetMean(), hSummary.GetRMS(), hSummary.GetEntries() , hSummary.GetRMS()/math.sqrt(hSummary.GetEntries())))


c.Clear()
c.cd()

hdummy.SetTitle("; th ; Slope tDiff VS x [ps/cm]  ")
hdummy.GetXaxis().SetRangeUser(0,30)
#hdummy.GetYaxis().SetRangeUser(-300, -100)
hdummy.GetYaxis().SetRangeUser(100, 300)
hdummy.Draw("")
#gVsTh.Print()

gVsTh.SetLineColor(ROOT.kBlack)
gVsTh.SetMarkerColor(ROOT.kBlack)
gVsTh.SetMarkerStyle(21)
#gVsTh.Fit("pol1", "Q")
gVsTh.Draw("pl")
c.SaveAs('%s/%s_vsTh.png'%(outdir,c.GetName() ))


