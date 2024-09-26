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


colors = { 
         'L' : ROOT.kBlue+2, 
         'R' : ROOT.kRed+2, 
         'L-R' : ROOT.kBlack, 
         'allChs' : ROOT.kGreen+1, 
}



parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("-l",  "--label",   required=True, type=str, help="label in the form: HPK_2E14_C25_LYSO815_Vov1.50_T-30C, HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C")
parser.add_argument("-i",  "--inputFolder",  required=True, type=str, help="input folder")
parser.add_argument("-o",  "--outFolder",   required=True, type=str, help="out folder")
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
  angle = 52 # Sept Module were at 49°
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
           #5: [0,1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14, 15], #not particularly good
           #7: [0,1, 2, 3, 4, 5,  6, 7, 8, 9, 10, 11, 12, 13, 14, 15], #only few bars have good MIP peak
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
  #offsetX = 5 #central bar
  offsetX = 0 #central bar


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
  #offsetX = 5 #central bar
  offsetX = 0 #central bar

elif args.label == 'HPK_nonIrr_C25_LYSO818_Vov1.00_T5C':
  goodBars = {
           5: [0,1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15],
           7: [0,1, 2, 3, 4, 5,  8, 9, 10, 11, 12, 13, 14, 15],
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
  #offsetX = 5 #central bar
  offsetX = 0 #central bar
  #6L refbar 7, 9, 10, 11 for th05/07 not fitted
  #refth = 15

barConversionFact = 0.312 / math.cos(angle*math.pi/180) # [cm]


label = args.label
outdir = '%s/%s'%(args.outFolder, label) 
#'/eos/user/f/fcetorel/www/MTD/TBSept23/TOFHIR2C/DPGstudies_4DigitizationModel/energy/Sept24/'%label
inputdir = args.inputFolder #'/afs/cern.ch/work/f/fcetorel/private/work2/dev_TB_CERN_Sept2023/'

#--- prepare outdir
if (os.path.isdir(outdir) == False):
    os.system('mkdir %s'%outdir)
print (outdir)

# Canvas things
hdummy = ROOT.TH2F('hdummy','',100,-30,30,5000,0.,2000)
hdummy.GetXaxis().SetTitle('x [cm]')
hdummy.GetYaxis().SetTitle('Energy [a.u.]')
hdummy.GetXaxis().SetRangeUser(1, 6 )

latex = ROOT.TLatex(0.65,0.84,'%s'%(irradiation))
latex.SetNDC()
latex.SetTextSize(0.038)
latex.SetTextFont(42)

c = ROOT.TCanvas("","",600,500)
#outfile   = ROOT.TFile.Open(outdir+'/uniformityCheck_%s.root'%label,'recreate')
gsVsTh = {}
for l in  ['L', 'R', 'L-R', "allChs"]: gsVsTh[l] = ROOT.TGraphErrors()
graphname = 'energy'
print ('--------------- Doing ', graphname)

for refth, goodbars in goodBars.items():
    print ("Ref is %02d, good bars are: "%refth, goodBars [refth])
    Escale = {} # to store MPV of L-R to normalize E peak
    g_energy_vs_x = OrderedDict()
    for bar in goodbars:
        for l in  ['L', 'R', 'L-R']:
            g_energy_vs_x ['%02d%s'%(bar,l)] = ROOT.TGraphErrors()
    #print (g_energy_vs_x)

    for refbar in range(refbarmin,refbarmax+1): #loop on coincidence bars--> REF bars
        filename = '%s/moduleCharacterization_step2_%s_refbar%d.root'%(inputdir, label, refbar)
        f = ROOT.TFile.Open(filename)
        if f == None: continue
        if args.debug: print (filename)
    
        for bar, gEn in  g_energy_vs_x.items(): #bar is in the form %02d%s 0-15, and s is L,R or L-R
        #for bar in goodbars:
            #print ("DOING BAR %02d"%bar)
            #g_energy_vs_x[bar].Print()
            gnameComplete = 'h1_energy_bar%s_Vov%.02f_th%02d'%(bar, vov, refth)
            if args.debug: print (gnameComplete)
            h1_energy = None
            h1_energy = f.Get(gnameComplete)
            fitFunc = None
            if h1_energy == None: 
                print ("%s not found"%gnameComplete)
                continue
            fitFunc = h1_energy.GetFunction('f_landau_bar%s_Vov%.02f_vth1_%02d'%(bar, vov, refth))
            if fitFunc == None: 
                print ("f_landau_bar%s_Vov%.02f_vth1_%02d not found"%(bar, vov, refth))
                continue
 
            
            gEn.SetPoint(gEn.GetN(), (refbar-offsetX)*barConversionFact, fitFunc.GetParameter(1))
            gEn.SetPointError(gEn.GetN()-1, 0, fitFunc.GetParError(1))
            
    c1 =  ROOT.TCanvas('c_%s_vs_x_bar'%graphname,'c_%s_vs_x_bar'%graphname,600,500)
    
    
    leg = ROOT.TLegend(0.65,0.7,0.92,0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetNColumns(2)
    
    # define 1 graph for th and each channel
    hSummary = {}
    if '815' in args.label: 
        for l in ['L', 'R', 'L-R', 'allChs']: hSummary[l] = ROOT.TH1F("h_%s_refth%02d"%(l,refth), "h_%s_refth%02d"%(l,refth), 45,-0.15, 0.15 )

    elif  '100056' in args.label: 
        for l in ['L', 'L-R', 'allChs']: hSummary[l] = ROOT.TH1F("h_%s_refth%02d"%(l,refth), "h_%s_refth%02d"%(l,refth), 45,-0.15, 0.15 )
        for l in [ 'R' ]: hSummary[l] = ROOT.TH1F("h_%s_refth%02d"%(l,refth), "h_%s_refth%02d"%(l,refth), 60,-0.15, 0.15 )
 
    else: 
        for l in ['L', 'R', 'L-R', 'allChs']: hSummary[l] = ROOT.TH1F("h_%s_refth%02d"%(l,refth), "h_%s_refth%02d"%(l,refth), 60,-0.15, 0.15 )
    tl = ROOT.TLatex() 
    # draw Plots of L-R and do fit to the get MPV for bar
    for bar in goodbars:
        c1.Clear() 
        c1.SetGridy()
        c1.SetGridx()
        c1.cd()
 
        lin0 = ROOT.TF1("pol0_%02d"%(bar), "pol0", -10, 10)        
        hdummy.Draw()
        hdummy.GetYaxis().SetRangeUser(0,1024)
        for l in ["L", "R", "L-R"]:
            #print ("DOING BAR %02d"%bar)
            #g.Print()
            g = g_energy_vs_x ['%02d%s'%(bar,l)]   
            g.SetMarkerStyle(20)
            g.SetMarkerColor(colors[l])
            g.SetMarkerSize(1)
            g.SetLineStyle(1)
            g.SetLineWidth(1)
            g.Draw('p same')
            if l == "L-R":
                lin0.SetRange(g.GetPointX(0)-0.1,g.GetPointX(g.GetN()-1)+0.1 )
                g.Fit(lin0, "QNR")
                Escale [bar] = lin0.GetParameter(0)
                
                tl.DrawLatex(1,900,"%s fit with pol0: %.1f"%(l, lin0.GetParameter(0)))
                lin0.SetLineStyle(ROOT.kDashed)
                #lin0.SetLineWidth(1)
                lin0.SetLineColor(ROOT.kBlack)
                lin0.Draw("same")
    
            c1.SaveAs('%s/%s%02d_ref%02d_beforeNorm.png'%(outdir,c1.GetName(), bar, refth ))
            #c1.SaveAs(outdir+c1.GetName()+'.pdf')
    
    # draw Plots of E normalized to L-R MPV versus X
    for bar in goodbars:
        c1.Clear() 
        hdummy.GetYaxis().SetRangeUser(0.4, 1.6 )
        c1.SetGridy()
        c1.SetGridx()
        c1.cd()
        hdummy.Draw()
        tl = ROOT.TLatex()
        tl.SetNDC()
        tl.SetTextFont(42)
        tl.SetTextSize(0.045)
 
        i = 0 # for the tlatex position
        for l in ["L", "R", "L-R"]:

            g = g_energy_vs_x ['%02d%s'%(bar,l)]   
            if args.debug:
                print ("DOING BAR %02d%s  "%(bar,l))
                print ("Before scaling")
                g.Print()
                print ("---- that's the scaleFactor: %.1f "%Escale [bar])
            lin1 = ROOT.TF1("pol1_%02d%s"%(bar,l), "pol1", -10, 10)        
            lin1.SetRange(g.GetPointX(0)-0.1,g.GetPointX(g.GetN()-1)+0.1 )
 
            if Escale [bar] != 0:  g.Scale(1./ Escale [bar] ) 
            else: print ("Sorry not able to scale, Escale = 0 for bar %02d"%bar)

            if args.debug: 
                print ("after scaling")
                g.Print()

            g.SetMarkerStyle(20)
            g.SetMarkerColor(colors[l])
            g.SetMarkerSize(1)
            g.SetLineStyle(1)
            g.SetLineWidth(1)
            g.Draw('p same')
            g.Fit(lin1, "QNR")

            
            #lin1.SetLineWidth(1)
            lin1.SetLineStyle(ROOT.kDashed)
            lin1.SetLineColor(colors[l])
            lin1.Draw("same")
            tl.DrawLatex(0.20,0.85-i,"%s fit %.4f x [a.u. / cm] + %.4f"%(l, lin1.GetParameter(1), lin1.GetParameter(0)))
            i = i +0.05
            hSummary[l].Fill(lin1.GetParameter(1))
            if l == "L" or l == "R": hSummary['allChs'].Fill(abs(lin1.GetParameter(1)))
        c1.SaveAs('%s/%s%02d_ref%02d.png'%(outdir,c1.GetName(), bar, refth ))
        #c1.SaveAs(outdir+c1.GetName()+'.pdf')
    


    

    #Finally drawing the summary plots of energy slope, and also average value for each th
    for l, hSum in hSummary.items():            
        c =  ROOT.TCanvas('Eslope_summary%s'%l,'Eslope_summary%s'%l,600,500)
        c.cd()
        hSum.Draw("histo")
        hSum.SetTitle("; Slope energy VS x [%/cm] ; ")
        gausF = ROOT.TF1("gaus", "gaus", -500, 500)
        gausF.SetRange(hSum.GetMean()-3*hSum.GetRMS(), hSum.GetMean()+3*hSum.GetRMS())
        hSum.Fit(gausF, "QSRN")
        gausF.SetLineColor(colors[l])
        gausF.SetLineStyle(colors[l])
        hSum.GetXaxis().SetRangeUser(hSum.GetMean()-5*hSum.GetRMS(), hSum.GetMean()+5*hSum.GetRMS())
        gausF.SetRange(gausF.GetParameter(1)-3*gausF.GetParameter(2), gausF.GetParameter(1)+3*gausF.GetParameter(2))
        hSum.Fit(gausF, "QSR")

        gausF.Draw("same")
        c.SaveAs('%s/%s_refTh%02d.png'%(outdir,c.GetName(), refth ))
    
        gVsTh = gsVsTh [l]
        gVsTh.SetPoint(gVsTh.GetN(), refth, gausF.GetParameter(1))
        gVsTh.SetPointError(gVsTh.GetN()-1, 0, gausF.GetParError(1))

# vs th plots
leg = ROOT.TLegend(0.4, 0.75, 0.6, 0.9)
leg.SetNColumns(2)
cTh =  ROOT.TCanvas('Eslope_summary_vsTh','Eslope_summary_vsTh',600,500)
cTh.SetGridx()
cTh.SetGridy()
# Canvas things
hdummy2 = ROOT.TH2F('hdummy','',100,-30,30, 100, -2, 2)
hdummy2.GetXaxis().SetTitle('bar')
hdummy2.GetYaxis().SetTitle('x [cm]')
hdummy2.GetYaxis().SetTitle('Energy [a.u.]')
hdummy2.GetXaxis().SetRangeUser(0, 30 )
hdummy2.GetYaxis().SetRangeUser(-0.1,  0.1)
hdummy2.Draw("")

for l in ["L", "R", "L-R", "allChs"]:
    # VS th
    cTh.cd()
    gVsTh = gsVsTh [l]
    #gVsTh.Print()
    gVsTh.SetLineColor(colors[l])
    gVsTh.SetMarkerColor(colors[l])
    gVsTh.SetMarkerStyle(21)
    gVsTh.Draw("pl")
    leg.AddEntry(gVsTh, l, "pl")
leg.Draw("same")
cTh.SaveAs('%s/%s.png'%(outdir,cTh.GetName()))


