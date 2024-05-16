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


# --- colors
cols = { 'refbar 11' : 632,  
         'refbar 2' : 451,  
         'refbar 3' : 888,  
         'refbar 4' : 50,  
         'refbar 5' : 46, 
         'refbar 7' : 51+44, 
         'refbar 6' : 48, 
         'refbar 8' : 51+20, 
         'refbar 9' : 636, 
         'refbar 10' : 52, 
         'refbar 12' : 800, 
         'refbar 1' : 666, 
}
# we can think about adding a parser
# very hardcoded, modify before running, select the module here

parser = argparse.ArgumentParser(description='Module characterization summary plots')
#parser.add_argument("-g",  "--gnames",   required=True, type=str, help="comma-separated list of graphs read from summaryPlots")
parser.add_argument("-l",  "--label",   required=True, type=str, help="label in the form: HPK_2E14_C25_LYSO815_Vov1.50_T-30C, HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C")
parser.add_argument("-b",  "--baseFolder",  required=True, type=str, help="base folder")
parser.add_argument("-o",  "--outFolder",   required=True, type=str, help="out folder")
args = parser.parse_args()

#label = 'HPK_2E14_C25_LYSO815_Vov1.50_T-30C' 
#label = 'HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C'
irradiation = ''
xmin = -0.5
xmax = 15.5
bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
if args.label == 'HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C':
  goodbars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
  vov = 1.00
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = ''
  refth = '05'
elif args.label == 'HPK_nonIrr_C25_LYSO813_Vov3.50_T-30C':
  goodbars = [1,3,4,5,6,8,9,10,11,12,13,14,15]
  vov = 3.50
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = ''
  refth = '05'
elif args.label == 'HPK_2E14_C25_LYSO815_Vov1.50_T-30C':
  goodbars = [1, 2, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 1.50
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  refth = '05'
elif args.label == 'HPK_2E14_C25_LYSO815_Vov1.50_T-35C':
  goodbars = [1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 1.50
  tresmin = 0
  tresmax = 120
  refbarmin = 2
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  refth = '05'
elif args.label == 'HPK_2E14_C25_LYSO815_Vov2.00_T-30C':
  goodbars = [1,2,4,5,7,9,10,11,12,13,14,15]
  vov = 2.00
  tresmin = 0
  tresmax = 120
  refbarmin = 3
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  refth = '05'
elif args.label == 'HPK_2E14_C25_LYSO100056_Vov1.50_T-35C':
  goodbars = [0,1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 1.50
  tresmin = 0
  tresmax = 120
  refbarmin = 3
  refbarmax = 13
  cellsize = '25 #mum'
  irradiation = '2 x 10^{14} 1 MeV n_{eq}/cm^{2}'
  refth = '15'
elif args.label == 'HPK_nonIrr_C25_LYSO818_Vov3.50_T5C':
  goodbars = [0,1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]
  vov = 3.50
  tresmin = 0
  tresmax = 120
  refbarmin = 6
  refbarmax = 12
  cellsize = '25 #mum'
  irradiation = 'non irradiated'
  refth = '15'


#decide the graph to run on here
#graphnames = ['deltaT_totRatioCorr_bestTh_vs_bar', 'deltaT_energyRatioCorr_bestTh_vs_bar' ]    
label = args.label
#graphnames = ['deltaT_energyRatioCorr_bestTh_vs_bar' , 'deltaT_totRatioCorr_bestTh_vs_bar' , 'deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar', 
graphnames = ['deltaT_totRatioCorr_bestTh_vs_bar',
             'energyL_vs_bar', 'energyR_vs_bar', 'energyL-R_vs_bar']    
outdir = args.outFolder +label + '/' #'/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Jan23/%s/'%label
basedir = args.baseFolder #'/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/'
print (outdir)

# Canvas things
hdummy = ROOT.TH2F('hdummy','',100,xmin,xmax,100,0,1000)
hdummy.GetXaxis().SetTitle('bar')
hdummy.GetYaxis().SetTitle('time resolution [ps]')


latex = ROOT.TLatex(0.65,0.84,'%s'%(irradiation))
latex.SetNDC()
latex.SetTextSize(0.038)
latex.SetTextFont(42)

gg_en = OrderedDict() 
g_tRes_average = OrderedDict()
tot = [] 
htotal = ROOT.TH1F("", "", 500,0, 100 )

c = ROOT.TCanvas("","",600,500)
outfile   = ROOT.TFile.Open(outdir+'/uniformityCheck_%s.root'%label,'recreate')

for graphname in graphnames: #main loop on the different graphs (eg. energy, tres...)

  if 'deltaT' in graphname: 
    g_tRes_average[graphname] = ROOT.TGraphErrors()
    g_tRes_average[graphname].SetName('g_tRes_average_'+graphname)
    
  gg = OrderedDict()

  print ('Doing ', graphname)

  for refbar in range(refbarmin,refbarmax): #loop on coincidence bars
    if ('3.50' in args.label and refbar == 7): continue 
    gnameComplete = graphname
    f = ROOT.TFile.Open('%s/%s_refbar%i/summaryPlots_%s_refbar%i.root'%(basedir,label, refbar,label, refbar))
    if f == None: continue
    print ('%s/%s_refbar%i/summaryPlots_%s_refbar%i.root'%(basedir,label, refbar,label, refbar))
    if 'deltaT' in graphname: 
      hdummy.GetYaxis().SetTitle('time resolution [ps]')
      hdummy.GetYaxis().SetRangeUser(tresmin, tresmax)
      gnameComplete = '%s_Vov%.02f_enBin01'%(graphname, vov)
    else: 
      hdummy.GetYaxis().SetTitle('Energy a.u.')
      hdummy.GetYaxis().SetRangeUser(0, 450)
      gnameComplete = '%s_Vov%.02f_th%s'%(graphname, vov,refth)
      print (gnameComplete)
    if f.Get('g_%s'%(gnameComplete)) == None: continue
    gg['refbar %i'%refbar] = f.Get('g_%s'%(gnameComplete))
    gg['refbar %i'%refbar].SetName('g_%s_refbar%i'%(gnameComplete, refbar))

    if 'deltaT' in gnameComplete:

      ## storing the average on all the bars of a module to produce the final plot: tresAve vs coincidence bar
      fitpol0 = ROOT.TF1('fitpol0','pol0',-1,16) 
      gg['refbar %i'%refbar].Fit(fitpol0,'QSN')
      #gg['refbar %i'%refbar].Print()
      ave, err = [fitpol0.GetParameter(0),fitpol0.GetParError(0)]
      print (ave, err)
      g_tRes_average[graphname].SetPoint(g_tRes_average[graphname].GetN(),refbar , ave)
      g_tRes_average[graphname].SetPointError(g_tRes_average[graphname].GetN()-1, 0, err)
    

      c.SetGridy()
      c.SetGridx()
      c.cd()
      hdummy.Draw()
      gg['refbar %i'%refbar].Draw("pl same")
      fitpol0.Draw("same") 

      #CMS_lumi.CMS_lumi(c,  1,  0)
      c.SaveAs(outdir+'g_%s_refbar%i.png'%(gnameComplete, refbar))


    f.Close()  
    
   
    
  # draw cumulative plot  
  hdummy.GetXaxis().SetTitle('bar')
  c1 =  ROOT.TCanvas('c_%s'%(gnameComplete),'c_%s'%(gnameComplete),600,500)
  c1.SetGridy()
  c1.SetGridx()
  c1.cd()
  hdummy.Draw()
  leg = ROOT.TLegend(0.15,0.74,0.80,0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.SetNColumns(2);
  
  for key, g in gg.items():
      g.SetMarkerStyle(20)
      g.SetMarkerSize(1)
      g.SetMarkerColor(cols[key])
      g.SetLineColor(cols[key])
      g.SetLineStyle(1)
      g.SetLineWidth(1)
      g.Draw('plsame')
      leg.AddEntry(g, key, 'PL')
  leg.Draw("same")
  
  c1.SaveAs(outdir+c1.GetName()+'.png')
  c1.SaveAs(outdir+c1.GetName()+'.pdf')


  if 'deltaT' in gnameComplete:
    c1.SetGridx(0)

    hdummy.GetXaxis().SetTitle('coincidence bar')
    hdummy.GetYaxis().SetRangeUser(0, 120)
    hdummy.Draw()
    #leg = ROOT.TLegend(0.15,0.74,0.80,0.92)
    #leg.SetBorderSize(0)
    #leg.SetFillStyle(0)
    #leg.SetNColumns(2);
    
    g_tRes_average[graphname].SetMarkerStyle(20) 
    g_tRes_average[graphname].SetMarkerSize(1)
    g_tRes_average[graphname].SetMarkerColor(601)
    g_tRes_average[graphname].Draw("p") 
    latex.Draw("")
    #CMS_lumi.CMS_lumi(c1,  1,  0)
    c1.SaveAs(outdir+graphname+'_tResAve_vs_coincBar.png')
    outfile.cd()
    g_tRes_average[graphname].Write(g_tRes_average[graphname].GetName())



 
  # this loop to make plots of single bar quantities (energy/tres VS coicidence bar)
  gPerBars = OrderedDict() #to make 1 graph for each DUT bar: time res/ en res vs reference bar
  hPerBars = OrderedDict() #to make 1 histo for each DUT bar: time res/ en res vs reference bar

  for idx,bar in enumerate(goodbars): #NB goodbars != all bars
    print ("index = ", idx, "goodbar = ", bar)
    #if (bar not in goodbars): continue
    #print ("... Is in goodbars :)")
    gdummy = ROOT.TGraphErrors()
    hdum = ROOT.TH1F("", "", 200,tresmin,tresmax)
    for key, g in gg.items():
      gdummy.SetName("g_%s_bar%i_Vs_refbar"%(gnameComplete, bar))

      #NB energy retains all bars, 
      #so to pick up correct value you can extract by bar number
      #while tres has already selection on goodbars hence select directly by idx
      gbar = idx #for tres
      if 'energy' in gnameComplete: gbar = bar
 
      yval = g.GetPointY(gbar) #for energy
      yerr = g.GetErrorY(gbar)
      hdum.SetName("h_%s_bar%i_Vs_refbar"%(gnameComplete, bar))

      gdummy.SetPoint(gdummy.GetN(), int(key[-2:]) , yval ) # since key is a string like this: refbar n
      gdummy.SetPointError(gdummy.GetN()-1, 0, yerr )

      hdum.Fill(yval)
      if ('deltaT_totRatioCorr_bestTh_vs_bar') in gnameComplete: #storing tres to perform final plots
        htotal.Fill(yval)
        tot.append(yval)

      print (gnameComplete, 'key = ', key, 'idx of loop = ', idx, 'bar = ', bar, 'gbar = ', gbar, 'tRes or En',  yval)
    gPerBars['bar%i'%bar] = gdummy
    hPerBars['bar%i'%bar] = hdum

    if (refth) in gnameComplete: #storing energies to produce summary plots later on
      gg_en[gnameComplete.replace('bar','refbar')+'_bar'+str(bar)] = gdummy 

    c =  ROOT.TCanvas('c_%s_bar%i'%(gnameComplete, bar),'c_%s_bar%i'%(gnameComplete, bar),600,500)
    c.SetGridy()
    c.SetGridx()
    c.cd()
    hdummy.Draw()
    gPerBars['bar%i'%bar].Draw('plsame')

    latex.Draw("")
    #CMS_lumi.CMS_lumi(c,  1,  0)
  
    c.SaveAs(outdir+c1.GetName()+'_bar%i.png'%bar)
    c.SaveAs(outdir+c1.GetName()+'_bar%i.pdf'%bar)

    outfile.cd()
    gPerBars['bar%i'%bar].Write(gPerBars['bar%i'%bar].GetName())



    c.Clear()
    c.cd()
    hPerBars['bar%i'%bar].Draw('histo')
    #print (b , hPerBars['bar%i'%b].GetMean(), hPerBars['bar%i'%b].GetRMS())
    gaus = ROOT.TF1("","gaus",tresmin,tresmax)
    gaus.SetLineColor(2)
    gaus.SetParameter(1,  hPerBars['bar%i'%bar].GetMean())
    hPerBars['bar%i'%bar].Fit(gaus, "qR")
    leg = ROOT.TLegend(0.15,0.74,0.80,0.92)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    gaus.Draw("same")    

    #c.SaveAs(outdir+c1.GetName()+'_bar%i_histo.png'%b)
    #c.SaveAs(outdir+c1.GetName()+'_bar%i_histo.pdf'%b)


# Summary of tot as tres - <tres> / <tres>
c =  ROOT.TCanvas('c_b','c_b',600,500)
c.SetGridy()
c.SetGridx()
c.cd()

 
hsummary = ROOT.TH1F("tResSpread", "", 100, -0.5, 0.5)
for t in tot:
  #print (t, htotal.GetMean())
  hsummary.Fill((t - htotal.GetMean() )/ htotal.GetMean() )
  
g = ROOT.TF1("", "gaus", -0.5, 0.5)
g.SetParameter(1, hsummary.GetMean())
hsummary.Fit(g, "q")
hsummary.GetXaxis().SetTitle("(#sigma_{t} - <#sigma_{t}>) / #sigma_{t}  ")
g.SetLineColor(3)  
 
hsummary.Draw('histo') 
g.Draw("same")

latex.Draw("")
#CMS_lumi.CMS_lumi(c,  1,  0)
 
c.SaveAs(outdir+hsummary.GetName()+'.png')
c.SaveAs(outdir+hsummary.GetName()+'.pdf')
outfile.cd()
hsummary.Write(hsummary.GetName())


 
# Summary plots for energy resolution
c.Clear()
c.SetGridy()
c.SetGridx()
c.cd()

 
hdummy.GetXaxis().SetTitle('coincidence bar')
hdummy.GetYaxis().SetTitle('energy [a.u.]')

fitEne = ROOT.TF1('fitEne','pol0',-1,16) 
for i,bar in enumerate(goodbars):
  en_LR = gg_en['energyL-R_vs_refbar_Vov%.02f_th%s_bar%i'%(vov,refth,  bar)]
  en_L = gg_en['energyL_vs_refbar_Vov%.02f_th%s_bar%i'%(vov, refth, bar)]
  en_R = gg_en['energyR_vs_refbar_Vov%.02f_th%s_bar%i'%(vov,refth, bar)]

  outfile.cd()
  en_LR.Write('energyL-R_vs_refbar_Vov%.02f_th%s_bar%i'%(vov, refth, bar))
  en_L.Write('energyL_vs_refbar_Vov%.02f_th%s_bar%i'%(vov, refth,bar))
  en_R.Write('energyR_vs_refbar_Vov%.02f_th%s_bar%i'%(vov,refth, bar))

  hdummy.GetYaxis().SetRangeUser(0, 450)
  hdummy.Draw()

  en_L.SetMarkerColor(633)
  en_R.SetMarkerColor(601)
  en_LR.SetMarkerColor(1)

  en_L.Draw('plsame')
  en_R.Draw('plsame')
  en_LR.Draw('plsame')


  leg = ROOT.TLegend(0.15,0.74,0.80,0.92)
  leg.SetBorderSize(0)
  leg.SetFillStyle(0)
  leg.AddEntry(en_L, 'Left', 'PL')
  leg.AddEntry(en_R, 'Right', 'PL')
  leg.AddEntry(en_LR, 'Average', 'PL')
 

  leg.Draw("same")
 
         
  latex.Draw("")
          
  c.SaveAs(outdir+'energySummary_bar%i.png'%bar)
  c.SaveAs(outdir+'energySummary_bar%i.pdf'%bar)

  c.Clear()
  ROOT.gStyle.SetOptFit(0)
  #Normalizing to mean of L-R
  en_LR.Fit("fitEne")
  aveEne = fitEne.GetParameter(0)

  en_LR.Scale(1./aveEne)
  en_L.Scale(1./aveEne)
  en_R.Scale(1./aveEne)

  en_LR.Write('energyL-RNorm_vs_refbar_Vov%.02f_th05_bar%i'%(vov, bar))
  en_L.Write('energyLNorm_vs_refbar_Vov%.02f_th05_bar%i'%(vov, bar))
  en_R.Write('energyRNorm_vs_refbar_Vov%.02f_th05_bar%i'%(vov,bar))



  hdummy.GetYaxis().SetRangeUser(0, 1.)
  hdummy.Draw()

  en_L.Draw('aplsame')
  en_L.GetXaxis().SetLimits(xmin, xmax)
  en_L.SetTitle("; coincidence bar; energy [a.u.]")
  en_L.GetYaxis().SetRangeUser(0.4, 1.6)
  en_R.Draw('plsame')
  en_LR.Draw('plsame')

  leg.Draw("same")
 
         
  latex.Draw("")
          
  c.SaveAs(outdir+'energySummaryNorm_bar%i.png'%bar)
  c.SaveAs(outdir+'energySummaryNorm_bar%i.pdf'%bar)

  c.Clear() 
