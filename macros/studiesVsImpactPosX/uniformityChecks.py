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
ROOT.gStyle.SetTitleOffset(0.8,'X')
ROOT.gStyle.SetTitleOffset(0.8,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning   


# --- colors
cols = { 
         0 : ROOT.kMagenta+1, 
         1 : ROOT.kOrange+8,  
         2 : ROOT.kOrange+4,  
         3: ROOT.kOrange+2,  
         4: ROOT.kRed+1,  
         5: ROOT.kPink+10, 
         6 : ROOT.kOrange,  
         7: ROOT.kMagenta+3, 
         8: ROOT.kViolet+10, 
         9: ROOT.kBlue+2, 
         10 : ROOT.kAzure+10, 
         11 :ROOT.kCyan, 
         12 :ROOT.kTeal+10, 
         13 :ROOT.kGreen-3, 
         14:ROOT.kGreen+2 , 
         15 :ROOT.kGreen+4,
         0.9: ROOT.kGreen+2 ,
         1.20: ROOT.kBlue+2 ,
}
# fancy palette
cols[4] = ROOT.kPink +1
cols[5] = ROOT.kOrange
cols[6] = ROOT.kRed + 1
cols[7] = ROOT.kMagenta + 2
cols[8] = ROOT.kBlue + 1
cols[9] = ROOT.kCyan
cols[10] = ROOT.kGreen + 2


parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("-g",  "--gname",   required=True, type=str, help="name of the deltaT graph from summaryPlots")
parser.add_argument("-l",  "--label",   required=True, type=str, help="label in the form: HPK_2E14_C25_LYSO815_Vov1.50_T-30C, HPK_nonIrr_C25_LYSO813_Vov1.00_T-30C")
parser.add_argument("-b",  "--baseFolder",  required=True, type=str, help="base folder")
parser.add_argument("-o",  "--outFolder",   required=True, type=str, help="out folder")
parser.add_argument( "--debug",   action='store_true', help="Debugging mode")
args = parser.parse_args()

irradiation = ''
bars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
if args.label == 'DM9001_SM691_Vov0.90_T18C':
  goodbars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
  vov = 0.90
  tresmin = 0
  tresmax = 120
  refbarmin = 4
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = ''
  angle = 52
  barConversionFact = 0.312 / math.cos(angle*math.pi/180)
if args.label == 'DM9001_SM691_Vov1.20_T18C':
  goodbars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
  vov = 1.20
  tresmin = 0
  tresmax = 120
  refbarmin = 4
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = ''
  angle = 52
  barConversionFact = 0.312 / math.cos(angle*math.pi/180)
if args.label == 'DM9003_SM800_Vov0.90_T18C':
  goodbars = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
  vov = 0.90
  tresmin = 0
  tresmax = 120
  refbarmin = 4
  refbarmax = 10
  cellsize = '25 #mum'
  irradiation = ''
  angle = 52
  barConversionFact = 0.312 / math.cos(angle*math.pi/180)


#decide the graph to run on here
label = args.label
graphname = args.gname #'deltaT_totRatioCorr_bestTh_vs_bar'
outdir = args.outFolder +label + '/' #'/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/uniformityStudy_TBpaper_Jan23/%s/'%label
basedir = args.baseFolder #'/eos/user/f/fcetorel/www/MTD/TBMay23/TOFHIR2C/ModuleCharacterization/'
print (outdir)
if args.debug: print( "bar conversion factor =", barConversionFact)
# Canvas things
# vs ref bar plots
hdummy = ROOT.TH2F('hdummy','',100,-0.5,15.5,100,0,120)
hdummy.GetXaxis().SetTitle('bar')
hdummy.GetYaxis().SetTitle('time resolution [ps]')
# vs x in cm plots
hdummyX = ROOT.TH2F('hdummy','',100,0,6,100,0,120)
hdummyX.GetXaxis().SetTitle('Hodoscope x [cm] ')
hdummyX.GetYaxis().SetTitle('time resolution [ps]')


latex = ROOT.TLatex(0.65,0.84,'%s'%(irradiation))
latex.SetNDC()
latex.SetTextSize(0.038)
latex.SetTextFont(42)

# tgraph with average tres (average over 16 bars of a module) vs x position
g_tRes_average = OrderedDict()
# list with all time resolution (per bar of DUT and per REF position)
tot = [] 
# histo to store all the time resolution ((per bar of DUT and per REF position))
htotal = ROOT.TH1F("", "", 200 ,0, 100 ) 

c = ROOT.TCanvas("","",600,500)
outfile   = ROOT.TFile.Open(outdir+'/uniformityCheck_%s.root'%label,'recreate')


g_tRes_average[graphname] = ROOT.TGraphErrors()
g_tRes_average[graphname].SetName('g_tRes_average_'+graphname)
  
gg = OrderedDict()

print ('Doing ', graphname)

for refbar in range(refbarmin,refbarmax+1): #loop on coincidence bars
  gnameComplete = graphname
  fname = '%s/rootFiles/summaryPlots_%s_refbar%i.root'%(basedir,label, refbar)
  f = ROOT.TFile.Open(fname)
  if f == None: 
      print ("File not found:: %s"%f.GetName())
      continue
  gnameComplete = '%s_Vov%.02f_enBin01'%(graphname, vov)
  if f.Get('g_%s'%(gnameComplete)) == None: 
      print ("Graph NOT found:: %s"%(gnameComplete))
      continue
  gg[refbar] = f.Get('g_%s'%(gnameComplete))
  gg[refbar].SetName('g_%s_refbar%i'%(gnameComplete, refbar))
  if args.debug: 
      print ("Doing ref bar %d"%refbar)
      print (fname)
      print (gnameComplete)
      for idx in range(0,len(goodbars)):
          print (idx, " bar: %d , tres: %f "%(gg[refbar].GetPointX(idx), gg[refbar].GetPointY(idx)))


  ## storing the average on all the bars of a module to produce the final plot: tresAve vs coincidence bar
  fitpol0 = ROOT.TF1('fitpol0','pol0',-1,16) 
  fitpol1 = ROOT.TF1('fitpol1','pol1',-1,16) 
  gg[refbar].Fit(fitpol0,'QSN')
  ave, err = [fitpol0.GetParameter(0),fitpol0.GetParError(0)]
  print ('RefBar %02d Average tRes = %.01f, spread (RMS) of tRes = %.01f'%(refbar, fitpol0.GetParameter(0), 100* gg[refbar].GetRMS(2)/ gg[refbar].GetMean(2)))
  g_tRes_average[graphname].SetPoint(g_tRes_average[graphname].GetN(),refbar , ave)
  g_tRes_average[graphname].SetPointError(g_tRes_average[graphname].GetN()-1, 0, err)
  

  c.SetGridy()
  c.cd()
  hdummy.Draw()

  gg[refbar].SetMarkerStyle(21)
  gg[refbar].SetMarkerColor(cols[vov])
  gg[refbar].Draw("pe same")
  fitpol0.SetLineColor(cols[vov]) 
  fitpol0.SetLineStyle(2) 
  fitpol0.Draw("same") 

  c.SaveAs(outdir+'g_%s_refbar%i.png'%(gnameComplete, refbar))

  f.Close()  
  
# draw cumulative plot  
hdummy.GetXaxis().SetTitle('bar')
c1 =  ROOT.TCanvas('c_%s'%(gnameComplete),'c_%s'%(gnameComplete),600,500)
c1.SetGridy()
c1.cd()
hdummy.Draw()
leg = ROOT.TLegend(0.15,0.74,0.80,0.92)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetNColumns(2);

for refbar, g in gg.items():
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1)
    g.SetMarkerColor(cols[refbar])
    g.SetLineColor(cols[refbar])
    g.SetLineStyle(1)
    g.SetLineWidth(1)
    g.Draw('pe same')
    leg.AddEntry(g, "REF bar %d"%refbar, 'PL')
leg.Draw("same")

c1.SaveAs(outdir+c1.GetName()+'.png')
c1.SaveAs(outdir+c1.GetName()+'.pdf')


hdummy.GetXaxis().SetTitle('coincidence bar')
hdummy.GetYaxis().SetRangeUser(0, 120)
hdummy.Draw()

g_tRes_average[graphname].SetMarkerStyle(20) 
g_tRes_average[graphname].SetMarkerSize(1)
g_tRes_average[graphname].SetMarkerColor(601)
g_tRes_average[graphname].Draw("pe") 
latex.Draw("")

c1.SaveAs(outdir+graphname+'_tResAve_vs_coincBar.png')
outfile.cd()
g_tRes_average[graphname].Write(g_tRes_average[graphname].GetName())




# this loop to make plots of single bar quantities (energy/tres VS coicidence bar)
gPerBars = OrderedDict() #to make 1 graph for each DUT bar: time res/ en res vs reference bar
hPerBars = OrderedDict() #to make 1 histo for each DUT bar: time res/ en res vs reference bar

for idx,bar in enumerate(goodbars): #NB goodbars != all bars
  if args.debug: print ("Now filling the tres vs refbar graph, reading from index = ", idx, "goodbar = ", bar)
  gdummy = ROOT.TGraphErrors()
  gdummyX = ROOT.TGraphErrors()
  hdum = ROOT.TH1F("", "", 200,tresmin,tresmax)
  for refbar, g in gg.items():
    gdummy.SetName("g_%s_bar%i_vs_refbar"%(gnameComplete, bar))
    gdummyX.SetName("g_%s_bar%i_vs_x"%(gnameComplete, bar))

    #tres has already selection on goodbars hence select directly by idx
    gbar = idx #for tres

    tres = g.GetPointY(gbar) 
    tres_err = g.GetErrorY(gbar)
    hdum.SetName("h_%s_bar%i_Vs_refbar"%(gnameComplete, bar))

    gdummy.SetPoint(gdummy.GetN(), refbar , tres ) 
    gdummy.SetPointError(gdummy.GetN()-1, 0, tres_err )

    gdummyX.SetPoint(gdummyX.GetN(), float(refbar)*barConversionFact, tres )
    gdummyX.SetPointError(gdummyX.GetN()-1, 0, tres_err )


    hdum.Fill(tres)
    htotal.Fill(tres)
    tot.append(tres)

    if args.debug: print (gnameComplete, 'refbar = ', refbar, ' idx of loop = ', idx, ' bar = ', bar, ' gbar = ', gbar, ' tRes =',  tres)
  gPerBars['bar%i'%bar] = gdummy
  gPerBars['bar%iX'%bar] = gdummyX
  hPerBars['bar%i'%bar] = hdum
  if args.debug:
       print ("VS ref bar graphs")
       for idx in range(0,gdummy.GetN()):
          print (idx, " bar: %d , tres: %f "%(gdummy.GetPointX(idx), gdummy.GetPointY(idx)))

   
  # plot tres of each bar vs ref bar
  c =  ROOT.TCanvas('c_%s_vsREFbar_bar%i'%(gnameComplete, bar),'c_%s_vsREFbar_bar%i'%(gnameComplete, bar),600,500)
  c.SetGridy()
  c.cd()
  hdummy.Draw()
  gPerBars['bar%i'%bar].SetMarkerColor(cols[vov])
  gPerBars['bar%i'%bar].Draw('pesame')
  fitpol0.SetRange(-0.5,15.5)
  gPerBars['bar%i'%bar].Fit(fitpol0, "Q")
  fitpol0.SetLineColor(cols[vov]) 
  fitpol0.SetLineStyle(2) 
  fitpol0.Draw("same") 

  print ('Bar %02d Average tRes = %.01f, spread (RMS) of tRes = %.01f'%(bar, fitpol0.GetParameter(0), 100* gPerBars['bar%i'%bar].GetRMS(2)/ gPerBars['bar%i'%bar].GetMean(2)))

  latex.Draw("")


  c.SaveAs(outdir+c.GetName()+'.png')
  outfile.cd()
  gPerBars['bar%i'%bar].Write(gPerBars['bar%i'%bar].GetName())

  # plot tres of each bar vs x in cm
  c =  ROOT.TCanvas('c_%s_vsX_bar%i'%(gnameComplete, bar),'c_%s_vsX_bar%i'%(gnameComplete, bar),600,500)
  c.SetGridy()
  c.cd()
  hdummyX.Draw()
  gPerBars['bar%iX'%bar].SetMarkerColor(cols[vov])
  gPerBars['bar%iX'%bar].Draw('pesame')
  fitpol1.SetRange(0,6)
  gPerBars['bar%iX'%bar].Fit(fitpol1, "Q")

  fitpol1.SetLineColor(cols[vov]) 
  fitpol1.SetLineStyle(2) 
  fitpol1.Draw("same") 


  latex.Draw("")

  c.SaveAs(outdir+c.GetName()+'.png')

  outfile.cd()
  gPerBars['bar%iX'%bar].Write(gPerBars['bar%iX'%bar].GetName())



# Summary of tRes for all bars and all refs
c =  ROOT.TCanvas('c_b','c_b',600,500)
c.Clear()
c.cd()
ROOT.gStyle.SetOptStat(1)
htotal.Draw("histo")
htotal.GetXaxis().SetTitle('#sigma_{t}')
htotal.GetXaxis().SetRangeUser(htotal.GetMean()-6*htotal.GetRMS(), htotal.GetMean()+6*htotal.GetRMS() )
c.SaveAs(outdir+'tRes_histoAllRefBars_%s.png'%graphname)
c.Clear()


# Summary of tot as tres - <tres> / <tres>
ROOT.gStyle.SetOptStat(1)
hsummary = ROOT.TH1F("tResSpread_"+graphname, "", 60, -0.3, 0.3)
for t in tot:
  # print (t, htotal.GetMean())
  hsummary.Fill((t - htotal.GetMean() )/ htotal.GetMean() )
  
g = ROOT.TF1("", "gaus", hsummary.GetMean()-3*hsummary.GetRMS(), hsummary.GetMean()+3*hsummary.GetRMS())
g.SetParameter(1, hsummary.GetMean())
#hsummary.Fit(g, "QRSN")
#hsummary.Fit(g,"QRS","", g.GetParameter(1)-3*g.GetParameter(2),g.GetParameter(1)+3*g.GetParameter(2) )
#hsummary.Fit(g,"QRS+","", g.GetParameter(1)-3*g.GetParameter(2),g.GetParameter(1)+3*g.GetParameter(2) )
hsummary.GetXaxis().SetTitle("(#sigma_{t} - <#sigma_{t}>) / <#sigma_{t}>  ")
g.SetLineColor(3)  
 
hsummary.Draw('histo') 
#g.Draw("same")
 
c.SaveAs(outdir+hsummary.GetName()+'.png')
c.SaveAs(outdir+hsummary.GetName()+'.pdf')
outfile.cd()
hsummary.Write(hsummary.GetName())
