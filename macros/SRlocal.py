#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import argparse

import ROOT                                                                                                                                         
import CMS_lumi, tdrstyle                                                                                                                           
                                                                                                                                                    
#set the tdr style                                                                                                                                  
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.05,'X')
ROOT.gStyle.SetLabelSize(0.05,'Y')
ROOT.gStyle.SetTitleSize(0.05,'X')
ROOT.gStyle.SetTitleSize(0.05,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')                                                                                                                
ROOT.gROOT.SetBatch(True)

downEnergyCut15 = {
"00" :"20",
"02" :"3.5",
"03" :"20",
"04" :"5",
"05" :"5",
"06" :"5",
"07" :"50",
"08" :"20",
"09" :"10",
"10" :"30",
"11" :"30",
"12" :"10" ,
"13" :"20" ,
"14" :"20",
}
upEnergyCut15 = {
"00" :"160",
"02" :"10",
"03" :"100",
"04" :"10",
"05" :"10",
"06" :"80",
"07" :"150",
"08" :"100",
"09" :"100",
"10" :"100",
"11" :"150",
"12" :"50" ,
"13" :"100" ,
"14" :"100",
}


downEnergyCut35 = {
"00" :"50",
"02" :"1",
"03" :"20",
"04" :"5",
"05" :"5",
"06" :"0",
"07" :"50",
"08" :"20",
"09" :"10",
"10" :"20",
"11" :"20",
"12" :"0" ,
"13" :"20" ,
"14" :"20",
}

upEnergyCut35 = {
"00" :"50",
"02" :"1",
"03" :"20",
"04" :"5",
"05" :"5",
"06" :"0",
"07" :"50",
"08" :"20",
"09" :"10",
"10" :"20",
"11" :"20",
"12" :"0" ,
"13" :"20" ,
"14" :"20",
}


upEnergyCut = upEnergyCut15
downEnergyCut = downEnergyCut15
vov = '1.50'
run = '5196'
thr = '05'
inputPulse = '/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/pulseShape_HPK_nonIrr_LYSO528_T10C_Vov'+vov+'.root'
inputFile  = '/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/moduleCharacterization_step1_TbToffFit_run'+run+'.root'
outFileName = '/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plots/myprovona.root'
outfile = ROOT.TFile(outFileName, 'RECREATE' )
outdir = '/eos/user/f/fcetorel/www/MTD/TBjune22/run'+run+'/SR_updated/updownEnCut_onlyGoodB_thr_'+thr
#outdir = '/eos/user/f/fcetorel/www/MTD/TBjune22/run'+run+'/updownEnCut_thr_'+thr
#outdir = '/afs/cern.ch/work/f/fcetorel/private/work2/TBJune22/Lab5015Analysis/plottoniProvona/'


infile = ROOT.TFile.Open(inputFile)
inpulse = ROOT.TFile.Open(inputPulse)
print (infile)

os.system("mkdir -p %s"%outdir)
os.system("cp index.php %s"%outdir)


#bars = ["00","02","03","04","05","06","07","08","09","10","11","12","13","14"]
bars = ["00","03","06","07","08","09","10","11","12","13","14"]
SRpulse = []
SRlocal = []
index = 0
gsumR= ROOT.TGraphErrors()
gsumL= ROOT.TGraphErrors()
hlocal = ROOT.TH1F("histo_local", "histo_local", 40, 0, 20)
hpulse = ROOT.TH1F("histo_pulse", "histo_pulse", 40, 0, 20)
for l in ["L","R"]:
  index = 0
  for b in bars:
    tree = infile.Get("data_bar"+b+"L-R_Vov"+vov+"_th"+thr)
    gr = inpulse.Get("g_pulseShape"+l+"_bar"+str(b)+"_Vov"+vov)
    
    c = ROOT.TCanvas()
    if vov == "1.50":    histo = ROOT.TH1F("histo_"+b+"_"+l, "histo_"+b+"_"+l, 60, 0, 30)
    else:    histo = ROOT.TH1F("histo_"+b+"_"+l, "histo_"+b+"_"+l, 50, 0, 50)
    histo.SetTitle("; Local Slew Rate " )
    #tree.Draw("(8*0.313)/tot"+l+">>histo_"+b+"_"+l, "(energyL+energyR)/2 > 10", "histo")
    tree.Draw("(8*0.313)/tot"+l+">>histo_"+b+"_"+l, "(energyL+energyR)/2 > "+downEnergyCut[b]+" && (energyL+energyR)/2 < "+upEnergyCut[b] , "histo")
    histo.SetLineColor(2)
    #gauss35 = ROOT.TF1("","gaus", 5,30)
    #gauss15 = ROOT.TF1("","gaus", 2,12)

    #if vov == "1.50": gauss = gauss15
    #else: gauss = gauss35  
    gauss = ROOT.TF1("", "gaus", histo.GetMean() - 2* histo.GetRMS(), histo.GetMean() + 2* histo.GetRMS())
    gauss.SetParameter(1, histo.GetMean())
    histo.Fit(gauss, "QR")
    SRlocal = gauss.GetParameter(1)
    SRlocalErr = gauss.GetParError(1)
    histo.Draw()
    c.SaveAs(outdir+"/run"+run+"_gaussFit_"+b+l+".png")
    
    c.Clear()

    gr.SetTitle("; time ; amp")
    if thr == "03": lin = ROOT.TF1("","pol1", gr.GetPointX(0), gr.GetPointX(3))
    else: lin = ROOT.TF1("","pol1", gr.GetPointX(1), gr.GetPointX(4))

    lin.SetParameter(1,gr.GetPointY(1)-gr.GetPointY(0) / gr.GetPointX(1)-gr.GetPointX(0) )
    gr.Fit(lin, "QR")
 
    print ( "First attempt :: ", lin.GetParameter(1) )    
    if lin.GetParameter(1) < 0: 
      lin.SetParameters(0.,gr.GetPointY(2)-gr.GetPointY(0) / gr.GetPointX(2)-gr.GetPointX(0) )    
      gr.Fit(lin, "QR")

      print ( "Second attempt :: ", lin.GetParameter(1) )    
      if lin.GetParameter(1) < 0:
        lin.SetParameters(0.,gr.GetPointY(3)-gr.GetPointY(0) / gr.GetPointX(3)-gr.GetPointX(0) )
        gr.Fit(lin, "QR")
        print ( "Third attempt :: ", lin.GetParameter(1) )    
    SRpulse = lin.GetParameter(1)
    SRpulseErr = lin.GetParError(1)
    gr.Draw("AP")


    
    c.SaveAs(outdir+"/run"+run+"_linearFit_"+b+l+".png")
    
    

    # graph with the results
    if ( SRpulse > 0 and SRlocal > 0 and  SRpulse < 40 and SRlocal < 40):
    	if l == "L": 
          gsumL.SetPoint(gsumL.GetN(), SRlocal, SRpulse)
    	  gsumL.SetPointError( gsumL.GetN()-1 , SRlocalErr, SRpulseErr)
    	else: 
          gsumR.SetPoint(gsumR.GetN(), SRlocal, SRpulse)
    	  gsumR.SetPointError( gsumR.GetN()-1 , SRlocalErr, SRpulseErr)

        hlocal.Fill(SRlocal)
        hpulse.Fill(SRpulse)
   
    print (index, " ",  b, " ", l," local: ", SRlocal , "Slew Rate: ", SRpulse)


gsumL.SetTitle(";  SR local ; SR pulse shape")
gsumR.SetMarkerStyle(22)
gsumL.SetMarkerStyle(22)
gsumR.SetMarkerColor(2)
gsumL.SetMarkerColor(4)

gsumL.Draw("AP")
gsumR.Draw("P same")
#if run == "5195":
#  gsum.GetYaxis().SetRangeUser(10,18)
#  gsum.GetXaxis().SetRangeUser(7,21)

gsumMean = ROOT.TGraphErrors()

gsumMean.SetTitle(";  SR local ; SR pulse shape")
for i in range(0,len(bars)):
  gsumMean.SetPoint(gsumMean.GetN(), (gsumL.GetPointX(gsumMean.GetN())+gsumR.GetPointX(gsumMean.GetN()))/2., (gsumL.GetPointY(gsumMean.GetN())+gsumR.GetPointY(gsumMean.GetN()))/2.)
  gsumMean.SetPointError(gsumMean.GetN()-1, (gsumL.GetErrorX(gsumMean.GetN()-1)+gsumR.GetErrorX(gsumMean.GetN()))/2., 
                                            (gsumL.GetErrorY(gsumMean.GetN()-1)+gsumR.GetErrorY(gsumMean.GetN()-1))/2.)
  #print (gsumL.GetPointX(i)  , "       ", gsumR.GetPointX(i) , "      ", (gsumL.GetPointX(i)+gsumR.GetPointX(i))/2. )

f = ROOT.TF1("bisect","x",0,20) 
f.SetLineColor(2)
f.SetLineWidth(2)

f.Draw("same")

c.SaveAs(outdir+"/SRlocalVSpulse_run"+run+".png")

gsumMean.Draw("AP")
f.Draw("same")
c.SaveAs(outdir+"/SRlocalVSpulse_mean_run"+run+".png")

hlocal.Draw("histo")
c.SaveAs(outdir+"/SRlocal_run"+run+".png")

hpulse.Draw("histo")
c.SaveAs(outdir+"/SRpulse_run"+run+".png")

