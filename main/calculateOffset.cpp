#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
//#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
//#include "interface/Na22SpectrumAnalyzerModule_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TSpectrum.h"







////---- find energy bins
//void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){
//
//  for(unsigned int i = 1; i < r->size(); i++){
//    TH1F *binHisto = new TH1F ( "binHisto", "binHisto", h -> FindBin(r->at(i)) - h->FindBin(r-> at(i-1)), r-> at(i-1), r->at(i));
//    int j = 1;
//    for (int bin = h->FindBin(r->at(i-1)) ; bin < h -> FindBin(r->at(i))+1 ; bin++){
//      binHisto -> SetBinContent( j, h->GetBinContent(bin));
//      j++;
//    }
//    b[i] = binHisto -> GetMean();
//    binHisto -> Delete();
//  }
//}



// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv)
{
  setTDRStyle();
  float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

  gErrorIgnoreLevel = kError;
  
  typedef std::numeric_limits<double> dbl;
        std::cout.precision(dbl::max_digits10);
  if( argc < 2 )
  {
    std::cout << ">>> calculateOffset::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  

  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  //system(Form("rm -r %s", plotDir.c_str())); // questo non va bene se stiamo lavorando in parallelo
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/tot_fitted_cut/",plotDir.c_str()));
  
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");

  std::string runs = opts.GetOpt<std::string>("Input.runs"); 
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
      map_energyMins[Vov[ii]] = energyMins[ii];
      map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  


  //--- open files
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");

  std::ofstream myfile; // out file
  //  std::ofstream myfileL; // out file


  myfile.open (Form("offset4thr_%s.txt",runs.c_str()));
  //myfileL.open (Form("exampleL_%s.txt",runs.c_str()));

  myfile << "mean,emean,bar,lab,vov,thr\n";
  //myfileL << "mean,emean,bar,lab,vov,thr\n";
  
  std::map<std::string,TTree*> trees;
  
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;  
  
  TList* list = inFile -> GetListOfKeys();
  TIter next(list);
  TObject* object = 0;
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::size_t found;
     
    found = name.find("data_");
    //tree
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    found = name.find("h1_energy_b");
    if( found!=std::string::npos )
    {
     //Vov e th
      std::string stepLabel = tokens[3]+"_"+tokens[4];
      VovLabels[tokens[3]] += 1;
      thLabels[tokens[4]] += 1;
      stepLabels.push_back(stepLabel);
      std::string string_Vov = tokens[3];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      std::string string_th = tokens[4];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  
  
 
  //--- get plot settings
  TCanvas* c;
  TCanvas* c2;
  float* vals = new float[6];
  TLatex* latex;
  TH1F* histo;
  TProfile* prof;
  TH2F* h2;

  
  //------------------
  //--- draw 1st plots
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string Na22 = "Na22";
  std::string Na22SingleBar = "Na22SingleBar";
  std::string Co60 = "Co60";
  std::string Co60SumPeak = "Co60SumPeak";	
  std::string Laser = "Laser";	
  std::string TB = "TB";
  std::string keepAll = "keepAll";
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg
  

  for(auto stepLabel : stepLabels) {
    

    TrackProcess(cpu, mem, vsz, rss);
    
    float Vov = map_Vovs[stepLabel];
    float vth1 = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.2f",Vov));
    std::string thLabel(Form("th%02.0f",vth1));

   
    //--------------------------------------------------------
    // --- loop over bars
    for(int iBar = 0; iBar < 16; ++iBar) {
      
      bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
      if (!barFound) continue;
      
      int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );
      
      // -- loop over L, R, LR
      for(auto LRLabel : LRLabels ) {
	
	//label histo
        std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));

        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
        if (LRLabel == "L-R") { 
          latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	}
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
          

        if (LRLabel == "R" || LRLabel == "L") {

	  // -- draw ToT
          c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
          gPad -> SetLogy();
        
          histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
	  if (histo){
	    histo -> SetTitle(";ToT [ns];entries");
	    histo -> SetLineColor(kRed);
	    histo -> Draw();
            float max = 0.3;
            TF1* fitFunc = new TF1 ( Form("func_%d",index), "gaus", max-10*histo->GetRMS(), max+10*histo->GetRMS() );
            fitFunc -> SetLineColor(kBlack);
            fitFunc -> SetLineWidth(2);
            histo -> Fit( fitFunc, "NQR");
            fitFunc -> SetRange(fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5, fitFunc->GetParameter(1) + fitFunc->GetParameter(2)*5 );
            histo -> Fit( fitFunc, "QRS+");
            fitFunc -> Draw("same");
            float mean = fitFunc -> GetParameter(1);
            float emean = fitFunc -> GetParError(1);
  

	    latex -> Draw("same");      

            myfile << Form("%.4f,%.4f,%02d,%s,%.2f,%d \n",mean,emean,iBar,LRLabel.c_str(),Vov,int(vth1));

	    c -> Print(Form("%s/tot_fitted_cut/c_tot__%s.png",plotDir.c_str(),label.c_str()));
	    c -> Print(Form("%s/tot_fitted_cut/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
	    delete c;
	  }
        }
        
     }// end loop over L, R, L-R labels
      

    }// -- end loop over bars
    
  } // -- end loop over stepLabels
  
  // ---  end 1st plots
  
  
  
  myfile.close();
  //myfileR.close();
}



