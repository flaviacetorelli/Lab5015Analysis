#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
 
#include "interface/SetTDRStyle.h"  

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TApplication.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>


using namespace std;
void ilTriangoloNo(float s12, float s23, float s13, float err_s12, float err_s23, float err_s13, std::vector<float> &result, std::vector<float> &err) {
  
  float s1 = sqrt( 0.5 * ( s12*s12 + s13*s13 - s23*s23) );
  float s2 = sqrt( 0.5 * ( s12*s12 + s23*s23 - s13*s13) );
  float s3 = sqrt( 0.5 * ( s13*s13 + s23*s23 - s12*s12) );

  float err_s1 = 1./2/s1 * sqrt( pow(s12*err_s12,2) + pow(s13*err_s13,2) + pow(s23*err_s23,2) );
  float err_s2 = 1./2/s2 * sqrt( pow(s12*err_s12,2) + pow(s23*err_s23,2) + pow(s13*err_s13,2) );
  float err_s3 = 1./2/s3 * sqrt( pow(s13*err_s13,2) + pow(s23*err_s23,2) + pow(s12*err_s12,2) );

  result.push_back(s1);
  result.push_back(s2);
  result.push_back(s3);

  err.push_back(err_s1);
  err.push_back(err_s2);
  err.push_back(err_s3);
}



void triangulation(std::vector<int> barList, map<int,TH1F*> hLR, map<std::string,TH1F*> hsingle, TGraphErrors *g_tRes_Corr_LR, map<std::string,TGraphErrors*> g_tRes_Corr, TGraphErrors *g_tRes_L, TGraphErrors *g_tRes_R, TGraphErrors *g_tRes_chRef) {
  std::cout << "Triangulation ..." << std::endl;
 
  for (int iBar = 0 ; iBar < 16; iBar++){
    bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
    if (!barFound) continue;
 
    if (hLR[iBar] -> GetEntries() == 0) continue;
    if (hsingle[ Form("bar%02dR", iBar)] -> GetEntries() == 0) continue;
    if (hsingle[ Form("bar%02dL", iBar)] -> GetEntries() == 0) continue;

    
    float tRes_L_R = 0;
    float tRes_L_chRef = 0;
    float tRes_R_chRef = 0;

    float err_tRes_L_R = 0;
    float err_tRes_L_chRef = 0;
    float err_tRes_R_chRef = 0;
    
    for ( int i = 0; i < g_tRes_Corr_LR-> GetN(); i++){
      if (  int(g_tRes_Corr_LR -> GetPointX(i)) == iBar) {
	tRes_L_R = g_tRes_Corr_LR-> GetPointY(i);
	err_tRes_L_R = g_tRes_Corr_LR-> GetErrorY(i);
	break;
      }
    }

    for ( int i = 0; i < g_tRes_Corr["L"]-> GetN(); i++){
      if (  g_tRes_Corr["L"]-> GetPointX(i) == iBar) {
	tRes_L_chRef = g_tRes_Corr["L"]-> GetPointY(i);
	err_tRes_L_chRef = g_tRes_Corr["L"]-> GetErrorY(i);
	break;
      }
    }

    for ( int i = 0; i < g_tRes_Corr["R"]-> GetN(); i++){
      if (  g_tRes_Corr["R"]-> GetPointX(i) == iBar) {
	tRes_R_chRef = g_tRes_Corr["R"]-> GetPointY(i);
	err_tRes_R_chRef = g_tRes_Corr["R"]-> GetErrorY(i);
	break;
      }
    }


    //cout <<  "  *** " << tRes_L_R << "   " << tRes_R_chRef  << "   " << tRes_L_chRef <<endl;
    if ( tRes_L_R==0 || tRes_R_chRef==0 || tRes_L_chRef==0) continue;    

    std::vector<float> tRes, err_tRes ;
    tRes.clear();
    err_tRes.clear();

    ilTriangoloNo(tRes_L_R,tRes_R_chRef,tRes_L_chRef, err_tRes_L_R, err_tRes_R_chRef, err_tRes_L_chRef, tRes, err_tRes);
    
    if ( !isnan(tRes[0]) ) g_tRes_L -> SetPoint( g_tRes_L->GetN(), iBar, tRes[0]);
    if ( !isnan(tRes[1]) ) g_tRes_R -> SetPoint( g_tRes_R->GetN(), iBar, tRes[1]);
    if ( !isnan(tRes[2]) ) g_tRes_chRef -> SetPoint( g_tRes_chRef->GetN(), iBar, tRes[2]);

    if ( err_tRes[0]>=0 ) g_tRes_L -> SetPointError( g_tRes_L->GetN()-1, 0, err_tRes[0]);
    if ( err_tRes[1]>=0 ) g_tRes_R -> SetPointError( g_tRes_R->GetN()-1, 0, err_tRes[1]);
    if ( err_tRes[2]>=0 ) g_tRes_chRef -> SetPointError( g_tRes_chRef->GetN()-1, 0, err_tRes[2]);

    //cout <<  "  ***  BELLA DEBUGGING ***  " << endl; 
    //cout <<  "  *** " << tRes[0] << "   " << tRes[1]  << "   " << tRes[2]  <<endl;
    //cout <<  "  *** " << g_tRes_L -> GetPointY( g_tRes_L->GetN()) << "   " << g_tRes_R -> GetPointY( g_tRes_R->GetN())  << "   " << g_tRes_chRef -> GetPointY( g_tRes_chRef->GetN()) <<endl;
   //cout <<  "                             " << endl; 
  }

}





void getTimeResolution(TH1F *histo, TF1 *& fitFunc){

  float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
  float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;

  fitFunc -> SetParameters(histo->GetMaximum(), histo->GetBinCenter(histo->GetMaximumBin()), histo->GetRMS());
  fitFunc -> SetRange(fitXMin, fitXMax);
  histo -> Fit(fitFunc,"QRNL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-1.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.0*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRNL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRSL+");
  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
  histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
 }

void plotTimeResolution(TGraphErrors *g_tRes_L,TGraphErrors *g_tRes_R, TGraphErrors *g_tRes_chRef, float ymin, float ymax, bool useTimeAverage,int color, string plotDir, string label ){
  //cout << "     DEBUG                         "  << label.c_str() << endl;
  //g_tRes_L -> Print();
  //g_tRes_R -> Print();
  //g_tRes_chRef -> Print();
  //cout << "                             "   << endl;
   
  TCanvas *c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  //if ( step1 == 1.50){
  //  ymax = 160;
  //}
  TH2F* hdummy = new TH2F("hdummy","",16,-0.5,15.5,100,ymin,ymax);
  hdummy->GetXaxis()->SetTitle("bar");
  hdummy->GetYaxis()->SetTitle("#sigma_{t} [ps]");
  hdummy->Draw();
  g_tRes_L ->SetName(Form("%s_L", label.c_str()));
  g_tRes_L ->SetMarkerStyle(20);
  g_tRes_L -> SetLineColor(color);
  g_tRes_L -> SetMarkerColor(color);
  g_tRes_R ->SetName(Form("%s_R", label.c_str()));
  g_tRes_R -> SetMarkerStyle(24);
  g_tRes_R -> SetLineColor(color);
  g_tRes_R -> SetMarkerColor(color);
  g_tRes_chRef ->SetName(Form("%s_chRef", label.c_str()));
  g_tRes_chRef ->SetMarkerStyle(21);
  g_tRes_chRef ->SetMarkerColor(12);
  g_tRes_chRef ->SetLineColor(12);
  g_tRes_L->Draw("psame") ;
  g_tRes_R->Draw("psame") ;
  g_tRes_chRef->Draw("psame") ;

  TLegend *leg = new TLegend(0.20, 0.78, 0.35, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g_tRes_L, "L", "P");
  leg->AddEntry(g_tRes_R, "R", "P");
  if (useTimeAverage) leg->AddEntry(g_tRes_chRef, "Ref (average)", "P");
  else leg->AddEntry(g_tRes_chRef, "Ref (single)", "P");
  leg->Draw("same");

  TF1 *myfitL = new TF1("myfitL","pol0", 0, 100);
  myfitL -> SetLineColor(color);
  myfitL -> SetParameter(0, g_tRes_L -> GetMean(2));
  g_tRes_L ->Fit("myfitL");
  TLatex *latexL = new TLatex(0.50,0.88,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitL->GetParameter(0), myfitL->GetParError(0)));
  latexL -> SetNDC();
  latexL -> SetTextFont(42);
  latexL -> SetTextSize(0.04);
  latexL->Draw("same");
  TF1 *myfitR = new TF1("myfitR","pol0", 0, 100);
  myfitR -> SetLineColor(color);
  myfitR -> SetLineStyle(2);
  myfitR -> SetParameter(0, g_tRes_R -> GetMean(2));
  g_tRes_R ->Fit("myfitR");
  TLatex *latexR = new TLatex(0.50,0.83,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitR->GetParameter(0), myfitR->GetParError(0)));
  latexR -> SetNDC();
  latexR -> SetTextFont(42);
  latexR -> SetTextSize(0.04);
  latexR->Draw("same");
  TF1 *myfitChRef = new TF1("myfitChRef","pol0", 0, 100);
  myfitChRef -> SetLineColor(12);
  myfitChRef -> SetParameter(0,g_tRes_chRef -> GetMean(2));
  g_tRes_chRef ->Fit("myfitChRef");
  TLatex *latexRef = new TLatex(0.50,0.78,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitChRef->GetParameter(0), myfitChRef->GetParError(0)));
  latexRef -> SetNDC();
  latexRef -> SetTextFont(42);
  latexRef -> SetTextSize(0.04);
  latexRef->Draw("same");
  
  std::cout << "L    : " <<  "   <sigma_t> = " << myfitL->GetParameter(0) << "+/-" << myfitL->GetParError(0) <<std::endl;
  std::cout << "R    : " <<  "   <sigma_t> = " << myfitR->GetParameter(0) << "+/-" << myfitR->GetParError(0) <<std::endl;
  std::cout << "ChRef: " <<  "   <sigma_t> = " << myfitChRef->GetParameter(0) << "+/-" << myfitChRef->GetParError(0) <<std::endl;

  c->Print(Form("%s/c_tRes_%s.png",plotDir.c_str(), label.c_str()));
  c->Print(Form("%s/c_tRes_%s.pdf",plotDir.c_str(), label.c_str()));
  c->SaveAs(Form("%s/c_tRes_%s.root",plotDir.c_str(), label.c_str()));
 // delete c;
 // delete myfitL;
 // delete myfitR;
 // delete myfitChRef;
 // delete latexL;
 // delete latexR;
 // delete latexRef;

}


// ================================================
int main(int argc, char** argv){
  
  setTDRStyle();
  gErrorIgnoreLevel = kError;  

  if( argc < 2 )
    {
      std::cout << ">>> drawPulseShape::usage:   " << argv[0] << " configFile.cfg" << std::endl;
      return -1;
    }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int chRef1 = opts.GetOpt<int>("Input.chRef1");
  int chRef2 = opts.GetOpt<int>("Input.chRef2");
  int useTimeAverage = opts.GetOpt<int>("Input.useTimeAverage"); 
  std::cout<< "Reference channels :  " <<  chRef1 << "  " << chRef2 <<std::endl;
  if ( useTimeAverage ) std::cout << "Using tAverage of the front module as reference"<<std::endl; 

  std::string outDir = opts.GetOpt<std::string>("Output.outDir");
  std::string mainPlotDir = opts.GetOpt<std::string>("Output.plotDir");


  float ymin = opts.GetOpt<float>("Plots.ymin");
  float ymax = opts.GetOpt<float>("Plots.ymax");
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList"); 

 
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  int array = opts.GetOpt<int>("Channels.array"); 

  std::vector<std::string> labelLR = {"L","R"};
  map<std::string, int> chID;  
  for(int iBar = 0; iBar < 16; ++iBar){
    for (auto label : labelLR ){
      if ( label == "L") chID[ Form("bar%02d%s", iBar, label.c_str()) ] = channelMapping[iBar*2+0]+array*64;
      if ( label == "R") chID[ Form("bar%02d%s", iBar, label.c_str()) ] = channelMapping[iBar*2+1]+array*64;
    }
  }

  
  for(int iBar = 0; iBar < 16; ++iBar){
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str()); 
      std::cout << iBar << "  " << label << "  " << chID[chLabel] <<   std::endl;
    }
  }

  int maxActiveChannels =  opts.GetOpt<int>("Cuts.maxActiveChannels");
  float minEnergy = opts.GetOpt<int>("Cuts.minEnergy");
  float maxEnergy = 850;
  int mystep2 = opts.GetOpt<int>("Cuts.step2"); 

  std::string offsetFileName = opts.GetOpt<std::string>("Input.offsetFileName");
  std::map < std::pair<int, std::string>, float> offsetCh;

  if (offsetFileName != "NOOFFSET") {
    std::cout<< " >>>>> Calculating TbT with the offset provided: \n " <<  offsetFileName << std::endl;
    std::ifstream offsetFile;
    offsetFile.open(offsetFileName);
    std::string line;
    int bar;
    float ov;
    std::string side;
    float value;
    while ( offsetFile.good() ){
      getline(offsetFile, line);
      std::istringstream ss(line);
      ss >> bar >> ov >> side >> value;
      offsetCh[std::make_pair(bar,side)] = value * 1000; // from nano to pico
      std::cout<< bar <<  "   " << side << "  " << offsetCh[std::make_pair(bar,side)] <<std::endl;
    }
  }
  else{
    std::cout<< "!!!!   NO offset provided" <<std::endl;
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
       for (auto label : labelLR ){
        offsetCh[std::make_pair(iBar, label)] = 0.;
      }
    }
  }

  

  
  // --- reading tree
  //------------------------------
  TChain* data = new TChain("data","data");
  
  std::stringstream ss(runs); 
  std::string token;
  while( std::getline(ss,token,',') )
    {
      std::stringstream ss2(token);
      std::string token2;
      int runMin = -1;
      int runMax = -1;
      while( std::getline(ss2,token2,'-') )
	{
	  if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
	  if( runMin == -1 ) runMin = atoi(token2.c_str());
	}
      if( runMax == -1 ) runMax = runMin;
      
      for(int run = runMin; run <= runMax; ++run) {
	//std::string inFileName = Form("/data1/cmsdaq/tofhir2/h8/reco/%04d/*_e.root",run);
	std::string inFileName = Form("%s/%04d/*_e.root",inputDir.c_str(), run);
	std::cout << ">>> Adding file " << inFileName << std::endl;
	data -> Add(inFileName.c_str());
      }
    }

  //--- define branches
  float step1, step2;
  int channelIdx[128];
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  std::vector<long long> *time = 0;
  std::vector<unsigned short>* t1fine = 0;
  
  data -> SetBranchStatus("*",0);
  data -> SetBranchStatus("step1",  1); data -> SetBranchAddress("step1",  &step1);
  data -> SetBranchStatus("step2",  1); data -> SetBranchAddress("step2",  &step2);
  data -> SetBranchStatus("channelIdx",  1); data -> SetBranchAddress("channelIdx",  channelIdx);
  data -> SetBranchStatus("tot",    1); data -> SetBranchAddress("tot",       &tot);
  data -> SetBranchStatus("energy", 1); data -> SetBranchAddress("energy", &energy);
  data -> SetBranchStatus("time",   1); data -> SetBranchAddress("time",     &time);
  data -> SetBranchStatus("t1fine",   1); data -> SetBranchAddress("t1fine",     &t1fine);

  int nEntries = data->GetEntries();
  cout << "Number of entries = " << nEntries << endl;
  //  int maxEntries = 200000;
  int maxEntries = nEntries;
  
  //int maxEntries = 2000000;


  // -- book histograms 
  TH1F *h_nActiveChannels0 = new TH1F("h_nActiveChannels0","h_nActiveChannels0",32,-0.5,31.5);
  TH1F *h_nActiveChannels1 = new TH1F("h_nActiveChannels1","h_nActiveChannels1",32,-0.5,31.5);
  
  TH1F *h_energy_chRef = new TH1F("h_energy_chRef","h_energy_chRef",512,0,1024);
  TH1F *h_tot_chRef = new TH1F("h_tot_chRef","h_tot_chRef",512,0,1024);

  map<std::string,TH1F*>      h_energy;
  map<std::string,TH1F*>      h_energyRatio;
  map<std::string,TH2F*>      h2_energyRatio_vs_totRatio;
  map<std::string,TH1F*>      h_deltaT;
  map<std::string,TH1F*>      h_deltaT_energyRatioCorr;
  map<std::string,TH1F*>      h_deltaT_energyRatioCorr_phaseCorr;
  map<std::string,TProfile*>  p_deltaT_vs_energyRatio;
  map<std::string,TH2F*>      h2_deltaT_vs_energyRatio;
  map<std::string,TProfile*>  p_deltaT_energyRatioCorr_vs_t1fine;
  map<std::string,TH2F*>      h2_deltaT_energyRatioCorr_vs_t1fine;
  map<std::string,TProfile*>  p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef;
  map<std::string,TH2F*>      h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef;

  map<std::string,TH1F*>      h_tot;
  map<std::string,TH1F*>      h_totRatio;
  map<std::string,TH1F*>      h_deltaT_energyRatioCorr_totRatioCorr;
  map<std::string,TH1F*>      h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr;
  map<std::string,TProfile*>  p_deltaT_energyRatioCorr_vs_totRatio;
  map<std::string,TH2F*>      h2_deltaT_energyRatioCorr_vs_totRatio;
  map<std::string,TProfile*>  p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine;
  map<std::string,TH2F*>      h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine;
  map<std::string,TProfile*>  p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef;
  map<std::string,TH2F*>      h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef;




  map<int,TH1F*>      h_energyRatio_LR;
  map<int,TH1F*>      h_deltaT_LR;
  map<int,TH1F*>      h_deltaT_LR_energyRatioCorr;
  map<int,TH1F*>      h_deltaT_LR_energyRatioCorr_phaseCorr;
  map<int,TProfile*>  p_deltaT_LR_vs_energyRatio;
  map<int,TProfile*>  p_deltaT_LR_energyRatioCorr_vs_t1fine;

  map<int,TH1F*>      h_totRatio_LR;
  map<int,TProfile*>  p_deltaT_LR_energyRatioCorr_vs_totRatio;
  map<int,TH1F*>      h_deltaT_LR_energyRatioCorr_totRatioCorr;
  map<int,TProfile*>  p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine;
  map<int,TH1F*>      h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr;

  for(int iBar = 0; iBar < 16; ++iBar){        
    
    // L-R
    // energy related
    h_energyRatio_LR[iBar] = new TH1F(Form("h_energyRatio_LR_bar%02d", iBar), Form("h_energyRatio_LR_bar%02d", iBar), 100, 0, 3);
    h_deltaT_LR[iBar] = new TH1F(Form("h_deltaT_LR_bar%02d", iBar), Form("h_deltaT_LR_bar%02d", iBar), 1000, -12000, 12000);
    h_deltaT_LR_energyRatioCorr[iBar] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_bar%02d", iBar), Form("h_deltaT_LR_energyRatioCorr_bar%02d", iBar), 1000, -12000, 12000);
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d", iBar), Form("h_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d", iBar), 1000, -12000, 12000);
    p_deltaT_LR_vs_energyRatio[iBar] = new TProfile(Form("p_deltaT_LR_vs_energyRatio_bar%02d",iBar), Form("p_deltaT_LR_vs_energyRatio_bar%02d",iBar), 60, 0, 3);
    p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar] = new TProfile(Form("p_deltaT_LR_energyRatioCorr_vs_t1fine_bar%02d",iBar), Form("p_deltaT_LR_energyRatioCorr_vs_t1fine_bar%02d",iBar), 50, 0, 1000);

    //tot related
    h_totRatio_LR[iBar] = new TH1F(Form("h_totRatio_LR_bar%02d", iBar), Form("h_totRatio_LR_bar%02d", iBar), 100, 0, 3);
    h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_totRatioCorr_bar%02d", iBar), Form("h_deltaT_LR_energyRatioCorr_totRatioCorr_bar%02d", iBar), 1000, -12000, 12000);
    p_deltaT_LR_energyRatioCorr_vs_totRatio[iBar] = new TProfile(Form("p_deltaT_LR_energyRatioCorr_vs_totRatio_bar%02d",iBar), Form("p_deltaT_LR_energyRatioCorr_vs_totRatio_bar%02d",iBar), 60, 0, 3);
    h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr_bar%02d", iBar), Form("h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr_bar%02d", iBar), 1000, -12000, 12000);
    p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar] = new TProfile(Form("p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine_bar%02d",iBar), Form("p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine_bar%02d",iBar), 50, 0, 1000);
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      
      // energy related
      h_energy[chLabel] = new TH1F(Form("h_energy_%s", chLabel.c_str()) , Form("h_energy_%s", chLabel.c_str()), 512, 0, 1024);
      h_energyRatio[chLabel] = new TH1F(Form("h_energyRatio_%s", chLabel.c_str()) , Form("h_energyRatio_%s", chLabel.c_str()), 100, 0, 3);
      h2_energyRatio_vs_totRatio[chLabel] = new TH2F(Form("h2_energyRatio_vs_totRatio_%s", chLabel.c_str()) , Form("h_energyRatio_vs_totRatio_%s", chLabel.c_str()), 100, 0, 3,100,0,3);
      h_deltaT[chLabel] = new TH1F(Form("h_deltaT_%s", chLabel.c_str()) , Form("h_deltaT_%s", chLabel.c_str()), 1000, -12000, 12000);
      h_deltaT_energyRatioCorr[chLabel] = new TH1F(Form("h_deltaT_energyRatioCorr_%s", chLabel.c_str()) , Form("h_deltaT_energyRatioCorr_%s", chLabel.c_str()), 1000, -12000, 12000);
      h_deltaT_energyRatioCorr_phaseCorr[chLabel] = new TH1F(Form("h_deltaT_energyRatioCorr_phaseCorr_%s", chLabel.c_str()) , Form("h_deltaT_energyRatioCorr_phaseCorr_%s", chLabel.c_str()), 1000, -12000, 12000);
      p_deltaT_vs_energyRatio[chLabel] = new TProfile(Form("p_deltaT_vs_energyRatio_%s", chLabel.c_str()) , Form("p_deltaT_vs_energyRatio_%s", chLabel.c_str()), 60, 0, 3);
      h2_deltaT_vs_energyRatio[chLabel] = new TH2F(Form("h2_deltaT_vs_energyRatio_%s", chLabel.c_str()) , Form("h2_deltaT_vs_energyRatio_%s", chLabel.c_str()), 60, 0, 3, 1000, -12000, 12000);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel] = new TProfile(Form("p_deltaT_energyRatioCorr_vs_t1fine_%s", chLabel.c_str()) , Form("p_deltaT_energyRatioCorr_vs_t1fine_%s", chLabel.c_str()), 50, 0, 1000);
      h2_deltaT_energyRatioCorr_vs_t1fine[chLabel] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_t1fine_%s", chLabel.c_str()) , Form("h2_deltaT_energyRatioCorr_vs_t1fine_%s", chLabel.c_str()), 50, 0, 1000, 1000, -12000, 12000);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel] = new TProfile(Form("p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()) , Form("p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()), 50, 0, 1000);
      h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel] = new TH2F(Form("h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()) , Form("h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()), 50, 0, 1000, 1000, -12000, 12000);


      // tot related !!! remember to set the right RANGES
      h_tot[chLabel] = new TH1F(Form("h_tot_%s", chLabel.c_str()) , Form("h_tot_%s", chLabel.c_str()), 100, 0, 1000);
      h_totRatio[chLabel] = new TH1F(Form("h_totRatio_%s", chLabel.c_str()) , Form("h_totRatio_%s", chLabel.c_str()), 100, 0, 3);
      h_deltaT_energyRatioCorr_totRatioCorr[chLabel] = new TH1F(Form("h_deltaT_energyRatioCorr_totRatioCorr_%s", chLabel.c_str()) , Form("h_deltaT_energyRatioCorr_totRatioCorr_%s", chLabel.c_str()), 1000, -12000, 12000);
      h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[chLabel] = new TH1F(Form("h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s", chLabel.c_str()) , Form("h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s", chLabel.c_str()), 1000, -12000, 12000);
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel] = new TProfile(Form("p_deltaT_energyRatioCorr_vs_totRatio_%s", chLabel.c_str()) , Form("p_deltaT_energyRatioCorr_vs_totRatio_%s", chLabel.c_str()), 60, 0, 3);
      h2_deltaT_energyRatioCorr_vs_totRatio[chLabel] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_totRatio_%s", chLabel.c_str()) , Form("h2_deltaT_energyRatioCorr_vs_totRatio_%s", chLabel.c_str()), 60, 0, 3, 1000, -12000, 12000);
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] = new TProfile(Form("p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s", chLabel.c_str()) , Form("p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s", chLabel.c_str()), 50, 0, 1000);
      h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] = new TH2F(Form("h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s", chLabel.c_str()) , Form("h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s", chLabel.c_str()), 50, 0, 1000, 1000, -12000, 12000);
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel] = new TProfile(Form("p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()) , Form("p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()), 50, 0, 1000);
      h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel] = new TH2F(Form("h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()) , Form("h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s", chLabel.c_str()), 50, 0, 1000, 1000, -12000, 12000);
 

    }
  }

  TH1F *h_deltaT_LR_chRef = new TH1F("h_deltaT_LR_chRef","h_deltaT_LR_chRef",  1000, -12000, 12000); 

  map <int, map<int, bool> > acceptEvent;
  map<int, bool>  acceptEvent_chRef;

  map<int, int>  nActiveChannels0;
  map<int, int>  nActiveChannels1;

  // -- first loop over events
  cout << "First loop over events to find the mip peak" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    data->GetEntry(entry);

    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;

    if (step2 != mystep2) continue;
    
    // -- count active channels in the two modules
    nActiveChannels0[entry] = 0;
    for (int ch = 0; ch < 32; ch++) {
      if ( channelIdx[ch] < 0 ) continue;
      if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
      if ((*energy)[channelIdx[ch]] > 0)
	nActiveChannels0[entry]+=1;
    }    
    
    nActiveChannels1[entry] = 0;
    for (int ch = 64; ch < 96 ; ch++) {
      if ( channelIdx[ch] < 0 ) continue;
      if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
      if ((*energy)[channelIdx[ch]] > 0)
	nActiveChannels1[entry]+=1;
    }

    h_nActiveChannels0 ->Fill(nActiveChannels0[entry]);
    h_nActiveChannels1 ->Fill(nActiveChannels1[entry]);

    if (nActiveChannels0[entry] >  maxActiveChannels) continue;
    if (nActiveChannels1[entry] >  maxActiveChannels) continue;


    //-- ref channel
    if ( channelIdx[chRef1] < 0 ) continue;
    if ( channelIdx[chRef2] < 0 ) continue;
    if ( (*tot)[channelIdx[chRef1]]/1000 < -10. || (*tot)[channelIdx[chRef1]]/1000 > 50. ) continue;      
    if ( (*tot)[channelIdx[chRef2]]/1000 < -10. || (*tot)[channelIdx[chRef1]]/1000 > 50. ) continue;      

    float energyRef = 0.5* ((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]] );
    if ( !useTimeAverage) energyRef  = (*energy)[channelIdx[chRef1]];

    h_energy_chRef -> Fill( energyRef );

    

    for(int iBar = 0; iBar < 16; ++iBar){        
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( channelIdx[ch] < 0 ) continue;
	if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
	
	h_energy[chLabel] -> Fill( (*energy)[channelIdx[ch]] );	
      }
    }

  }// -- end first loop over entries


  // -- find min energy
  TF1 *fitLandau_chRef = new TF1("fitLandau_chRef","landau", 0, 1000);
  h_energy_chRef->GetXaxis()->SetRangeUser(200,800);
  int maxbin = h_energy_chRef->GetMaximumBin();
  float peak = h_energy_chRef->GetBinCenter(maxbin);
  fitLandau_chRef->SetRange(peak*0.8, peak*1.2);
  fitLandau_chRef->SetParameter(1,peak);
  fitLandau_chRef->SetParameter(2,0.1*peak);
  h_energy_chRef->Fit("fitLandau_chRef","QR");
  float energyMin_chRef = fitLandau_chRef->GetParameter(1)*0.8;
  h_energy_chRef->GetXaxis()->SetRangeUser(0,1000);

  map<std::string, float> energyMin;
  map<std::string, TF1*> fitLandau;

  for(int iBar = 0; iBar < 16; ++iBar){        
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      
      fitLandau[chLabel] = new TF1(Form("fitLandau_%s",chLabel.c_str()),"landau", 0, 1000);    
      h_energy[chLabel]->GetXaxis()->SetRangeUser(minEnergy,800);
      if (chLabel == "bar13L") h_energy[chLabel]->GetXaxis()->SetRangeUser(minEnergy/5,800);
      if (chLabel == "bar13L" && runs == "5352") h_energy[chLabel]->GetXaxis()->SetRangeUser(60,800);
      if (chLabel == "bar06L" && runs == "5299") h_energy[chLabel]->GetXaxis()->SetRangeUser(40,800);
      int maxbin = h_energy[chLabel]->GetMaximumBin();
      float peak = h_energy[chLabel]->GetBinCenter(maxbin);
      fitLandau[chLabel]->SetRange(peak*0.85, peak*1.2);
      fitLandau[chLabel]->SetParameter(1,peak);
      fitLandau[chLabel]->SetParameter(2,0.1*peak);
      h_energy[chLabel]->Fit(fitLandau[chLabel],"QR");
      fitLandau[chLabel]->SetRange(fitLandau[chLabel]->GetParameter(1)*0.85, fitLandau[chLabel]->GetParameter(1)*1.2);  
      energyMin[chLabel] = fitLandau[chLabel]->GetParameter(1)*0.75;
      h_energy[chLabel]->GetXaxis()->SetRangeUser(0,1000);
    }
  }


  // -- second loop over events to get amp walk corrections
  cout << "Second loop over events to get amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry);

    acceptEvent_chRef[entry] = false;

    if (step2 != mystep2) continue;
   
    // -- remove showering events
    if (nActiveChannels0[entry] > maxActiveChannels) continue;
    if (nActiveChannels1[entry] > maxActiveChannels) continue;

    //-- ref channel
    if ( channelIdx[chRef1] < 0 ) continue;
    if ( channelIdx[chRef2] < 0 ) continue;
    if ( (*tot)[channelIdx[chRef1]]/1000 < -10. || (*tot)[channelIdx[chRef1]]/1000 > 50. ) continue;      
    if ( (*tot)[channelIdx[chRef2]]/1000 < -10. || (*tot)[channelIdx[chRef2]]/1000 > 50. ) continue;      

    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    float totRef = 0.5*((*tot)[channelIdx[chRef1]]+(*tot)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);

    if ( !useTimeAverage){
      energyRef  = (*energy)[channelIdx[chRef1]];
      totRef  = (*tot)[channelIdx[chRef1]];
      tRef  = (*time)[channelIdx[chRef1]];
    }

    if ( energyRef < energyMin_chRef || energyRef > 850) continue;

    acceptEvent_chRef[entry] = true;

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 

	acceptEvent[entry][ch] = false;

	if ( runs == "5352" && chLabel == "bar00R")  maxEnergy = 600.;
	else if ( runs == "5352" && chLabel == "bar01L")  maxEnergy = 700.;
	else maxEnergy = 850;

	if ( channelIdx[ch] < 0 ) continue;
	if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
	if ( (*energy)[channelIdx[ch]] < energyMin[chLabel] || (*energy)[channelIdx[ch]] > maxEnergy ) continue;
	
	acceptEvent[entry][ch] = true;
	
	long long deltaT = (*time)[channelIdx[ch]] - tRef;
      
	if ( fabs(deltaT)>10000) continue;
	
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef ;
	float totRatio = ((*tot)[channelIdx[ch]]-offsetCh[std::make_pair(iBar , label.c_str())])/totRef ;
        // cout << iBar << " " << label.c_str() << " " << (*tot)[channelIdx[ch]] << "  - " << offsetCh[std::make_pair(iBar , label.c_str())] <<   " = " << totRatio << endl; 
	
	h_deltaT[chLabel]   -> Fill( deltaT );	
	h_energyRatio[chLabel] -> Fill( energyRatio );	
	h_totRatio[chLabel] -> Fill( totRatio );	
	h2_energyRatio_vs_totRatio[chLabel] -> Fill( totRatio, energyRatio );	
	p_deltaT_vs_energyRatio[chLabel] -> Fill( energyRatio , deltaT );	
	h2_deltaT_vs_energyRatio[chLabel] -> Fill( energyRatio , deltaT );
      }

      // tDiff
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      if ( acceptEvent[entry][chL]  && acceptEvent[entry][chR]  && acceptEvent_chRef[entry])    {      


	float totRatio = ((*tot)[channelIdx[chL]]-offsetCh[std::make_pair(iBar , "L")])/((*tot)[channelIdx[chR]]-offsetCh[std::make_pair(iBar , "R")]);
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]]; 
        h_energyRatio_LR[iBar]-> Fill( energyRatio);
        h_totRatio_LR[iBar]-> Fill( totRatio );
	
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	if ( fabs(deltaT)<10000  ){
	  h_deltaT_LR[iBar]-> Fill( deltaT );
	  p_deltaT_LR_vs_energyRatio[iBar]-> Fill((*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]], deltaT );
	}
      }
    
    }// end loop over bars

    // fill deltaT L-R for ref bar
    int ch1 =  0+array*64;
    int ch2 = 31+array*64;
    if ( acceptEvent[entry][ch1]  && acceptEvent[entry][ch2]  && acceptEvent_chRef[entry])    {      
      h_deltaT_LR_chRef ->Fill(  (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] );
    }

  }// -- end second loop over entries
  
  
  // ---  amp walk corr
  map<std::string,TF1*> fitFun_energyRatio;
  map<std::string,TF1*> fitFun_energyRatioCorr;

  map<int,TF1*> fitFun_energyRatio_LR;
  map<int,TF1*> fitFun_energyRatioCorr_LR;

  for(int iBar = 0; iBar < 16; ++iBar){  

    // L-R
    fitFun_energyRatio_LR[iBar] = new TF1(Form("fitFun_energyRatio_LR_%02d", iBar), "gaus", 0,10);  
    h_energyRatio_LR[iBar] -> Fit(fitFun_energyRatio_LR[iBar],"QR");
    
    fitFun_energyRatioCorr_LR[iBar] = new TF1(Form("fitFun_energyRatioCorr_LR_bar%02d", iBar), "pol3", 0,10);
    fitFun_energyRatioCorr_LR[iBar]->SetRange( fitFun_energyRatio_LR[iBar]->GetParameter(1)-5*fitFun_energyRatio_LR[iBar]->GetParameter(2), fitFun_energyRatio_LR[iBar]->GetParameter(1)+5*fitFun_energyRatio_LR[iBar]->GetParameter(2));
    p_deltaT_LR_vs_energyRatio[iBar] -> Fit(fitFun_energyRatioCorr_LR[iBar],"QRS"); 


    // - single channels 
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      int ch = chID[chLabel]; 

      fitFun_energyRatio[chLabel] = new TF1(Form("fitFun_energyRatio_%s",chLabel.c_str()), "gaus", 0,10);  
      h_energyRatio[chLabel] -> Fit(fitFun_energyRatio[chLabel],"QR");

      fitFun_energyRatioCorr[chLabel] = new TF1(Form("fitFun_energyRatioCorr_ch%02d",ch), "pol3", 0,10);
      //fitFun_energyRatioCorr[chLabel]->SetRange( fitFun_energyRatio[chLabel]->GetParameter(1) - 3*fitFun_energyRatio[chLabel]->GetParameter(2), fitFun_energyRatio[chLabel]->GetParameter(1) + 3*fitFun_energyRatio[chLabel]->GetParameter(2) );
      fitFun_energyRatioCorr[chLabel]->SetRange( fitFun_energyRatio[chLabel]->GetParameter(1) - 5*fitFun_energyRatio[chLabel]->GetParameter(2), fitFun_energyRatio[chLabel]->GetParameter(1) + 5*fitFun_energyRatio[chLabel]->GetParameter(2) );
      
      p_deltaT_vs_energyRatio[chLabel] -> Fit(fitFun_energyRatioCorr[chLabel],"QRS");
    }
  }
  
  // -- third loop over events to apply amp walk corrections
  cout << "Third loop over events to apply amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry);

    if ( !acceptEvent_chRef[entry]) continue;

    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    float totRef = 0.5*((*tot)[channelIdx[chRef1]]+(*tot)[channelIdx[chRef2]]);
    if ( !useTimeAverage) {
      energyRef  = (*energy)[channelIdx[chRef1]];
      tRef       = (*time)[channelIdx[chRef1]];
      totRef       = (*tot)[channelIdx[chRef1]];
    }

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( !acceptEvent[entry][ch] ) continue;
	
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef;
	float energyRatioCorr = fitFun_energyRatioCorr[chLabel] -> Eval( energyRatio ) - fitFun_energyRatioCorr[chLabel] -> Eval( fitFun_energyRatio[chLabel]->GetParameter(1) ); 
	long long deltaT = (*time)[channelIdx[ch]] - tRef;

	//float totRatio = (*tot)[channelIdx[ch]]/totRef;
	float totRatio = ((*tot)[channelIdx[ch]]-offsetCh[std::make_pair(iBar , label.c_str())])/totRef ;
        //cout << " tot =  " << (*tot)[channelIdx[ch]] << " - offset " << offsetCh[std::make_pair(iBar , label.c_str())] << " =  " << (*tot)[channelIdx[ch]]-offsetCh[std::make_pair(iBar , label.c_str())];
	if ( fabs(deltaT)>10000) continue;   
	if ( fabs(deltaT-energyRatioCorr)>10000) continue;   


  //map<std::string,TH1F*>      h_totRatio;

	
	h_deltaT_energyRatioCorr[chLabel] -> Fill( deltaT - energyRatioCorr);
        p_deltaT_energyRatioCorr_vs_totRatio[chLabel] -> Fill ( totRatio ,deltaT - energyRatioCorr );
        h2_deltaT_energyRatioCorr_vs_totRatio[chLabel] -> Fill ( totRatio ,deltaT - energyRatioCorr );
	p_deltaT_energyRatioCorr_vs_t1fine[chLabel]  -> Fill( (*t1fine)[channelIdx[ch]] , deltaT - energyRatioCorr );
	h2_deltaT_energyRatioCorr_vs_t1fine[chLabel] -> Fill( (*t1fine)[channelIdx[ch]] , deltaT - energyRatioCorr );
      } // end loop L,R


      // tDiff L-R
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      if ( acceptEvent[entry][chL]  && acceptEvent[entry][chR]  && acceptEvent_chRef[entry])    {      
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
	float totRatio = ((*tot)[channelIdx[chL]]-offsetCh[std::make_pair(iBar , "L")])/((*tot)[channelIdx[chR]]-offsetCh[std::make_pair(iBar , "R")]);
	float energyRatioCorr = fitFun_energyRatioCorr_LR[iBar] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR[iBar] -> Eval( fitFun_energyRatio_LR[iBar]->GetParameter(1) ); 
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	
	if ( fabs(deltaT) < 10000 && fabs(deltaT - energyRatioCorr) < 10000){   // aggiungere quelli L-R
	  h_deltaT_LR_energyRatioCorr[iBar]-> Fill( deltaT - energyRatioCorr);
	  p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar]-> Fill( 0.5*( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]), deltaT - energyRatioCorr);
	}
      }


    }// end loop over bars
  }    

  // ---  amp walk corr : tot
  map<std::string,TF1*> fitFun_totRatio;
  map<std::string,TF1*> fitFun_energyRatioCorr_totRatioCorr;

  map<int,TF1*> fitFun_totRatio_LR;
  map<int,TF1*> fitFun_energyRatioCorr_totRatioCorr_LR;

  for(int iBar = 0; iBar < 16; ++iBar){  

    // L-R
    fitFun_totRatio_LR[iBar] = new TF1(Form("fitFun_totRatio_LR_bar%02d", iBar), "gaus", 0,10);  
    h_totRatio_LR[iBar] -> Fit(fitFun_totRatio_LR[iBar],"QR");
    
    fitFun_energyRatioCorr_totRatioCorr_LR[iBar] = new TF1(Form("fitFun_energyRatioCorr_totRatioCorr_LR_bar%02d", iBar), "pol3", 0,10);
    fitFun_energyRatioCorr_totRatioCorr_LR[iBar]->SetRange( fitFun_totRatio_LR[iBar]->GetParameter(1)-5*fitFun_totRatio_LR[iBar]->GetParameter(2), fitFun_totRatio_LR[iBar]->GetParameter(1)+5*fitFun_totRatio_LR[iBar]->GetParameter(2));
    p_deltaT_LR_energyRatioCorr_vs_totRatio[iBar] -> Fit(fitFun_energyRatioCorr_totRatioCorr_LR[iBar],"QRS"); 


    // - single channels 
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      int ch = chID[chLabel]; 

      fitFun_totRatio[chLabel] = new TF1(Form("fitFun_totRatio_%s",chLabel.c_str()), "gaus", 0,10);  
      h_totRatio[chLabel] -> Fit(fitFun_totRatio[chLabel],"QR");

      fitFun_energyRatioCorr_totRatioCorr[chLabel] = new TF1(Form("fitFun_energyRatioCorr_totRatioCorr_ch%02d",ch), "pol3", 0,10);
      //fitFun_energyRatioCorr[chLabel]->SetRange( fitFun_energyRatio[chLabel]->GetParameter(1) - 3*fitFun_energyRatio[chLabel]->GetParameter(2), fitFun_energyRatio[chLabel]->GetParameter(1) + 3*fitFun_energyRatio[chLabel]->GetParameter(2) );
      fitFun_energyRatioCorr_totRatioCorr[chLabel]->SetRange( fitFun_totRatio[chLabel]->GetParameter(1) - 5*fitFun_totRatio[chLabel]->GetParameter(2), fitFun_totRatio[chLabel]->GetParameter(1) + 5*fitFun_totRatio[chLabel]->GetParameter(2) );
       
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel] -> Fit(fitFun_energyRatioCorr_totRatioCorr[chLabel],"QRS");
    }
  }
 


  // -- fourth loop over events to apply amp walk corrections
  cout << "Fourth loop over events to apply phase corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data ->GetEntry(entry);

    if ( !acceptEvent_chRef[entry] ) continue;
    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    float  totRef = 0.5*((*tot)[channelIdx[chRef1]]+(*tot)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    if ( !useTimeAverage) {
      energyRef  = (*energy)[channelIdx[chRef1]];
      totRef  = (*tot)[channelIdx[chRef1]];
      tRef       = (*time)[channelIdx[chRef1]];
    }

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( !acceptEvent[entry][ch] ) continue;
	
	long long deltaT = (*time)[channelIdx[ch]] - tRef;
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef;
	//float totRatio = (*tot)[channelIdx[ch]]/totRef;
	float totRatio = ((*tot)[channelIdx[ch]]-offsetCh[std::make_pair(iBar , label.c_str())])/totRef ;
	float energyRatioCorr = fitFun_energyRatioCorr[chLabel] -> Eval( energyRatio ) - fitFun_energyRatioCorr[chLabel] -> Eval( fitFun_energyRatio[chLabel]->GetParameter(1) ); 
	int bin1  = p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->FindBin( (*t1fine)[channelIdx[ch]]) ; 
	int bin2 = p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->FindBin( p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->GetMean() );
	float phaseCorr = p_deltaT_energyRatioCorr_vs_t1fine[chLabel] -> GetBinContent(bin1) -  p_deltaT_energyRatioCorr_vs_t1fine[chLabel] -> GetBinContent(bin2);

	float energyRatioCorr_totRatioCorr = fitFun_energyRatioCorr_totRatioCorr[chLabel] -> Eval( totRatio ) - fitFun_energyRatioCorr_totRatioCorr[chLabel] -> Eval( fitFun_totRatio[chLabel]->GetParameter(1) ); 
	
      
	if ( fabs(deltaT)>10000) continue;   
	if ( fabs(deltaT-energyRatioCorr)>10000) continue;   
	if ( fabs(deltaT-energyRatioCorr- energyRatioCorr_totRatioCorr)>10000) continue;   
	
	h_deltaT_energyRatioCorr_phaseCorr[chLabel] -> Fill( deltaT - energyRatioCorr - phaseCorr);
	h_deltaT_energyRatioCorr_totRatioCorr[chLabel] -> Fill( deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr);

        p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] -> Fill( (*t1fine)[channelIdx[ch]] , deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr);
        h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] -> Fill( (*t1fine)[channelIdx[ch]], deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr);

	p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]  -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr - phaseCorr);
	h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel] -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr - phaseCorr);

      } // end loop L,R

      
      // tDiff L-R
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      if ( acceptEvent[entry][chL]  && acceptEvent[entry][chR]  && acceptEvent_chRef[entry])    {      
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
	//float totRatio = (*tot)[channelIdx[chL]]/(*tot)[channelIdx[chR]];
	float totRatio = ((*tot)[channelIdx[chL]]-offsetCh[std::make_pair(iBar , "L")])/((*tot)[channelIdx[chR]]-offsetCh[std::make_pair(iBar , "R")]);

	float energyRatioCorr = fitFun_energyRatioCorr_LR[iBar] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR[iBar] -> Eval( fitFun_energyRatio_LR[iBar]->GetParameter(1) ); 
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	int bin1  = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar]->FindBin( 0.5 * ( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]) ); 
	int bin2 = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar]->FindBin( p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar]->GetMean() );
	float phaseCorr = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar] -> GetBinContent(bin1) -  p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar] -> GetBinContent(bin2);


	float energyRatioCorr_totRatioCorr = fitFun_energyRatioCorr_totRatioCorr_LR[iBar] -> Eval(totRatio) - fitFun_energyRatioCorr_totRatioCorr_LR[iBar] -> Eval( fitFun_totRatio_LR[iBar]->GetParameter(1) ); 

	if ( fabs(deltaT)<10000 && fabs(deltaT-energyRatioCorr)<10000 && fabs(deltaT-energyRatioCorr- energyRatioCorr_totRatioCorr)<10000){
	  h_deltaT_LR_energyRatioCorr_phaseCorr[iBar] -> Fill( deltaT - energyRatioCorr - phaseCorr);
	  h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar] -> Fill( deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr);
          p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar] -> Fill(0.5*( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]) ,deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr);

	}
      }
      
    }    // end loop over bars
  }  
  


  cout << "Fifth loop over events to apply phase corrections : on top of tot" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data ->GetEntry(entry);

    if ( !acceptEvent_chRef[entry] ) continue;
    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    float  totRef = 0.5*((*tot)[channelIdx[chRef1]]+(*tot)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    if ( !useTimeAverage) {
      energyRef  = (*energy)[channelIdx[chRef1]];
      totRef  = (*tot)[channelIdx[chRef1]];
      tRef       = (*time)[channelIdx[chRef1]];
    }

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( !acceptEvent[entry][ch] ) continue;
	
	long long deltaT = (*time)[channelIdx[ch]] - tRef;
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef;
	//float totRatio = (*tot)[channelIdx[ch]]/totRef;
	float totRatio = ((*tot)[channelIdx[ch]]-offsetCh[std::make_pair(iBar , label.c_str())])/totRef ;
	float energyRatioCorr = fitFun_energyRatioCorr[chLabel] -> Eval( energyRatio ) - fitFun_energyRatioCorr[chLabel] -> Eval( fitFun_energyRatio[chLabel]->GetParameter(1) ); 
	int bin1 = p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->FindBin( (*t1fine)[channelIdx[ch]]) ; 
	int bin2 = p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->FindBin( p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->GetMean() );
	float phaseCorr = p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] -> GetBinContent(bin1) -  p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel] -> GetBinContent(bin2);

	float energyRatioCorr_totRatioCorr = fitFun_energyRatioCorr_totRatioCorr[chLabel] -> Eval( totRatio ) - fitFun_energyRatioCorr_totRatioCorr[chLabel] -> Eval( fitFun_totRatio[chLabel]->GetParameter(1) ); 
	
      
	if ( fabs(deltaT)>10000) continue;   
	if ( fabs(deltaT-energyRatioCorr)>10000) continue;   
	if ( fabs(deltaT-energyRatioCorr- energyRatioCorr_totRatioCorr)>10000) continue;   
	
	h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[chLabel] -> Fill( deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr -  phaseCorr);

	p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]  -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr - phaseCorr);
	h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel] -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr -  energyRatioCorr_totRatioCorr - phaseCorr);

      } // end loop L,R

      
      // tDiff L-R
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      if ( acceptEvent[entry][chL]  && acceptEvent[entry][chR]  && acceptEvent_chRef[entry])    {      
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
	//float totRatio = (*tot)[channelIdx[chL]]/(*tot)[channelIdx[chR]];
	float totRatio = ((*tot)[channelIdx[chL]]-offsetCh[std::make_pair(iBar , "L")])/((*tot)[channelIdx[chR]]-offsetCh[std::make_pair(iBar , "R")]);

	float energyRatioCorr = fitFun_energyRatioCorr_LR[iBar] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR[iBar] -> Eval( fitFun_energyRatio_LR[iBar]->GetParameter(1) ); 
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	int bin1 = p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar]->FindBin( 0.5 * ( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]) ); 
	int bin2 = p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar]->FindBin( p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar]->GetMean() );
	float phaseCorr = p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar] -> GetBinContent(bin1) -  p_deltaT_LR_energyRatioCorr_totRatioCorr_vs_t1fine[iBar] -> GetBinContent(bin2);


	float energyRatioCorr_totRatioCorr = fitFun_energyRatioCorr_totRatioCorr_LR[iBar] -> Eval(totRatio) - fitFun_energyRatioCorr_totRatioCorr_LR[iBar] -> Eval( fitFun_totRatio_LR[iBar]->GetParameter(1) ); 

	if ( fabs(deltaT)<10000 && fabs(deltaT-energyRatioCorr)<10000 && fabs(deltaT-energyRatioCorr- energyRatioCorr_totRatioCorr)<10000){
	  h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar] -> Fill( deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr -  phaseCorr);

	}
      }
      
    }    // end loop over bars
  }  
 

  // -- gaus fit deltaT for each channel
  map<std::string,TGraphErrors*> g_tRes;
  map<std::string,TGraphErrors*> g_tRes_energyRatioCorr;
  map<std::string,TGraphErrors*> g_tRes_energyRatioCorr_totRatioCorr;
  map<std::string,TGraphErrors*> g_tRes_energyRatioCorr_phaseCorr;
  map<std::string,TGraphErrors*> g_tRes_energyRatioCorr_totRatioCorr_phaseCorr;

  map<std::string,TF1*> fitGaus;
  map<std::string,TF1*> fitGaus_energyRatioCorr;
  map<std::string,TF1*> fitGaus_energyRatioCorr_totRatioCorr;
  map<std::string,TF1*> fitGaus_energyRatioCorr_phaseCorr;
  map<std::string,TF1*> fitGaus_energyRatioCorr_totRatioCorr_phaseCorr;

  TGraphErrors *g_tRes_LR = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_LR = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_phaseCorr_LR = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_LR = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR = new TGraphErrors();

  map<int,TF1*> fitGaus_LR;
  map<int,TF1*> fitGaus_energyRatioCorr_LR;
  map<int,TF1*> fitGaus_energyRatioCorr_phaseCorr_LR;
  map<int,TF1*> fitGaus_energyRatioCorr_totRatioCorr_LR;
  map<int,TF1*> fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR;
  

  // delta T L-R
  std::cout << "Time resolution from tDiff(L-R)"<<std::endl;
  for(int iBar = 0; iBar < 16; ++iBar){
    
    if ( h_deltaT_LR[iBar] -> GetEntries() == 0) continue;

    if (runs == "5352" && iBar == 13) continue; // bad energy in one channel of bar13, skip!

    // -- no corr                                                                                                                                                         
    fitGaus_LR[iBar] = new TF1(Form("fitGaus_LR_bar%02d", iBar), "gaus",-10000,10000);                                                                                  
    getTimeResolution(h_deltaT_LR[iBar], fitGaus_LR[iBar]);
    g_tRes_LR-> SetPoint( g_tRes_LR->GetN(), iBar, fitGaus_LR[iBar]->GetParameter(2));
    g_tRes_LR-> SetPointError( g_tRes_LR->GetN()-1, 0, fitGaus_LR[iBar]->GetParError(2));

    // -- energy corr
    fitGaus_energyRatioCorr_LR[iBar] = new TF1(Form("fitGaus_energyRatioCorr_LR_bar%02d", iBar), "gaus",-10000,10000);
    getTimeResolution(h_deltaT_LR_energyRatioCorr[iBar], fitGaus_energyRatioCorr_LR[iBar]); 
    g_tRes_energyRatioCorr_LR-> SetPoint( g_tRes_energyRatioCorr_LR->GetN(), iBar, fitGaus_energyRatioCorr_LR[iBar]->GetParameter(2));
    g_tRes_energyRatioCorr_LR-> SetPointError( g_tRes_energyRatioCorr_LR->GetN()-1, 0, fitGaus_energyRatioCorr_LR[iBar]->GetParError(2));

    // -- energy + phase corr
    fitGaus_energyRatioCorr_phaseCorr_LR[iBar] = new TF1(Form("fitGaus_energyRatioCorr_phaseCorr_LR_bar%02d", iBar), "gaus",-10000,10000);
    getTimeResolution(h_deltaT_LR_energyRatioCorr_phaseCorr[iBar], fitGaus_energyRatioCorr_phaseCorr_LR[iBar]); 
    g_tRes_energyRatioCorr_phaseCorr_LR-> SetPoint( g_tRes_energyRatioCorr_phaseCorr_LR->GetN(), iBar, fitGaus_energyRatioCorr_phaseCorr_LR[iBar]->GetParameter(2));
    g_tRes_energyRatioCorr_phaseCorr_LR-> SetPointError( g_tRes_energyRatioCorr_phaseCorr_LR->GetN()-1, 0, fitGaus_energyRatioCorr_phaseCorr_LR[iBar]->GetParError(2));
    

    // -- energy + tot corr
    fitGaus_energyRatioCorr_totRatioCorr_LR[iBar] = new TF1(Form("fitGaus_energyRatioCorr_totRatioCorr_LR_bar%02d", iBar), "gaus",-10000,10000);
    getTimeResolution(h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar], fitGaus_energyRatioCorr_totRatioCorr_LR[iBar]); 
    g_tRes_energyRatioCorr_totRatioCorr_LR-> SetPoint( g_tRes_energyRatioCorr_totRatioCorr_LR->GetN(), iBar, fitGaus_energyRatioCorr_totRatioCorr_LR[iBar]->GetParameter(2));
    g_tRes_energyRatioCorr_totRatioCorr_LR-> SetPointError( g_tRes_energyRatioCorr_totRatioCorr_LR->GetN()-1, 0, fitGaus_energyRatioCorr_totRatioCorr_LR[iBar]->GetParError(2));

     // -- energy + tot + phase  corr
    fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR[iBar] = new TF1(Form("fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR_bar%02d", iBar), "gaus",-10000,10000);
    getTimeResolution(h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar], fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR[iBar]); 
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR-> SetPoint( g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR->GetN(), iBar, fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR[iBar]->GetParameter(2));
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR-> SetPointError( g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR->GetN()-1, 0, fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_LR[iBar]->GetParError(2));
 
      
  }



  // deltaT ch - Ref
  std::cout << "time resolution from tDiff(ch-Ref)"<<std::endl;
  for (auto label : labelLR ){
    g_tRes[label] = new TGraphErrors();
    g_tRes_energyRatioCorr[label] = new TGraphErrors();
    g_tRes_energyRatioCorr_totRatioCorr[label] = new TGraphErrors();
    g_tRes_energyRatioCorr_phaseCorr[label] = new TGraphErrors();
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label] = new TGraphErrors();
    
    for(int iBar = 0; iBar < 16; ++iBar){        
      bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
      if (!barFound) continue;
 
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());

      if ( h_deltaT[chLabel] -> GetEntries() == 0) continue;

      if (runs == "5352" && chLabel == "bar13L") continue; // bad energy in one channel of bar13, skip!

      // -- no corr
      fitGaus[chLabel] = new TF1(Form("fitGaus_%s",chLabel.c_str()), "gaus",-10000,10000);
      getTimeResolution(h_deltaT[chLabel], fitGaus[chLabel]);
      g_tRes[label]-> SetPoint( g_tRes[label]->GetN(), iBar, fitGaus[chLabel]->GetParameter(2));
      g_tRes[label]-> SetPointError( g_tRes[label]->GetN()-1, 0, fitGaus[chLabel]->GetParError(2));
      
      // -- energy corr
      fitGaus_energyRatioCorr[chLabel] = new TF1(Form("fitGaus_energyRatioCorr_%s",chLabel.c_str()), "gaus",-10000,10000);
      getTimeResolution(h_deltaT_energyRatioCorr[chLabel], fitGaus_energyRatioCorr[chLabel]); 
      g_tRes_energyRatioCorr[label]-> SetPoint( g_tRes_energyRatioCorr[label]->GetN(), iBar, fitGaus_energyRatioCorr[chLabel]->GetParameter(2));
      g_tRes_energyRatioCorr[label]-> SetPointError( g_tRes_energyRatioCorr[label]->GetN()-1, 0, fitGaus_energyRatioCorr[chLabel]->GetParError(2));

       // -- energy + tot corr
      fitGaus_energyRatioCorr_totRatioCorr[chLabel] = new TF1(Form("fitGaus_energyRatioCorr_totRatioCorr_%s",chLabel.c_str()), "gaus",-10000,10000);
      getTimeResolution(h_deltaT_energyRatioCorr_totRatioCorr[chLabel], fitGaus_energyRatioCorr_totRatioCorr[chLabel]); 
      g_tRes_energyRatioCorr_totRatioCorr[label]-> SetPoint( g_tRes_energyRatioCorr_totRatioCorr[label]->GetN(), iBar, fitGaus_energyRatioCorr_totRatioCorr[chLabel]->GetParameter(2));
      g_tRes_energyRatioCorr_totRatioCorr[label]-> SetPointError( g_tRes_energyRatioCorr_totRatioCorr[label]->GetN()-1, 0, fitGaus_energyRatioCorr_totRatioCorr[chLabel]->GetParError(2));
      
      // -- energy + phase corr
      fitGaus_energyRatioCorr_phaseCorr[chLabel] = new TF1(Form("fitGaus_energyRatioCorr_phaseCorr_%s",chLabel.c_str()), "gaus",-10000,10000);
      getTimeResolution(h_deltaT_energyRatioCorr_phaseCorr[chLabel], fitGaus_energyRatioCorr_phaseCorr[chLabel]);
      g_tRes_energyRatioCorr_phaseCorr[label]-> SetPoint( g_tRes_energyRatioCorr_phaseCorr[label]->GetN(), iBar, fitGaus_energyRatioCorr_phaseCorr[chLabel]->GetParameter(2));
      g_tRes_energyRatioCorr_phaseCorr[label]-> SetPointError( g_tRes_energyRatioCorr_phaseCorr[label]->GetN()-1, 0, fitGaus_energyRatioCorr_phaseCorr[chLabel]->GetParError(2));

       // -- energy + tot + phase corr
      fitGaus_energyRatioCorr_totRatioCorr_phaseCorr[chLabel] = new TF1(Form("fitGaus_energyRatioCorr_totRatioCorr_phaseCorr_%s",chLabel.c_str()), "gaus",-10000,10000);
      getTimeResolution(h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[chLabel], fitGaus_energyRatioCorr_totRatioCorr_phaseCorr[chLabel]); 
      g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]-> SetPoint( g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->GetN(), iBar, fitGaus_energyRatioCorr_totRatioCorr_phaseCorr[chLabel]->GetParameter(2));
      g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]-> SetPointError( g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->GetN()-1, 0, fitGaus_energyRatioCorr_totRatioCorr_phaseCorr[chLabel]->GetParError(2));
 
  
    }
  }


  
  std::cout << "Triangulation ..." << std::endl;
  std::cout << " >> Energy  ..." << std::endl;
   // triangulation energy 
  TGraphErrors *g_tRes_energyRatioCorr_L = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_R = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_chRef = new TGraphErrors();
  triangulation(barList, 
                h_deltaT_LR_energyRatioCorr, 
                h_deltaT_energyRatioCorr, 
                g_tRes_energyRatioCorr_LR ,
                g_tRes_energyRatioCorr, 
                g_tRes_energyRatioCorr_L,
                g_tRes_energyRatioCorr_R, 
                g_tRes_energyRatioCorr_chRef);

  
  // triangulation energy + phase
  std::cout << " >> Energy + phase ..." << std::endl;
  TGraphErrors *g_tRes_energyRatioCorr_phaseCorr_L = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_phaseCorr_R = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_phaseCorr_chRef = new TGraphErrors();
  triangulation(barList, 
                h_deltaT_LR_energyRatioCorr_phaseCorr, 
                h_deltaT_energyRatioCorr_phaseCorr,
                g_tRes_energyRatioCorr_phaseCorr_LR ,
                g_tRes_energyRatioCorr_phaseCorr, 
                g_tRes_energyRatioCorr_phaseCorr_L, 
                g_tRes_energyRatioCorr_phaseCorr_R, 
                g_tRes_energyRatioCorr_phaseCorr_chRef);

 

  // energy + tot corr
  std::cout << " >> Energy + tot ..." << std::endl;
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_L = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_R = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_chRef = new TGraphErrors();
  triangulation(barList,
                h_deltaT_LR_energyRatioCorr_totRatioCorr, 
                h_deltaT_energyRatioCorr_totRatioCorr, 
                g_tRes_energyRatioCorr_totRatioCorr_LR , 
                g_tRes_energyRatioCorr_totRatioCorr, 
                g_tRes_energyRatioCorr_totRatioCorr_L, 
                g_tRes_energyRatioCorr_totRatioCorr_R, 
                g_tRes_energyRatioCorr_totRatioCorr_chRef);

 // triangulation energy + tot + phase
  std::cout << " >> Energy + tot + phase ..." << std::endl;
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_L = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_R = new TGraphErrors();
  TGraphErrors *g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_chRef = new TGraphErrors();
  triangulation(barList, 
                h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr, 
                h_deltaT_energyRatioCorr_totRatioCorr_phaseCorr, 
                g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_LR,
                g_tRes_energyRatioCorr_totRatioCorr_phaseCorr, 
                g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_L, 
                g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_R,
                g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_chRef);

 


  // ======  save histograms in a file
//  string foutN = Form("%s/analysisSingleChannel_runs%s",outDir.c_str(),runs.c_str());
//  if (useTimeAverage) foutN = Form("%s/analysisSingleChannel_runs%s_timeAverage",outDir.c_str(),runs.c_str());
//  if (offsetFileName != "NOOFFSET") foutN = Form("%s_offset",foutN.c_str());
//  string foutName = Form ("%s.root", foutN.c_str());
//  TFile *fout = new TFile(foutName.c_str(),"recreate");
//
//  for (auto label : labelLR ){
//    for(int iBar = 0; iBar < 16; ++iBar){        
//      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
//      h_energyRatio[chLabel]->Write();
//      h2_energyRatio_vs_totRatio[chLabel]->Write();
//      h_energy[chLabel]->Write();
//      h_deltaT[chLabel]->Write();
//      h_deltaT_energyRatioCorr[chLabel]->Write();
//      p_deltaT_vs_energyRatio[chLabel]->Write();
//    }
//    
//    g_tRes[label]->Write(Form("g_tRes_%s",label.c_str()));
//    g_tRes_energyRatioCorr[label]->Write(Form("g_tRes_energyRatioCorr_%s",label.c_str()));
//    g_tRes_energyRatioCorr_phaseCorr[label]->Write( Form("g_tRes_energyRatioCorr_phaseCorr_%s",label.c_str()));
//    g_tRes_energyRatioCorr_totRatioCorr[label]->Write( Form("g_tRes_energyRatioCorr_totRatioCorr_%s",label.c_str()));
//  }
//  
//  fout->Close();


  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  // ======== PLOT --> INSERT IN CONFIG FILE 
  std::string plotDir = (Form("%s/analysisSingleChannel/%s/",mainPlotDir.c_str(),runs.c_str() ));
  if (useTimeAverage) plotDir = Form("%s/analysisSingleChannel/%s_timeAverage/",mainPlotDir.c_str(),runs.c_str() );
  if (offsetFileName != "NOOFFSET") plotDir = Form("%s/offset/",plotDir.c_str());
  system(Form("mkdir -p %s",plotDir.c_str()));  

  cout<< "Printing plots ..."<<endl;

  TCanvas *c;

  c = new TCanvas("c","c", 800, 600);
  h_nActiveChannels0 -> GetXaxis()->SetTitle("nActiveChannels");
  h_nActiveChannels0->Draw(); 
  TLine* line = new TLine(maxActiveChannels,0.,maxActiveChannels+0.5,h_nActiveChannels0 ->GetMaximum());
  line -> SetLineWidth(1);
  line -> SetLineStyle(7);
  line -> Draw("same");     
  c->Print(Form("%s/c_nActiveChannels_0.png",plotDir.c_str()));
  c->Print(Form("%s/c_nActiveChannels_0.pdf",plotDir.c_str()));
  delete c;

  c = new TCanvas("c","c", 800, 600);
  h_nActiveChannels1 -> GetXaxis()->SetTitle("nActiveChannels");
  h_nActiveChannels1->Draw();
  line -> Draw("same");           
  c->Print(Form("%s/c_nActiveChannels_1.png",plotDir.c_str()));
  c->Print(Form("%s/c_nActiveChannels_1.pdf",plotDir.c_str()));
  delete c;


  c = new TCanvas("c","c", 800, 600);
  //c->SetLogy();
  h_energy_chRef-> GetXaxis()->SetTitle("energy [ADC]");
  h_energy_chRef->Draw();      
  c->Print(Form("%s/c_energy_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_energy_chRef.pdf",plotDir.c_str()));
  delete c;


  c = new TCanvas("c","c", 800, 600);
  //c->SetLogy();
  h_tot_chRef-> GetXaxis()->SetTitle("tot [ps]");
  h_tot_chRef->Draw();      
  c->Print(Form("%s/c_tot_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_tot_chRef.pdf",plotDir.c_str()));
  delete c;



  c = new TCanvas("c","c", 800, 600);
  h_deltaT_LR_chRef-> GetXaxis()-> SetTitle("#Deltat [ps] ");
  TF1 *ff = new TF1("ff","gaus",-10000,10000);
  h_deltaT_LR_chRef-> Fit(ff,"QRS");
  ff->SetRange( ff->GetParameter(1)-2*ff->GetParameter(2), ff->GetParameter(1)+2*ff->GetParameter(2));
  h_deltaT_LR_chRef-> Fit(ff,"QRS");
  h_deltaT_LR_chRef-> SetMarkerStyle(20);
  h_deltaT_LR_chRef-> GetXaxis()-> SetRangeUser( ff->GetParameter(1)-7*ff->GetParameter(2),ff->GetParameter(1)+7*ff->GetParameter(2) );
  h_deltaT_LR_chRef-> Draw("e");
  c->Print(Form("%s/c_deltaT_LR_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_deltaT_LR_chRef.pdf",plotDir.c_str()));
  delete c;


  // Delta T Corr Histos for tDiff L-R
  for(int iBar = 0; iBar < 16; ++iBar){
    c = new TCanvas("c","c", 800, 600);
    h_deltaT_LR_energyRatioCorr[iBar]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
    h_deltaT_LR_energyRatioCorr[iBar]-> SetMarkerStyle(20);
    h_deltaT_LR_energyRatioCorr[iBar]-> Draw("e");
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;

    c = new TCanvas("c","c", 800, 600);
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar]-> SetMarkerStyle(20);
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar]-> Draw("e");
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;

    c = new TCanvas("c","c", 800, 600);
    h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
    h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar]-> SetMarkerStyle(20);
    h_deltaT_LR_energyRatioCorr_totRatioCorr[iBar]-> Draw("e");
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_totRatioCorr_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_totRatioCorr_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;

    c = new TCanvas("c","c", 800, 600);
    h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
    h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar]-> SetMarkerStyle(20);
    h_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr[iBar]-> Draw("e");
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_totRatioCorr_phaseCorr_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;




  }

  
  
  for (auto label : labelLR ){
    for(int iBar = 0; iBar < 16; ++iBar){        
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      
      // -- energy
      c = new TCanvas("c","c", 800, 600);
      //c->SetLogy();
      h_energy[chLabel]-> GetXaxis()->SetTitle("energy [ADC]");
      h_energy[chLabel]->Draw();
      TLine* line = new TLine(energyMin[chLabel],0.,energyMin[chLabel],h_energy[chLabel]->GetMaximum());
      line -> SetLineWidth(1);
      line -> SetLineStyle(7);
      line -> Draw("same");
      TLine* line2 = new TLine(maxEnergy,0.,maxEnergy,h_energy[chLabel]->GetMaximum());
      line2 -> SetLineWidth(1);
      line2 -> SetLineStyle(7);
      line2 -> Draw("same");
      c->Print(Form("%s/c_energy_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_energy_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- deltaT
      c = new TCanvas("c","c", 800, 600);
      h_deltaT[chLabel]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT[chLabel]->SetMarkerStyle(20);
      h_deltaT[chLabel]->Draw("e");
      c->Print(Form("%s/c_deltaT_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- tot Ratio
      c = new TCanvas("c","c", 800, 600);
      h_totRatio[chLabel]-> GetXaxis()->SetTitle("tot_{ch}/tot_{chRef}");
      h_totRatio[chLabel]->Draw();
      c->Print(Form("%s/c_totRatio_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_totRatio_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

      // -- energy Ratio
      c = new TCanvas("c","c", 800, 600);
      h_energyRatio[chLabel]-> GetXaxis()->SetTitle("E_{ch}/E_{chRef}");
      h_energyRatio[chLabel]->Draw();
      c->Print(Form("%s/c_energyRatio_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_energyRatio_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
 

      
      // -- deltaT vs energyRatio 
      c = new TCanvas("c","c", 800, 600);
      gStyle->SetOptFit(0);  
      p_deltaT_vs_energyRatio[chLabel]-> SetMarkerStyle(20);
      p_deltaT_vs_energyRatio[chLabel]-> SetMarkerSize(1);
      p_deltaT_vs_energyRatio[chLabel]-> GetXaxis()-> SetRangeUser( p_deltaT_vs_energyRatio[chLabel]->GetMean(1) - 3*p_deltaT_vs_energyRatio[chLabel]->GetRMS(1), p_deltaT_vs_energyRatio[chLabel]->GetMean(1) + 3*p_deltaT_vs_energyRatio[chLabel]->GetRMS(1));
      p_deltaT_vs_energyRatio[chLabel]-> GetYaxis()-> SetRangeUser(p_deltaT_vs_energyRatio[chLabel]->GetMean(2) - 3*p_deltaT_vs_energyRatio[chLabel]->GetRMS(2), p_deltaT_vs_energyRatio[chLabel]->GetMean(2) + 3*p_deltaT_vs_energyRatio[chLabel]->GetRMS(2));
      p_deltaT_vs_energyRatio[chLabel]-> GetYaxis()->SetTitle("E_{ch}/E_{chRef}");
      p_deltaT_vs_energyRatio[chLabel]-> GetYaxis()->SetTitle("#DeltaT [ps]");
      p_deltaT_vs_energyRatio[chLabel]->Draw();
      h2_deltaT_vs_energyRatio[chLabel]->Draw("colz same");
      p_deltaT_vs_energyRatio[chLabel]->Draw("same");
      c->Print(Form("%s/c_deltaT_vs_energyRatio_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_vs_energyRatio_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      gStyle->SetOptFit(1111);  

      // -- deltaT Energy Corr vs totRatio 
      c = new TCanvas("c","c", 800, 600);
      gStyle->SetOptFit(0);  
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> SetMarkerSize(1);
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> GetXaxis()-> SetRangeUser(p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetMean(1) - 3*p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetRMS(1), p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetMean(1) + 3*p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetRMS(1));
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> GetYaxis()-> SetRangeUser(p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetMean(2) - 3*p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetRMS(2), p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetMean(2) + 3*p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->GetRMS(2));
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> GetYaxis()->SetTitle("tot_{ch}/tot_{chRef}");
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]-> GetYaxis()->SetTitle("#DeltaT [ps]");
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->Draw();
      h2_deltaT_energyRatioCorr_vs_totRatio[chLabel]->Draw("colz same");
      p_deltaT_energyRatioCorr_vs_totRatio[chLabel]->Draw("same");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_totRatio%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_totRatio%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      gStyle->SetOptFit(1111);  
 

      
      // -- deltaT corr
      c = new TCanvas("c","c", 800, 600);
      h_deltaT_energyRatioCorr[chLabel]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT_energyRatioCorr[chLabel]-> SetMarkerStyle(20);
      h_deltaT_energyRatioCorr[chLabel]-> Draw("e");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

      // -- deltaT energy + tot corr
      c = new TCanvas("c","c", 800, 600);
      h_deltaT_energyRatioCorr_totRatioCorr[chLabel]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT_energyRatioCorr_totRatioCorr[chLabel]-> SetMarkerStyle(20);
      h_deltaT_energyRatioCorr_totRatioCorr[chLabel]-> Draw("e");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

      
      // -- deltaT energy+phase corr
      c = new TCanvas("c","c", 800, 600);
      h_deltaT_energyRatioCorr_phaseCorr[chLabel]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT_energyRatioCorr_phaseCorr[chLabel]-> SetMarkerStyle(20);
      h_deltaT_energyRatioCorr_phaseCorr[chLabel]-> Draw("e");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

      // -- deltaT corr vs t1fine
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]-> GetXaxis()->SetTitle("t1fine");
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]-> GetYaxis()->SetRangeUser( h2_deltaT_energyRatioCorr_vs_t1fine[chLabel]->GetMean(2) - 300, h2_deltaT_energyRatioCorr_vs_t1fine[chLabel]->GetMean(2) + 300);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->Draw();
      h2_deltaT_energyRatioCorr_vs_t1fine[chLabel]->Draw("colz same");
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_t1fine_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_t1fine_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- deltaT corr vs t1fine
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]-> GetXaxis()->SetTitle("t1fine");
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]-> GetYaxis()->SetRangeUser( h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->GetMean(2) - 300, h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->GetMean(2) + 300);
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->Draw();
      h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->Draw("colz same");
      p_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine[chLabel]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_vs_t1fine_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
 
      
      // -- deltaT corr vs t1fine ref channel
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]-> GetXaxis()->SetTitle("t1fine_{chRef}");
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]-> GetYaxis()->SetRangeUser( h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->GetMean(2) - 300, h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->GetMean(2) + 300);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw();
      h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw("colz same");
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

       // -- deltaT corr vs t1fine ref channel
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]-> GetXaxis()->SetTitle("t1fine_{chRef}");
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]-> GetYaxis()->SetRangeUser( h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->GetMean(2) - 300, h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->GetMean(2) + 300);
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw();
      h2_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw("colz same");
      p_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef[chLabel]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_vs_t1fineRef_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
     

    }
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  
  // -- no corr
  //float ymin = 150;
  //float ymax = 350;
  if (step1 == 1.50) ymax = 240;

  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  TH2F* hdummy1 = new TH2F("hdummy1","",16,-0.5,15.5,100,ymin,ymax);
  hdummy1->GetXaxis()->SetTitle("bar");
  hdummy1->GetYaxis()->SetTitle("#sigma(t_{ch} - t_{ref}) [ps]");
  hdummy1->Draw();
  for (auto label : labelLR ){   
    g_tRes[label]->SetMarkerStyle(20);
    g_tRes[label]->SetName(Form("g_tRes_%s", label.c_str()));
    if (label == "R") g_tRes[label]->SetMarkerStyle(24);
    g_tRes[label]->SetMarkerSize(1);
    g_tRes[label]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_vs_bar.png",plotDir.c_str()));                                                                                                     
  c->Print(Form("%s/c_tRes_vs_bar.pdf",plotDir.c_str()));                                                                                                     
  delete c;
  
  
  // -- energy corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  TH2F* hdummy2 = new TH2F("hdummy2","",16,-0.5,15.5,100,ymin,ymax);
  hdummy2->GetXaxis()->SetTitle("bar");
  hdummy2->GetYaxis()->SetTitle("#sigma(t_{ch} - t_{ref}) [ps]");
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr[label]->SetName(Form("g_tRes_energyRatioCorr_%s", label.c_str()));
    g_tRes_energyRatioCorr[label]->SetMarkerColor(613); // kMagenta - 3
    g_tRes_energyRatioCorr[label]->SetLineColor(613);
    g_tRes_energyRatioCorr[label]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr[label]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr[label]->SetMarkerSize(1);
    g_tRes_energyRatioCorr[label]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_vs_bar.pdf",plotDir.c_str()));
  c->SaveAs(Form("%s/c_tRes_energyRatioCorr_vs_bar.root",plotDir.c_str()));
  delete c;


  // -- energy + phase corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_phaseCorr[label]->SetName(Form("g_tRes_energyRatioCorr_phaseCorr_%s", label.c_str()));
    g_tRes_energyRatioCorr_phaseCorr[label]->SetMarkerColor(867); // kAzure + 7
    g_tRes_energyRatioCorr_phaseCorr[label]->SetLineColor(867);
    g_tRes_energyRatioCorr_phaseCorr[label]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_phaseCorr[label]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_phaseCorr[label]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_phaseCorr[label]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_bar.pdf",plotDir.c_str()));
  c->SaveAs(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_bar.root",plotDir.c_str()));
  delete c;





   // --tot corr +  energy corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_totRatioCorr[label]->SetName(Form("g_tRes_energyRatioCorr_totRatioCorr_%s", label.c_str()));
    g_tRes_energyRatioCorr_totRatioCorr[label]->SetMarkerColor(887); // kViolet + 7
    g_tRes_energyRatioCorr_totRatioCorr[label]->SetLineColor(887);
    g_tRes_energyRatioCorr_totRatioCorr[label]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_totRatioCorr[label]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_totRatioCorr[label]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_totRatioCorr[label]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_vs_bar.pdf",plotDir.c_str()));
  c->SaveAs(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_vs_bar.root",plotDir.c_str()));
  delete c;

  // -- energy + tot + phase corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetName(Form("g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_%s", label.c_str()));
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetMarkerColor(418); // kGreen + 2
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetLineColor(418);
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr[label]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_phaseCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_phaseCorr_vs_bar.pdf",plotDir.c_str()));
  c->SaveAs(Form("%s/c_tRes_energyRatioCorr_totRatioCorr_phaseCorr_vs_bar.root",plotDir.c_str()));
  delete c;


//   // --tot corr +  energy corr
//  c = new TCanvas("c","c", 800, 600);
//  c->SetGridx();
//  c->SetGridy();
//  hdummy2->Draw();
//  for (auto label : labelLR ){   
//    g_tRes_energyRatioCorr[label]->Draw("psame") ;
//    g_tRes_energyRatioCorr_totRatioCorr[label]->Draw("psame") ;
//  }
//
//  TLegend *lege = new TLegend(0.20, 0.78, 0.65, 0.88);
//  lege->SetBorderSize(0);
//  lege->SetFillStyle(0);
//  lege->AddEntry(g_tRes_energyRatioCorr["L"], "energy ratio corr", "P");
//  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr["L"], "energy & tot ratio corr", "P");
//  lege->Draw("same");
//
//
//
//  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr_vs_bar.png",plotDir.c_str()));
//  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr_vs_bar.pdf",plotDir.c_str()));
//  delete c;
 
  
    
  //-- tRes of L, R, chRef after triangulation
  //
 
  //std::cout << " >>>>>> Energy  ..." << std::endl;
  plotTimeResolution(g_tRes_energyRatioCorr_L,
                     g_tRes_energyRatioCorr_R, 
                     g_tRes_energyRatioCorr_chRef,
                     ymin, ymax, useTimeAverage, 613, plotDir,   
                     "energyRatioCorr" );

  //std::cout << " >>>>>> Energy +  phase ..." << std::endl;
  plotTimeResolution(g_tRes_energyRatioCorr_phaseCorr_L,
                     g_tRes_energyRatioCorr_phaseCorr_R, 
                     g_tRes_energyRatioCorr_phaseCorr_chRef, 
                     ymin, ymax, useTimeAverage, 867,  plotDir,   
                     "energyRatioCorr_phaseCorr" );


  //std::cout << " >>>>>> Energy + tot ..." << std::endl;
  plotTimeResolution(g_tRes_energyRatioCorr_totRatioCorr_L,
                     g_tRes_energyRatioCorr_totRatioCorr_R, 
                     g_tRes_energyRatioCorr_totRatioCorr_chRef, 
                     ymin, ymax, useTimeAverage, 887,  plotDir,   
                     "energyRatioCorr_totRatioCorr" );

  //std::cout << " >>>>>> Energy + tot + phase ..." << std::endl;
  plotTimeResolution(g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_L,
                     g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_R, 
                     g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_chRef, 
                     ymin, ymax, useTimeAverage, 418,  plotDir,   
                     "energyRatioCorr_totRatioCorr_phaseCorr" );

   // --tot corr +  energy corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  hdummy2->GetYaxis()->SetTitle("#sigma(t_{ch}) [ps]");
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_L->Draw("psame") ;
    g_tRes_energyRatioCorr_R->Draw("psame") ;
    g_tRes_energyRatioCorr_totRatioCorr_L->Draw("psame") ;
    g_tRes_energyRatioCorr_totRatioCorr_R->Draw("psame") ;
  }
  g_tRes_energyRatioCorr_chRef-> SetLineColor(921);  // kGray +1
  g_tRes_energyRatioCorr_chRef-> SetMarkerColor(921);
  g_tRes_energyRatioCorr_chRef->Draw("psame") ;
  g_tRes_energyRatioCorr_totRatioCorr_chRef->Draw("psame") ;

  TLegend *lege = new TLegend(0.20, 0.7, 0.8, 0.88);
  lege->SetBorderSize(0);
  lege->SetFillStyle(0);
  lege->SetNColumns(2);
  lege->AddEntry(g_tRes_energyRatioCorr_L, "E ratio corr, L", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_L, "E & tot ratio corr L", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_R, "E ratio corr, R", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_R, "E & tot ratio corr R", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_chRef, "E ratio corr, chRef", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_chRef, "E & tot ratio corr chRef", "P");
  lege->Draw("same");



  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr.pdf",plotDir.c_str()));
  delete c;
  delete lege; 

   // --tot corr +  energy corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_phaseCorr_L->Draw("psame") ;
    g_tRes_energyRatioCorr_phaseCorr_R->Draw("psame") ;
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_L->Draw("psame") ;
    g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_R->Draw("psame") ;
  }

  g_tRes_energyRatioCorr_phaseCorr_chRef->Draw("psame") ;
  g_tRes_energyRatioCorr_phaseCorr_chRef->SetLineColor(921) ;
  g_tRes_energyRatioCorr_phaseCorr_chRef->SetMarkerColor(921) ;
  g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_chRef->Draw("psame") ;
  lege = new TLegend(0.20, 0.7, 0.8, 0.88);
  lege->SetBorderSize(0);
  lege->SetFillStyle(0);
  lege->SetNColumns(2);
  lege->AddEntry(g_tRes_energyRatioCorr_phaseCorr_L, "E ratio + phase corr, L", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_L, "E & tot ratio + phase corr L", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_phaseCorr_R, "E ratio + phase corr, R", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_R, "E & tot ratio + phase corr R", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_phaseCorr_chRef, "E ratio + phase corr, chRef", "P");
  lege->AddEntry(g_tRes_energyRatioCorr_totRatioCorr_phaseCorr_chRef, "E & tot ratio + phase corr chRef", "P");
  lege->Draw("same");



  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr_phaseCorr.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_CFR_energyRatioCorr_totRatioCorr_phaseCorr.pdf",plotDir.c_str()));
  delete c;
  delete lege;  

  
}
 
