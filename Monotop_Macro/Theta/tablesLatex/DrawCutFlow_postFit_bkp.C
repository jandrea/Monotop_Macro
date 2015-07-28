//******************************************************
// Compute CutFlowTables from Proof output
//******************************************************

#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <fstream>
#include <vector>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <math.h>

#include <iostream>
#include <map>

using namespace std;

void CutFlow(TString signal_type){


  vector<TString> channels;
  channels.push_back("mujets");

  vector<TString> regions;
  regions.push_back("mWTmujetsSelectedSignalregion");
  regions.push_back("mWTmujetsWregionHighpt");
  regions.push_back("mWTmujetsttbarregionHighpt");

  vector<TString > samplelist_MC;
  samplelist_MC.push_back("TTMSDecays");
  samplelist_MC.push_back("WExclb"    );
  samplelist_MC.push_back("WExclc"    );
  samplelist_MC.push_back("WExcll"    );
  samplelist_MC.push_back("DY"        );
  samplelist_MC.push_back("VV"        );
  samplelist_MC.push_back("QCD"       );
  samplelist_MC.push_back("SingleTop" );
  samplelist_MC.push_back("SingleTopW");



  vector<TString > samplelist_signal;
  if(signal_type == "FCNC"){
/*    samplelist_signal.push_back("S4Inv100");
    samplelist_signal.push_back("S4Inv200");
    samplelist_signal.push_back("S4Inv300");
    samplelist_signal.push_back("S4Inv400");
*/    samplelist_signal.push_back("S4Inv500");
/*    samplelist_signal.push_back("S4Inv600");
    samplelist_signal.push_back("S4Inv700");
    samplelist_signal.push_back("S4Inv800");
    samplelist_signal.push_back("S4Inv900");
     samplelist_signal.push_back("S4Inv1000");
*/  }
  if(signal_type == "RES"){
    //samplelist_signal.push_back("S1Res300Inv100");
    samplelist_signal.push_back("S1Res500Inv100");
    samplelist_signal.push_back("S1Res700Inv100");
    samplelist_signal.push_back("S1Res900Inv100");
    //samplelist_signal.push_back("S1Res1300Inv100");
    //samplelist_signal.push_back("S1Res1500Inv100");
  }

  vector<TString > samplelist_data;
  samplelist_data.push_back("DATA");

  map<TString, TH1D* > histo_map;

  //*********************************
  //get histograms
  //*********************************
  TFile *f_data  = new TFile("../thetaInOut/inputTheta_merged_AllRegions.root");
  f_data->cd();

  for(unsigned int ireg=0; ireg<regions.size(); ireg++)
  {
    for(unsigned int idata=0; idata<samplelist_data.size(); idata++)
    {
        histo_map[ regions[ireg]+"__"+samplelist_data[idata]] = (TH1D*)f_data->Get( (regions[ireg]+"__"+samplelist_data[idata]).Data() )->Clone();
        if(histo_map[ regions[ireg]+"__"+samplelist_data[idata]] == 0) cout << " non existing histogram " << regions[ireg]+"__"+samplelist_data[idata] << endl;
    }
  }

  TFile *f_mc  = new TFile("../thetaInOut/outputTheta_merged_AllRegions.root");
  f_mc->cd();

  for(unsigned int ireg=0; ireg<regions.size(); ireg++)
  {
    for(unsigned int iMC=0; iMC<samplelist_MC.size(); iMC++)
    {
        if(samplelist_MC[iMC] == "QCD" && regions[ireg] == "mWTmujetsttbarregionHighpt") continue;
        histo_map[ regions[ireg]+"__"+samplelist_MC[iMC]] = (TH1D*)f_mc->Get( (regions[ireg]+"__"+samplelist_MC[iMC]).Data() )->Clone();
        if(histo_map[ regions[ireg]+"__"+samplelist_MC[iMC]] == 0) cout << "non existing histogram " << regions[ireg]+"__"+samplelist_MC[iMC] << endl;
    }
  }


  TFile *f_signal  = new TFile("../thetaInOut/inputTheta_merged_AllSignals_AllRegions.root");
  f_signal->cd();

  for(unsigned int ireg=0; ireg<regions.size(); ireg++)
  {
    for(unsigned int iMC=0; iMC<samplelist_signal.size(); iMC++)
    {
        histo_map[ regions[ireg]+"__"+samplelist_signal[iMC]] = (TH1D*)f_signal->Get( (regions[ireg]+"__"+samplelist_signal[iMC]).Data() )->Clone();
    }
  }


  string ofilenametex = "CutFlow_postFit.tex";
  ofstream ofile(ofilenametex.c_str());

  ofile<<"\\documentclass[10pt]{article}"<<endl;
  ofile<<"\\usepackage{amsthm}"<<endl;
  ofile<<"\\usepackage{color}"<<endl;
  //ofile<<"\\documentclass[amsmath,amssymb]{revtex4}"<<endl;

  ofile<<"\\begin{document}"<<endl;

  ofile.setf(ios::fixed);
  ofile.precision(1);
  vector<string> CutName;
  //CutName.push_back("Presel.");
  CutName.push_back("TTbar-enriched CR");
  CutName.push_back("W+Jets-enriched CR");
  CutName.push_back("Signal Region");
  //CutName.push_back("ATLAS FCNC SR");
  //CutName.push_back("ATLAS RES SR");

  //ofile.precision(3);
  //ofile << "\\clearpage" << endl;
  //ofile << "\\begin{landscape}" << endl;
  ofile << "\\begin{table}[p]" << endl;

  ofile << "\\begin{tabular}{|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << "\\hline" << endl;

  for(unsigned int isel = 0; isel <CutName.size(); isel++)
  {
     ofile << "  & "  << CutName[isel] ;
  }

  for(unsigned int iback = 0; iback < samplelist_MC.size(); iback++)
  {

    ofile << "\\hline " << endl;
        double total_background_0 = 0;
        double total_background_error_0 = 0;
        double total_background_1 = 0;
        double total_background_error_1 = 0;
        double total_background_2 = 0;
        double total_background_error_2 = 0;
        double total_SingleTop_0 = 0;
        double total_SingleTop_error_0 = 0;
        double total_SingleTop_1 = 0;
        double total_SingleTop_error_1 = 0;
        double total_SingleTop_2 = 0;
        double total_SingleTop_error_2 = 0;


    ofile << samplelist_MC[iback] <<" & ";
    for(unsigned int ireg = 0; ireg < regions.size(); ireg++)
    {

        //if(samplelist_MC[iback] == "QCD" && regions[ireg] == "mWTmujetsttbarregionHighpt") continue;
        TH1D* tmp = histo_map[ regions[ireg]+"__"+samplelist_MC[iback]];
        double err = 0;
        if(iback < 7)
        {
            cout << "sample = " << samplelist_MC[iback] << " | and Yield =" << tmp->Integral() << endl;
            ofile << tmp->IntegralAndError(1, tmp->GetNbinsX(), err) << "$\\pm$" << err << " & " ;
        }
        else if(iback == 7 || iback == 8)
        {
            total_SingleTop       += tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
            total_SingleTop_error += pow(err,2);
            if(iback == 8) ofile << total_SingleTop << "$\\pm$" << pow(total_SingleTop_error, 0.5) << " & " ;
        }

	    if(ireg == 0)
        {
            total_background_0        +=  tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
	        total_background_error_0  +=  pow(err, 2);
        }
        else if(ireg == 1)
        {
            total_background_1        +=  tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
	        total_background_error_1  +=  pow(err, 2);
        }
        else if(ireg == 2)
        {
            total_background_2        +=  tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
	        total_background_error_2  +=  pow(err, 2);
        }

cout << "BKD= " << total_background << endl;

      ofile << total_background << "$\\pm$" << pow(total_background_error, 0.5) << " & " ;


      for(unsigned int isign = 0;  isign<samplelist_signal.size(); isign++)
      {
         TH1D* tmp = histo_map[ regions[ireg]+"__"+samplelist_signal[isign]];
         double err = 0;
         ofile << tmp->IntegralAndError(1, tmp->GetNbinsX(), err) << "$\\pm$" << err << " & " ;
      }

      double total_DATA = 0;
      double total_DATA_error = 0;
      for(unsigned int idata = 0; idata <samplelist_data.size(); idata++)
      {
         TH1D* tmp = histo_map[ regions[ireg]+"__"+samplelist_data[idata]];
         double err = 0;
         ofile << tmp->IntegralAndError(1, tmp->GetNbinsX(), err) << "$\\pm$" << err << " & " ;
      }

      ofile << "\\\\ " << endl;
      ofile << "\\hline " << endl;
    }
    ofile << "\\hline " << endl;
    ofile << "\\hline " << endl;
  }

    ofile << "\\end{tabular} " << endl;
    ofile << "\\end{table} " << endl;
    //ofile << "\\end{landscape}" << endl;
    ofile << "\\end{document} " << endl;


}

void DrawCutFlow_postFit(){

  //CutFlow("FCNC");
  CutFlow("RES");

}
