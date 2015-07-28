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


  vector<TString > samplelistInput_MC;
  vector<TString > samplelist_MC;
  samplelistInput_MC.push_back("TTbar_Madgraph"       );     samplelist_MC.push_back("TTbar"     );
  samplelistInput_MC.push_back("WExclb"               );     samplelist_MC.push_back("WExclb"    );
  samplelistInput_MC.push_back("WExclc"               );     samplelist_MC.push_back("WExclc"    );
  samplelistInput_MC.push_back("WExcll"               );     samplelist_MC.push_back("WExcll"    );
  samplelistInput_MC.push_back("DYJetsToLL_M-10To50"  );     samplelist_MC.push_back("DY"        );
  samplelistInput_MC.push_back("DYJetsToLL_M-50"      );     samplelist_MC.push_back("DY"        );
  samplelistInput_MC.push_back("T_s"                  );     samplelist_MC.push_back("SingleTop" );
  samplelistInput_MC.push_back("T_t"                  );     samplelist_MC.push_back("SingleTop" );
  samplelistInput_MC.push_back("T_tW"                 );     samplelist_MC.push_back("SingleTop" );
  samplelistInput_MC.push_back("Tbar_t"               );     samplelist_MC.push_back("SingleTop" );
  samplelistInput_MC.push_back("Tbar_tW"              );     samplelist_MC.push_back("SingleTop" );
  samplelistInput_MC.push_back("WZ"                   );     samplelist_MC.push_back("VV"        );
  samplelistInput_MC.push_back("WW"                   );     samplelist_MC.push_back("VV"        );
  samplelistInput_MC.push_back("ZZ"                   );     samplelist_MC.push_back("VV"        );
  samplelistInput_MC.push_back("QCD"                  );     samplelist_MC.push_back("QCD"       );



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
    samplelist_signal.push_back("S1Res300Inv100");
    samplelist_signal.push_back("S1Res500Inv100");
    samplelist_signal.push_back("S1Res700Inv100");
    samplelist_signal.push_back("S1Res900Inv100");
    samplelist_signal.push_back("S1Res1300Inv100");
    samplelist_signal.push_back("S1Res1500Inv100");
  }

  vector<TString > samplelistInput_data;
  vector<TString > samplelist_data;
  samplelistInput_data.push_back("SingleMuA");                  samplelist_data.push_back("Data");
  samplelistInput_data.push_back("SingleMuB");                  samplelist_data.push_back("Data");
  samplelistInput_data.push_back("SingleMuC");                  samplelist_data.push_back("Data");
  samplelistInput_data.push_back("SingleMuD");                  samplelist_data.push_back("Data");

  map<TString, TH1D* > histo_map;

  //*********************************
  //get histograms
  //*********************************

  TFile *f_data  = new TFile("../../TreeReader/outputroot_withSyst/histo_merged.root");
  f_data->cd();

  for(unsigned int ichan=0; ichan<channels.size(); ichan++)
  {
    for(unsigned int idata=0; idata<samplelistInput_data.size(); idata++)
    {
        histo_map[ channels[ichan]+"_"+samplelistInput_data[idata]] = (TH1D*)gROOT->FindObject( ("CutFlow_"+channels[ichan]+"___"+samplelistInput_data[idata]).Data() );
        if(histo_map[ channels[ichan]+"_"+samplelistInput_data[idata]] == 0) cout << " non existing histogram " << ("CutFlow_"+channels[ichan]+"___"+samplelistInput_data[idata]) << endl;
    }
  }


  for(unsigned int ichan=0; ichan<channels.size(); ichan++)
  {
    for(unsigned int iMC=0; iMC<samplelistInput_MC.size(); iMC++)
    {
        histo_map[ channels[ichan]+"_"+samplelistInput_MC[iMC]] = (TH1D*)gROOT->FindObject( ("CutFlow_"+channels[ichan]+"___"+samplelistInput_MC[iMC]).Data() );
        if(histo_map[ channels[ichan]+"_"+samplelistInput_MC[iMC]]  == 0) cout << "non existing histogram " << "CutFlow_"+channels[ichan]+"___"+samplelistInput_MC[iMC] << endl;
    }
  }


  for(unsigned int ichan=0; ichan<channels.size(); ichan++)
  {
    for(unsigned int iMC=0; iMC<samplelistInput_signal.size(); iMC++)
    {
        histo_map[ channels[ichan]+"_"+samplelistInput_signal[iMC]] = (TH1D*)gROOT->FindObject( ("CutFlow_"+channels[ichan]+"___"+samplelistInput_signal[iMC]).Data() );
    }
  }


  string ofilenametex = "CutFlow_light.tex";
  ofstream ofile(ofilenametex.c_str());

  ofile<<"\\documentclass[10pt]{article}"<<endl;
  ofile<<"\\usepackage{amsthm}"<<endl;
  ofile<<"\\usepackage{color}"<<endl;
  //ofile<<"\\documentclass[amsmath,amssymb]{revtex4}"<<endl;

  ofile<<"\\begin{document}"<<endl;

  ofile.setf(ios::fixed);
  ofile.precision(1);
  vector<string> CutName;
  CutName.push_back("Presel.");
  CutName.push_back("Signal Region");
  CutName.push_back("TTbar-enriched CR");
  CutName.push_back("W+Jets-enriched CR");
  //CutName.push_back("ATLAS FCNC SR");
  //CutName.push_back("ATLAS RES SR");

  //ofile.precision(3);
  //ofile << "\\clearpage" << endl;
  //ofile << "\\begin{landscape}" << endl;
  ofile << "\\begin{table}[p]" << endl;

  ofile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  //ofile << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|}" << endl;
  ofile << "\\hline" << endl;
  ofile << "\\hline" << endl;



  for(unsigned int ichan = 0; ichan < channels.size(); ichan++)
  {
     ofile << "% for channel   " <<  channels[ichan] << endl;
     ofile << "      & "  ;
     for(unsigned int iback = 0; iback <samplelist_MC.size(); iback++)
     {
        if(iback == 0)                                                  ofile << samplelist_MC[iback] << "  & " ;
        if(iback > 0 && samplelist_MC[iback] != samplelist_MC[iback-1]) ofile << samplelist_MC[iback] << "  & " ;
     }



    ofile << "total back. & " ;
    //cout total background
    for(unsigned int isign = 0;  isign<samplelist_signal.size(); isign++)
    {
        ofile << samplelist_signal[isign] << "  & " ;
    }

    for(unsigned int idata = 0; idata <samplelist_data.size(); idata++)
    {
        if(idata == 0)                                                      ofile << samplelist_data[idata] << "  \\\\ "  << endl;;
        if(idata > 0 && samplelist_data[idata] != samplelist_data[idata-1]) ofile << samplelist_data[idata] << "  \\\\ "  << endl;;
    }


    ofile << "\\hline " << endl;
    for(unsigned int isel = 0; isel < CutName.size(); isel++)
    {
        unsigned short int bin = 0;
        if(isel == 0) bin = 5;   //Presel
        if(isel == 1) bin = 11;  //SR Prefit
        if(isel == 2) bin = 12;  //TT-CR Prefit
        if(isel == 3) bin = 13;  //W-CR Prefit
        //if(isel == 4) bin = 14;  //ATLAS FCNC SR
        //if(isel == 5) bin = 15;  //ATLAS RES SR

        if(CutName[isel] == "TTbar-enriched CR" || CutName[isel] == "W+Jets-enriched CR" || CutName[isel] == "Signal Region" ||  CutName[isel] == "ATLAS FCNC SR" || CutName[isel] == "ATLAS RES CR" ) ofile << "\\hline " << endl;
        //cout << "isel " << isel << endl;
        ofile <<CutName[isel] <<" & ";
        double total_background = 0;
        double total_background_error = 0;
        double total_DY = 0;
        double total_DY_error = 0;
        double total_SingleTop = 0;
        double total_SingleTop_error = 0;
        double total_VV = 0;
        double total_VV_error = 0;

        for(unsigned int iback = 0; iback <samplelistInput_MC.size(); iback++)
        {
            TH1D* tmp = histo_map[ channels[ichan]+"_"+samplelistInput_MC[iback]];
            if(iback < 4)
            {
                cout << "sample = " << samplelistInput_MC[iback] << " | and BinContent(" << bin << ")=" << tmp->GetBinContent(bin) << endl;
                ofile << tmp->GetBinContent(bin) << "$\\pm$" << tmp->GetBinError(bin)  << " & " ;
            }
            else if(iback == 4 || iback == 5 )
            {
                total_DY +=           tmp->GetBinContent(bin);
                total_DY_error += pow(tmp->GetBinError(bin),2);
                if(iback == 5) ofile << total_DY << "$\\pm$" << pow(total_DY_error, 0.5) << " & " ;
            }
            else if(iback >= 6 && iback <= 10 )
            {
                cout << "sample = " << samplelistInput_MC[iback] << " | and BinContent(" << bin << ")=" << tmp->GetBinContent(bin) << endl;
                total_SingleTop +=           tmp->GetBinContent(bin);
                total_SingleTop_error += pow(tmp->GetBinError(bin),2);
                if(iback == 10) ofile << total_SingleTop << "$\\pm$" << pow(total_SingleTop_error, 0.5) << " & " ;
            }
            else if(iback >= 11 && iback <= 13 )
            {
                cout << "sample = " << samplelistInput_MC[iback] << " | and BinContent(" << bin << ")=" << tmp->GetBinContent(bin) << endl;
                total_VV +=           tmp->GetBinContent(bin);
                total_VV_error += pow(tmp->GetBinError(bin),2);
                if(iback == 13) ofile << total_VV << "$\\pm$" << pow(total_VV_error, 0.5) << " & " ;
            }
            else if(iback == 14)
            {
                if(tmp->GetBinContent(bin) < 0)  ofile << 0.0 << "$\\pm$" << 0.0  << " & " ;
                else                                ofile << tmp->GetBinContent(bin) << "$\\pm$" << tmp->GetBinError(bin)  << " & " ;
            }

	        total_background        +=      tmp->GetBinContent(bin);
	        total_background_error  +=  pow(tmp->GetBinError(bin), 2);
      }

cout << "BKD= " << total_background << endl;

      ofile << total_background << "$\\pm$" << pow(total_background_error, 0.5) << " & " ;


      for(unsigned int isign = 0;  isign<samplelistInput_signal.size(); isign++)
      {
         TH1D* tmp = histo_map[ channels[ichan]+"_"+samplelistInput_signal[isign]];
         ofile << tmp->GetBinContent(bin) << "$\\pm$" << tmp->GetBinError(bin)  << " & " ;
      }

      double total_DATA = 0;
      double total_DATA_error = 0;
      for(unsigned int idata = 0; idata <samplelistInput_data.size(); idata++)
      {
         TH1D* tmp = histo_map[ channels[ichan]+"_"+samplelistInput_data[idata]];
	     //cout << "isel " << isel  << " " << channels[ichan]+"_"+samplelistInput_data[idata] << " " << tmp->GetBinContent(isel+1) << endl;
         //if(idata <= (samplelistInput_data.size()-1) ) ofile << tmp->GetBinContent(isel+1) << "$\\pm$" << tmp->GetBinError(isel+1)  << " & " ;
	     //else  ofile << tmp->GetBinContent(isel+1) << "$\\pm$" << tmp->GetBinError(isel+1)   ;
	     //tmp->Draw();
	     total_DATA      +=      tmp->GetBinContent(bin);
	     total_DATA_error+=  pow(tmp->GetBinError(bin), 2);
      }

      ofile << total_DATA << "$\\pm$" << pow(total_DATA_error, 0.5)  ;
cout << "DATA= " << total_DATA << endl;
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

void DrawCutFlow_light(){

  //CutFlow("FCNC");
  CutFlow("RES");

}
