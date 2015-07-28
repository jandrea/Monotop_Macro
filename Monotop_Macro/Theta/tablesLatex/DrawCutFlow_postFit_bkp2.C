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
#include <sstream>
#include <vector>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <math.h>

#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

string floatToString(float f)
{
     stringstream stream;
     stream.str("");
     stream.precision(3);
     stream << f;
     return stream.str();
}

string replaceString(string subject, const string& search, const string& replace)
{
     size_t pos = 0;
     while((pos = subject.find(search, pos)) != string::npos)
     {
         subject.replace(pos, search.length(), replace);
         pos += replace.length();
     }
     return subject;
}

string simplifyWriting(float myString)
{
     string number = floatToString(myString) ;
     number = replaceString(number, "e+0", "$\\times 10^{");
     if (number.find("times") != string::npos) number += "}$";
     return number;
}

void CutFlow(TString signal_type){


  //cout << setprecision(4) << endl;
  vector<TString> channels;
  channels.push_back("mujets");

  vector<TString> regions;
  regions.push_back("mWTmujetsWregionHighpt");
  regions.push_back("mWTmujetsttbarregionHighpt");
  regions.push_back("mWTmujetsSelectedSignalregion");

  vector<TString > samplelist_MC;
  vector<TString > samplename_MC;
  samplelist_MC.push_back("TTMSDecays");     samplename_MC.push_back("\\bm{$t\\bar{t}$}");
  samplelist_MC.push_back("WExclb"    );     samplename_MC.push_back("\\bm{$W(b\\text{\\textbf{-flav.}})$}" );
  samplelist_MC.push_back("WExclc"    );     samplename_MC.push_back("\\bm{$W(c\\text{\\textbf{-flav.}})$}" );
  samplelist_MC.push_back("WExcll"    );     samplename_MC.push_back("\\bm{$W(l\\text{\\textbf{-flav.}})$}" );
  samplelist_MC.push_back("DY"        );     samplename_MC.push_back("\\textbf{Drell-Yan}"   );
  samplelist_MC.push_back("VV"        );     samplename_MC.push_back("\\textbf{Diboson}"    );
  samplelist_MC.push_back("QCD"       );     samplename_MC.push_back("\\textbf{QCD}"        );
  samplelist_MC.push_back("SingleTopW");     samplename_MC.push_back("\\textbf{SingleTopW}" );
  samplelist_MC.push_back("SingleTop" );     samplename_MC.push_back("\\textbf{SingleTop}"  );



  vector<TString > samplelist_signal;
  vector<TString > samplename_signal;
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
    //samplelist_signal.push_back("S1Res500Inv100");  samplename_signal.push_back("\\bm{$S1(m_{X} = 500, m_{\\chi} = 100)$}");
    samplelist_signal.push_back("S1Res500Inv100");  samplename_signal.push_back("\\bm{$S1(500, 100)$}");
    samplelist_signal.push_back("S1Res700Inv100");  samplename_signal.push_back("\\bm{$S1(700, 100)$}");
    samplelist_signal.push_back("S1Res900Inv100");  samplename_signal.push_back("\\bm{$S1(900, 100)$}");
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
        if(samplelist_MC[iMC] == "QCD" && regions[ireg] == "mWTmujetsttbarregionHighpt") histo_map[ regions[ireg]+"__"+samplelist_MC[iMC]] = 0 ;
        else histo_map[ regions[ireg]+"__"+samplelist_MC[iMC]] = (TH1D*)f_mc->Get( (regions[ireg]+"__"+samplelist_MC[iMC]).Data() )->Clone();
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

  ofile << setprecision(4) << endl;
  ofile<<"\\documentclass[10pt]{article}"<<endl;
  ofile<<"\\usepackage{amsmath}"<<endl;
  ofile<<"\\usepackage{amsthm}"<<endl;
  ofile<<"\\usepackage{color}"<<endl;
  ofile<<"\\usepackage{bm}"<<endl;
  //ofile<<"\\documentclass[amsmath,amssymb]{revtex4}"<<endl;

  ofile<<"\\begin{document}"<<endl;

  //ofile.setf(ios::fixed);
  //ofile.precision(1);
  vector<string> CutName;
  //CutName.push_back("Presel.");
  //CutName.push_back("$t\\bar{t}$");
  CutName.push_back("\\bm{$W\\text{\\textbf{-enriched CR}}$}");
  CutName.push_back("\\bm{$t\\bar{t}\\text{\\textbf{-enriched CR}}$}");
  CutName.push_back("\\textbf{Signal Region}");
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

  ofile << "\\\\" << endl;

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


  for(unsigned int iback = 0; iback < samplelist_MC.size(); iback++)
  {

    if(samplelist_MC[iback] != "SingleTopW") ofile << "\\hline " << endl;
    if(samplelist_MC[iback] != "SingleTopW") ofile << samplename_MC[iback] <<" & ";

    for(unsigned int ireg = 0; ireg < regions.size(); ireg++)
    {
        TH1D* tmp = 0;
        if(samplelist_MC[iback] != "QCD" || regions[ireg] != "mWTmujetsttbarregionHighpt") tmp = histo_map[ regions[ireg]+"__"+samplelist_MC[iback]];

        double err = 0;
        if(iback < 7)
        {
            if (tmp == 0) { ofile << "0 $\\pm$ 0 & " ; continue; }
            else
            {
                ofile << simplifyWriting(tmp->IntegralAndError(1, tmp->GetNbinsX(), err)) ;
                ofile << "$\\pm$" << simplifyWriting(err) ;
            }
            if(ireg == 0 || ireg == 1) ofile << " & " ;
            else                       ofile << "\\\\" << endl;
        }
        else if(iback == 7 || iback == 8)
        {
            if(ireg == 0)
            {
                total_SingleTop_0       += tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
                total_SingleTop_error_0 += pow(err,2);
                if(iback == 8) ofile << simplifyWriting(total_SingleTop_0) << "$\\pm$" << simplifyWriting(pow(total_SingleTop_error_0, 0.5)) << " & " ;
            }
            else if(ireg == 1)
            {
                total_SingleTop_1       += tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
                total_SingleTop_error_1 += pow(err,2);
                if(iback == 8) ofile << simplifyWriting(total_SingleTop_1) << "$\\pm$" << simplifyWriting(pow(total_SingleTop_error_1, 0.5)) << " & " ;
            }
            else if(ireg == 2)
            {
                total_SingleTop_2       += tmp->IntegralAndError(1, tmp->GetNbinsX(), err);
                total_SingleTop_error_2 += pow(err,2);
                if(iback == 8) ofile << simplifyWriting(total_SingleTop_2) << "$\\pm$" << simplifyWriting(pow(total_SingleTop_error_2, 0.5)) << "\\\\" << endl;
            }

        }

        //cout << "sample = " << samplelist_MC[iback] << " | and Yield =" << tmp->Integral() <<  " | err = " << err << endl;

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
    }
   }
      ofile << "\\hline " << endl;
      ofile << "\\hline " << endl;


      ofile << "\\textbf{total SM} & " << simplifyWriting(total_background_0) << "$\\pm$" << simplifyWriting(pow(total_background_error_0, 0.5)) << " & " << simplifyWriting(total_background_1) << "$\\pm$" << simplifyWriting(pow(total_background_error_1, 0.5)) << " & " << simplifyWriting(total_background_2) << "$\\pm$" << simplifyWriting(pow(total_background_error_2, 0.5)) << "\\\\" << endl;

      ofile << "\\hline " << endl;
      ofile << "\\hline " << endl;

      for(unsigned int idata = 0; idata <samplelist_data.size(); idata++)
      {
        ofile << "\\textbf{data}  &";
        for(unsigned int ireg = 0; ireg < regions.size(); ireg++)
        {
         TH1D* tmp = histo_map[ regions[ireg]+"__"+samplelist_data[idata]];
         double err = 0;
         ofile << simplifyWriting(tmp->IntegralAndError(1, tmp->GetNbinsX(), err));
         ofile << "$\\pm$" << simplifyWriting(err);
         if(ireg == 0 || ireg == 1) ofile << " & " ;
         else                       ofile << "\\\\" << endl;
        }

      }


      ofile << "\\hline " << endl;
      ofile << "\\hline " << endl;

      for(unsigned int isign = 0;  isign<samplelist_signal.size(); isign++)
      {
        ofile << samplename_signal[isign] << " & ";
        for(unsigned int ireg = 0; ireg < regions.size(); ireg++)
        {
         TH1D* tmp = histo_map[ regions[ireg]+"__"+samplelist_signal[isign]];
         double err = 0;
         ofile << simplifyWriting(tmp->IntegralAndError(1, tmp->GetNbinsX(), err));
         ofile << "$\\pm$" << simplifyWriting(err);
         if(ireg == 0 || ireg == 1) ofile << " & " ;
         else                       ofile << "\\\\" << endl;
        }
      }

    ofile << "\\hline " << endl;

    ofile << "\\end{tabular} " << endl;
    ofile << "\\end{table} " << endl;
    //ofile << "\\end{landscape}" << endl;
    ofile << "\\end{document} " << endl;


}

void DrawCutFlow_postFit(){

  //CutFlow("FCNC");
  CutFlow("RES");

}
