#include "TString.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>
#include <vector>

using namespace std;

double getSFfit( TString thetasample, TString region, bool useAllRegions, bool useUpFit)
{
   TFile*               inputfile_prefit  = 0;
   //if(useAllRegions)    inputfile_prefit  = new TFile("inputTheta_merged_AllRegions_nosyst.root");
   if(useAllRegions)    inputfile_prefit  = new TFile("inputTheta_merged_AllRegions.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_merged_CRsOnly.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_mujets_merged_Selectedsignalregion.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_merged_SignalRegion.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_mujets_merged_ttbarregion_2j2b.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_WTTinterRegionSR.root");
   else                 inputfile_prefit  = new TFile("inputTest_noSignal.root");
   //else                 inputfile_prefit  = new TFile("interRegionTest/inputTheta_interRegionSR.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_interRegionOnly.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_merged_InterRegion.root");
   //else                 inputfile_prefit  = new TFile("inputTheta_mujets_merged_Wregion.root");

   TFile*                               inputfile_postfit = 0;
   if(useUpFit && !useAllRegions)       inputfile_postfit = new TFile( ("interRegionTest/histos_postFit_SignalAndInterRegion_"+thetasample+"_rate.root").Data() );
   //if(useUpFit && !useAllRegions)       inputfile_postfit = new TFile( ("testsWCRonly/histos_postFit_WCRonly_"+thetasample+"_rate.root").Data() );
   //if(useUpFit && !useAllRegions)       inputfile_postfit = new TFile( ("testsTTCRonly/histos_postFit_TTCRonly_"+thetasample+"_rate.root").Data() );
   //if(useUpFit && !useAllRegions)       inputfile_postfit = new TFile( ("testsSRonly/histos_postFit_SRonly_test_"+thetasample+"_rate.root").Data() );
   //if(useUpFit && !useAllRegions)       inputfile_postfit = new TFile( ("thetaInOut/histos_postFit_CRsOnly_"+thetasample+"_rate.root").Data() );
   else if(useUpFit && useAllRegions)   inputfile_postfit = new TFile( ("thetaInOut/histos_postFit_AllRegions_"+thetasample+"_rate.root").Data() );
   //else if(useUpFit && useAllRegions)   inputfile_postfit = new TFile( ("noSystTheta/histos_postFit_AllRegions_"+thetasample+"_rate.root").Data() );
   //else if(useAllRegions)               inputfile_postfit = new TFile(  "outputTheta_merged_AllRegions_2j2b_nosyst.root");
   else if(useAllRegions)               inputfile_postfit = new TFile(  "outputTheta_merged_AllRegions.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_merged_SignalRegionTest.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_merged_TTCRonly.root");
   else                                 inputfile_postfit = new TFile(  "outputTest_noSignal.root");
   //else                                 inputfile_postfit = new TFile(  "interRegionTest/outputTheta_interRegionSR.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_WTTinterRegionSR.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_interRegionOnly.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_merged_InterRegion.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_merged_WCRonly.root");
   //else                                 inputfile_postfit = new TFile(  "outputTheta_merged_CRsOnly.root");

   TString thetadistrib = "";
   bool isW = false;
   if( thetasample == "WExclc" || thetasample == "WExcll" || thetasample == "WExclb") isW = true;

   if(region == "Selectedsignalregion" && useAllRegions)                        thetadistrib = "mWTmujetsSelectedSignalregion";
   else if(region == "Wregion" )                                                thetadistrib = "mWTmujetsWregionHighpt";
   else if(region == "TTbarregion")                                             thetadistrib = "mWTmujetsttbarregionHighpt";
   //else if(region == "Selectedsignalregion" && thetasample == "TTMSDecays")     thetadistrib = "mWTmujetsttbarregionHighpt";
   //else if(region == "Selectedsignalregion" && isW )                            thetadistrib = "mWTmujetsWregionHighpt";
   else if(region == "Selectedsignalregion")                                    thetadistrib = "mWTmujetsSelectedSignalregion";
   //else if(region == "Selectedsignalregion")                                    thetadistrib = "mWTmujetsinterRegion";
   else if(region == "interRegion")                                             thetadistrib = "mWTmujetsinterRegion";

   TString whistoname_postfit = thetadistrib+"__"+thetasample;
   inputfile_postfit->cd();

   TH1D*       whisto_postfit = 0;
   if(thetadistrib != "") whisto_postfit = (TH1D*)inputfile_postfit->Get(whistoname_postfit)->Clone();

   TString whistoname_prefit  = thetadistrib+"__"+thetasample;

   inputfile_prefit->cd();
   TH1D*       whisto_prefit  = 0;
   if(thetadistrib != "") whisto_prefit = (TH1D*)inputfile_prefit->Get(whistoname_prefit)->Clone();

   double SF = 0;

   if( whisto_postfit != 0 && whisto_prefit != 0 && thetadistrib != "")            SF = whisto_postfit->Integral()/whisto_prefit->Integral();
   else if ( (whisto_postfit == 0 || whisto_prefit == 0) && thetadistrib != "")    SF = 0;
   else if ( thetadistrib == "")                                                   SF = 1;
   else cout << "Please check the 'getSFfit' function!" << endl;

   inputfile_prefit->Close();
   inputfile_postfit->Close();
   return SF;
}


void Macro_CheckFit(vector<TString> signalSample_list, vector<TString> mcSample_list, vector<TString> thetaSample_list, vector<int> colorVector, bool displaySignal, TString distrib, TString region, TString var, bool usePostFit, bool useAllRegions, bool useUpFit)
{

  bool setlogy = 0;
  if (displaySignal) setlogy = 1;

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);


  // For the Pad:
  gStyle->SetPadBorderMode(0);
  // ROOT . gStyle . SetPadBorderSize(Width_t size = 1);
  gStyle->SetPadColor(0); // kWhite
  gStyle->SetPadGridX(0); //false
  gStyle->SetPadGridY(0); //false
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);

  // For the frame:
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameFillStyle(0);
  gStyle->SetFrameLineColor(1);
  gStyle->SetFrameLineStyle(1);
  gStyle->SetFrameLineWidth(1);

  // For the histo:rebin
  // ROOT . gStyle . SetHistFillColor(1);
  // ROOT . gStyle . SetHistFillStyle(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  // ROOT . gStyle . SetLegoInnerR(Float_t rad = 0.5);
  // ROOT . gStyle . SetNumberContours(Int_t number = 20);

  gStyle->SetEndErrorSize(2);
  //ROOT . gStyle . SetErrorMarker(20);   /// I COMMENTED THIS OUT
  //ROOT . gStyle . SetErrorX(0.);
  //ROOT . gStyle . SetMarkerStyle(20);

  //For the fit/function:
  gStyle->SetOptFit(1011);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  //For the date:
  gStyle->SetOptDate(0);
  // ROOT . gStyle . SetDateX(Float_t x = 0.01);
  // ROOT . gStyle . SetDateY(Float_t y = 0.01);

  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(0); // kWhite
  gStyle->SetStatFont(42);
  //ROOT . gStyle . SetStatFontSize(0.025);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);
  // ROOT . gStyle . SetStatStyle(Style_t style = 1001);
  // ROOT . gStyle . SetStatX(Float_t x = 0);
  // ROOT . gStyle . SetStatY(Float_t y = 0);

  // Margins:
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  //ROOT . gStyle . SetPadRightMargin(0.12);
  gStyle->SetPadRightMargin(0.03);

  // For the Global title:

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
  // ROOT . gStyle . SetTitleH(0); // Set the height of the title box
  // ROOT . gStyle . SetTitleW(0); // Set the width of the title box
  // ROOT . gStyle . SetTitleX(0); // Set the position of the title box
  // ROOT . gStyle . SetTitleY(0.985); // Set the position of the title box
  // ROOT . gStyle . SetTitleStyle(Style_t style = 1001);
  // ROOT . gStyle . SetTitleBorderSize(2);

  // For the axis titles:

  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  // ROOT . gStyle . SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // ROOT . gStyle . SetTitleYSize(Float_t size = 0.02);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);
  // ROOT . gStyle . SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:

  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:

  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(1); // kTRUE
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  // Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);

  // Postscript options:
  gStyle->SetPaperSize(20.,20.);
  // ROOT . gStyle . SetLineScalePS(Float_t scale = 3);
  // ROOT . gStyle . SetLineStyleString(Int_t i, const char* text);
  // ROOT . gStyle . SetHeaderPS(const char* header);
  // ROOT . gStyle . SetTitlePS(const char* pstitle);

  // ROOT . gStyle . SetBarOffset(Float_t baroff = 0.5);
  // ROOT . gStyle . SetBarWidth(Float_t barwidth = 0.5);
  // ROOT . gStyle . SetPaintTextFormat(const char* format = "g");
  // ROOT . gStyle . SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // ROOT . gStyle . SetTimeOffset(Double_t toffset);
  // ROOT . gStyle . SetHistMinimumZero(kTRUE);


  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  c1->SetBottomMargin(0.3);
  c1->SetLogy(setlogy);
  c1->cd();


  vector<double> SF_fit;
  for (short unsigned int imc = 0; imc < thetaSample_list.size(); imc++)
  {
      double SF_tmp = 0;
      if(usePostFit)
      {
          if (region == "TTbarregion" && thetaSample_list[imc] == "QCD") SF_tmp = 0;
          else SF_tmp = getSFfit( thetaSample_list[imc], region, useAllRegions, useUpFit);

          if (!isfinite(SF_tmp) || isnan(SF_tmp)) SF_tmp = 0;
      }
      else                                        SF_tmp = 1;

      SF_fit.push_back(SF_tmp);
      cout << "SF_fit(" << thetaSample_list[imc] << ") = " << SF_tmp << endl;
  }
/*
  TFile * inputfile_toRescale ;
  inputfile_toRescale  = new TFile("../TreeReader/outputroot_withSyst/histo_merged.root");

  TString dataAname         = distrib+"__SingleMuA";
  TString dataBname         = distrib+"__SingleMuB";
  TString dataCname         = distrib+"__SingleMuC";
  TString dataDname         = distrib+"__SingleMuD";

  TH1D * distrib__RUNA  	= (TH1D*)inputfile_toRescale->Get(dataAname)->Clone();
  TH1D * distrib__RUNB  	= (TH1D*)inputfile_toRescale->Get(dataBname)->Clone();
  TH1D * distrib__RUNC  	= (TH1D*)inputfile_toRescale->Get(dataCname)->Clone();
  TH1D * distrib__RUND  	= (TH1D*)inputfile_toRescale->Get(dataDname)->Clone();

  TH1D * histo_data  	    = (TH1D*)distrib__RUNA;
         histo_data->Add(distrib__RUNB);
         histo_data->Add(distrib__RUNC);
         histo_data->Add(distrib__RUND);

  vector<TH1D *> histo_mcSamples;
  vector<TH1D *> signalSamples;

  for(unsigned int imc = 0; imc < mcSample_list.size(); imc++)
  {
      TString histo_mc_name   = distrib+"__"+mcSample_list[imc];

      TH1D * histo_tmp = (TH1D*)inputfile_toRescale->Get(histo_mc_name)->Clone();
      histo_tmp->Scale(SF_fit[imc]);
      histo_mcSamples.push_back(histo_tmp);
  }

  THStack *the_stack_histo= new THStack();
  for(unsigned int imc = 0; imc< histo_mcSamples.size(); imc++)
  {
      histo_mcSamples[imc]->SetFillStyle(1001);
      histo_mcSamples[imc]->SetFillColor(colorVector[imc]);
      histo_mcSamples[imc]->SetLineColor(colorVector[imc]);

      if(imc < histo_mcSamples.size() && colorVector[imc] != colorVector[imc+1] ) histo_mcSamples[imc]->SetLineColor(1);
      if(imc ==  histo_mcSamples.size())                                          histo_mcSamples[imc]->SetLineColor(1);
      the_stack_histo->Add(histo_mcSamples[imc]);
  }

  the_stack_histo->Draw("h");
  the_stack_histo->GetXaxis()->SetLabelSize(0.0);
  the_stack_histo->GetYaxis()->SetLabelSize(0.04);
  if(histo_data->GetMaximum() > the_stack_histo->GetMaximum() ) the_stack_histo->SetMaximum(     histo_data->GetMaximum()+0.1*     histo_data->GetMaximum());
  else                                                          the_stack_histo->SetMaximum(the_stack_histo->GetMaximum()+0.1*the_stack_histo->GetMaximum());

  if(setlogy) {the_stack_histo->SetMinimum(1.); the_stack_histo->SetMaximum(1000000);}
  if(var == "DeltaPhiLJ") the_stack_histo->GetXaxis()->SetRangeUser(0, 3.14);

  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(1.2);
  histo_data->SetLineColor(1);
  histo_data->Draw("epsame");
  if (displaySignal)
  {
    for( short unsigned int isig = 0; isig < signalSample_list.size(); isig++)
    {
        TString signalname      = distrib+"__"+signalSample_list[isig];
        TH1D * distrib__signal  = (TH1D*)inputfile_toRescale->Get(signalname)->Clone();
        if(distrib__signal != 0)
        {
            if(isig == 0)      distrib__signal->SetLineColor(kBlue+1);
            else if(isig == 1) distrib__signal->SetLineColor(kBlack);
            else if(isig == 2) distrib__signal->SetLineColor(kMagenta-5);
            else               distrib__signal->SetLineColor((int) isig);

            distrib__signal->SetLineStyle(2);
            distrib__signal->SetLineWidth(3);
            distrib__signal->Draw("hsame");
            signalSamples.push_back(distrib__signal);
        }
    }
  }

  //--------------------------
  //MC systematic plot
  //--------------------------
  TH1D * histo_syst_MC   = (TH1D*)(the_stack_histo->GetHistogram() )->Clone();
  for(unsigned int imc=0; imc< histo_mcSamples.size(); imc++)
  {
      histo_syst_MC->Add(histo_mcSamples[imc]);
  }

  TGraphErrors *thegraph = new TGraphErrors(histo_syst_MC);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  thegraph->Draw("e2same");


  //-------------------
  //legend and captions
  //-------------------
  TLatex *latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  latex->SetTextAlign(31);
  latex->DrawLatex(0.45, 0.95, "CMS Preliminary");


  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(31);
  latex2->DrawLatex(0.87, 0.95, "19.7 fb^{-1} at #sqrt{s} = 8 TeV");

  TString info_data = "#mu channel";


  TLatex * text2 = new TLatex(0.45,0.98, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.18);
  if (var == "DeltaPhiLJ" && region != "Selectedsignalregion") text2->SetX(0.70);
  text2->SetY(0.92);
  //text2->SetLineWidth(2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.046);
  //    text2->SetTextSizePixels(24);// dflt=28
  text2->Draw();

  TLegend*  qw;
  if( var == "DeltaPhiLJ" && region != "Selectedsignalregion")  qw = new TLegend(.35,.60,.50,.90);
  else                                                          qw = new TLegend(.70,.60,.90,.90);

  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  qw->AddEntry(histo_data,         "Data" ,                "ep");

  for(unsigned int i=0; i<mcSample_list.size(); i++)
  {
      if( mcSample_list[i] == "TTMSDecays_central"              )   qw->AddEntry( histo_mcSamples[i],  "t#bar{t}"	,"f");
      if( mcSample_list[i] == "WExclb"                          )   qw->AddEntry( histo_mcSamples[i],  "W+bjets"	,"f");
      if( mcSample_list[i] == "WExclc"                          )   qw->AddEntry( histo_mcSamples[i],  "W+cjets"	,"f");
      if( mcSample_list[i] == "WExcll"                          )   qw->AddEntry( histo_mcSamples[i],  "W+ljets"	,"f");
      if( mcSample_list[i] == "DYJetsToLL_M-50"                 )   qw->AddEntry( histo_mcSamples[i],  "DY"		    ,"f");
      if( mcSample_list[i] == "T_s" 		                    )   qw->AddEntry( histo_mcSamples[i+4],"Single top" ,"f");
      if( mcSample_list[i] == "WZ"                              )   qw->AddEntry( histo_mcSamples[i+2],"VV"	        ,"f");
      if( mcSample_list[i] == "QCD"  && region != "TTbarregion" )   qw->AddEntry( histo_mcSamples[i],  "QCD"	    ,"f");
   }


  if (displaySignal)
  {
    for( short unsigned int isig = 0; isig < signalSamples.size(); isig++)
    {
        qw->AddEntry( signalSamples[isig],signalSample_list[isig]	,"l");
    }
  }

  qw->Draw();

  //--------------------------
  //Data over background ratio
  //--------------------------
  TPad *canvas_2 = new TPad("canvas_2", "canvas_2", 0.0, 0.0, 1.0, 0.9);
  canvas_2->SetTopMargin(0.7);
  canvas_2->SetFillColor(0);
  canvas_2->SetFillStyle(0);
  canvas_2->SetGridy(1);
  canvas_2->Draw();
  canvas_2->cd(0);

  TH1D * histo_ratio_data = (TH1D*)histo_data->Clone();
  histo_ratio_data->Divide(histo_syst_MC);

  TGraphErrors *thegraph_tmp = (TGraphErrors*) thegraph->Clone();

  double *theErrorX  = thegraph_tmp->GetEX();
  double *theErrorY  = thegraph_tmp->GetEY();
  double *theY       = thegraph_tmp->GetY() ;
  for(short int i=0; i<thegraph_tmp->GetN(); i++)
  {
      if(theY[i]!=0) theErrorY[i] = theErrorY[i]/theY[i];
      theY[i]=1;
  }

  TGraphErrors *thegraph_ratio = new TGraphErrors(thegraph_tmp->GetN(), thegraph_tmp->GetX() , theY ,  thegraph_tmp->GetEX(),     thegraph_tmp->GetEY() );
  thegraph_ratio->SetFillStyle(3005);
  thegraph_ratio->SetFillColor(1);

  histo_ratio_data->SetMinimum(0.5);
  histo_ratio_data->SetMaximum(1.5);
  histo_ratio_data->GetXaxis()->SetTitleOffset(1.2);

  if     (var == "mWT")        histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else if(var == "MET")        histo_ratio_data->GetXaxis()->SetTitle("MET [GeV]");
  else if(var == "ptW")        histo_ratio_data->GetXaxis()->SetTitle("p(W)_{T} [GeV]");
  else if(var == "DeltaPhiLJ")
  {
      histo_ratio_data->GetXaxis()->SetTitle("#Delta#phi(lep - lead. jet) ");
      histo_ratio_data->GetXaxis()->SetRangeUser(0, 3.14);
  }
  else if(var == "NJet")       histo_ratio_data->GetXaxis()->SetTitle("jet mult.");
  else if(var == "NBJet")      histo_ratio_data->GetXaxis()->SetTitle("b-tagged jet mult.");
  else if(var == "mWT_full")   histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else if(var == "mWTplusMET") histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} + missing E_{T} [GeV]");
  else if(var == "mW")	       histo_ratio_data->GetXaxis()->SetTitle("m(W) [GeV]");
  else if(var == "etaW")	   histo_ratio_data->GetXaxis()->SetTitle("#eta(W) ");
  else if(var == "JetPt")	   histo_ratio_data->GetXaxis()->SetTitle("p_{T}(jets) [GeV]");
  else if(var == "LeadJetPt")	       histo_ratio_data->GetXaxis()->SetTitle("p_{T}(lead. jet) [GeV]");
  else if(var == "SubLeadJetPt")	   histo_ratio_data->GetXaxis()->SetTitle("p_{T}(sublead. jet) [GeV]");
  else if(var == "JetEta")	   histo_ratio_data->GetXaxis()->SetTitle("#eta(jets)");
  else if(var == "LeadJetEta")	       histo_ratio_data->GetXaxis()->SetTitle("#eta(lead. jet)");
  else if(var == "SubLeadJetEta")	   histo_ratio_data->GetXaxis()->SetTitle("#eta(sublead. jet)");
  else if(var == "LeptPt")	   histo_ratio_data->GetXaxis()->SetTitle("p_{T}(lepton) [GeV]");
  else if(var == "LeptEta")	   histo_ratio_data->GetXaxis()->SetTitle("#eta(lepton)");
  else if(var == "LeptIso")	   histo_ratio_data->GetXaxis()->SetTitle("iso(lepton)");
  else if(var == "BJetCSV")    histo_ratio_data->GetXaxis()->SetTitle("CSV discriminator");
  else if(var == "LeadJetBtagDiscr")histo_ratio_data->GetXaxis()->SetTitle("Btag discriminator(lead jet)");


  histo_ratio_data->GetXaxis()->SetLabelSize(0.04);
  histo_ratio_data->GetYaxis()->SetLabelSize(0.03);
  histo_ratio_data->GetYaxis()->SetNdivisions(6);
  histo_ratio_data->GetYaxis()->SetTitleSize(0.03);
  histo_ratio_data->Draw("E1X0");
  thegraph_ratio->Draw("e2same");

  TString outputname;
  if(usePostFit && useAllRegions)       outputname = "PostFitbySF_finalVersion_AllRegions_"+var+"_mujets";
  else if(usePostFit && !useAllRegions) outputname = "interRegionTest/PostFitbySF_interRegionSR_"+var+"_mujets";
  //else if(usePostFit && !useAllRegions) outputname = "PostFitbySF_finalVersion_CRsOnly_"+var+"_mujets";
  else                                  outputname = "PreFitbySF_finalVersion_" +var+"_mujets";

  TString endname_png;
  TString endname_pdf;

  if(!displaySignal)        outputname += "_noSignal";

  if (setlogy) { endname_png = "_logY.png"; endname_pdf = "_logY.pdf"; }
  else         { endname_png = ".png";      endname_pdf = ".pdf"; }

  //if(!useUpFit) c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+endname_png).Data());
  //if(!useUpFit) c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+endname_pdf).Data());
  if(!useUpFit) c1->SaveAs( (outputname+"_"+region+endname_png).Data());
  if(!useUpFit) c1->SaveAs( (outputname+"_"+region+endname_pdf).Data());
*/
}


void Macro_CheckFit(){

  bool displaySignal = false;
  bool usePostFit    = true;
  bool useAllRegions = false;
  bool useUpFit      = false;

  //-------------------------
  //define list of signal samples
  vector<TString> signalSample_list;

  if( displaySignal )
  {
    signalSample_list.push_back("S1Res700Inv100");
    signalSample_list.push_back("S4Inv500");
  }


  //-------------------------
  //define list of MC samples
  vector<TString> mcSample_list;
  vector<TString> thetaSample_list;
  vector<int> colorVector;
  mcSample_list.push_back("TTMSDecays_central");  colorVector.push_back(kRed+1);         thetaSample_list.push_back("TTMSDecays");
  mcSample_list.push_back("WExclb");              colorVector.push_back(kGreen+4);       thetaSample_list.push_back("WExclb");
  mcSample_list.push_back("WExclc");              colorVector.push_back(kGreen);         thetaSample_list.push_back("WExclc");
  mcSample_list.push_back("WExcll");              colorVector.push_back(kGreen-2);       thetaSample_list.push_back("WExcll");
  mcSample_list.push_back("DYJetsToLL_M-10To50"); colorVector.push_back(kAzure-2);       thetaSample_list.push_back("DY");
  mcSample_list.push_back("DYJetsToLL_M-50");     colorVector.push_back(kAzure-2);       thetaSample_list.push_back("DY");
  mcSample_list.push_back("T_s");                 colorVector.push_back(13);             thetaSample_list.push_back("SingleTop");
  mcSample_list.push_back("T_t");                 colorVector.push_back(13);             thetaSample_list.push_back("SingleTop");
  mcSample_list.push_back("Tbar_t");              colorVector.push_back(13);             thetaSample_list.push_back("SingleTop");
  mcSample_list.push_back("T_tW");                colorVector.push_back(13);             thetaSample_list.push_back("SingleTopW");
  mcSample_list.push_back("Tbar_tW");             colorVector.push_back(13);             thetaSample_list.push_back("SingleTopW");
  mcSample_list.push_back("WZ");                  colorVector.push_back(kRed+2);         thetaSample_list.push_back("VV");
  mcSample_list.push_back("WW");                  colorVector.push_back(kRed+2);         thetaSample_list.push_back("VV");
  mcSample_list.push_back("ZZ");                  colorVector.push_back(kRed+2);         thetaSample_list.push_back("VV");
  mcSample_list.push_back("QCD");                 colorVector.push_back(kYellow+1);      thetaSample_list.push_back("QCD");



//----------------------------------------------------------------------------//
//---------- mWT MET DeltaPhiLJ JetEta JetPt LeptPt LeptEta ptW NVtx----------//
//----------------------------------------------------------------------------//
  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "mWT_mujets_Selectedsignalregion", "Selectedsignalregion", "mWT", usePostFit, useAllRegions, useUpFit);
  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "mWT_mujets_ttbarregion_2j2b",            "TTbarregion",   "mWT", usePostFit, useAllRegions, useUpFit);
  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "mWT_mujets_Wregion_highpt",                    "Wregion", "mWT", usePostFit, useAllRegions, useUpFit);
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "mWT_mujets_interRegion",                    "interRegion", "mWT", usePostFit, useAllRegions, useUpFit);
//  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "LeadJetPt_mujets_ttbarregion_2j2b",            "TTbarregion",   "LeadJetPt", usePostFit, useAllRegions, useUpFit);
//  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, displaySignal, "LeadJetPt_mujets_ttbarregion_2j2b",            "TTbarregion",   "LeadJetPt", usePostFit, useAllRegions, useUpFit);


}
