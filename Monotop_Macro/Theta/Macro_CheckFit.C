#include "TString.h"
#include "TH1F.h"
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

double getSFfit( TFile* inputfile_prefit, TFile* inputfile_postfit, TString inputdistrib, TString thetainputdistrib, TString sample, TString thetasample)
{
   TString whistoname_prefit  = inputdistrib+"__"+sample;
   TH1F*       whisto_prefit  = (TH1F*)inputfile_prefit->Get(whistoname_prefit)->Clone();

   TString whistoname_postfit = thetainputdistrib+"__"+thetasample;
   TH1F*       whisto_postfit = (TH1F*)inputfile_postfit->Get(whistoname_postfit)->Clone();

   double SF = 0;
   cout << "Yield ( " << inputdistrib << " ) = " << whisto_prefit->Integral() << " and yield ( " << thetainputdistrib << " ) = " << whisto_postfit->Integral() << endl;
   if( whisto_postfit != 0 && whisto_prefit != 0) SF = whisto_postfit->Integral()/whisto_prefit->Integral();
   else cout << "Please check the 'getSFfit' function!" << endl;

   return SF;
}


void Macro_CheckFit(std::vector<TString> signalSample_list, std::vector<TString> mcSample_list, std::vector<TString> thetaSample_list, std::vector<int> colorVector, bool usePostFit, bool displaySignal, TString inputdistrib , TString thetainputdistrib, TString region, TString var, TString outputdistrib){

  bool setlogy = 1;
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

  TFile * inputfile_prefit ;
  TFile * inputfile_postfit ;
  inputfile_prefit  = new TFile("../TreeReader/outputroot_withSyst/histo_merged.root");
  //inputfile_postfit = new TFile("histos-mle.root");
  inputfile_postfit = new TFile("histos-mle_bkp18_03_15.root");

  vector<TString> TMPlist = mcSample_list;
  if (usePostFit) TMPlist = thetaSample_list;
  vector<double> SF_fit;
  for (short unsigned int imc = 0; imc < TMPlist.size(); imc++)
  {
    if (usePostFit)
    {
        double SF_tmp = getSFfit( inputfile_prefit, inputfile_postfit, inputdistrib, thetainputdistrib, mcSample_list[imc], thetaSample_list[imc]);
        if (isfinite(SF_tmp) && !isnan(SF_tmp)) SF_fit.push_back(SF_tmp);
        else                                    SF_fit.push_back(0);
        cout << "SF_fit(" << thetaSample_list[imc] << ") = " << SF_tmp << endl;
    }
    else SF_fit.push_back(1);
  }

  TString dataAname         = outputdistrib+"__SingleMuA";
  TString dataBname         = outputdistrib+"__SingleMuB";
  TString dataCname         = outputdistrib+"__SingleMuC";
  TString dataDname         = outputdistrib+"__SingleMuD";

  TH1F * distrib__RUNA  	= (TH1F*)inputfile_prefit->Get(dataAname)->Clone();
  TH1F * distrib__RUNB  	= (TH1F*)inputfile_prefit->Get(dataBname)->Clone();
  TH1F * distrib__RUNC  	= (TH1F*)inputfile_prefit->Get(dataCname)->Clone();
  TH1F * distrib__RUND  	= (TH1F*)inputfile_prefit->Get(dataDname)->Clone();

  TH1F * histo_data  	    = (TH1F*)distrib__RUNA;
         histo_data->Add(distrib__RUNB);
         histo_data->Add(distrib__RUNC);
         histo_data->Add(distrib__RUND);

  //histo_data->Rebin(5);
  std::vector<TH1F *> histo_mcSamples;
  std::vector<TH1F *> signalSamples;

  for(unsigned int imc = 0; imc < mcSample_list.size(); imc++){

    TString histo_mc_name   = outputdistrib+"__"+mcSample_list[imc];
    TH1F * histo_tmp = (TH1F*)inputfile_prefit->Get(histo_mc_name);
    histo_tmp->Scale(SF_fit[imc]);
    histo_mcSamples.push_back(histo_tmp);
  }

  THStack *the_stack_histo= new THStack();
  for(unsigned int imc = 0; imc< histo_mcSamples.size(); imc++){
    histo_mcSamples[imc]->SetFillStyle(1001);
    histo_mcSamples[imc]->SetFillColor(colorVector[imc]);
    histo_mcSamples[imc]->SetLineColor(colorVector[imc]);
    if(imc < histo_mcSamples.size() && colorVector[imc] != colorVector[imc+1] ) histo_mcSamples[imc]->SetLineColor(1);
    if(imc ==  histo_mcSamples.size())                                          histo_mcSamples[imc]->SetLineColor(1);
    //histo_mcSamples[imc]->Rebin(5);
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
        TString signalname      = outputdistrib+"__"+signalSample_list[isig];
        TH1F * distrib__signal  = (TH1F*)inputfile_prefit->Get(signalname)->Clone();
        if(distrib__signal != 0)
        {
            //cout << "histo: " << signalname << " | integral = " << distrib__signal->Integral() << endl;
            distrib__signal->SetLineColor((int) isig + 2);
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
  TH1F * histo_syst_MC   = (TH1F*)(the_stack_histo->GetHistogram() )->Clone();
  for(unsigned int imc=0; imc< histo_mcSamples.size(); imc++){
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

  TString info_data = "mu channel";


  TLatex * text2 = new TLatex(0.45,0.98, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.18);
  if (var == "DeltaPhiLJ") text2->SetX(0.70);
  text2->SetY(0.92);
  //text2->SetLineWidth(2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.0610687);
  //    text2->SetTextSizePixels(24);// dflt=28
  text2->Draw();

  TLegend*  qw;
  if( var == "DeltaPhiLJ")  qw = new TLegend(.35,.60,.50,.90);
  else                      qw = new TLegend(.75,.60,.90,.90);

  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);

  qw->AddEntry(histo_data,         "Data" ,                "ep");
  for(unsigned int i=0; i<mcSample_list.size(); i++)
  {
    if( (mcSample_list[i] == "TTbar_Madgraph"  ) || ( mcSample_list[i] == "TTbarMadgraph") ) qw->AddEntry( histo_mcSamples[i],	"t#bar{t}"	,"f");
    if(  mcSample_list[i] == "WExclb"          )                                             qw->AddEntry( histo_mcSamples[i],	"W+bjets"	,"f");
    if(  mcSample_list[i] == "WExclc"          )                                             qw->AddEntry( histo_mcSamples[i],	"W+cjets"	,"f");
    if(  mcSample_list[i] == "WExcll"          )                                             qw->AddEntry( histo_mcSamples[i],	"W+ljets"	,"f");
    if(  mcSample_list[i] == "WExcl"           )                                             qw->AddEntry( histo_mcSamples[i],	"W+jets"	,"f");
    if( (mcSample_list[i] == "DYJetsToLL_M-50" ) || ( mcSample_list[i] == "DY50"         ) ) qw->AddEntry( histo_mcSamples[i],  "DY"		,"f");
    if( (mcSample_list[i] == "T_s" 		       ) || ( mcSample_list[i] == "Ts" 			 ) ) qw->AddEntry( histo_mcSamples[i+4],"single top","f");
    if( (mcSample_list[i] == "WZ"              ) || ( mcSample_list[i] == "WZ"           ) ) qw->AddEntry( histo_mcSamples[i+2],"VV"	    ,"f");
    if( (mcSample_list[i] == "QCD_A"           )                                           ) qw->AddEntry( histo_mcSamples[i+3],"QCD"	    ,"f");
    if( (mcSample_list[i] == "QCD"             )                                           ) qw->AddEntry( histo_mcSamples[i],  "QCD"	    ,"f");
  }

  if (displaySignal)
  {
    for( short unsigned int isig = 0; isig < signalSamples.size(); isig++)
    {
        qw->AddEntry( signalSamples[isig],signalSample_list[isig]	    ,"l");
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


  TH1F * histo_ratio_data = (TH1F*)histo_data->Clone();
  histo_ratio_data->Divide(histo_syst_MC);

  TGraphErrors *thegraph_tmp = (TGraphErrors*) thegraph->Clone();

  double *theErrorX  = thegraph_tmp->GetEX();
  double *theErrorY  = thegraph_tmp->GetEY();
  double *theY       = thegraph_tmp->GetY() ;
  for(short int i=0; i<thegraph_tmp->GetN(); i++){
    theErrorY[i] = theErrorY[i]/theY[i];
    theY[i]=1;
  }

  TGraphErrors *thegraph_ratio = new TGraphErrors(thegraph_tmp->GetN(), thegraph_tmp->GetX() , theY ,  thegraph_tmp->GetEX(),     thegraph_tmp->GetEY() );
  thegraph_ratio->SetFillStyle(3005);
  thegraph_ratio->SetFillColor(1);

  histo_ratio_data->SetMinimum(0.5);
  histo_ratio_data->SetMaximum(1.5);
  histo_ratio_data->GetXaxis()->SetTitleOffset(1.2);
  if      (var == "mWT") histo_ratio_data->GetXaxis()->SetTitle("m(W)_{T} [GeV]");
  else if (var == "MET") histo_ratio_data->GetXaxis()->SetTitle("MET [GeV]");
  else if (var == "ptW") histo_ratio_data->GetXaxis()->SetTitle("p(W)_{T} [GeV]");
  else if (var == "JetPt") histo_ratio_data->GetXaxis()->SetTitle("p_{T} (jets) [GeV]");
  else if (var == "ptLept") histo_ratio_data->GetXaxis()->SetTitle("p_{T} (#mu) [GeV]");
  else if (var == "DeltaPhiLJ") { histo_ratio_data->GetXaxis()->SetTitle("#Delta #phi (lep - lead. jet) "); histo_ratio_data->GetXaxis()->SetRangeUser(0, 3.14); }
  histo_ratio_data->GetXaxis()->SetLabelSize(0.04);
  histo_ratio_data->GetYaxis()->SetLabelSize(0.03);
  histo_ratio_data->GetYaxis()->SetNdivisions(6);
  histo_ratio_data->GetYaxis()->SetTitleSize(0.03);
  histo_ratio_data->Draw("E1X0");
  thegraph_ratio->Draw("e2same");

  TString outputname;
  if (!usePostFit) outputname = "PreFit_"+var+"_mujets";
  else             outputname = "PostFit_"+var+"_mujets";

  TString endname_png;
  TString endname_pdf;
  //if (setlogy) endname = "_logY.pdf";
  if (setlogy) { endname_png = "_logY_.png"; endname_pdf = "_logY_.pdf"; }
  else         { endname_png = "_.png";      endname_pdf = "_.pdf"; }

  c1->SaveAs( ("plots_fits/"+outputname+"_"+region+endname_png).Data());
  c1->SaveAs( ("plots_fits/"+outputname+"_"+region+endname_pdf).Data());

}


void Macro_CheckFit(){

  bool usePostFit = true;
  bool displaySignal = false;

  //-------------------------
  //define list of signal samples
  std::vector<TString> signalSample_list;

  if( displaySignal )
  {
    //signalSample_list.push_back("S1_1000_100");
    //signalSample_list.push_back("S1_1000_800");
    //signalSample_list.push_back("S1_500_100");
   // signalSample_list.push_back("S4_400");
   // signalSample_list.push_back("S4_600");
   // signalSample_list.push_back("S4_700");
  }


  //-------------------------
  //define list of MC samples
  std::vector<TString> thetaSample_list;
  std::vector<TString> mcSample_list;
  std::vector<int> colorVector;
  mcSample_list.push_back("TTbar_Madgraph");     thetaSample_list.push_back("TTbarMadgraph");   colorVector.push_back(kRed+1);
  mcSample_list.push_back("WExclb");             thetaSample_list.push_back("WExclb");          colorVector.push_back(kGreen-2);
  mcSample_list.push_back("WExclc");             thetaSample_list.push_back("WExclc");          colorVector.push_back(kGreen);
  mcSample_list.push_back("WExcll");             thetaSample_list.push_back("WExcll");          colorVector.push_back(kGreen+4);
  mcSample_list.push_back("DYJetsToLL_M-10To50");thetaSample_list.push_back("DY10To50");        colorVector.push_back(kAzure-2);
  mcSample_list.push_back("DYJetsToLL_M-50");    thetaSample_list.push_back("DY50");            colorVector.push_back(kAzure-2);
  mcSample_list.push_back("T_s");                thetaSample_list.push_back("Ts");              colorVector.push_back(kRed+2);
  mcSample_list.push_back("T_t");                thetaSample_list.push_back("Tt");              colorVector.push_back(kRed+2);
  mcSample_list.push_back("T_tW");               thetaSample_list.push_back("TtW");             colorVector.push_back(kRed+2);
  mcSample_list.push_back("Tbar_t");             thetaSample_list.push_back("Tbart");           colorVector.push_back(kRed+2);
  mcSample_list.push_back("Tbar_tW");            thetaSample_list.push_back("TbartW");          colorVector.push_back(kRed+2);
  mcSample_list.push_back("WZ");                 thetaSample_list.push_back("WZ");              colorVector.push_back(13);
  mcSample_list.push_back("WW");                 thetaSample_list.push_back("WW");              colorVector.push_back(13);
  mcSample_list.push_back("ZZ");                 thetaSample_list.push_back("ZZ");              colorVector.push_back(13);

  mcSample_list.push_back("QCD_A");             thetaSample_list.push_back("QCD");             colorVector.push_back(kYellow+1);
  mcSample_list.push_back("QCD_B");                                                            colorVector.push_back(kYellow+1);
  mcSample_list.push_back("QCD_C");                                                            colorVector.push_back(kYellow+1);
  mcSample_list.push_back("QCD_D");                                                            colorVector.push_back(kYellow+1);


//--------------------------//
//------- Delta Phi --------//
//--------------------------//
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_signalregion", "mWTmujetsSignalregion", "signalregion", "DeltaPhiLJ", "DeltaPhiLJ_mujets_signalregion");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt", "TTbarregion", "DeltaPhiLJ", "DeltaPhiLJ_mujets_ttbarregion_highpt");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "DeltaPhiLJ", "DeltaPhiLJ_mujets_Wregion_highpt");


//--------------------------//
//---------- MET -----------//
//--------------------------//
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Selectedsignalregion", "mWTmujetsSelectedSignalregion", "Selectedsignalregion", "MET", "MET_mujets_Selectedsignalregion");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt", "TTbarregion", "MET", "MET_mujets_ttbarregion_highpt");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "MET", "MET_mujets_Wregion_highpt");

//--------------------------//
//---------- pWT -----------//
//--------------------------//
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_signalregion", "mWTmujetsSignalregion", "signalregion", "ptW", "ptW_mujets_signalregion");
  ////Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt", "TTbarregion", "ptW", "ptW_mujets_ttbarregion_highpt");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "ptW", "ptW_mujets_Wregion_highpt");


  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_signalregion", "mWTmujetsSignalregion", "signalregion", "JetPt", "JetPt_mujets_afterleptsel");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_signalregion", "mWTmujetsSignalregion", "signalregion", "mWT", "mWT_mujets_signalregion");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "ptLept", "LeptPt_mujets_Wregion_highpt");
  //Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "mWT", "mWT_mujets_Wregion_highpt");
  Macro_CheckFit(signalSample_list, mcSample_list, thetaSample_list, colorVector, usePostFit, displaySignal, "mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", "Wregion", "DeltaPhiLMet", "DeltaPhiLMet_mujets_Wregion_highpt");
}
