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

#include "myCMS_lumi.h"

using namespace std;

void PlotFitResults(  vector<TString> signalSample_list, vector<TString> mcSample_list, vector<int> colorVector, bool usePostFit, bool displaySignal, TString inputdistrib , TString thetainputdistrib , TString inputfilename, TString region, TString variable, bool useAllRegions, bool useMergedThetaSamples, bool useElectronChannel, TString ttCRJetMult)
{

  Int_t stati=0;
  Bool_t  fit=1;
  Bool_t logy=0;

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
  TFile * inputfile_fit ;
  //if(usePostFit && useAllRegions)    inputfile_fit = new TFile("thetaInOut/outputTheta_merged_S1Res700Inv50_electrons_AllRegions.root");
  //if(usePostFit && useAllRegions)    inputfile_fit = new TFile("thetaInOut/outputTheta_merged_AllRegions.root");
//if(usePostFit && useAllRegions)    inputfile_fit = new TFile("outputTheta_merged_AllRegions.root");
  if(usePostFit && useAllRegions)    inputfile_fit = new TFile("outputTheta_merged_testInterRegion.root");
  else if(usePostFit)                inputfile_fit = new TFile("interRegionTest/outputTheta_interRegionOnly.root");
  //else if(usePostFit)                inputfile_fit = new TFile("outputTheta_merged_InterRegion.root");
  //else if(usePostFit)                inputfile_fit = new TFile("outputTheta_merged_CRsOnly.root");
  else                               inputfile_fit = new TFile(inputfilename);

  TFile * inputfile         = new TFile(inputfilename );

  TString dataAname         = inputdistrib+"__SingleMuA";
  TString dataBname         = inputdistrib+"__SingleMuB";
  TString dataCname         = inputdistrib+"__SingleMuC";
  TString dataDname         = inputdistrib+"__SingleMuD";

  if(useElectronChannel)
  {
    dataAname         = inputdistrib+"__SingleElA";
    dataBname         = inputdistrib+"__SingleElB";
    dataCname         = inputdistrib+"__SingleElC";
    dataDname         = inputdistrib+"__SingleElD";
  }

  TH1F * distrib__RUNA  	= (TH1F*)inputfile->Get(dataAname)->Clone();
  TH1F * distrib__RUNB  	= (TH1F*)inputfile->Get(dataBname)->Clone();
  TH1F * distrib__RUNC  	= (TH1F*)inputfile->Get(dataCname)->Clone();
  TH1F * distrib__RUND  	= (TH1F*)inputfile->Get(dataDname)->Clone();

  TH1F * histo_data  	    = (TH1F*)distrib__RUNA;
         histo_data->Add(distrib__RUNB);
         histo_data->Add(distrib__RUNC);
         histo_data->Add(distrib__RUND);


  vector<TH1F *> histo_mcSamples;
  vector<TH1F *> signalSamples;

  TString   inputdistrib_signal = inputdistrib;
  if(usePostFit)  inputdistrib  = thetainputdistrib;

  for(unsigned int imc = 0; imc < mcSample_list.size(); imc++)
  {
    if (usePostFit && mcSample_list[imc] == "QCD" && region == "TTbarregion") continue;
    if (usePostFit && mcSample_list[imc] == "DY10To50" && region == "Selectedsignalregion") continue;
    TString histo_mc_name   = inputdistrib+"__"+mcSample_list[imc];
    TH1F * histo_tmp = (TH1F*)inputfile_fit->Get(histo_mc_name)->Clone();
    cout << "Processing " << histo_mc_name << " ..." << endl;
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
  the_stack_histo->GetYaxis()->SetTitleSize(0.04);
  the_stack_histo->GetYaxis()->SetLabelSize(0.04);
  //the_stack_histo->GetYaxis()->SetTitleOffset(0.012);
  the_stack_histo->GetYaxis()->SetTitle("Events / 20 GeV");
  the_stack_histo->SetMaximum(300.);
  if(histo_data->GetMaximum() > the_stack_histo->GetMaximum() ) the_stack_histo->SetMaximum(     histo_data->GetMaximum()+0.3*     histo_data->GetMaximum());
  else                                                          the_stack_histo->SetMaximum(the_stack_histo->GetMaximum()+0.3*the_stack_histo->GetMaximum());

  if(displaySignal) {the_stack_histo->SetMinimum(1.); }
  else if(setlogy)  {the_stack_histo->SetMinimum(1.); }
  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(1.2);
  histo_data->SetLineColor(1);
  histo_data->Draw("epsame");

  if (displaySignal)
  {
    for( short unsigned int isig = 0; isig < signalSample_list.size(); isig++)
    {
        TString signalname      = inputdistrib_signal+"__"+signalSample_list[isig];
        TH1F * distrib__signal  = (TH1F*)inputfile->Get(signalname)->Clone();
        if(distrib__signal != 0)
        {
            if(isig == 0) distrib__signal->SetLineColor(kBlue+1);
            else if(isig == 1) distrib__signal->SetLineColor(kBlack);
            else if(isig == 2) distrib__signal->SetLineColor(kMagenta-5);
            else if(isig == 4) distrib__signal->SetLineColor(kOrange+1);
            else               distrib__signal->SetLineColor((int) isig);
            distrib__signal->SetLineStyle(2);
            distrib__signal->SetLineWidth(3);
            distrib__signal->Draw("hsame");
            signalSamples.push_back(distrib__signal);
            if(distrib__signal->GetMaximum() > the_stack_histo->GetMaximum() ) the_stack_histo->SetMaximum(1.3*distrib__signal->GetMaximum());
        }
    }
  }

  //--------------------------
  //MC systematic plot
  //--------------------------
  TH1F * histo_syst_MC   = (TH1F*)(the_stack_histo->GetHistogram() )->Clone();
  for(unsigned int imc=0; imc< histo_mcSamples.size(); imc++)
  {
    histo_syst_MC->Add(histo_mcSamples[imc]);
  }

  c1->cd();
  TGraphErrors *thegraph = new TGraphErrors(histo_syst_MC);
  thegraph->SetFillStyle(3005);
  thegraph->SetFillColor(1);
  thegraph->Draw("e2same");

  //-------------------
  //legend and captions
  //-------------------

  myCMS_lumi(c1);

  TString info_data = "#mu channel";

  TLatex * text2 = new TLatex(0.55,0.96, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetTextAngle(0);
  text2->SetTextFont(42);
  text2->SetX(0.45);
  text2->SetY(0.88);
  //text2->SetLineWidth(2);
  text2->SetTextFont(42);
  text2->SetTextSize(0.040);
  //    text2->SetTextSizePixels(24);// dflt=28
  text2->Draw();

  TLegend*  qw = new TLegend(.75,.60,.90,.90);

  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);

  qw->AddEntry(histo_data,         "Data" ,                "ep");
  for(unsigned int i=0; i<mcSample_list.size(); i++)
  {
    if(!usePostFit)
    {
      if( mcSample_list[i] == "TTMSDecays_central"  )               qw->AddEntry( histo_mcSamples[i],  "t#bar{t}"	,"f");
      if( mcSample_list[i] == "WExclb"              )               qw->AddEntry( histo_mcSamples[i],  "W+bjets"	,"f");
      if( mcSample_list[i] == "WExclc"              )               qw->AddEntry( histo_mcSamples[i],  "W+cjets"	,"f");
      if( mcSample_list[i] == "WExcll"              )               qw->AddEntry( histo_mcSamples[i],  "W+ljets"	,"f");
      if( mcSample_list[i] == "DYJetsToLL_M-50"     )               qw->AddEntry( histo_mcSamples[i],  "DY"		    ,"f");
      if( mcSample_list[i] == "T_s" 		        )               qw->AddEntry( histo_mcSamples[i+4],"Single top" ,"f");
      if( mcSample_list[i] == "WZ"                  )               qw->AddEntry( histo_mcSamples[i+2],"VV"	        ,"f");
      if( mcSample_list[i] == "QCD" && region != "TTbarregion" )    qw->AddEntry( histo_mcSamples[i],  "QCD"	    ,"f");
    }
    else
    {
      if( mcSample_list[i] == "TTMSDecays"  )       qw->AddEntry( histo_mcSamples[i],  "t#bar{t}"	,"f");
      if( mcSample_list[i] == "WExclb"      )       qw->AddEntry( histo_mcSamples[i],  "W+bjets"	,"f");
      if( mcSample_list[i] == "WExclc"      )       qw->AddEntry( histo_mcSamples[i],  "W+cjets"	,"f");
      if( mcSample_list[i] == "WExcll"      )       qw->AddEntry( histo_mcSamples[i],  "W+ljets"	,"f");

      if(!useMergedThetaSamples)
      {
          if( mcSample_list[i] == "DY50"    )       qw->AddEntry( histo_mcSamples[i],  "DY"		    ,"f");
          if( mcSample_list[i] == "Ts" 		)       qw->AddEntry( histo_mcSamples[i+4],"Single top" ,"f");
          if( mcSample_list[i] == "WZ"      )       qw->AddEntry( histo_mcSamples[i+2],"VV"	        ,"f");
      }
      else
      {
          if( mcSample_list[i] == "DY"          )   qw->AddEntry( histo_mcSamples[i],  "DY"		            ,"f");
          if( mcSample_list[i] == "SingleTop"   )   qw->AddEntry( histo_mcSamples[i+1],"Single top"         ,"f");
          if( mcSample_list[i] == "VV"          )   qw->AddEntry( histo_mcSamples[i],  "VV"	                ,"f");
      }
      if( mcSample_list[i] == "QCD" && region != "TTbarregion" )   qw->AddEntry( histo_mcSamples[i],  "QCD"	,"f");
    }
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
  histo_ratio_data->GetXaxis()->SetTitleOffset(1.1);
  if (variable == "DeltaPhiLJ")     histo_ratio_data->GetXaxis()->SetTitle("#Delta #phi (lept - leading jet) [GeV]");
  else if(variable == "MET")        histo_ratio_data->GetXaxis()->SetTitle("MET [GeV]");
  else if(variable == "mWT")        histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else if(variable == "ptW")        histo_ratio_data->GetXaxis()->SetTitle("p_{T}^{W} [GeV]");

  histo_ratio_data->GetXaxis()->SetLabelSize(0.04);
  histo_ratio_data->GetXaxis()->SetTitleSize(0.054);
  histo_ratio_data->GetYaxis()->SetLabelSize(0.03);
  histo_ratio_data->GetYaxis()->SetNdivisions(6);
  histo_ratio_data->GetYaxis()->SetTitleSize(0.03);
  histo_ratio_data->GetYaxis()->SetTitle("Data / MC");
  histo_ratio_data->Draw("E1X0");

  thegraph_ratio->Draw("e2same");

  TString outputname;
  TString                channel = "mujets";
  if(useElectronChannel) channel = "eljets";
  if (!usePostFit)          outputname = "PreFit_finalVersion_"+variable+"_"+channel;
  else if (useAllRegions)   outputname = "PostFit_finalVersion_AllRegions_";
  //else if (useAllRegions)   outputname = "PostFit_finalVersion_AllRegions_"+variable+"_"+channel;
  else                      outputname = "PostFit_interRegionOnly_";
  //else                      outputname = "PostFit_SR_"+variable+"_"+channel;
  //else                      outputname = "PostFit_finalVersion_InterRegionOnly_"+variable+"_"+channel;
  //else                      outputname = "PostFit_finalVersion_CRsOnly_"+variable+"_"+channel;

  if(!displaySignal)                                    outputname += "_noSignal";
  else if(signalSample_list[0].SubString("S4") == "S4") outputname += "_S4";
  else if(signalSample_list[0].SubString("S1") == "S1") outputname += "_S1";

  TString endname_root;
  TString endname_pdf;
  TString endname_png;
  if (setlogy) { endname_root = "_logY.root"; endname_pdf = "_logY.pdf"; endname_png = "_logY.png";}
  else         { endname_root = ".root";      endname_pdf = ".pdf";      endname_png = ".png"; }


  //c1->SaveAs( ("newfinalPlots_AN/"+outputname+"_"+region+endname_root).Data());
  //c1->SaveAs( ("newfinalPlots_AN/"+outputname+"_"+region+endname_pdf).Data());
/*
  if( !usePostFit && region != "TTbarregion") c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+endname_root).Data());
  else                                        c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+ttCRJetMult+endname_root).Data());
  if( !usePostFit && region != "TTbarregion") c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+endname_pdf).Data());
  else                                        c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+ttCRJetMult+endname_pdf).Data());
  if( !usePostFit && region != "TTbarregion") c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+endname_png).Data());
  else                                        c1->SaveAs( ("plots_PAS_AN/"+outputname+"_"+region+ttCRJetMult+endname_png).Data());
*/
  c1->SaveAs( ("interRegionTest/"+outputname+"_"+region+endname_root).Data());
  c1->SaveAs( ("interRegionTest/"+outputname+"_"+region+endname_pdf).Data());
  c1->SaveAs( ("interRegionTest/"+outputname+"_"+region+endname_png).Data());

}


void PlotFitResults()
{

  bool usePostFit = true;
  bool displaySignal = false;
  bool useAllRegions = false;
  bool useMergedThetaSamples = true;
  bool useElectronChannel = false;
  TString ttCRJetMult = "_2j2b";

  TString                channel = "mujets";
  if(useElectronChannel) channel = "eljets";
  TString                inputfilename = "../TreeReader/outputroot_withSyst/histo_merged.root";
  //TString                inputfilename = "../TreeReader/outputroot_withSyst/2j2b/histo_merged.root";
  if(useElectronChannel) inputfilename = "../TreeReader/outputroot_withSyst/histo_merged_electron_syst.root";


  //-------------------------
  //define list of signal samples
  vector<TString> signalSample_list;

  if( displaySignal )
  {

 //   signalSample_list.push_back("S4Inv100");
 //   signalSample_list.push_back("S4Inv500");
 //   signalSample_list.push_back("S4Inv1000");


    signalSample_list.push_back("S1Res700Inv100");
    signalSample_list.push_back("S1Res1300Inv100");
    signalSample_list.push_back("S1Res2100Inv100");

/*
    signalSample_list.push_back("S1Res500Inv10");
    signalSample_list.push_back("S1Res500Inv50");
    signalSample_list.push_back("S1Res500Inv100");
    signalSample_list.push_back("S1Res500Inv150");
    signalSample_list.push_back("S1Res500Inv200");
*/
    }


   //-------------------------
  //define list of MC samples
  vector<TString> mcSample_list;
  vector<int> colorVector;

  if(!usePostFit)
  {
      mcSample_list.push_back("TTMSDecays_central");          colorVector.push_back(kRed+1);
      mcSample_list.push_back("WExclb");                      colorVector.push_back(kGreen+4);
      //mcSample_list.push_back("WExclb");                      colorVector.push_back(kGreen-2);
      mcSample_list.push_back("WExclc");                      colorVector.push_back(kGreen);
      mcSample_list.push_back("WExcll");                      colorVector.push_back(kGreen-2);
      //mcSample_list.push_back("WExcll");                      colorVector.push_back(kGreen+4);
      mcSample_list.push_back("DYJetsToLL_M-10To50");         colorVector.push_back(kAzure-2);
      mcSample_list.push_back("DYJetsToLL_M-50");             colorVector.push_back(kAzure-2);
      mcSample_list.push_back("T_s");                         colorVector.push_back(13);
      mcSample_list.push_back("T_t");                         colorVector.push_back(13);
      mcSample_list.push_back("T_tW");                        colorVector.push_back(13);
      mcSample_list.push_back("Tbar_t");                      colorVector.push_back(13);
      mcSample_list.push_back("Tbar_tW");                     colorVector.push_back(13);
      mcSample_list.push_back("WZ");                          colorVector.push_back(kRed+2);
      mcSample_list.push_back("WW");                          colorVector.push_back(kRed+2);
      mcSample_list.push_back("ZZ");                          colorVector.push_back(kRed+2);
      mcSample_list.push_back("QCD");                         colorVector.push_back(kYellow+1);
  }
  else
  {
      mcSample_list.push_back("TTMSDecays");    colorVector.push_back(kRed+1);
      mcSample_list.push_back("WExclb");        colorVector.push_back(kGreen+4);
      mcSample_list.push_back("WExclc");        colorVector.push_back(kGreen);
      mcSample_list.push_back("WExcll");        colorVector.push_back(kGreen-2);

      if(!useMergedThetaSamples)
      {
          mcSample_list.push_back("DY10To50");      colorVector.push_back(kAzure-2);
          mcSample_list.push_back("DY50");          colorVector.push_back(kAzure-2);
          mcSample_list.push_back("Ts");            colorVector.push_back(13);
          mcSample_list.push_back("Tt");            colorVector.push_back(13);
          mcSample_list.push_back("TtW");           colorVector.push_back(13);
          mcSample_list.push_back("Tbart");         colorVector.push_back(13);
          mcSample_list.push_back("TbartW");        colorVector.push_back(13);
          mcSample_list.push_back("WZ");            colorVector.push_back(kRed+2);
          mcSample_list.push_back("WW");            colorVector.push_back(kRed+2);
          mcSample_list.push_back("ZZ");            colorVector.push_back(kRed+2);
      }
      else
      {
          mcSample_list.push_back("DY");            colorVector.push_back(kAzure-2);
          mcSample_list.push_back("SingleTop");     colorVector.push_back(13);
          mcSample_list.push_back("SingleTopW");    colorVector.push_back(13);
          mcSample_list.push_back("VV");            colorVector.push_back(kRed+2);
      }

      mcSample_list.push_back("QCD");           colorVector.push_back(kYellow+1);
  }


  //PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_mujets_ATLASRESsignalregion", "mWTmujetsATLASRESSignalregion", "../TreeReader/outputroot_withSyst/histo_merged.root","ATLASsignalregion", "mWT", useAllRegions, useMergedThetaSamples);


  //if (!usePostFit || useAllRegions) PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_"+channel+"_Selectedsignalregion", "mWT"+channel+"SelectedSignalregion", inputfilename,"Selectedsignalregion", "mWT", useAllRegions, useMergedThetaSamples, useElectronChannel, ttCRJetMult);
  //PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_"+channel+"_Selectedsignalregion", "mWT"+channel+"SelectedSignalregion", inputfilename,"Selectedsignalregion", "mWT", useAllRegions, useMergedThetaSamples, useElectronChannel, ttCRJetMult);
  PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_"+channel+"_interRegion", "mWT"+channel+"interRegion", inputfilename,"interRegion", "mWT", useAllRegions, useMergedThetaSamples, useElectronChannel, ttCRJetMult);
  //PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_"+channel+"_ttbarregion"+ttCRJetMult, "mWT"+channel+"ttbarregionHighpt", inputfilename, "TTbarregion", "mWT", useAllRegions, useMergedThetaSamples, useElectronChannel, ttCRJetMult);
  //PlotFitResults(signalSample_list, mcSample_list, colorVector, usePostFit,  displaySignal, "mWT_"+channel+"_Wregion_highpt", "mWT"+channel+"WregionHighpt", inputfilename, "Wregion", "mWT", useAllRegions, useMergedThetaSamples, useElectronChannel, ttCRJetMult);


}
