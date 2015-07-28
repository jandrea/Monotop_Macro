#include "TString.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include <iostream>

using namespace std;

bool PlotStack(TString varname, TString namechan, TString selection, bool setlogy,
     vector<TString> dataSample_list,
     vector<TString> channel_list, vector<TString> mcSample_list, vector<TString > signalSample_list, vector<int> colorVector, short int QCDCorr, bool sumChannels, bool useElectronChannel)
{

  bool QCDregion = false;
  TString channel = "";
  if(!sumChannels ) channel = namechan;

  TString filename;

  //if( QCDCorr == 0) filename = "../TreeReader/outputroot_QCDcorr_iso0p6/histo_merged.root";
  if( QCDCorr == 0) filename = "../TreeReader/outputroot_QCDcorr_iso0p5/backup_final_muons/histo_merged.root";
  //if( QCDCorr == 0) filename = "../TreeReader/outputroot_QCDcorr/histo_merged.root";
  if( QCDCorr == 1) filename = "../TreeReader/outputroot_woQCDcorr/backup_final_muons/histo_merged.root";
  if( QCDCorr == 2) filename = "../TreeReader/outputroot_withSyst_beforeJuly15/histo_merged_nosyst.root";
  if( QCDCorr == 3 && !useElectronChannel) filename = "../TreeReader/outputroot_withSyst/histo_merged.root";
  else if( QCDCorr == 3)                   filename = "../TreeReader/outputroot_withSyst_beforeJuly15/histo_merged_electron_syst.root";

  Int_t stati=0;
  Bool_t  fit=1;

  //-----------------------------------
  //define the canvas
  //-----------------------------------

  // For the canvas:
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetCanvasDefH(600); //Height of canvas
  gStyle->SetCanvasDefW(600); //Width of canvas
  gStyle->SetCanvasDefX(0);   //POsition on screen
  gStyle->SetCanvasDefY(0);


  // For the Pad:
  gStyle->SetPadBorderMode(0);
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

  // For the histo:
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineStyle(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetEndErrorSize(2);


  //For the fit/function:
  gStyle->SetOptFit(1011);
  gStyle->SetFitFormat("5.4g");
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);

  //For the date:
  gStyle->SetOptDate(0);

  // For the statistics box:
  gStyle->SetOptFile(0);
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  gStyle->SetStatColor(0); // kWhite
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatTextColor(1);
  gStyle->SetStatFormat("6.4g");
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.15);

  // Margins:
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);

  // For the Global title:
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.06, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.25);

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

  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  c1->SetBottomMargin(0.3);
  c1->SetLogy(setlogy);
  c1->cd();

  TFile * filechannel = new TFile(filename);

  TH1D *histo_data = 0;
  vector<TH1D *> histo_mcSamples;
  vector<TH1D *> histo_mcSignal;
  //---------------------------------------------
  //retrieve and add histograms
  //---------------------------------------------
  int niter_data = 0;
  int niter_chan  = 0;
  for(unsigned int ichan = 0; ichan < channel_list.size() ; ichan++)
  {
      if( !sumChannels  &&  channel_list[ichan] != channel ) continue;


      //--------------------
      //loop over datasamples
      //--------------------
      for(unsigned int idatasampl = 0; idatasampl< dataSample_list.size() ; idatasampl++)
      {
          TString histo_Data_name   = varname+"_"+channel_list[ichan]+"_"+selection+ "__"+dataSample_list[idatasampl];
          if( niter_data == 0 && niter_chan == 0)
          {
              histo_data = (TH1D*)filechannel->Get(histo_Data_name)->Clone();
              if( histo_data== 0)  cout << "  no existing histo data with name  " << histo_Data_name  << endl;
          }
          else if(niter_data!=0)
          {
              TH1D * histo_tmp = (TH1D*)filechannel->Get(histo_Data_name)->Clone();
              if(histo_tmp == 0)  cout << "  no existing histo data with name  " << histo_Data_name  << endl;
              histo_data->Add(histo_data, histo_tmp);
              cout << "DataBin(1)= " << histo_data->GetBinContent(1) << endl;
          }
          niter_data++;
      }

      if(QCDregion) histo_data->GetXaxis()->SetRangeUser(0,80);

      //--------------------
      //loop over MC samples
      //--------------------
      for(unsigned int imc = 0; imc < mcSample_list.size(); imc++)
      {
          TString histo_mc_name   = varname+"_"+channel_list[ichan]+"_"+selection+ "__"+ mcSample_list[imc];
          TH1D * histo_tmp = (TH1D*)filechannel->Get(histo_mc_name)->Clone();
          int numchan = -1;

          if(channel_list[ichan] == "mujets" || channel_list[ichan] == "eljets") numchan = 0;

          if(histo_tmp == 0)  cout << "  no existing histo with name  " << histo_mc_name << endl;
          if(sumChannels)
          {
              if(niter_chan == 0)
              {
                  histo_mcSamples.push_back(histo_tmp);
              }
              else
              {
                  histo_mcSamples[imc]->Add(histo_mcSamples[imc], histo_tmp);
              }
          }
          else
          {
              histo_mcSamples.push_back(histo_tmp);
          }
      }

      niter_chan++;

  }// loop over the channels

  THStack *the_stack_histo= new THStack();
  for(unsigned int imc = 0; imc< histo_mcSamples.size(); imc++)
  {
    histo_mcSamples[imc]->SetFillStyle(1001);
    histo_mcSamples[imc]->SetFillColor(colorVector[imc]);
    histo_mcSamples[imc]->SetLineColor(colorVector[imc]);
    if(QCDregion) histo_mcSamples[imc]->GetXaxis()->SetRangeUser(0,80);
    if(imc < histo_mcSamples.size() && colorVector[imc] != colorVector[imc+1] )  histo_mcSamples[imc]->SetLineColor(1);
    if(imc ==  histo_mcSamples.size()) histo_mcSamples[imc]->SetLineColor(1);
    the_stack_histo->Add(histo_mcSamples[imc]);
  }

  the_stack_histo->Draw("h");
  if(QCDregion) the_stack_histo->GetXaxis()->SetRangeUser(0,80);
  the_stack_histo->GetXaxis()->SetLabelSize(0.0);
  the_stack_histo->GetYaxis()->SetTitle("Entries");
  the_stack_histo->GetYaxis()->SetTitleSize(0.05);
  the_stack_histo->GetYaxis()->SetTitleOffset(1.5);
  if(histo_data->GetMaximum() > the_stack_histo->GetMaximum() ) the_stack_histo->SetMaximum(histo_data->GetMaximum()+0.3*histo_data->GetMaximum());
  else the_stack_histo->SetMaximum(the_stack_histo->GetMaximum()+0.3*the_stack_histo->GetMaximum());

  if(setlogy) the_stack_histo->SetMinimum(1);


  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(1.2);
  histo_data->SetLineColor(1);
  histo_data->Draw("epsame");


  //--------------------------
  //Add signal plots
  //--------------------------
  for(unsigned int isign = 0; isign < signalSample_list.size(); isign++)
  {
    TString histo_signal_name   = varname+"_"+namechan+"_"+selection+ "__"+ signalSample_list[isign];
    TH1D * histo_signal = (TH1D*)filechannel->Get(histo_signal_name)->Clone();
    histo_mcSignal.push_back(histo_signal);
    if(isign == 0)
    {
        histo_signal->SetLineColor(kRed);
    }
    else if(isign == 1)
    {
        histo_signal->SetLineColor(kBlue);
    }
    else if(isign == 2)
    {
        histo_signal->SetLineColor(kGreen);
    }
    else histo_signal->SetLineColor(kOrange);

    histo_signal->SetLineStyle(2);
    histo_signal->SetLineWidth(2);
    histo_signal->Draw("hsame");
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
  latex->DrawLatex(0.5, 0.95, "CMS Preliminary");

  TLatex *latex2 = new TLatex();
  latex2->SetNDC();
  latex2->SetTextSize(0.04);
  latex2->SetTextAlign(31);
  latex2->DrawLatex(0.96, 0.94, "#sqrt{s} = 8 TeV, L = 19.7 fb^{-1}");

  TString info_data;
  if (namechan=="mujets") info_data = "#mu channel";
  if (namechan=="eljets") info_data = "e- channel";

  TLatex * text2 = new TLatex(0.45,0.98, info_data);
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.21);
  text2->SetY(0.90);
  text2->SetTextFont(42);
  text2->SetTextSize(0.050687);
  text2->Draw();


  TLegend* qw = 0;
  if(varname == "DeltaPhiLJ" && selection != "Selectedsignalregion")  qw = new TLegend(.50,.60,.65,.90);
  else                                                                qw = new TLegend(.80,.60,.95,.90);
  if(QCDregion)                                                       qw = new TLegend(.80,.45,.95,.70);

  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);

  qw->AddEntry(histo_data,         "Data" ,                "ep");

  if(QCDCorr == 0)
  {
        for(unsigned int i=0; i<mcSample_list.size(); i++)
        {
            if(mcSample_list[i] == "NTuple_53X_TTJetsMadgraphZ2"		 ) qw->AddEntry( histo_mcSamples[0],	"t#bar{t}"	 ,"f");
            if(mcSample_list[i] == "NTuple_53X_T_s-channel" 			 ) qw->AddEntry( histo_mcSamples[6],	"single top" ,"f");
            if(mcSample_list[i] == "NTuple_53X_DYJetsToLL_M-50" 	     ) qw->AddEntry( histo_mcSamples[3],	"DY"		 ,"f");
            if(mcSample_list[i] == "NTuple_53X_WJetsToLNu"           	 ) qw->AddEntry( histo_mcSamples[1],	"W+jets"	 ,"f");
            if(mcSample_list[i] == "NTuple_53X_WWJetsIncl"             	 ) qw->AddEntry( histo_mcSamples[11],	"VV"	     ,"f");
            if(mcSample_list[i] == "QCD_Pt-120to170_MuEnrichedPt5"       ) qw->AddEntry( histo_mcSamples[16],	"QCD"	     ,"f");
        }
  }
  else if (QCDCorr == 1)
  {
        for(unsigned int i=0; i<mcSample_list.size(); i++)
        {
            if(mcSample_list[i] == "QCD_Pt-120to170"            ) qw->AddEntry( histo_mcSamples[20],	"QCD"	    ,"f");
            if(mcSample_list[i] == "WW"             		    ) qw->AddEntry( histo_mcSamples[13],	"VV"	    ,"f");
            if(mcSample_list[i] == "T_s" 			            ) qw->AddEntry( histo_mcSamples[10],    "single top","f");
            if(mcSample_list[i] == "DYJetsToLL_M-50" 	        ) qw->AddEntry( histo_mcSamples[5],	    "DY"		,"f");
            if(mcSample_list[i] == "WExcll"           	  	    ) qw->AddEntry( histo_mcSamples[3],	    "W+ljets"	,"f");
            if(mcSample_list[i] == "WExclc"           	  	    ) qw->AddEntry( histo_mcSamples[2],	    "W+cjets"	,"f");
            if(mcSample_list[i] == "WExclb"           	  	    ) qw->AddEntry( histo_mcSamples[1],	    "W+bjets"	,"f");
            if(mcSample_list[i] == "TTMSDecays_central"         ) qw->AddEntry( histo_mcSamples[0],	    "t#bar{t}"	,"f");
        }
  }
  else if (QCDCorr == 2 || QCDCorr == 3)
  {
        for(unsigned int i=0; i<mcSample_list.size(); i++)
        {
            if(mcSample_list[i] == "QCD"                        ) qw->AddEntry( histo_mcSamples[14],	"QCD"	    ,"f");
            if(mcSample_list[i] == "WW"             		    ) qw->AddEntry( histo_mcSamples[13],	"VV"	    ,"f");
            if(mcSample_list[i] == "T_s" 			            ) qw->AddEntry( histo_mcSamples[10],    "single top","f");
            if(mcSample_list[i] == "DYJetsToLL_M-50" 	        ) qw->AddEntry( histo_mcSamples[5],	    "DY"		,"f");
            if(mcSample_list[i] == "WExcll"           	  	    ) qw->AddEntry( histo_mcSamples[3],	    "W+ljets"	,"f");
            if(mcSample_list[i] == "WExclc"           	  	    ) qw->AddEntry( histo_mcSamples[2],	    "W+cjets"	,"f");
            if(mcSample_list[i] == "WExclb"           	  	    ) qw->AddEntry( histo_mcSamples[1],	    "W+bjets"	,"f");
            if(mcSample_list[i] == "TTMSDecays_central"         ) qw->AddEntry( histo_mcSamples[0],	    "t#bar{t}"	,"f");
        }
  }


  if (QCDCorr == 1 || QCDCorr == 2 || QCDCorr == 3)
  {
        for(unsigned int i=0; i<histo_mcSignal.size(); i++)
        {
            qw->AddEntry( histo_mcSignal[i],	    signalSample_list[i]	,"l");
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
  for(int i=0; i<thegraph_tmp->GetN(); i++)
  {
    if(theY[i]!=0) theErrorY[i] = theErrorY[i]/theY[i];
    theY[i]=1;
  }


  TGraphErrors *thegraph_ratio = new TGraphErrors(thegraph_tmp->GetN(), thegraph_tmp->GetX() , theY ,  thegraph_tmp->GetEX(),     thegraph_tmp->GetEY() );
  thegraph_ratio->SetFillStyle(3005);
  thegraph_ratio->SetFillColor(1);


  if(varname == "InvM_ll") 	        histo_ratio_data->GetXaxis()->SetTitle("M_{ll} [GeV/c^{-1}]");
  else if(varname == "NJet")        histo_ratio_data->GetXaxis()->SetTitle("jet mult.");
  else if(varname == "NBJet")       histo_ratio_data->GetXaxis()->SetTitle("b-tagged jet mult.");
  else if(varname == "MET")	        histo_ratio_data->GetXaxis()->SetTitle("MET [GeV]");
  else if(varname == "NVtx")	    histo_ratio_data->GetXaxis()->SetTitle("NVtx");
  else if(varname == "mWT")	        histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else if(varname == "mWT_full")	histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else if(varname == "mWTplusMET")  histo_ratio_data->GetXaxis()->SetTitle("m_{T}^{W} + missing E_{T} [GeV]");
  else if(varname == "mW")	        histo_ratio_data->GetXaxis()->SetTitle("m(W) [GeV]");
  else if(varname == "ptW")	        histo_ratio_data->GetXaxis()->SetTitle("p_{T}(W) [GeV]");
  else if(varname == "etaW")	    histo_ratio_data->GetXaxis()->SetTitle("#eta(W) ");
  else if(varname == "JetPt")	    histo_ratio_data->GetXaxis()->SetTitle("p_{T}(jets) [GeV]");
  else if(varname == "DeltaPhiLJ")	histo_ratio_data->GetXaxis()->SetTitle("#Delta#phi (lep. jet)");
  else if(varname == "JetEta")	    histo_ratio_data->GetXaxis()->SetTitle("#eta(jets)");
  else if(varname == "LeptPt")	    histo_ratio_data->GetXaxis()->SetTitle("p_{T}(lepton) [GeV]");
  else if(varname == "LeptEta")	    histo_ratio_data->GetXaxis()->SetTitle("#eta(lepton)");
  else if(varname == "LeptIso")	    histo_ratio_data->GetXaxis()->SetTitle("iso(lepton)");
  else if(varname == "BJetCSV")     histo_ratio_data->GetXaxis()->SetTitle("CSV discriminator");
  else if(varname == "LeadJetBtagDiscr")     histo_ratio_data->GetXaxis()->SetTitle("Btag discriminator(lead jet)");


  histo_ratio_data->SetMinimum(0.5);
  histo_ratio_data->SetMaximum(1.5);
  histo_ratio_data->SetTitle(0);
  histo_ratio_data->SetStats(0);
  histo_ratio_data->GetXaxis()->SetTitleOffset(1.0);
  histo_ratio_data->GetXaxis()->SetLabelSize(0.04);
  histo_ratio_data->GetXaxis()->SetTitleSize(0.05);
  histo_ratio_data->GetYaxis()->SetTitle("data/MC");
  histo_ratio_data->GetYaxis()->SetLabelSize(0.04);
  histo_ratio_data->GetYaxis()->SetNdivisions(6);
  histo_ratio_data->GetYaxis()->SetTitleSize(0.045);
  histo_ratio_data->GetYaxis()->SetTitleOffset(1.1);
  histo_ratio_data->Draw("E1X0");

  thegraph_ratio->Draw("e2same");

  TString end_name_png=".png", end_name_pdf=".pdf", end_name_C=".C", end_name_root=".root", end_name_ps=".ps", end_name_eps=".eps";

  TString outputname;
  //if(QCDCorr == 0)      outputname = "plots_QCDcorr/"    +varname+"_"+namechan+"_"+selection;
  if(QCDCorr == 0)      outputname = "plots_AN/isoInv0p5/new"    +varname+"_"+namechan+"_"+selection;
  //if(QCDCorr == 0)      outputname = "plots_QCDcorr_iso0p6/"    +varname+"_"+namechan+"_"+selection;
  else if(QCDCorr == 1) outputname = "plots_AN/withoutQCDCorr/"  +varname+"_"+namechan+"_"+selection;
  else if(QCDCorr == 2) outputname = "plots_AN/withQCDCorr/new_"+varname+"_"+namechan+"_"+selection;
  else if(QCDCorr == 3 && !useElectronChannel) outputname = "plots_AN_final/PreFit_"   +varname+"_"+namechan+"_"+selection;
  else if(QCDCorr == 3)                        outputname = "plots_electrons/PreFit_"   +varname+"_"+namechan+"_"+selection;


  if     (QCDregion && QCDCorr == 1) outputname = "plots_testQCD/beforeQCDCorr"+varname+"_"+namechan+"_"+selection;
  else if(QCDregion && QCDCorr == 2) outputname = "plots_testQCD/afterQCDCorr" +varname+"_"+namechan+"_"+selection;
  else if(QCDregion && QCDCorr == 3) outputname = "plots_testQCD/afterQCDCorr" +varname+"_"+namechan+"_"+selection;
  else if(QCDregion)                 outputname = "plots_testQCD/"             +varname+"_"+namechan+"_"+selection;

  if(setlogy) c1->SaveAs((outputname+"_Logy").Data());
  else
  {
              c1->SaveAs((outputname+end_name_pdf).Data());
              c1->SaveAs((outputname+end_name_png).Data());
              c1->SaveAs((outputname+end_name_ps).Data());
              c1->SaveAs((outputname+end_name_eps).Data());
  //            c1->SaveAs((outputname+end_name_C).Data());
              c1->SaveAs((outputname+end_name_root).Data());
  }

}
