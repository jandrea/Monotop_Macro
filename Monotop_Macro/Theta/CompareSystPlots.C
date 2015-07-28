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



void CompareSystPlots(TString var, TString syst, TString process, TString variable, TString channel, TString region)
{
  bool isLogY = false;

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");


  TFile * inputfile_fit = new TFile("../TreeReader/outputroot_withSyst/histo_merged.root");

  TString process_tmp = process;
  if(process == "TTMSDecays") process_tmp = "TTMSDecays_central";
  TString histoname         = var+"__"+process_tmp;
  TString histoname_minus;
  TString histoname_plus;

  cout << endl;
  cout << "Variable = " << variable << " | Process = " << process << " | Syst = " << syst << endl;

  if( syst != "matching" && syst != "scale") histoname_minus   = var+"__"+process_tmp+"__"+syst+"__minus";
  else if(process == "TTMSDecays")           histoname_minus   = var+"__"+process+"_"+syst+"down";
  else {cout << "ERROR: process '" << process << "' don't have syst '" << syst << "'!" << endl; return;}
  if( syst != "matching" && syst != "scale") histoname_plus    = var+"__"+process_tmp+"__"+syst+"__plus";
  else                                       histoname_plus    = var+"__"+process+"_"+syst+"up";

  TH1F * hist_comp_         = (TH1F*)inputfile_fit->Get(histoname)->Clone();
  TH1F * hist_comp___minus  = (TH1F*)inputfile_fit->Get(histoname_minus)->Clone();
  TH1F * hist_comp___plus   = (TH1F*)inputfile_fit->Get(histoname_plus)->Clone();


  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  if(isLogY) c1->SetLogy(1);
  else       c1->SetLogy(0);
  c1->cd();

  hist_comp___plus->SetTitle("");


  hist_comp_->       SetLineColor(1);
  hist_comp___minus->SetLineColor(2);
  hist_comp___plus-> SetLineColor(4);


  hist_comp_->       SetLineWidth(3);
  hist_comp___minus->SetLineWidth(3);
  hist_comp___plus-> SetLineWidth(3);

  if( variable == "mWT") hist_comp___plus->GetXaxis()->SetTitle("m_{T}^{W} [GeV]");
  else                   hist_comp___plus->GetXaxis()->SetTitle((var+" [GeV]").Data());
  hist_comp___plus->GetYaxis()->SetTitle("Normalized events");
  hist_comp___plus->GetYaxis()->SetTitleOffset(1.3);
  if(!isLogY) hist_comp___plus->SetMinimum(0);
  hist_comp___plus->DrawNormalized("");
  hist_comp___minus->DrawNormalized("same");
  hist_comp_->DrawNormalized("same");

  TLatex * text2 = new TLatex(0.60,0.88, (region+"    "+channel).Data());
  text2->SetNDC();
  text2->SetTextAlign(13);
  text2->SetX(0.45);
  text2->SetY(0.88);
  text2->SetTextFont(42);
  text2->SetTextSize(0.046);
  text2->Draw();


  TLegend* qw = new TLegend(.55,.60,.85,.80);
  //TLegend* qw = new TLegend(.45,.55,.85,.80);
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  qw->AddEntry(hist_comp___plus,    (process+"__"+syst+"__plus").Data()  , "l");
  qw->AddEntry(hist_comp_,          (process).Data() , "l");
  qw->AddEntry(hist_comp___minus,   (process+"__"+syst+"__minus").Data()  , "l");
  //qw->AddEntry(hist_comp___plus,    (var+"__"+process+"__"+syst+"__plus").Data()  , "l");
  //qw->AddEntry(hist_comp_,          (var+"__"+process).Data() , "l");
  //qw->AddEntry(hist_comp___minus,   (var+"__"+process+"__"+syst+"__minus").Data()  , "l");

  qw->Draw();
  if(isLogY) c1->SaveAs( ("plots_compareSyst/"+var+"__"+process+"__"+syst+"_logY.png").Data());
  else       c1->SaveAs( ("plots_compareSyst/"+var+"__"+process+"__"+syst+".png").Data());
  if(isLogY) c1->SaveAs( ("plots_compareSyst/"+var+"__"+process+"__"+syst+"_logY.pdf").Data());
  else       c1->SaveAs( ("plots_compareSyst/"+var+"__"+process+"__"+syst+".pdf").Data());
}

void CompareSystPlots(){

    std::vector<TString > systlist;
    systlist.push_back("lept"               );
    systlist.push_back("trig"               );
    systlist.push_back("PDF"                );
    systlist.push_back("PU"                 );
    systlist.push_back("toppt"              );
    systlist.push_back("btag"               );
    systlist.push_back("mistag"             );
    systlist.push_back("jes"                );
    systlist.push_back("jer"                );
    systlist.push_back("metuncls"           );
    systlist.push_back("matching"           );
    systlist.push_back("scale"              );

    std::vector<TString > processlist;
    //processlist.push_back("S1Res900Inv100"       );
    //processlist.push_back("S4Inv500"             );
    //processlist.push_back("DATA"                 );
    processlist.push_back("TTMSDecays"           );
    //processlist.push_back("WExclb"                );
    //processlist.push_back("WExclc"                );
    //processlist.push_back("WExcll"                );
    //processlist.push_back("WExcl"                );
/*    processlist.push_back("DYJetsToLL_M-10To50"  );
    processlist.push_back("DYJetsToLL_M-50"      );
    processlist.push_back("T_s"                  );
    processlist.push_back("T_t"                  );
    processlist.push_back("T_tW"                 );
    processlist.push_back("Tbar_t"               );
    processlist.push_back("Tbar_tW"              );
    processlist.push_back("WZ"                   );
    processlist.push_back("WW"                   );
    processlist.push_back("ZZ"                   );
    processlist.push_back("QCD"                  );
*/
    std::vector<TString > varlist;
    //varlist.push_back("mWT_mujets_Selectedsignalregion"     );
    //varlist.push_back("NVtx_mujets_Wregion_highpt"   );
    varlist.push_back("mWT_mujets_ttbarregion_2j2b"   );

    TString variable, channel, region;

    for(unsigned int ivar = 0; ivar < varlist.size(); ivar++)
    {
        for(unsigned int isyst = 0; isyst < systlist.size(); isyst++)
        {
            for(unsigned int iprocess = 0; iprocess < processlist.size(); iprocess++)
            {
                if(varlist[ivar] == "mWT_mujets_Selectedsignalregion")    { variable = "mWT"; channel = "#mu channel"; region = "Signal Region";}
                else if(varlist[ivar] == "mWT_mujets_Wregion_highpt")     { variable = "mWT"; channel = "#mu channel"; region = "W-enriched CR";}
                else if(varlist[ivar] == "mWT_mujets_ttbarregion_2j2b")   { variable = "mWT"; channel = "#mu channel"; region = "t#bar{t}-enriched CR";}
                CompareSystPlots(varlist[ivar], systlist[isyst], processlist[iprocess], variable, channel, region);
            }
        }
    }
}
