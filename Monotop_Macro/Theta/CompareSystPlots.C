
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



void CompareSystPlots(TString var, TString syst, TString process){

  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0); // must be kWhite but I dunno how to do that in PyROOT
  gStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");


  TFile * inputfile_fit = new TFile("../TreeReader/outputroot_withSyst/histo_merged_woWCorr.root");

  TString histoname         = var+"__"+process;
  TString histoname_minus   = var+"__"+process+"__"+syst+"__minus";
  TString histoname_plus    = var+"__"+process+"__"+syst+"__plus";

  TH1F * hist_comp_         = (TH1F*)inputfile_fit->Get(histoname)->Clone();
  TH1F * hist_comp___minus  = (TH1F*)inputfile_fit->Get(histoname_minus)->Clone();
  TH1F * hist_comp___plus   = (TH1F*)inputfile_fit->Get(histoname_plus)->Clone();


  TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
  c1->cd();

  hist_comp_->SetTitle("");


  hist_comp_->       SetLineColor(1);
  hist_comp___minus->SetLineColor(2);
  hist_comp___plus-> SetLineColor(4);


  hist_comp_->       SetLineWidth(3);
  hist_comp___minus->SetLineWidth(3);
  hist_comp___plus-> SetLineWidth(3);


  hist_comp___plus->Draw("");
  hist_comp___minus->Draw("same");
  hist_comp_->Draw("same");


  TLegend* qw = new TLegend(.45,.60,.85,.85);
  qw->SetShadowColor(0);
  qw->SetFillColor(0);
  qw->SetLineColor(0);
  qw->AddEntry(hist_comp_,          (var+"__"+process).Data() , "l");
  qw->AddEntry(hist_comp___minus,   (var+"__"+process+"__"+syst+"__minus").Data()  , "l");
  qw->AddEntry(hist_comp___plus,    (var+"__"+process+"__"+syst+"__plus").Data()  , "l");

  qw->Draw();
  c1->SaveAs( ("plots/"+var+"__"+process+"__"+syst+".eps").Data());
}

void CompareSystPlots(){

    std::vector<TString > systlist;
    //systlist.push_back("W"                  );
    systlist.push_back("lept"               );
    systlist.push_back("trig"               );
    //systlist.push_back("PDF"              );
    systlist.push_back("PU"                 );
    systlist.push_back("toppt"              );
    systlist.push_back("btag"               );
    systlist.push_back("mistag"             );
    systlist.push_back("jes"                );
    systlist.push_back("jer"                );
    systlist.push_back("metuncls"           );
/*    systlist.push_back("btag__JES"          );
    systlist.push_back("btag__CSVLF"        );
    systlist.push_back("btag__CSVHFStats1"  );
    systlist.push_back("btag__CSVHFStats2"  );
    systlist.push_back("btag__CSVCErr1"     );
    systlist.push_back("btag__CSVCErr2"     );
    systlist.push_back("btag__CSVHF"        );
    systlist.push_back("btag__CSVLFStats1"  );
    systlist.push_back("btag__CSVLFStats2"  );
*/
    std::vector<TString > processlist;
/*    processlist.push_back("S1"                   );
    processlist.push_back("SingleMuA"            );
    processlist.push_back("SingleMuB"            );
    processlist.push_back("SingleMuC"            );
    processlist.push_back("SingleMuD"            );
    processlist.push_back("TTbar_Madgraph"       );
*/    processlist.push_back("WExclb"                );
    processlist.push_back("WExclc"                );
    processlist.push_back("WExcll"                );
    processlist.push_back("WExcl"                );
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
    processlist.push_back("QCD_A"                );
    processlist.push_back("QCD_B"                );
    processlist.push_back("QCD_C"                );
    processlist.push_back("QCD_D"                );
*/
    std::vector<TString > varlist;
    varlist.push_back("mWT_mujets_signalregion"     );
    varlist.push_back("MET_mujets_signalregion"     );
    varlist.push_back("mWT_mujets_Wregion_highpt"   );
    varlist.push_back("MET_mujets_Wregion_highpt"   );

    for(unsigned int ivar = 0; ivar < varlist.size(); ivar++)
    {
        for(unsigned int isyst = 0; isyst < systlist.size(); isyst++)
        {
            for(unsigned int iprocess = 0; iprocess < processlist.size(); iprocess++)
            {
                CompareSystPlots(varlist[ivar], systlist[isyst], processlist[iprocess]);
            }
        }
    }
}
