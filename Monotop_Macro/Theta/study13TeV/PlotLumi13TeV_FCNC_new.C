#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "THStack.h"
#include <iostream>
#include <vector>

using namespace std;



double getExactLimit(double xA, double yA, double xB, double yB)
{
    //limit = 1/a * (1 - b) because limit is when signal strength equals 1
    return ((xB-xA)/(yB-yA))*( 1 - 0.5*(yA+yB - (xA+xB)*(yB-yA)/(xB-xA) ) );
}

void PlotLumi13TeV_FCNC_new()
{

  const unsigned int n=10;   // number of mass points

  Double_t x[n]           = {100,   200,  300,  400,  500,  600,  700,  800,  900,  1000};
  Double_t x_sec[n]       = {599,   174,  7.33, 3.55, 1.74, 0.98, 0.58, 0.35, 0.22, 0.15};

     //Double_t lumi0p1_exp[n]       = {0.055,  0.139, 2.006, 3.403, 4.227, 7.11, 10.38, 15.5, 23.15, 31.3   };
     //Double_t lumi0p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     Double_t lumi1p0_exp[n]       = {0.0194, 0.052, 0.734, 1.073, 1.227, 1.98, 2.734, 4.13, 6.067, 7.413  };
     //Double_t lumi1p5_exp[n]       = {0.0163, 0.051, 0.613, 0.973, 1.057, 1.67, 2.324, 3.42, 4.947, 6.35   };
     Double_t lumi2p0_exp[n]       = {0.0153, 0.048, 0.540, 0.813, 0.944, 1.55, 2.094, 2.92, 4.317, 5.08  };
     //Double_t lumi2p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     //Double_t lumi3p0_exp[n]       = {0.0133, 0.043, 0.497, 0.733, 0.837, 1.37, 1.834, 2.39, 3.667, 4.43  };
     //Double_t lumi3p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     Double_t lumi4p0_exp[n]       = {0.0113, 0.040, 0.463, 0.693, 0.687, 1.18, 1.684, 2.05, 2.977, 3.93  };
     //Double_t lumi5p0_exp[n]       = {0.0123, 0.038, 0.423, 0.633, 0.677, 1.11, 1.464, 2.05, 3.047, 3.48  };
     //Double_t lumi6p0_exp[n]       = {0.0104, 0.040, 0.401, 0.622, 0.629, 1.09, 1.414, 1.83, 2.807, 3.21  };
     //Double_t lumi7p0_exp[n]       = {0.0100, 0.039, 0.400, 0.582, 0.627, 0.96, 1.384, 1.69, 2.757, 2.99  };
     Double_t lumi8p0_exp[n]       = {0.0100, 0.038, 0.391, 0.522, 0.569, 0.90, 1.274, 1.57, 2.687, 2.84  };
     Double_t lumi16p0_exp[n]      = {0.0090, 0.036, 0.311, 0.472, 0.499, 0.66, 0.914, 1.15, 1.947, 2.06  };

     vector<Double_t*> lumis;
     lumis.push_back(lumi1p0_exp);
     lumis.push_back(lumi2p0_exp);
     lumis.push_back(lumi4p0_exp);
     lumis.push_back(lumi8p0_exp);
     lumis.push_back(lumi16p0_exp);

     int lumiSize = lumis.size();
     Double_t myLumis[5] = {1., 2., 4., 8., 16.};
     Double_t myLimits[lumiSize];

     for(unsigned short int lumi = 0; lumi<lumis.size(); lumi++)
     {
         for(unsigned short int ns = 0; ns<n; ns++)
         {
             double mInvminus  = 100*(ns+1);
             double mInvplus   = 100*(ns+2);

             if(ns != n-1 && lumis[lumi][ns] < 1 && lumis[lumi][ns+1] > 1)
             {
                 cout << "Limit is between: " << mInvminus << " and " << mInvplus << endl;
                 cout << getExactLimit(mInvminus, lumis[lumi][ns], mInvplus, lumis[lumi][ns+1]) << endl;
                 myLimits[lumi] = getExactLimit(mInvminus, lumis[lumi][ns], mInvplus, lumis[lumi][ns+1]);
             }
         }
     }

  TCanvas *c1 = new TCanvas();
  //c1->SetLogy(1);
  c1->cd();


  TGraph* limits = new TGraph(lumis.size(), myLumis, myLimits);
  limits->SetLineColor(kBlack);
  limits->SetLineWidth(2);
  limits->SetLineStyle(1);


  limits->Draw();
  limits->SetTitle("");
  limits->GetXaxis()->SetTitleFont(62);
  limits->GetYaxis()->SetTitleFont(62);
  limits->GetXaxis()->SetTitleSize(0.045);
  limits->GetYaxis()->SetTitleSize(0.055);
  limits->GetXaxis()->SetLabelSize(0.04);
  limits->GetYaxis()->SetLabelSize(0.05);
  limits->GetXaxis()->SetTitle("luminosity (fb^{-1})");
  limits->GetXaxis()->SetTitleOffset(0.8);
  limits->GetYaxis()->SetTitle("limit on m_{#chi} [GeV]   ");
  limits->GetYaxis()->SetTitleOffset(0.8);
  limits->SetMinimum(x[0]);
  limits->SetMaximum(x[n-1]);
  limits->GetXaxis()->SetRangeUser(myLumis[0],myLumis[lumiSize-1]);

  TPaveText* text_2;
  text_2 = new TPaveText(0.5, 0.35, 0.85, 0.55, "NDC");
  text_2->AddText("Non-resonant model");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(22);
  text_2->Draw();


  TLine* line_h = new TLine(1, 617, 7.2, 617);
  line_h->SetLineWidth(2);
  line_h->SetLineColor(kRed);
  line_h->SetLineStyle(2);
  line_h->Draw("same");

  TLine* line_v = new TLine(7.2, 100, 7.2, 617);
  line_v->SetLineWidth(2);
  line_v->SetLineColor(kRed);
  line_v->SetLineStyle(2);
  line_v->Draw("same");

  TLegend *leg = new TLegend(0.55,0.70,0.85,0.85);
  leg->AddEntry(limits  ,"13 TeV, a_{res} = 0.1","l");
  leg->AddEntry(line_h ,"8 TeV current limit, a_{res} = 0.1","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("same");



  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, #sqrt{s} = 13 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw("same");

  c1->SaveAs("Comparison_lumi13TeV_FCNC.png");
  c1->SaveAs("Comparison_lumi13TeV_FCNC.pdf");
  c1->SaveAs("Comparison_lumi13TeV_FCNC.eps");

}
