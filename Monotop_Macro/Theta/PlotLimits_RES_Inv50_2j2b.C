#include "TString.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>
#include <vector>

#include "writeLimits_RES.C"

using namespace std;

void PlotLimits()
{

  bool isExtern = true;
  //gROOT->ProcessLine(".L writeLimits.C+");
  const unsigned int n=9;   // number of mass points

  Double_t y_obs[n], y_obs_err[n], y_exp2up[n], y_expup[n], y_exp[n], y_expdn[n], y_exp2dn[n] ;
  Double_t x[n] = {500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100};

  writeLimits("RES", y_obs, y_obs_err, y_exp2up, y_expup, y_exp, y_expdn, y_exp2dn, isExtern);

  bool useATLASvsCMS = false;
  bool displayExclXsec = true;

  Double_t x_sec[n]       = {11.8  ,2.77, 0.87, 0.31, 0.11, 0.05, 0.022, 0.0097, 0.0045};
  //Double_t x_sec[n]       = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};

  TString name = "ShapeAnalysis";
/*
  Double_t y_obs_err[n]= {0.00006, 0.00005, 0.00011, 0.00029, 0.0007, 0.0020, 0.0052, 0.0169, 0.042 };
  Double_t y_obs[n]    = {0.0148,  0.0219,  0.0385,  0.0827,  0.2398, 0.560,  1.404,  4.548,  9.77  };
  Double_t y_exp2up[n] = {0.0303,  0.0531,  0.0854,  0.1794,  0.5193, 1.300,  2.750,  8.848, 21.95  };
  Double_t y_expup[n]  = {0.0182,  0.0306,  0.0511,  0.1128,  0.3171, 0.786,  1.758,  5.424, 13.25  };
  Double_t y_exp[n]    = {0.0114,  0.0189,  0.0324,  0.0703,  0.2015, 0.491,  1.114,  3.463, 8.464  };
  Double_t y_expdn[n]  = {0.0073,  0.0126,  0.0217,  0.0470,  0.1353, 0.330,  0.754,  2.312, 5.717  };
  Double_t y_exp2dn[n] = {0.0049,  0.0088,  0.0156,  0.0348,  0.1015, 0.244,  0.577,  1.728, 4.225  };
*/
/*
  Double_t y_obs[n]    = {0.0126, 0.0206, 0.0362, 0.0772, 0.2378, 0.562, 1.464, 4.418 , 10.20  };
  Double_t y_obs_err[n]= {0.00008, 0.00012, 0.00018, 0.00049, 0.0011, 0.0042, 0.0066, 0.0098, 0.015 };
  Double_t y_exp2up[n] = {0.0293, 0.0473, 0.0784, 0.1729, 0.5203, 1.270, 3.300, 8.848 , 22.75  };
  Double_t y_expup[n]  = {0.0185, 0.0296, 0.0490, 0.1078, 0.3128, 0.748, 1.918, 5.494 , 13.81  };
  Double_t y_exp[n]    = {0.0117, 0.0188, 0.0319, 0.0690, 0.1955, 0.469, 1.194, 3.463 , 8.444  };
  Double_t y_expdn[n]  = {0.0076, 0.0125, 0.0213, 0.0460, 0.1313, 0.314, 0.803, 2.302 , 5.677  };
  Double_t y_exp2dn[n] = {0.0051, 0.0088, 0.0157, 0.0342, 0.0985, 0.238, 0.599, 1.728 , 4.275  };
*/

  if(displayExclXsec)
  {
      for(unsigned short int ind = 0; ind < n; ind++)
      {
          y_obs[ind]      *= x_sec[ind];
          y_obs_err[ind]  *= x_sec[ind];
          y_exp2up[ind]   *= x_sec[ind];
          y_expup[ind]    *= x_sec[ind];
          y_exp[ind]      *= x_sec[ind];
          y_expdn[ind]    *= x_sec[ind];
          y_exp2dn[ind]   *= x_sec[ind];
      }
  }

  TGraph* xsec_th_0p1 = new TGraph(n, x, x_sec);
  xsec_th_0p1->SetLineColor(kRed);
  xsec_th_0p1->SetLineWidth(2);
  xsec_th_0p1->SetLineStyle(2);


  Double_t x_sec_0p05[n]       = {0};
  for(unsigned short int ii = 0; ii < n; ii++)
  {
    x_sec_0p05[ii] = x_sec[ii]/4. ;
  }

  TGraph* xsec_th_0p05 = new TGraph(n, x, x_sec_0p05);
  xsec_th_0p05->SetLineColor(kBlue);
  xsec_th_0p05->SetLineWidth(2);
  xsec_th_0p05->SetLineStyle(2);


  Double_t y_sigmaup[n], y_sigmadn[n], y_sigma2up[n], y_sigma2dn[n];

  for(int i = 0; i<n; i++){
     y_sigmadn[i]  = y_exp[i] - y_expdn[i];
     y_sigmaup[i]  = y_expup[i] - y_exp[i];
     y_sigma2dn[i] = y_exp[i] - y_exp2dn[i];
     y_sigma2up[i] = y_exp2up[i] - y_exp[i];
  }


  TCanvas *c1 = new TCanvas();
  if(displayExclXsec) { c1->SetLogy(1);}
  c1->cd();

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetName("");

  TGraphErrors* gr_obs = new TGraphErrors(n, x, y_obs, 0, y_obs_err);
  gr_obs->SetLineWidth(2);
  gr_obs->SetMarkerStyle(20);

  TGraphErrors* gr_exp = new TGraphErrors(n, x, y_exp, 0, 0);
  gr_exp->SetLineWidth(2);
  gr_exp->SetLineStyle(2);
  gr_exp->SetMarkerStyle(22);

  TGraphAsymmErrors* sigma2 = new TGraphAsymmErrors(n, x, y_exp, 0, 0, y_sigma2dn, y_sigma2up);
  sigma2->SetFillColor(kYellow);
  sigma2->SetLineColor(kYellow);

  TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(n, x, y_exp, 0, 0, y_sigmadn, y_sigmaup);
  sigma1->SetFillColor(kGreen);
  sigma1->SetLineColor(kGreen);

  TGraphErrors* gr_expup = new TGraphErrors(n, x, y_expup, 0, 0);
  gr_expup->SetLineColor(kRed);
  gr_expup->SetLineWidth(2);
  gr_expup->SetLineStyle(2);

  TGraphErrors* gr_expdn = new TGraphErrors(n, x, y_expdn, 0, 0);
  gr_expdn->SetLineColor(kRed);
  gr_expdn->SetLineWidth(2);
  gr_expdn->SetLineStyle(2);

  TGraphErrors* gr_exp2up = new TGraphErrors(n, x, y_exp2up, 0, 0);
  gr_exp2up->SetLineColor(kBlue);
  gr_exp2up->SetLineWidth(2);
  gr_exp2up->SetLineStyle(4);

  TGraphErrors* gr_exp2dn = new TGraphErrors(n, x, y_exp2dn, 0, 0);
  gr_exp2dn->SetLineColor(kBlue);
  gr_exp2dn->SetLineWidth(2);
  gr_exp2dn->SetLineStyle(4);

  mg->Add(sigma2);
  mg->Add(sigma1);
  mg->Draw("3AC");
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetTitleOffset(.8);
  mg->GetYaxis()->SetTitleOffset(.7);
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.055);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetXaxis()->SetTitle("m_{Res.} [GeV]");
  mg->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  if(displayExclXsec)
  {
      mg->GetYaxis()->SetTitle("#sigma(p p #rightarrow t f_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
      mg->GetYaxis()->SetTitleSize(0.04);
      mg->GetYaxis()->SetTitleOffset(1.2);
      mg->SetMinimum(0.008);
      mg->SetMaximum(x_sec[0]);
  }
  else
  {
      mg->GetYaxis()->SetTitle("signal strength");
      mg->GetYaxis()->SetRangeUser(0.,y_exp2up[n-1]);
  }

  gr_exp  ->Draw("samelp");
  gr_obs  ->Draw("samelp");

  TLine* line = new TLine(x[0], 1., x[n-1], 1.);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  if(!displayExclXsec) line->Draw("same");
  else                 {xsec_th_0p1->Draw("same");xsec_th_0p05->Draw("same");}

  TLegend *leg;
  if(!displayExclXsec) leg = new TLegend(0.3,0.65,0.55,0.85);
  else                 leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry(gr_obs  ," Observed 95\% CL limit","lpe");
  leg->AddEntry(gr_exp  ," Expected 95\% CL limit","lp");
  leg->AddEntry(sigma1  ," Expected #pm1#sigma","f");
  leg->AddEntry(sigma2  ," Expected #pm2#sigma","f");
  if(displayExclXsec) {leg->AddEntry(xsec_th_0p1  ," Theory (LO), a_{res} = 0.1","l"); leg->AddEntry(xsec_th_0p05  ," Theory (LO), a_{res} = 0.05","l");}
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 19.7 fb^{-1}, #sqrt{s} = 8 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw();

  TPaveText* text_2;
  if(!displayExclXsec) text_2 = new TPaveText(0.3, 0.45, 0.55, 0.55, "NDC");
  else                 text_2 = new TPaveText(0.6, 0.45, 0.85, 0.55, "NDC");
  text_2->AddText("Resonant model");
  text_2->AddText("m(Inv.) = 50 GeV");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(22);
  text_2->Draw();

  TString             Y_axis = "SignalStrength";
  if(displayExclXsec) Y_axis = "ExcludedXsection";

  TString externOfSyst;
  if(isExtern) externOfSyst = "_Extern";
  else         externOfSyst = "_noExtern";

  c1->SaveAs("limitPlots/finalLimits_RES"+externOfSyst+"_DM50_"+name+"_"+Y_axis+".png");
  c1->SaveAs("limitPlots/finalLimits_RES"+externOfSyst+"_DM50_"+name+"_"+Y_axis+".pdf");
  c1->SaveAs("limitPlots/finalLimits_RES"+externOfSyst+"_DM50_"+name+"_"+Y_axis+".eps");

}

void PlotLimits_RES_Inv50_2j2b()
{
    PlotLimits();
}
