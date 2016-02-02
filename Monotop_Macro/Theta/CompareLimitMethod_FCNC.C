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

#include "writeLimits_FCNC.C"

using namespace std;

void CompareLimitMethod()
{

  short int option = 0;
  //gROOT->ProcessLine(".L writeLimits.C+");
  const unsigned int n=10;   // number of mass points

  Double_t y_obs[n], y_obs_err[n], y_exp2up[n], SA_expup[n], SA_exp[n], SA_expdn[n], y_exp2dn[n] ;
  Double_t x[n] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

  writeLimits("FCNC", y_obs, y_obs_err, y_exp2up, SA_expup, SA_exp, SA_expdn, y_exp2dn, option);

  bool useATLASvsCMS   = false;
  bool displayExclXsec = true;

  Double_t x_sec[n]       = {28.0  ,8.1  ,2.94, 1.22, 0.57, 0.29, 0.15, 0.08, 0.05, 0.03};

  TString name = "CCvsSA";
/*
     Double_t y_obs[n]    = {0.0641  , 0.173  , 0.245  , 0.418  , 0.739  , 1.129 , 1.553 , 2.439 , 4.31 , 6.02  };
     Double_t y_obs_err[n]= {0.0003  , 0.0008 , 0.0009 , 0.0013 , 0.0026 , 0.0035, 0.0048, 0.0071, 0.015, 0.022 };
     Double_t y_exp2up[n] = {0.342   , 0.712  , 0.877  , 1.145  , 2.13   , 2.81  , 3.820 , 5.475 , 10.0 , 15.16 };
     Double_t SA_expup[n]  = {0.148   , 0.33   , 0.46   , 0.631  , 1.14   , 1.60  , 2.141 , 3.311 , 5.614, 8.22  };
     Double_t SA_exp[n]    = {0.081   , 0.185  , 0.267  , 0.393  , 0.674  , 0.975 , 1.342 , 2.083 , 3.435, 4.94  };
     Double_t SA_expdn[n]  = {0.0485  , 0.113  , 0.168  , 0.254  , 0.429  , 0.637 , 0.889 , 1.372 , 2.241, 3.24  };
     Double_t y_exp2dn[n] = {0.0317  , 0.076  , 0.116  , 0.178  , 0.298  , 0.451 , 0.643 , 0.994 , 1.621, 2.27  };
*/

     Double_t CC_CMS_expup[n]   = {0.345   , 0.952  , 1.953  , 3.723  , 6.556  , 12.14 , 20.61 , 35.96 , 59.35, 96.43 };
     Double_t CC_CMS_exp[n]     = {0.207   , 0.579  , 1.162  , 2.225  , 4.013  , 7.371 , 12.11 , 21.70 , 35.54, 59.82 };
     Double_t CC_CMS_expdn[n]   = {0.126   , 0.361  , 0.696  , 1.363  , 2.453  , 4.401 , 7.433 , 13.67 , 21.46, 35.49 };

     Double_t CC_ATLAS_expup[n] = {0.400   , 0.983  , 1.219  , 2.049  , 3.338  , 5.536 , 7.704 , 11.89 , 19.25, 30.07 };
     Double_t CC_ATLAS_exp[n]   = {0.254   , 0.590  , 0.728  , 1.251  , 2.036  , 3.402 , 4.723 , 7.192 , 11.75, 18.63 };
     Double_t CC_ATLAS_expdn[n] = {0.147   , 0.372  , 0.451  , 0.794  , 1.256  , 2.111 , 2.976 , 4.566 , 7.341, 11.72 };

/*
     Double_t y_obs[n]    = {0.043   , 0.108  , 0.162  , 0.372  , 0.402  , 0.646 , 1.245 , 1.814 , 3.66 , 3.58  };
     Double_t y_obs_err[n]= {0.0003  , 0.0008 , 0.0013 , 0.0031 , 0.0026 , 0.0037, 0.0060, 0.0114, 0.023, 0.022};
     Double_t y_exp2up[n] = {0.332   , 1.112  , 1.23   , 3.108  , 2.64   , 3.97  , 5.770 , 9.295 , 16.3 , 15.73 };
     Double_t y_expup[n]  = {0.172   , 0.44   , 0.57   , 1.108  , 1.15   , 1.69  , 2.475 , 4.531 , 7.224, 7.23  };
     Double_t y_exp[n]    = {0.068   , 0.156  , 0.247  , 0.444  , 0.53   , 0.719 , 1.130 , 2.083 , 3.115, 3.166 };
     Double_t y_expdn[n]  = {0.0317  , 0.073  , 0.127  , 0.225  , 0.275  , 0.404 , 0.641 , 1.125 , 1.695, 1.8   };
     Double_t y_exp2dn[n] = {0.0193  , 0.044  , 0.076  , 0.144  , 0.172  , 0.267 , 0.438 , 0.724 , 1.131, 1.27  };
*/



/*
     Double_t y_obs[n]    = {0.0697  , 0.1990 , 0.2098 , 0.373  , 0.423  , 0.7118, 1.102 , 1.860 , 2.87 , 4.07  };
     Double_t y_obs_err[n]= {0.001086, 0.00265, 0.00230, 0.00283, 0.00268, 0.0048, 0.0070, 0.0148, 0.020, 0.0286};
     Double_t y_exp2up[n] = {0.454   , 0.828  , 1.34   , 1.718  , 2.190  , 3.655 , 5.068 , 8.920 , 10.3 , 12.63 };
     Double_t y_expup[n]  = {0.2330  , 0.34   , 0.63   , 0.880  , 1.162  , 2.132 , 3.019 , 4.241 , 6.224, 7.43  };
     Double_t y_exp[n]    = {0.0985  , 0.156  , 0.312  , 0.454  , 0.661  , 1.222 , 1.782 , 2.483 , 3.755, 4.466 };
     Double_t y_expdn[n]  = {0.0454  , 0.079  , 0.166  , 0.258  , 0.373  , 0.718 , 1.039 , 1.501 , 2.285, 2.7   };
     Double_t y_exp2dn[n] = {0.0312  , 0.047  , 0.0923 , 0.132  , 0.223  , 0.428 , 0.596 , 1.008 , 1.431, 1.75  };
*/

/*

     // Safety limits
     Double_t y_obs[n]    = {0.0503  , 0.1346 , 0.1698 , 0.322  , 0.379  , 0.568 , 1.028 , 1.530 , 2.287, 3.50  };
     Double_t y_obs_err[n]= {0.000796, 0.00189, 0.00168, 0.00287, 0.00268, 0.0049, 0.0073, 0.0122, 0.015, 0.0266};
     Double_t y_exp2up[n] = {0.399   , 0.76   , 1.63   , 2.128  , 2.190  , 3.805 , 5.218 , 7.080 , 12.93, 15.37 };
     Double_t y_expup[n]  = {0.1920  , 0.36   , 0.73   , 1.007  , 1.252  , 2.163 , 3.199 , 4.281 , 8.294, 9.72  };
     Double_t y_exp[n]    = {0.0855  , 0.152  , 0.336  , 0.521  , 0.665  , 1.172 , 1.769 , 2.438 , 4.430, 5.485 };
     Double_t y_expdn[n]  = {0.0422  , 0.078  , 0.167  , 0.269  , 0.367  , 0.653 , 0.989 , 1.461 , 2.515, 3.189 };
     Double_t y_exp2dn[n] = {0.0222  , 0.044  , 0.0833 , 0.161  , 0.221  , 0.366 , 0.629 , 0.908 , 1.391, 2.006 };
*/

/*
     // With all the 4 distr in SR and mTW in CRs
     Double_t y_obs[n]    = {0.0353  , 0.0846 , 0.1006 , 0.118  , 0.114  , 0.112 , 0.177 , 0.247 , 0.422, 0.61  };
     Double_t y_obs_err[n]= {0.000426, 0.00115, 0.00083, 0.00087, 0.00097, 0.0010, 0.0018, 0.0025, 0.005, 0.0077};
     Double_t y_exp2up[n] = {0.209   , 0.43   , 0.464  , 0.497  , 0.512  , 0.597 , 0.781 , 0.922 , 1.81 , 2.56  };
     Double_t y_expup[n]  = {0.0946  , 0.23   , 0.27   , 0.289  , 0.312  , 0.376 , 0.455 , 0.575 , 0.99 , 1.49  };
     Double_t y_exp[n]    = {0.0447  , 0.111  , 0.143  , 0.172  , 0.180  , 0.211 , 0.259 , 0.340 , 0.596, 0.883 };
     Double_t y_expdn[n]  = {0.0224  , 0.057  , 0.077  , 0.100  , 0.102  , 0.124 , 0.153 , 0.213 , 0.342, 0.532 };
     Double_t y_exp2dn[n] = {0.0125  , 0.033  , 0.0466 , 0.055  , 0.0611 , 0.070 , 0.096 , 0.125 , 0.203, 0.321 };
*/

  if(displayExclXsec)
  {
      for(unsigned short int ind = 0; ind < n; ind++)
      {
          y_obs[ind]      *= x_sec[ind];
          y_obs_err[ind]  *= x_sec[ind];

          SA_expup[ind]   *= x_sec[ind];
          SA_exp[ind]     *= x_sec[ind];
          SA_expdn[ind]   *= x_sec[ind];

          CC_CMS_expup[ind]   *= x_sec[ind];
          CC_CMS_exp[ind]     *= x_sec[ind];
          CC_CMS_expdn[ind]   *= x_sec[ind];

          CC_ATLAS_expup[ind] *= x_sec[ind];
          CC_ATLAS_exp[ind]   *= x_sec[ind];
          CC_ATLAS_expdn[ind] *= x_sec[ind];
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

  Double_t SA_sigmaup[n], SA_sigmadn[n], CC_CMS_sigmaup[n], CC_CMS_sigmadn[n], CC_ATLAS_sigmaup[n], CC_ATLAS_sigmadn[n], y_sigma2up[n], y_sigma2dn[n];

  for(int i = 0; i<n; i++)
  {
     SA_sigmadn[i]  = SA_exp[i] - SA_expdn[i];
     SA_sigmaup[i]  = SA_expup[i] - SA_exp[i];
     CC_CMS_sigmadn[i]  = CC_CMS_exp[i] - CC_CMS_expdn[i];
     CC_CMS_sigmaup[i]  = CC_CMS_expup[i] - CC_CMS_exp[i];
     CC_ATLAS_sigmadn[i]  = CC_ATLAS_exp[i] - CC_ATLAS_expdn[i];
     CC_ATLAS_sigmaup[i]  = CC_ATLAS_expup[i] - CC_ATLAS_exp[i];

     y_sigma2dn[i] = SA_exp[i] - y_exp2dn[i];
     y_sigma2up[i] = y_exp2up[i] - SA_exp[i];
  }


  TCanvas *c1 = new TCanvas();
  if(displayExclXsec) { c1->SetLogy(1);}
  c1->cd();

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetName("");

  TGraphErrors* gr_obs = new TGraphErrors(n, x, y_obs, 0, y_obs_err);
  gr_obs->SetLineWidth(2);
  gr_obs->SetMarkerStyle(20);

  TGraphErrors* gr_SA_exp = new TGraphErrors(n, x, SA_exp, 0, 0);
  gr_SA_exp->SetLineWidth(2);
  gr_SA_exp->SetLineStyle(2);
  gr_SA_exp->SetMarkerStyle(22);

  TGraph* gr_CC_CMS_exp = new TGraph(n, x, CC_CMS_exp);
  gr_CC_CMS_exp->SetLineColor(kBlack);
  gr_CC_CMS_exp->SetLineWidth(2);
  gr_CC_CMS_exp->SetMarkerStyle(20);

  TGraph* gr_CC_ATLAS_exp = new TGraph(n, x, CC_ATLAS_exp);
  gr_CC_ATLAS_exp->SetLineColor(kBlack);
  gr_CC_ATLAS_exp->SetLineWidth(2);
  gr_CC_ATLAS_exp->SetMarkerStyle(21);

  TGraphAsymmErrors* sigma2 = new TGraphAsymmErrors(n, x, SA_exp, 0, 0, y_sigma2dn, y_sigma2up);
  sigma2->SetFillColor(kYellow-4);
  sigma2->SetLineColor(kYellow-4);

  TGraphAsymmErrors* sigma1_SA = new TGraphAsymmErrors(n, x, SA_exp, 0, 0, SA_sigmadn, SA_sigmaup);
  sigma1_SA->SetFillColor(kGreen-3);
  sigma1_SA->SetLineColor(kGreen-3);

  TGraphAsymmErrors* sigma1_CC_CMS = new TGraphAsymmErrors(n, x, CC_CMS_exp, 0, 0, CC_CMS_sigmadn, CC_CMS_sigmaup);
  sigma1_CC_CMS->SetFillColor(kYellow-4);
  sigma1_CC_CMS->SetLineColor(kYellow-4);

  TGraphAsymmErrors* sigma1_CC_ATLAS = new TGraphAsymmErrors(n, x, CC_ATLAS_exp, 0, 0, CC_ATLAS_sigmadn, CC_ATLAS_sigmaup);
  sigma1_CC_ATLAS->SetFillColor(kOrange);
  sigma1_CC_ATLAS->SetLineColor(kOrange);

  TGraphErrors* gr_SA_expup = new TGraphErrors(n, x, SA_expup, 0, 0);
  gr_SA_expup->SetLineColor(kRed);
  gr_SA_expup->SetLineWidth(2);
  gr_SA_expup->SetLineStyle(2);

  TGraphErrors* gr_SA_expdn = new TGraphErrors(n, x, SA_expdn, 0, 0);
  gr_SA_expdn->SetLineColor(kRed);
  gr_SA_expdn->SetLineWidth(2);
  gr_SA_expdn->SetLineStyle(2);

  TGraphErrors* gr_exp2up = new TGraphErrors(n, x, y_exp2up, 0, 0);
  gr_exp2up->SetLineColor(kBlue);
  gr_exp2up->SetLineWidth(2);
  gr_exp2up->SetLineStyle(4);

  TGraphErrors* gr_exp2dn = new TGraphErrors(n, x, y_exp2dn, 0, 0);
  gr_exp2dn->SetLineColor(kBlue);
  gr_exp2dn->SetLineWidth(2);
  gr_exp2dn->SetLineStyle(4);

  //mg->Add(sigma2);
  mg->Add(sigma1_SA);
  mg->Add(sigma1_CC_CMS);
  mg->Add(sigma1_CC_ATLAS);
  mg->Draw("3AC");
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetTitleOffset(.8);
  mg->GetYaxis()->SetTitleOffset(.7);
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.055);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetXaxis()->SetTitle("m_{Inv.} [GeV]");
  mg->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  if(displayExclXsec)
  {
      mg->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
      mg->GetYaxis()->SetTitleSize(0.04);
      mg->GetYaxis()->SetTitleOffset(1.2);
      //mg->SetMinimum(0.01);
      mg->SetMinimum(SA_expdn[n-1]);
      mg->SetMaximum(x_sec[0]);
  }
  else
  {
      mg->GetYaxis()->SetTitle("signal strength");
      mg->GetYaxis()->SetRangeUser(0.,y_exp2up[n-1]);
  }

  gr_SA_exp  ->Draw("samelp");
  gr_CC_CMS_exp  ->Draw("samelp");
  gr_CC_ATLAS_exp  ->Draw("samelp");
  //gr_obs  ->Draw("samelp");

  TLine* line = new TLine(x[0], 1., x[n-1], 1.);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  if(!displayExclXsec) line->Draw("same");
  else                 {xsec_th_0p1->Draw("same");xsec_th_0p05->Draw("same");}

  TLegend *leg;
  if(!displayExclXsec) leg = new TLegend(0.3,0.65,0.55,0.85);
  else                 leg = new TLegend(0.55,0.60,0.85,0.85);
  leg->AddEntry(gr_SA_exp    ,"Shape Ana. CMS"    ,"lp");
  leg->AddEntry(sigma1_SA  ," Expected #pm1#sigma (Shape Ana.)","f");
  //leg->AddEntry(sigma2  ," Expected #pm2#sigma","f");
  leg->AddEntry(gr_CC_CMS_exp ,"Cut&Count CMS Sel." ,"lp");
  leg->AddEntry(sigma1_CC_CMS  ," Expected #pm1#sigma (Cut&Count CMS)","f");
  leg->AddEntry(gr_CC_ATLAS_exp ,"Cut&Count ATLAS Sel." ,"lp");
  leg->AddEntry(sigma1_CC_ATLAS  ," Expected #pm1#sigma (Cut&Count ATLAS)","f");
  if(displayExclXsec) {leg->AddEntry(xsec_th_0p1  ," Theory (LO), a_{non-res} = 0.1","l"); leg->AddEntry(xsec_th_0p05  ," Theory (LO), a_{non-res} = 0.05","l");}
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
  else                 text_2 = new TPaveText(0.6, 0.48, 0.85, 0.58, "NDC");
  text_2->AddText("Non-resonant model");
  text_2->AddText("Expected 95\% CL limits");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(18);
  text_2->Draw();

  TString             Y_axis = "SignalStrength";
  if(displayExclXsec) Y_axis = "ExcludedXsection";

  TString externOfSyst;
  if(isExtern) externOfSyst = "_Extern";
  else         externOfSyst = "_noExtern";

  c1->SaveAs("limitPlots/finalLimits_FCNC"+externOfSyst+"_methodComparison_"+name+".png");
  c1->SaveAs("limitPlots/finalLimits_FCNC"+externOfSyst+"_methodComparison_"+name+".pdf");
  c1->SaveAs("limitPlots/finalLimits_FCNC"+externOfSyst+"_methodComparison_"+name+".eps");

}

void CompareLimitMethod_FCNC()
{
    CompareLimitMethod();
}
