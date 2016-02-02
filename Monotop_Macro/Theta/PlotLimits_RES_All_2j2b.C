#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TString.h>
#include <sstream>
#include <TLatex.h>


#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TMultiGraph.h"

  const unsigned int nRes = 9;   // number of resonant mass points
  const unsigned int nInv = 5;   // number of invisible mass points

void getExclusionMass(Double_t * x,  Double_t* x_sec, Double_t* y_obs, Double_t* limits)
{

  TGraphErrors* gr_obs  = new TGraphErrors(nRes, x, y_obs, 0, 0);
  TGraphErrors* gr_pred = new TGraphErrors(nRes, x, x_sec, 0, 0);

  gr_obs->Draw();
  gr_pred->Draw("same");


  int ncross_obs = 0;

  for(unsigned int i=1; i<(nRes-1); i++)
  {
    if( x_sec[i] > y_obs[i] && x_sec[i+1] < y_obs[i+1] )
    {
      ncross_obs = i;
      break;
    }
  }


  TF1 * line_obs = new TF1("line_obs", "pol1", x[ncross_obs], x[ncross_obs+1]);
  gr_obs->Fit(line_obs, "R");


  TF1 * line_sec = new TF1("line_sec", "pol1", x[ncross_obs], x[ncross_obs+1]);
  gr_pred->Fit(line_sec, "R");

  float a1 = line_obs->GetParameter(1);
  float b1 = line_obs->GetParameter(0);

  float a2 = line_sec->GetParameter(1);
  float b2 = line_sec->GetParameter(0);

  float x_val = (b2-b1)/(a1-a2);
  float y_val = a1*x_val+b1;

  limits[0] = x_val;
  limits[1] = y_val;

  delete line_obs, line_sec;

}


void PlotLimits_RES_All_2j2b()
{

  Double_t x[nRes]     = {500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100};

  //---------------------------------------------------------------------
  //for invariant mass of 10 GeV
  Double_t x_sec_minv10[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv10[nRes]    = {0.0147,  0.0207,  0.0369,  0.0815,  0.2333, 0.561,  1.422,  3.854,  10.78 };
  Double_t y_exp2up_minv10[nRes] = {0.0293,  0.0496,  0.0789,  0.1707,  0.4831, 1.159,  2.986,  7.912,  22.47 };
  Double_t y_expup_minv10[nRes]  = {0.0181,  0.0294,  0.0498,  0.1078,  0.3048, 0.721,  1.853,  4.944,  13.60 };
  Double_t y_exp_minv10[nRes]    = {0.0114,  0.0182,  0.0320,  0.0685,  0.1945, 0.470,  1.178,  3.184,  8.724 };
  Double_t y_expdn_minv10[nRes]  = {0.0073,  0.0118,  0.0211,  0.0466,  0.1314, 0.318,  0.794,  2.122,  5.863 };
  Double_t y_exp2dn_minv10[nRes] = {0.0049,  0.0082,  0.0152,  0.0347,  0.0973, 0.237,  0.592,  1.586,  4.332 };

  Double_t limits_obs_minv10[2]     = {0,0};
  Double_t limits_exp2up_minv10[2]  = {0,0};
  Double_t limits_expup_minv10[2]   = {0,0};
  Double_t limits_exp_minv10[2]     = {0,0};
  Double_t limits_expdn_minv10[2]   = {0,0};
  Double_t limits_exp2dn_minv10[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 50 GeV
  Double_t x_sec_minv50[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv50[nRes]    = {0.0148,  0.0219,  0.0385,  0.0827,  0.2398, 0.560,  1.404,  4.548,  9.77  };
  Double_t y_exp2up_minv50[nRes] = {0.0303,  0.0531,  0.0854,  0.1794,  0.5193, 1.300,  2.750,  8.848, 21.95  };
  Double_t y_expup_minv50[nRes]  = {0.0182,  0.0306,  0.0511,  0.1128,  0.3171, 0.786,  1.758,  5.424, 13.25  };
  Double_t y_exp_minv50[nRes]    = {0.0114,  0.0189,  0.0324,  0.0703,  0.2015, 0.491,  1.114,  3.463, 8.464  };
  Double_t y_expdn_minv50[nRes]  = {0.0073,  0.0126,  0.0217,  0.0470,  0.1353, 0.330,  0.754,  2.312, 5.717  };
  Double_t y_exp2dn_minv50[nRes] = {0.0049,  0.0088,  0.0156,  0.0348,  0.1015, 0.244,  0.577,  1.728, 4.225  };

  Double_t limits_obs_minv50[2]     = {0,0};
  Double_t limits_exp2up_minv50[2]  = {0,0};
  Double_t limits_expup_minv50[2]   = {0,0};
  Double_t limits_exp_minv50[2]     = {0,0};
  Double_t limits_expdn_minv50[2]   = {0,0};
  Double_t limits_exp2dn_minv50[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 100 GeV
  Double_t x_sec_minv100[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv100[nRes]    = {0.0169,  0.0221,  0.0382,  0.0855,  0.3648, 0.560,  1.510,  4.110,  11.09 };
  Double_t y_exp2up_minv100[nRes] = {0.0323,  0.0506,  0.0809,  0.1841,  0.7841, 1.154,  3.214,  8.852,  23.76 };
  Double_t y_expup_minv100[nRes]  = {0.0198,  0.0303,  0.0493,  0.1139,  0.5028, 0.706,  1.978,  5.364,  14.02 };
  Double_t y_exp_minv100[nRes]    = {0.0125,  0.0192,  0.0315,  0.0722,  0.3169, 0.462,  1.258,  3.352,  8.878 };
  Double_t y_expdn_minv100[nRes]  = {0.0082,  0.0127,  0.0211,  0.0487,  0.2084, 0.318,  0.841,  2.275,  5.963 };
  Double_t y_exp2dn_minv100[nRes] = {0.0056,  0.0089,  0.0152,  0.0358,  0.1491, 0.239,  0.638,  1.701,  4.452 };

  Double_t limits_obs_minv100[2]    = {0,0};
  Double_t limits_exp2up_minv100[2] = {0,0};
  Double_t limits_expup_minv100[2]  = {0,0};
  Double_t limits_exp_minv100[2]    = {0,0};
  Double_t limits_expdn_minv100[2]  = {0,0};
  Double_t limits_exp2dn_minv100[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 150 GeV
  Double_t x_sec_minv150[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv150[nRes]    = {0.0211,  0.0245,  0.0409,  0.0847,  0.2323, 0.542,  1.495,  4.195,  10.92 };
  Double_t y_exp2up_minv150[nRes] = {0.0415,  0.0542,  0.0914,  0.1908,  0.4591, 1.138,  2.994,  9.037,  21.55 };
  Double_t y_expup_minv150[nRes]  = {0.0253,  0.0331,  0.0569,  0.1158,  0.2988, 0.719,  1.874,  5.604,  13.74 };
  Double_t y_exp_minv150[nRes]    = {0.0158,  0.0209,  0.0356,  0.0728,  0.1924, 0.458,  1.215,  3.514,  8.844 };
  Double_t y_expdn_minv150[nRes]  = {0.0102,  0.0135,  0.0232,  0.0482,  0.1294, 0.310,  0.825,  2.342,  6.003 };
  Double_t y_exp2dn_minv150[nRes] = {0.0071,  0.0097,  0.0167,  0.0357,  0.0952, 0.232,  0.619,  1.758,  4.491 };

  Double_t limits_obs_minv150[2]    = {0,0};
  Double_t limits_exp2up_minv150[2] = {0,0};
  Double_t limits_expup_minv150[2]  = {0,0};
  Double_t limits_exp_minv150[2]    = {0,0};
  Double_t limits_expdn_minv150[2]  = {0,0};
  Double_t limits_exp2dn_minv150[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 200 GeV
  Double_t x_sec_minv200[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv200[nRes]    = {0.0272,  0.0288,  0.0413,  0.0911,  0.2273, 0.551,  1.430,  4.098,  11.90 };
  Double_t y_exp2up_minv200[nRes] = {0.0557,  0.0671,  0.0925,  0.2013,  0.4768, 1.146,  2.917,  8.812,  22.94 };
  Double_t y_expup_minv200[nRes]  = {0.0346,  0.0398,  0.0553,  0.1218,  0.3008, 0.717,  1.833,  5.435,  14.54 };
  Double_t y_exp_minv200[nRes]    = {0.0214,  0.0249,  0.0349,  0.0773,  0.1902, 0.457,  1.190,  3.464,  9.244 };
  Double_t y_expdn_minv200[nRes]  = {0.0135,  0.0161,  0.0234,  0.0516,  0.1280, 0.314,  0.814,  2.342,  6.286 };
  Double_t y_exp2dn_minv200[nRes] = {0.0093,  0.0114,  0.0172,  0.0377,  0.0962, 0.235,  0.610,  1.716,  4.728 };

  Double_t limits_obs_minv200[2]    = {0,0};
  Double_t limits_exp2up_minv200[2] = {0,0};
  Double_t limits_expup_minv200[2]  = {0,0};
  Double_t limits_exp_minv200[2]    = {0,0};
  Double_t limits_expdn_minv200[2]  = {0,0};
  Double_t limits_exp2dn_minv200[2] = {0,0};


  for(unsigned short int ind = 0; ind < nRes; ind++)
  {
    y_obs_minv10[ind]   *= x_sec_minv10[ind];
    y_exp2up_minv10[ind]*= x_sec_minv10[ind];
    y_expup_minv10[ind] *= x_sec_minv10[ind];
    y_exp_minv10[ind]   *= x_sec_minv10[ind];
    y_expdn_minv10[ind] *= x_sec_minv10[ind];
    y_exp2dn_minv10[ind]*= x_sec_minv10[ind];

    y_obs_minv50[ind]   *= x_sec_minv50[ind];
    y_exp2up_minv50[ind]*= x_sec_minv50[ind];
    y_expup_minv50[ind] *= x_sec_minv50[ind];
    y_exp_minv50[ind]   *= x_sec_minv50[ind];
    y_expdn_minv50[ind] *= x_sec_minv50[ind];
    y_exp2dn_minv50[ind]*= x_sec_minv50[ind];

    y_obs_minv100[ind]   *= x_sec_minv100[ind];
    y_exp2up_minv100[ind]*= x_sec_minv100[ind];
    y_expup_minv100[ind] *= x_sec_minv100[ind];
    y_exp_minv100[ind]   *= x_sec_minv100[ind];
    y_expdn_minv100[ind] *= x_sec_minv100[ind];
    y_exp2dn_minv100[ind]*= x_sec_minv100[ind];

    y_obs_minv150[ind]   *= x_sec_minv150[ind];
    y_exp2up_minv150[ind]*= x_sec_minv150[ind];
    y_expup_minv150[ind] *= x_sec_minv150[ind];
    y_exp_minv150[ind]   *= x_sec_minv150[ind];
    y_expdn_minv150[ind] *= x_sec_minv150[ind];
    y_exp2dn_minv150[ind]*= x_sec_minv150[ind];

    y_obs_minv200[ind]   *= x_sec_minv200[ind];
    y_exp2up_minv200[ind]*= x_sec_minv200[ind];
    y_expup_minv200[ind] *= x_sec_minv200[ind];
    y_exp_minv200[ind]   *= x_sec_minv200[ind];
    y_expdn_minv200[ind] *= x_sec_minv200[ind];
    y_exp2dn_minv200[ind]*= x_sec_minv200[ind];


  }

  getExclusionMass(x, x_sec_minv10,  y_obs_minv10,  limits_obs_minv10);
  getExclusionMass(x, x_sec_minv50,  y_obs_minv50,  limits_obs_minv50);
  getExclusionMass(x, x_sec_minv100, y_obs_minv100, limits_obs_minv100);
  getExclusionMass(x, x_sec_minv150, y_obs_minv150, limits_obs_minv150);
  getExclusionMass(x, x_sec_minv200, y_obs_minv200, limits_obs_minv200);


  getExclusionMass(x, x_sec_minv10,  y_exp2up_minv10,  limits_exp2up_minv10);
  getExclusionMass(x, x_sec_minv50,  y_exp2up_minv50,  limits_exp2up_minv50);
  getExclusionMass(x, x_sec_minv100, y_exp2up_minv100, limits_exp2up_minv100);
  getExclusionMass(x, x_sec_minv150, y_exp2up_minv150, limits_exp2up_minv150);
  getExclusionMass(x, x_sec_minv200, y_exp2up_minv200, limits_exp2up_minv200);


  getExclusionMass(x, x_sec_minv10,  y_expup_minv10,  limits_expup_minv10);
  getExclusionMass(x, x_sec_minv50,  y_expup_minv50,  limits_expup_minv50);
  getExclusionMass(x, x_sec_minv100, y_expup_minv100, limits_expup_minv100);
  getExclusionMass(x, x_sec_minv150, y_expup_minv150, limits_expup_minv150);
  getExclusionMass(x, x_sec_minv200, y_expup_minv200, limits_expup_minv200);


  getExclusionMass(x, x_sec_minv10,  y_exp_minv10,  limits_exp_minv10);
  getExclusionMass(x, x_sec_minv50,  y_exp_minv50,  limits_exp_minv50);
  getExclusionMass(x, x_sec_minv100, y_exp_minv100, limits_exp_minv100);
  getExclusionMass(x, x_sec_minv150, y_exp_minv150, limits_exp_minv150);
  getExclusionMass(x, x_sec_minv200, y_exp_minv200, limits_exp_minv200);


  getExclusionMass(x, x_sec_minv10,  y_expdn_minv10,  limits_expdn_minv10);
  getExclusionMass(x, x_sec_minv50,  y_expdn_minv50,  limits_expdn_minv50);
  getExclusionMass(x, x_sec_minv100, y_expdn_minv100, limits_expdn_minv100);
  getExclusionMass(x, x_sec_minv150, y_expdn_minv150, limits_expdn_minv150);
  getExclusionMass(x, x_sec_minv200, y_expdn_minv200, limits_expdn_minv200);


  getExclusionMass(x, x_sec_minv10,  y_exp2dn_minv10,  limits_exp2dn_minv10);
  getExclusionMass(x, x_sec_minv50,  y_exp2dn_minv50,  limits_exp2dn_minv50);
  getExclusionMass(x, x_sec_minv100, y_exp2dn_minv100, limits_exp2dn_minv100);
  getExclusionMass(x, x_sec_minv150, y_exp2dn_minv150, limits_exp2dn_minv150);
  getExclusionMass(x, x_sec_minv200, y_exp2dn_minv200, limits_exp2dn_minv200);



  double excluded_obs[nInv]     = { 0, 0, 0, 0, 0};
  double excluded_exp2up[nInv]  = { 0, 0, 0, 0, 0};
  double excluded_expup[nInv]   = { 0, 0, 0, 0, 0};
  double excluded_exp[nInv]     = { 0, 0, 0, 0, 0};
  double excluded_expdn[nInv]   = { 0, 0, 0, 0, 0};
  double excluded_exp2dn[nInv]  = { 0, 0, 0, 0, 0};

  double inv_mass[nInv] = {10, 50, 100, 150, 200};


  for(unsigned int i=0; i<nInv; i++)
  {
    if(i==0)
    {
      excluded_obs[0]    = limits_obs_minv10[0];
      excluded_exp2up[0] = limits_exp2up_minv10[0];
      excluded_expup[0]  = limits_expup_minv10[0];
      excluded_exp[0]    = limits_exp_minv10[0];
      excluded_expdn[0]  = limits_expdn_minv10[0];
      excluded_exp2dn[0] = limits_exp2dn_minv10[0];
    }
    if(i==1)
    {
      excluded_obs[1]    = limits_obs_minv50[0];
      excluded_exp2up[1] = limits_exp2up_minv50[0];
      excluded_expup[1]  = limits_expup_minv50[0];
      excluded_exp[1]    = limits_exp_minv50[0];
      excluded_expdn[1]  = limits_expdn_minv50[0];
      excluded_exp2dn[1] = limits_exp2dn_minv50[0];
    }
    if(i==2)
    {
      excluded_obs[2]    = limits_obs_minv100[0];
      excluded_exp2up[2] = limits_exp2up_minv100[0];
      excluded_expup[2]  = limits_expup_minv100[0];
      excluded_exp[2]    = limits_exp_minv100[0];
      excluded_expdn[2]  = limits_expdn_minv100[0];
      excluded_exp2dn[2] = limits_exp2dn_minv100[0];
    }
    if(i==3)
    {
      excluded_obs[3]    = limits_obs_minv150[0];
      excluded_exp2up[3] = limits_exp2up_minv150[0];
      excluded_expup[3]  = limits_expup_minv150[0];
      excluded_exp[3]    = limits_exp_minv150[0];
      excluded_expdn[3]  = limits_expdn_minv150[0];
      excluded_exp2dn[3] = limits_exp2dn_minv150[0];
    }
    if(i==4)
    {
      excluded_obs[4]    = limits_obs_minv200[0];
      excluded_exp2up[4] = limits_exp2up_minv200[0];
      excluded_expup[4]  = limits_expup_minv200[0];
      excluded_exp[4]    = limits_exp_minv200[0];
      excluded_expdn[4]  = limits_expdn_minv200[0];
      excluded_exp2dn[4] = limits_exp2dn_minv200[0];
    }
  }

  //for(unsigned int i=0; i<nInv; i++) cout << "mass limit obs "<< i << "  "  << excluded_obs[i] << endl;
  ////for(unsigned int i=0; i<nInv; i++) cout << "mass limit exp "<< i << "  "  << excluded_exp[i] << endl;
  ////for(unsigned int i=0; i<nInv; i++) cout << "mass limit exp2up "<< i << "  "  << excluded_exp2up[i] << endl;
  //for(unsigned int i=0; i<nInv; i++) cout << "mass limit expup "<< i << "  "  << excluded_expup[i] << endl;
  //for(unsigned int i=0; i<nInv; i++) cout << "mass limit expdn "<< i << "  "  << excluded_expdn[i] << endl;
  ////for(unsigned int i=0; i<nInv; i++) cout << "mass limit exp2dn "<< i << "  "  << excluded_exp2dn[i] << endl;

  Double_t y_sigmaup[nInv], y_sigmadn[nInv], y_sigma2up[nInv], y_sigma2dn[nInv];

  for(unsigned int i = 0; i<nInv; i++)
  {
     y_sigmadn[i]  = excluded_exp[i] - excluded_expdn[i];
     y_sigmaup[i]  = excluded_expup[i] - excluded_exp[i];
     y_sigma2dn[i] = excluded_exp[i] - excluded_exp2dn[i];
     y_sigma2up[i] = excluded_exp2up[i] - excluded_exp[i];
  }



  TGraphErrors* limits_obs     = new TGraphErrors(nInv, inv_mass, excluded_obs, 0, 0);
  limits_obs->SetDrawOption(0);
  limits_obs->SetLineWidth(2);
  limits_obs->SetMarkerStyle(20);

  TGraphErrors* limits_exp     = new TGraphErrors(nInv, inv_mass, excluded_exp, 0, 0);
  limits_exp->SetLineWidth(2);
  limits_exp->SetLineStyle(2);
  limits_exp->SetMarkerStyle(22);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  TGraphErrors* limits_exp2up  = new TGraphErrors(nInv, inv_mass, excluded_exp2up, 0, 0);
  limits_exp2up->SetLineColor(kBlue);
  limits_exp2up->SetLineWidth(2);
  limits_exp2up->SetLineStyle(4);


  TGraphErrors* limits_expup   = new TGraphErrors(nInv, inv_mass, excluded_expup, 0, 0);
  limits_expup->SetLineColor(kRed);
  limits_expup->SetLineWidth(2);
  limits_expup->SetLineStyle(2);

  TGraphErrors* limits_expdn   = new TGraphErrors(nInv, inv_mass, excluded_expdn, 0, 0);
  limits_expdn->SetLineColor(kRed);
  limits_expdn->SetLineWidth(2);
  limits_expdn->SetLineStyle(2);


  TGraphErrors* limits_exp2dn  = new TGraphErrors(nInv, inv_mass, excluded_exp2dn, 0, 0);
  limits_exp2dn->SetLineColor(kBlue);
  limits_exp2dn->SetLineWidth(2);
  limits_exp2dn->SetLineStyle(4);


  double zeros[nInv] = {0, 0, 0, 0, 0};

  TGraphAsymmErrors* sigma2 = new TGraphAsymmErrors(nInv, inv_mass, excluded_exp, 0, 0 , y_sigma2dn , y_sigma2up);
  sigma2->SetFillColor(kYellow);
  sigma2->SetLineColor(kYellow);

  //TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(nInv, inv_mass, excluded_exp, zeros, zeros, excluded_expup, excluded_expdn);
  TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(nInv, inv_mass, excluded_exp, 0, 0 , y_sigmadn, y_sigmaup);
  sigma1->SetFillColor(kGreen);
  sigma1->SetLineColor(kGreen);

  for(unsigned int i=0; i<nInv; i++)  cout << "mass limit "<< inv_mass[i] << "  "  << excluded_obs[i] << "  exp " << excluded_exp[i] << " ["<< excluded_expup[i]<< "," << excluded_expdn[i] << "]"  << " ["<< excluded_exp2up[i]<< "," << excluded_exp2dn[i] << "]"  << endl;

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(sigma2);
  mg->Add(sigma1);
  mg->Draw("3AC");
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetTitleOffset(.8);
  mg->GetYaxis()->SetTitleOffset(1.0);
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.050);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetXaxis()->SetTitle("m_{Inv.} [GeV]");
  mg->GetYaxis()->SetTitle("exclusion limit on m_{Res.} [GeV]");
  mg->GetXaxis()->SetRangeUser(x[0],x[nRes-1]);
  mg->SetMaximum(2100);
  mg->SetMinimum(1350);


  limits_obs->Draw("samelp");
  limits_exp ->Draw("samelp");


  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry( limits_obs ," Observed 95\% CL limit","lpe");
  leg->AddEntry( limits_exp ," Expected 95\% CL limit","lp");
  leg->AddEntry( sigma1 ," Expected #pm 1 #sigma ","f");
  leg->AddEntry( sigma2 ," Expected #pm 2 #sigma ","f");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();


  TLatex* text = new TLatex(0.17, 0.92, "CMS Preliminary, L = 19.7 fb^{-1}, #sqrt{s} = 8 TeV, a=0.1");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw();

  c1->SaveAs("limitPlots/finalLimits_RES_All_ExcludedXsection_2j2b.png");
  c1->SaveAs("limitPlots/finalLimits_RES_All_ExcludedXsection_2j2b.pdf");
  c1->SaveAs("limitPlots/finalLimits_RES_All_ExcludedXsection_2j2b.eps");

}
