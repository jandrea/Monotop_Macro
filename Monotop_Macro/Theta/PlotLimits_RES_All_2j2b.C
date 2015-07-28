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
  Double_t y_obs_minv10[nRes]    = {0.0124, 0.0196, 0.0365, 0.0779, 0.2333, 0.557, 1.465, 3.824 , 11.005 };
  Double_t y_exp2up_minv10[nRes] = {0.0307, 0.0451, 0.0842, 0.1807, 0.5301, 1.239, 3.186, 8.432 , 24.07  };
  Double_t y_expup_minv10[nRes]  = {0.0195, 0.0283, 0.0512, 0.1118, 0.3208, 0.752, 1.923, 5.034 , 14.52  };
  Double_t y_exp_minv10[nRes]    = {0.0121, 0.0179, 0.0327, 0.0707, 0.2035, 0.469, 1.221, 3.144 , 8.984  };
  Double_t y_expdn_minv10[nRes]  = {0.0078, 0.0119, 0.0217, 0.0476, 0.1374, 0.317, 0.826, 2.122 , 6.043  };
  Double_t y_exp2dn_minv10[nRes] = {0.0052, 0.0086, 0.0156, 0.0354, 0.1032, 0.239, 0.619, 1.576 , 4.512  };

  Double_t limits_obs_minv10[2]     = {0,0};
  Double_t limits_exp2up_minv10[2]  = {0,0};
  Double_t limits_expup_minv10[2]   = {0,0};
  Double_t limits_exp_minv10[2]     = {0,0};
  Double_t limits_expdn_minv10[2]   = {0,0};
  Double_t limits_exp2dn_minv10[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 50 GeV
  Double_t x_sec_minv50[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv50[nRes]    = {0.0126, 0.0206, 0.0362, 0.0772, 0.2378, 0.562, 1.464, 4.418 , 10.20  };
  Double_t y_exp2up_minv50[nRes] = {0.0293, 0.0473, 0.0784, 0.1729, 0.5203, 1.270, 3.300, 8.848 , 22.75  };
  Double_t y_expup_minv50[nRes]  = {0.0185, 0.0296, 0.0490, 0.1078, 0.3128, 0.748, 1.918, 5.494 , 13.81  };
  Double_t y_exp_minv50[nRes]    = {0.0117, 0.0188, 0.0319, 0.0690, 0.1955, 0.469, 1.194, 3.463 , 8.444  };
  Double_t y_expdn_minv50[nRes]  = {0.0076, 0.0125, 0.0213, 0.0460, 0.1313, 0.314, 0.803, 2.302 , 5.677  };
  Double_t y_exp2dn_minv50[nRes] = {0.0051, 0.0088, 0.0157, 0.0342, 0.0985, 0.238, 0.599, 1.728 , 4.275  };

  Double_t limits_obs_minv50[2]     = {0,0};
  Double_t limits_exp2up_minv50[2]  = {0,0};
  Double_t limits_expup_minv50[2]   = {0,0};
  Double_t limits_exp_minv50[2]     = {0,0};
  Double_t limits_expdn_minv50[2]   = {0,0};
  Double_t limits_exp2dn_minv50[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 100 GeV
  Double_t x_sec_minv100[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv100[nRes]    = {0.0138, 0.0203, 0.0362, 0.0820, 0.3518, 0.562, 1.537, 4.080 , 11.095 };
  Double_t y_exp2up_minv100[nRes] = {0.0339, 0.0481, 0.0812, 0.1841, 0.7671, 1.244, 3.424, 8.832 , 23.76  };
  Double_t y_expup_minv100[nRes]  = {0.0216, 0.0303, 0.0510, 0.1129, 0.4738, 0.757, 2.038, 5.254 , 14.02  };
  Double_t y_exp_minv100[nRes]    = {0.0136, 0.0193, 0.0327, 0.0724, 0.2969, 0.478, 1.288, 3.282 , 8.878  };
  Double_t y_expdn_minv100[nRes]  = {0.0087, 0.0128, 0.0218, 0.0489, 0.2004, 0.323, 0.859, 2.215 , 5.963  };
  Double_t y_exp2dn_minv100[nRes] = {0.0059, 0.0091, 0.0159, 0.0361, 0.1491, 0.242, 0.649, 1.661 , 4.452  };

  Double_t limits_obs_minv100[2]    = {0,0};
  Double_t limits_exp2up_minv100[2] = {0,0};
  Double_t limits_expup_minv100[2]  = {0,0};
  Double_t limits_exp_minv100[2]    = {0,0};
  Double_t limits_expdn_minv100[2]  = {0,0};
  Double_t limits_exp2dn_minv100[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 150 GeV
  Double_t x_sec_minv150[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv150[nRes]    = {0.0171, 0.0238, 0.0367, 0.0789, 0.2293, 0.536, 1.495, 4.055 , 11.419 };
  Double_t y_exp2up_minv150[nRes] = {0.0425, 0.0531, 0.0862, 0.1818, 0.5171, 1.138, 3.244, 9.037 , 26.55  };
  Double_t y_expup_minv150[nRes]  = {0.0263, 0.0344, 0.0529, 0.1088, 0.3178, 0.703, 1.954, 5.434 , 15.61  };
  Double_t y_exp_minv150[nRes]    = {0.0165, 0.0219, 0.0340, 0.0699, 0.2014, 0.445, 1.227, 3.384 , 9.724  };
  Double_t y_expdn_minv150[nRes]  = {0.0106, 0.0145, 0.0228, 0.0473, 0.1364, 0.303, 0.813, 2.272 , 6.513  };
  Double_t y_exp2dn_minv150[nRes] = {0.0071, 0.0100, 0.0164, 0.0350, 0.1016, 0.228, 0.619, 1.708 , 4.811  };

  Double_t limits_obs_minv150[2]    = {0,0};
  Double_t limits_exp2up_minv150[2] = {0,0};
  Double_t limits_expup_minv150[2]  = {0,0};
  Double_t limits_exp_minv150[2]    = {0,0};
  Double_t limits_expdn_minv150[2]  = {0,0};
  Double_t limits_exp2dn_minv150[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 200 GeV
  Double_t x_sec_minv200[nRes]    = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t y_obs_minv200[nRes]    = {0.0203, 0.0266, 0.0380, 0.0903, 0.2263, 0.541, 1.440, 3.838 , 11.504 };
  Double_t y_exp2up_minv200[nRes] = {0.0584, 0.0631, 0.0844, 0.2013, 0.5138, 1.200, 3.006, 9.022 , 25.14  };
  Double_t y_expup_minv200[nRes]  = {0.0365, 0.0396, 0.0530, 0.1258, 0.3158, 0.727, 1.893, 5.265 , 15.19  };
  Double_t y_exp_minv200[nRes]    = {0.0224, 0.0249, 0.0340, 0.0779, 0.1985, 0.462, 1.198, 3.264 , 9.544  };
  Double_t y_expdn_minv200[nRes]  = {0.0142, 0.0163, 0.0229, 0.0527, 0.1330, 0.310, 0.809, 2.202 , 6.366  };
  Double_t y_exp2dn_minv200[nRes] = {0.0095, 0.0113, 0.0167, 0.0384, 0.0981, 0.233, 0.614, 1.656 , 4.718  };

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
