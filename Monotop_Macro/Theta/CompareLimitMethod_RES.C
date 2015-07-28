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

void getExclusionMass(Double_t * x,  Double_t* x_sec, Double_t* y_exp, Double_t* limits)
{

  TGraphErrors* gr_exp  = new TGraphErrors(nRes, x, y_exp, 0, 0);
  TGraphErrors* gr_pred = new TGraphErrors(nRes, x, x_sec, 0, 0);

  gr_exp->Draw();
  gr_pred->Draw("same");


  int ncross_exp = 0;

  for(unsigned int i=1; i<(nRes-1); i++)
  {
    if( x_sec[i] > y_exp[i] && x_sec[i+1] < y_exp[i+1] )
    {
      ncross_exp = i;
      break;
    }
  }


  TF1 * line_exp = new TF1("line_exp", "pol1", x[ncross_exp], x[ncross_exp+1]);
  gr_exp->Fit(line_exp, "R");


  TF1 * line_sec = new TF1("line_sec", "pol1", x[ncross_exp], x[ncross_exp+1]);
  gr_pred->Fit(line_sec, "R");

  float a1 = line_exp->GetParameter(1);
  float b1 = line_exp->GetParameter(0);

  float a2 = line_sec->GetParameter(1);
  float b2 = line_sec->GetParameter(0);

  float x_val = (b2-b1)/(a1-a2);
  float y_val = a1*x_val+b1;

  cout << "x_val= " << x_val << " | y_val= " << y_val << endl;

  limits[0] = x_val;
  limits[1] = y_val;

  delete line_exp, line_sec;

}


void CompareLimitMethod_RES()
{

  Double_t x[nRes]     = {500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100};

  TString name = "CCvsSA";




//--------------------------------------------------------
//---------------- Shape Analysis Part  ------------------
//--------------------------------------------------------

  //---------------------------------------------------------------------
  //for invariant mass of 10 GeV
  Double_t x_sec_minv10[nRes]     = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t SA_expup_minv10[nRes]  = {0.0195, 0.0283, 0.0512, 0.1118, 0.3208, 0.752, 1.923, 5.034 , 14.52  };
  Double_t SA_exp_minv10[nRes]    = {0.0121, 0.0179, 0.0327, 0.0707, 0.2035, 0.469, 1.221, 3.144 , 8.984  };
  Double_t SA_expdn_minv10[nRes]  = {0.0078, 0.0119, 0.0217, 0.0476, 0.1374, 0.317, 0.826, 2.122 , 6.043  };

  Double_t limits_expup_SA_minv10[2]   = {0,0};
  Double_t limits_exp_SA_minv10[2]     = {0,0};
  Double_t limits_expdn_SA_minv10[2]   = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 50 GeV
  Double_t x_sec_minv50[nRes]     = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t SA_expup_minv50[nRes]  = {0.0185, 0.0296, 0.0490, 0.1078, 0.3128, 0.748, 1.918, 5.494 , 13.81  };
  Double_t SA_exp_minv50[nRes]    = {0.0117, 0.0188, 0.0319, 0.0690, 0.1955, 0.469, 1.194, 3.463 , 8.444  };
  Double_t SA_expdn_minv50[nRes]  = {0.0076, 0.0125, 0.0213, 0.0460, 0.1313, 0.314, 0.803, 2.302 , 5.677  };

  Double_t limits_expup_SA_minv50[2]   = {0,0};
  Double_t limits_exp_SA_minv50[2]     = {0,0};
  Double_t limits_expdn_SA_minv50[2]   = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 100 GeV
  Double_t x_sec_minv100[nRes]     = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t SA_expup_minv100[nRes]  = {0.0216, 0.0303, 0.0510, 0.1129, 0.4738, 0.757, 2.038, 5.254 , 14.02  };
  Double_t SA_exp_minv100[nRes]    = {0.0136, 0.0193, 0.0327, 0.0724, 0.2969, 0.478, 1.288, 3.282 , 8.878  };
  Double_t SA_expdn_minv100[nRes]  = {0.0087, 0.0128, 0.0218, 0.0489, 0.2004, 0.323, 0.859, 2.215 , 5.963  };

  Double_t limits_expup_SA_minv100[2]  = {0,0};
  Double_t limits_exp_SA_minv100[2]    = {0,0};
  Double_t limits_expdn_SA_minv100[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 150 GeV
  Double_t x_sec_minv150[nRes]     = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t SA_expup_minv150[nRes]  = {0.0263, 0.0344, 0.0529, 0.1088, 0.3178, 0.703, 1.954, 5.434 , 15.61  };
  Double_t SA_exp_minv150[nRes]    = {0.0165, 0.0219, 0.0340, 0.0699, 0.2014, 0.445, 1.227, 3.384 , 9.724  };
  Double_t SA_expdn_minv150[nRes]  = {0.0106, 0.0145, 0.0228, 0.0473, 0.1364, 0.303, 0.813, 2.272 , 6.513  };

  Double_t limits_expup_SA_minv150[2]  = {0,0};
  Double_t limits_exp_SA_minv150[2]    = {0,0};
  Double_t limits_expdn_SA_minv150[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 200 GeV
  Double_t x_sec_minv200[nRes]     = {11.8  , 2.77  , 0.87  , 0.31  , 0.11  , 0.05 , 0.022, 0.0097, 0.0045 };
  Double_t SA_expup_minv200[nRes]  = {0.0365, 0.0396, 0.0530, 0.1258, 0.3158, 0.727, 1.893, 5.265 , 15.19  };
  Double_t SA_exp_minv200[nRes]    = {0.0224, 0.0249, 0.0340, 0.0779, 0.1985, 0.462, 1.198, 3.264 , 9.544  };
  Double_t SA_expdn_minv200[nRes]  = {0.0142, 0.0163, 0.0229, 0.0527, 0.1330, 0.310, 0.809, 2.202 , 6.366  };

  Double_t limits_expup_SA_minv200[2]  = {0,0};
  Double_t limits_exp_SA_minv200[2]    = {0,0};
  Double_t limits_expdn_SA_minv200[2]  = {0,0};


  for(unsigned short int ind = 0; ind < nRes; ind++)
  {
    SA_expup_minv10[ind] *= x_sec_minv10[ind];
    SA_exp_minv10[ind]   *= x_sec_minv10[ind];
    SA_expdn_minv10[ind] *= x_sec_minv10[ind];

    SA_expup_minv50[ind] *= x_sec_minv50[ind];
    SA_exp_minv50[ind]   *= x_sec_minv50[ind];
    SA_expdn_minv50[ind] *= x_sec_minv50[ind];

    SA_expup_minv100[ind] *= x_sec_minv100[ind];
    SA_exp_minv100[ind]   *= x_sec_minv100[ind];
    SA_expdn_minv100[ind] *= x_sec_minv100[ind];

    SA_expup_minv150[ind] *= x_sec_minv150[ind];
    SA_exp_minv150[ind]   *= x_sec_minv150[ind];
    SA_expdn_minv150[ind] *= x_sec_minv150[ind];

    SA_expup_minv200[ind] *= x_sec_minv200[ind];
    SA_exp_minv200[ind]   *= x_sec_minv200[ind];
    SA_expdn_minv200[ind] *= x_sec_minv200[ind];
  }

  getExclusionMass(x, x_sec_minv10,  SA_expup_minv10,  limits_expup_SA_minv10);
  getExclusionMass(x, x_sec_minv50,  SA_expup_minv50,  limits_expup_SA_minv50);
  getExclusionMass(x, x_sec_minv100, SA_expup_minv100, limits_expup_SA_minv100);
  getExclusionMass(x, x_sec_minv150, SA_expup_minv150, limits_expup_SA_minv150);
  getExclusionMass(x, x_sec_minv200, SA_expup_minv200, limits_expup_SA_minv200);

  getExclusionMass(x, x_sec_minv10,  SA_exp_minv10,  limits_exp_SA_minv10);
  getExclusionMass(x, x_sec_minv50,  SA_exp_minv50,  limits_exp_SA_minv50);
  getExclusionMass(x, x_sec_minv100, SA_exp_minv100, limits_exp_SA_minv100);
  getExclusionMass(x, x_sec_minv150, SA_exp_minv150, limits_exp_SA_minv150);
  getExclusionMass(x, x_sec_minv200, SA_exp_minv200, limits_exp_SA_minv200);

  getExclusionMass(x, x_sec_minv10,  SA_expdn_minv10,  limits_expdn_SA_minv10);
  getExclusionMass(x, x_sec_minv50,  SA_expdn_minv50,  limits_expdn_SA_minv50);
  getExclusionMass(x, x_sec_minv100, SA_expdn_minv100, limits_expdn_SA_minv100);
  getExclusionMass(x, x_sec_minv150, SA_expdn_minv150, limits_expdn_SA_minv150);
  getExclusionMass(x, x_sec_minv200, SA_expdn_minv200, limits_expdn_SA_minv200);


  double excluded_SA_expup[nInv]   = { 0, 0, 0, 0, 0};
  double excluded_SA_exp[nInv]     = { 0, 0, 0, 0, 0};
  double excluded_SA_expdn[nInv]   = { 0, 0, 0, 0, 0};


//--------------------------------------------------------
//------------------- Cut&Count Part  --------------------
//--------------------------------------------------------

  //---------------------------------------------------------------------
  //for invariant mass of 10 GeV
  Double_t CC_expup_minv10[nRes]  = {0.0638, 0.1010, 0.2510, 0.6753, 2.2142, 5.462, 15.67, 40.34 , 112.6  };
  Double_t CC_exp_minv10[nRes]    = {0.0407, 0.0649, 0.1610, 0.4316, 1.3935, 3.496, 10.12, 25.63 , 71.53  };
  Double_t CC_expdn_minv10[nRes]  = {0.0274, 0.0441, 0.1080, 0.2912, 0.9266, 2.339, 6.995, 17.08 , 47.55  };

  Double_t limits_expup_CC_minv10[2]   = {0,0};
  Double_t limits_exp_CC_minv10[2]     = {0,0};
  Double_t limits_expdn_CC_minv10[2]   = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 50 GeV
  Double_t CC_expup_minv50[nRes]  = {0.0635, 0.1054, 0.2652, 0.6868, 2.1178, 5.464, 16.16, 41.14 , 110.5  };
  Double_t CC_exp_minv50[nRes]    = {0.0401, 0.0671, 0.1672, 0.4360, 1.3685, 3.456, 10.17, 26.26 , 70.49  };
  Double_t CC_expdn_minv50[nRes]  = {0.0269, 0.0450, 0.1129, 0.2930, 0.9213, 2.319, 7.060, 17.62 , 47.03  };

  Double_t limits_expup_CC_minv50[2]   = {0,0};
  Double_t limits_exp_CC_minv50[2]     = {0,0};
  Double_t limits_expdn_CC_minv50[2]   = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 100 GeV
  Double_t CC_expup_minv100[nRes]  = {0.0791, 0.1074, 0.2410, 0.7109, 3.2158, 5.452, 15.92, 44.86 , 109.9  };
  Double_t CC_exp_minv100[nRes]    = {0.0506, 0.0695, 0.1574, 0.4534, 2.0359, 3.549, 10.19, 28.12 , 70.54  };
  Double_t CC_expdn_minv100[nRes]  = {0.0337, 0.0461, 0.1056, 0.3029, 1.3664, 2.401, 7.114, 18.63 , 46.71  };

  Double_t limits_expup_CC_minv100[2]  = {0,0};
  Double_t limits_exp_CC_minv100[2]    = {0,0};
  Double_t limits_expdn_CC_minv100[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 150 GeV
  Double_t CC_expup_minv150[nRes]  = {0.1069, 0.1141, 0.2520, 0.6718, 2.1578, 5.261, 15.03, 41.89 , 114.2  };
  Double_t CC_exp_minv150[nRes]    = {0.0688, 0.0729, 0.1624, 0.4299, 1.3945, 3.424, 9.808, 27.04 , 72.83  };
  Double_t CC_expdn_minv150[nRes]  = {0.0461, 0.0494, 0.1089, 0.2843, 0.9224, 2.275, 6.828, 17.93 , 48.66  };

  Double_t limits_expup_CC_minv150[2]  = {0,0};
  Double_t limits_exp_CC_minv150[2]    = {0,0};
  Double_t limits_expdn_CC_minv150[2]  = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 200 GeV
  Double_t CC_expup_minv200[nRes]  = {0.2163, 0.1266, 0.2549, 0.7048, 2.1258, 5.228, 15.15, 38.61 , 123.2  };
  Double_t CC_exp_minv200[nRes]    = {0.1370, 0.0826, 0.1634, 0.4525, 1.3635, 3.397, 9.700, 24.87 , 80.19  };
  Double_t CC_expdn_minv200[nRes]  = {0.0910, 0.0549, 0.1108, 0.2977, 0.9005, 2.285, 6.679, 16.70 , 54.85  };

  Double_t limits_expup_CC_minv200[2]  = {0,0};
  Double_t limits_exp_CC_minv200[2]    = {0,0};
  Double_t limits_expdn_CC_minv200[2]  = {0,0};


  for(unsigned short int ind = 0; ind < nRes; ind++)
  {
    CC_expup_minv10[ind] *= x_sec_minv10[ind];
    CC_exp_minv10[ind]   *= x_sec_minv10[ind];
    CC_expdn_minv10[ind] *= x_sec_minv10[ind];

    CC_expup_minv50[ind] *= x_sec_minv50[ind];
    CC_exp_minv50[ind]   *= x_sec_minv50[ind];
    CC_expdn_minv50[ind] *= x_sec_minv50[ind];

    CC_expup_minv100[ind] *= x_sec_minv100[ind];
    CC_exp_minv100[ind]   *= x_sec_minv100[ind];
    CC_expdn_minv100[ind] *= x_sec_minv100[ind];

    CC_expup_minv150[ind] *= x_sec_minv150[ind];
    CC_exp_minv150[ind]   *= x_sec_minv150[ind];
    CC_expdn_minv150[ind] *= x_sec_minv150[ind];

    CC_expup_minv200[ind] *= x_sec_minv200[ind];
    CC_exp_minv200[ind]   *= x_sec_minv200[ind];
    CC_expdn_minv200[ind] *= x_sec_minv200[ind];
  }

  getExclusionMass(x, x_sec_minv10,  CC_expup_minv10,  limits_expup_CC_minv10);
  getExclusionMass(x, x_sec_minv50,  CC_expup_minv50,  limits_expup_CC_minv50);
  getExclusionMass(x, x_sec_minv100, CC_expup_minv100, limits_expup_CC_minv100);
  getExclusionMass(x, x_sec_minv150, CC_expup_minv150, limits_expup_CC_minv150);
  getExclusionMass(x, x_sec_minv200, CC_expup_minv200, limits_expup_CC_minv200);

  getExclusionMass(x, x_sec_minv10,  CC_exp_minv10,  limits_exp_CC_minv10);
  getExclusionMass(x, x_sec_minv50,  CC_exp_minv50,  limits_exp_CC_minv50);
  getExclusionMass(x, x_sec_minv100, CC_exp_minv100, limits_exp_CC_minv100);
  getExclusionMass(x, x_sec_minv150, CC_exp_minv150, limits_exp_CC_minv150);
  getExclusionMass(x, x_sec_minv200, CC_exp_minv200, limits_exp_CC_minv200);

  getExclusionMass(x, x_sec_minv10,  CC_expdn_minv10,  limits_expdn_CC_minv10);
  getExclusionMass(x, x_sec_minv50,  CC_expdn_minv50,  limits_expdn_CC_minv50);
  getExclusionMass(x, x_sec_minv100, CC_expdn_minv100, limits_expdn_CC_minv100);
  getExclusionMass(x, x_sec_minv150, CC_expdn_minv150, limits_expdn_CC_minv150);
  getExclusionMass(x, x_sec_minv200, CC_expdn_minv200, limits_expdn_CC_minv200);


  double excluded_CC_expup[nInv]   = { 0, 0, 0, 0, 0};
  double excluded_CC_exp[nInv]     = { 0, 0, 0, 0, 0};
  double excluded_CC_expdn[nInv]   = { 0, 0, 0, 0, 0};





//--------------------------------------------------------
//-------------------- Common Part  ----------------------
//--------------------------------------------------------


  double inv_mass[nInv] = {10, 50, 100, 150, 200};

  for(unsigned int i=0; i<nInv; i++)
  {
    if(i==0)
    {
      excluded_SA_expup[0]  = limits_expup_SA_minv10[0];
      excluded_SA_exp[0]    = limits_exp_SA_minv10[0];
      excluded_SA_expdn[0]  = limits_expdn_SA_minv10[0];

      excluded_CC_expup[0]  = limits_expup_CC_minv10[0];
      excluded_CC_exp[0]    = limits_exp_CC_minv10[0];
      excluded_CC_expdn[0]  = limits_expdn_CC_minv10[0];
    }
    if(i==1)
    {
      excluded_SA_expup[1]  = limits_expup_SA_minv50[0];
      excluded_SA_exp[1]    = limits_exp_SA_minv50[0];
      excluded_SA_expdn[1]  = limits_expdn_SA_minv50[0];

      excluded_CC_expup[1]  = limits_expup_CC_minv50[0];
      excluded_CC_exp[1]    = limits_exp_CC_minv50[0];
      excluded_CC_expdn[1]  = limits_expdn_CC_minv50[0];
    }
    if(i==2)
    {
      excluded_SA_expup[2]  = limits_expup_SA_minv100[0];
      excluded_SA_exp[2]    = limits_exp_SA_minv100[0];
      excluded_SA_expdn[2]  = limits_expdn_SA_minv100[0];

      excluded_CC_expup[2]  = limits_expup_CC_minv100[0];
      excluded_CC_exp[2]    = limits_exp_CC_minv100[0];
      excluded_CC_expdn[2]  = limits_expdn_CC_minv100[0];
    }
    if(i==3)
    {
      excluded_SA_expup[3]  = limits_expup_SA_minv150[0];
      excluded_SA_exp[3]    = limits_exp_SA_minv150[0];
      excluded_SA_expdn[3]  = limits_expdn_SA_minv150[0];

      excluded_CC_expup[3]  = limits_expup_CC_minv150[0];
      excluded_CC_exp[3]    = limits_exp_CC_minv150[0];
      excluded_CC_expdn[3]  = limits_expdn_CC_minv150[0];
    }
    if(i==4)
    {
      excluded_SA_expup[4]  = limits_expup_SA_minv200[0];
      excluded_SA_exp[4]    = limits_exp_SA_minv200[0];
      excluded_SA_expdn[4]  = limits_expdn_SA_minv200[0];

      excluded_CC_expup[4]  = limits_expup_CC_minv200[0];
      excluded_CC_exp[4]    = limits_exp_CC_minv200[0];
      excluded_CC_expdn[4]  = limits_expdn_CC_minv200[0];
    }
  }

  Double_t y_SA_sigmaup[nInv], y_SA_sigmadn[nInv], y_CC_sigmaup[nInv], y_CC_sigmadn[nInv];

  for(unsigned int i = 0; i<nInv; i++)
  {
     y_SA_sigmadn[i]  = excluded_SA_exp[i] - excluded_SA_expdn[i];
     y_SA_sigmaup[i]  = excluded_SA_expup[i] - excluded_SA_exp[i];
     y_CC_sigmadn[i]  = excluded_CC_exp[i] - excluded_CC_expdn[i];
     y_CC_sigmaup[i]  = excluded_CC_expup[i] - excluded_CC_exp[i];
  }


  TGraphErrors* limits_SA_exp     = new TGraphErrors(nInv, inv_mass, excluded_SA_exp, 0, 0);
  limits_SA_exp->SetLineWidth(2);
  limits_SA_exp->SetLineStyle(2);
  limits_SA_exp->SetMarkerStyle(22);

  TGraphErrors* limits_CC_exp     = new TGraphErrors(nInv, inv_mass, excluded_CC_exp, 0, 0);
  limits_CC_exp->SetLineWidth(2);
  limits_CC_exp->SetLineStyle(2);
  limits_CC_exp->SetMarkerStyle(20);

  TCanvas *c1 = new TCanvas();
  c1->cd();

  double zeros[nInv] = {0, 0, 0, 0, 0};

  TGraphAsymmErrors* sigma_CC = new TGraphAsymmErrors(nInv, inv_mass, excluded_CC_exp, 0, 0 , y_CC_sigmadn , y_CC_sigmaup);
  sigma_CC->SetFillColor(kYellow);
  sigma_CC->SetLineColor(kYellow);

  TGraphAsymmErrors* sigma_SA = new TGraphAsymmErrors(nInv, inv_mass, excluded_SA_exp, 0, 0 , y_SA_sigmadn, y_SA_sigmaup);
  sigma_SA->SetFillColor(kGreen);
  sigma_SA->SetLineColor(kGreen);

  for(unsigned int i=0; i<nInv; i++)  cout << "mass limit "<< inv_mass[i] << "  "  << "  exp " << excluded_SA_exp[i] << " ["<< excluded_SA_expup[i]<< "," << excluded_SA_expdn[i] << "]"  << endl;

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(sigma_CC);
  mg->Add(sigma_SA);
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
  mg->GetXaxis()->SetRangeUser(inv_mass[0],inv_mass[nInv-1]);
  mg->SetMaximum(2100);
  mg->SetMinimum(1100);


  limits_SA_exp ->Draw("samelp");
  limits_CC_exp ->Draw("samelp");


  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry(limits_SA_exp    ,"Shape Ana. CMS"    ,"lp");
  leg->AddEntry(sigma_SA  ," Expected #pm1#sigma (Shape Ana.)","f");
  leg->AddEntry(limits_CC_exp ,"Cut&Count ATLAS Sel." ,"lp");
  leg->AddEntry(sigma_CC  ," Expected #pm1#sigma (Cut&Count)","f");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();

  TLatex* text = new TLatex(0.17, 0.92, "CMS Preliminary, L = 19.7 fb^{-1}, #sqrt{s} = 8 TeV, a=0.1");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw();

  TPaveText* text_2 = new TPaveText(0.2, 0.65, 0.45, 0.75, "NDC");
  text_2->AddText("Resonant model");
  text_2->AddText("Expected 95\% CL limits");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(18);
  text_2->Draw();

  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_RES.png");
  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_RES.pdf");
  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_RES.eps");

}
