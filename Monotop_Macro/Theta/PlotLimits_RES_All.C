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

  int n=7;   // number of mass points

void getExclusionMass(Double_t * x,  Double_t* x_sec, Double_t* y_obs, Double_t* limits) {


  TGraphErrors* gr_obs = new TGraphErrors(n, x, y_obs, 0, 0);

  TGraphErrors* gr_pred = new TGraphErrors(n, x, x_sec, 0, 0);

  gr_obs->Draw();
  gr_pred->Draw("same");


  int ncross_obs = 0;

  for(unsigned int i=1; i<(n-1); i++){
    if( x_sec[i] > y_obs[i] && x_sec[i+1] < y_obs[i+1] ) {
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

  cout << "a1 " << a1  << "  b1 " << b1 << endl;
  cout << "a2 " << a2  << "  b2 " << b2 << endl;

  float x_val = (b2-b1)/(a1-a2);
  float y_val = a1*x_val+b1;

  cout << "x_val " << x_val << "  y_val " << y_val << endl;
  limits[0] = x_val;
  limits[1] = y_val;



  delete line_obs, line_sec;

}


void PlotLimits_RES_All() {

  Double_t x[7]     = {300, 500, 700, 900, 1100, 1300, 1500};

  //---------------------------------------------------------------------
  //for invariant mass of 10 GeV
  Double_t x_sec_minv10[7] = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t y_obs_minv10[7] = {0.0120 ,0.0065 ,0.01155,0.02190,0.0515 ,0.1443 ,0.3871 };
  Double_t y_exp2up_minv10[7] = {0.109  ,0.0536 ,0.0625 ,0.1136 ,0.2643 ,0.745  ,1.928  };
  Double_t y_expup_minv10[7]  = {0.060  ,0.0264 ,0.0312 ,0.0556 ,0.1374 ,0.349  ,0.876  };
  Double_t y_exp_minv10[7]    = {0.0275 ,0.0113 ,0.0145 ,0.0269 ,0.0619 ,0.1577 ,0.3771 };
  Double_t y_expdn_minv10[7]  = {0.013  ,0.0054 ,0.0084 ,0.0149 ,0.0341 ,0.0901 ,0.2111 };
  Double_t y_exp2dn_minv10[7] = {0.0077 ,0.0029 ,0.0053 ,0.010  ,0.0225 ,0.0622 ,0.1482 };

  Double_t limits_obs_minv10[2] = {0,0};
  Double_t limits_exp2up_minv10[2] = {0,0};
  Double_t limits_expup_minv10[2] = {0,0};
  Double_t limits_exp_minv10[2] = {0,0};
  Double_t limits_expdn_minv10[2] = {0,0};
  Double_t limits_exp2dn_minv10[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 50 GeV
  Double_t x_sec_minv50[7]  = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t y_obs_minv50[7]  = {0.0581 ,0.0151 ,0.02465,0.0469 ,0.0991 ,0.3223 ,0.6681 };
  Double_t y_exp2up_minv50[7] = {0.354  ,0.0975 ,0.1339 ,0.2555 ,0.4813 ,1.544  ,3.195  };
  Double_t y_expup_minv50[7]  = {0.200  ,0.0490 ,0.0670 ,0.1208 ,0.2164 ,0.727  ,1.616  };
  Double_t y_exp_minv50[7]    = {0.100  ,0.0213 ,0.0323 ,0.0559 ,0.1059 ,0.3411 ,0.7371 };
  Double_t y_expdn_minv50[7]  = {0.053  ,0.0097 ,0.0175 ,0.0312 ,0.0604 ,0.1851 ,0.4141 };
  Double_t y_exp2dn_minv50[7] = {0.032  ,0.0054 ,0.0112 ,0.019  ,0.0415 ,0.1287 ,0.2772 };

  Double_t limits_obs_minv50[2] = {0,0};
  Double_t limits_exp2up_minv50[2] = {0,0};
  Double_t limits_expup_minv50[2] = {0,0};
  Double_t limits_exp_minv50[2] = {0,0};
  Double_t limits_expdn_minv50[2] = {0,0};
  Double_t limits_exp2dn_minv50[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 100 GeV
  Double_t x_sec_minv100[7]    = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t y_obs_minv100[7]    = {0.0724 ,0.0075 ,0.01178,0.02330,0.0523 ,0.3973 ,0.3651 };
  Double_t y_exp2up_minv100[7] = {0.635  ,0.0540 ,0.0577 ,0.1145 ,0.2243 ,1.700  ,1.585  };
  Double_t y_expup_minv100[7]  = {0.346  ,0.0270 ,0.0337 ,0.0558 ,0.1154 ,0.818  ,0.786  };
  Double_t y_exp_minv100[7]    = {0.167  ,0.0123 ,0.0167 ,0.0258 ,0.0559 ,0.3991 ,0.3741 };
  Double_t y_expdn_minv100[7]  = {0.063  ,0.0058 ,0.0091 ,0.0142 ,0.0329 ,0.2261 ,0.2111 };
  Double_t y_exp2dn_minv100[7] = {0.043  ,0.0036 ,0.0059 ,0.009  ,0.0225 ,0.1557 ,0.1442 };

  Double_t limits_obs_minv100[2] = {0,0};
  Double_t limits_exp2up_minv100[2] = {0,0};
  Double_t limits_expup_minv100[2] = {0,0};
  Double_t limits_exp_minv100[2] = {0,0};
  Double_t limits_expdn_minv100[2] = {0,0};
  Double_t limits_exp2dn_minv100[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 150 GeV
  Double_t x_sec_minv150[7]    = {0, 11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t y_obs_minv150[7]    = {0, 0.0210 ,0.02665,0.05250,0.1133 ,0.3173 ,0.7121 };
  Double_t y_exp2up_minv150[7] = {0, 0.1335 ,0.1279 ,0.2515 ,0.5543 ,1.381  ,3.485  };
  Double_t y_expup_minv150[7]  = {0, 0.0690 ,0.0650 ,0.1298 ,0.2624 ,0.649  ,1.599  };
  Double_t y_exp_minv150[7]    = {0, 0.0286 ,0.0307 ,0.0629 ,0.1299 ,0.3211 ,0.7341 };
  Double_t y_expdn_minv150[7]  = {0, 0.0130 ,0.0177 ,0.0349 ,0.0705 ,0.1831 ,0.4241 };
  Double_t y_exp2dn_minv150[7] = {0, 0.0078 ,0.0110 ,0.023  ,0.0478 ,0.1237 ,0.2852 };

  Double_t limits_obs_minv150[2] = {0,0};
  Double_t limits_exp2up_minv150[2] = {0,0};
  Double_t limits_expup_minv150[2] = {0,0};
  Double_t limits_exp_minv150[2] = {0,0};
  Double_t limits_expdn_minv150[2] = {0,0};
  Double_t limits_exp2dn_minv150[2] = {0,0};

  //---------------------------------------------------------------------
  //for invariant mass of 200 GeV
  Double_t x_sec_minv200[7]    = {0, 11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t y_obs_minv200[7]    = {0, 0.0262 ,0.03165,0.04720,0.1173 ,0.3203 ,0.6971 };
  Double_t y_exp2up_minv200[7] = {0,0.1915 ,0.2009 ,0.2265 ,0.4923 ,1.428  ,3.525  };
  Double_t y_expup_minv200[7]  = {0,0.1030 ,0.1030 ,0.1168 ,0.2584 ,0.700  ,1.677  };
  Double_t y_exp_minv200[7]    = {0,0.0470 ,0.0485 ,0.0552 ,0.1209 ,0.3341 ,0.7411 };
  Double_t y_expdn_minv200[7]  = {0,0.0232 ,0.0254 ,0.0306 ,0.0704 ,0.1931 ,0.4271 };
  Double_t y_exp2dn_minv200[7] = {0,0.0144 ,0.0162 ,0.020  ,0.0479 ,0.1257 ,0.2992 };

  Double_t limits_obs_minv200[2] = {0,0};
  Double_t limits_exp2up_minv200[2] = {0,0};
  Double_t limits_expup_minv200[2] = {0,0};
  Double_t limits_exp_minv200[2] = {0,0};
  Double_t limits_expdn_minv200[2] = {0,0};
  Double_t limits_exp2dn_minv200[2] = {0,0};




  for(unsigned short int ind = 0; ind < n; ind++){
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


  for(unsigned short int ii = 0; ii < n; ii++)
  {
    x_sec_minv10[ii]  = x_sec_minv10[ii]/5.;
    x_sec_minv50[ii]  = x_sec_minv50[ii]/5. ;
    x_sec_minv100[ii] = x_sec_minv100[ii]/5. ;
    x_sec_minv150[ii] = x_sec_minv150[ii]/5. ;
    x_sec_minv200[ii] = x_sec_minv200[ii]/5. ;
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



  double excluded_obs[5] ={ 0, 0, 0, 0, 0};
  double excluded_exp2up[5] ={ 0, 0, 0, 0, 0};
  double excluded_expup[5] ={ 0, 0, 0, 0, 0};
  double excluded_exp[5] ={ 0, 0, 0, 0, 0};
  double excluded_expdn[5] ={ 0, 0, 0, 0, 0};
  double excluded_exp2dn[5] ={ 0, 0, 0, 0, 0};

  double inv_mass[5] ={ 10, 50, 100, 150, 200};




  for(int i=0; i<5; i++){
    if(i==0){
      excluded_obs[0]    = limits_obs_minv10[0];
      excluded_exp2up[0] = limits_exp2up_minv10[0];
      excluded_expup[0]  = limits_expup_minv10[0];
      excluded_exp[0]    = limits_exp_minv10[0];
      excluded_expdn[0]  = limits_expdn_minv10[0];
      excluded_exp2dn[0] = limits_exp2dn_minv10[0];
    }
    if(i==1){
      excluded_obs[1]    = limits_obs_minv50[0];
      excluded_exp2up[1] = limits_exp2up_minv50[0];
      excluded_expup[1]  = limits_expup_minv50[0];
      excluded_exp[1]    = limits_exp_minv50[0];
      excluded_expdn[1]  = limits_expdn_minv50[0];
      excluded_exp2dn[1] = limits_exp2dn_minv50[0];
    }
    if(i==2){
      excluded_obs[2]    = limits_obs_minv100[0];
      excluded_exp2up[2] = limits_exp2up_minv100[0];
      excluded_expup[2]  = limits_expup_minv100[0];
      excluded_exp[2]    = limits_exp_minv100[0];
      excluded_expdn[2]  = limits_expdn_minv100[0];
      excluded_exp2dn[2] = limits_exp2dn_minv100[0];
    }
    if(i==3){
      excluded_obs[3]    = limits_obs_minv150[0];
      excluded_exp2up[3] = limits_exp2up_minv150[0];
      excluded_expup[3]  = limits_expup_minv150[0];
      excluded_exp[3]    = limits_exp_minv150[0];
      excluded_expdn[3]  = limits_expdn_minv150[0];
      excluded_exp2dn[3] = limits_exp2dn_minv150[0];
    }
    if(i==4){
      cout << limits_obs_minv200[0] << endl;
      excluded_obs[4]    = limits_obs_minv200[0];
      cout << "excluded_obs[4]  " << excluded_obs[4] << endl;
      excluded_exp2up[4]  = limits_exp2up_minv200[0];
      excluded_expup[4]   = limits_expup_minv200[0];
      excluded_exp[4]     = limits_exp_minv200[0];
      excluded_expdn[4]   = limits_expdn_minv200[0];
      excluded_exp2dn[4]  = limits_exp2dn_minv200[0];
    }
  }

  for(int i=0; i<5; i++) cout << "mass limit obs "<< i << "  "  << excluded_obs[i] << endl;
  //for(int i=0; i<5; i++) cout << "mass limit exp "<< i << "  "  << excluded_exp[i] << endl;
  //for(int i=0; i<5; i++) cout << "mass limit exp2up "<< i << "  "  << excluded_exp2up[i] << endl;
  for(int i=0; i<5; i++) cout << "mass limit expup "<< i << "  "  << excluded_expup[i] << endl;
  for(int i=0; i<5; i++) cout << "mass limit expdn "<< i << "  "  << excluded_expdn[i] << endl;
  //for(int i=0; i<5; i++) cout << "mass limit exp2dn "<< i << "  "  << excluded_exp2dn[i] << endl;


  TGraphErrors* limits_obs     = new TGraphErrors(5, inv_mass, excluded_obs, 0, 0);
  limits_obs->SetDrawOption(0);
  limits_obs->SetLineWidth(2);
  limits_obs->SetMarkerStyle(20);

  TGraphErrors* limits_exp     = new TGraphErrors(5, inv_mass, excluded_exp, 0, 0);
  limits_exp->SetLineWidth(2);
  limits_exp->SetLineStyle(2);
  limits_exp->SetMarkerStyle(22);

  TCanvas *c1 = new TCanvas();
  c1->cd();


  TH1F * histo = new TH1F("histo", "histo", 100, 0, 210);
  histo->SetMaximum(1500);
  histo->SetMinimum(1000);
  histo->SetTitle("");
  histo->SetStats(0);

  histo->GetXaxis()->SetTitleFont(62);
  histo->GetYaxis()->SetTitleFont(62);
  histo->GetXaxis()->SetTitleOffset(.8);
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetXaxis()->SetTitleSize(0.055);
  histo->GetYaxis()->SetTitleSize(0.055);
  histo->GetXaxis()->SetLabelSize(0.04);
  histo->GetYaxis()->SetLabelSize(0.04);
  histo->GetXaxis()->SetTitle("m_{Inv.} (GeV)");
  histo->GetYaxis()->SetTitle("exclusion limit on m_{Res.} (GeV)");

  histo->Draw();
  limits_obs->Draw("samelp");
  limits_exp ->Draw("samelp");


  TGraphErrors* limits_exp2up  = new TGraphErrors(5, inv_mass, excluded_exp2up, 0, 0);
  limits_exp2up->SetLineColor(kBlue);
  limits_exp2up->SetLineWidth(2);
  limits_exp2up->SetLineStyle(4);


  TGraphErrors* limits_expup   = new TGraphErrors(5, inv_mass, excluded_expup, 0, 0);
  limits_expup->SetLineColor(kRed);
  limits_expup->SetLineWidth(2);
  limits_expup->SetLineStyle(2);

  TGraphErrors* limits_expdn   = new TGraphErrors(5, inv_mass, excluded_expdn, 0, 0);
  limits_expdn->SetLineColor(kRed);
  limits_expdn->SetLineWidth(2);
  limits_expdn->SetLineStyle(2);


  TGraphErrors* limits_exp2dn  = new TGraphErrors(5, inv_mass, excluded_exp2dn, 0, 0);
  limits_exp2dn->SetLineColor(kBlue);
  limits_exp2dn->SetLineWidth(2);
  limits_exp2dn->SetLineStyle(4);


  double zeros[5] = {0, 0, 0, 0, 0};

  TGraphAsymmErrors* sigma2 = new TGraphAsymmErrors(5, inv_mass, excluded_exp, 0, 0 , excluded_exp2up , excluded_exp2dn);
  sigma2->SetFillColor(kYellow);
  sigma2->SetLineColor(kYellow);

  //TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(n, inv_mass, excluded_exp, zeros, zeros, excluded_expup, excluded_expdn);
  TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(5, inv_mass, excluded_exp, 0, 0 , excluded_expup, excluded_expdn);
  sigma1->SetFillColor(kGreen);
  sigma1->SetLineColor(kGreen);

  for(int i=0; i<5; i++)  cout << "mass limit "<< inv_mass[i] << "  "  << excluded_obs[i] << "  exp " << excluded_exp[i] << " ["<< excluded_expup[i]<< "," << excluded_expdn[i] << "]"
  << " ["<< excluded_exp2up[i]<< "," << excluded_exp2dn[i] << "]"  << endl;

  limits_obs->Draw("samelp");
  limits_exp ->Draw("samelp");
  limits_expup->Draw("samelp");
  limits_expdn->Draw("samelp");
  //limits_exp2up->Draw("samelp");
  //limits_exp2dn->Draw("samelp");

  /*TMultiGraph *mg = new TMultiGraph();
  //mg->Add(sigma2);
  mg->Add(sigma1);
  mg->Draw("3AC");
  mg->GetXaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetTitleOffset(.8);
  mg->GetYaxis()->SetTitleOffset(1.2);
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleSize(0.055);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetXaxis()->SetTitle("m_{DM} (GeV)");
  mg->GetYaxis()->SetTitle("excludsion limite on m_{res} (GeV)");
  mg->GetXaxis()->SetRangeUser(0,210);


  limits_obs->Draw("samelp");
  limits_exp ->Draw("samelp");

  */
  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry( limits_obs ," Observed 95\% CL limit","lpe");
  leg->AddEntry( limits_exp ," Expected 95\% CL limit","lp");
  leg->AddEntry( limits_expup ," Expected #pm 1 #sigma ","lp");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw();


  TLatex* text = new TLatex(0.17, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV, a=0.05");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw();




}
