{

  const unsigned int n=9;   // number of mass points

  Double_t x[n] = {500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100};
  //Double_t x[n] = {300, 500, 700, 900, 1100, 1300, 1500};

  bool useATLASvsCMS = false;
  bool displayExclXsec = true;

  Double_t x_sec[n]       = {11.8  ,2.77, 0.87, 0.31, 0.11, 0.05, 0.022, 0.0097, 0.0045};
  //Double_t x_sec[n]       = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};

  if(!useATLASvsCMS){
     TString name = "ShapeAnalysis";

  Double_t y_obs_err[n]= {0.00006, 0.00006, 0.00011, 0.00029, 0.0012, 0.0017, 0.0053, 0.0162, 0.015 };
  Double_t y_obs[n]    = {0.0169,  0.0221,  0.0382,  0.0855,  0.3648, 0.560,  1.510,  4.110,  11.09 };
  Double_t y_exp2up[n] = {0.0323,  0.0506,  0.0809,  0.1841,  0.7841, 1.154,  3.214,  8.852,  23.76 };
  Double_t y_expup[n]  = {0.0198,  0.0303,  0.0493,  0.1139,  0.5028, 0.706,  1.978,  5.364,  14.02 };
  Double_t y_exp[n]    = {0.0125,  0.0192,  0.0315,  0.0722,  0.3169, 0.462,  1.258,  3.352,  8.878 };
  Double_t y_expdn[n]  = {0.0082,  0.0127,  0.0211,  0.0487,  0.2084, 0.318,  0.841,  2.275,  5.963 };
  Double_t y_exp2dn[n] = {0.0056,  0.0089,  0.0152,  0.0358,  0.1491, 0.239,  0.638,  1.701,  4.452 };

/*
  Double_t y_obs_err[n]= {0.00008, 0.00012, 0.00018, 0.00049, 0.0011, 0.0042, 0.0066, 0.0098, 0.015 };
  Double_t y_obs[n]    = {0.0138, 0.0203, 0.0362, 0.0820, 0.3518, 0.562, 1.537, 4.080 , 11.095 };
  Double_t y_exp2up[n] = {0.0339, 0.0481, 0.0812, 0.1841, 0.7671, 1.244, 3.424, 8.832 , 23.76  };
  Double_t y_expup[n]  = {0.0216, 0.0303, 0.0510, 0.1129, 0.4738, 0.757, 2.038, 5.254 , 14.02  };
  Double_t y_exp[n]    = {0.0136, 0.0193, 0.0327, 0.0724, 0.2969, 0.478, 1.288, 3.282 , 8.878  };
  Double_t y_expdn[n]  = {0.0087, 0.0128, 0.0218, 0.0489, 0.2004, 0.323, 0.859, 2.215 , 5.963  };
  Double_t y_exp2dn[n] = {0.0059, 0.0091, 0.0159, 0.0361, 0.1491, 0.242, 0.649, 1.661 , 4.452  };
*/
     /*
     Double_t y_obs[n]    = {0.0120 ,0.0065 ,0.01155,0.02190,0.0515 ,0.1443 ,0.3871 };
     Double_t y_obs_err[n]= {0.00007,0.00006,0.00006,0.00010,0.00024,0.0007 ,0.0022 };
     Double_t y_exp2up[n] = {0.109  ,0.0536 ,0.0625 ,0.1136 ,0.2643 ,0.745  ,1.928  };
     Double_t y_expup[n]  = {0.060  ,0.0264 ,0.0312 ,0.0556 ,0.1374 ,0.349  ,0.876  };
     Double_t y_exp[n]    = {0.0275 ,0.0113 ,0.0145 ,0.0269 ,0.0619 ,0.1577 ,0.3771 };
     Double_t y_expdn[n]  = {0.013  ,0.0054 ,0.0084 ,0.0149 ,0.0341 ,0.0901 ,0.2111 };
     Double_t y_exp2dn[n] = {0.0077 ,0.0029 ,0.0053 ,0.010  ,0.0225 ,0.0622 ,0.1482 };
*/

/*
     // with the 4 distributions in SR and mTW only in CRs
     Double_t y_obs[n]    = {0.0522 ,0.0081 ,0.00918,0.01566,0.0315 ,0.0792 };
     Double_t y_obs_err[n]= {0.00038,0.00010,0.00010,0.00019,0.00025,0.0006 };
     Double_t y_exp2up[n] = {0.257  ,0.0507 ,0.0508 ,0.0656 ,0.1547 ,0.354  };
     Double_t y_expup[n]  = {0.143  ,0.0279 ,0.0329 ,0.0422 ,0.0971 ,0.219  };
     Double_t y_exp[n]    = {0.063  ,0.0113 ,0.0163 ,0.0249 ,0.0565 ,0.1274 };
     Double_t y_expdn[n]  = {0.031  ,0.0048 ,0.0083 ,0.0141 ,0.0330 ,0.0764 };
     Double_t y_exp2dn[n] = {0.016  ,0.0029 ,0.0052 ,0.0088 ,0.0217 ,0.0475 };
*/
   }else{
  // Delta phi
     TString name = "ATLASvsCMS";
     Double_t y_obs[n]    = {0.236  ,0.367  ,0.405  };
     Double_t y_obs_err[n]= {0.00196,0.00477,0.00454};
     Double_t y_exp2up[n] = {1.79   ,1.4    ,1.31   };
     Double_t y_expup[n]  = {0.755  ,0.808  ,0.833  };
     Double_t y_exp[n]    = {0.329  ,0.447  ,0.474  };
     Double_t y_expdn[n]  = {0.166  ,0.254  ,0.295  };
     Double_t y_exp2dn[n] = {0.105  ,0.166  ,0.19   };

  }

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
  mg->GetXaxis()->SetTitle("m_{Res.} (GeV)");
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

  TLine line (x[0], 1., x[n-1], 1.);
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
  text_2->AddText("m(Inv.) = 100 GeV");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(22);
  text_2->Draw();

  TString             Y_axis = "SignalStrength";
  if(displayExclXsec) Y_axis = "ExcludedXsection";

  c1->SaveAs("limitPlots/finalLimits_RES_DM100_"+name+"_"+Y_axis+".png");
  c1->SaveAs("limitPlots/finalLimits_RES_DM100_"+name+"_"+Y_axis+".pdf");
  c1->SaveAs("limitPlots/finalLimits_RES_DM100_"+name+"_"+Y_axis+".eps");

}
