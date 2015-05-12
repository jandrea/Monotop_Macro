{

  const unsigned int n=6;   // number of mass points

  //Double_t x[n] = {300, 500, 700, 900, 1100, 1300, 1500};
  Double_t x[n] = {500, 700, 900, 1100, 1300, 1500};

  bool useATLASvsCMS = false;
  bool displayExclXsec = true;

  //Double_t x_sec[n]       = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
  Double_t x_sec[n]       = {11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};

  if(!useATLASvsCMS){
     TString name = "ShapeAnalysis";

     Double_t y_obs[n]    = {0.0072 ,0.01388,0.02900,0.0661 ,0.4853 ,0.4771 };
     Double_t y_obs_err[n]= {0.00005,0.00010,0.00014,0.00044,0.0030 ,0.0028 };
     Double_t y_exp2up[n] = {0.0540 ,0.0591 ,0.0875 ,0.1893 ,1.233  ,1.202  };
     Double_t y_expup[n]  = {0.0319 ,0.0346 ,0.0491 ,0.1144 ,0.742  ,0.721  };
     Double_t y_exp[n]    = {0.0165 ,0.0210 ,0.0313 ,0.0698 ,0.4731 ,0.4621 };
     Double_t y_expdn[n]  = {0.0088 ,0.0129 ,0.0205 ,0.0494 ,0.3241 ,0.3211 };
     Double_t y_exp2dn[n] = {0.0053 ,0.0087 ,0.015  ,0.0367 ,0.2467 ,0.2392 };

/*
     Double_t y_obs[n]    = {0.0724 ,0.0075 ,0.01178,0.02330,0.0523 ,0.3973 ,0.3651 };
     Double_t y_obs_err[n]= {0.00044,0.00007,0.00006,0.00013,0.00030,0.0019 ,0.0019 };
     Double_t y_exp2up[n] = {0.635  ,0.0540 ,0.0577 ,0.1145 ,0.2243 ,1.700  ,1.585  };
     Double_t y_expup[n]  = {0.346  ,0.0270 ,0.0337 ,0.0558 ,0.1154 ,0.818  ,0.786  };
     Double_t y_exp[n]    = {0.167  ,0.0123 ,0.0167 ,0.0258 ,0.0559 ,0.3991 ,0.3741 };
     Double_t y_expdn[n]  = {0.063  ,0.0058 ,0.0091 ,0.0142 ,0.0329 ,0.2261 ,0.2111 };
     Double_t y_exp2dn[n] = {0.043  ,0.0036 ,0.0059 ,0.009  ,0.0225 ,0.1557 ,0.1442 };
*/

/*
     Double_t y_obs[n]    = {0.0581 ,0.0156 ,0.02165,0.03860,0.0783 ,0.2003 ,0.4891 };
     Double_t y_obs_err[n]= {0.00062,0.00022,0.00021,0.00033,0.00060,0.0016 ,0.0032 };
     Double_t y_exp2up[n] = {0.354  ,0.0675 ,0.0699 ,0.1065 ,0.2153 ,0.554  ,1.405  };
     Double_t y_expup[n]  = {0.200  ,0.0360 ,0.0410 ,0.0638 ,0.1444 ,0.349  ,0.866  };
     Double_t y_exp[n]    = {0.100  ,0.0173 ,0.0243 ,0.0389 ,0.0829 ,0.2121 ,0.5081 };
     Double_t y_expdn[n]  = {0.053  ,0.0086 ,0.0145 ,0.0242 ,0.0504 ,0.1251 ,0.3341 };
     Double_t y_exp2dn[n] = {0.032  ,0.0054 ,0.0095 ,0.016  ,0.0335 ,0.0867 ,0.2462 };
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
  mg->GetXaxis()->SetTitle("m_{mediator} (GeV)");
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

  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw();

  TPaveText* text_2;
  if(!displayExclXsec) text_2 = new TPaveText(0.3, 0.45, 0.55, 0.55, "NDC");
  else                 text_2 = new TPaveText(0.6, 0.45, 0.85, 0.55, "NDC");
  text_2->AddText("Resonant model");
  text_2->AddText("m(DM) = 100 GeV");
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
