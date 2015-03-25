{

  const unsigned int n=10;   // number of mass points

  Double_t x[n] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

  bool useATLASvsCMS   = false;
  bool displayExclXsec = true;

  Double_t x_sec[n]       = {28.0  ,8.1  ,2.94, 1.22, 0.57, 0.29, 0.15, 0.08, 0.05, 0.03};

  if(!useATLASvsCMS){
     TString name = "ShapeAnalysis";
     Double_t y_obs[n]    = {0.0503  , 0.1346 , 0.1698 , 0.322  , 0.114  , 0.112 , 0.177 , 0.247 , 0.422, 0.61  };
     Double_t y_obs_err[n]= {0.000796, 0.00189, 0.00168, 0.00287, 0.00097, 0.0010, 0.0018, 0.0025, 0.005, 0.0077};
     Double_t y_exp2up[n] = {0.399   , 0.76   , 1.63   , 2.128  , 0.512  , 0.597 , 0.781 , 0.922 , 1.81 , 2.56  };
     Double_t y_expup[n]  = {0.1920  , 0.36   , 0.73   , 1.007  , 0.312  , 0.376 , 0.455 , 0.575 , 0.99 , 1.49  };
     Double_t y_exp[n]    = {0.0855  , 0.152  , 0.336  , 0.521  , 0.180  , 0.211 , 0.259 , 0.340 , 0.596, 0.883 };
     Double_t y_expdn[n]  = {0.0422  , 0.078  , 0.167  , 0.269  , 0.102  , 0.124 , 0.153 , 0.213 , 0.342, 0.532 };
     Double_t y_exp2dn[n] = {0.0222  , 0.044  , 0.0833 , 0.161  , 0.0611 , 0.070 , 0.096 , 0.125 , 0.203, 0.321 };


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

  TGraph* xsec_th = new TGraph(n, x, x_sec);
  xsec_th->SetLineColor(kRed);
  xsec_th->SetLineWidth(2);
  xsec_th->SetLineStyle(2);


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
  mg->GetXaxis()->SetTitle("m_{DM} (GeV)");
  mg->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  if(displayExclXsec)
  {
      mg->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
      mg->GetYaxis()->SetTitleSize(0.04);
      mg->GetYaxis()->SetTitleOffset(1.2);
      mg->SetMinimum(0.01);
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
  else                 xsec_th->Draw("same");

  TLegend *leg;
  if(!displayExclXsec) leg = new TLegend(0.3,0.65,0.55,0.85);
  else                 leg = new TLegend(0.6,0.65,0.85,0.85);
  leg->AddEntry(gr_obs  ," Observed 95\% CL limit","lpe");
  leg->AddEntry(gr_exp  ," Expected 95\% CL limit","lp");
  leg->AddEntry(sigma1  ," Expected #pm1#sigma","f");
  leg->AddEntry(sigma2  ," Expected #pm2#sigma","f");
  if(displayExclXsec) leg->AddEntry(xsec_th  ," Theory (LO), a_{non-res} = 0.1","l");
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
  text_2->AddText("Non-resonant model");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(22);
  text_2->Draw();

  TString             Y_axis = "SignalStrength";
  if(displayExclXsec) Y_axis = "ExcludedXsection";

  c1->SaveAs("limitPlots/limits_FCNC_"+name+"_"+Y_axis+".png");
  c1->SaveAs("limitPlots/limits_FCNC_"+name+"_"+Y_axis+".pdf");
  c1->SaveAs("limitPlots/limits_FCNC_"+name+"_"+Y_axis+".eps");

}
