{

  const unsigned int n=10;   // number of mass points

  Double_t x[n] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

  bool useATLASvsCMS   = false;
  bool displayExclXsec = true;

  Double_t x_sec[n]       = {28.0  ,8.1  ,2.94, 1.22, 0.57, 0.29, 0.15, 0.08, 0.05, 0.03};

  if(!useATLASvsCMS){
     TString name = "ShapeAnalysis";

     Double_t y_obs[n]    = {0.0641  , 0.173  , 0.245  , 0.418  , 0.739  , 1.129 , 1.553 , 2.439 , 4.31 , 6.02  };
     Double_t y_obs_err[n]= {0.0003  , 0.0008 , 0.0009 , 0.0013 , 0.0026 , 0.0035, 0.0048, 0.0071, 0.015, 0.022 };
     Double_t y_exp2up[n] = {0.342   , 0.712  , 0.877  , 1.145  , 2.13   , 2.81  , 3.820 , 5.475 , 10.0 , 15.16 };
     Double_t y_expup[n]  = {0.148   , 0.33   , 0.46   , 0.631  , 1.14   , 1.60  , 2.141 , 3.311 , 5.614, 8.22  };
     Double_t y_exp[n]    = {0.081   , 0.185  , 0.267  , 0.393  , 0.674  , 0.975 , 1.342 , 2.083 , 3.435, 4.94  };
     Double_t y_expdn[n]  = {0.0485  , 0.113  , 0.168  , 0.254  , 0.429  , 0.637 , 0.889 , 1.372 , 2.241, 3.24  };
     Double_t y_exp2dn[n] = {0.0317  , 0.076  , 0.116  , 0.178  , 0.298  , 0.451 , 0.643 , 0.994 , 1.621, 2.27  };

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
  mg->GetXaxis()->SetTitle("m_{Inv.} (GeV)");
  mg->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  if(displayExclXsec)
  {
      mg->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
      mg->GetYaxis()->SetTitleSize(0.04);
      mg->GetYaxis()->SetTitleOffset(1.2);
      //mg->SetMinimum(0.01);
      mg->SetMinimum(y_exp2dn[n-1]);
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

  c1->SaveAs("limitPlots/finalLimits_FCNC_"+name+"_"+Y_axis+"_2j2b.png");
  c1->SaveAs("limitPlots/finalLimits_FCNC_"+name+"_"+Y_axis+"_2j2b.pdf");
  c1->SaveAs("limitPlots/finalLimits_FCNC_"+name+"_"+Y_axis+"_2j2b.eps");

}
