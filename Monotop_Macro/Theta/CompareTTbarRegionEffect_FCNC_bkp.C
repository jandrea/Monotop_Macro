{

  const unsigned int n=10;   // number of mass points

  Double_t x[n] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

     TString name = "TTregionComparison";
     //Double_t CC_NoSel_exp[n]       = {,,,,,,,,,};
     //Double_t CC_Sel_exp[n]         = {,,,,,,,,,};
     //Double_t CC_ATLASSel_exp[n]    = {0.26,0.74,1.01,1.95,3.28,5.34,9.67,13.1,27.1,31.9};
     //Double_t SA_SR_exp[n]          = {,,,,,,,,,};
     //Double_t SA_mTWonly_exp[n]       = {0.0697  , 0.1990 , 0.2098 , 0.373  , 0.423  , 0.7118, 1.102 , 1.860 , 2.87 , 4.07  };  // with all CR
     Double_t SA_tt2j2b_expup[n]     = {0.148   , 0.33   , 0.46   , 0.631  , 1.14   , 1.60  , 2.141 , 3.311 , 5.614, 8.22  };
     Double_t SA_tt2j2b_exp[n]       = {0.081   , 0.185  , 0.267  , 0.393  , 0.674  , 0.975 , 1.342 , 2.083 , 3.435, 4.94  };
     Double_t SA_tt2j2b_expdn[n]     = {0.0485  , 0.113  , 0.168  , 0.254  , 0.429  , 0.637 , 0.889 , 1.372 , 2.241, 3.24  };
     Double_t SA_tt3j2b_exp[n]       = {0.084   , 0.189  , 0.282  , 0.411  , 0.692  , 1.015 , 1.402 , 2.123 , 3.515, 5.14  };
     Double_t SA_tt4j2b_exp[n]       = {0.073   , 0.299  , 0.356  , 0.603  , 0.707  , 1.212 , 1.812 , 2.733 , 4.235, 5.45  };
     //Double_t SA_AllCR_exp[n]       = {0.038,0.087,0.097,0.111,0.104,0.108,0.161,0.23,0.37,0.53};  // with all CR
     //Double_t SA_ATLASSR_exp[n]     = {,,,,,,,,,};
     //Double_t CC_ATLAS_exp[n]       = {0.030,0.043,0.063,0.090,0.13,0.17,0.23,0.30,0.40,0.50};   //ex limit
     //Double_t CC_ATLAS_exp[n]         = {0.1773,0.403,0.637,1.1790,1.427,3.07,4.108,6.65,11.58,14.49};


     Double_t x_sec[n]                = {28.0  ,8.1  ,2.94, 1.22, 0.57, 0.29, 0.15, 0.08, 0.05, 0.03};

     for(unsigned short int ind = 0; ind < n; ind++)
     {
          SA_tt2j2b_expup[ind]    *= x_sec[ind];
          SA_tt2j2b_exp[ind]    *= x_sec[ind];
          SA_tt2j2b_expdn[ind]    *= x_sec[ind];
          SA_tt3j2b_exp[ind]    *= x_sec[ind];
          SA_tt4j2b_exp[ind]    *= x_sec[ind];
     }

  TCanvas *c1 = new TCanvas();
  c1->SetLogy(1);
  c1->cd();


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

  TGraph* gr_SA_tt2j2b_exp = new TGraph(n, x, SA_tt2j2b_exp);
  gr_SA_tt2j2b_exp->SetLineWidth(2);
  gr_SA_tt2j2b_exp->SetLineColor(kBlack);
  gr_SA_tt2j2b_exp->SetLineStyle(1);

  TGraph* gr_SA_tt3j2b_exp = new TGraph(n, x, SA_tt3j2b_exp);
  gr_SA_tt3j2b_exp->SetLineWidth(2);
  gr_SA_tt3j2b_exp->SetLineColor(kGreen);
  gr_SA_tt3j2b_exp->SetLineStyle(1);

  TGraph* gr_SA_tt4j2b_exp = new TGraph(n, x, SA_tt4j2b_exp);
  gr_SA_tt4j2b_exp->SetLineWidth(2);
  gr_SA_tt4j2b_exp->SetLineColor(kOrange);
  gr_SA_tt4j2b_exp->SetLineStyle(1);

  Double_t y_sigmaup[n], y_sigmadn[n];

  for(int i = 0; i<n; i++)
  {
     y_sigmadn[i]  = SA_tt2j2b_exp[i] - SA_tt2j2b_expdn[i];
     y_sigmaup[i]  = SA_tt2j2b_expup[i] - SA_tt2j2b_exp[i];
  }

  TGraphAsymmErrors* sigma1 = new TGraphAsymmErrors(n, x, SA_tt2j2b_exp, 0, 0, y_sigmadn, y_sigmaup);
  sigma1->SetFillColor(kGreen);
  sigma1->SetLineColor(kGreen);




  gr_SA_tt2j2b_exp->Draw();
  gr_SA_tt2j2b_exp->SetTitle("");
  gr_SA_tt2j2b_exp->GetXaxis()->SetTitleFont(62);
  gr_SA_tt2j2b_exp->GetYaxis()->SetTitleFont(62);
  gr_SA_tt2j2b_exp->GetXaxis()->SetTitleSize(0.055);
  gr_SA_tt2j2b_exp->GetYaxis()->SetTitleSize(0.045);
  gr_SA_tt2j2b_exp->GetXaxis()->SetLabelSize(0.04);
  gr_SA_tt2j2b_exp->GetYaxis()->SetLabelSize(0.04);
  gr_SA_tt2j2b_exp->GetXaxis()->SetTitle("m_{Inv.} (GeV)");
  gr_SA_tt2j2b_exp->GetXaxis()->SetTitleOffset(0.8);
  gr_SA_tt2j2b_exp->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
  gr_SA_tt2j2b_exp->GetYaxis()->SetTitleOffset(0.7);
  gr_SA_tt2j2b_exp->SetMinimum(SA_tt2j2b_exp[n-1]);
  gr_SA_tt2j2b_exp->SetMaximum(x_sec[0]);
  gr_SA_tt2j2b_exp->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  sigma1->Draw("same");
  gr_SA_tt3j2b_exp->Draw("same");
  gr_SA_tt4j2b_exp->Draw("same");

/*
  TLine line (x[0], 1.,x[n-1], 1.);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw("same");
*/
  xsec_th_0p1->Draw("same");
  xsec_th_0p05->Draw("same");


  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  leg->AddEntry(gr_SA_tt2j2b_exp    ,"TT-CR ( == 2j , == 2b )"    ,"l");
  leg->AddEntry(gr_SA_tt3j2b_exp    ,"TT-CR ( == 3j , == 2b )"    ,"l");
  leg->AddEntry(gr_SA_tt4j2b_exp    ,"TT-CR ( >= 4j , == 2b )"    ,"l");
  leg->AddEntry(xsec_th_0p1  ," Theory (LO), a_{res} = 0.1","l");
  leg->AddEntry(xsec_th_0p05  ," Theory (LO), a_{res} = 0.05","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("same");

  TPaveText* text_2;
  text_2 = new TPaveText(0.6, 0.45, 0.85, 0.55, "NDC");
  text_2->AddText("Non-resonant model");
  text_2->SetFillColor(kWhite);
  text_2->SetFillStyle(0);
  text_2->SetBorderSize(0);
  text_2->SetTextFont(43);
  text_2->SetTextSize(22);
  text_2->Draw();



  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw("same");

  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_FCNC.png");
  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_FCNC.pdf");
  c1->SaveAs("limitPlots/finalLimits_methodComparison_"+name+"_FCNC.eps");

}
