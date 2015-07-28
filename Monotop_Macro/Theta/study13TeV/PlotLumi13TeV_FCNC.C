{

  const unsigned int n=10;   // number of mass points

  Double_t x[n]           = {100,   200,  300,  400,  500,  600,  700,  800,  900,  1000};
  Double_t x_sec[n]       = {599,   174,  7.33, 3.55, 1.74, 0.98, 0.58, 0.35, 0.22, 0.15};

     //Double_t lumi0p1_exp[n]       = {0.055,  0.139, 2.006, 3.403, 4.227, 7.11, 10.38, 15.5, 23.15, 31.3   };
     //Double_t lumi0p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     Double_t lumi1p0_exp[n]       = {0.0194, 0.052, 0.734, 1.073, 1.227, 1.98, 2.734, 4.13, 6.067, 7.413  };
     //Double_t lumi1p5_exp[n]       = {0.0163, 0.051, 0.613, 0.973, 1.057, 1.67, 2.324, 3.42, 4.947, 6.35   };
     Double_t lumi2p0_exp[n]       = {0.0153, 0.048, 0.540, 0.813, 0.944, 1.55, 2.094, 2.92, 4.317, 5.08  };
     //Double_t lumi2p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     //Double_t lumi3p0_exp[n]       = {0.0133, 0.043, 0.497, 0.733, 0.837, 1.37, 1.834, 2.39, 3.667, 4.43  };
     //Double_t lumi3p5_exp[n]       = {0.0263, 0.071, 0.983, 1.373, 1.807, 3.14, 4.224, 5.93, 8.667, 11.43  };
     Double_t lumi4p0_exp[n]       = {0.0113, 0.040, 0.463, 0.693, 0.687, 1.18, 1.684, 2.05, 2.977, 3.93  };
     //Double_t lumi5p0_exp[n]       = {0.0123, 0.038, 0.423, 0.633, 0.677, 1.11, 1.464, 2.05, 3.047, 3.48  };
     //Double_t lumi6p0_exp[n]       = {0.0104, 0.040, 0.401, 0.622, 0.629, 1.09, 1.414, 1.83, 2.807, 3.21  };
     //Double_t lumi7p0_exp[n]       = {0.0100, 0.039, 0.400, 0.582, 0.627, 0.96, 1.384, 1.69, 2.757, 2.99  };
     Double_t lumi8p0_exp[n]       = {0.0100, 0.038, 0.391, 0.522, 0.569, 0.90, 1.274, 1.57, 2.687, 2.84  };
     Double_t lumi16p0_exp[n]      = {0.0090, 0.036, 0.311, 0.472, 0.499, 0.66, 0.914, 1.15, 1.947, 2.06  };

     for(unsigned short int ind = 0; ind < n; ind++)
     {
          //lumi0p1_exp[ind]    *= x_sec[ind];
          //lumi0p5_exp[ind]    *= x_sec[ind];
          lumi1p0_exp[ind]    *= x_sec[ind];
          //lumi1p5_exp[ind]    *= x_sec[ind];
          lumi2p0_exp[ind]    *= x_sec[ind];
          //lumi3p0_exp[ind]    *= x_sec[ind];
          lumi4p0_exp[ind]    *= x_sec[ind];
          //lumi5p0_exp[ind]    *= x_sec[ind];
          //lumi6p0_exp[ind]    *= x_sec[ind];
          //lumi7p0_exp[ind]    *= x_sec[ind];
          lumi8p0_exp[ind]    *= x_sec[ind];
          lumi16p0_exp[ind]   *= x_sec[ind];
     }

  TCanvas *c1 = new TCanvas();
  c1->SetLogy(1);
  c1->cd();


  TGraph* xsec_th_0p1 = new TGraph(n, x, x_sec);
  xsec_th_0p1->SetLineColor(kRed);
  xsec_th_0p1->SetLineWidth(2);
  xsec_th_0p1->SetLineStyle(2);
/*
  TGraph* SA_lumi0p1 = new TGraph(n, x, lumi0p1_exp);
  SA_lumi0p1->SetLineWidth(2);
  SA_lumi0p1->SetLineColor(kBlack);
  SA_lumi0p1->SetLineStyle(1);

  TGraph* SA_lumi0p5 = new TGraph(n, x, lumi0p5_exp);
  SA_lumi0p5->SetLineWidth(2);
  SA_lumi0p5->SetLineColor(kGreen);
  SA_lumi0p5->SetLineStyle(1);
*/
  TGraph* SA_lumi1p0 = new TGraph(n, x, lumi1p0_exp);
  SA_lumi1p0->SetLineWidth(2);
  SA_lumi1p0->SetLineColor(kBlue);
  SA_lumi1p0->SetLineStyle(1);
/*
  TGraph* SA_lumi1p5 = new TGraph(n, x, lumi1p5_exp);
  SA_lumi1p5->SetLineWidth(2);
  SA_lumi1p5->SetLineColor(kOrange-2);
  SA_lumi1p5->SetLineStyle(1);
*/
  TGraph* SA_lumi2p0 = new TGraph(n, x, lumi2p0_exp);
  SA_lumi2p0->SetLineWidth(2);
  SA_lumi2p0->SetLineColor(kGray);
  SA_lumi2p0->SetLineStyle(1);
/*
  TGraph* SA_lumi3p0 = new TGraph(n, x, lumi3p0_exp);
  SA_lumi3p0->SetLineWidth(2);
  SA_lumi3p0->SetLineColor(kRed-4);
  SA_lumi3p0->SetLineStyle(1);
*/
  TGraph* SA_lumi4p0 = new TGraph(n, x, lumi4p0_exp);
  SA_lumi4p0->SetLineWidth(2);
  SA_lumi4p0->SetLineColor(kGreen);
  SA_lumi4p0->SetLineStyle(1);
/*
  TGraph* SA_lumi5p0 = new TGraph(n, x, lumi5p0_exp);
  SA_lumi5p0->SetLineWidth(2);
  SA_lumi5p0->SetLineColor(kOrange-2);
  SA_lumi5p0->SetLineStyle(1);

  TGraph* SA_lumi6p0 = new TGraph(n, x, lumi6p0_exp);
  SA_lumi6p0->SetLineWidth(2);
  SA_lumi6p0->SetLineColor(kBlack);
  SA_lumi6p0->SetLineStyle(1);
*//*
  TGraph* SA_lumi7p0 = new TGraph(n, x, lumi7p0_exp);
  SA_lumi7p0->SetLineWidth(2);
  SA_lumi7p0->SetLineColor(kBlack);
  SA_lumi7p0->SetLineStyle(1);
*/
  TGraph* SA_lumi8p0 = new TGraph(n, x, lumi8p0_exp);
  SA_lumi8p0->SetLineWidth(2);
  SA_lumi8p0->SetLineColor(kOrange-2);
  SA_lumi8p0->SetLineStyle(1);

  TGraph* SA_lumi16p0 = new TGraph(n, x, lumi16p0_exp);
  SA_lumi16p0->SetLineWidth(2);
  SA_lumi16p0->SetLineColor(kBlack);
  SA_lumi16p0->SetLineStyle(1);



  SA_lumi1p0->Draw();
  SA_lumi1p0->SetTitle("");
  SA_lumi1p0->GetXaxis()->SetTitleFont(62);
  SA_lumi1p0->GetYaxis()->SetTitleFont(62);
  SA_lumi1p0->GetXaxis()->SetTitleSize(0.055);
  SA_lumi1p0->GetYaxis()->SetTitleSize(0.045);
  SA_lumi1p0->GetXaxis()->SetLabelSize(0.04);
  SA_lumi1p0->GetYaxis()->SetLabelSize(0.04);
  SA_lumi1p0->GetXaxis()->SetTitle("m_{DM} (GeV)");
  SA_lumi1p0->GetXaxis()->SetTitleOffset(0.8);
  SA_lumi1p0->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
  SA_lumi1p0->GetYaxis()->SetTitleOffset(0.7);
  SA_lumi1p0->SetMinimum(lumi16p0_exp[n-1]);
  SA_lumi1p0->SetMaximum(x_sec[0]);
  SA_lumi1p0->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  //SA_lumi0p1->Draw("same");
  SA_lumi1p0->Draw("same");
  //SA_lumi1p5->Draw("same");
  SA_lumi2p0->Draw("same");
  //SA_lumi3p0->Draw("same");
  SA_lumi4p0->Draw("same");
  //SA_lumi5p0->Draw("same");
  //SA_lumi6p0->Draw("same");
  //SA_lumi7p0->Draw("same");
  SA_lumi8p0->Draw("same");
  SA_lumi16p0->Draw("same");

  xsec_th_0p1->Draw("same");


  TLegend *leg = new TLegend(0.55,0.65,0.85,0.85);
  //leg->AddEntry(SA_lumi0p1,     "L = 0.1 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi0p5,     "L = 0.5 fb^{-1}" ,"l");
  leg->AddEntry(SA_lumi1p0,     "L = 1.0 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi1p5,     "L = 1.5 fb^{-1}" ,"l");
  leg->AddEntry(SA_lumi2p0,     "L = 2.0 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi3p0,     "L = 3.0 fb^{-1}" ,"l");
  leg->AddEntry(SA_lumi4p0,     "L = 4.0 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi5p0,     "L = 5.0 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi6p0,     "L = 6.0 fb^{-1}" ,"l");
  //leg->AddEntry(SA_lumi7p0,     "L = 7.0 fb^{-1}" ,"l");
  leg->AddEntry(SA_lumi8p0,     "L = 8.0 fb^{-1}" ,"l");
  leg->AddEntry(SA_lumi16p0,    "L =  16 fb^{-1}" ,"l");
  leg->AddEntry(xsec_th_0p1  ," Theory (LO), a_{res} = 0.1","l");
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



  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, #sqrt{s} = 13 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw("same");

  c1->SaveAs("Comparison_lumi13TeV_FCNC.png");
  c1->SaveAs("Comparison_lumi13TeV_FCNC.pdf");
  c1->SaveAs("Comparison_lumi13TeV_FCNC.eps");

}
