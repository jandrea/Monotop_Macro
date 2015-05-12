{

  const unsigned int n=6;   // number of mass points

  Double_t x[n] = {500, 700, 900, 1100, 1300, 1500};
  //Double_t x[n] = {300, 500, 700, 900, 1100, 1300, 1500};

     TString name = "CCvsSA";
     //Double_t CC_NoSel_exp[n]       = {,,,,,,,,,};
     //Double_t CC_Sel_exp[n]         = {,,,,,,,,,};
     //Double_t CC_ATLASSel_exp[n]    = {0.26,0.74,1.01,1.95,3.28,5.34,9.67,13.1,27.1,31.9};
     //Double_t SA_SR_exp[n]          = {,,,,,,,,,};
     Double_t SA_mTWonly_exp[n]       = {0.0165 ,0.0210 ,0.0313 ,0.0698 ,0.4731 ,0.4621 };
     //Double_t SA_mTWonly_exp[n]       = {0.0581 ,0.0156 ,0.02165,0.03860,0.0783 ,0.2003 ,0.4891 };
     //Double_t SA_AllCR_exp[n]       = {0.038,0.087,0.097,0.111,0.104,0.108,0.161,0.23,0.37,0.53};  // with all CR
     //Double_t SA_ATLASSR_exp[n]     = {,,,,,,,,,};
     //Double_t CC_ATLAS_exp[n]       = {0.030,0.043,0.063,0.090,0.13,0.17,0.23,0.30,0.40,0.50};   //ex limit
     Double_t CC_ATLAS_exp[n]         = {0.057 ,0.0718 ,0.1640 ,0.486 ,3.358 ,3.59};
     //Double_t CC_ATLAS_exp[n]         = {100 ,0.050 ,0.0655 ,0.1540 ,0.414 ,1.263 ,3.41};


     Double_t x_sec[n]       = {11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};
     //Double_t x_sec[n]       = {77.2  ,11.8  ,2.77, 0.87, 0.31, 0.11, 0.05};

     for(unsigned short int ind = 0; ind < n; ind++)
     {
          SA_mTWonly_exp[ind]    *= x_sec[ind];
          CC_ATLAS_exp[ind]      *= x_sec[ind];
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

/*
  TGraph* CC_NoSel = new TGraph(n, x, CC_NoSel_exp);
  CC_NoSel->SetTitle("");
  CC_NoSel->SetLineWidth(2);
  CC_NoSel->SetLineColor(kBlue);
  CC_NoSel->SetLineStyle(2);

  TGraph* CC_Sel = new TGraph(n, x, CC_Sel_exp);
  CC_Sel->SetLineWidth(2);
  CC_Sel->SetLineColor(kBlack);
  CC_Sel->SetLineStyle(2);

  TGraph* CC_ATLASSel = new TGraph(n, x, CC_ATLASSel_exp);
  CC_ATLASSel->SetLineWidth(2);
  CC_ATLASSel->SetLineColor(kGreen);
  CC_ATLASSel->SetLineStyle(2);

  TGraph* SA_SR = new TGraph(n, x, SA_SR_exp);
  SA_SR->SetTitle("");
  SA_SR->SetLineWidth(2);
  SA_SR->SetLineColor(kBlack);
  SA_SR->SetLineStyle(1);
*/
  TGraph* SA_mTWonly = new TGraph(n, x, SA_mTWonly_exp);
  SA_mTWonly->SetLineWidth(2);
  SA_mTWonly->SetLineColor(kBlack);
  SA_mTWonly->SetLineStyle(1);
/*
  TGraph* SA_AllCR = new TGraph(n, x, SA_AllCR_exp);
  SA_AllCR->SetLineWidth(2);
  SA_AllCR->SetLineColor(kRed);
  SA_AllCR->SetLineStyle(1);

  TGraph* SA_ATLASSR = new TGraph(n, x, SA_ATLASSR_exp);
  SA_ATLASSR->SetLineWidth(2);
  SA_ATLASSR->SetLineColor(kGreen);
  SA_ATLASSR->SetLineStyle(1);
*/
  TGraph* CC_ATLAS = new TGraph(n, x, CC_ATLAS_exp);
  CC_ATLAS->SetLineWidth(2);
  CC_ATLAS->SetLineColor(kGreen);
  CC_ATLAS->SetLineStyle(1);

 // CC_NoSel->Draw();
 // CC_Sel->Draw("same");
  CC_ATLAS->Draw();
  CC_ATLAS->SetTitle("");
  CC_ATLAS->GetXaxis()->SetTitleFont(62);
  CC_ATLAS->GetYaxis()->SetTitleFont(62);
  CC_ATLAS->GetXaxis()->SetTitleSize(0.055);
  CC_ATLAS->GetYaxis()->SetTitleSize(0.045);
  CC_ATLAS->GetXaxis()->SetLabelSize(0.04);
  CC_ATLAS->GetYaxis()->SetLabelSize(0.04);
  CC_ATLAS->GetXaxis()->SetTitle("m_{mediator} (GeV)");
  CC_ATLAS->GetXaxis()->SetTitleOffset(0.8);
  CC_ATLAS->GetYaxis()->SetTitle("#sigma(p p #rightarrow t v_{met} #times BR(t #rightarrow b l #nu) [pb]   ");
  CC_ATLAS->GetYaxis()->SetTitleOffset(0.7);
  CC_ATLAS->SetMinimum(0.008);
  CC_ATLAS->SetMaximum(x_sec[0]);
  //CC_ATLAS->GetYaxis()->SetTitle("95% CL limit on signal strength");
  CC_ATLAS->GetXaxis()->SetRangeUser(x[0],x[n-1]);
 // SA_SR->Draw("same");
  SA_mTWonly->Draw("same");
 // SA_mTWonlyAllCR->Draw("same");
 // SA_ATLASSR->Draw("same");

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
  //leg->AddEntry(SA_SR       ,"Shape Ana. SR-only"   ,"l");
  //leg->AddEntry(SA_AllCR    ,"Shape Ana. All-CR"    ,"l");
  leg->AddEntry(SA_mTWonly    ,"Shape Ana. CMS"    ,"l");
  //leg->AddEntry(SA_AllCR    ,"Shape Ana. CMS"    ,"l");
  //leg->AddEntry(SA_ATLASSR  ,"Shape Ana. ATLAS-SR"  ,"l");
  //leg->AddEntry(CC_NoSel    ,"Cut&Count PreSel"     ,"l");
  //leg->AddEntry(CC_Sel      ,"Cut&Count CMS Sel."   ,"l");
  //leg->AddEntry(CC_ATLASSel ,"Cut&Count ATLAS Sel." ,"l");
  leg->AddEntry(CC_ATLAS ,"Cut&Count ATLAS" ,"l");
  leg->AddEntry(xsec_th_0p1  ," Theory (LO), a_{res} = 0.1","l");
  leg->AddEntry(xsec_th_0p05  ," Theory (LO), a_{res} = 0.05","l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("same");

  TPaveText* text_2;
  text_2 = new TPaveText(0.6, 0.45, 0.85, 0.55, "NDC");
  text_2->AddText("Resonant model");
  text_2->AddText("m(DM) = 100 GeV");
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

  c1->SaveAs("limitPlots/new_limits_comparison_"+name+"_RES.png");
  c1->SaveAs("limitPlots/new_limits_comparison_"+name+"_RES.pdf");
  c1->SaveAs("limitPlots/new_limits_comparison_"+name+"_RES.eps");

}
