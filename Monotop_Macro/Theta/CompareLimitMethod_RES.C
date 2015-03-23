{

  const unsigned int n=6;   // number of mass points

  Double_t x[n] = {300, 500, 700, 900, 1100, 1300};

     TString name = "CCvsSA";
     Double_t CC_NoSel_obs[n]       = {0.39 ,0.95   ,3.16   ,10.6   ,31.7   ,98.4 };
     Double_t CC_Sel_obs[n]         = {0.21 ,0.04   ,0.11   ,0.34   ,1.06   ,3.14 };
     Double_t CC_ATLASSel_obs[n]    = {100  ,0.05   ,0.07   ,0.17   ,0.49   ,1.45 };
     Double_t SA_SR_obs[n]          = {0.09 ,0.037  ,0.0429 ,0.05   ,0.10   ,0.26 };
     Double_t SA_AllCR_obs[n]       = {0.04 ,0.011  ,0.009  ,0.015  ,0.029  ,0.075};
     Double_t SA_ATLASSR_obs[n]     = {100  ,0.053  ,0.041  ,0.056  ,0.14   ,0.32 };
     //Double_t y_exp2dn[n] = {0.016  ,0.0029 ,0.0052 ,0.0088 ,0.0217 ,0.0475 };

  TCanvas *c1 = new TCanvas();
  c1->SetLogy(1);
  c1->cd();

  TGraph* CC_NoSel = new TGraph(n, x, CC_NoSel_obs);
  CC_NoSel->SetTitle("");
  CC_NoSel->SetLineWidth(2);
  CC_NoSel->SetLineColor(kBlue);
  CC_NoSel->SetLineStyle(2);

  TGraph* CC_Sel = new TGraph(n, x, CC_Sel_obs);
  CC_Sel->SetLineWidth(2);
  CC_Sel->SetLineColor(kBlack);
  CC_Sel->SetLineStyle(2);

  TGraph* CC_ATLASSel = new TGraph(n, x, CC_ATLASSel_obs);
  CC_ATLASSel->SetLineWidth(2);
  CC_ATLASSel->SetLineColor(kGreen);
  CC_ATLASSel->SetLineStyle(2);

  TGraph* SA_SR = new TGraph(n, x, SA_SR_obs);
  SA_SR->SetTitle("");
  SA_SR->SetLineWidth(2);
  SA_SR->SetLineColor(kBlack);
  SA_SR->SetLineStyle(1);

  TGraph* SA_AllCR = new TGraph(n, x, SA_AllCR_obs);
  SA_AllCR->SetLineWidth(2);
  SA_AllCR->SetLineColor(kRed);
  SA_AllCR->SetLineStyle(1);

  TGraph* SA_ATLASSR = new TGraph(n, x, SA_ATLASSR_obs);
  SA_ATLASSR->SetLineWidth(2);
  SA_ATLASSR->SetLineColor(kGreen);
  SA_ATLASSR->SetLineStyle(1);


  CC_NoSel->Draw();
  CC_Sel->Draw("same");
  CC_ATLASSel->Draw("same");
  CC_NoSel->GetXaxis()->SetTitle("m_{mediator} (GeV)");
  CC_NoSel->GetXaxis()->SetTitleSize(0.05);
  CC_NoSel->GetXaxis()->SetTitleOffset(0.9);
  CC_NoSel->GetYaxis()->SetTitle("95% CL limit on signal strength");
  CC_NoSel->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  CC_NoSel->SetMinimum(0.005);
  CC_NoSel->SetMaximum(100);
  SA_SR->Draw("same");
  SA_AllCR->Draw("same");
  SA_ATLASSR->Draw("same");

  TLine line (300, 1.,1300., 1.);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw("same");

  TLegend *leg = new TLegend(0.2,0.65,0.45,0.89);
  leg->AddEntry(SA_SR       ,"Shape Ana. SR-only"   ,"l");
  leg->AddEntry(SA_AllCR    ,"Shape Ana. All-CR"    ,"l");
  leg->AddEntry(SA_ATLASSR  ,"Shape Ana. ATLAS-SR"  ,"l");
  leg->AddEntry(CC_NoSel    ,"Cut&Count PreSel"     ,"l");
  leg->AddEntry(CC_Sel      ,"Cut&Count CMS Sel."   ,"l");
  leg->AddEntry(CC_ATLASSel ,"Cut&Count ATLAS Sel." ,"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("same");


  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw("same");

  c1->SaveAs("limitPlots/limits"+name+".png");
  c1->SaveAs("limitPlots/limits"+name+".pdf");
  c1->SaveAs("limitPlots/limits"+name+".eps");

}
