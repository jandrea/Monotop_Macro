{

  const unsigned int n=10;   // number of mass points

  Double_t x[n] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};

     TString name = "CCvsSA";
     //Double_t CC_NoSel_obs[n]       = {,,,,,,,,,};
     //Double_t CC_Sel_obs[n]         = {,,,,,,,,,};
     //Double_t CC_ATLASSel_obs[n]    = {0.26,0.74,1.01,1.95,3.28,5.34,9.67,13.1,27.1,31.9};
     //Double_t SA_SR_obs[n]          = {,,,,,,,,,};
     Double_t SA_AllCR_obs[n]       = {0.038,0.087,0.097,0.111,0.104,0.108,0.161,0.23,0.37,0.53};
     //Double_t SA_ATLASSR_obs[n]     = {,,,,,,,,,};
     Double_t CC_ATLAS_obs[n]       = {0.030,0.043,0.063,0.090,0.13,0.17,0.23,0.30,0.40,0.50};

  TCanvas *c1 = new TCanvas();
  c1->SetLogy(1);
  c1->cd();
/*
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
*/
  TGraph* SA_AllCR = new TGraph(n, x, SA_AllCR_obs);
  SA_AllCR->SetLineWidth(2);
  SA_AllCR->SetLineColor(kRed);
  SA_AllCR->SetLineStyle(1);
/*
  TGraph* SA_ATLASSR = new TGraph(n, x, SA_ATLASSR_obs);
  SA_ATLASSR->SetLineWidth(2);
  SA_ATLASSR->SetLineColor(kGreen);
  SA_ATLASSR->SetLineStyle(1);
*/
  TGraph* CC_ATLAS = new TGraph(n, x, CC_ATLAS_obs);
  CC_ATLAS->SetLineWidth(2);
  CC_ATLAS->SetLineColor(kGreen);
  CC_ATLAS->SetLineStyle(1);

 // CC_NoSel->Draw();
 // CC_Sel->Draw("same");
  CC_ATLAS->Draw();
  CC_ATLAS->SetTitle("");
  CC_ATLAS->GetXaxis()->SetTitle("m_{DM} (GeV)");
  CC_ATLAS->GetXaxis()->SetTitleSize(0.05);
  CC_ATLAS->GetXaxis()->SetTitleOffset(0.9);
  CC_ATLAS->GetYaxis()->SetTitle("95% CL limit on signal strength");
  CC_ATLAS->GetXaxis()->SetRangeUser(x[0],x[n-1]);
  CC_ATLAS->SetMinimum(0.005);
  CC_ATLAS->SetMaximum(100);
 // SA_SR->Draw("same");
  SA_AllCR->Draw("same");
 // SA_ATLASSR->Draw("same");

  TLine line (x[0], 1.,x[n-1], 1.);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->SetLineStyle(2);
  line->Draw("same");

  TLegend *leg = new TLegend(0.2,0.65,0.45,0.89);
  //leg->AddEntry(SA_SR       ,"Shape Ana. SR-only"   ,"l");
  //leg->AddEntry(SA_AllCR    ,"Shape Ana. All-CR"    ,"l");
  leg->AddEntry(SA_AllCR    ,"Shape Ana. CMS"    ,"l");
  //leg->AddEntry(SA_ATLASSR  ,"Shape Ana. ATLAS-SR"  ,"l");
  //leg->AddEntry(CC_NoSel    ,"Cut&Count PreSel"     ,"l");
  //leg->AddEntry(CC_Sel      ,"Cut&Count CMS Sel."   ,"l");
  //leg->AddEntry(CC_ATLASSel ,"Cut&Count ATLAS Sel." ,"l");
  leg->AddEntry(CC_ATLAS ,"Cut&Count ATLAS" ,"l");
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->Draw("same");


  TLatex* text = new TLatex(0.37, 0.92, "CMS Preliminary, L = 20 fb^{-1}, #sqrt{s} = 8 TeV");
  //text->SetTextAlign(33);  // 33-left, 22-center, 11-right
  text->SetNDC(true);
  text->SetTextFont(62);
  text->SetTextSize(0.045);
  text->Draw("same");

  c1->SaveAs("limitPlots/limits"+name+"_FCNC.png");
  c1->SaveAs("limitPlots/limits"+name+"_FCNC.pdf");
  c1->SaveAs("limitPlots/limits"+name+"_FCNC.eps");

}
