#include "TString.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>
#include <vector>


TH1D* getSystematic(TString plusORminus, TString inputdistrib, TString outputdistrib, TString syst, TString sample, TFile* inputfile)
{
    TString fullSystName;

    if (plusORminus != "" && syst != "")
    {
        fullSystName   = inputdistrib+"__"+sample+"__"+syst+"__"+plusORminus;
    }
    else
    {
        fullSystName   = inputdistrib+"__"+sample;
    }
    TH1D * distrib__sample  	= (TH1D*)inputfile->Get(fullSystName)->Clone();

    return distrib__sample;
}


void ProdTemplate(TString inputdistrib, TString outputdistrib, std::vector<TString> signalList, std::vector<TString> thetaSignalList, std::vector<TString> sampleList, std::vector<TString> thetaSampleList, std::vector<TString> stytList, std::vector<TString> thetaStytList, TString intputfilename, std::vector<double> scaleCMStoATLAS, bool rescaleToATLAS, bool doCutnCount ){

  TFile * inputfile	  = new TFile( intputfilename.Data() );

  TString dataAname         = inputdistrib+"__SingleMuA";
  TString dataBname         = inputdistrib+"__SingleMuB";
  TString dataCname         = inputdistrib+"__SingleMuC";
  TString dataDname         = inputdistrib+"__SingleMuD";

  TH1D * distrib__RUNA  	= (TH1D*)inputfile->Get(dataAname)->Clone();
  TH1D * distrib__RUNB  	= (TH1D*)inputfile->Get(dataBname)->Clone();
  TH1D * distrib__RUNC  	= (TH1D*)inputfile->Get(dataCname)->Clone();
  TH1D * distrib__RUND  	= (TH1D*)inputfile->Get(dataDname)->Clone();

  TH1D * distrib__DATA  	= (TH1D*)distrib__RUNA;
         distrib__DATA->Add(distrib__RUNB);
         distrib__DATA->Add(distrib__RUNC);
         distrib__DATA->Add(distrib__RUND);


  std::vector< TH1D* > distrib_MC;
  std::vector< TH1D* > distrib_MC_sys;
  std::vector< TH1D* > distrib_signal_sys;
  std::vector< TH1D* > distrib_signal;
  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(i <  signalList.size() )
      {
          TString signalname          = inputdistrib+"__"+signalList[i];
          TH1D * distrib__Monotop     = (TH1D*)inputfile->Get(signalname)->Clone();
          if(rescaleToATLAS) distrib__Monotop->Scale(scaleCMStoATLAS[i]);
          TString outputsignalname    = outputdistrib+"__"+thetaSignalList[i];
          distrib__Monotop->SetName(outputsignalname );
          distrib_signal.push_back( (TH1D*)distrib__Monotop);
      }
      TH1D* distrib__nominale = 0;
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_highpt") continue;
      if(sampleList[i] == "DYJetsToLL_M-10To50" && inputdistrib == "mWT_mujets_Selectedsignalregion") continue;

      distrib__nominale    = getSystematic(""     , inputdistrib, outputdistrib, ""       , sampleList[i], inputfile);
      distrib__nominale->SetName((outputdistrib+"__"+thetaSampleList[i]).Data());
      distrib_MC.push_back( (TH1D*)distrib__nominale);
  }

  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_highpt") continue;
      if(sampleList[i] == "DYJetsToLL_M-10To50" && inputdistrib == "mWT_mujets_Selectedsignalregion") continue;

      for(unsigned int j=0; j<stytList.size(); j++)
      {
          if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching" || stytList[j]=="toppt") && sampleList[i] != "TTMSDecays_central") continue;
          if( (stytList[j] =="Iso" || stytList[j] == "BgdContam") && sampleList[i] != "QCD" ) continue;
          if( sampleList[i] == "QCD" && (stytList[j] !="Iso" && stytList[j] != "BgdContam") ) continue;

          if (i < signalList.size() && stytList[j] != "mass" && stytList[j] != "scale" && stytList[j] != "matching" && stytList[j] != "toppt" && stytList[j] != "Iso" && stytList[j] != "BgdContam")
          {

	      TH1D* distribsignal__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              TH1D* distribsignal__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              distribsignal__plus->SetName( (outputdistrib+"__"+thetaSignalList[i]+"__"+thetaStytList[j]+"__plus").Data());
              distribsignal__minus->SetName((outputdistrib+"__"+thetaSignalList[i]+"__"+thetaStytList[j]+"__minus").Data());
              distrib_signal_sys.push_back( (TH1D*)distribsignal__plus );
              distrib_signal_sys.push_back( (TH1D*)distribsignal__minus);
          }

          TH1D* distribmc__plus = 0;
          TH1D* distribmc__minus = 0;


          if( sampleList[i] == "TTMSDecays_central" && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
          {
              TString inputdistribnameplus;
              TString inputdistribnameminus;
              if(stytList[j] == "mass")
              {
                  inputdistribnameplus    = inputdistrib+"__TTMSDecays_mass173_5";
                  inputdistribnameminus   = inputdistrib+"__TTMSDecays_mass171_5";
              }
              else if(stytList[j] == "scale")
              {
                  inputdistribnameplus    = inputdistrib+"__TTMSDecays_scaleup";
                  inputdistribnameminus   = inputdistrib+"__TTMSDecays_scaledown";
              }
              else if(stytList[j] == "matching")
              {
                  inputdistribnameplus    = inputdistrib+"__TTMSDecays_matchingup";
                  inputdistribnameminus   = inputdistrib+"__TTMSDecays_matchingdown";
              }
              else cout << "There is a problem at line " << __LINE__ << " : Please check TTbar systematics" << endl;

              TH1D* tmp_inputdistribplus      = (TH1D*)inputfile->Get(inputdistribnameplus)->Clone() ;
              TH1D* tmp_inputdistribminus     = (TH1D*)inputfile->Get(inputdistribnameminus)->Clone() ;
              TString outputdistribnameplus   = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__plus";
              TString outputdistribnameminus  = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__minus";
              tmp_inputdistribminus->SetName(outputdistribnameminus );
              tmp_inputdistribplus->SetName(outputdistribnameplus   );
              distrib_MC_sys.push_back( (TH1D*)tmp_inputdistribplus );
              distrib_MC_sys.push_back( (TH1D*)tmp_inputdistribminus);
          }
          else
          {
              distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
              distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
              distribmc__plus->SetName( (outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__plus").Data());
              distribmc__minus->SetName((outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__minus").Data());
              distrib_MC_sys.push_back( (TH1D*)distribmc__plus );
              distrib_MC_sys.push_back( (TH1D*)distribmc__minus);
          }
      }
  }

  TString outputfilename;
  if      (inputdistrib == "mWT_mujets_Wregion_highpt")         outputfilename = "inputTheta_splitted_Wregion";
  else if (inputdistrib == "mWT_mujets_ttbarregion_highpt")     outputfilename = "inputTheta_splitted_ttbarregion";
  else if (inputdistrib == "mWT_mujets_Selectedsignalregion")   outputfilename = "inputTheta_splitted_Selectedsignalregion";
  else                                                          outputfilename = "inputTheta_splitted_"+inputdistrib;

  outputfilename+=".root";

  TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );
  outputfile->cd();

  if(doCutnCount)   distrib__DATA->Rebin(distrib__DATA->GetNbinsX());

  distrib__DATA->Write((outputdistrib+"__DATA").Data());

  for(unsigned int i=0; i<distrib_MC.size(); i++)
  {
      if(doCutnCount) distrib_MC[i]->Rebin(distrib_MC[i]->GetNbinsX());
      distrib_MC[i]->Write();
  }
  for(unsigned int i=0; i<distrib_signal.size(); i++)
  {
      if(doCutnCount) distrib_signal[i]->Rebin(distrib_signal[i]->GetNbinsX());
      distrib_signal[i]->Write();
  }
  for(unsigned int i=0; i<distrib_signal_sys.size(); i++)
  {
      if(doCutnCount) distrib_signal_sys[i]->Rebin(distrib_signal_sys[i]->GetNbinsX());
      distrib_signal_sys[i]->Write();
  }
  for(unsigned int i=0; i<distrib_MC_sys.size(); i++)
  {
      if(doCutnCount) distrib_MC_sys[i]->Rebin(distrib_MC_sys[i]->GetNbinsX());
      distrib_MC_sys[i]->Write();
  }

}


void ProdTemplate_splittedSamples(){

  bool rescaleToATLAS       = false;
  bool doCutnCount          = false;

  std::vector<TString> signalList;
  std::vector<TString> thetaSignalList;
  std::vector<double>  scaleCMStoATLAS;
  signalList.push_back("S4Inv700");  thetaSignalList.push_back("S4Inv700");
  scaleCMStoATLAS.push_back(1.001/11.79);
  std::vector<TString> sampleList;
  std::vector<TString> thetaSampleList;
  sampleList.push_back("TTMSDecays_central" );      thetaSampleList.push_back("TTMSDecays"      );
  sampleList.push_back("WExclb"             );      thetaSampleList.push_back("WExclb"          );
  sampleList.push_back("WExclc"             );      thetaSampleList.push_back("WExclc"          );
  sampleList.push_back("WExcll"             );      thetaSampleList.push_back("WExcll"          );
  sampleList.push_back("DYJetsToLL_M-10To50");      thetaSampleList.push_back("DY10To50"        );
  sampleList.push_back("DYJetsToLL_M-50"    );      thetaSampleList.push_back("DY50"            );
  sampleList.push_back("T_s"                );      thetaSampleList.push_back("Ts"              );
  sampleList.push_back("T_t"                );      thetaSampleList.push_back("Tt"              );
  sampleList.push_back("T_tW"               );      thetaSampleList.push_back("TtW"             );
  sampleList.push_back("Tbar_t"             );      thetaSampleList.push_back("Tbart"           );
  sampleList.push_back("Tbar_tW"            );      thetaSampleList.push_back("TbartW"          );
  sampleList.push_back("WZ"                 );      thetaSampleList.push_back("WZ"              );
  sampleList.push_back("ZZ"                 );      thetaSampleList.push_back("ZZ"              );
  sampleList.push_back("WW"                 );      thetaSampleList.push_back("WW"              );
  sampleList.push_back("QCD"                );      thetaSampleList.push_back("QCD"             );

  std::vector<TString> systlist;
  std::vector<TString> thetaSystlist;

  systlist.push_back("lept"             );          thetaSystlist.push_back("lept"           );
  systlist.push_back("trig"             );          thetaSystlist.push_back("trig"           );
  //systlist.push_back("PDF"              );          thetaSystlist.push_back("PDF"            );
  systlist.push_back("PU"               );          thetaSystlist.push_back("PU"             );
  systlist.push_back("jes"              );          thetaSystlist.push_back("jes"            );
  systlist.push_back("jer"              );          thetaSystlist.push_back("jer"            );
  systlist.push_back("metuncls"         );          thetaSystlist.push_back("metuncls"       );

  //systlist.push_back("mass"             );          thetaSystlist.push_back("mass"           );
  systlist.push_back("scale"            );          thetaSystlist.push_back("scale"          );
  systlist.push_back("matching"         );          thetaSystlist.push_back("matching"       );
  systlist.push_back("toppt"            );          thetaSystlist.push_back("toppt"          );

  systlist.push_back("btag"             );          thetaSystlist.push_back("btag"           );
  systlist.push_back("mistag"           );          thetaSystlist.push_back("mistag"         );


  // ------------------------------------------------------- //
  // "signalregion" means no cut on (MET, DeltaPhiLJ)        //
  // "Selectedsignalregion" means (MET>100, DeltaPhiLJ<1.8)  //
  // ------------------------------------------------------- //


  //ProdTemplate("mWT_mujets_ATLASRESsignalregion", "mWTmujetsATLASRESSignalregion",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

  ProdTemplate("mWT_mujets_Selectedsignalregion", "mWTmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

}
