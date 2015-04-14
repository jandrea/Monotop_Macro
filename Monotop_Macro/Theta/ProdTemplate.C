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


void ProdTemplate(TString inputdistrib, TString outputdistrib, std::vector<TString> signalList, std::vector<TString> thetaSignalList, std::vector<TString> sampleList, std::vector<TString> stytList, TString intputfilename, std::vector<double> scaleCMStoATLAS, bool rescaleToATLAS, bool doCutnCount ){

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
      if(sampleList[i] == "TTMSDecays")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "TTMSDecays_central", inputfile);
      }
      else if(sampleList[i] == "WExcll" ||sampleList[i] == "WExclc" || sampleList[i] == "WExclb" || sampleList[i] == "QCD")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , sampleList[i], inputfile);
      }
     else if(sampleList[i] == "SingleTop")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_s", inputfile);
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_t", inputfile));
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_tW", inputfile));
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_t", inputfile));
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_tW", inputfile));
      }
      else if(sampleList[i] == "DY")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "DYJetsToLL_M-10To50", inputfile);
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "DYJetsToLL_M-50", inputfile));
      }
      else if(sampleList[i] == "VV")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "WW", inputfile);
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "ZZ", inputfile));
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "WZ", inputfile));
      }
      else
      {
                 cout << "WARNING: sample: " << sampleList[i] << " not treated. Please check the config file." << endl;
      }
      distrib__nominale->SetName((outputdistrib+"__"+sampleList[i]).Data());
      distrib_MC.push_back( (TH1D*)distrib__nominale);
  }
  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_highpt") continue;

      for(unsigned int j=0; j<stytList.size(); j++)
      {
          if( stytList[j]=="toppt" && sampleList[i] != "TTMSDecays") continue;
          if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching") && (sampleList[i] != "WExcll" && sampleList[i] != "WExclb" && sampleList[i] != "WExclc" && sampleList[i] != "TTMSDecays") ) continue;
          if(  stytList[j]=="mTWtail" && (sampleList[i] != "WExclc" && sampleList[i] != "WExclb" && sampleList[i] != "WExcll" )) continue;
          if( (stytList[j] =="Iso" || stytList[j] == "BgdContam") && sampleList[i] != "QCD" ) continue;
          if( sampleList[i] == "QCD" && (stytList[j] !="Iso" && stytList[j] != "BgdContam") ) continue;

          if (i < signalList.size() && stytList[j] != "mass" && stytList[j] != "scale" && stytList[j] != "matching" && stytList[j] != "toppt" && stytList[j]!="mTWtail" && stytList[j] != "Iso" && stytList[j] != "BgdContam" && stytList[j] != "PDF")
          {
              TH1D* distribsignal__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              TH1D* distribsignal__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              //if(rescaleToATLAS) tmp_inputsignaldistribminus->Scale(scaleCMStoATLAS[i]);
              //if(rescaleToATLAS) tmp_inputsignaldistribplus->Scale(scaleCMStoATLAS[i]);
              distribsignal__plus->SetName( (outputdistrib+"__"+thetaSignalList[i]+"__"+stytList[j]+"__plus").Data());
              distribsignal__minus->SetName((outputdistrib+"__"+thetaSignalList[i]+"__"+stytList[j]+"__minus").Data());
              distrib_signal_sys.push_back( (TH1D*)distribsignal__plus );
              distrib_signal_sys.push_back( (TH1D*)distribsignal__minus);
          }

          TH1D* distribmc__plus = 0;
          TH1D* distribmc__minus = 0;
          TH1D* distribWtail__nominale = 0;
          TH1D* distribWtail__up = 0;
          TH1D* distribWtail__down = 0;

          if( (sampleList[i] == "TTMSDecays" || sampleList[i] == "WExcll" || sampleList[i] == "WExclc" || sampleList[i] == "WExclb")  && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
          //if( (sampleList[i] == "TTMSDecays" )  && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
          {
              TString inputdistribnameplus;
              TString inputdistribnameminus;
              if(stytList[j] == "mass")
              {
                  inputdistribnameplus    = inputdistrib+"__"+sampleList[i]+"_mass173_5";
                  inputdistribnameminus   = inputdistrib+"__"+sampleList[i]+"_mass171_5";
              }
              else if(stytList[j] == "scale")
              {
                  inputdistribnameplus    = inputdistrib+"__"+sampleList[i]+"_scaleup";
                  inputdistribnameminus   = inputdistrib+"__"+sampleList[i]+"_scaledown";
              }
              else if(stytList[j] == "matching")
              {
                  inputdistribnameplus    = inputdistrib+"__"+sampleList[i]+"_matchingup";
                  inputdistribnameminus   = inputdistrib+"__"+sampleList[i]+"_matchingdown";
              }

              TH1D* tmp_inputdistribplus      = (TH1D*)inputfile->Get(inputdistribnameplus)->Clone() ;
              TH1D* tmp_inputdistribminus     = (TH1D*)inputfile->Get(inputdistribnameminus)->Clone() ;
              TString outputdistribnameplus   = outputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__plus";
              TString outputdistribnameminus  = outputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__minus";
              tmp_inputdistribminus->SetName(outputdistribnameminus );
              tmp_inputdistribplus->SetName(outputdistribnameplus   );
              distrib_MC_sys.push_back( (TH1D*)tmp_inputdistribplus );
              distrib_MC_sys.push_back( (TH1D*)tmp_inputdistribminus);
          }
          else
          {
                if(sampleList[i] == "TTMSDecays")
                {
                    distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "TTMSDecays_central", inputfile);
                    distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "TTMSDecays_central", inputfile);
                }
                else if(sampleList[i] == "WExcll" ||sampleList[i] == "WExclc" || sampleList[i] == "WExclb" || sampleList[i] == "QCD")
                {
                    distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
                    distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
                }
                else if(sampleList[i] == "SingleTop")
                {
                    distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_s", inputfile);
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_t", inputfile));
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile));
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                    distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_s", inputfile);
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_t", inputfile));
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile));
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                }
                else if(sampleList[i] == "DY")
                {
                    distribmc__plus           = getSystematic("plus"    , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-10To50", inputfile);
                    distribmc__plus->Add(   	getSystematic("plus"    , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-50", inputfile));

                    distribmc__minus          = getSystematic("minus"   , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-10To50", inputfile);
                    distribmc__minus->Add(   	getSystematic("minus"   , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-50", inputfile));
                }
                else if(sampleList[i] == "VV")
                {
                    distribmc__plus         	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "WW", inputfile);
                    distribmc__plus->Add(   	  getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "ZZ", inputfile));
                    distribmc__plus->Add(   	  getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "WZ", inputfile));

                    distribmc__minus         	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "WW", inputfile);
                    distribmc__minus->Add(   	  getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "ZZ", inputfile));
                    distribmc__minus->Add(   	  getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "WZ", inputfile));
                }
                else
                {
                    cout << "WARNING: sample: " << sampleList[i] << " not treated. Please check the config file." << endl;
                }

                distribmc__plus->SetName( (outputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__plus").Data());
                distribmc__minus->SetName((outputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__minus").Data());
                distrib_MC_sys.push_back( (TH1D*)distribmc__plus );
                distrib_MC_sys.push_back( (TH1D*)distribmc__minus);

          }
      }
  }

  TString outputfilename;
  if      (inputdistrib == "mWT_mujets_Wregion_highpt")         outputfilename = "inputTheta_Wregion";
  else if (inputdistrib == "mWT_mujets_ttbarregion_highpt")     outputfilename = "inputTheta_ttbarregion";
  else if (inputdistrib == "mWT_mujets_Selectedsignalregion")   outputfilename = "inputTheta_Selectedsignalregion_mWT";
  else if (inputdistrib == "MET_mujets_Selectedsignalregion")   outputfilename = "inputTheta_Selectedsignalregion_MET";
  else if (inputdistrib == "DeltaPhiLJ_mujets_signalregion")    outputfilename = "inputTheta_signalregion_DeltaPhiLJ";
  else                                                          outputfilename = "inputTheta_"+inputdistrib;

  outputfilename+=".root";

  TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );

  if(doCutnCount)  distrib__DATA->Rebin(distrib__DATA->GetNbinsX());

  outputfile->cd();
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


void ProdTemplate(){

  bool rescaleToATLAS       = false;
  bool doCutnCount          = false;

  std::vector<TString> signalList;
  std::vector<TString> thetaSignalList;
  std::vector<double>  scaleCMStoATLAS;
  signalList.push_back("S4_500_fastSim");  thetaSignalList.push_back("S4inv500");
  scaleCMStoATLAS.push_back(1.001/11.79);
  std::vector<TString> sampleList;
  sampleList.push_back("TTMSDecays" );
  sampleList.push_back("WExclb"         );
  sampleList.push_back("WExclc"         );
  sampleList.push_back("WExcll"         );
  sampleList.push_back("DY"             );
  sampleList.push_back("SingleTop"      );
  sampleList.push_back("VV"             );
  sampleList.push_back("QCD"            );

  std::vector<TString> systlist;

  systlist.push_back("lept"             );
  systlist.push_back("trig"             );
  systlist.push_back("PDF"              );
  systlist.push_back("PU"               );
  systlist.push_back("jes"              );
  systlist.push_back("jer"              );
  systlist.push_back("metuncls"         );

  //systlist.push_back("mass"             );
  systlist.push_back("scale"            );
  systlist.push_back("matching"         );
  systlist.push_back("toppt"            );

  systlist.push_back("btag"             );
  systlist.push_back("mistag"           );

  systlist.push_back("Iso"              );
  systlist.push_back("BgdContam"        );


  // ------------------------------------------------------- //
  // "signalregion" means no cut on (MET, DeltaPhiLJ)        //
  // "Selectedsignalregion" means (MET>100, DeltaPhiLJ<1.8)  //
  // ------------------------------------------------------- //


  //ProdTemplate("mWT_mujets_ATLASRESsignalregion", "mWTmujetsATLASRESSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("mWT_mujets_ATLASFCNCsignalregion", "mWTmujetsATLASFCNCSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_Selectedsignalregion", "mWTmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("MET_mujets_Selectedsignalregion", "METmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("DeltaPhiLJ_mujets_Selectedsignalregion", "DeltaPhiLJmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("ptW_mujets_Selectedsignalregion", "ptWmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

  //ProdTemplate("MET_mujets_signalregion", "METmujetsSignalregion", signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("DeltaPhiLJ_mujets_signalregion", "DeltaPhiLJmujetsSignalregion", signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("mWT_mujets_signalregion", "mWTmujetsSignalregion", signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt",signalList, thetaSignalList, sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

}
