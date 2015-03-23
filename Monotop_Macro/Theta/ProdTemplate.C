#include "TString.h"
#include "TH1F.h"
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

TH1F* getQCDsystematic(TString plusORminus, TString inputdistrib, TString outputdistrib)
{
    TFile* qcdfile;
    qcdfile	  = new TFile("../TreeReader/outputroot_withSyst/histo_merged.root");
    TString QCDAname, QCDBname, QCDCname, QCDDname;

    if (plusORminus != "")
    {
        QCDAname   = inputdistrib+"__QCD_A__"+plusORminus;
        QCDBname   = inputdistrib+"__QCD_B__"+plusORminus;
        QCDCname   = inputdistrib+"__QCD_C__"+plusORminus;
        QCDDname   = inputdistrib+"__QCD_D__"+plusORminus;
    }
    else
    {
        QCDAname   = inputdistrib+"__QCD_A";
        QCDBname   = inputdistrib+"__QCD_B";
        QCDCname   = inputdistrib+"__QCD_C";
        QCDDname   = inputdistrib+"__QCD_D";
    }

    TH1F * distrib__QCDA  	= (TH1F*)qcdfile->Get(QCDAname)->Clone();
    //TH1F * distrib__QCDB  	= (TH1F*)qcdfile->Get(QCDBname)->Clone();
    //TH1F * distrib__QCDC  	= (TH1F*)qcdfile->Get(QCDCname)->Clone();
    //TH1F * distrib__QCDD  	= (TH1F*)qcdfile->Get(QCDDname)->Clone();

    TH1F * distrib__QCD  	= (TH1F*)distrib__QCDA;
    //       distrib__QCD->Add(distrib__QCDB);
    //       distrib__QCD->Add(distrib__QCDC);
    //       distrib__QCD->Add(distrib__QCDD);

    if (plusORminus != "")  distrib__QCD->SetName((outputdistrib+"__QCD__iso__"+plusORminus).Data());
    else                    distrib__QCD->SetName((outputdistrib+"__QCD").Data());

    return distrib__QCD;
}


void ProdTemplate(TString inputdistrib, TString outputdistrib, std::vector<TString> signalList, std::vector<TString> thetaSignalList, std::vector<TString> sampleList, std::vector<TString> thetaSampleList, std::vector<TString> stytList, std::vector<TString> thetaStytList, TString intputfilename, std::vector<double> scaleCMStoATLAS, bool rescaleToATLAS, bool doCutnCount ){


  TFile * inputfile	  = new TFile( intputfilename.Data() );

  TString dataAname         = inputdistrib+"__SingleMuA";
  TString dataBname         = inputdistrib+"__SingleMuB";
  TString dataCname         = inputdistrib+"__SingleMuC";
  TString dataDname         = inputdistrib+"__SingleMuD";

  TH1F * distrib__RUNA  	= (TH1F*)inputfile->Get(dataAname)->Clone();
  TH1F * distrib__RUNB  	= (TH1F*)inputfile->Get(dataBname)->Clone();
  TH1F * distrib__RUNC  	= (TH1F*)inputfile->Get(dataCname)->Clone();
  TH1F * distrib__RUND  	= (TH1F*)inputfile->Get(dataDname)->Clone();

  TH1F * distrib__DATA  	= (TH1F*)distrib__RUNA;
         distrib__DATA->Add(distrib__RUNB);
         distrib__DATA->Add(distrib__RUNC);
         distrib__DATA->Add(distrib__RUND);

  TH1F * distrib__QCD         	= getQCDsystematic(""     , inputdistrib, outputdistrib);
  //TH1F * distrib__QCD__plus  	= getQCDsystematic("plus" , inputdistrib, outputdistrib);
  //TH1F * distrib__QCD__minus  	= getQCDsystematic("minus", inputdistrib, outputdistrib);

  std::vector< TH1F* > distrib_MC;
  std::vector< TH1F* > distrib_MC_sys;
  std::vector< TH1F* > distrib_signal_sys;
  std::vector< TH1F* > distrib_signal;


  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(i <  signalList.size() )
      {
          TString signalname          = inputdistrib+"__"+signalList[i];
          TH1F * distrib__Monotop     = (TH1F*)inputfile->Get(signalname)->Clone();
          if(rescaleToATLAS) distrib__Monotop->Scale(scaleCMStoATLAS[i]);
          TString outputsignalname    = outputdistrib+"__"+thetaSignalList[i];
          distrib__Monotop->SetName(outputsignalname );
          distrib_signal.push_back( (TH1F*)distrib__Monotop);
      }
      TString inputdistribname    = inputdistrib+"__"+sampleList[i];
      TH1F* tmp_inputdistrib      = (TH1F*)inputfile->Get(inputdistribname)->Clone() ;
      TString outputdistribname   = outputdistrib+"__"+thetaSampleList[i];
      tmp_inputdistrib->SetName(outputdistribname );
      distrib_MC.push_back( (TH1F*)tmp_inputdistrib);
  }

  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      for(unsigned int j=0; j<stytList.size(); j++)
      {
          if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching" || stytList[j]=="toppt") && sampleList[i] != "TTbar_Madgraph") continue;

          if (i < signalList.size() && stytList[j] != "mass" && stytList[j] != "scale" && stytList[j] != "matching" && stytList[j] != "toppt" )
          {
              TString signalnameplus                = inputdistrib+"__"+signalList[i]+"__"+stytList[j]+"__plus";
              TString signalnameminus               = inputdistrib+"__"+signalList[i]+"__"+stytList[j]+"__minus";
              TH1F* tmp_inputsignaldistribminus     = (TH1F*)inputfile->Get(signalnameplus)->Clone() ;
              TH1F* tmp_inputsignaldistribplus      = (TH1F*)inputfile->Get(signalnameminus)->Clone() ;
              if(rescaleToATLAS) tmp_inputsignaldistribminus->Scale(scaleCMStoATLAS[i]);
              if(rescaleToATLAS) tmp_inputsignaldistribplus->Scale(scaleCMStoATLAS[i]);
              TString outputsignaldistribnameplus   = outputdistrib+"__"+thetaSignalList[i]+"__"+thetaStytList[j]+"__plus";
              TString outputsignaldistribnameminus  = outputdistrib+"__"+thetaSignalList[i]+"__"+thetaStytList[j]+"__minus";
              tmp_inputsignaldistribminus->SetName(outputsignaldistribnameminus );
              tmp_inputsignaldistribplus->SetName(outputsignaldistribnameplus   );
              distrib_signal_sys.push_back( (TH1F*)tmp_inputsignaldistribplus );
              distrib_signal_sys.push_back( (TH1F*)tmp_inputsignaldistribminus);
          }

          if( sampleList[i] == "TTbar_Madgraph" && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
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

              TH1F* tmp_inputdistribplus      = (TH1F*)inputfile->Get(inputdistribnameplus)->Clone() ;
              TH1F* tmp_inputdistribminus     = (TH1F*)inputfile->Get(inputdistribnameminus)->Clone() ;
              TString outputdistribnameplus   = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__plus";
              TString outputdistribnameminus  = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__minus";
              tmp_inputdistribminus->SetName(outputdistribnameminus );
              tmp_inputdistribplus->SetName(outputdistribnameplus   );
              distrib_MC_sys.push_back( (TH1F*)tmp_inputdistribplus );
              distrib_MC_sys.push_back( (TH1F*)tmp_inputdistribminus);
          }
          else
          {
              TString inputdistribnameplus    = inputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__plus";
              TString inputdistribnameminus   = inputdistrib+"__"+sampleList[i]+"__"+stytList[j]+"__minus";
              TH1F* tmp_inputdistribplus      = (TH1F*)inputfile->Get(inputdistribnameplus)->Clone() ;
              TH1F* tmp_inputdistribminus     = (TH1F*)inputfile->Get(inputdistribnameminus)->Clone() ;
              TString outputdistribnameplus   = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__plus";
              TString outputdistribnameminus  = outputdistrib+"__"+thetaSampleList[i]+"__"+thetaStytList[j]+"__minus";
              tmp_inputdistribminus->SetName(outputdistribnameminus );
              tmp_inputdistribplus->SetName(outputdistribnameplus   );
              distrib_MC_sys.push_back( (TH1F*)tmp_inputdistribplus );
              distrib_MC_sys.push_back( (TH1F*)tmp_inputdistribminus);
          }
      }
  }

  TString outputfilename;
  if      (inputdistrib == "mWT_mujets_Wregion_highpt")     outputfilename = "inputTheta_Wregion";
  else if (inputdistrib == "mWT_mujets_ttbarregion_highpt") outputfilename = "inputTheta_ttbarregion";
  else if (inputdistrib == "mWT_mujets_signalregion")       outputfilename = "inputTheta_signalregion_mWT";
  else if (inputdistrib == "MET_mujets_signalregion")       outputfilename = "inputTheta_signalregion_MET";
  else if (inputdistrib == "DeltaPhiLJ_mujets_signalregion")outputfilename = "inputTheta_signalregion_DeltaPhiLJ";
  else                                                      outputfilename = "inputTheta_"+inputdistrib;

  outputfilename+=".root";

  TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );
  outputfile->cd();

  if(doCutnCount)
  {
      distrib__DATA->Rebin(distrib__DATA->GetNbinsX());
      distrib__QCD->Rebin(distrib__QCD->GetNbinsX());
     // distrib__QCD__plus->Rebin(distrib__QCD__plus->GetNbinsX());
     // distrib__QCD__minus->Rebin(distrib__QCD__minus->GetNbinsX());
  }

  distrib__DATA->Write((outputdistrib+"__DATA").Data());
  distrib__QCD->Write();
 // distrib__QCD__plus->Write();
 // distrib__QCD__minus->Write();

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

  bool useWFlavourSplitting = true;
  bool rescaleToATLAS       = false;
  bool doCutnCount          = false;

  std::vector<TString> signalList;
  std::vector<TString> thetaSignalList;
  std::vector<double>  scaleCMStoATLAS;
  signalList.push_back("S4_100_fastSim");  thetaSignalList.push_back("S4inv100");
  scaleCMStoATLAS.push_back(1.001/11.79);
  std::vector<TString> sampleList;
  std::vector<TString> thetaSampleList;
  sampleList.push_back("TTbar_Madgraph"     );      thetaSampleList.push_back("TTbarMadgraph"   );

  if( !useWFlavourSplitting)
  {
      sampleList.push_back("WExcl"          );      thetaSampleList.push_back("WExcl"           );
  }
  else
  {
      sampleList.push_back("WExclb"         );      thetaSampleList.push_back("WExclb"          );
      sampleList.push_back("WExclc"         );      thetaSampleList.push_back("WExclc"          );
      sampleList.push_back("WExcll"         );      thetaSampleList.push_back("WExcll"          );
  }

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

  std::vector<TString> systlist;
  std::vector<TString> thetaSystlist;

  systlist.push_back("lept"             );          thetaSystlist.push_back("lept"           );
  systlist.push_back("trig"             );          thetaSystlist.push_back("trig"           );
  //systlist.push_back("PDF"            );
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
  ProdTemplate("MET_mujets_Selectedsignalregion", "METmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("DeltaPhiLJ_mujets_Selectedsignalregion", "DeltaPhiLJmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("ptW_mujets_Selectedsignalregion", "ptWmujetsSelectedSignalregion",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

  //ProdTemplate("MET_mujets_signalregion", "METmujetsSignalregion", signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("DeltaPhiLJ_mujets_signalregion", "DeltaPhiLJmujetsSignalregion", signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  //ProdTemplate("mWT_mujets_signalregion", "mWTmujetsSignalregion", signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt", signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );
  ProdTemplate("mWT_mujets_ttbarregion_highpt", "mWTmujetsttbarregionHighpt",signalList, thetaSignalList, sampleList, thetaSampleList, systlist, thetaSystlist, "../TreeReader/outputroot_withSyst/histo_merged.root", scaleCMStoATLAS, rescaleToATLAS, doCutnCount );

}
