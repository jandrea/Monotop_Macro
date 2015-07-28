#include "TString.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include <iostream>
#include <vector>

using namespace std;

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


void ProdTemplate(TString inputdistrib, TString outputdistrib, vector<TString> signalList, vector<double> signalXsecs, vector<TString> sampleList, vector<double> sampleXsecs, vector<TString> stytList, TString intputfilename, vector<double> scaleCMStoATLAS, bool rescaleToATLAS, bool doCutnCount, bool useElectronChannel, bool useExtrapol13TeV )
{

  double lumiscale;
  lumiscale = 16.0;
  lumiscale/=19.7;

  TString                channel = "mujets";
  if(useElectronChannel) channel = "eljets";
  TFile * inputfile	  = new TFile( intputfilename.Data() );

  TString dataAname         = inputdistrib+"__SingleMuA";
  TString dataBname         = inputdistrib+"__SingleMuB";
  TString dataCname         = inputdistrib+"__SingleMuC";
  TString dataDname         = inputdistrib+"__SingleMuD";

  if(useElectronChannel)
  {
    dataAname         = inputdistrib+"__SingleElA";
    dataBname         = inputdistrib+"__SingleElB";
    dataCname         = inputdistrib+"__SingleElC";
    dataDname         = inputdistrib+"__SingleElD";
  }

  TH1D * distrib__RUNA  	= (TH1D*)inputfile->Get(dataAname)->Clone();
  TH1D * distrib__RUNB  	= (TH1D*)inputfile->Get(dataBname)->Clone();
  TH1D * distrib__RUNC  	= (TH1D*)inputfile->Get(dataCname)->Clone();
  TH1D * distrib__RUND  	= (TH1D*)inputfile->Get(dataDname)->Clone();

  TH1D * distrib__DATA  	= (TH1D*)distrib__RUNA;
         distrib__DATA->Add(distrib__RUNB);
         distrib__DATA->Add(distrib__RUNC);
         distrib__DATA->Add(distrib__RUND);
  vector< TH1D* > distrib_MC;
  vector< TH1D* > distrib_MC_sys;
  vector< TH1D* > distrib_signal_sys;
  vector< TH1D* > distrib_signal;


  TH1D* distrib__Pseudo_DATA = new TH1D("Pseudo_DATA", "Pseudo_DATA", distrib__DATA->GetNbinsX(), distrib__DATA->GetXaxis()->GetXmin(), distrib__DATA->GetXaxis()->GetXmax());

  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(i <  signalList.size() )
      {
          TString signalname          = inputdistrib+"__"+signalList[i];
          TH1D * distrib__Monotop     = (TH1D*)inputfile->Get(signalname)->Clone();
          if(rescaleToATLAS) distrib__Monotop->Scale(scaleCMStoATLAS[i]);
          TString outputsignalname    = outputdistrib+"__"+signalList[i];
          distrib__Monotop->SetName(outputsignalname );
          if(useExtrapol13TeV) distrib__Monotop->Scale(double(lumiscale)*double(signalXsecs[i]));
          distrib_signal.push_back( (TH1D*)distrib__Monotop);
      }
      TH1D* distrib__nominale = 0;
      //if(sampleList[i] == "QCD" && inputdistrib == "mWT_"+channel+"_ttbarregion_highpt") continue;
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
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_t", inputfile));
                 //distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_tW", inputfile));
                 //distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_tW", inputfile));
      }
     else if(sampleList[i] == "SingleTopW")
      {
                 distrib__nominale          = getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_tW", inputfile);
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_tW", inputfile));
      }
      else if(sampleList[i] == "DY")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "DYJetsToLL_M-50", inputfile);
                 if(inputdistrib != "mWT_mujets_Selectedsignalregion") distrib__nominale->Add(getSystematic("" , inputdistrib, outputdistrib, "" , "DYJetsToLL_M-10To50", inputfile));
                 //if(inputdistrib != "mWT_"+channel+"_Selectedsignalregion") distrib__nominale->Add(getSystematic("" , inputdistrib, outputdistrib, "" , "DYJetsToLL_M-10To50", inputfile));
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
      if(useExtrapol13TeV) distrib__nominale->Scale(double(lumiscale)*double(sampleXsecs[i]));
      if(useExtrapol13TeV) distrib__Pseudo_DATA->Add(distrib__nominale);
      distrib_MC.push_back( (TH1D*)distrib__nominale);
  }
  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_highpt") continue;
      //if(sampleList[i] == "QCD" && inputdistrib == "mWT_"+channel+"_ttbarregion_highpt") continue;

      for(unsigned int j=0; j<stytList.size(); j++)
      {
          if( stytList[j]=="toppt" && sampleList[i] != "TTMSDecays") continue;
          //if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching") && (sampleList[i] != "WExcll" && sampleList[i] != "WExclb" && sampleList[i] != "WExclc" && sampleList[i] != "TTMSDecays") ) continue;
          if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching") && sampleList[i] != "TTMSDecays" ) continue;
          if( (stytList[j] =="Iso" || stytList[j] == "BgdContam") && sampleList[i] != "QCD" ) continue;
          if( sampleList[i] == "QCD" && (stytList[j] !="Iso" && stytList[j] != "BgdContam") ) continue;

          if (i < signalList.size() && stytList[j] != "mass" && stytList[j] != "scale" && stytList[j] != "matching" && stytList[j] != "toppt" && stytList[j] != "Iso" && stytList[j] != "BgdContam" && stytList[j] != "PDF")
          {
              TH1D* distribsignal__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              TH1D* distribsignal__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              //if(rescaleToATLAS) tmp_inputsignaldistribminus->Scale(scaleCMStoATLAS[i]);
              //if(rescaleToATLAS) tmp_inputsignaldistribplus->Scale(scaleCMStoATLAS[i]);
              distribsignal__plus->SetName( (outputdistrib+"__"+signalList[i]+"__"+stytList[j]+"__plus").Data());
              distribsignal__minus->SetName((outputdistrib+"__"+signalList[i]+"__"+stytList[j]+"__minus").Data());
              if(useExtrapol13TeV) distribsignal__plus->Scale(double(lumiscale)*double(signalXsecs[i]));
              if(useExtrapol13TeV) distribsignal__minus->Scale(double(lumiscale)*double(signalXsecs[i]));
              distrib_signal_sys.push_back( (TH1D*)distribsignal__plus );
              distrib_signal_sys.push_back( (TH1D*)distribsignal__minus);
          }

          TH1D* distribmc__plus = 0;
          TH1D* distribmc__minus = 0;

          //if( (sampleList[i] == "TTMSDecays" || sampleList[i] == "WExcll" || sampleList[i] == "WExclc" || sampleList[i] == "WExclb")  && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
          if( (sampleList[i] == "TTMSDecays" )  && (stytList[j] == "mass" || stytList[j] == "scale" || stytList[j] == "matching") )
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
              if(useExtrapol13TeV) tmp_inputdistribplus->Scale(double(lumiscale)*double(sampleXsecs[i]));
              if(useExtrapol13TeV) tmp_inputdistribminus->Scale(double(lumiscale)*double(sampleXsecs[i]));
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
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));
                    //distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile));
                    //distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                    distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_s", inputfile);
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_t", inputfile));
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));
                    //distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile));
                    //distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                }
                else if(sampleList[i] == "SingleTopW")
                {
                    distribmc__plus         = getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile);
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                    distribmc__minus        = getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_tW", inputfile);
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_tW", inputfile));

                }
                else if(sampleList[i] == "DY")
                {
                    distribmc__plus           = getSystematic("plus"    , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-50", inputfile);
                    if(inputdistrib != "mWT_"+channel+"_Selectedsignalregion") distribmc__plus->Add(   	getSystematic("plus", inputdistrib, outputdistrib, stytList[j],"DYJetsToLL_M-10To50", inputfile));

                    distribmc__minus          = getSystematic("minus"   , inputdistrib, outputdistrib, stytList[j]      , "DYJetsToLL_M-50", inputfile);
                    if(inputdistrib != "mWT_"+channel+"_Selectedsignalregion") distribmc__minus->Add(   	getSystematic("minus", inputdistrib, outputdistrib, stytList[j],"DYJetsToLL_M-10To50", inputfile));
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
                if(useExtrapol13TeV) distribmc__plus->Scale(double(lumiscale)*double(sampleXsecs[i]));
                if(useExtrapol13TeV) distribmc__minus->Scale(double(lumiscale)*double(sampleXsecs[i]));
                distrib_MC_sys.push_back( (TH1D*)distribmc__plus );
                distrib_MC_sys.push_back( (TH1D*)distribmc__minus);

          }
      }
  }

  TString outputfilename;
  if      (inputdistrib == "mWT_"+channel+"_Wregion_highpt")         outputfilename = "inputTheta_"+channel+"_merged_Wregion";
  else if (inputdistrib == "mWT_"+channel+"_ttbarregion_highpt")     outputfilename = "inputTheta_"+channel+"_merged_ttbarregion";
  else if (inputdistrib == "mWT_"+channel+"_Selectedsignalregion")   outputfilename = "inputTheta_"+channel+"_merged_Selectedsignalregion";
  else                                                               outputfilename = "inputTheta_merged_"+inputdistrib;

  outputfilename+=".root";

  TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );

  if(doCutnCount)  distrib__DATA->Rebin(distrib__DATA->GetNbinsX());

  outputfile->cd();
  if( !useExtrapol13TeV) distrib__DATA->Write((outputdistrib+"__DATA").Data());
  else
  {
      TRandom theRand;
      for(int i=1; i<= distrib__Pseudo_DATA->GetNbinsX(); i++)
      {
          distrib__Pseudo_DATA->SetBinContent(i, theRand.PoissonD( distrib__Pseudo_DATA->GetBinContent(i))  );
          distrib__Pseudo_DATA->SetBinError(i, pow(distrib__Pseudo_DATA->GetBinContent(i), 0.5) );
      }

      distrib__Pseudo_DATA->Write((outputdistrib+"__DATA").Data());
  }

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


void ProdTemplate_PseudoData()
{

  bool rescaleToATLAS       = false;
  bool doCutnCount          = false;
  bool useElectronChannel   = false;
  bool useExtrapol13TeV     = true;

  TString                channel = "mujets";
  if(useElectronChannel) channel = "eljets";
  TString                inputfilename = "../../TreeReader/outputroot_withSyst/histo_merged.root";
  if(useElectronChannel) inputfilename = "../../TreeReader/outputroot_withSyst/histo_merged_electron_syst.root";

  vector<TString> signalList;
  vector<double> signalXsecs;
  vector<double>  scaleCMStoATLAS;
  signalList.push_back("S1Res1500Inv100");  signalXsecs.push_back(0.28/0.05);
  scaleCMStoATLAS.push_back(1.001/11.79);
  vector<TString> sampleList;
  vector<double> sampleXsecs;

  sampleList.push_back("TTMSDecays"     );  sampleXsecs.push_back(831.76/245.8    );
  sampleList.push_back("WExclb"         );  sampleXsecs.push_back(61464./36864.3  );
  sampleList.push_back("WExclc"         );  sampleXsecs.push_back(61464./36864.3  );
  sampleList.push_back("WExcll"         );  sampleXsecs.push_back(61464./36864.3  );
  sampleList.push_back("DY"             );  sampleXsecs.push_back(6025./3531.     );
  sampleList.push_back("SingleTop"      );  sampleXsecs.push_back(224.17/90.89    );
  sampleList.push_back("SingleTopW"     );  sampleXsecs.push_back(71.2/22.2       );
  sampleList.push_back("VV"             );  sampleXsecs.push_back(182./92.2       );
  sampleList.push_back("QCD"            );  sampleXsecs.push_back(831.76/245.8    );



  vector<TString> systlist;

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


  //ProdTemplate("mWT_"+channel+"_ATLASRESsignalregion", "mWT"+channel+"ATLASRESSignalregion",signalList, signalXsecs, sampleList,  sampleXsecs, systlist, inputfilename, scaleCMStoATLAS, rescaleToATLAS, doCutnCount, useElectronChannel, useExtrapol13TeV );
  //ProdTemplate("mWT_"+channel+"_ATLASFCNCsignalregion", "mWT"+channel+"ATLASFCNCSignalregion",signalList, signalXsecs, sampleList,  sampleXsecs, systlist, inputfilename, scaleCMStoATLAS, rescaleToATLAS, doCutnCount, useElectronChannel, useExtrapol13TeV );
  ProdTemplate("mWT_"+channel+"_Selectedsignalregion", "mWT"+channel+"SelectedSignalregion",signalList, signalXsecs, sampleList,  sampleXsecs, systlist,  inputfilename, scaleCMStoATLAS, rescaleToATLAS, doCutnCount, useElectronChannel, useExtrapol13TeV );
  ProdTemplate("mWT_"+channel+"_Wregion_highpt", "mWT"+channel+"WregionHighpt", signalList, signalXsecs, sampleList,  sampleXsecs, systlist,  inputfilename, scaleCMStoATLAS, rescaleToATLAS, doCutnCount, useElectronChannel, useExtrapol13TeV  );
  ProdTemplate("mWT_"+channel+"_ttbarregion_highpt", "mWT"+channel+"ttbarregionHighpt",signalList, signalXsecs, sampleList,  sampleXsecs, systlist,  inputfilename, scaleCMStoATLAS, rescaleToATLAS, doCutnCount, useElectronChannel, useExtrapol13TeV  );

}
