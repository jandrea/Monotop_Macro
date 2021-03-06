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

using namespace std;

TH1D* getSystematic(TString plusORminus, TString inputdistrib, TString outputdistrib, TString syst, TString sample, TFile* inputfile)
{
    TString fullSystName;
    TString nominalPDFName = inputdistrib+"__"+sample;

    if (plusORminus != "" && syst != "")
    {
        fullSystName   = inputdistrib+"__"+sample+"__"+syst+"__"+plusORminus;
    }
    else
    {
        fullSystName   = inputdistrib+"__"+sample;
    }
    TH1D * distrib__sample  	= (TH1D*)inputfile->Get(fullSystName)->Clone();
    //TH1D * distrib__nominal  	= (TH1D*)inputfile->Get(fullSystName)->Clone();

    return distrib__sample;
}

TH1D* mergeFiveLastBins(TH1D* initialHisto)
{
    int nInitialBins = initialHisto->GetNbinsX();
    int nNewBins     = initialHisto->GetNbinsX() - 5;
    Float_t binsLowerBand[nNewBins + 1];
    binsLowerBand[0] = initialHisto->GetXaxis()->GetXmin();
    int binWidth = (initialHisto->GetXaxis()->GetXmax() - initialHisto->GetXaxis()->GetXmin() )/nInitialBins;

    for( unsigned short int j = 1; j<= nNewBins; j++)
    {
        binsLowerBand[j] = binsLowerBand[j-1] + binWidth;
    }
    //TH1D *newHisto = new TH1D(initialHisto->GetName(), initialHisto->GetName(), nNewBins, initialHisto->GetXaxis()->GetMinimum(), initialHisto->GetXaxis()->GetMaximum() );
    //TH1D *newHisto = new TH1D(initialHisto->GetName(), initialHisto->GetName(), nNewBins, initialHisto->GetXaxis()->GetXmin(), initialHisto->GetXaxis()->GetXmax() );
    TH1D *newHisto = new TH1D(initialHisto->GetName(), initialHisto->GetName(), nNewBins, binsLowerBand );
//cout << "Histo: " << initialHisto->GetName() << " | Xmin = " << initialHisto->GetXaxis()->GetXmin()<< " | Xmax= " << initialHisto->GetXaxis()->GetXmax()  << endl;
    newHisto->Sumw2();
    double errorTmp = 0;
    double contentTmp = 0;

    for ( unsigned short int bin = 1; bin <= nInitialBins; bin++)
    {
        if(bin < nNewBins)
        {
            newHisto->SetBinContent(bin, initialHisto->GetBinContent(bin) );
            newHisto->SetBinError(bin, initialHisto->GetBinError(bin) );
        }
        else
        {
            contentTmp += initialHisto->GetBinContent(bin);
            errorTmp += pow(initialHisto->GetBinError(bin), 2);
        }
    }

    newHisto->SetBinContent(newHisto->GetNbinsX(), contentTmp );
    newHisto->SetBinError(newHisto->GetNbinsX(), sqrt(errorTmp) );

    return newHisto;

}

void ProdTemplate(TString inputdistrib, TString outputdistrib, vector<TString> signalList, vector<TString> thetaSignalList, vector<TString> sampleList, vector<TString> stytList, TString intputfilename, vector<double> scaleCMStoATLAS, bool doCutnCount, bool useElectronChannel, bool mergeLastBins , bool mergeWsamples )
{

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


  for(unsigned int i=0; i<sampleList.size(); i++)
  {
      if(i <  signalList.size() )
      {
          TString signalname          = inputdistrib+"__"+signalList[i];
          TH1D * distrib__Monotop     = (TH1D*)inputfile->Get(signalname)->Clone();
          TString outputsignalname    = outputdistrib+"__"+thetaSignalList[i];
          distrib__Monotop->SetName(outputsignalname );
          distrib_signal.push_back( (TH1D*)distrib__Monotop);
      }
      TH1D* distrib__nominale = 0;
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_2j2b") continue;
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_3j2b") continue;
      if(sampleList[i] == "QCD" && inputdistrib == "mWT_mujets_ttbarregion_4j2b") continue;

      if(sampleList[i] == "TTMSDecays")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "TTMSDecays_central", inputfile);
      }
      else if(sampleList[i] == "WExcll" ||sampleList[i] == "WExclc" || sampleList[i] == "WExclb" || sampleList[i] == "QCD")
      {
                 if (!mergeWsamples || sampleList[i] == "QCD") distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , sampleList[i], inputfile);
                 else
                 {
                     distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "WExclb", inputfile);
                     distrib__nominale->Add(      getSystematic(""     , inputdistrib, outputdistrib, ""       , "WExclc", inputfile));
                     distrib__nominale->Add(      getSystematic(""     , inputdistrib, outputdistrib, ""       , "WExcll", inputfile));
                 }
      }
      else if(sampleList[i] == "SingleTop")
      {
                 distrib__nominale         	= getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_s", inputfile);
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "T_t", inputfile));
                 distrib__nominale->Add(   	  getSystematic(""     , inputdistrib, outputdistrib, ""       , "Tbar_t", inputfile));
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


      // ---------------------------------------------------------------
      // --------------- deal with systematics -------------------------
      // ---------------------------------------------------------------
      for(unsigned int j=0; j<stytList.size(); j++)
      {
          if( stytList[j]=="toppt" && sampleList[i] != "TTMSDecays") continue;
          if( (stytList[j]=="mass" || stytList[j]=="scale" || stytList[j]=="matching") && sampleList[i] != "TTMSDecays" ) continue;
          if( (stytList[j] =="Iso" || stytList[j] == "BgdContam") && sampleList[i] != "QCD" ) continue;
          if( sampleList[i] == "QCD" && (stytList[j] !="Iso" && stytList[j] != "BgdContam") ) continue;

          if (i < signalList.size() && stytList[j] != "mass" && stytList[j] != "scale" && stytList[j] != "matching" && stytList[j] != "toppt" && stytList[j] != "Iso" && stytList[j] != "BgdContam")
          {
              TH1D* distribsignal__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              TH1D* distribsignal__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , signalList[i], inputfile);
              distribsignal__plus->SetName( (outputdistrib+"__"+thetaSignalList[i]+"__"+stytList[j]+"__plus").Data());
              distribsignal__minus->SetName((outputdistrib+"__"+thetaSignalList[i]+"__"+stytList[j]+"__minus").Data());
              distrib_signal_sys.push_back( (TH1D*)distribsignal__plus );
              distrib_signal_sys.push_back( (TH1D*)distribsignal__minus);
          }

          TH1D* distribmc__plus = 0;
          TH1D* distribmc__minus = 0;

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
                    if(!mergeWsamples || sampleList[i] == "QCD") distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
                    else
                    {
                        distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "WExclb", inputfile);
                        distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "WExclc", inputfile));
                        distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]       , "WExcll", inputfile));
                    }
                    if(!mergeWsamples || sampleList[i] == "QCD") distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , sampleList[i], inputfile);
                    else
                    {
                        distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "WExclb", inputfile);
                        distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "WExclc", inputfile));
                        distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]       , "WExcll", inputfile));
                    }
                }
                else if(sampleList[i] == "SingleTop")
                {
                    distribmc__plus     	= getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_s", inputfile);
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "T_t", inputfile));
                    distribmc__plus->Add(     getSystematic("plus"     , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));

                    distribmc__minus     	= getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_s", inputfile);
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "T_t", inputfile));
                    distribmc__minus->Add(    getSystematic("minus"    , inputdistrib, outputdistrib, stytList[j]      , "Tbar_t", inputfile));

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
                distrib_MC_sys.push_back( (TH1D*)distribmc__plus );
                distrib_MC_sys.push_back( (TH1D*)distribmc__minus);

          }
      }
  }

  TString outputfilename;
  if      (inputdistrib == "mWT_"+channel+"_Wregion_highpt")         outputfilename = "inputTheta_"+channel+"_merged_Wregion";
  else if (inputdistrib == "mWT_"+channel+"_ttbarregion_2j2b")       outputfilename = "inputTheta_"+channel+"_merged_ttbarregion_2j2b";
  else if (inputdistrib == "mWT_"+channel+"_ttbarregion_3j2b")       outputfilename = "inputTheta_"+channel+"_merged_ttbarregion_3j2b";
  else if (inputdistrib == "mWT_"+channel+"_ttbarregion_4j2b")       outputfilename = "inputTheta_"+channel+"_merged_ttbarregion_4j2b";
  else if (inputdistrib == "mWT_"+channel+"_Selectedsignalregion")   outputfilename = "inputTheta_"+channel+"_merged_Selectedsignalregion";
  else if (inputdistrib == "mWT_"+channel+"_interRegion")            outputfilename = "inputTheta_"+channel+"_merged_interRegion";
  else                                                               outputfilename = "inputTheta_merged_"+inputdistrib;

  outputfilename+=".root";

  TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );

  if(mergeLastBins) distrib__DATA = mergeFiveLastBins( distrib__DATA);
  if(doCutnCount)  distrib__DATA->Rebin(distrib__DATA->GetNbinsX());

  outputfile->cd();
  distrib__DATA->Write((outputdistrib+"__DATA").Data());

  for(unsigned int i=0; i<distrib_MC.size(); i++)
  {
      if(mergeLastBins) distrib_MC[i] = mergeFiveLastBins( distrib_MC[i]);
      if(doCutnCount) distrib_MC[i]->Rebin(distrib_MC[i]->GetNbinsX());
      distrib_MC[i]->Write();
  }
  for(unsigned int i=0; i<distrib_signal.size(); i++)
  {
      if(mergeLastBins) distrib_signal[i] = mergeFiveLastBins( distrib_signal[i]);
      if(doCutnCount) distrib_signal[i]->Rebin(distrib_signal[i]->GetNbinsX());
  //    distrib_signal[i]->Scale(0.0000000001/distrib_signal[i]->Integral());
      distrib_signal[i]->Write();
  }
  for(unsigned int i=0; i<distrib_signal_sys.size(); i++)
  {
      if(mergeLastBins) distrib_signal_sys[i] = mergeFiveLastBins( distrib_signal_sys[i]);
      if(doCutnCount) distrib_signal_sys[i]->Rebin(distrib_signal_sys[i]->GetNbinsX());
  //    distrib_signal_sys[i]->Scale(0.0000000001/distrib_signal_sys[i]->Integral());
      distrib_signal_sys[i]->Write();
  }
  for(unsigned int i=0; i<distrib_MC_sys.size(); i++)
  {
      if(mergeLastBins) distrib_MC_sys[i] = mergeFiveLastBins( distrib_MC_sys[i]);
      if(doCutnCount) distrib_MC_sys[i]->Rebin(distrib_MC_sys[i]->GetNbinsX());
      distrib_MC_sys[i]->Write();
  }
}


void ProdTemplate_mergedSamples(){

  bool mergeLastBins        = false;
  bool doCutnCount          = false;
  bool useElectronChannel   = false;
  bool mergeWsamples        = false;

  TString                channel = "mujets";
  if(useElectronChannel) channel = "eljets";
  TString                inputfilename = "../TreeReader/outputroot_withSyst/histo_merged.root";
  if(useElectronChannel) inputfilename = "../TreeReader/outputroot_withSyst/histo_merged_electron_syst.root";

  vector<TString> signalList;
  vector<TString> thetaSignalList;
  vector<double>  scaleCMStoATLAS;
  signalList.push_back("S4Inv1000");  thetaSignalList.push_back("S4Inv1000");
  //signalList.push_back("S1Res2100Inv200");  thetaSignalList.push_back("S1Res2100Inv200");
  //signalList.push_back("S4Inv1000");        thetaSignalList.push_back("S4Inv1000");
  //signalList.push_back("S1Res2100Inv50");  thetaSignalList.push_back("S1Res2100Inv50");
  //signalList.push_back("S1Res1300Inv100");  thetaSignalList.push_back("S1Res1300Inv100");
  //signalList.push_back("S1Res2100Inv100");  thetaSignalList.push_back("S1Res2100Inv100");
  scaleCMStoATLAS.push_back(1.001/11.79);
  vector<TString> sampleList;

  sampleList.push_back("TTMSDecays"     );
  sampleList.push_back("WExclb"         );
  if(!mergeWsamples) sampleList.push_back("WExclc"         );
  if(!mergeWsamples) sampleList.push_back("WExcll"         );
  sampleList.push_back("DY"             );
  sampleList.push_back("SingleTop"      );
  sampleList.push_back("SingleTopW"     );
  sampleList.push_back("VV"             );
  sampleList.push_back("QCD"            );



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


//  ProdTemplate("mWT_"+channel+"_ATLASRESsignalregion", "mWT"+channel+"ATLASRESSignalregion",signalList, thetaSignalList, sampleList,  systlist, inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples );
  //ProdTemplate("mWT_"+channel+"_ATLASFCNCsignalregion", "mWT"+channel+"ATLASFCNCSignalregion",signalList, thetaSignalList, sampleList,  systlist, inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples  );
  ProdTemplate("mWT_"+channel+"_Selectedsignalregion", "mWT"+channel+"SelectedSignalregion",signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples  );
  ProdTemplate("mWT_"+channel+"_interRegion", "mWT"+channel+"interRegion",signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples  );
  ProdTemplate("mWT_"+channel+"_Wregion_highpt", "mWT"+channel+"WregionHighpt", signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples    );
  ProdTemplate("mWT_"+channel+"_ttbarregion_2j2b", "mWT"+channel+"ttbarregionHighpt",signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples   );
  //ProdTemplate("mWT_"+channel+"_ttbarregion_3j2b", "mWT"+channel+"ttbarregionHighpt",signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples   );
  //ProdTemplate("mWT_"+channel+"_ttbarregion_4j2b", "mWT"+channel+"ttbarregionHighpt",signalList, thetaSignalList, sampleList,  systlist,  inputfilename, scaleCMStoATLAS, doCutnCount, useElectronChannel, mergeLastBins, mergeWsamples  );

}
