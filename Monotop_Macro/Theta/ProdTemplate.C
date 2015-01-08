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



void setNjetSel(TH1F * thehisto, int njetsel){

  if(njetsel!=0) {
    for(int i=1; i<=njetsel; i++){
      thehisto->SetBinContent(i, 0.0000000001);
      thehisto->SetBinError(i, 0.0000000001);
    }
  }

}


void ProdTemplate(TString inputdistrib, TString outputdistrib, std::vector<TString> sampleList, std::vector<TString> thetaSampleList, std::vector<TString> stytList, std::vector<TString> thetaStytList, TString intputfilename){

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

  std::vector< TH1F* > distrib_MC;
  std::vector< TH1F* > distrib_MC_sys;
  std::vector< TH1F* > distrib_signal_sys;




  for(unsigned int i=0; i<sampleList.size(); i++){
    //cout << inputdistrib+"__"+sampleList[i]  << endl;
    if(i == 0)
    {
        TString signalname          = inputdistrib+"__S1";
        TH1F * distrib__Monotop     = (TH1F*)inputfile->Get(signalname)->Clone();
        TString outputsignalname    = outputdistrib+"__Monotop";
        distrib__Monotop->SetName(outputsignalname );
    }
    TString inputdistribname    = inputdistrib+"__"+sampleList[i];
    TH1F* tmp_inputdistrib      = (TH1F*)inputfile->Get(inputdistribname)->Clone() ;
    TString outputdistribname   = outputdistrib+"__"+thetaSampleList[i];
    tmp_inputdistrib->SetName(outputdistribname );
    distrib_MC.push_back( (TH1F*)tmp_inputdistrib);
  }


  for(unsigned int i=0; i<sampleList.size(); i++){
    for(unsigned int j=0; j<stytList.size(); j++){
      if (i == 0)
      {
        TString signalnameplus                = inputdistrib+"__S1__"+stytList[j]+"__plus";
        TString signalnameminus               = inputdistrib+"__S1__"+stytList[j]+"__minus";
        TH1F* tmp_inputsignaldistribminus     = (TH1F*)inputfile->Get(signalnameplus)->Clone() ;
        TH1F* tmp_inputsignaldistribplus      = (TH1F*)inputfile->Get(signalnameminus)->Clone() ;
        TString outputsignaldistribnameplus   = outputdistrib+"__Monotop__"+thetaStytList[j]+"__plus";
        TString outputsignaldistribnameminus  = outputdistrib+"__Monotop__"+thetaStytList[j]+"__minus";
        tmp_inputsignaldistribminus->SetName(outputsignaldistribnameminus );
        tmp_inputsignaldistribplus->SetName(outputsignaldistribnameplus   );
        distrib_signal_sys.push_back( (TH1F*)tmp_inputsignaldistribplus );
        distrib_signal_sys.push_back( (TH1F*)tmp_inputsignaldistribminus);
      }
      if(stytList[j] == "scale"    && sampleList[i] !="TT") continue;
      if(stytList[j] == "matching" && sampleList[i] !="TT") continue;
      //cout << "  " << (inputdistrib+"__"+sampleList[i]+"__"+stytList[j]) << endl;


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

  //for(unsigned int i=0; i<stytList.size(); i++) outputfilename+="_"+stytList[i];
  //outputfilename+=".root";

  //TFile * outputfile = new TFile(  outputfilename.Data(), "recreate" );
  TFile * outputfile = new TFile(  "inputTheta.root", "recreate" );
  outputfile->cd();

  distrib__DATA->Write((outputdistrib+"__DATA").Data());

  for(unsigned int i=0; i<distrib_MC.size(); i++)
  {
      if (i == 0) distrib__Monotop->Write();
      distrib_MC[i]->Write();
  }
  for(unsigned int i=0; i<distrib_signal_sys.size(); i++)
  {
      distrib_signal_sys[i]->Write();
  }
  for(unsigned int i=0; i<distrib_MC_sys.size(); i++)
  {
      distrib_MC_sys[i]->Write();
  }

}


void ProdTemplate(){

  std::vector<TString> sampleList;
  std::vector<TString> thetaSampleList;
  sampleList.push_back("TTbar_Madgraph"     );      thetaSampleList.push_back("TTbarMadgraph"   );
  sampleList.push_back("WJets"              );      thetaSampleList.push_back("WJets"           );
  sampleList.push_back("DYJetsToLL_M-10To50");      thetaSampleList.push_back("DY10To50"        );
  sampleList.push_back("DYJetsToLL_M-50"    );      thetaSampleList.push_back("DY50"            );
  sampleList.push_back("T_s"                );      thetaSampleList.push_back("Ts"              );
  sampleList.push_back("T_t"                );      thetaSampleList.push_back("Tt"              );
  sampleList.push_back("T_tW"               );      thetaSampleList.push_back("TtW"             );
  sampleList.push_back("Tbar_t"             );      thetaSampleList.push_back("Tbart"           );
  sampleList.push_back("Tbar_tW"            );      thetaSampleList.push_back("TbartW"          );
  //sampleList.push_back("TTZ"              );
  //sampleList.push_back("TTW"              );
  sampleList.push_back("WZ"                 );      thetaSampleList.push_back("WZ"              );
  sampleList.push_back("ZZ"                 );      thetaSampleList.push_back("ZZ"              );
  sampleList.push_back("WW"                 );      thetaSampleList.push_back("WW"              );
  //sampleList.push_back("tZq"              );
  sampleList.push_back("QCD_A"              );      thetaSampleList.push_back("QCDA"            );
  sampleList.push_back("QCD_B"              );      thetaSampleList.push_back("QCDB"            );
  sampleList.push_back("QCD_C"              );      thetaSampleList.push_back("QCDC"            );
  sampleList.push_back("QCD_D"              );      thetaSampleList.push_back("QCDD"            );


  std::vector<TString> systlist;
  std::vector<TString> thetaSystlist;

  systlist.push_back("lept"             );          thetaSystlist.push_back("lept"           );
  //systlist.push_back("trig"           );
  //systlist.push_back("PDF"            );
  //systlist.push_back("PU"             );
  systlist.push_back("jes"              );          thetaSystlist.push_back("jes"            );
  systlist.push_back("jer"              );          thetaSystlist.push_back("jer"            );
  //systlist.push_back("metuncls"       );

  //systlist.push_back("scale"          );
  //systlist.push_back("matching"       );
  //systlist.push_back("toppt"          );


  systlist.push_back("btag__JES"        );          thetaSystlist.push_back("btagJES"        );
  systlist.push_back("btag__CSVLF"      );          thetaSystlist.push_back("btagCSVLF"      );
  systlist.push_back("btag__CSVHFStats1");          thetaSystlist.push_back("btagCSVHFStats1");
  systlist.push_back("btag__CSVHFStats2");          thetaSystlist.push_back("btagCSVHFStats2");
  systlist.push_back("btag__CSVCErr1"   );          thetaSystlist.push_back("btagCSVCErr1"   );
  systlist.push_back("btag__CSVCErr2"   );          thetaSystlist.push_back("btagCSVCErr2"   );
  systlist.push_back("btag__CSVHF"      );          thetaSystlist.push_back("btagCSVHF"      );
  systlist.push_back("btag__CSVLFStats1");          thetaSystlist.push_back("btagCSVLFStats1");
  systlist.push_back("btag__CSVLFStats2");          thetaSystlist.push_back("btagCSVLFStats2");


  ProdTemplate("mWT_mujets_signalregion", "mWTmujetsSignalregion",        sampleList,  thetaSampleList, systlist, thetaSystlist,  "../TreeReader/outputroot_withSyst/histo_merged.root" );
//  ProdTemplate("mWT_mujets_Wregion_highpt", "mWTmujetsWregionHighpt",     sampleList,  thetaSampleList, systlist, thetaSystlist,  "../TreeReader/outputroot_withSyst/histo_merged.root" );
//  ProdTemplate("mWT_mujets_1bjetregion", "mWTmujets1bjetregion",     sampleList,  thetaSampleList, systlist, thetaSystlist,  "../TreeReader/outputroot_withSyst/histo_merged.root" );
//  ProdTemplate("mWT_mujets_ttbarregion_highpt", sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", "mWT_mujets_ttbarregion_highpt"  );
//  ProdTemplate("mWT_mujets_ttbarregion_lowpt",  sampleList,  systlist,  "../TreeReader/outputroot_withSyst/histo_merged.root", "mWT_mujets_ttbarregion_lowpt"   );

}
