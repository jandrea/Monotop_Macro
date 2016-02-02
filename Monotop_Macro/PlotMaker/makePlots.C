

{

 bool useElectronChannel = false;

 gROOT->ProcessLine(".L PlotStack.C+");
 gROOT->SetStyle("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptStat(0);

 //----------------
 //list of channels
 std::vector<TString> channel_list;
 if(!useElectronChannel) channel_list.push_back("mujets");
 else                    channel_list.push_back("eljets");


    /////////////////////////////////////////
    //////   QCDCorr possible values   //////
    ////   0 : QCD region (iso > 0.5)    ////
    ///    1 : Analysis wo QCD correction ///
    //     2 : Analysis with QCD correction /
    ///    3 : Analysis with systematics  ///
    /////////////////////////////////////////
    short int QCDCorr = 3;


  //---------------------------
  //define list of data samples
  std::vector<TString> dataSample_list;
  if(QCDCorr == 0)
  {
        dataSample_list.push_back("NTuple_53X_SingleMuRun2012A");
        dataSample_list.push_back("NTuple_53X_SingleMuRun2012B");
        dataSample_list.push_back("NTuple_53X_SingleMuRun2012C");
        dataSample_list.push_back("NTuple_53X_SingleMuRun2012D");
  }
  else if(QCDCorr == 1 || QCDCorr == 2 || QCDCorr == 3)
  {
        if(!useElectronChannel)
        {
            dataSample_list.push_back("SingleMuA");
            dataSample_list.push_back("SingleMuB");
            dataSample_list.push_back("SingleMuC");
            dataSample_list.push_back("SingleMuD");
        }
        else
        {
            dataSample_list.push_back("SingleElA");
            dataSample_list.push_back("SingleElB");
            dataSample_list.push_back("SingleElC");
            dataSample_list.push_back("SingleElD");
        }
  }


 //-------------------------
 //define list of MC samples
  std::vector<TString> mcSample_list;
  std::vector<int> colorVector;

  if(QCDCorr == 0)
  {
        mcSample_list.push_back("NTuple_53X_TTJetsMadgraphZ2");             colorVector.push_back(kRed+1);
        //mcSample_list.push_back("NTuple_53X_TTWJets_8TeVmadgraph");       colorVector.push_back(kRed+3);
        //mcSample_list.push_back("NTuple_53X_TTZJets_8TeVmadgraph_v2");    colorVector.push_back(kRed+3);
        //mcSample_list.push_back("NTuple_53X_TZJetsTo3LNuB");              colorVector.push_back(kMagenta);
        mcSample_list.push_back("NTuple_53X_WJetsToLNu");                   colorVector.push_back(kGreen+2);
        mcSample_list.push_back("NTuple_53X_DYJetsToLL_M-10To50");          colorVector.push_back(kAzure-2);
        mcSample_list.push_back("NTuple_53X_DYJetsToLL_M-50");              colorVector.push_back(kAzure-2);
        mcSample_list.push_back("NTuple_53X_T_s-channel");                  colorVector.push_back(kRed+2);
        mcSample_list.push_back("NTuple_53X_T_t-channel");                  colorVector.push_back(kRed+2);
        mcSample_list.push_back("NTuple_53X_T_tW-channel");                 colorVector.push_back(kRed+2);
        mcSample_list.push_back("NTuple_53X_Tbar_t-channel");               colorVector.push_back(kRed+2);
        mcSample_list.push_back("NTuple_53X_Tbar_tW-channel");              colorVector.push_back(kRed+2);
        //mcSample_list.push_back("NTuple_53X_Tbar_s-channel");             colorVector.push_back();
        mcSample_list.push_back("NTuple_53X_WZJetsIncl");                   colorVector.push_back(13);
        mcSample_list.push_back("NTuple_53X_WWJetsIncl");                   colorVector.push_back(13);
        mcSample_list.push_back("NTuple_53X_ZZJetsIncl");                   colorVector.push_back(13);
        mcSample_list.push_back("QCD_Pt-20to30_MuEnrichedPt5");             colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-30to50_MuEnrichedPt5");             colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-50to80_MuEnrichedPt5");             colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-80to120_MuEnrichedPt5");            colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-120to170_MuEnrichedPt5");           colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-170to300_MuEnrichedPt5");           colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-300to470_MuEnrichedPt5");           colorVector.push_back(kYellow+1);
  }
  else if(QCDCorr == 1 || QCDCorr == 2 || QCDCorr == 3)
  {
        //mcSample_list.push_back("TTbar_Madgraph");              colorVector.push_back(kRed+1);
        mcSample_list.push_back("TTMSDecays_central");              colorVector.push_back(kRed+1);
        //mcSample_list.push_back("TTWJets_8TeVmadgraph");      colorVector.push_back(kRed+3);
        //mcSample_list.push_back("TTZJets_8TeVmadgraph_v2");   colorVector.push_back(kRed+3);
        //mcSample_list.push_back("TZJetsTo3LNuB");             colorVector.push_back(kMagenta);
        //mcSample_list.push_back("WJets");                       colorVector.push_back(kGreen+2);
        //mcSample_list.push_back("W1Jets");                       colorVector.push_back(kGreen);
        //mcSample_list.push_back("W2Jets");                       colorVector.push_back(kGreen+2);
        //mcSample_list.push_back("W3Jets");                       colorVector.push_back(kGreen+4);
        //mcSample_list.push_back("W4Jets");                       colorVector.push_back(kGreen-2);
        //mcSample_list.push_back("WExcl");                       colorVector.push_back(kGreen-2);
        mcSample_list.push_back("WExclb");                       colorVector.push_back(kGreen+4);
        mcSample_list.push_back("WExclc");                       colorVector.push_back(kGreen);
        mcSample_list.push_back("WExcll");                       colorVector.push_back(kGreen-2);
        //mcSample_list.push_back("Wminus_Powheg");             colorVector.push_back(kGreen+2);
        //mcSample_list.push_back("Wplus_Powheg");              colorVector.push_back(kGreen+2);
        mcSample_list.push_back("DYJetsToLL_M-10To50");         colorVector.push_back(kAzure-2);
        mcSample_list.push_back("DYJetsToLL_M-50");             colorVector.push_back(kAzure-2);
        mcSample_list.push_back("T_s");                         colorVector.push_back(13);
        mcSample_list.push_back("T_t");                         colorVector.push_back(13);
        mcSample_list.push_back("T_tW");                        colorVector.push_back(13);
        mcSample_list.push_back("Tbar_t");                      colorVector.push_back(13);
        mcSample_list.push_back("Tbar_tW");                     colorVector.push_back(13);
        //mcSample_list.push_back("Tbar_s-channel");            colorVector.push_back();
        mcSample_list.push_back("WZ");                          colorVector.push_back(kRed-3);
        mcSample_list.push_back("WW");                          colorVector.push_back(kRed-3);
        mcSample_list.push_back("ZZ");                          colorVector.push_back(kRed-3);
  }

  if(QCDCorr == 1)
  {
        mcSample_list.push_back("QCD_Pt-20to30");       colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-30to50");       colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-50to80");       colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-80to120");      colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-120to170");     colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-170to300");     colorVector.push_back(kYellow+1);
        mcSample_list.push_back("QCD_Pt-300to470");     colorVector.push_back(kYellow+1);
  }

  if(QCDCorr == 2 || QCDCorr == 3)
  {
        mcSample_list.push_back("QCD");     colorVector.push_back(kYellow+1);
  }



 //-----------------------------
 //define list of signal samples
  std::vector<TString> signalSample_list;
/*
  if(QCDCorr == 0)
  {
        signalSample_list.push_back("S1_1000_100");
  }
*/  if(QCDCorr == 1 || QCDCorr == 2 || QCDCorr == 3)
  {
//        signalSample_list.push_back("S1_700_100_fullSim");
//        signalSample_list.push_back("S4_500_fastSim");
  }



 //--------------------------
 //define list of systematics
 std::vector<TString> syst_list;
/*
 syst_list.push_back("leptup");
 syst_list.push_back("leptdown");
 syst_list.push_back("trigup");
 syst_list.push_back("trigdown");
 //syst_list.push_back("PDFup");
 //syst_list.push_back("PDFdown");
 syst_list.push_back("PUup");
 syst_list.push_back("PUdown");
 syst_list.push_back("jesup");
 syst_list.push_back("jesdown");
 syst_list.push_back("jerup");
 syst_list.push_back("jerdown");
 syst_list.push_back("metunclsup");
 syst_list.push_back("metunclsdown");

 syst_list.push_back("scaleup");
 syst_list.push_back("scaledown");
 syst_list.push_back("matchingup");
 syst_list.push_back("matchingdown");
*/


 //--------------------------
 //define list of systematics
 std::vector<TString> selectionStep_list;
 //selectionStep_list.push_back("");
 //selectionStep_list.push_back("afterleptsel");
 //selectionStep_list.push_back("1bjetregion");
 //selectionStep_list.push_back("qcdLregion");
 //selectionStep_list.push_back("qcdWregion");
 //selectionStep_list.push_back("qcdTTregion");
 //selectionStep_list.push_back("qcdSregion");
 //selectionStep_list.push_back("afterbjetsel");
 //selectionStep_list.push_back("afterMWsel");
 //selectionStep_list.push_back("ttbarregion_lowjetpt");
 //selectionStep_list.push_back("QCDnormLregion");
 //selectionStep_list.push_back("QCDnormWregion");
 //selectionStep_list.push_back("QCDnormSregion");
 //selectionStep_list.push_back("QCDnormTTregion");
 //selectionStep_list.push_back("Selectedsignalregion");
 selectionStep_list.push_back("ttbarregion_2j2b");
 //selectionStep_list.push_back("ttbarregion_3j2b");
 //selectionStep_list.push_back("ttbarregion_4j2b");
 //selectionStep_list.push_back("Wregion_highpt");


 //------------------------
 //define list of variables
 std::vector<TString> variables_list;
 //variables_list.push_back("NVtx");
 //variables_list.push_back("NJet");
 //variables_list.push_back("NBJet");
 variables_list.push_back("mWT");
 //variables_list.push_back("DeltaPhiLJ");
 //variables_list.push_back("mWT_full");
 variables_list.push_back("MET");
 //variables_list.push_back("ptW");
 //variables_list.push_back("etaW");
 variables_list.push_back("LeptPt");
 //variables_list.push_back("LeptEta");
 variables_list.push_back("JetPt");
 //variables_list.push_back("JetEta");
 //variables_list.push_back("LeadJetBtagDiscr");
// variables_list.push_back("mW");
// variables_list.push_back("mWTplusMET");
// variables_list.push_back("LeptIso");
 //variables_list.push_back("DeltaPhiLJ");
 //variables_list.push_back("DeltaRLJ");
/* variables_list.push_back("JetPt1");
 variables_list.push_back("JetEta1");
 variables_list.push_back("JetPt2");
 variables_list.push_back("JetEta2");
 variables_list.push_back("JetPt3");
 variables_list.push_back("JetEta3");
 variables_list.push_back("JetPt4");
 variables_list.push_back("JetEta4");
 variables_list.push_back("DeltaPhiLJ");
 variables_list.push_back("DeltaRLJ");
*/
 bool dology = false;

  for(int iselstep=0; iselstep < selectionStep_list.size(); iselstep++)
  {
      for(int ivar=0; ivar < variables_list.size(); ivar++)
      {
          for (int ichan=0; ichan<channel_list.size(); ichan++)
          {
              PlotStack(variables_list[ivar],  channel_list[ichan],  selectionStep_list[iselstep] , dology, dataSample_list,  channel_list, mcSample_list, signalSample_list, colorVector, QCDCorr, false, useElectronChannel);
          }
      }
  } // end loop jchan

}
