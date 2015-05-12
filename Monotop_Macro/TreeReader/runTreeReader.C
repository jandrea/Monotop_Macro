{
    gSystem->Load("../BTagSFNotReshaping/BTagSFUtil_C.so");
    gROOT->ProcessLine(".L TreeReader.C+");

    vector<TString > emptysystlist;
    emptysystlist.push_back(""                      );

    vector<TString > systlist;
    systlist.push_back(""                           );
    systlist.push_back("lept__plus"                 );
    systlist.push_back("lept__minus"                );
    systlist.push_back("trig__plus"                 );
    systlist.push_back("trig__minus"                );
    //systlist.push_back("PDF__plus"                  );
    //systlist.push_back("PDF__minus"                 );
    systlist.push_back("PU__plus"                   );
    systlist.push_back("PU__minus"                  );
    systlist.push_back("jes__plus"                  );
    systlist.push_back("jes__minus"                 );
    systlist.push_back("jer__plus"                  );
    systlist.push_back("jer__minus"                 );
    systlist.push_back("metuncls__plus"             );
    systlist.push_back("metuncls__minus"            );
    systlist.push_back("toppt__plus"                );
    systlist.push_back("toppt__minus"               );
    systlist.push_back("btag__plus"                 );
    systlist.push_back("btag__minus"                );
    systlist.push_back("mistag__plus"               );
    systlist.push_back("mistag__minus"              );


    /////////////////////////////////////////
    //////   CorrOption possible values /////
    ////   0 : QCD region (iso > 0.4)    ////
    ///    1 : Analysis wo QCD correction ///
    //     2 : Analysis with QCD correction /
    ///    3 : Analysis with systematics  ///
    /////////////////////////////////////////
    short int CorrOption = 3;


    vector<TString > datalist;
    if(CorrOption == 0)
    {
        datalist.push_back("NTuple_53X_SingleMuRun2012A");
        datalist.push_back("NTuple_53X_SingleMuRun2012B");
        datalist.push_back("NTuple_53X_SingleMuRun2012C");
        datalist.push_back("NTuple_53X_SingleMuRun2012D");
    }
    else
    {

        datalist.push_back("SingleMuA");
        datalist.push_back("SingleMuB");
        datalist.push_back("SingleMuC");
        datalist.push_back("SingleMuD");
/*
        datalist.push_back("SingleElA");
        datalist.push_back("SingleElB");
        datalist.push_back("SingleElC");
        datalist.push_back("SingleElD");
*/
    }

    vector<TString > mclist;
    if(CorrOption == 0)
    {
        mclist.push_back("NTuple_53X_TTJetsMadgraphZ2"            );
        mclist.push_back("NTuple_53X_WJetsToLNu"                  );
        mclist.push_back("NTuple_53X_DYJetsToLL_M-10To50"         );
        mclist.push_back("NTuple_53X_DYJetsToLL_M-50"             );
        mclist.push_back("NTuple_53X_T_s-channel"                 );
        mclist.push_back("NTuple_53X_T_t-channel"                 );
        mclist.push_back("NTuple_53X_T_tW-channel"                );
        mclist.push_back("NTuple_53X_Tbar_t-channel"              );
        mclist.push_back("NTuple_53X_Tbar_tW-channel"             );
        mclist.push_back("NTuple_53X_WZJetsIncl"                  );
        mclist.push_back("NTuple_53X_WWJetsIncl"                  );
        mclist.push_back("NTuple_53X_ZZJetsIncl"                  );
        mclist.push_back("QCD_Pt-20to30_MuEnrichedPt5"            );
        mclist.push_back("QCD_Pt-30to50_MuEnrichedPt5"            );
        mclist.push_back("QCD_Pt-50to80_MuEnrichedPt5"            );
        mclist.push_back("QCD_Pt-80to120_MuEnrichedPt5"           );
        mclist.push_back("QCD_Pt-120to170_MuEnrichedPt5"          );
        mclist.push_back("QCD_Pt-170to300_MuEnrichedPt5"          );
        mclist.push_back("QCD_Pt-300to470_MuEnrichedPt5"          );
    }
    else
    {
        mclist.push_back("TTMSDecays_central"   );
        mclist.push_back("W0Jets"               );
        mclist.push_back("W1Jets"               );
        mclist.push_back("W2Jets"               );
        mclist.push_back("W3Jets"               );
        mclist.push_back("W4Jets"               );
        mclist.push_back("DYJetsToLL_M-10To50"  );
        mclist.push_back("DYJetsToLL_M-50"      );
        mclist.push_back("T_s"                  );
        mclist.push_back("T_t"                  );
        mclist.push_back("T_tW"                 );
        mclist.push_back("Tbar_t"               );
        mclist.push_back("Tbar_tW"              );
        mclist.push_back("WZ"                   );
        mclist.push_back("WW"                   );
        mclist.push_back("ZZ"                   );
    }

    vector<TString > qcdlist;
    qcdlist.push_back("QCD_Pt-20to30"   );
    qcdlist.push_back("QCD_Pt-30to50"   );
    qcdlist.push_back("QCD_Pt-50to80"   );
    qcdlist.push_back("QCD_Pt-80to120"  );
    qcdlist.push_back("QCD_Pt-120to170" );
    qcdlist.push_back("QCD_Pt-170to300" );
    qcdlist.push_back("QCD_Pt-300to470" );

    vector<TString > qcdcorrectedlist;
    qcdcorrectedlist.push_back("QCD");

    vector<TString > signallist;

    signallist.push_back("S1Res1500Inv10");
    signallist.push_back("S1Res1300Inv10");
    signallist.push_back("S1Res1100Inv10");
    signallist.push_back("S1Res900Inv10");
    signallist.push_back("S1Res700Inv10");
    signallist.push_back("S1Res500Inv10");
    //signallist.push_back("S1Res300Inv10");

    signallist.push_back("S1Res1500Inv50");
    signallist.push_back("S1Res1300Inv50");
    signallist.push_back("S1Res1100Inv50");
    signallist.push_back("S1Res900Inv50");
    signallist.push_back("S1Res700Inv50");
    signallist.push_back("S1Res500Inv50");
    //signallist.push_back("S1Res300Inv50");

    signallist.push_back("S1Res1500Inv100");
    signallist.push_back("S1Res1300Inv100");
    signallist.push_back("S1Res1100Inv100");
    signallist.push_back("S1Res900Inv100");
    signallist.push_back("S1Res700Inv100");
    signallist.push_back("S1Res500Inv100");
    signallist.push_back("S1Res300Inv100");

    signallist.push_back("S1Res1500Inv150");
    signallist.push_back("S1Res1300Inv150");
    signallist.push_back("S1Res1100Inv150");
    signallist.push_back("S1Res900Inv150");
    signallist.push_back("S1Res700Inv150");
    signallist.push_back("S1Res500Inv150");

    signallist.push_back("S1Res1500Inv200");
    signallist.push_back("S1Res1300Inv200");
    signallist.push_back("S1Res1100Inv200");
    signallist.push_back("S1Res900Inv200");
    signallist.push_back("S1Res700Inv200");
    signallist.push_back("S1Res500Inv200");


    signallist.push_back("S4Inv100");
    signallist.push_back("S4Inv200");
    signallist.push_back("S4Inv300");
    signallist.push_back("S4Inv400");
    signallist.push_back("S4Inv500");
    signallist.push_back("S4Inv600");
    signallist.push_back("S4Inv700");
    signallist.push_back("S4Inv800");
    signallist.push_back("S4Inv900");
    signallist.push_back("S4Inv1000");

    vector<TString > TTWsystlist;
    //TTsystlist.push_back("TTMSDecays_mass171_5");
    //TTsystlist.push_back("TTMSDecays_mass173_5");
    TTWsystlist.push_back("TTMSDecays_matchingdown");
    TTWsystlist.push_back("TTMSDecays_matchingup");
    TTWsystlist.push_back("TTMSDecays_scaledown");
    TTWsystlist.push_back("TTMSDecays_scaleup");
/*
    TTWsystlist.push_back("WJets_matchingdown");
    TTWsystlist.push_back("WJets_matchingup");
    TTWsystlist.push_back("WJets_scaledown");
    TTWsystlist.push_back("WJets_scaleup");
*/

  TTree* tree=0;
  //if(CorrOption == 3 || CorrOption == 2) systlist = emptysystlist;

  if(CorrOption == 0)
  {
    TreeReader * tree_ = new TreeReader(CorrOption, tree, "QCDdatadriven", 0);
    tree_.Loop(CorrOption, datalist, mclist, "QCDdatadriven", systlist, "noflav", 0);
    delete tree_;

    for (unsigned int imc = 0; imc < mclist.size(); imc++)
    {
        TreeReader * tree_ = new TreeReader(CorrOption, tree, mclist[imc], 0);
        tree_.Loop(CorrOption, datalist, mclist, mclist[imc], systlist, "noflav", 0);
        delete tree_;
    }

    for (unsigned int idata = 0; idata < datalist.size(); idata++)
    {
        TreeReader * tree_ = new TreeReader(CorrOption, tree, datalist[idata], 0);
        tree_.Loop(CorrOption, datalist, mclist, datalist[idata], systlist, "noflav", 0);
        delete tree_;
    }

  }

  if(CorrOption == 1 || CorrOption == 2 || CorrOption == 3)
  {

    for (unsigned int isig = 0; isig < signallist.size(); isig++)
    {
        TreeReader * tree_ = new TreeReader(CorrOption, tree, signallist[isig], 0);
        tree_.Loop(CorrOption, datalist,  mclist, signallist[isig], systlist, "noflav", 0);
        delete tree_;
    }

/*    for (unsigned int imc = 0; imc < mclist.size(); imc++)
    {
        if (mclist[imc] != "WJets" && mclist[imc] != "W0Jets" && mclist[imc] != "W1Jets" && mclist[imc] != "W2Jets" && mclist[imc] != "W3Jets" && mclist[imc] != "W4Jets")
        {
            TreeReader * tree_ = new TreeReader(CorrOption, tree, mclist[imc] , 0);
            tree_.Loop(CorrOption, datalist,  mclist, mclist[imc], systlist, "noflav", 0);
            delete tree_;
        }
        else
        {

            TreeReader * tree_ = new TreeReader(CorrOption, tree, mclist[imc], 0);
            tree_.Loop(CorrOption, datalist,  mclist, mclist[imc], systlist, "b", 0);
            delete tree_;

            TreeReader * tree_ = new TreeReader(CorrOption, tree, mclist[imc], 0);
            tree_.Loop(CorrOption, datalist,  mclist, mclist[imc], systlist, "c", 0);
            delete tree_;

            TreeReader * tree_ = new TreeReader(CorrOption, tree, mclist[imc], 0);
            tree_.Loop(CorrOption, datalist,  mclist, mclist[imc], systlist, "l", 0);
            delete tree_;

        }
    }
*//*
    for (unsigned int idata = 0; idata < datalist.size(); idata++)
    {
        TreeReader * tree_ = new TreeReader(CorrOption, tree, datalist[idata], 0);
        tree_.Loop(CorrOption, datalist,  mclist, datalist[idata], emptysystlist, "noflav", 0);
        delete tree_;
    }
*/
  }

  if(CorrOption == 1)
  {
    for (unsigned int iqcd = 0; iqcd < qcdlist.size(); iqcd++)
    {
        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdlist[iqcd], 0);
        tree_.Loop(CorrOption, datalist,  mclist, qcdlist[iqcd], systlist, "noflav", 0);
        delete tree_;
    }
  }
/*
  if(CorrOption == 2 || CorrOption == 3)
  {

    for (unsigned int iqcd = 0; iqcd < qcdcorrectedlist.size(); iqcd++)
    {

//        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdcorrectedlist[iqcd], 0);
//        tree_.Loop(CorrOption, datalist,  mclist, qcdcorrectedlist[iqcd], emptysystlist, "noflav", 0);
//        delete tree_;

        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdcorrectedlist[iqcd], -2);
        tree_.Loop(CorrOption, datalist,  mclist, qcdcorrectedlist[iqcd], emptysystlist, "noflav", -2);
        delete tree_;

        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdcorrectedlist[iqcd], -1);
        tree_.Loop(CorrOption, datalist,  mclist, qcdcorrectedlist[iqcd], emptysystlist, "noflav", -1);
        delete tree_;

        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdcorrectedlist[iqcd], 1);
        tree_.Loop(CorrOption, datalist,  mclist, qcdcorrectedlist[iqcd], emptysystlist, "noflav", 1);
        delete tree_;

        TreeReader * tree_ = new TreeReader(CorrOption, tree, qcdcorrectedlist[iqcd], 2);
        tree_.Loop(CorrOption, datalist,  mclist, qcdcorrectedlist[iqcd], emptysystlist, "noflav", 2);
        delete tree_;

    }

  }
*/
  if(CorrOption == 3 && systlist.size() > 1)
  {
/*
      for(unsigned short int itt = 0; itt < TTWsystlist.size(); itt++)
      {

          if (TTWsystlist[itt] == "TTMSDecays_scaledown" || TTWsystlist[itt] == "TTMSDecays_scaleup" || TTWsystlist[itt] == "TTMSDecays_matchingdown" || TTWsystlist[itt] == "TTMSDecays_matchingup")
          {
              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 0);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "noflav", 0);
              delete tree_;
          }
          else if (TTWsystlist[itt] == "WJets_scaledown")
          {
              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "b", -2);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "c", -2);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "l", -2);
              delete tree_;
          }
          else if (TTWsystlist[itt] == "WJets_matchingdown")
          {
              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "b", -1);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "c", -1);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], -1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "l", -1);
              delete tree_;
          }
          else if (TTWsystlist[itt] == "WJets_matchingup")
          {
              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "b", 1);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "c", 1);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 1);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "l", 1);
              delete tree_;
          }
          else if (TTWsystlist[itt] == "WJets_scaleup")
          {
              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "b", 2);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "c", 2);
              delete tree_;

              TreeReader * tree_ = new TreeReader(CorrOption, tree, TTWsystlist[itt], 2);
              tree_.Loop(CorrOption, datalist,  mclist, TTWsystlist[itt], emptysystlist, "l", 2);
              delete tree_;
          }



      }
*/
  }

}
