#define TreeReader_cxx
#include "TreeReader.h"
#include "common.h"
#include <string>
#include <TGraphAsymmErrors.h>



void TreeReader::Loop(short int CorrOption, vector<TString> datalist, vector<TString> mclist, TString sample, vector<TString> systlist, TString flavtag, short int systType)
{
   TString thesample  = sample;
   TFile * theoutputfile = 0;
   isData       = (sample == "SingleElA" || sample == "SingleElB" || sample == "SingleElC" || sample == "SingleElD" || sample == "SingleMuA" || sample == "SingleMuB" || sample == "SingleMuC" || sample == "SingleMuD" || sample == "NTuple_53X_SingleMuRun2012A" || sample == "NTuple_53X_SingleMuRun2012B" || sample == "NTuple_53X_SingleMuRun2012C" || sample == "NTuple_53X_SingleMuRun2012D" );
   isQCD        = (sample == "QCD" || sample == "QCDdatadriven" );
   isW          = (sample == "WJets"  || sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets" || sample == "WJets_matchingdown" || sample == "WJets_matchingup" || sample == "WJets_scaledown" ||  sample == "WJets_scaleup" );
   isWExcl      = (sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets" || sample == "WJets_matchingdown" || sample == "WJets_matchingup" || sample == "WJets_scaledown" ||  sample == "WJets_scaleup" );
   isWIncl      = (sample == "WJets");
   isTTbarSyst  = (sample == "TTMSDecays_matchingup" || sample == "TTMSDecays_matchingdown" || sample == "TTMSDecays_scaledown" ||  sample == "TTMSDecays_scaleup" );
   isSignal     = ( (sample.SubString("S1Res") == "S1Res") || (sample.SubString("S4Inv") == "S4Inv") );

   TString sampleroot;
   if      ( isW && flavtag == "b") sampleroot = thesample+"_bflav";
   else if ( isW && flavtag == "c") sampleroot = thesample+"_cflav";
   else if ( isW && flavtag == "l") sampleroot = thesample+"_lflav";
   else                             sampleroot = thesample;
   double qcdIsoCut = 0.5;

   if(CorrOption == 3 && !isData && !isQCD && !isTTbarSyst) loadPDFSumWeight(thechannel,  sample, flavtag);
   //cout << "thechannel = " << thechannel << " | sample = " << sample << " | flavtag = " << flavtag << endl;

   if      (CorrOption == 0)  theoutputfile = new TFile( ("outputroot_QCDcorr_iso0p5/histofile_"     +sampleroot+".root" ).Data() , "recreate");
   else if (CorrOption == 1)  theoutputfile = new TFile( ("outputroot_woQCDcorr/histofile_"   +sampleroot+".root" ).Data() , "recreate");
   else if (CorrOption == 2)
   {
       if (systType == -2 && isQCD)      theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_"    +sampleroot+"__BgdContam__minus.root" ).Data() , "recreate");
       else if (systType == -1 && isQCD)  { theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_"    +sampleroot+"__Iso__minus.root" ).Data() , "recreate"); qcdIsoCut = 0.6; }
       else if (systType == 1 && isQCD)   { theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_"    +sampleroot+"__Iso__plus.root"  ).Data() , "recreate"); qcdIsoCut = 0.4; }
       else if (systType == 2 && isQCD)  theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_"    +sampleroot+"__BgdContam__plus.root"  ).Data() , "recreate");
       else                       theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_"    +sampleroot+".root"        ).Data() , "recreate");
   }
   else if (CorrOption == 3)
   {
       if (systType == -2 && isQCD)      theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__BgdContam__minus.root" ).Data() , "recreate");
       else if (systType == -1 && isQCD) { theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__Iso__minus.root" ).Data() , "recreate"); qcdIsoCut = 0.6; }
       else if (systType == 1 && isQCD)  { theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__Iso__plus.root"  ).Data() , "recreate"); qcdIsoCut = 0.4; }
       else if (systType == 2 && isQCD)  theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__BgdContam__plus.root"  ).Data() , "recreate");
       else                       theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+".root"        ).Data() , "recreate");
   }
   else cout << "ERROR: Wrong value of CorrOption or flavtag! Allowed values: 0,1,2,3 for corrOption and Incl, b, c, l for flavtag (if CorrOption == 3)" << endl;

   for(unsigned int i=0; i< systlist.size(); i++)
   {
     TString samplename = "";
     if(      systlist[i] == "" && isQCD && systType == -2) samplename = thesample+"__BgdContam__minus";
     else if( systlist[i] == "" && isQCD && systType == -1) samplename = thesample+"__Iso__minus";
     else if( systlist[i] == "" && isQCD && systType ==  1) samplename = thesample+"__Iso__plus";
     else if( systlist[i] == "" && isQCD && systType ==  2) samplename = thesample+"__BgdContam__plus";
     else if( systlist[i] == "")                            samplename = thesample;
     else                                                   samplename = thesample+"__"+systlist[i];

     bool     firstinit = false;
     if(i==0) firstinit = true;
     initializeHisto(samplename, systlist[i], firstinit, flavtag, systType);
   }


   SetUpCSVreweighting();

   //-------------------------------
   //determine the btag scale factor
   //-------------------------------

   //fBTagSF = new BTagSFUtil("CSVL","ABCD");//ReReco
   fBTagSF = new BTagSFUtil("CSVM","ABCD");//ReReco
   //fBTagSF = new BTagSFUtil("CSVT","ABCD");//ReReco

   if (fChain == 0) return;

   cout << endl;
   cout << "Processing the events... " << endl;
   cout << endl;

   Long64_t nentries = fChain->GetEntriesFast();

   double SF_QCD_W = 0, SF_QCD_S = 0, SF_QCD_TT = 0, SF_QCD_L = 0;
   if ((CorrOption == 2 || CorrOption == 3) && isQCD)
   {
       SF_QCD_L  = getQCDscalefactor(datalist, mclist, "mWT_"+thechannel+"_QCDnormLregion" );
       SF_QCD_W  = getQCDscalefactor(datalist, mclist, "mWT_"+thechannel+"_QCDnormWregion" );
       SF_QCD_S  = getQCDscalefactor(datalist, mclist, "mWT_"+thechannel+"_QCDnormSregion" );
       SF_QCD_TT = getQCDscalefactor(datalist, mclist, "mWT_"+thechannel+"_QCDnormTTregion");
   }


   TGraphAsymmErrors* ratio_eta_0to0p9;
   TGraphAsymmErrors* ratio_eta_0p9to1p2;
   TGraphAsymmErrors* ratio_eta_1p2to2p1;


   TFile* ftrigger = new TFile(string("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root").c_str());
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD",ratio_eta_0to0p9);
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD",ratio_eta_0p9to1p2);
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD",ratio_eta_1p2to2p1);
   ftrigger->Close();

   TRandom3* rand_ = new TRandom3(12345);

   Long64_t nbytes = 0, nb = 0;
   double n0b = 0, n1b = 0, n2b = 0, n3b = 0, dilepton = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (nentries > 10000 && (jentry % (nentries / 10000) == 0)) 			            printProgressBar(jentry,nentries,sample);
      if (nentries > 100   && nentries <= 10000 && (jentry % (nentries / 100) == 0)) 	printProgressBar(jentry,nentries,sample);

      bool hasBjet  = false, hasCjet  = false, isBflavour = false, isCflavour = false, isLflavour = false, match = false;
      short int nbflav = 0;
      if(CorrOption != 0 && isW && (flavtag != "Incl" && flavtag != "allflav") )
      {
         int     iter_jets   = smalltree_njets;
         float * jet_pt	     = smalltree_jet_pt;
         float * jet_eta     = smalltree_jet_eta;
         int   * jet_flav    = smalltree_jet_flav;

         for(short int ijet=0; ijet<iter_jets; ijet++)
          {
              if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4)  continue;
              if(abs(jet_flav[ijet]) == 5)                       {hasBjet = true; nbflav++;}
              if(abs(jet_flav[ijet]) == 4)                        hasCjet = true;
          }

          //---------------
          //W flavour study
          //---------------
          if      ( hasBjet)             isBflavour = true;
          else if (!hasBjet &&  hasCjet) isCflavour = true;
          else if (!hasBjet && !hasCjet) isLflavour = true;
          else cout << "WARNING: This event has been rejected!" << endl;

          //---------------------------------------
          //Test the match with the output filename
          //---------------------------------------
          if (flavtag == "b" && isBflavour) match = true;
          if (flavtag == "c" && isCflavour) match = true;
          if (flavtag == "l" && isLflavour) match = true;

      }
      else if (CorrOption != 0 && isW && flavtag == "allflav") match = true;

      //---------------------------------------------------------------------------
      //If the event content differs from the output filename then reject the event
      //---------------------------------------------------------------------------
      if(CorrOption != 0 && isW && !match) continue;
      if(CorrOption == 0 && smalltree_lept_iso[0] < 0.5) continue;
      if(CorrOption == 3 && isQCD && smalltree_lept_iso[0] < qcdIsoCut) continue;

      if ( nbflav == 0) n0b++;
      if ( nbflav == 1) n1b++;
      if ( nbflav == 2) n2b++;
      if ( nbflav >= 3) n3b++;

      double SFtrigger      = 1;
      double SFtriggerError = 0;
      if(!useElectronChannel)
      {
        if (abs(smalltree_lept_eta[0]) < 0.9)
        {
          SFtrigger      = getSFtrigger(ratio_eta_0to0p9,  smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_0to0p9,  smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
        }
        else if (abs(smalltree_lept_eta[0]) >= 0.9 && abs(smalltree_lept_eta[0]) < 1.2)
        {
          SFtrigger      = getSFtrigger(ratio_eta_0p9to1p2,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_0p9to1p2,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
        }
        else if (abs(smalltree_lept_eta[0]) >= 1.2 && abs(smalltree_lept_eta[0]) <= 2.1)
        {
          SFtrigger      = getSFtrigger(ratio_eta_1p2to2p1,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_1p2to2p1,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
        }
      }

      //----------------------------------------------------------------------
      //apply event selection
      //fifth argument is for the systematics
      //"" means no systematics
      //else, enter the systematic strings that appears in the histograms name
      //systematic names are :
      //"leptup", "leptdown", "trigup", "trigdown", "PDFup", "PDFdown"
      //"jesup", "jesdown", "jerup", "jerdown", "metunclsup", "metunclsdown"
      //----------------------------------------------------------------------

      double randomNumber = rand_->Uniform(1.);

      if      (CorrOption == 0)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, systlist[0], thesample, flavtag, systType, randomNumber);
      else if (CorrOption == 1)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, systlist[0], thesample, flavtag, systType, randomNumber);
      else if (CorrOption == 2)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, systlist[0], thesample, flavtag, systType, randomNumber);
      else if (CorrOption == 3)
      {
          for(unsigned int isyst=0; isyst< systlist.size(); isyst++)
          {
              if( (systlist[isyst] == "toppt__plus" || systlist[isyst] == "toppt__minus") && thesample != "TTMSDecays_central" ) continue;
              if( isData && systlist[isyst] != "") continue;
              applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, systlist[isyst], thesample, flavtag, systType, randomNumber);
          }
      }
      else cout << "ERROR: Wrong value of CorrOption! Allowed valued: 0,1,2,3" << endl;

   }

   manageOverflows();

   cout << endl;
   cout << "N0b = " << n0b << " | N1b = " << n1b << " | N2b = " << n2b << " | N3b = " << n3b << endl;
   theoutputfile->Write();
   deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;

}



bool TreeReader::applyEventSel(short int CorrOption, double SF_QCD_L, double SF_QCD_W, double SF_QCD_S, double SF_QCD_TT, double SFtrigger, double SFtriggerError, TString systtype, TString sample, TString flavtag, short int systType, double randomNumber){

      bool applyCSV_reshape = 0;
      bool applyCSV         = 1;

      int  btagSys          = 0;
      int  mistagSys        = 0;

      double met_pt    = smalltree_met_pt;
      double met_phi   = smalltree_met_phi;
      double evtweight = -1;

      TString thesample = sample;
      if(systtype != "")               thesample = thesample + "__"+systtype;
      if(isQCD && systType == -2)      thesample = thesample + "__BgdContam__minus";
      else if(isQCD && systType == -1) thesample = thesample + "__Iso__minus";
      else if(isQCD && systType ==  1) thesample = thesample + "__Iso__plus";
      else if(isQCD && systType ==  2) thesample = thesample + "__BgdContam__plus";

      TString WtagForPDF = sample;
      if      ( isW && flavtag == "b") WtagForPDF = "WExclb";
      else if ( isW && flavtag == "c") WtagForPDF = "WExclc";
      else if ( isW && flavtag == "l") WtagForPDF = "WExcll";

      int iter_jets          = 0;
      float * jet_pt         = 0;
      float * jet_eta        = 0;
      float * jet_phi        = 0;
      float * jet_btagdiscri = 0;
      int   * jet_flav       = 0;

      double wCSV = 1.;
      wCSV = GetCSVweight(0,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      evtweight = smalltree_evtweight;
      iter_jets      = smalltree_njets;
      jet_pt	     = smalltree_jet_pt;
      jet_eta	     = smalltree_jet_eta;
      jet_phi	     = smalltree_jet_phi;
      jet_btagdiscri = smalltree_jet_btagdiscri;
      jet_flav       = smalltree_jet_flav;


      if(systtype == "" ||
         systtype == "lept__plus"  ||   systtype == "lept__minus"  ||
         systtype == "trig__plus"  ||   systtype == "trig__minus"  ||
         systtype == "PU__plus"    ||   systtype == "PU__minus"    ||
         systtype == "toppt__plus" ||   systtype == "toppt__minus" ||
         systtype == "PDF__plus"   ||   systtype == "PDF__minus"

      ){
	iter_jets      = smalltree_njets;
        jet_pt         = smalltree_jet_pt;
        jet_eta        = smalltree_jet_eta;
        jet_phi        = smalltree_jet_phi;
        jet_btagdiscri = smalltree_jet_btagdiscri;
        jet_flav       = smalltree_jet_flav;

	if(!isQCD) evtweight = smalltree_evtweight;
    else       evtweight = smalltree_evtweight_nominal;


    if(     systtype == "lept__plus")                           evtweight = smalltree_weight_leptup;
	else if(systtype == "lept__minus")                          evtweight = smalltree_weight_leptdown;
	else if(systtype == "trig__plus" && useElectronChannel)     evtweight = smalltree_weight_trigup;
	else if(systtype == "trig__plus")                           evtweight = (SFtrigger + SFtriggerError)*smalltree_evtweight/SFtrigger;
	else if(systtype == "trig__minus" && useElectronChannel)    evtweight = smalltree_weight_trigdown;
	else if(systtype == "trig__minus")                          evtweight = (SFtrigger - SFtriggerError)*smalltree_evtweight/SFtrigger;
	else if(systtype == "PU__plus")                             evtweight = smalltree_weight_PUup;
	else if(systtype == "PU__minus")                            evtweight = smalltree_weight_PUdown;
	else if(systtype == "PDF__plus")                            evtweight = smalltree_weight_PDFup;
	//else if(systtype == "PDF__plus")                            evtweight = smalltree_weight_PDFup*(PDFSumWeight["pdfWeight_"+thechannel+"___"+WtagForPDF+"__"+systtype]);
	else if(systtype == "PDF__minus")                           evtweight = smalltree_weight_PDFdown;
	//else if(systtype == "PDF__minus")                           evtweight = smalltree_weight_PDFdown*(PDFSumWeight["pdfWeight_"+thechannel+"___"+WtagForPDF+"__"+systtype]);
	else if(systtype == "toppt__plus")                          evtweight = smalltree_weight_toppt*smalltree_evtweight;
	else if(systtype == "toppt__minus")                         evtweight = smalltree_evtweight;

      }else if(systtype == "jes__plus"){

	met_pt         = smalltree_met_jesup_pt;
	met_phi        = smalltree_met_jesup_phi;
	iter_jets      = smalltree_jesup_njets;
        jet_pt         = smalltree_jet_jesup_pt;
        jet_eta        = smalltree_jet_jesup_eta;
        jet_phi        = smalltree_jet_jesup_phi;
        jet_btagdiscri = smalltree_jet_jesup_btagdiscri;
        jet_flav       = smalltree_jet_jesup_flav;

      }else if(systtype == "jes__minus"){
	met_pt         = smalltree_met_jesdown_pt;
	met_phi        = smalltree_met_jesdown_phi;
	iter_jets      = smalltree_jesdown_njets;
        jet_pt         = smalltree_jet_jesdown_pt;
        jet_eta        = smalltree_jet_jesdown_eta;
        jet_phi        = smalltree_jet_jesdown_phi;
        jet_btagdiscri = smalltree_jet_jesdown_btagdiscri;
        jet_flav       = smalltree_jet_jesdown_flav;

      }else if(systtype == "jer__plus"){
	met_pt         = smalltree_met_jerup_pt;
	met_phi        = smalltree_met_jerup_phi;
	iter_jets      = smalltree_jerup_njets;
        jet_pt         = smalltree_jet_jerup_pt;
        jet_eta        = smalltree_jet_jerup_eta;
        jet_phi        = smalltree_jet_jerup_phi;
        jet_btagdiscri = smalltree_jet_jerup_btagdiscri;
        jet_flav       = smalltree_jet_jerup_flav;

      }else if(systtype == "jer__minus"){
	met_pt         = smalltree_met_jerdown_pt;
	met_phi        = smalltree_met_jerdown_phi;
	iter_jets      = smalltree_jerdown_njets;
        jet_pt         = smalltree_jet_jerdown_pt;
        jet_eta        = smalltree_jet_jerdown_eta;
        jet_phi        = smalltree_jet_jerdown_phi;
        jet_btagdiscri = smalltree_jet_jerdown_btagdiscri;
        jet_flav       = smalltree_jet_jerdown_flav;

      }else if(systtype == "metuncls__plus"){
	met_pt         = smalltree_met_unclsup_pt;
	met_phi        = smalltree_met_unclsup_phi;

      }else if(systtype == "metuncls__minus"){
	met_pt         = smalltree_met_unclsdown_pt;
	met_phi        = smalltree_met_unclsdown_phi;

      }else if(systtype == "btag__plus"   ){btagSys =  1;}
      else if(systtype  == "btag__minus"  ){btagSys = -1;}
      else if(systtype  == "mistag__plus" ){mistagSys =  1;}
      else if(systtype  == "mistag__minus"){mistagSys = -1;}
      else{
        cout << "WARNING syst type " << systtype << " not recognized !! " << endl;
	    cout << "correct syst types are " << endl;
	    cout << " \"\",  \"lept__plus\", \"lept__minus\", \"trig__plus\", \"trig__minus\", \"PDF__plus\", \"PDF__minus\"  " << endl;
	    cout << "\"jes__plus\", \"jes__minus\", \"jer__plus\", \"jer__minus\", \"metuncls__plus\", \"metuncls__minus\" " << endl;
      }

      if( applyCSV_reshape && !applyCSV ) evtweight *= wCSV;

      bool isBsyst = false;
      //determine the btag systematic
      if(systtype == "btag__plus" )
      {
        btagSys   = 1;
        evtweight = smalltree_evtweight;
        isBsyst = true;
      }
      else if(systtype == "btag__minus" )
      {
        btagSys   = -1;
        evtweight = smalltree_evtweight;
        isBsyst = true;
      }
      else if(systtype == "mistag__plus" )
      {
        mistagSys = 1;
        evtweight = smalltree_evtweight;
        isBsyst = true;
      }
      else if(systtype == "mistag__minus")
      {
        mistagSys = -1;
        evtweight = smalltree_evtweight;
        isBsyst = true;
      }

     // Add the trigger efficiency
     if(isData)      evtweight  = 1;
     else if(!isQCD) evtweight *= SFtrigger;


     double QCDContamweight = 1;
     if  (CorrOption == 2 || CorrOption == 3 )
     {
         if(sample == "QCD" && systType == 0)         QCDContamweight = smalltree_evtweight_nominal;
         else if(sample == "QCD" && systType == -2)   QCDContamweight = smalltree_evtweight_minus;
         else if(sample == "QCD" && systType ==  2)   QCDContamweight = smalltree_evtweight_plus;
         if(sample == "QCD" )                         evtweight = QCDContamweight;
     }


     if (systtype == "PDF__plus"  && (!isfinite(smalltree_weight_PDFup)   || isnan(smalltree_weight_PDFup)   )) return false;
     if (systtype == "PDF__minus" && (!isfinite(smalltree_weight_PDFdown) || isnan(smalltree_weight_PDFdown) )) return false;

     if (!isfinite(evtweight) || isnan(evtweight))
     {
         cout << "WARNING: Infinite/NaN value of evtweight !! (syst = " << systtype << ")" << endl;
         return false;
     }

     fillHisto("CutFlow",     "",  thesample, systtype, 1 ,  evtweight, flavtag, systType);

     fillHisto("sumWeight",    "", thesample, systtype, 0.5 , evtweight , flavtag, systType);
     if( systtype == "PDF__plus" || systtype == "PDF__minus") fillHisto("pdfWeight",   "", thesample, systtype, 0.5 , evtweight, flavtag, systType);


     if(smalltree_lept_pt[0] > 33 && fabs(smalltree_lept_eta[0]) < 2.1 )
     {
       fillHisto("CutFlow",     "",  thesample, systtype, 2 ,  evtweight, flavtag, systType);
       TLorentzVector lept, met, leadingJet;
       lept.SetPtEtaPhiM(smalltree_lept_pt[0], smalltree_lept_eta[0],  smalltree_lept_phi[0], 0.);
       met.SetPtEtaPhiM(met_pt, 0, met_phi , 0.);
       double mTW = pow( 2*smalltree_lept_pt[0]*met_pt*(1-cos(smalltree_lept_phi[0] -  met_phi)) ,0.5);

       if(iter_jets>0) leadingJet.SetPtEtaPhiM(jet_pt[0], jet_eta[0], jet_phi[0], 0);
       if(CorrOption != 0 && (sample == "QCD"))          evtweight = QCDContamweight*SF_QCD_L;

       int njets=0;
       int nbjets = 0;
       bool jetsup70 = false;
       vector<int > btagged_jet_idx;
       fillHisto("LeadJetBtagDiscr",     "afterleptsel",  thesample, systtype,  jet_btagdiscri[0] ,  evtweight, flavtag, systType);
       for(int ijet=0; ijet<iter_jets; ijet++)
       {
         if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
         njets++;
         fillHisto("JetPt",     "afterleptsel",  thesample, systtype,  jet_pt[ijet] ,  evtweight, flavtag, systType);
         fillHisto("JetEta",    "afterleptsel",  thesample, systtype,  jet_eta[ijet] , evtweight, flavtag, systType);
         //if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.244)  {nbjets++;   btagged_jet_idx.push_back(ijet); }
         if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.679)  {nbjets++;   btagged_jet_idx.push_back(ijet); }
         //if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.898)  {nbjets++;   btagged_jet_idx.push_back(ijet); }

         if( applyCSV)
         {
           bool isbtag = 0;
           if(isData) isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], -999999 , jet_pt[ijet] , jet_eta[ijet], btagSys, randomNumber); // for nominal sample
           else
           {
               if(abs(jet_flav[ijet]) == 5 || abs(jet_flav[ijet]) == 4)
               {
                   isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], jet_flav[ijet] , jet_pt[ijet] , jet_eta[ijet], btagSys, randomNumber);
               }
               else
               {
                   isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], jet_flav[ijet] , jet_pt[ijet] , jet_eta[ijet], mistagSys, randomNumber);
               }
           }
           fillHisto("BTagProba", "afterleptsel", thesample, systtype, isbtag , evtweight, flavtag, systType);
           if(isbtag) {nbjets++; btagged_jet_idx.push_back(ijet); }
         }

	     if(jet_pt[ijet] > 70 ) jetsup70 = true;
       }

       if (njets > 0) fillHisto("CutFlow",     "",  thesample, systtype, 3 ,  evtweight, flavtag, systType);

       double term1 = lept.Pz()*(lept.Px()*met.Px() + lept.Py()*met.Py() + (80.399)*(80.399)/2.);
       double det = lept.Px()*met.Px() + lept.Py()*met.Py() + (80.399)*(80.399)/2.
       		  - met.Pt()*met.Pt()*(lept.E()*lept.E() - lept.Pz()*lept.Pz());
       if(det<0) det=0;
       double term2 = lept.E()*pow(det, 0.5);
       double denom = lept.E()*lept.E() - lept.Pz()*lept.Pz();

       double sol1 = (term1 - term2)/denom;

       double neutrE = pow( pow(met.Px(),2) + pow(met.Py(),2) + sol1*sol1, 0.5);//neglecting neut mass
       double mWTplusMET = mTW + met_pt;
       TLorentzVector neutrino;
       neutrino.SetPxPyPzE( met.Px(), met.Py(), sol1, neutrE);

       TLorentzVector theWcand = neutrino   + lept;
       TLorentzVector jetLept  = leadingJet + lept;

       fillHisto("NJet",       "afterleptsel",  thesample, systtype,   njets,   		        evtweight, flavtag, systType);
       fillHisto("NBJet",      "afterleptsel",  thesample, systtype,   nbjets, 		        evtweight, flavtag, systType);

       fillHisto("ptW",        "afterleptsel",  thesample, systtype,   theWcand.Pt(),       	evtweight, flavtag, systType);
       fillHisto("etaW",       "afterleptsel",  thesample, systtype,   theWcand.Eta(),       	evtweight, flavtag, systType);

       fillHisto("mWT",        "afterleptsel",  thesample, systtype,   mTW,    		        evtweight, flavtag, systType);
       fillHisto("mW",         "afterleptsel",  thesample, systtype,   theWcand.M(),  	        evtweight, flavtag, systType);
       fillHisto("MET",        "afterleptsel",  thesample, systtype,   met_pt, 		        evtweight, flavtag, systType);
       fillHisto("mWTplusMET", "afterleptsel",  thesample, systtype,   mWTplusMET, 		        evtweight, flavtag, systType);
       fillHisto("LeptPt",     "afterleptsel",  thesample, systtype,   smalltree_lept_pt[0], 	evtweight, flavtag, systType);
       fillHisto("LeptEta",    "afterleptsel",  thesample, systtype,   smalltree_lept_eta[0],       evtweight, flavtag, systType);

        //***********************************
        //qcd region
        //***********************************
    if(CorrOption == 0)
    {

        fillHisto("mWT",          "qcdLregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
        fillHisto("MET",          "qcdLregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);

        if(jetsup70 && njets == 1 && nbjets == 0)
        {
          fillHisto("mWT",        "qcdWregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto("MET",        "qcdWregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        else if(njets == 1 && nbjets == 1)
        {
          fillHisto("mWT",        "qcdBregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto("MET",        "qcdBregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        if(jetsup70 && njets == 1 && nbjets == 1)
        {
          fillHisto("mWT",        "qcdSregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto("MET",        "qcdSregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        else if(jetsup70 && njets >= 4 && nbjets ==2 )
        {
          fillHisto("mWT",        "qcdTTregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto("MET",        "qcdTTregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
	}

        //***********************************
        //QCD normalization region
        //***********************************
    if(CorrOption == 0 || CorrOption == 1)
    {

        if(mTW <= 40) fillHisto("mWT",        "QCDnormLregion",   thesample, systtype,   mTW,    		evtweight, flavtag, systType);

        if(jetsup70 && njets == 1 && nbjets == 0 && mTW <= 40)
        {
          fillHisto("mWT",        "QCDnormWregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        else if(njets == 1 && nbjets == 1 && mTW <= 40)
        {
          fillHisto("mWT",        "QCDnormBregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        if(jetsup70 && njets == 1 && nbjets == 1 && mTW <= 40)
        {
          fillHisto("mWT",        "QCDnormSregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        else if(jetsup70 && njets >= 4 && nbjets == 2 && mTW <= 40)
        {
          fillHisto("mWT",        "QCDnormTTregion",  thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
	}

    if(CorrOption == 1 || CorrOption == 2 || CorrOption == 3)
    {
        //***********************************
        //W enriched region
        //***********************************
	    if(jetsup70 && njets == 1 && nbjets == 0)
        {

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_W;
          fillHisto("mWT_full",   "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 12 ,  evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",  "Wregion_highpt",  thesample, systtype,  jet_pt[ijet] ,                 evtweight, flavtag, systType);
                fillHisto("JetEta", "Wregion_highpt",  thesample, systtype,  jet_eta[ijet] ,                evtweight, flavtag, systType);
            }
            fillHisto("NVtx",       "Wregion_highpt",  thesample, systtype,   smalltree_nvertex,                evtweight, flavtag, systType);
            fillHisto("NJet",       "Wregion_highpt",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto("NBJet",      "Wregion_highpt",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto("LeadJetBtagDiscr","Wregion_highpt",thesample, systtype,jet_btagdiscri[0],                evtweight, flavtag, systType);
            fillHisto("ptW",        "Wregion_highpt",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto("etaW",       "Wregion_highpt",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);
            fillHisto("mWT",        "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("mWT_high",   "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("MET",        "Wregion_highpt",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ", "Wregion_highpt",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),   evtweight, flavtag, systType);
            fillHisto("DeltaPhiWptMet",  "Wregion_highpt",  thesample, systtype,  abs(theWcand.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto("DeltaPhiLMet","Wregion_highpt", thesample, systtype,   abs(lept.DeltaPhi(met)),          evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",   "Wregion_highpt",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);
            fillHisto("LeptPt",     "Wregion_highpt",  thesample, systtype,   smalltree_lept_pt[0], 	        evtweight, flavtag, systType);
            fillHisto("LeptEta",    "Wregion_highpt",  thesample, systtype,   smalltree_lept_eta[0], 	        evtweight, flavtag, systType);
          }
	    }

        //***********************************
        //Intermediate region
        //***********************************
	    if(jetsup70 && njets == 1 && nbjets == 1 && met_pt <= 100)
        {

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_S;
          fillHisto("mWT_full",   "interRegion",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 15 ,  evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",  "interRegion",  thesample, systtype,  jet_pt[ijet] ,                 evtweight, flavtag, systType);
                fillHisto("JetEta", "interRegion",  thesample, systtype,  jet_eta[ijet] ,                evtweight, flavtag, systType);
            }
            fillHisto("NVtx",       "interRegion",  thesample, systtype,   smalltree_nvertex,                evtweight, flavtag, systType);
            fillHisto("NJet",       "interRegion",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto("NBJet",      "interRegion",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto("LeadJetBtagDiscr","interRegion",thesample, systtype,jet_btagdiscri[0],                evtweight, flavtag, systType);
            fillHisto("ptW",        "interRegion",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto("etaW",       "interRegion",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);
            fillHisto("mWT",        "interRegion",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("mWT_high",   "interRegion",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("MET",        "interRegion",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ", "interRegion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),   evtweight, flavtag, systType);
            fillHisto("DeltaPhiWptMet",  "interRegion",  thesample, systtype,  abs(theWcand.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto("DeltaPhiLMet","interRegion", thesample, systtype,   abs(lept.DeltaPhi(met)),          evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",   "interRegion",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);
            fillHisto("LeptPt",     "interRegion",  thesample, systtype,   smalltree_lept_pt[0], 	        evtweight, flavtag, systType);
            fillHisto("LeptEta",    "interRegion",  thesample, systtype,   smalltree_lept_eta[0], 	        evtweight, flavtag, systType);
          }
	    }


        //***********************************
        //FCNC signal enriched region
        //***********************************
        if(mTW > 40 && njets > 0)
        {
            fillHisto("CutFlow",     "",  thesample, systtype, 4 ,  evtweight, flavtag, systType);

            if(nbjets == 1)
            {
                fillHisto("CutFlow",     "",  thesample, systtype, 5 ,  evtweight, flavtag, systType);
                if(njets == 1)
                {
                    fillHisto("CutFlow",     "",  thesample, systtype, 6 ,  evtweight, flavtag, systType);
                    if(jetsup70)
                    {
                        fillHisto("ptW",            "Signalregion",  thesample, systtype,   theWcand.Pt(),       evtweight, flavtag, systType);
                        fillHisto("MET",            "Signalregion",  thesample, systtype,   met_pt, 			evtweight, flavtag, systType);
                        fillHisto("DeltaPhiLJ",     "Signalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
                        fillHisto("mWT",            "Signalregion",  thesample, systtype,   mTW,    			evtweight, flavtag, systType);

                        fillHisto("CutFlow",     "",  thesample, systtype, 7 ,  evtweight, flavtag, systType);
                        if(met_pt > 100)
                        {
                            fillHisto("CutFlow",     "",  thesample, systtype, 8 ,  evtweight, flavtag, systType);
                            if(theWcand.Pt() > 50)
                            {
                                fillHisto("CutFlow",     "",  thesample, systtype, 9 ,  evtweight, flavtag, systType);
                                if(abs(lept.DeltaPhi(leadingJet)) < 1.7)
                                {
                                    fillHisto("CutFlow",     "",  thesample, systtype, 10 ,  evtweight, flavtag, systType);
                                    if(sample == "QCD")          evtweight = QCDContamweight*SF_QCD_S;

                                    for(int ijet=0; ijet<iter_jets; ijet++)
                                    {
                                        if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                                        fillHisto("JetPt",      "Selectedsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                                        fillHisto("JetEta",     "Selectedsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
                                    }

                                        tree_EvtWeight  = evtweight;
                                        tree_wMass      = mTW;
                                        tree_met        = met_pt;
                                        tree_wPt        = theWcand.Pt();
                                        tree_deltaPhilb = abs(lept.DeltaPhi(leadingJet));

                                        if(theTree_map[thesample] != 0)  theTree_map[thesample]->Fill();

                                        //if tmeme is dilepton
                                        //if(smalltree_tmeme != 0 || smalltree_tmeme != 1 || smalltree_tmeme != 10 || smalltree_tmeme != 11000 || smalltree_tmeme != 10100 ) dilepton++;

                                        fillHisto("NVtx",           "Selectedsignalregion",  thesample, systtype,   smalltree_nvertex,   evtweight, flavtag, systType);
                                        fillHisto("LeadJetBtagDiscr","Selectedsignalregion", thesample, systtype,   jet_btagdiscri[0] ,  evtweight, flavtag, systType);
                                        fillHisto("NJet",           "Selectedsignalregion",  thesample, systtype,   njets, 	    	     evtweight, flavtag, systType);
                                        fillHisto("NBJet",          "Selectedsignalregion",  thesample, systtype,   nbjets, 		     evtweight, flavtag, systType);
                                        fillHisto("LeptPt",         "Selectedsignalregion",  thesample, systtype,   smalltree_lept_pt[0],evtweight, flavtag, systType);
                                        fillHisto("LeptEta",        "Selectedsignalregion",  thesample, systtype,   smalltree_lept_eta[0],evtweight, flavtag, systType);
                                        fillHisto("ptW",            "Selectedsignalregion",  thesample, systtype,   theWcand.Pt(),       evtweight, flavtag, systType);
                                        fillHisto("etaW",           "Selectedsignalregion",  thesample, systtype,   theWcand.Eta(),      evtweight, flavtag, systType);
                                        fillHisto("mWT",            "Selectedsignalregion",  thesample, systtype,   mTW,    			evtweight, flavtag, systType);
                                        fillHisto("mWT_high",       "Selectedsignalregion",  thesample, systtype,   mTW,    			evtweight, flavtag, systType);
                                        fillHisto("MET",            "Selectedsignalregion",  thesample, systtype,   met_pt, 			evtweight, flavtag, systType);
                                        fillHisto("DeltaPhiLJ",     "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
                                        fillHisto("DeltaPhiLMet",   "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
                                        fillHisto("DeltaPhiJMet",   "Selectedsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
                                        fillHisto("DeltaPhiLJMet",  "Selectedsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
                                        fillHisto("DeltaPhiWptMet",  "Selectedsignalregion",  thesample, systtype,  abs(theWcand.DeltaPhi(met)),     evtweight, flavtag, systType);
                                        fillHisto("DeltaRLJ",       "Selectedsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	evtweight, flavtag, systType);
                                }
                            }
                        }
                    }
                }
            }
        }

        //***********************************
        //ATLAS signal enriched region
        //***********************************
        if(njets == 1 && nbjets == 1 ){

          if(sample == "QCD")          evtweight = QCDContamweight*SF_QCD_S;

          if(mWTplusMET > 60 && mTW > 250 && abs(lept.DeltaPhi(leadingJet)) < 1.4)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 13 , evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",      "ATLASFCNCsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                fillHisto("JetEta",     "ATLASFCNCsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
            }
            fillHisto("ptW",            "ATLASFCNCsignalregion",  thesample, systtype,   theWcand.Pt(),       	     evtweight, flavtag, systType);
            fillHisto("mWT",            "ATLASFCNCsignalregion",  thesample, systtype,   mTW,    			     evtweight, flavtag, systType);
            fillHisto("mWT_high",       "ATLASFCNCsignalregion",  thesample, systtype,   mTW,    			     evtweight, flavtag, systType);
            fillHisto("MET",            "ATLASFCNCsignalregion",  thesample, systtype,   met_pt, 			     evtweight, flavtag, systType);
            fillHisto("mWTplusMET",     "ATLASFCNCsignalregion",  thesample, systtype,   mWTplusMET, 		     evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ",     "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
            fillHisto("DeltaPhiLMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
            fillHisto("DeltaPhiJMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJMet",  "ATLASFCNCsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",       "ATLASFCNCsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	     evtweight, flavtag, systType);
          }

          if(mWTplusMET > 60 && mTW > 210 && abs(lept.DeltaPhi(leadingJet)) < 1.2)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 14 , evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",      "ATLASRESsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                fillHisto("JetEta",     "ATLASRESsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
            }
            fillHisto("ptW",            "ATLASRESsignalregion",  thesample, systtype,   theWcand.Pt(),       	    evtweight, flavtag, systType);
            fillHisto("mWT",            "ATLASRESsignalregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
            fillHisto("mWT_high",       "ATLASRESsignalregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
            fillHisto("MET",            "ATLASRESsignalregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
            fillHisto("mWTplusMET",     "ATLASRESsignalregion",  thesample, systtype,   mWTplusMET, 		    evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ",     "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
            fillHisto("DeltaPhiLMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
            fillHisto("DeltaPhiJMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJMet",  "ATLASRESsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",       "ATLASRESsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag, systType);
          }

        }


        //***********************************
        //ttbar (2j2b) enriched region
        //***********************************
        if(jetsup70 && njets == 2 && nbjets == 2){

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_TT;
          fillHisto("mWT_full",     "ttbarregion_2j2b",  thesample, systtype,   mTW,    			        evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 11 ,  evtweight, flavtag, systType);
  	        for(int ijet=0; ijet<iter_jets; ijet++){
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",     "ttbarregion_2j2b",  thesample, systtype,  jet_pt[ijet], 		evtweight, flavtag, systType);
                fillHisto("JetEta",    "ttbarregion_2j2b",  thesample, systtype,  jet_eta[ijet], 		evtweight, flavtag, systType);
            }

     //cout << "evtweight (syst : " << systtype << " ) = " << evtweight << endl;
            fillHisto("NVtx",       "ttbarregion_2j2b",  thesample, systtype,   smalltree_nvertex,                evtweight, flavtag, systType);
            fillHisto("LeadJetBtagDiscr","ttbarregion_2j2b",thesample, systtype,  jet_btagdiscri[0] ,             evtweight, flavtag, systType);
            fillHisto("NJet",       "ttbarregion_2j2b",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto("NBJet",      "ttbarregion_2j2b",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto("ptW",        "ttbarregion_2j2b",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto("etaW",       "ttbarregion_2j2b",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);

            fillHisto("mWT",        "ttbarregion_2j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("mWT_high",   "ttbarregion_2j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("MET",        "ttbarregion_2j2b",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto("LeptPt",     "ttbarregion_2j2b",  thesample, systtype,   smalltree_lept_pt[0], 	        evtweight, flavtag, systType);
            fillHisto("LeptEta",    "ttbarregion_2j2b",  thesample, systtype,   smalltree_lept_eta[0], 	        evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ", "ttbarregion_2j2b",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),   evtweight, flavtag, systType);
            fillHisto("DeltaPhiWptMet",  "ttbarregion_2j2b",  thesample, systtype,  abs(theWcand.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",   "ttbarregion_2j2b",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);

            fillHisto("LeadJetPt",     "ttbarregion_2j2b",  thesample, systtype,  jet_pt[0], 			evtweight, flavtag, systType);
            fillHisto("LeadJetEta",    "ttbarregion_2j2b",  thesample, systtype,  jet_eta[0], 			evtweight, flavtag, systType);

            fillHisto("SubLeadJetPt",  "ttbarregion_2j2b",  thesample, systtype,  jet_pt[1], 			evtweight, flavtag, systType);
            fillHisto("SubLeadJetEta", "ttbarregion_2j2b",  thesample, systtype,  jet_eta[1], 			evtweight, flavtag, systType);


          }
        }

        //***********************************
        //ttbar (3j2b) enriched region
        //***********************************
        if(jetsup70 && njets == 3 && nbjets == 2){

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_TT;
          fillHisto("mWT_full",     "ttbarregion_3j2b",  thesample, systtype,   mTW,    			        evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 11 ,  evtweight, flavtag, systType);
  	        for(int ijet=0; ijet<iter_jets; ijet++){
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",     "ttbarregion_3j2b",  thesample, systtype,  jet_pt[ijet], 		evtweight, flavtag, systType);
                fillHisto("JetEta",    "ttbarregion_3j2b",  thesample, systtype,  jet_eta[ijet], 		evtweight, flavtag, systType);
            }

            fillHisto("NVtx",       "ttbarregion_3j2b",  thesample, systtype,   smalltree_nvertex,                evtweight, flavtag, systType);
            fillHisto("LeadJetBtagDiscr","ttbarregion_3j2b",thesample, systtype,  jet_btagdiscri[0] ,             evtweight, flavtag, systType);
            fillHisto("NJet",       "ttbarregion_3j2b",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto("NBJet",      "ttbarregion_3j2b",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto("ptW",        "ttbarregion_3j2b",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto("etaW",       "ttbarregion_3j2b",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);

            fillHisto("mWT",        "ttbarregion_3j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("mWT_high",   "ttbarregion_3j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("MET",        "ttbarregion_3j2b",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto("LeptPt",     "ttbarregion_3j2b",  thesample, systtype,   smalltree_lept_pt[0], 	        evtweight, flavtag, systType);
            fillHisto("LeptEta",    "ttbarregion_3j2b",  thesample, systtype,   smalltree_lept_eta[0], 	        evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ", "ttbarregion_3j2b",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),   evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",   "ttbarregion_3j2b",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);
          }
        }


        //***********************************
        //ttbar (4j2b) enriched region
        //***********************************
        if(jetsup70 && njets >= 4 && nbjets == 2){

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_TT;
          fillHisto("mWT_full",     "ttbarregion_4j2b",  thesample, systtype,   mTW,    			        evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto("CutFlow",     "",  thesample, systtype, 11 ,  evtweight, flavtag, systType);
  	        for(int ijet=0; ijet<iter_jets; ijet++){
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto("JetPt",     "ttbarregion_4j2b",  thesample, systtype,  jet_pt[ijet], 		evtweight, flavtag, systType);
                fillHisto("JetEta",    "ttbarregion_4j2b",  thesample, systtype,  jet_eta[ijet], 		evtweight, flavtag, systType);
            }

            fillHisto("NVtx",       "ttbarregion_4j2b",  thesample, systtype,   smalltree_nvertex,                evtweight, flavtag, systType);
            fillHisto("LeadJetBtagDiscr","ttbarregion_4j2b",thesample, systtype,  jet_btagdiscri[0] ,             evtweight, flavtag, systType);
            fillHisto("NJet",       "ttbarregion_4j2b",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto("NBJet",      "ttbarregion_4j2b",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto("ptW",        "ttbarregion_4j2b",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto("etaW",       "ttbarregion_4j2b",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);

            fillHisto("mWT",        "ttbarregion_4j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("mWT_high",   "ttbarregion_4j2b",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto("MET",        "ttbarregion_4j2b",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto("LeptPt",     "ttbarregion_4j2b",  thesample, systtype,   smalltree_lept_pt[0], 	        evtweight, flavtag, systType);
            fillHisto("LeptEta",    "ttbarregion_4j2b",  thesample, systtype,   smalltree_lept_eta[0], 	        evtweight, flavtag, systType);
            fillHisto("DeltaPhiLJ", "ttbarregion_4j2b",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),   evtweight, flavtag, systType);
            fillHisto("DeltaRLJ",   "ttbarregion_4j2b",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);

          }
        }


    }


  } //lepton selection

  return true;
}


//------------------------------------------------------
//initialize the historams for the analysis
//------------------------------------------------------


void TreeReader::initializeHisto(TString sample, TString syst, bool isfirstset, TString flavtag, short int systType)
{

//  cout << endl;
//  cout << endl;
//  cout << "Initializing histograms for sample: " << sample << " ... " << endl;

  if(isfirstset)
  {
    numb_histo = 0;
    TH1D * first_emptyHisto = new TH1D("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
    histo_list.push_back(first_emptyHisto);

    numb_histo++;
  }

  addHisto("CutFlow",               "",     sample.Data(), syst.Data(),    20,    0,    20, flavtag, systType);
  addHisto("sumWeight",             "",     sample.Data(), syst.Data(),     1,    0,     1, flavtag, systType);
  addHisto("pdfWeight",             "",     sample.Data(), syst.Data(),     1,    0,     1, flavtag, systType);

  //after lepton selection
  addHisto("NJet",      "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("mWT",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mW",        "afterleptsel", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "afterleptsel", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "afterleptsel", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("MET",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET","afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("JetPt",     "afterleptsel",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",    "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "afterleptsel",  	sample.Data(), syst.Data(),   30,   0.,   300, flavtag, systType);
  addHisto("LeptEta",   "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("BTagProba", "afterleptsel",  	sample.Data(), syst.Data(),   2,  -0.5,   1.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "afterleptsel",  sample.Data(), syst.Data(),  20,     0,     1, flavtag, systType);

  addHisto("NVtx",      "Wregion_highpt",   sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",      "Wregion_highpt",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "Wregion_highpt",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",     "Wregion_highpt",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",    "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",       "Wregion_highpt", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "Wregion_highpt", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",       "Wregion_highpt",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",  "Wregion_highpt",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("LeptPt",    "Wregion_highpt",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeptEta",   "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet","Wregion_highpt",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiWptMet","Wregion_highpt",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "Wregion_highpt",sample.Data(), syst.Data(),   20,    0,     1, flavtag, systType);

  addHisto("NVtx",      "interRegion",   sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",      "interRegion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "interRegion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",     "interRegion",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",    "interRegion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",       "interRegion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "interRegion", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",       "interRegion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",  "interRegion",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",  "interRegion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "interRegion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("LeptPt",    "interRegion",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeptEta",   "interRegion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","interRegion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet","interRegion",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiWptMet","interRegion",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "interRegion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "interRegion",sample.Data(), syst.Data(),   20,    0,     1, flavtag, systType);



  addHisto("mWT",       "QCDnormLregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "QCDnormWregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "QCDnormBregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "QCDnormSregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "QCDnormTTregion", 	sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);

  addHisto("mWT",       "qcdLregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("MET",       "qcdLregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "qcdWregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("MET",       "qcdWregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "qcdBregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("MET",       "qcdBregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "qcdSregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("MET",       "qcdSregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("mWT",       "qcdTTregion", 	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("MET",       "qcdTTregion", 	    sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);

  addHisto("ptW",           "Signalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWT",           "Signalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("MET",           "Signalregion",  	sample.Data(), syst.Data(),   20,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "Signalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

  addHisto("NVtx",          "Selectedsignalregion",     sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",          "Selectedsignalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",         "Selectedsignalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",         "Selectedsignalregion",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",        "Selectedsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",        "Selectedsignalregion",  	sample.Data(), syst.Data(),   30,   0,    300, flavtag, systType);
  addHisto("LeptEta",       "Selectedsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "Selectedsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",          "Selectedsignalregion", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   20,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiWptMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "Selectedsignalregion",  sample.Data(), syst.Data(),   20,    0,     1, flavtag, systType);

  addHisto("JetPt",         "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",        "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "ATLASFCNCsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWT",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",      "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("MET",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

  addHisto("JetPt",         "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",        "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "ATLASRESsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWT",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",      "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("MET",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

  addHisto("NVtx",      "ttbarregion_2j2b",     sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",      "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("NBJet",     "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("mWT",       "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",  "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",  "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "ttbarregion_2j2b", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "ttbarregion_2j2b", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeptEta",   "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","ttbarregion_2j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiWptMet","ttbarregion_2j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "ttbarregion_2j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("JetPt",         "ttbarregion_2j2b",     sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeadJetPt",     "ttbarregion_2j2b",     sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("SubLeadJetPt",  "ttbarregion_2j2b",     sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",        "ttbarregion_2j2b",     sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeadJetEta",    "ttbarregion_2j2b",     sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("SubLeadJetEta", "ttbarregion_2j2b",     sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "ttbarregion_2j2b", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("NVtx",      "ttbarregion_3j2b",     sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",      "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("NBJet",     "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("mWT",       "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",  "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",  "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "ttbarregion_3j2b", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "ttbarregion_3j2b", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeptEta",   "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","ttbarregion_3j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "ttbarregion_3j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("JetPt",     "ttbarregion_3j2b",     sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",    "ttbarregion_3j2b",     sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "ttbarregion_3j2b", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("NVtx",      "ttbarregion_4j2b",     sample.Data(), syst.Data(),   60,    0,    60, flavtag, systType);
  addHisto("NJet",      "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("NBJet",     "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("mWT",       "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_high",  "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   28,   40,   600, flavtag, systType);
  addHisto("mWT_full",  "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "ttbarregion_4j2b", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "ttbarregion_4j2b", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("LeptEta",   "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","ttbarregion_4j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "ttbarregion_4j2b",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("JetPt",     "ttbarregion_4j2b",     sample.Data(), syst.Data(),   30,    0,   300, flavtag, systType);
  addHisto("JetEta",    "ttbarregion_4j2b",     sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "ttbarregion_4j2b", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

//  cout << "      Histograms properly initialized!   " << endl;
//  cout << endl;

  // Prepare a TTree for the BdT
  TString treename = "Ttree_"+sample;
  TTree * TheTree = new TTree(treename.Data(),treename.Data());
  TheTree->Branch("tree_EvtWeight" ,&tree_EvtWeight  ,"tree_EvtWeight/F"    );
  TheTree->Branch("tree_wMass"     ,&tree_wMass      ,"tree_wMass/F"        );
  TheTree->Branch("tree_wPt"       ,&tree_wPt        ,"tree_wPt/F"          );
  TheTree->Branch("tree_met"       ,&tree_met        ,"tree_met/F"          );
  TheTree->Branch("tree_deltaPhilb",&tree_deltaPhilb ,"tree_deltaPhilb/F"   );

  theTree_list.push_back(TheTree);
  theTree_map[sample] = theTree_list.back();

  tree_EvtWeight  = -10000;
  tree_wMass      = -10000;
  tree_wPt        = -10000;
  tree_met        = -10000;
  tree_deltaPhilb = -10000;


}


//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1D binning
//creates one histograms per channel
//-------------------------------------------------------------

void TreeReader::addHisto(TString var, TString selstep, TString sample, TString syst, int nbins, float min, float max, TString flavtag, short int systType)
{


  TString name;

  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(     syst == "" && systType == -1) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_matchingdown";
      else if(syst == "" && systType ==  1) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_matchingup";
      else if(syst == "" && systType == -2) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_scaledown";
      else if(syst == "" && systType ==  2) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_scaleup";
      else if(syst == "")                   name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag;
      else                                  name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"__"+syst;
  }
  else                                      name =  var+"_"+thechannel+"_"+selstep+"__"+sample;


  TH1D * thehisto = new TH1D(name,name,nbins,min,max);
  thehisto->Sumw2();
  histo_list.push_back(thehisto);
  histo_map[name.Data()] = numb_histo;

  numb_histo++;


}

//-------------------------------------------------------------
//fill histograms
//first parameter is the channel,
//second parameter is the variable name,
//third parameter is the selection step (like "afterleptsel")
//forth parameter is the sample name (like "Z)
//others are value and weight
//-------------------------------------------------------------
void TreeReader::fillHisto(TString var, TString selstep, TString sample, TString syst, float val, float weight, TString flavtag, short int systType)
{

  TString name;

  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(     syst == "" && systType == -1) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_matchingdown";
      else if(syst == "" && systType == 1)  name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_matchingup";
      else if(syst == "" && systType == -2) name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_scaledown";
      else if(syst == "" && systType == 2)  name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"_scaleup";
      else if(syst == "")                   name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag;
      else                                  name =  var+"_"+thechannel+"_"+selstep+"__WExcl"+flavtag+"__"+syst;
  }
  else                                      name =  var+"_"+thechannel+"_"+selstep+"__"+sample;


  if(histo_map[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histogram " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }

  if(thechannel == "mujets" || thechannel == "eljets")     histo_list[histo_map[name.Data()]]->Fill(val, weight);

}

void TreeReader::manageOverflows()
{

  if(histo_list.size() == 0)
  {
    cout << "   WARNING trying to manage a non existing histogram at line: " << __LINE__ << endl;
  }

  for ( unsigned short int hist = 0; hist < histo_list.size(); hist++)
  {
      double lastBinContent  = histo_list[hist]->GetBinContent(histo_list[hist]->GetNbinsX()    );
      double overflowContent = histo_list[hist]->GetBinContent(histo_list[hist]->GetNbinsX() +1 );
      double lastBinError    = histo_list[hist]->GetBinError(  histo_list[hist]->GetNbinsX()    );
      double overflowError   = histo_list[hist]->GetBinError(  histo_list[hist]->GetNbinsX() +1 );
      double newLastBinError = sqrt(lastBinError*lastBinError + overflowError*overflowError     );
      histo_list[hist]->SetBinContent(histo_list[hist]->GetNbinsX()   , lastBinContent + overflowContent);
      histo_list[hist]->SetBinError(  histo_list[hist]->GetNbinsX()   , newLastBinError                 );
      histo_list[hist]->SetBinContent(histo_list[hist]->GetNbinsX()+1 , 0                               );
      histo_list[hist]->SetBinError(  histo_list[hist]->GetNbinsX()+1 , 0                               );
  }
}



void TreeReader::deleteHisto()
{
}

//SetUp CSV reweighting
void TreeReader::SetUpCSVreweighting()
{

  // Do not set it up if we're running on collision data
  if(isData){ return; }
  f_CSVwgt_HF = new TFile ("../BTagCSV/CSVRW/csv_rwt_hf.root");
  f_CSVwgt_LF = new TFile ("../BTagCSV/CSVRW/csv_rwt_lf.root");


  // CSV reweighting
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c  = "final";
    TString syst_csv_suffix_lf = "final";

    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp";
      syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown";
      syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp";
      syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown";
      syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up";
      syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down";
      syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up";
      syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down";
      syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<6; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<6; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }

    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)f_CSVwgt_LF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

}

// Get event weight
double TreeReader::GetCSVweight(const int iSys, int jet_n, float *jet_pt,float *jet_eta,float *jet_btagdiscri,int *jet_flav)
{
  if (isData) return 1.0;

  int iSysHF = 0;
  int iSysC = 0;
  int iSysLF = 0;

  // systematic variations
  // 1 HF/LF JESup iSysHF=1 iSysC=0 iSysLF=1
  // 2 HF/LF JESdown iSysHF=2 iSysC=0 iSysLF=2
  // 3 HF CSVLFup iSysHF=3 iSysC=0 iSysLF=0
  // 4 HF CSVLFdown iSysHF=4 iSysC=0 iSysLF=0
  // 5 HF CSVHFStats1up iSysHF=5 iSysC=0 iSysLF=0
  // 6 HF CSVHFStats1down iSysHF=6 iSysC=0 iSysLF=0
  // 7 HF CSVHFStats2up iSysHF=7 iSysC=0 iSysLF=0
  // 8 HF CSVHFStats2down iSysHF=8 iSysC=0 iSysLF=0
  // 9 C CSVCErr1up iSysHF=0 iSysC=1 iSysLF=0
  // 10 C CSVCErr1down iSysHF=0 iSysC=2 iSysLF=0
  // 11 C CSVCErr2up iSysHF=0 iSysC=3 iSysLF=0
  // 12 C CSVCErr2down iSysHF=0 iSysC=4 iSysLF=0
  // 13 LF CSVHFup iSysHF=0 iSysC=0 iSysLF=3
  // 14 LF CSVHFdown iSysHF=0 iSysC=0 iSysLF=4
  // 15 LF CSVLFStats1up iSysHF=0 iSysC=0 iSysLF=5
  // 16 LF CSVLFStats1down iSysHF=0 iSysC=0 iSysLF=6
  // 17 LF CSVLFStats2up iSysHF=0 iSysC=0 iSysLF=7
  // 18 LF CSVLFStats2down iSysHF=0 iSysC=0 iSysLF=8
  if( iSys == 1 ) {iSysHF=1; iSysC=0; iSysLF=1;}
  else if( iSys == 2 ) {iSysHF=2; iSysC=0; iSysLF=2;}
  else if( iSys == 3 ) {iSysHF=3; iSysC=0; iSysLF=0;}
  else if( iSys == 4 ) {iSysHF=4; iSysC=0; iSysLF=0;}
  else if( iSys == 5 ) {iSysHF=5; iSysC=0; iSysLF=0;}
  else if( iSys == 6 ) {iSysHF=6; iSysC=0; iSysLF=0;}
  else if( iSys == 7 ) {iSysHF=7; iSysC=0; iSysLF=0;}
  else if( iSys == 8 ) {iSysHF=8; iSysC=0; iSysLF=0;}
  else if( iSys == 9 ) {iSysHF=0; iSysC=1; iSysLF=0;}
  else if( iSys == 10 ) {iSysHF=0; iSysC=2; iSysLF=0;}
  else if( iSys == 11 ) {iSysHF=0; iSysC=3; iSysLF=0;}
  else if( iSys == 12 ) {iSysHF=0; iSysC=4; iSysLF=0;}
  else if( iSys == 13 ) {iSysHF=0; iSysC=0; iSysLF=3;}
  else if( iSys == 14 ) {iSysHF=0; iSysC=0; iSysLF=4;}
  else if( iSys == 15 ) {iSysHF=0; iSysC=0; iSysLF=5;}
  else if( iSys == 16 ) {iSysHF=0; iSysC=0; iSysLF=6;}
  else if( iSys == 17 ) {iSysHF=0; iSysC=0; iSysLF=7;}
  else if( iSys == 18 ) {iSysHF=0; iSysC=0; iSysLF=8;}

//  iSysHF=1; // JESup
//  iSysHF=2; // JESdown
//  iSysHF=3; // CSVLFup
//  iSysHF=4; // CSVLFdown
//  iSysHF=5; // CSVHFStats1up
//  iSysHF=6; // CSVHFStats1down
//  iSysHF=7; // CSVHFStats2up
//  iSysHF=8; // CSVHFStats2down

//  iSysC=1; // CSVCErr1up
//  iSysC=2; // CSVCErr1down
//  iSysC=3; // CSVCErr2up
//  iSysC=4; // CSVCErr2down

//  iSysLF=1; // JESup
//  iSysLF=2; // JESdown
//  iSysLF=3; // CSVHFup
//  iSysLF=4; // CSVHFdown
//  iSysLF=5; // CSVLFStats1up
//  iSysLF=6; // CSVLFStats1down
//  iSysLF=7; // CSVLFStats2up
//  iSysLF=8; // CSVLFStats2down


  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int ij=0;ij<jet_n;ij++){

    double csv = jet_btagdiscri[ij];
    double jetPt = jet_pt[ij];
    double jetAbsEta = fabs( jet_eta[ij] );
    int flavor = abs( jet_flav[ij] );

    if( jetPt < 20. || jetAbsEta > 2.5 ) continue;

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100 && jetPt<160) iPt = 4;
    else if (jetPt >=160 && jetPt<10000) iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;
    // kskovpen hack - should change eta cut to 2.4 !
    else if ( jetAbsEta>=2.41 && jetAbsEta<=2.5) iEta = 2;

    if (iPt < 0 || iEta < 0) cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;

    }
    else {
      if (iPt >=3) iPt=3;       /// [>60]
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;


    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  return csvWgtTotal;
}


double TreeReader::getQCDscalefactor(vector<TString> datalist, vector<TString> mclist, TString normregion)
{
  double SF = 0;

    //--------------------
    //load QCDdata mTW file
    //--------------------
    double yieldQCD = 0;
    TString fileToGetQCDData_name;
    TFile* fileToGetQCDData = 0;
    TString histo_qcddatadriven_name;
    TH1D* histo_qcddatadriven = 0;
    if(!useElectronChannel) fileToGetQCDData_name     = "outputroot_QCDcorr_iso0p5/backup_final_muons/histofile_QCDdatadriven.root";
    else                    fileToGetQCDData_name     = "outputroot_QCDcorr_iso0p5/backup_final_electrons/histofile_QCDdatadriven.root";
    fileToGetQCDData          = new TFile(fileToGetQCDData_name);
    histo_qcddatadriven_name  = normregion+"__QCDdatadriven";
    histo_qcddatadriven       = (TH1D*)fileToGetQCDData->Get(histo_qcddatadriven_name);

    if (histo_qcddatadriven != 0)
    {
        yieldQCD += histo_qcddatadriven->Integral();
    }
    else cout << "histo_qcddatadriven: " << histo_qcddatadriven_name << " not found! " << endl;
    fileToGetQCDData->Close();
    //--------------------
    //loop over Datasamples
    //--------------------
    double yieldData = 0;
    TString fileToGetData_name;
    TFile* fileToGetData = 0;
    TString histo_data_name;
    TH1D* histo_data = 0;
    for(unsigned int dataSample = 0; dataSample < datalist.size() ; dataSample++)
    {
        if(!useElectronChannel) fileToGetData_name  = "outputroot_woQCDcorr/backup_final_muons/histofile_"+datalist[dataSample]+".root";
        else                    fileToGetData_name  = "outputroot_woQCDcorr/backup_final_electrons/histofile_"+datalist[dataSample]+".root";
        fileToGetData       = new TFile(fileToGetData_name);
        histo_data_name     = normregion+"__"+datalist[dataSample];
        histo_data          = (TH1D*)fileToGetData->Get(histo_data_name);

        if (histo_data != 0)
        {
            yieldData += histo_data->Integral();
        }
        else cout << "histo_data: " << histo_data_name << " not found! " << endl;
        fileToGetData->Close();
    }

    //--------------------
    //loop over MCsamples
    //--------------------
    double yieldMC = 0;
    TString fileToGetMC_name;
    TFile* fileToGetMC = 0;
    TString histo_mc_name;
    TH1D* histo_mc = 0;
    for(unsigned int mcSample = 0; mcSample < mclist.size() ; mcSample++)
    {
        if(!useElectronChannel) fileToGetMC_name    = "outputroot_woQCDcorr/backup_final_muons/histofile_"+mclist[mcSample]+".root";
        else                    fileToGetMC_name    = "outputroot_woQCDcorr/backup_final_electrons/histofile_"+mclist[mcSample]+".root";
        fileToGetMC         = new TFile(fileToGetMC_name);

        if( mclist[mcSample] == "W0Jets" ||  mclist[mcSample] == "W1Jets" ||  mclist[mcSample] == "W2Jets" ||  mclist[mcSample] == "W3Jets" ||  mclist[mcSample] == "W4Jets" )
        {
             histo_mc    =     (TH1D*)fileToGetMC->Get( (normregion+"__WExclb").Data() )->Clone();
             histo_mc->Add(    (TH1D*)fileToGetMC->Get( (normregion+"__WExclc").Data() )->Clone() );
             histo_mc->Add(    (TH1D*)fileToGetMC->Get( (normregion+"__WExcll").Data() )->Clone() );
        }
        else
        {
            histo_mc_name       = normregion+"__"+mclist[mcSample];
            histo_mc            = (TH1D*)fileToGetMC->Get(histo_mc_name)->Clone();
        }

        if (histo_mc != 0)
        {
    	    yieldMC += histo_mc->Integral();
        }
        else cout << "histo_mc: " << histo_mc_name << " not found! " << endl;

        if (fileToGetMC != 0) fileToGetMC->Close();
    }
    if (yieldQCD != 0)  SF = (yieldData - yieldMC)/yieldQCD;
    else 	            SF = 0;
    return SF;
}

vector<double> TreeReader::getY( TGraphAsymmErrors* graph, double Xvalue)
{
  double Yvalue = 1;
  double Yerror = 0;
  vector<double> Youtput;
  double* X = graph->GetX();
  double* Y = graph->GetY();

  for ( int i = 0; X[i] <= 500 ; i++)
  {
       if( i != 0 && X[i] > Xvalue && Xvalue <= 170)
       {
            // Y = a*X + b with b equal to Y[i] - a*X[i] ~interpolation between two points of the graph
            Yvalue = ((Y[i] - Y[i-1])/(X[i] - X[i-1]))* Xvalue + Y[i] - ((Y[i] - Y[i-1])/(X[i] - X[i-1]))*X[i];
            Yerror = (graph->GetErrorY(i) + graph->GetErrorY(i-1))/2;
            break;
        }
       else if( Xvalue > 170)
       {
            Yvalue = Y[7];
            Yvalue = Y[7];
            Yerror = graph->GetErrorY(7);
            break;
       }
  }

  Youtput.push_back(Yvalue);
  Youtput.push_back(Yerror);
  return Youtput;
}

vector<double> TreeReader::getSFtrigger( TGraphAsymmErrors* ratioPlot,  double pT, double eta)
{
   vector<double> SF;
   double SFvalue = 1;
   double SFerror = 0;

   if (abs(eta) <= 2.1 && pT <= 500)
   {
       SFvalue = getY( ratioPlot, pT )[0];
       SFerror = getY( ratioPlot, pT )[1];
   }
   else
   {
       SFvalue = getY( ratioPlot, 499)[0];
       SFerror = getY( ratioPlot, 499)[1];
   }

   SF.push_back(SFvalue);
   SF.push_back(SFerror);

   return SF;
}

void TreeReader::loadPDFSumWeight(TString channel, TString tmp_sample, TString flavtag)
{

    TString sample = tmp_sample;
    if      ( isW && flavtag == "b") sample+="_bflav";
    else if ( isW && flavtag == "c") sample+="_cflav";
    else if ( isW && flavtag == "l") sample+="_lflav";

    TFile *theinput_file = new TFile(("InputPDF_Weights/histofile_"+sample+".root").Data());

    if      ( isW && flavtag == "b") sample = "WExclb";
    else if ( isW && flavtag == "c") sample = "WExclc";
    else if ( isW && flavtag == "l") sample = "WExcll";

    TString histoname_pdfUp    = "pdfWeight_"+channel+"___"+sample+"__PDF__plus";
    TH1D* histo_pdfUp          = (TH1D*)theinput_file->Get( histoname_pdfUp.Data() )->Clone();

    TString histoname_pdfDown  = "pdfWeight_"+channel+"___"+sample+"__PDF__minus";
    TH1D* histo_pdfDown        = (TH1D*)theinput_file->Get( histoname_pdfDown.Data() )->Clone();


    TString histoname_sum  = "sumWeight_"+channel+"___"+sample;
    TH1D* histo_sum        = (TH1D*)theinput_file->Get( histoname_sum.Data() )->Clone();


    //PDFSumWeight[histoname_pdfUp]   = histo_sum->GetBinContent(1);
    //PDFSumWeight[histoname_pdfDown] = histo_sum->GetBinContent(1);
    PDFSumWeight[histoname_pdfUp]   = histo_sum->GetBinContent(1)/histo_pdfUp->GetBinContent(1);
    PDFSumWeight[histoname_pdfDown] = histo_sum->GetBinContent(1)/histo_pdfDown->GetBinContent(1);
    cout << "PDF_Sum = " << histo_sum->GetBinContent(1) << " | PDF_down = " << histo_pdfDown->GetBinContent(1) << endl;

    theinput_file->Close();
    delete theinput_file;

}

