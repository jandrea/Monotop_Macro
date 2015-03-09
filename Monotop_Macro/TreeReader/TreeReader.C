#define TreeReader_cxx
#include "TreeReader.h"
#include "common.h"
#include <string>
#include <TGraphAsymmErrors.h>

using namespace std;


void TreeReader::Loop(short int CorrOption, std::vector<TString> datalist, std::vector<TString> datalist_longnames, std::vector<TString> mclist, TString sample, std::vector<TString> thesystlist, TString flavtag, short int QCDsyst)
{
   TString thesample  = sample;
   TFile * theoutputfile = 0;
   isQCD    = (sample == "QCD_A"  || sample == "QCD_B"  || sample == "QCD_C"  || sample == "QCD_D" );
   isW      = (sample == "WJets"  || sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets");
   isWExcl  = (sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets");
   isWIncl  = (sample == "WJets");
   TString sampleroot;
   if      ( isW && flavtag == "b") sampleroot = thesample+"_bflav";
   else if ( isW && flavtag == "c") sampleroot = thesample+"_cflav";
   else if ( isW && flavtag == "l") sampleroot = thesample+"_lflav";
   else                             sampleroot = thesample;

   double qcdIsoCut = 0.5;

   //if      (CorrOption == 0)  theoutputfile = new TFile( ("outputroot_QCDcorr/histofile_"     +sampleroot+".root" ).Data() , "recreate");
   if      (CorrOption == 0)  theoutputfile = new TFile( ("outputroot_QCDcorr_iso0p5/histofile_"     +sampleroot+".root" ).Data() , "recreate");
   //if      (CorrOption == 0)  theoutputfile = new TFile( ("outputroot_QCDcorr_iso0p6/histofile_"     +sampleroot+".root" ).Data() , "recreate");
   else if (CorrOption == 1)  theoutputfile = new TFile( ("outputroot_woQCDcorr/histofile_"   +sampleroot+".root" ).Data() , "recreate");
   else if (CorrOption == 2)  theoutputfile = new TFile( ("outputroot_withQCDcorr/histofile_" +sampleroot+".root" ).Data() , "recreate");
   else if (CorrOption == 3)
   {
       if (QCDsyst == -1)       { theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__minus.root" ).Data() , "recreate"); qcdIsoCut = 0.6; }
       else if (QCDsyst == 1)   { theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+"__plus.root"  ).Data() , "recreate"); qcdIsoCut = 0.4; }
       else                       theoutputfile = new TFile( ("outputroot_withSyst/histofile_"    +sampleroot+".root"        ).Data() , "recreate");
   }
   else cout << "ERROR: Wrong value of CorrOption or flavtag! Allowed values: 0,1,2,3 for corrOption and Incl, b, c, l for flavtag (if CorrOption == 3)" << endl;

   for(unsigned int i=0; i< systlist.size(); i++){
     TString samplename = "";
     if(      systlist[i] == "" && isQCD && QCDsyst == -1) samplename = thesample+"__minus";
     else if( systlist[i] == "" && isQCD && QCDsyst == 1 ) samplename = thesample+"__plus";
     else if( systlist[i] == "")                           samplename = thesample;
     else                                                  samplename = thesample+"__"+systlist[i];

     bool     firstinit = false;
     if(i==0) firstinit = true;
     initializeHisto(samplename, systlist[i], firstinit, flavtag);
   }

   isData = (sample == "SingleMuA" || sample == "SingleMuB" || sample == "SingleMuC" || sample == "SingleMuD" || sample == "NTuple_53X_SingleMuRun2012A" || sample == "NTuple_53X_SingleMuRun2012B" || sample == "NTuple_53X_SingleMuRun2012C" || sample == "NTuple_53X_SingleMuRun2012D" );

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

   double SF_QCD_W = 0, SF_QCD_B = 0, SF_QCD_S = 0, SF_QCD_TT = 0;
   if (CorrOption == 2 || CorrOption == 3)
   {
       SF_QCD_W  = getQCDscalefactor(datalist, datalist_longnames, mclist, "mWT_mujets_QCDnormWregion" , flavtag );
       SF_QCD_B  = getQCDscalefactor(datalist, datalist_longnames, mclist, "mWT_mujets_QCDnormBregion" , flavtag );
       SF_QCD_S  = getQCDscalefactor(datalist, datalist_longnames, mclist, "mWT_mujets_QCDnormSregion" , flavtag );
       SF_QCD_TT = getQCDscalefactor(datalist, datalist_longnames, mclist, "mWT_mujets_QCDnormTTregion", flavtag );
   }

   TGraphAsymmErrors* ratio_eta_0to0p9;
   TGraphAsymmErrors* ratio_eta_0p9to1p2;
   TGraphAsymmErrors* ratio_eta_1p2to2p1;

   TFile* ftrigger = new TFile(string("SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root").c_str());
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Barrel_0to0p9_pt25-500_2012ABCD",ratio_eta_0to0p9);
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Transition_0p9to1p2_pt25-500_2012ABCD",ratio_eta_0p9to1p2);
   ftrigger->GetObject("IsoMu24_eta2p1_DATA_over_MC_TightID_IsodB_PT_ABSETA_Endcaps_1p2to2p1_pt25-500_2012ABCD",ratio_eta_1p2to2p1);
   ftrigger->Close();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if (nentries > 10000 && (jentry % (nentries / 10000) == 0)) 			            printProgressBar(jentry,nentries,sample);
      if (nentries > 100   && nentries <= 10000 && (jentry % (nentries / 100) == 0)) 	printProgressBar(jentry,nentries,sample);

      TString thechannel = "mujets";

      bool hasBjet  = false, hasCjet  = false, isBflavour = false, isCflavour = false, isLflavour = false, match = false;
      if(CorrOption != 0 && isW && (flavtag != "Incl" && flavtag != "allflav") )
      {
         int     iter_jets   = smalltree_njets;
         float * jet_pt	     = smalltree_jet_pt;
         float * jet_eta	 = smalltree_jet_eta;
         int   * jet_flav    = smalltree_jet_flav;

         for(short int ijet=0; ijet<iter_jets; ijet++)
          {
              if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4)  continue;
              if(abs(jet_flav[ijet]) == 5)                        hasBjet = true;
              if(abs(jet_flav[ijet]) == 4)                        hasCjet = true;
          }

          //---------------
          //W flavour study
          //---------------
          if      ( hasBjet && !hasCjet) isBflavour = true;
          else if (!hasBjet &&  hasCjet) isCflavour = true;
          else if (!hasBjet && !hasCjet) isLflavour = true;

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


      double SFtrigger      = 1;
      double SFtriggerError = 0;

      if (abs(smalltree_lept_eta[0]) < 0.9)
      {
          SFtrigger      = getSFtrigger(ratio_eta_0to0p9,  smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_0to0p9,  smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
      }
      else if (abs(smalltree_lept_eta[0]) > 0.9 && abs(smalltree_lept_eta[0]) < 1.2)
      {
          SFtrigger      = getSFtrigger(ratio_eta_0p9to1p2,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_0p9to1p2,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
      }
      else if (abs(smalltree_lept_eta[0]) > 1.2 && abs(smalltree_lept_eta[0]) < 2.1)
      {
          SFtrigger      = getSFtrigger(ratio_eta_1p2to2p1,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[0];
          SFtriggerError = getSFtrigger(ratio_eta_1p2to2p1,smalltree_lept_pt[0] ,smalltree_lept_eta[0] )[1];
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
      if      (CorrOption == 0)         applyEventSel(CorrOption, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0],     thesample, flavtag);
      else if (CorrOption == 1)         applyEventSel(CorrOption, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0],     thesample, flavtag);
      else if (CorrOption == 2)         applyEventSel(CorrOption, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0],     thesample, flavtag);
      else if (CorrOption == 3)
      {
          for(unsigned int isyst=0; isyst< systlist.size(); isyst++)
          {
              //if( (systlist[isyst] == "toppt__plus" || systlist[isyst] == "toppt__minus") && thesample != "TTbar_Madgraph" ) continue;
              applyEventSel(CorrOption, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[isyst], thesample, flavtag);
          }
      }
      else cout << "ERROR: Wrong value of CorrOption! Allowed valued: 0,1,2,3" << endl;


   }

/*
   if ( CorrOption == 3 && isW )
   {
       TString sampleTag;
       if (isWExcl) sampleTag = "WExcl";
       else         sampleTag = "WIncl";
       TH1F* histoWCorrWeights = getWCorrWeights(datalist, mclist, "mWT_mujets_Wregion_highpt", flavtag );
       TString nosyst;
       if (flavtag == "b" || flavtag == "c" || flavtag == "l") nosyst = sampleTag+flavtag;
       else                                                    nosyst = sampleTag;

       scaleHisto("mujets", nosyst.Data(),     histoWCorrWeights);
   }
*/

   cout << endl;
   theoutputfile->Write();
   deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;

}



bool TreeReader::applyEventSel(short int CorrOption, double SF_QCD_W, double SF_QCD_B, double SF_QCD_S, double SF_QCD_TT, double SFtrigger, double SFtriggerError, TString thechannel, TString systtype, TString sample, TString flavtag){

      bool applyCSV_reshape = 0;
      bool applyCSV         = 1;

      int  btagSys          = 1;
      int  mistagSys        = 1;

      double met_pt    = smalltree_met_pt;
      double met_phi   = smalltree_met_phi;
      double evtweight = -1;
      TString thesample = sample;
      if(systtype != "") thesample = thesample + "__"+systtype;

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
         systtype == "PDF__plus"   ||   systtype == "PDF__minus"   ||
     ( ( systtype == "W__plus"     ||   systtype == "W__minus"  )  && isW )

      ){
	    iter_jets      = smalltree_njets;
        jet_pt         = smalltree_jet_pt;
        jet_eta        = smalltree_jet_eta;
        jet_phi        = smalltree_jet_phi;
        jet_btagdiscri = smalltree_jet_btagdiscri;
        jet_flav       = smalltree_jet_flav;

		evtweight = smalltree_evtweight;

		if(     systtype == "lept__plus")    evtweight = smalltree_weight_leptup;
		else if(systtype == "lept__minus")   evtweight = smalltree_weight_leptdown;
		else if(systtype == "trig__plus")    evtweight = (SFtrigger + SFtriggerError)*smalltree_evtweight/SFtrigger;
		else if(systtype == "trig__minus")   evtweight = (SFtrigger - SFtriggerError)*smalltree_evtweight/SFtrigger;
		else if(systtype == "PU__plus")      evtweight = smalltree_weight_PUup;
		else if(systtype == "PU__minus")     evtweight = smalltree_weight_PUdown;
		else if(systtype == "PDF__plus")     evtweight = smalltree_weight_PDFup;
		else if(systtype == "PDF__minus")    evtweight = smalltree_weight_PDFdown;
		else if(systtype == "toppt__plus")   evtweight = smalltree_weight_toppt*smalltree_evtweight;//FIXME check the toppt systs
		else if(systtype == "toppt__minus")  evtweight = smalltree_evtweight;
		else if(systtype == "W__plus")       evtweight = smalltree_evtweight;
		else if(systtype == "W__minus")      evtweight = smalltree_evtweight;

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

      }else if(systtype == "btag__JES__plus"){
	    wCSV           = GetCSVweight(1,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__JES__minus"){
	    wCSV           = GetCSVweight(2,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLF__plus"){
	    wCSV           = GetCSVweight(3,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLF__minus"){
	    wCSV           = GetCSVweight(4,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHFStats1__plus"){
	    wCSV           = GetCSVweight(5,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHFStats1__minus"){
	    wCSV           = GetCSVweight(6,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHFStats2__plus"){
	    wCSV           = GetCSVweight(7,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHFStats2__minus"){
	    wCSV           = GetCSVweight(8,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVCErr1__plus"){
	    wCSV           = GetCSVweight(9,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVCErr1__minus"){
	    wCSV           = GetCSVweight(10,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVCErr2__plus"){
	    wCSV           = GetCSVweight(11,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVCErr2__minus"){
	    wCSV           = GetCSVweight(12,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHF__plus"){
	    wCSV           = GetCSVweight(13,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVHF__minus"){
	    wCSV           = GetCSVweight(14,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLFStats1__plus"){
	    wCSV           = GetCSVweight(15,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLFStats1__minus"){
	    wCSV           = GetCSVweight(16,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLFStats2__plus"){
	    wCSV           = GetCSVweight(17,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__CSVLFStats2__minus"){
	    wCSV           = GetCSVweight(18,smalltree_njets,smalltree_jet_pt,smalltree_jet_eta,smalltree_jet_btagdiscri,smalltree_jet_flav);

      }else if(systtype == "btag__plus"   ){btagSys =  1;}
      else if(systtype  == "btag__minus"  ){btagSys = -1;}
      else if(systtype  == "mistag__plus" ){btagSys =  1;}
      else if(systtype  == "mistag__minus"){btagSys = -1;}
      else if(systtype  == "W__minus" || systtype == "W__plus") ;
      else{
        cout << "WARNING syst type " << systtype << " not recognized !! " << endl;
	    cout << "correct syst types are " << endl;
	    cout << " \"\",  \"lept__plus\", \"lept__minus\", \"trig__plus\", \"trig__minus\", \"PDF__plus\", \"PDF__minus\"  " << endl;
	    cout << "\"jes__plus\", \"jes__minus\", \"jer__plus\", \"jer__minus\", \"metuncls__plus\", \"metuncls__minus\" " << endl;
      }

      if( applyCSV_reshape && !applyCSV ) evtweight *= wCSV;

      //determine the btag systematic
      if(systtype == "btag__plus" ){
        btagSys   = 1;
        evtweight = smalltree_evtweight;
      }
      if(systtype == "btag__minus" ){
        btagSys   = -1;
        evtweight = smalltree_evtweight;
      }
      if(systtype == "mistag__plus" ){
        mistagSys = 1;
        evtweight = smalltree_evtweight;
      }
      if(systtype == "mistag__minus"){
        mistagSys = -1;
        evtweight = smalltree_evtweight;
      }

     // Add the trigger efficiency
     if(sample == "SingleMuA" || sample == "SingleMuB" || sample == "SingleMuC" || sample == "SingleMuD" || sample == "NTuple_53X_SingleMuRun2012A" || sample == "NTuple_53X_SingleMuRun2012B" || sample == "NTuple_53X_SingleMuRun2012C" || sample == "NTuple_53X_SingleMuRun2012D" ) evtweight = 1;
     else evtweight *= SFtrigger;

     if ((CorrOption == 1 || CorrOption == 2 || CorrOption == 3) && thesample == "Wminus_Pohweg")   evtweight *= 5074.7/6200;
     if ((CorrOption == 1 || CorrOption == 2 || CorrOption == 3) && thesample == "Wplus_Pohweg")    evtweight *= 7213.4/6200;
     if  (CorrOption == 2 || CorrOption == 3 )
     {
         if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")       evtweight = 1;
     }

     if (!isfinite(evtweight)) return false;

     if(smalltree_lept_pt[0] > 33 && fabs(smalltree_lept_eta[0]) < 2.1 ){

       TLorentzVector lept, met, leadingJet;
       lept.SetPtEtaPhiM(smalltree_lept_pt[0], smalltree_lept_eta[0],  smalltree_lept_phi[0], 0.);
       met.SetPtEtaPhiM(met_pt, 0, met_phi , 0.);

       if(iter_jets>0) leadingJet.SetPtEtaPhiM(jet_pt[0], jet_eta[0], jet_phi[0], 0);

       int njets=0;
       int nbjets = 0;
       bool jetsup70 = false;
       std::vector<int > btagged_jet_idx;
       for(int ijet=0; ijet<iter_jets; ijet++){
         if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
         njets++;
         fillHisto(thechannel, "JetPt",     "afterleptsel",  thesample, systtype,  jet_pt[ijet] ,  evtweight, flavtag);
         fillHisto(thechannel, "JetEta",    "afterleptsel",  thesample, systtype,  jet_eta[ijet] , evtweight, flavtag);
         //if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.244)  {nbjets++;   btagged_jet_idx.push_back(ijet); }
         if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.679)  {nbjets++;   btagged_jet_idx.push_back(ijet); }
         //if(applyCSV_reshape && !applyCSV && jet_btagdiscri[ijet] > 0.898)  {nbjets++;   btagged_jet_idx.push_back(ijet); }

         if( applyCSV){
           bool isbtag = 0;
           if(isData) isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], -999999 , jet_pt[ijet] , jet_eta[ijet], btagSys); // for nominal sample
           else{
               if(abs(jet_flav[ijet]) == 5 || abs(jet_flav[ijet])== 4){
                   isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], jet_flav[ijet] , jet_pt[ijet] , jet_eta[ijet], btagSys);
               }else{
                   isbtag = fBTagSF->IsTagged( jet_btagdiscri[ijet], jet_flav[ijet] , jet_pt[ijet] , jet_eta[ijet], mistagSys);
               }
               fillHisto(thechannel, "BTagProba", "afterleptsel", thesample, systtype, isbtag , evtweight, flavtag);
           }
           if(applyCSV && isbtag) {nbjets++; btagged_jet_idx.push_back(ijet); }
         }

	     if(jet_pt[ijet] > 70 ) jetsup70 = true;
       }



       double mTW = pow( 2*smalltree_lept_pt[0]*met_pt*(1-cos(smalltree_lept_phi[0] -  met_phi)) ,0.5);

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

       fillHisto(thechannel, "NJet",       "afterleptsel",  thesample, systtype,   iter_jets, 		        evtweight, flavtag);
       fillHisto(thechannel, "NBJet",      "afterleptsel",  thesample, systtype,   nbjets, 		            evtweight, flavtag);


       fillHisto(thechannel, "mWT",        "afterleptsel",  thesample, systtype,   mTW,    		            evtweight, flavtag);
       fillHisto(thechannel, "mW",         "afterleptsel",  thesample, systtype,   theWcand.M(),  	        evtweight, flavtag);
       fillHisto(thechannel, "MET",        "afterleptsel",  thesample, systtype,   met_pt, 		            evtweight, flavtag);
       fillHisto(thechannel, "mWTplusMET", "afterleptsel",  thesample, systtype,   mWTplusMET, 		        evtweight, flavtag);
       fillHisto(thechannel, "LeptPt",     "afterleptsel",  thesample, systtype,   smalltree_lept_pt[0], 	evtweight, flavtag);
       fillHisto(thechannel, "LeptEta",    "afterleptsel",  thesample, systtype,   smalltree_lept_eta[0],   evtweight, flavtag);

        //***********************************
        //qcd region
        //***********************************
    if(CorrOption == 0){

        if(jetsup70 && njets == 1 && nbjets == 0){
          fillHisto(thechannel, "mWT",        "qcdWregion",  thesample, systtype,   mTW,    			    evtweight, flavtag);
          fillHisto(thechannel, "MET",        "qcdWregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag);
        }
        else if(njets == 1 && nbjets == 1){
          fillHisto(thechannel, "mWT",        "qcdBregion",  thesample, systtype,   mTW,    			    evtweight, flavtag);
          fillHisto(thechannel, "MET",        "qcdBregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag);
        }
        if(jetsup70 && njets == 1 && nbjets == 1){
          fillHisto(thechannel, "mWT",        "qcdSregion",  thesample, systtype,   mTW,    			    evtweight, flavtag);
          fillHisto(thechannel, "MET",        "qcdSregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag);
        }
        else if(jetsup70 && njets >= 4 && nbjets ==2 ){
          fillHisto(thechannel, "mWT",        "qcdTTregion",  thesample, systtype,   mTW,    			    evtweight, flavtag);
          fillHisto(thechannel, "MET",        "qcdTTregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag);
        }
	}

        //***********************************
        //QCD normalization region
        //***********************************
    if(CorrOption == 0 || CorrOption == 1){

        if(jetsup70 && njets == 1 && nbjets == 0 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormWregion",   thesample, systtype,   mTW,    			evtweight, flavtag);
        }
        else if(njets == 1 && nbjets == 1 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormBregion",   thesample, systtype,   mTW,    			evtweight, flavtag);
        }
        if(jetsup70 && njets == 1 && nbjets == 1 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormSregion",   thesample, systtype,   mTW,    			evtweight, flavtag);
        }
        else if(jetsup70 && njets >= 4 && nbjets == 2 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormTTregion",  thesample, systtype,   mTW,    			evtweight, flavtag);
        }
	}

    if(CorrOption == 1 || CorrOption == 2 || CorrOption == 3)
    {
        //***********************************
        //1bjet region
        //***********************************
        if(nbjets == 1 && mTW > 40 ){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")                      evtweight = SF_QCD_B;
          fillHisto(thechannel, "mWT",        "1bjetregion",  thesample, systtype,   mTW,    			            evtweight, flavtag);
          fillHisto(thechannel, "MET",        "1bjetregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag);
          fillHisto(thechannel, "DeltaPhiLJ", "1bjetregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),evtweight, flavtag);
          fillHisto(thechannel, "DeltaRLJ",   "1bjetregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag);
          fillHisto(thechannel, "LeptPt",     "1bjetregion",  thesample, systtype,   smalltree_lept_pt[0], 		    evtweight, flavtag);
          fillHisto(thechannel, "LeptEta",    "1bjetregion",  thesample, systtype,   smalltree_lept_eta[0], 	    evtweight, flavtag);
	    }

        //***********************************
        //W enriched region
        //***********************************
	    if(jetsup70 && njets == 1 && nbjets == 0){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")                          evtweight = SF_QCD_W;
          if(mTW > 40)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",  "Wregion_highpt",  thesample, systtype,  jet_pt[ijet] ,                 evtweight, flavtag);
                fillHisto(thechannel, "JetEta", "Wregion_highpt",  thesample, systtype,  jet_eta[ijet] ,                evtweight, flavtag);
            }
            fillHisto(thechannel, "ptW",        "Wregion_highpt",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag);
            fillHisto(thechannel, "mWT_full",   "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag);
            fillHisto(thechannel, "mWT",        "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag);
            fillHisto(thechannel, "MET",        "Wregion_highpt",  thesample, systtype,   met_pt, 			                evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ", "Wregion_highpt",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), 	evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLMet","Wregion_highpt", thesample, systtype,   abs(lept.DeltaPhi(met)),          evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",   "Wregion_highpt",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag);
            fillHisto(thechannel, "LeptPt",     "Wregion_highpt",  thesample, systtype,   smalltree_lept_pt[0], 		    evtweight, flavtag);
            fillHisto(thechannel, "LeptEta",    "Wregion_highpt",  thesample, systtype,   smalltree_lept_eta[0], 		    evtweight, flavtag);
          }
	    }

        //***********************************
        //signal enriched region
        //***********************************
        if(jetsup70 && njets == 1 && nbjets == 1 ){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")          evtweight = SF_QCD_S;
          fillHisto(thechannel, "mWT_full",     "signalregion",  thesample, systtype,   mTW,    		evtweight, flavtag);
          // Manage overflows in SR - put higher values in the last bin by hands
          if(mTW    >= 500) mTW    = 499.95;
          if(met_pt >= 500) met_pt = 499.95;
          if(mTW > 40)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",  "signalregion",  thesample, systtype,  jet_pt[ijet] ,             evtweight, flavtag);
                fillHisto(thechannel, "JetEta", "signalregion",  thesample, systtype,  jet_eta[ijet] ,            evtweight, flavtag);
            }
            fillHisto(thechannel, "ptW",        "signalregion",  thesample, systtype,   theWcand.Pt(),       	         evtweight, flavtag);
            fillHisto(thechannel, "mWT",        "signalregion",  thesample, systtype,   mTW,    			             evtweight, flavtag);
            fillHisto(thechannel, "MET",        "signalregion",  thesample, systtype,   met_pt, 			             evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ", "signalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),  evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",   "signalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	     evtweight, flavtag);
          }
        }


        //***********************************
        //FCNC signal enriched region
        //***********************************
        if(jetsup70 && njets == 1 && nbjets == 1 ){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")          evtweight = SF_QCD_S;
          fillHisto(thechannel, "mWT_full",     "Selectedsignalregion",  thesample, systtype,   mTW,    		evtweight, flavtag);
          // Manage overflows in SR - put higher values in the last bin by hands
          if(mTW    >= 500) mTW    = 499.95;
          if(met_pt >= 500) met_pt = 499.95;
          if(mTW > 40 && abs(lept.DeltaPhi(leadingJet)) < 1.7 && met_pt > 100 && theWcand.Pt() > 50)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",      "Selectedsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag);
                fillHisto(thechannel, "JetEta",     "Selectedsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag);
            }
            fillHisto(thechannel, "ptW",            "Selectedsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag);
            fillHisto(thechannel, "mWT",            "Selectedsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag);
            fillHisto(thechannel, "MET",            "Selectedsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ",     "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLMet",   "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiJMet",   "Selectedsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJMet",  "Selectedsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",       "Selectedsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag);
          }
        }

        //***********************************
        //ATLAS signal enriched region
        //***********************************
        if(njets == 1 && nbjets == 1 ){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")          evtweight = SF_QCD_S;
          // Manage overflows in SR - put higher values in the last bin by hands
          if(mTW    >= 500) mTW    = 499.95;
          if(met_pt >= 500) met_pt = 499.95;
          if(mWTplusMET > 60 && mTW > 250 && abs(lept.DeltaPhi(leadingJet)) < 1.4)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",      "ATLASFCNCsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag);
                fillHisto(thechannel, "JetEta",     "ATLASFCNCsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag);
            }
            fillHisto(thechannel, "ptW",            "ATLASFCNCsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag);
            fillHisto(thechannel, "mWT",            "ATLASFCNCsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag);
            fillHisto(thechannel, "MET",            "ATLASFCNCsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag);
            fillHisto(thechannel, "mWTplusMET",     "ATLASFCNCsignalregion",  thesample, systtype,   mWTplusMET, 		            evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ",     "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiJMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJMet",  "ATLASFCNCsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",       "ATLASFCNCsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag);
          }

          if(mWTplusMET > 60 && mTW > 210 && abs(lept.DeltaPhi(leadingJet)) < 1.2)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",      "ATLASRESsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag);
                fillHisto(thechannel, "JetEta",     "ATLASRESsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag);
            }
            fillHisto(thechannel, "ptW",            "ATLASRESsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag);
            fillHisto(thechannel, "mWT",            "ATLASRESsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag);
            fillHisto(thechannel, "MET",            "ATLASRESsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag);
            fillHisto(thechannel, "mWTplusMET",     "ATLASRESsignalregion",  thesample, systtype,   mWTplusMET, 		            evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ",     "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiJMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJMet",  "ATLASRESsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",       "ATLASRESsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag);
          }

        }


        //***********************************
        //ttbar enriched region
        //***********************************
        if(jetsup70 && njets >= 4 && nbjets ==2){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")                          evtweight = SF_QCD_TT;
          fillHisto(thechannel, "mWT_full",     "ttbarregion_highpt",  thesample, systtype,   mTW,    			        evtweight, flavtag);
          if(mTW > 40)
          {
            fillHisto(thechannel, "mWT",        "ttbarregion_highpt",  thesample, systtype,   mTW,    			             evtweight, flavtag);
            fillHisto(thechannel, "MET",        "ttbarregion_highpt",  thesample, systtype,   met_pt, 			             evtweight, flavtag);
            fillHisto(thechannel, "LeptPt",     "ttbarregion_highpt",  thesample, systtype,   smalltree_lept_pt[0], 	     evtweight, flavtag);
            fillHisto(thechannel, "LeptEta",    "ttbarregion_highpt",  thesample, systtype,   smalltree_lept_eta[0], 	     evtweight, flavtag);
            fillHisto(thechannel, "DeltaPhiLJ", "ttbarregion_highpt",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),evtweight, flavtag);
            fillHisto(thechannel, "DeltaRLJ",   "ttbarregion_highpt",  thesample, systtype,   lept.DeltaR(leadingJet), 	     evtweight, flavtag);
          }
        }

        //***********************************
        //ttbar enriched region lowjetpt
        //***********************************
        if( njets >= 4 && nbjets ==2 && mTW > 40 ){

          if(sample == "QCD_A" || sample == "QCD_B" || sample == "QCD_C" || sample == "QCD_D")                          evtweight = SF_QCD_TT;
          fillHisto(thechannel, "mWT",        "ttbarregion_lowjetpt",  thesample, systtype,   mTW,    			        evtweight, flavtag);
          fillHisto(thechannel, "MET",        "ttbarregion_lowjetpt",  thesample, systtype,   met_pt, 			        evtweight, flavtag);
          fillHisto(thechannel, "LeptPt",     "ttbarregion_lowjetpt",  thesample, systtype,   smalltree_lept_pt[0], 	evtweight, flavtag);
          fillHisto(thechannel, "LeptEta",    "ttbarregion_lowjetpt",  thesample, systtype,   smalltree_lept_eta[0], 	evtweight, flavtag);

  	        for(int ijet=0; ijet<iter_jets; ijet++){
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",     "ttbarregion_lowjetpt",  thesample, systtype,  jet_pt[ijet], 		evtweight, flavtag);
                fillHisto(thechannel, "JetEta",    "ttbarregion_lowjetpt",  thesample, systtype,  jet_eta[ijet], 		evtweight, flavtag);
            }

          fillHisto(thechannel, "JetPt1",     "ttbarregion_lowjetpt",  thesample, systtype,  jet_pt[0], 			evtweight, flavtag);
          fillHisto(thechannel, "JetEta1",    "ttbarregion_lowjetpt",  thesample, systtype,  jet_eta[0], 			evtweight, flavtag);

          fillHisto(thechannel, "JetPt2",     "ttbarregion_lowjetpt",  thesample, systtype,  jet_pt[1], 			evtweight, flavtag);
          fillHisto(thechannel, "JetEta2",    "ttbarregion_lowjetpt",  thesample, systtype,  jet_eta[1], 			evtweight, flavtag);

          fillHisto(thechannel, "JetPt3",     "ttbarregion_lowjetpt",  thesample, systtype,  jet_pt[2], 			evtweight, flavtag);
          fillHisto(thechannel, "JetEta3",    "ttbarregion_lowjetpt",  thesample, systtype,  jet_eta[2], 			evtweight, flavtag);

          fillHisto(thechannel, "JetPt4",     "ttbarregion_lowjetpt",  thesample, systtype,  jet_pt[3], 			evtweight, flavtag);
          fillHisto(thechannel, "JetEta4",    "ttbarregion_lowjetpt",  thesample, systtype,  jet_eta[3], 			evtweight, flavtag);

        }
    }


  } //lepton selection
}




//------------------------------------------------------
//initialize the historams for the analysis
//------------------------------------------------------


void TreeReader::initializeHisto(TString sample, TString syst, bool isfirstset, TString flavtag){

  cout << endl;
  cout << endl;
  cout << "Initializing histograms for sample: " << sample << " ... " << endl;

  if(isfirstset){
    numb_histo = 0;
    TH1F * first_emptyHisto = new TH1F("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
    histo_list_mujets.push_back(first_emptyHisto);

    numb_histo++;
  }
  addHisto("NVtx",                  "",     sample.Data(), syst.Data(),    60,    0,    60, flavtag);

  //after lepton selection
  addHisto("NJet",      "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag);
  addHisto("NBJet",     "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag);
  addHisto("mWT",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("mW",        "afterleptsel", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("MET",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("mWTplusMET","afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("JetPt",     "afterleptsel",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",    "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("LeptPt",    "afterleptsel",  	sample.Data(), syst.Data(),   100,  0.,   200, flavtag);
  addHisto("LeptEta",   "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("BTagProba", "afterleptsel",  	sample.Data(), syst.Data(),   2,  -0.5,   1.5, flavtag);

  addHisto("JetPt",     "Wregion_highpt",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",    "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("ptW",       "Wregion_highpt", 	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("mWT",       "Wregion_highpt",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("mWT_full",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("MET",       "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("LeptPt",    "Wregion_highpt",  	sample.Data(), syst.Data(),   250,  0.,   500, flavtag);
  addHisto("LeptEta",   "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("DeltaPhiLJ","Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLMet","Wregion_highpt",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("mWT",       "1bjetregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("MET",       "1bjetregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("LeptPt",    "1bjetregion",  	sample.Data(), syst.Data(),   100,  0.,   200, flavtag);
  addHisto("LeptEta",   "1bjetregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("DeltaPhiLJ","1bjetregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",  "1bjetregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("mWT",       "QCDnormWregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "QCDnormBregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "QCDnormSregion",  	sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "QCDnormTTregion", 	sample.Data(), syst.Data(),   50,    0,   300, flavtag);

  addHisto("mWT",       "qcdWregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("MET",       "qcdWregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "qcdBregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("MET",       "qcdBregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "qcdSregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("MET",       "qcdSregion",  	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("mWT",       "qcdTTregion", 	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);
  addHisto("MET",       "qcdTTregion", 	    sample.Data(), syst.Data(),   50,    0,   300, flavtag);

  addHisto("JetPt",     "signalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",    "signalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("ptW",       "signalregion", 	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("mWT",       "signalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("mWT_full",  "signalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("MET",       "signalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("DeltaPhiLJ","signalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",  "signalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("JetPt",         "Selectedsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",        "Selectedsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("ptW",           "Selectedsignalregion", 	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("mWT",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("mWT_full",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("MET",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("DeltaPhiLJ",    "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiJMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLJMet", "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("JetPt",         "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",        "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("ptW",           "ATLASFCNCsignalregion", 	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("mWT",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("MET",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("mWTplusMET",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("DeltaPhiLJ",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiJMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLJMet", "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",      "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("JetPt",         "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("JetEta",        "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("ptW",           "ATLASRESsignalregion", 	sample.Data(), syst.Data(),   10,    0,   300, flavtag);
  addHisto("mWT",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("MET",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("mWTplusMET",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("DeltaPhiLJ",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiJMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaPhiLJMet", "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",      "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("mWT",       "ttbarregion_highpt",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("mWT_full",  "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("MET",       "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("LeptPt",    "ttbarregion_highpt",  	sample.Data(), syst.Data(),   100,  0.,   200, flavtag);
  addHisto("LeptEta",   "ttbarregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("DeltaPhiLJ","ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);
  addHisto("DeltaRLJ",  "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag);

  addHisto("mWT",       "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   23,   40,   500, flavtag);
  addHisto("MET",       "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   25,    0,   500, flavtag);
  addHisto("JetPt",     "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,   0,   300, flavtag);
  addHisto("JetEta",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);
  addHisto("LeptPt",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,  0.,   200, flavtag);
  addHisto("LeptEta",   "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);


  addHisto("JetPt1",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,   0,   450, flavtag);
  addHisto("JetEta1",   "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);

  addHisto("JetPt2",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,   0,   450, flavtag);
  addHisto("JetEta2",   "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);

  addHisto("JetPt3",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,   0,   450, flavtag);
  addHisto("JetEta3",   "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);

  addHisto("JetPt4",    "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   100,   0,   450, flavtag);
  addHisto("JetEta4",   "ttbarregion_lowjetpt", sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag);


  cout << "      Histograms properly initialized!   " << endl;
  cout << endl;
}


//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------

void TreeReader::addHisto(TString var, TString selstep, TString sample, TString syst, int nbins, float min, float max, TString flavtag){


  TString name_mujets;
  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(syst == "") name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag;
      else           name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"__"+syst;
  }
  else if (isWExcl)
  {
      if(syst == "") name_mujets =  var+"_mujets_"+selstep+"__WExcl";
      else           name_mujets =  var+"_mujets_"+selstep+"__WExcl"+"__"+syst;
  }
  else if (isWIncl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(syst == "") name_mujets =  var+"_mujets_"+selstep+"__WIncl"+flavtag;
      else           name_mujets =  var+"_mujets_"+selstep+"__WIncl"+flavtag+"__"+syst;
  }
  else if (isWIncl)
  {
      if(syst == "") name_mujets =  var+"_mujets_"+selstep+"__WIncl";
      else           name_mujets =  var+"_mujets_"+selstep+"__WIncl"+"__"+syst;
  }
  else               name_mujets =  var+"_mujets_"+selstep+"__"+sample;


  TH1F * thehisto_mujets = new TH1F(name_mujets,name_mujets,nbins,min,max);
  //FIXME add if(mWT) thehisto_mujets->Overflow...
  thehisto_mujets->Sumw2();
  histo_list_mujets.push_back(thehisto_mujets);
  histo_map_mujets[name_mujets.Data()] = numb_histo;

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
void TreeReader::fillHisto(TString channel, TString var, TString selstep, TString sample, TString syst, float val, float weight, TString flavtag){

  TString name;
  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(syst == "") name =  var+"_mujets_"+selstep+"__WExcl"+flavtag;
      else           name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"__"+syst;
  }
  else if (isWExcl)
  {
      if(syst == "") name =  var+"_mujets_"+selstep+"__WExcl";
      else           name =  var+"_mujets_"+selstep+"__WExcl"+"__"+syst;
  }
  else if (isWIncl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(syst == "") name =  var+"_mujets_"+selstep+"__WIncl"+flavtag;
      else           name =  var+"_mujets_"+selstep+"__WIncl"+flavtag+"__"+syst;
  }
  else if (isWIncl)
  {
      if(syst == "") name =  var+"_mujets_"+selstep+"__WIncl";
      else           name =  var+"_mujets_"+selstep+"__WIncl"+"__"+syst;
  }
  else               name =  var+"_mujets_"+selstep+"__"+sample;


  if(channel == "mujets" && histo_map_mujets[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histogram " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }

  if(channel == "mujets")     histo_list_mujets[histo_map_mujets[name.Data()]]->Fill(val, weight);

}


//normalize the W histos with the weights
void TreeReader::scaleHisto(TString channel, TString sample, TH1F* histoWCorrWeights){

  TString name_TT = "mWT_"+channel+"_ttbarregion_highpt__"+sample;
  TString name_W  = "mWT_"+channel+"_Wregion_highpt__"+sample;
  TString name_S  = "mWT_"+channel+"_signalregion__"+sample;
  TString name_B  = "mWT_"+channel+"_1bjetregion__"+sample;

  if(channel == "mujets" && (histo_map_mujets[name_TT.Data()] == 0  || histo_map_mujets[name_W.Data()] == 0  || histo_map_mujets[name_S.Data()] == 0  || histo_map_mujets[name_B.Data()] == 0  )) {
    cout << "   warning trying to scale a non existing histogram " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name_TT << endl;
  }

      histo_list_mujets[histo_map_mujets[ name_TT.Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ name_W.Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ name_S.Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ name_B.Data()]]->Multiply(histoWCorrWeights);

      histo_list_mujets[histo_map_mujets[ ( name_TT+ "__W__minus" ).Data()]];
      histo_list_mujets[histo_map_mujets[ ( name_W + "__W__minus" ).Data()]];
      histo_list_mujets[histo_map_mujets[ ( name_S + "__W__minus" ).Data()]];
      histo_list_mujets[histo_map_mujets[ ( name_B + "__W__minus" ).Data()]];

      histoWCorrWeights->Scale(2);
      histo_list_mujets[histo_map_mujets[ ( name_TT+ "__W__plus" ).Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ ( name_W + "__W__plus" ).Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ ( name_S + "__W__plus" ).Data()]]->Multiply(histoWCorrWeights);
      histo_list_mujets[histo_map_mujets[ ( name_B + "__W__plus" ).Data()]]->Multiply(histoWCorrWeights);

}


void TreeReader::deleteHisto(){
}

//SetUp CSV reweighting
void TreeReader::SetUpCSVreweighting(){

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
double TreeReader::GetCSVweight(const int iSys, int jet_n,
				float *jet_pt,float *jet_eta,float *jet_btagdiscri,int *jet_flav){
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

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

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

TH1F* TreeReader::getWCorrWeights(vector<TString> datalist, vector<TString> mclist, TString normregion, TString flavtag )
{

    //--------------------
    //loop over Datasamples
    //--------------------
    TString fileToGetData_name;
    TFile* fileToGetData = 0;
    TString histo_data_name;
    TH1F* histo_data = 0;
    for(unsigned int dataSample = 0; dataSample < datalist.size() ; dataSample++)
    {
        fileToGetData_name  = "outputroot_withQCDcorr/histofile_"+datalist[dataSample]+".root";
        fileToGetData       = new TFile(fileToGetData_name);
        histo_data_name     = normregion+"__"+datalist[dataSample];
        TH1F* histo_tmp     = (TH1F*)fileToGetData->Get(histo_data_name.Data());

        if (dataSample == 0) histo_data = (TH1F*) histo_tmp->Clone();
        else                 histo_data->Add( histo_tmp );
    }

    //--------------------
    //loop over MCsamples
    //--------------------
    TString fileToGetMC_name;
    TFile* fileToGetMC = 0;
    TString histo_mc_name;
    TH1F* histo_W = 0;
    TFile* fileToGetW = 0;
    for(unsigned int mcSample = 0; mcSample < mclist.size() ; mcSample++)
    {
        if     ( mclist[mcSample] == "WJets" )
        {
            fileToGetW    = new TFile("outputroot_withQCDcorr/histofile_WIncl.root");
            if( flavtag == "allflav")   histo_W    =     (TH1F*)fileToGetW->Get( (normregion+"__WIncl").Data() )->Clone();
            else
            {
                                        histo_W    =     (TH1F*)fileToGetW->Get( (normregion+"__WInclb").Data() )->Clone();
                                        histo_W->Add(    (TH1F*)fileToGetW->Get( (normregion+"__WInclc").Data() )->Clone() );
                                        histo_W->Add(    (TH1F*)fileToGetW->Get( (normregion+"__WIncll").Data() )->Clone() );
            }
        }
        else if( mclist[mcSample] == "W0Jets")
        {
            fileToGetW    = new TFile("outputroot_withQCDcorr/histofile_WExcl.root");
            if( flavtag == "allflav")   histo_W    =     (TH1F*)fileToGetW->Get( (normregion+"__WExcl").Data() )->Clone();
            else
            {
                                        histo_W    =     (TH1F*)fileToGetW->Get( (normregion+"__WExclb").Data() )->Clone();
                                        histo_W->Add(    (TH1F*)fileToGetW->Get( (normregion+"__WExclc").Data() )->Clone() );
                                        histo_W->Add(    (TH1F*)fileToGetW->Get( (normregion+"__WExcll").Data() )->Clone() );
            }
        }
        else if (mclist[mcSample] == "W1Jets" || mclist[mcSample] == "W2Jets" || mclist[mcSample] == "W3Jets" || mclist[mcSample] == "W4Jets") {}
        else
        {
            fileToGetMC_name    = "outputroot_withQCDcorr/histofile_"+mclist[mcSample]+".root";
            fileToGetMC         = new TFile(fileToGetMC_name);
            histo_mc_name       = normregion+"__"+mclist[mcSample];
            histo_data->Add( (TH1F*)fileToGetMC->Get( (histo_mc_name.Data()) )->Clone(), -1);
        }

   }

    //--------------------
    //loop over QCDsamples
    //--------------------
    std::vector<TString> qcdcorrectedlist;
    qcdcorrectedlist.push_back("QCD_A");
    qcdcorrectedlist.push_back("QCD_B");
    qcdcorrectedlist.push_back("QCD_C");
    qcdcorrectedlist.push_back("QCD_D");

    TString fileToGetQCD_name;
    TFile* fileToGetQCD = 0;
    TString histo_qcd_name;
    for(unsigned int qcdSample = 0; qcdSample < qcdcorrectedlist.size() ; qcdSample++)
    {
        fileToGetQCD_name    = "outputroot_withQCDcorr/histofile_"+qcdcorrectedlist[qcdSample]+".root";
        fileToGetQCD         = new TFile(fileToGetQCD_name);
        histo_qcd_name       = normregion+"__"+qcdcorrectedlist[qcdSample];
        histo_data->Add( (TH1F*)fileToGetQCD->Get( (histo_qcd_name.Data()) )->Clone(), -1);
    }

    if( histo_data != 0 && histo_W != 0) histo_data->Divide(histo_W);
    else                                 cout << "There is a problem in the W-normalization!" << endl;
    return histo_data;
}

double TreeReader::getQCDscalefactor(vector<TString> datalist, vector<TString> datalist_longnames, vector<TString> mclist, TString normregion, TString flavtag )
{
  double SF = 0;

    //--------------------
    //load QCDdata mTW file
    //--------------------
    double yieldQCD = 0;
    TString fileToGetQCDData_name;
    TFile* fileToGetQCDData = 0;
    TString histo_qcddatadriven_name;
    TH1F* histo_qcddatadriven = 0;
    for(unsigned int dataSample = 0; dataSample < datalist_longnames.size() ; dataSample++)
    {
        //fileToGetQCDData_name     = "outputroot_QCDcorr/histofile_"+datalist_longnames[dataSample]+".root";
        fileToGetQCDData_name     = "outputroot_QCDcorr_iso0p5/histofile_"+datalist_longnames[dataSample]+".root";
        fileToGetQCDData          = new TFile(fileToGetQCDData_name);
        histo_qcddatadriven_name  = normregion+"__"+datalist_longnames[dataSample];
        histo_qcddatadriven       = (TH1F*)fileToGetQCDData->Get(histo_qcddatadriven_name);

        if (histo_qcddatadriven != 0)
        {
            yieldQCD += histo_qcddatadriven->Integral();
        }
        else cout << "histo_qcddatadriven: " << histo_qcddatadriven_name << " not found! " << endl;
        fileToGetQCDData->Close();
    }

    //--------------------
    //loop over Datasamples
    //--------------------
    double yieldData = 0;
    TString fileToGetData_name;
    TFile* fileToGetData = 0;
    TString histo_data_name;
    TH1F* histo_data = 0;
    for(unsigned int dataSample = 0; dataSample < datalist.size() ; dataSample++)
    {
        fileToGetData_name  = "outputroot_woQCDcorr/histofile_"+datalist[dataSample]+".root";
        fileToGetData       = new TFile(fileToGetData_name);
        histo_data_name     = normregion+"__"+datalist[dataSample];
        histo_data          = (TH1F*)fileToGetData->Get(histo_data_name);

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
    TH1F* histo_mc = 0;
    for(unsigned int mcSample = 0; mcSample < mclist.size() ; mcSample++)
    {
        // Patch to use TTMSDecays instead of TTbar_Madgraph to make use of the right systematics (match, scale, mass)
        //if(mclist[0] != "TTbar_Madgraph") mclist[0] = "TTbar_Madgraph";
        //cout << "sample= " << mclist[mcSample] << " | flavtag= " << flavtag << endl;
        if     ( mclist[mcSample] == "WJets" )
        {
            fileToGetMC    = new TFile("outputroot_woQCDcorr/histofile_WIncl.root");
            if( flavtag == "allflav")   histo_mc    =     (TH1F*)fileToGetMC->Get( (normregion+"__WIncl").Data() )->Clone();
            else
            {
                                        histo_mc    =     (TH1F*)fileToGetMC->Get( (normregion+"__WInclb").Data() )->Clone();
                                        histo_mc->Add(    (TH1F*)fileToGetMC->Get( (normregion+"__WInclc").Data() )->Clone() );
                                        histo_mc->Add(    (TH1F*)fileToGetMC->Get( (normregion+"__WIncll").Data() )->Clone() );
            }
        }
        else if( mclist[mcSample] == "W0Jets")
        {
            fileToGetMC    = new TFile("outputroot_woQCDcorr/histofile_WExcl.root");
            if( flavtag == "allflav")   histo_mc    =     (TH1F*)fileToGetMC->Get( (normregion+"__WExcl").Data() )->Clone();
            else
            {
                                        histo_mc    =     (TH1F*)fileToGetMC->Get( (normregion+"__WExclb").Data() )->Clone();
                                        histo_mc->Add(    (TH1F*)fileToGetMC->Get( (normregion+"__WExclc").Data() )->Clone() );
                                        histo_mc->Add(    (TH1F*)fileToGetMC->Get( (normregion+"__WExcll").Data() )->Clone() );
            }
        }
        else if (mclist[mcSample] == "W1Jets" || mclist[mcSample] == "W2Jets" || mclist[mcSample] == "W3Jets" || mclist[mcSample] == "W4Jets") continue;
        else
        {
            fileToGetMC_name    = "outputroot_woQCDcorr/histofile_"+mclist[mcSample]+".root";
            fileToGetMC         = new TFile(fileToGetMC_name);
            histo_mc_name       = normregion+"__"+mclist[mcSample];
            histo_mc            = (TH1F*)fileToGetMC->Get(histo_mc_name);
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
       SFvalue = 1;
       SFerror = 0;
   }

   SF.push_back(SFvalue);
   SF.push_back(SFerror);
   return SF;

}

