#define TreeReader_cxx
#include "TreeReader.h"
#include "common.h"
#include <string>
#include <TGraphAsymmErrors.h>

using namespace std;


void TreeReader::Loop(short int CorrOption, std::vector<TString> datalist, std::vector<TString> datalist_longnames, std::vector<TString> mclist, TString sample, std::vector<TString> thesystlist, TString flavtag, short int systType)
{
   TString thesample  = sample;
   TFile * theoutputfile = 0;
   isQCD    = (sample == "QCD" );
   isW      = (sample == "WJets"  || sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets" || sample == "WJets_matchingdown" || sample == "WJets_matchingup" || sample == "WJets_scaledown" ||  sample == "WJets_scaleup" );
   isWExcl  = (sample == "W0Jets" || sample == "W1Jets" || sample == "W2Jets" || sample == "W3Jets" || sample == "W4Jets" || sample == "WJets_matchingdown" || sample == "WJets_matchingup" || sample == "WJets_scaledown" ||  sample == "WJets_scaleup" );
   isWIncl  = (sample == "WJets");
   TString sampleroot;
   if      ( isW && flavtag == "b") sampleroot = thesample+"_bflav";
   else if ( isW && flavtag == "c") sampleroot = thesample+"_cflav";
   else if ( isW && flavtag == "l") sampleroot = thesample+"_lflav";
   else                             sampleroot = thesample;
   double qcdIsoCut = 0.5;

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

   for(unsigned int i=0; i< systlist.size(); i++){
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

   double SF_QCD_W = 0, SF_QCD_B = 0, SF_QCD_S = 0, SF_QCD_TT = 0, SF_QCD_L = 0;
   if ((CorrOption == 2 || CorrOption == 3) && isQCD)
   {
       SF_QCD_L  = getQCDscalefactor(datalist, datalist_longnames, mclist, "mWT_mujets_QCDnormLregion" , flavtag );
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
      if      (CorrOption == 0)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0], thesample, flavtag, systType);
      else if (CorrOption == 1)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0], thesample, flavtag, systType);
      else if (CorrOption == 2)         applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[0], thesample, flavtag, systType);
      else if (CorrOption == 3)
      {
          for(unsigned int isyst=0; isyst< systlist.size(); isyst++)
          {
             // if( (systlist[isyst] == "toppt__plus" || systlist[isyst] == "toppt__minus") && thesample != "TTMSDecays_central" ) continue;
              applyEventSel(CorrOption, SF_QCD_L, SF_QCD_W, SF_QCD_B, SF_QCD_S, SF_QCD_TT, SFtrigger, SFtriggerError, thechannel, systlist[isyst], thesample, flavtag, systType);
          }
      }
      else cout << "ERROR: Wrong value of CorrOption! Allowed valued: 0,1,2,3" << endl;


   }

/*
   if ( CorrOption == 3 )
   {
       TH1D* histoWCorrWeights = getWCorrWeights();

       TCanvas *c1 = new TCanvas("c1","c1", 1000, 800);
       c1->cd();
       histoWCorrWeights->Draw();
       c1->SaveAs("WtailCorrWeights.root");
       c1->Close();
   }
*/

   manageOverflows();

   cout << endl;
   theoutputfile->Write();
   deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;

}



bool TreeReader::applyEventSel(short int CorrOption, double SF_QCD_L, double SF_QCD_W, double SF_QCD_B, double SF_QCD_S, double SF_QCD_TT, double SFtrigger, double SFtriggerError, TString thechannel, TString systtype, TString sample, TString flavtag, short int systType){

      bool applyCSV_reshape = 0;
      bool applyCSV         = 1;

      int  btagSys          = 1;
      int  mistagSys        = 1;

      double met_pt    = smalltree_met_pt;
      double met_phi   = smalltree_met_phi;
      double evtweight = -1;

      TString thesample = sample;
      if(systtype != "")               thesample = thesample + "__"+systtype;
      if(isQCD && systType == -2)      thesample = thesample + "__BgdContam__minus";
      else if(isQCD && systType == -1) thesample = thesample + "__Iso__minus";
      else if(isQCD && systType ==  1) thesample = thesample + "__Iso__plus";
      else if(isQCD && systType ==  2) thesample = thesample + "__BgdContam__plus";

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
      else if(systtype  == "mistag__plus" ){btagSys =  1;}
      else if(systtype  == "mistag__minus"){btagSys = -1;}
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

     double QCDContamweight = 1;
     if  (CorrOption == 2 || CorrOption == 3 )
     {
         if(sample == "QCD" && systType == 0)         QCDContamweight = smalltree_evtweight_nominal;
         else if(sample == "QCD" && systType == -2)   QCDContamweight = smalltree_evtweight_minus;
         else if(sample == "QCD" && systType ==  2)   QCDContamweight = smalltree_evtweight_plus;
         if(sample == "QCD" )                         evtweight = QCDContamweight;
 //if (QCDContamweight != 1) cout << "Evtweight= " << QCDContamweight << endl;
     }
     if (!isfinite(evtweight)) return false;

     fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 1 ,  evtweight, flavtag, systType);


     if(smalltree_lept_pt[0] > 33 && fabs(smalltree_lept_eta[0]) < 2.1 ){

       fillHisto(thechannel, "NVtx",        "",  thesample, systtype, smalltree_nvertex ,  evtweight, flavtag, systType);
       fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 2 ,  evtweight, flavtag, systType);
       TLorentzVector lept, met, leadingJet;
       lept.SetPtEtaPhiM(smalltree_lept_pt[0], smalltree_lept_eta[0],  smalltree_lept_phi[0], 0.);
       met.SetPtEtaPhiM(met_pt, 0, met_phi , 0.);
       double mTW = pow( 2*smalltree_lept_pt[0]*met_pt*(1-cos(smalltree_lept_phi[0] -  met_phi)) ,0.5);

       if(iter_jets>0) leadingJet.SetPtEtaPhiM(jet_pt[0], jet_eta[0], jet_phi[0], 0);
       if(CorrOption != 0 && (sample == "QCD"))          evtweight = QCDContamweight*SF_QCD_L;
       //if(CorrOption == 3 && isW && mTW > 375){if(mTW > 500) mTW = 499.5; evtweight *= getW_SF(mTW, systType);cout << "SF= " <<  getW_SF(mTW, systType) << " | evtweight= " << evtweight << endl; }
       if(CorrOption == 3 && isW && mTW >= 280 && systtype == "")                evtweight *= getW_SF(mTW, systType);

       int njets=0;
       int nbjets = 0;
       bool jetsup70 = false;
       std::vector<int > btagged_jet_idx;
       fillHisto(thechannel, "LeadJetBtagDiscr",     "afterleptsel",  thesample, systtype,  jet_btagdiscri[0] ,  evtweight, flavtag, systType);
       for(int ijet=0; ijet<iter_jets; ijet++){
         if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
         njets++;
         fillHisto(thechannel, "JetPt",     "afterleptsel",  thesample, systtype,  jet_pt[ijet] ,  evtweight, flavtag, systType);
         fillHisto(thechannel, "JetEta",    "afterleptsel",  thesample, systtype,  jet_eta[ijet] , evtweight, flavtag, systType);
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
               fillHisto(thechannel, "BTagProba", "afterleptsel", thesample, systtype, isbtag , evtweight, flavtag, systType);
           }
           if(applyCSV && isbtag) {nbjets++; btagged_jet_idx.push_back(ijet); }
         }

	     if(jet_pt[ijet] > 70 ) jetsup70 = true;
       }

       if (njets > 0) fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 3 ,  evtweight, flavtag, systType);
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

       fillHisto(thechannel, "NJet",       "afterleptsel",  thesample, systtype,   njets,   		        evtweight, flavtag, systType);
       fillHisto(thechannel, "NBJet",      "afterleptsel",  thesample, systtype,   nbjets, 		            evtweight, flavtag, systType);

       fillHisto(thechannel, "ptW",        "afterleptsel",  thesample, systtype,   theWcand.Pt(),       	evtweight, flavtag, systType);
       fillHisto(thechannel, "etaW",       "afterleptsel",  thesample, systtype,   theWcand.Eta(),       	evtweight, flavtag, systType);

       fillHisto(thechannel, "mWT",        "afterleptsel",  thesample, systtype,   mTW,    		            evtweight, flavtag, systType);
       fillHisto(thechannel, "mW",         "afterleptsel",  thesample, systtype,   theWcand.M(),  	        evtweight, flavtag, systType);
       fillHisto(thechannel, "MET",        "afterleptsel",  thesample, systtype,   met_pt, 		            evtweight, flavtag, systType);
       fillHisto(thechannel, "mWTplusMET", "afterleptsel",  thesample, systtype,   mWTplusMET, 		        evtweight, flavtag, systType);
       fillHisto(thechannel, "LeptPt",     "afterleptsel",  thesample, systtype,   smalltree_lept_pt[0], 	evtweight, flavtag, systType);
       fillHisto(thechannel, "LeptEta",    "afterleptsel",  thesample, systtype,   smalltree_lept_eta[0],   evtweight, flavtag, systType);

        //***********************************
        //qcd region
        //***********************************
    if(CorrOption == 0){

        fillHisto(thechannel, "mWT",          "qcdLregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
        fillHisto(thechannel, "MET",          "qcdLregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);

        if(jetsup70 && njets == 1 && nbjets == 0){
          fillHisto(thechannel, "mWT",        "qcdWregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto(thechannel, "MET",        "qcdWregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        else if(njets == 1 && nbjets == 1){
          fillHisto(thechannel, "mWT",        "qcdBregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto(thechannel, "MET",        "qcdBregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        if(jetsup70 && njets == 1 && nbjets == 1){
          fillHisto(thechannel, "mWT",        "qcdSregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto(thechannel, "MET",        "qcdSregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
        else if(jetsup70 && njets >= 4 && nbjets ==2 ){
          fillHisto(thechannel, "mWT",        "qcdTTregion",  thesample, systtype,   mTW,    			    evtweight, flavtag, systType);
          fillHisto(thechannel, "MET",        "qcdTTregion",  thesample, systtype,   met_pt, 			    evtweight, flavtag, systType);
        }
	}

        //***********************************
        //QCD normalization region
        //***********************************
    if(CorrOption == 0 || CorrOption == 1){

        if(mTW <= 40) fillHisto(thechannel, "mWT",        "QCDnormLregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);

        if(jetsup70 && njets == 1 && nbjets == 0 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormWregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        else if(njets == 1 && nbjets == 1 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormBregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        if(jetsup70 && njets == 1 && nbjets == 1 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormSregion",   thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
        else if(jetsup70 && njets >= 4 && nbjets == 2 && mTW <= 40){
          fillHisto(thechannel, "mWT",        "QCDnormTTregion",  thesample, systtype,   mTW,    			evtweight, flavtag, systType);
        }
	}

    if(CorrOption == 1 || CorrOption == 2 || CorrOption == 3)
    {
        //***********************************
        //1bjet region
        //***********************************
        if(nbjets == 1 && mTW > 40 ){

          if(sample == "QCD")                      evtweight = QCDContamweight*SF_QCD_B;
          fillHisto(thechannel, "mWT",        "1bjetregion",  thesample, systtype,   mTW,    			            evtweight, flavtag, systType);
          fillHisto(thechannel, "MET",        "1bjetregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag, systType);
          fillHisto(thechannel, "DeltaPhiLJ", "1bjetregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),evtweight, flavtag, systType);
          fillHisto(thechannel, "DeltaRLJ",   "1bjetregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag, systType);
          fillHisto(thechannel, "LeptPt",     "1bjetregion",  thesample, systtype,   smalltree_lept_pt[0], 		    evtweight, flavtag, systType);
          fillHisto(thechannel, "LeptEta",    "1bjetregion",  thesample, systtype,   smalltree_lept_eta[0], 	    evtweight, flavtag, systType);
	    }

        //***********************************
        //W enriched region
        //***********************************
	    if(jetsup70 && njets == 1 && nbjets == 0){

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_W;
          fillHisto(thechannel, "mWT_full",   "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 12 ,  evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",  "Wregion_highpt",  thesample, systtype,  jet_pt[ijet] ,                 evtweight, flavtag, systType);
                fillHisto(thechannel, "JetEta", "Wregion_highpt",  thesample, systtype,  jet_eta[ijet] ,                evtweight, flavtag, systType);
            }
            fillHisto(thechannel, "NJet",       "Wregion_highpt",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto(thechannel, "NBJet",      "Wregion_highpt",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto(thechannel, "LeadJetBtagDiscr","Wregion_highpt",thesample, systtype,jet_btagdiscri[0] ,               evtweight, flavtag, systType);
            fillHisto(thechannel, "ptW",        "Wregion_highpt",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto(thechannel, "etaW",       "Wregion_highpt",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);
            //if(CorrOption == 3 && isW && mTW > 375 && systtype ==""){if(mTW > 500) mTW = 499.5; cout << "SF= " <<  getW_SF(mTW, systType) << " | evtweight= " << evtweight << endl; }
            fillHisto(thechannel, "mWT",        "Wregion_highpt",  thesample, systtype,   mTW,    			                evtweight, flavtag, systType);
            fillHisto(thechannel, "MET",        "Wregion_highpt",  thesample, systtype,   met_pt, 			                evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJ", "Wregion_highpt",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), 	evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLMet","Wregion_highpt", thesample, systtype,   abs(lept.DeltaPhi(met)),          evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaRLJ",   "Wregion_highpt",  thesample, systtype,   lept.DeltaR(leadingJet), 	        evtweight, flavtag, systType);
            fillHisto(thechannel, "LeptPt",     "Wregion_highpt",  thesample, systtype,   smalltree_lept_pt[0], 		    evtweight, flavtag, systType);
            fillHisto(thechannel, "LeptEta",    "Wregion_highpt",  thesample, systtype,   smalltree_lept_eta[0], 		    evtweight, flavtag, systType);
          }
	    }

        //***********************************
        //signal enriched region
        //***********************************
        if(jetsup70 && njets == 1 && nbjets == 1 ){

          if(sample == "QCD")          evtweight = QCDContamweight*SF_QCD_S;
          fillHisto(thechannel, "mWT_full",     "signalregion",  thesample, systtype,   mTW,    		evtweight, flavtag, systType);
          // Manage overflows in SR - put higher values in the last bin by hands
 //         if(mTW    >= 500) mTW    = 499.95;
 //         if(met_pt >= 500) met_pt = 499.95;
          if(mTW > 40)
          {
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",  "signalregion",  thesample, systtype,  jet_pt[ijet] ,             evtweight, flavtag, systType);
                fillHisto(thechannel, "JetEta", "signalregion",  thesample, systtype,  jet_eta[ijet] ,            evtweight, flavtag, systType);
            }
            fillHisto(thechannel, "NJet",       "signalregion",  thesample, systtype,   njets,  		          evtweight, flavtag, systType);
            fillHisto(thechannel, "NBJet",      "signalregion",  thesample, systtype,   nbjets, 		          evtweight, flavtag, systType);

            fillHisto(thechannel, "LeptPt",     "signalregion",  thesample, systtype,   smalltree_lept_pt[0], 	         evtweight, flavtag, systType);
            fillHisto(thechannel, "LeptEta",    "signalregion",  thesample, systtype,   smalltree_lept_eta[0], 		     evtweight, flavtag, systType);
            fillHisto(thechannel, "LeadJetBtagDiscr","signalregion",thesample, systtype,  jet_btagdiscri[0] ,            evtweight, flavtag, systType);
            fillHisto(thechannel, "ptW",        "signalregion",  thesample, systtype,   theWcand.Pt(),       	         evtweight, flavtag, systType);
            fillHisto(thechannel, "etaW",       "signalregion",  thesample, systtype,   theWcand.Eta(),       	         evtweight, flavtag, systType);
            fillHisto(thechannel, "mWT",        "signalregion",  thesample, systtype,   mTW,    			             evtweight, flavtag, systType);
            fillHisto(thechannel, "MET",        "signalregion",  thesample, systtype,   met_pt, 			             evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJ", "signalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),  evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaRLJ",   "signalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	     evtweight, flavtag, systType);
          }
        }


        //***********************************
        //FCNC signal enriched region
        //***********************************
        if(mTW > 40 && njets > 0)
        {
            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 4 ,  evtweight, flavtag, systType);
            // Manage overflows in SR - put higher values in the last bin by hands
//            if(mTW    >= 500) mTW    = 499.95;
//            if(met_pt >= 500) met_pt = 499.95;

            if(nbjets == 1)
            {
                fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 5 ,  evtweight, flavtag, systType);
                if(njets == 1)
                {
                    fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 6 ,  evtweight, flavtag, systType);
                    if(jetsup70)
                    {
                        fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 7 ,  evtweight, flavtag, systType);
                        if(met_pt > 100)
                        {
                            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 8 ,  evtweight, flavtag, systType);
                            if(theWcand.Pt() > 50)
                            {
                                fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 9 ,  evtweight, flavtag, systType);
                                if(abs(lept.DeltaPhi(leadingJet)) < 1.7)
                                {
                                    fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 10 ,  evtweight, flavtag, systType);
                                    if(sample == "QCD")          evtweight = QCDContamweight*SF_QCD_S;

                                    for(int ijet=0; ijet<iter_jets; ijet++)
                                    {
                                        if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                                        fillHisto(thechannel, "JetPt",      "Selectedsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                                        fillHisto(thechannel, "JetEta",     "Selectedsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
                                    }

                                        tree_EvtWeight  = evtweight;
                                        tree_wMass      = mTW;
                                        tree_met        = met_pt;
                                        tree_wPt        = theWcand.Pt();
                                        tree_deltaPhilb = abs(lept.DeltaPhi(leadingJet));

                                        if(theTree_map[thesample] != 0)  theTree_map[thesample]->Fill();

                                        fillHisto(thechannel, "LeadJetBtagDiscr","Selectedsignalregion", thesample, systtype,   jet_btagdiscri[0] ,             evtweight, flavtag, systType);
                                        fillHisto(thechannel, "NJet",           "Selectedsignalregion",  thesample, systtype,   njets, 	    	          evtweight, flavtag, systType);
                                        fillHisto(thechannel, "NBJet",          "Selectedsignalregion",  thesample, systtype,   nbjets, 		          evtweight, flavtag, systType);
                                        fillHisto(thechannel, "LeptPt",         "Selectedsignalregion",  thesample, systtype,   smalltree_lept_pt[0], 		    evtweight, flavtag, systType);
                                        fillHisto(thechannel, "LeptEta",        "Selectedsignalregion",  thesample, systtype,   smalltree_lept_eta[0], 		    evtweight, flavtag, systType);
                                        fillHisto(thechannel, "ptW",            "Selectedsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag, systType);
                                        fillHisto(thechannel, "etaW",           "Selectedsignalregion",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);
                                        fillHisto(thechannel, "mWT",            "Selectedsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag, systType);
                                        fillHisto(thechannel, "MET",            "Selectedsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag, systType);
                                        fillHisto(thechannel, "DeltaPhiLJ",     "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
                                        fillHisto(thechannel, "DeltaPhiLMet",   "Selectedsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
                                        fillHisto(thechannel, "DeltaPhiJMet",   "Selectedsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
                                        fillHisto(thechannel, "DeltaPhiLJMet",  "Selectedsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
                                        fillHisto(thechannel, "DeltaRLJ",       "Selectedsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag, systType);
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
          // Manage overflows in SR - put higher values in the last bin by hands
//          if(mTW    >= 500) mTW    = 499.95;
//          if(met_pt >= 500) met_pt = 499.95;
          if(mWTplusMET > 60 && mTW > 250 && abs(lept.DeltaPhi(leadingJet)) < 1.4)
          {
            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 13 , evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",      "ATLASFCNCsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                fillHisto(thechannel, "JetEta",     "ATLASFCNCsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
            }
            fillHisto(thechannel, "ptW",            "ATLASFCNCsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag, systType);
            fillHisto(thechannel, "mWT",            "ATLASFCNCsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag, systType);
            fillHisto(thechannel, "MET",            "ATLASFCNCsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag, systType);
            fillHisto(thechannel, "mWTplusMET",     "ATLASFCNCsignalregion",  thesample, systtype,   mWTplusMET, 		            evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJ",     "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiJMet",   "ATLASFCNCsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJMet",  "ATLASFCNCsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaRLJ",       "ATLASFCNCsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag, systType);
          }

          if(mWTplusMET > 60 && mTW > 210 && abs(lept.DeltaPhi(leadingJet)) < 1.2)
          {
            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 14 , evtweight, flavtag, systType);
            for(int ijet=0; ijet<iter_jets; ijet++)
            {
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",      "ATLASRESsignalregion",  thesample, systtype,  jet_pt[ijet] ,                   evtweight, flavtag, systType);
                fillHisto(thechannel, "JetEta",     "ATLASRESsignalregion",  thesample, systtype,  jet_eta[ijet] ,                  evtweight, flavtag, systType);
            }
            fillHisto(thechannel, "ptW",            "ATLASRESsignalregion",  thesample, systtype,   theWcand.Pt(),       	        evtweight, flavtag, systType);
            fillHisto(thechannel, "mWT",            "ATLASRESsignalregion",  thesample, systtype,   mTW,    			            evtweight, flavtag, systType);
            fillHisto(thechannel, "MET",            "ATLASRESsignalregion",  thesample, systtype,   met_pt, 			            evtweight, flavtag, systType);
            fillHisto(thechannel, "mWTplusMET",     "ATLASRESsignalregion",  thesample, systtype,   mWTplusMET, 		            evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJ",     "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)), evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(lept.DeltaPhi(met)),        evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiJMet",   "ATLASRESsignalregion",  thesample, systtype,   abs(leadingJet.DeltaPhi(met)),  evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJMet",  "ATLASRESsignalregion",  thesample, systtype,   abs(jetLept.DeltaPhi(met)),     evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaRLJ",       "ATLASRESsignalregion",  thesample, systtype,   lept.DeltaR(leadingJet), 	    evtweight, flavtag, systType);
          }

        }


        //***********************************
        //ttbar enriched region
        //***********************************
        if(jetsup70 && njets >= 4 && nbjets ==2){

          if(sample == "QCD")                          evtweight = QCDContamweight*SF_QCD_TT;
          fillHisto(thechannel, "mWT_full",     "ttbarregion_highpt",  thesample, systtype,   mTW,    			        evtweight, flavtag, systType);
          if(mTW > 40)
          {
            fillHisto(thechannel, "CutFlow",     "",  thesample, systtype, 11 ,  evtweight, flavtag, systType);
  	        for(int ijet=0; ijet<iter_jets; ijet++){
                if(jet_pt[ijet] < 30 || fabs(jet_eta[ijet]) > 2.4) continue;
                fillHisto(thechannel, "JetPt",     "ttbarregion_highpt",  thesample, systtype,  jet_pt[ijet], 		evtweight, flavtag, systType);
                fillHisto(thechannel, "JetEta",    "ttbarregion_highpt",  thesample, systtype,  jet_eta[ijet], 		evtweight, flavtag, systType);
            }

            fillHisto(thechannel, "LeadJetBtagDiscr","ttbarregion_highpt",thesample, systtype,  jet_btagdiscri[0] ,             evtweight, flavtag, systType);
            fillHisto(thechannel, "NJet",       "ttbarregion_highpt",  thesample, systtype,   njets, 		                    evtweight, flavtag, systType);
            fillHisto(thechannel, "NBJet",      "ttbarregion_highpt",  thesample, systtype,   nbjets, 		                    evtweight, flavtag, systType);
            fillHisto(thechannel, "ptW",        "ttbarregion_highpt",  thesample, systtype,   theWcand.Pt(),       	            evtweight, flavtag, systType);
            fillHisto(thechannel, "etaW",       "ttbarregion_highpt",  thesample, systtype,   theWcand.Eta(),       	        evtweight, flavtag, systType);

            fillHisto(thechannel, "mWT",        "ttbarregion_highpt",  thesample, systtype,   mTW,    			             evtweight, flavtag, systType);
            fillHisto(thechannel, "MET",        "ttbarregion_highpt",  thesample, systtype,   met_pt, 			             evtweight, flavtag, systType);
            fillHisto(thechannel, "LeptPt",     "ttbarregion_highpt",  thesample, systtype,   smalltree_lept_pt[0], 	     evtweight, flavtag, systType);
            fillHisto(thechannel, "LeptEta",    "ttbarregion_highpt",  thesample, systtype,   smalltree_lept_eta[0], 	     evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaPhiLJ", "ttbarregion_highpt",  thesample, systtype,   abs(lept.DeltaPhi(leadingJet)),evtweight, flavtag, systType);
            fillHisto(thechannel, "DeltaRLJ",   "ttbarregion_highpt",  thesample, systtype,   lept.DeltaR(leadingJet), 	     evtweight, flavtag, systType);
          }
        }

    }


  } //lepton selection
}




//------------------------------------------------------
//initialize the historams for the analysis
//------------------------------------------------------


void TreeReader::initializeHisto(TString sample, TString syst, bool isfirstset, TString flavtag, short int systType){

  cout << endl;
  cout << endl;
  cout << "Initializing histograms for sample: " << sample << " ... " << endl;

  if(isfirstset){
    numb_histo = 0;
    TH1D * first_emptyHisto = new TH1D("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
    histo_list_mujets.push_back(first_emptyHisto);

    numb_histo++;
  }
  addHisto("CutFlow",               "",     sample.Data(), syst.Data(),    15,    0,    15, flavtag, systType);
  addHisto("NVtx",                  "",     sample.Data(), syst.Data(),    60,    0,    60, flavtag, systType);

  //after lepton selection
  addHisto("NJet",      "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "afterleptsel",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("mWT",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mW",        "afterleptsel", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "afterleptsel", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "afterleptsel", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("MET",       "afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET","afterleptsel",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("JetPt",     "afterleptsel",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",    "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "afterleptsel",  	sample.Data(), syst.Data(),   40,   0.,   200, flavtag, systType);
  addHisto("LeptEta",   "afterleptsel",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("BTagProba", "afterleptsel",  	sample.Data(), syst.Data(),   2,  -0.5,   1.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "afterleptsel", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("NJet",      "Wregion_highpt",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "Wregion_highpt",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",     "Wregion_highpt",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",    "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",       "Wregion_highpt", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "Wregion_highpt", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",       "Wregion_highpt",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_full",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("LeptPt",    "Wregion_highpt",  	sample.Data(), syst.Data(),   25,   0.,   500, flavtag, systType);
  addHisto("LeptEta",   "Wregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet","Wregion_highpt",	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "Wregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "Wregion_highpt", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("mWT",       "1bjetregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("MET",       "1bjetregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("LeptPt",    "1bjetregion",  	sample.Data(), syst.Data(),   40,   0.,   200, flavtag, systType);
  addHisto("LeptEta",   "1bjetregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","1bjetregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "1bjetregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

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

  addHisto("NJet",      "signalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",     "signalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",     "signalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",    "signalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "signalregion",  	sample.Data(), syst.Data(),   50,  0.,    500, flavtag, systType);
  addHisto("LeptEta",   "signalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",       "signalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "signalregion", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",       "signalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_full",  "signalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "signalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ","signalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "signalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "signalregion", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("NJet",          "Selectedsignalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("NBJet",         "Selectedsignalregion",  	sample.Data(), syst.Data(),   5,  -0.5,   4.5, flavtag, systType);
  addHisto("JetPt",         "Selectedsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",        "Selectedsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",        "Selectedsignalregion",  	sample.Data(), syst.Data(),   50,  0.,    500, flavtag, systType);
  addHisto("LeptEta",       "Selectedsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "Selectedsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",          "Selectedsignalregion", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("mWT",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_full",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",           "Selectedsignalregion",  	sample.Data(), syst.Data(),   20,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "Selectedsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "Selectedsignalregion", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  addHisto("JetPt",         "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",        "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "ATLASFCNCsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWT",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("MET",           "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "ATLASFCNCsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

  addHisto("JetPt",         "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   10,    0,   300, flavtag, systType);
  addHisto("JetEta",        "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("ptW",           "ATLASRESsignalregion", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWT",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("MET",           "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("mWTplusMET",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("DeltaPhiLJ",    "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiJMet",  "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaPhiLJMet", "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",      "ATLASRESsignalregion",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);

  addHisto("NJet",      "ttbarregion_highpt",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("NBJet",     "ttbarregion_highpt",  	sample.Data(), syst.Data(),   8,  -0.5,   7.5, flavtag, systType);
  addHisto("mWT",       "ttbarregion_highpt",  	sample.Data(), syst.Data(),   23,   40,   500, flavtag, systType);
  addHisto("mWT_full",  "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("MET",       "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("ptW",       "ttbarregion_highpt", 	sample.Data(), syst.Data(),   25,    0,   500, flavtag, systType);
  addHisto("etaW",      "ttbarregion_highpt", 	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeptPt",    "ttbarregion_highpt",  	sample.Data(), syst.Data(),   40,   0.,   200, flavtag, systType);
  addHisto("LeptEta",   "ttbarregion_highpt",  	sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("DeltaPhiLJ","ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("DeltaRLJ",  "ttbarregion_highpt",  	sample.Data(), syst.Data(),   25,    0,  3.14, flavtag, systType);
  addHisto("JetPt",     "ttbarregion_highpt",   sample.Data(), syst.Data(),   50,    0,   300, flavtag, systType);
  addHisto("JetEta",    "ttbarregion_highpt",   sample.Data(), syst.Data(),   26, -2.5,   2.5, flavtag, systType);
  addHisto("LeadJetBtagDiscr", "ttbarregion_highpt", sample.Data(), syst.Data(),  20,  0,     1, flavtag, systType);

  cout << "      Histograms properly initialized!   " << endl;
  cout << endl;

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

void TreeReader::addHisto(TString var, TString selstep, TString sample, TString syst, int nbins, float min, float max, TString flavtag, short int systType){


  TString name_mujets;
  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(     syst == "" && systType == -1) name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_matchingdown";
      else if(syst == "" && systType == 1)  name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_matchingup";
      else if(syst == "" && systType == -2) name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_scaledown";
      else if(syst == "" && systType == 2)  name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_scaleup";
      else if(syst == "")                   name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag;
      else                                  name_mujets =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"__"+syst;
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


  TH1D * thehisto_mujets = new TH1D(name_mujets,name_mujets,nbins,min,max);
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
void TreeReader::fillHisto(TString channel, TString var, TString selstep, TString sample, TString syst, float val, float weight, TString flavtag, short int systType){

  TString name;
  if      (isWExcl && (flavtag == "b" || flavtag == "c" || flavtag == "l") )
  {
      if(     syst == "" && systType == -1) name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_matchingdown";
      else if(syst == "" && systType == 1)  name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_matchingup";
      else if(syst == "" && systType == -2) name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_scaledown";
      else if(syst == "" && systType == 2)  name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"_scaleup";
      else if(syst == "")                   name =  var+"_mujets_"+selstep+"__WExcl"+flavtag;
      else                                  name =  var+"_mujets_"+selstep+"__WExcl"+flavtag+"__"+syst;
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
  //if(var == "mWT" && val > 375 && selstep == "Wregion_highpt") cout << "histo = " << name << " val= " << val << " | weight= " << weight << endl;
  if(channel == "mujets")     histo_list_mujets[histo_map_mujets[name.Data()]]->Fill(val, weight);

}

void TreeReader::manageOverflows()
{

  if(histo_list_mujets.size() == 0)
  {
    cout << "   WARNING trying to manage a non existing histogram at line: " << __LINE__ << endl;
  }

  for ( unsigned short int hist = 0; hist < histo_list_mujets.size(); hist++)
  {
      double lastBinContent  = histo_list_mujets[hist]->GetBinContent(histo_list_mujets[hist]->GetNbinsX() );
      double overflowContent = histo_list_mujets[hist]->GetBinContent(histo_list_mujets[hist]->GetNbinsX() +1);
      histo_list_mujets[hist]->SetBinContent(histo_list_mujets[hist]->GetNbinsX() , lastBinContent + overflowContent);
  }
}



//normalize the W histos with the weights
void TreeReader::scaleHisto(TString channel, TString sample, TH1D* histoWCorrWeights)
{

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

TH1D* TreeReader::getWCorrWeights()
{

    //------------------
    //loop over samples
    //------------------
    TFile* fileToGetData      = new TFile("../Theta/inputTheta_merged.root");
    TH1D*  histo_data         = (TH1D*)fileToGetData->Get("mWTmujetsWregionHighpt__DATA")->Clone();

    //--------------------
    //loop over MCsamples
    //--------------------
    //TFile* fileToGetPostFitMC = new TFile("../Theta/histos-mle.root");
    TString histo_mc_name;
    TH1D* histo_W =     (TH1D*)fileToGetData->Get("mWTmujetsWregionHighpt__WExclb")->Clone();
          histo_W->Add( (TH1D*)fileToGetData->Get("mWTmujetsWregionHighpt__WExclc")->Clone() );
          histo_W->Add( (TH1D*)fileToGetData->Get("mWTmujetsWregionHighpt__WExcll")->Clone() );

    histo_data->Add( (TH1D*)fileToGetData->Get( "mWTmujetsWregionHighpt__TTMSDecays"   )->Clone(), -1);
    histo_data->Add( (TH1D*)fileToGetData->Get( "mWTmujetsWregionHighpt__DY"           )->Clone(), -1);
    histo_data->Add( (TH1D*)fileToGetData->Get( "mWTmujetsWregionHighpt__SingleTop"    )->Clone(), -1);
    histo_data->Add( (TH1D*)fileToGetData->Get( "mWTmujetsWregionHighpt__VV"           )->Clone(), -1);
    histo_data->Add( (TH1D*)fileToGetData->Get( "mWTmujetsWregionHighpt__QCD"          )->Clone(), -1);

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
    TH1D* histo_qcddatadriven = 0;
    fileToGetQCDData_name     = "outputroot_QCDcorr_iso0p5/histofile_QCDdatadriven.root";
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
        fileToGetData_name  = "outputroot_woQCDcorr/histofile_"+datalist[dataSample]+".root";
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
        fileToGetMC_name    = "outputroot_woQCDcorr/histofile_"+mclist[mcSample]+".root";
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
//cout << "Yield_data= " << yieldData << " | Yield_MC= " << yieldMC << " | Yield_QCD= " << yieldQCD << endl;
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

double TreeReader::getW_SF(double mTW, short int systType)
{
    if(systType == -1)      return 0.01;
    else if(systType == 0)  return 1;
    else if(systType == 1)  return 100;

    //if(systType == -1)      return 1;
    //else if(systType == 0)  return 53.4689 - 0.0624099*mTW;
    //else if(systType == 1)  return 2*(53.4689 - 0.0624099*mTW);
    //if(systType == -1)      return (2.23818-0.0041937*mTW)/2.;
    //else if(systType == 0)  return 2.23818-0.0041937*mTW;
    //else if(systType == 1)       return 1;
}
