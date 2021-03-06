//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 24 16:03:04 2014 by ROOT version 5.34/21
// from TTree SmallTree_WZ/SmallTree_WZ
// found on file: proof_WZ.root
//////////////////////////////////////////////////////////

#ifndef TreeReader_h
#define TreeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TGraphAsymmErrors.h>

#include "TLorentzVector.h"

#include <iostream>

#include "../BTagSFNotReshaping/BTagSFUtil.h"

// Header file for the classes stored in the TTree if any.

using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile *f_CSVwgt_HF;
   TFile *f_CSVwgt_LF;

   TH1D *h_csv_wgt_hf[100][100];
   TH1D *c_csv_wgt_hf[100][100];
   TH1D *h_csv_wgt_lf[100][100][100];

   bool useElectronChannel;
   bool isData;
   bool isW;
   bool isWExcl;
   bool isWIncl;
   bool isQCD;
   bool isTTbarSyst;
   bool isSignal;
   TString thechannel;

   // Declaration of leaf types
   Int_t           smalltree_nlepton;
   Float_t         smalltree_lept_pt[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_eta[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_phi[3];   //[smalltree_nlepton]
   Float_t         smalltree_lept_iso[3];   //[smalltree_nlepton]
   Int_t           smalltree_lept_flav[3];   //[smalltree_nlepton]
   Int_t           smalltree_njets;
   Float_t         smalltree_jet_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_up[100];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_down[100];   //[smalltree_njets]
   Int_t           smalltree_jet_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jesup_njets;
   Float_t         smalltree_jet_jesup_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jesup_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jesdown_njets;
   Float_t         smalltree_jet_jesdown_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jesdown_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jerup_njets;
   Float_t         smalltree_jet_jerup_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jerup_flav[100];   //[smalltree_njets]
   Int_t           smalltree_jerdown_njets;
   Float_t         smalltree_jet_jerdown_pt[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_eta[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_phi[100];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_btagdiscri[100];   //[smalltree_njets]
   Int_t           smalltree_jet_jerdown_flav[100];   //[smalltree_njets]
   Float_t         smalltree_met_jesup_pt;
   Float_t         smalltree_met_jesup_phi;
   Float_t         smalltree_met_jesdown_pt;
   Float_t         smalltree_met_jesdown_phi;
   Float_t         smalltree_met_jerup_pt;
   Float_t         smalltree_met_jerup_phi;
   Float_t         smalltree_met_jerdown_pt;
   Float_t         smalltree_met_jerdown_phi;
   Float_t         smalltree_met_unclsup_pt;
   Float_t         smalltree_met_unclsup_phi;
   Float_t         smalltree_met_unclsdown_pt;
   Float_t         smalltree_met_unclsdown_phi;
   Float_t         smalltree_met_pt;
   Float_t         smalltree_metnosmear_phi;
   Float_t         smalltree_metnosmear_pt;
   Float_t         smalltree_met_phi;
   Float_t         smalltree_weight_trigup;
   Float_t         smalltree_weight_trigdown;
   Float_t         smalltree_weight_leptup;
   Float_t         smalltree_weight_leptdown;
   Float_t         smalltree_weight_PDFup;
   Float_t         smalltree_weight_PDFdown;
   Float_t         smalltree_weight_PUup;
   Float_t         smalltree_weight_PUdown;
   Float_t         smalltree_evtweight;
   Float_t         smalltree_evtweight_nominal;
   Float_t         smalltree_evtweight_minus;
   Float_t         smalltree_evtweight_plus;
   Float_t         smalltree_weight_toppt;
   //Int_t           smalltree_IChannel;
   Int_t           smalltree_nvertex;
   Int_t           smalltree_tmeme;

   // List of branches
   TBranch        *b_smalltree_nlepton;   //!
   TBranch        *b_smalltree_lept_pt;   //!
   TBranch        *b_smalltree_lept_eta;   //!
   TBranch        *b_smalltree_lept_phi;   //!
   TBranch        *b_smalltree_lept_iso;   //!
   TBranch        *b_smalltree_lept_flav;   //!
   TBranch        *b_smalltree_njets;   //!
   TBranch        *b_smalltree_jet_pt;   //!
   TBranch        *b_smalltree_jet_eta;   //!
   TBranch        *b_smalltree_jet_phi;   //!
   TBranch        *b_smalltree_jet_btagdiscri;   //!
   TBranch        *b_smalltree_jet_btagdiscri_up;   //!
   TBranch        *b_smalltree_jet_btagdiscri_down;   //!
   TBranch        *b_smalltree_jet_flav;   //!
   TBranch        *b_smalltree_jesup_njets;   //!
   TBranch        *b_smalltree_jet_jesup_pt;   //!
   TBranch        *b_smalltree_jet_jesup_eta;   //!
   TBranch        *b_smalltree_jet_jesup_phi;   //!
   TBranch        *b_smalltree_jet_jesup_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jesup_flav;   //!
   TBranch        *b_smalltree_jesdown_njets;   //!
   TBranch        *b_smalltree_jet_jesdown_pt;   //!
   TBranch        *b_smalltree_jet_jesdown_eta;   //!
   TBranch        *b_smalltree_jet_jesdown_phi;   //!
   TBranch        *b_smalltree_jet_jesdown_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jesdown_flav;   //!
   TBranch        *b_smalltree_jerup_njets;   //!
   TBranch        *b_smalltree_jet_jerup_pt;   //!
   TBranch        *b_smalltree_jet_jerup_eta;   //!
   TBranch        *b_smalltree_jet_jerup_phi;   //!
   TBranch        *b_smalltree_jet_jerup_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jerup_flav;   //!
   TBranch        *b_smalltree_jerdown_njets;   //!
   TBranch        *b_smalltree_jet_jerdown_pt;   //!
   TBranch        *b_smalltree_jet_jerdown_eta;   //!
   TBranch        *b_smalltree_jet_jerdown_phi;   //!
   TBranch        *b_smalltree_jet_jerdown_btagdiscri;   //!
   TBranch        *b_smalltree_jet_jerdown_flav;   //!
   TBranch        *b_smalltree_met_jesup_pt;   //!
   TBranch        *b_smalltree_met_jesup_phi;   //!
   TBranch        *b_smalltree_met_jesdown_pt;   //!
   TBranch        *b_smalltree_met_jesdown_phi;   //!
   TBranch        *b_smalltree_met_jerup_pt;   //!
   TBranch        *b_smalltree_met_jerup_phi;   //!
   TBranch        *b_smalltree_met_jerdown_pt;   //!
   TBranch        *b_smalltree_met_jerdown_phi;   //!
   TBranch        *b_smalltree_met_unclsup_pt;   //!
   TBranch        *b_smalltree_met_unclsup_phi;   //!
   TBranch        *b_smalltree_met_unclsdown_pt;   //!
   TBranch        *b_smalltree_met_unclsdown_phi;   //!
   TBranch        *b_smalltree_met_pt;   //!
   TBranch        *b_smalltree_met_phi;   //!
   TBranch        *b_smalltree_metnosmear_pt;   //!
   TBranch        *b_smalltree_metnosmear_phi;   //!
   TBranch        *b_smalltree_weight_trigup;   //!
   TBranch        *b_smalltree_weight_trigdown;   //!
   TBranch        *b_smalltree_weight_leptup;   //!
   TBranch        *b_smalltree_weight_leptdown;   //!
   TBranch        *b_smalltree_weight_PDFup;   //!
   TBranch        *b_smalltree_weight_PDFdown;   //!
   TBranch        *b_smalltree_weight_PUup;   //!
   TBranch        *b_smalltree_weight_PUdown;   //!
   TBranch        *b_smalltree_evtweight;   //!
   TBranch        *b_smalltree_evtweight_nominal;   //!
   TBranch        *b_smalltree_evtweight_minus;   //!
   TBranch        *b_smalltree_evtweight_plus;   //!
   TBranch        *b_smalltree_weight_toppt;   //!
   //TBranch        *b_smalltree_IChannel;   //!
   TBranch        *b_smalltree_nvertex;   //!
   TBranch        *b_smalltree_tmeme;   //!




   TreeReader(short int CorrOption, TTree *tree=0, TString sample="", short int systType = 0);
   virtual ~TreeReader();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(short int CorrOption, TString sample, TTree *tree );
   virtual void     Loop(short int CorrOption, vector<TString> datalist, vector<TString> mclist, TString sample, vector<TString> thesystlist, TString flavtag, short int systType);
   virtual Bool_t   Notify();

   void             initializeHisto(TString sample, TString syst, bool isfirstset, TString flavtag, short int systType);
   void             addHisto( TString var, TString selstep, TString sample, TString syst, int nbins, float min, float max, TString flavtag, short int systType);
   void             fillHisto(TString var, TString selstep, TString sample, TString syst, float val, float weight, TString flavtag, short int systType);
   void             manageOverflows();

   bool             applyEventSel(short int CorrOption, double SF_QCD_L, double SF_QCD_W, double SF_QCD_S, double SF_QCD_TT, double SFtrigger, double SFtriggerError, TString systtype, TString sample, TString flavtag, short int systType, double randomNumber);
   void             SetUpCSVreweighting();
   double           GetCSVweight(const int iSys, int jet_n, float *jet_pt,float *jet_eta,float *jet_btagdiscri,int *jet_flav);
   double           getQCDscalefactor(vector<TString> datalist, vector<TString> mclist, TString normregion);
   vector<double>   getY( TGraphAsymmErrors* graph, double Wvalue);
   vector<double>   getSFtrigger( TGraphAsymmErrors* ratioPlot, double pT, double eta);

   vector<TTree*> theTree_list;
   map<TString,TTree*> theTree_map;


   vector<TH1D*> histo_list;
   map<string, int> histo_map;
   int numb_histo;

   void deleteHisto();

   int   tree_SampleType;
   int   tree_Channel;

   void loadPDFSumWeight(TString channel, TString tmp_sample, TString flavtag);
   map<TString, float> PDFSumWeight;

   float tree_EvtWeight;
   float tree_met;
   float tree_wMass;
   float tree_deltaPhilb;
   float tree_wPt;

   BTagSFUtil *fBTagSF;

};

#endif

#ifdef TreeReader_cxx
TreeReader::TreeReader(short int CorrOption, TTree *tree, TString sample, short int systType) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   useElectronChannel = false;
   if (tree == 0) {
      TFile *f = 0;
      if(CorrOption == 0 && !useElectronChannel){
          f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p5.root");
      }else if(CorrOption == 0){
          f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p5_electron.root");
      }else if(CorrOption == 1 && !useElectronChannel){
          f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles/proof_merged_monotop_muon.root");
      }else if(CorrOption == 1){
          f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles/proof_merged_monotop_electron.root");
      }else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType == -1){
          if(!useElectronChannel) f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p6.root");
          else                    f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p6_electron.root");
      }else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType == 1){
          if(!useElectronChannel) f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p4.root");
          else                    f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p4_electron.root");
      }else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD"){
          if(!useElectronChannel) f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p5.root");
          else                    f = (TFile*)gROOT->GetListOfFiles()->FindObject("proof_QCDdatadriven_iso0p5_electron.root");
      }else if (CorrOption == 2 || CorrOption == 3){
                                 f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles/proof_merged_monotop_muon.root" );
          if(useElectronChannel) f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles/proof_merged_monotop_electron.root");
      }else cout << "ERROR: Wrong value of CorrOption! Allowed values: 0,1,2" << endl;

      if (!f || !f->IsOpen())
      {
         if(CorrOption == 0 && !useElectronChannel)                           f = new TFile("proof_QCDdatadriven_iso0p5.root");
         else if(CorrOption == 0)                                             f = new TFile("proof_QCDdatadriven_iso0p5_electron.root");
         else if(CorrOption == 1 && !useElectronChannel)                      f = new TFile("../InputFiles/proof_merged_monotop_muon.root");
         else if(CorrOption == 1)                                             f = new TFile("../InputFiles/proof_merged_monotop_electron.root");
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType == -1 && !useElectronChannel)
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p6.root");
         }
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType == -1)
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p6_electron.root");
         }
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType ==  1 && !useElectronChannel)
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p4.root");
         }
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && systType ==  1)
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p4_electron.root");
         }
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD" && !useElectronChannel)
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p5.root");
         }
         else if((CorrOption == 2 || CorrOption == 3) && sample == "QCD")
         {
                                                                              f = new TFile("proof_QCDdatadriven_iso0p5_electron.root");
         }
         else if( CorrOption == 2 || CorrOption == 3)
         {
                                                f = new TFile("../InputFiles/proof_merged_monotop_muon.root");
                         if(useElectronChannel) f = new TFile("../InputFiles/proof_merged_monotop_electron.root");
         }
         else cout << "ERROR: Wrong value of CorrOption! Allowed values: 0,1,2,3" << endl;
      }
   if     ((CorrOption == 2 || CorrOption == 3) && sample == "QCD")     f->GetObject( "SmallTree_QCDdatadriven" ,tree);
   else                                                                 f->GetObject(("SmallTree_"+sample).Data()             ,tree);
cout << "f_name= " << f->GetName() << endl;
   }

   Init(CorrOption,sample, tree);

   isData      = 1;
   isQCD       = 1;
   isTTbarSyst = 1;
   isSignal    = 1;
   if(useElectronChannel) thechannel = "eljets";
   else                   thechannel = "mujets";
}

TreeReader::~TreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeReader::Init(short int CorrOption, TString sample,  TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   if(sample == "QCD" || sample == "QCDdatadriven")
   {
       fChain->SetBranchAddress("smalltree_evtweight_nominal",  &smalltree_evtweight_nominal,  &b_smalltree_evtweight_nominal);
       fChain->SetBranchAddress("smalltree_evtweight_minus",    &smalltree_evtweight_minus,    &b_smalltree_evtweight_minus);
       fChain->SetBranchAddress("smalltree_evtweight_plus",     &smalltree_evtweight_plus,     &b_smalltree_evtweight_plus);
   }
   else fChain->SetBranchAddress("smalltree_evtweight",&smalltree_evtweight,        &b_smalltree_evtweight);

   fChain->SetBranchAddress("smalltree_nlepton",       &smalltree_nlepton,          &b_smalltree_nlepton);
   fChain->SetBranchAddress("smalltree_lept_pt",        smalltree_lept_pt,          &b_smalltree_lept_pt);
   fChain->SetBranchAddress("smalltree_lept_eta",       smalltree_lept_eta,         &b_smalltree_lept_eta);
   fChain->SetBranchAddress("smalltree_lept_phi",       smalltree_lept_phi,         &b_smalltree_lept_phi);
   fChain->SetBranchAddress("smalltree_lept_iso",       smalltree_lept_iso,         &b_smalltree_lept_iso);
   fChain->SetBranchAddress("smalltree_lept_flav",      smalltree_lept_flav,        &b_smalltree_lept_flav);
   fChain->SetBranchAddress("smalltree_njets",         &smalltree_njets,            &b_smalltree_njets);
   fChain->SetBranchAddress("smalltree_jet_pt",         smalltree_jet_pt,           &b_smalltree_jet_pt);
   fChain->SetBranchAddress("smalltree_jet_eta",        smalltree_jet_eta,          &b_smalltree_jet_eta);
   fChain->SetBranchAddress("smalltree_jet_phi",        smalltree_jet_phi,          &b_smalltree_jet_phi);
   fChain->SetBranchAddress("smalltree_jet_flav",       smalltree_jet_flav,         &b_smalltree_jet_flav);
   fChain->SetBranchAddress("smalltree_jet_btagdiscri", smalltree_jet_btagdiscri,   &b_smalltree_jet_btagdiscri);
   fChain->SetBranchAddress("smalltree_met_pt",        &smalltree_met_pt,           &b_smalltree_met_pt);
   fChain->SetBranchAddress("smalltree_met_phi",       &smalltree_met_phi,          &b_smalltree_met_phi);

   if(CorrOption != 0 && sample != "QCD")
   {

        fChain->SetBranchAddress("smalltree_jet_btagdiscri_up",     smalltree_jet_btagdiscri_up,        &b_smalltree_jet_btagdiscri_up);
        fChain->SetBranchAddress("smalltree_jet_btagdiscri_down",   smalltree_jet_btagdiscri_down,      &b_smalltree_jet_btagdiscri_down);
        fChain->SetBranchAddress("smalltree_jesup_njets",          &smalltree_jesup_njets,              &b_smalltree_jesup_njets);
        fChain->SetBranchAddress("smalltree_jet_jesup_pt",          smalltree_jet_jesup_pt,             &b_smalltree_jet_jesup_pt);
        fChain->SetBranchAddress("smalltree_jet_jesup_eta",         smalltree_jet_jesup_eta,            &b_smalltree_jet_jesup_eta);
        fChain->SetBranchAddress("smalltree_jet_jesup_phi",         smalltree_jet_jesup_phi,            &b_smalltree_jet_jesup_phi);
        fChain->SetBranchAddress("smalltree_jet_jesup_btagdiscri",  smalltree_jet_jesup_btagdiscri,     &b_smalltree_jet_jesup_btagdiscri);
        fChain->SetBranchAddress("smalltree_jet_jesup_flav",        smalltree_jet_jesup_flav,           &b_smalltree_jet_jesup_flav);
        fChain->SetBranchAddress("smalltree_jesdown_njets",        &smalltree_jesdown_njets,            &b_smalltree_jesdown_njets);
        fChain->SetBranchAddress("smalltree_jet_jesdown_pt",        smalltree_jet_jesdown_pt,           &b_smalltree_jet_jesdown_pt);
        fChain->SetBranchAddress("smalltree_jet_jesdown_eta",       smalltree_jet_jesdown_eta,          &b_smalltree_jet_jesdown_eta);
        fChain->SetBranchAddress("smalltree_jet_jesdown_phi",       smalltree_jet_jesdown_phi,          &b_smalltree_jet_jesdown_phi);
        fChain->SetBranchAddress("smalltree_jet_jesdown_btagdiscri",smalltree_jet_jesdown_btagdiscri,   &b_smalltree_jet_jesdown_btagdiscri);
        fChain->SetBranchAddress("smalltree_jet_jesdown_flav",      smalltree_jet_jesdown_flav,         &b_smalltree_jet_jesdown_flav);
        fChain->SetBranchAddress("smalltree_jerup_njets",          &smalltree_jerup_njets,              &b_smalltree_jerup_njets);
        fChain->SetBranchAddress("smalltree_jet_jerup_pt",          smalltree_jet_jerup_pt,             &b_smalltree_jet_jerup_pt);
        fChain->SetBranchAddress("smalltree_jet_jerup_eta",         smalltree_jet_jerup_eta,            &b_smalltree_jet_jerup_eta);
        fChain->SetBranchAddress("smalltree_jet_jerup_phi",         smalltree_jet_jerup_phi,            &b_smalltree_jet_jerup_phi);
        fChain->SetBranchAddress("smalltree_jet_jerup_btagdiscri",  smalltree_jet_jerup_btagdiscri,     &b_smalltree_jet_jerup_btagdiscri);
        fChain->SetBranchAddress("smalltree_jet_jerup_flav",        smalltree_jet_jerup_flav,           &b_smalltree_jet_jerup_flav);
        fChain->SetBranchAddress("smalltree_jerdown_njets",        &smalltree_jerdown_njets,            &b_smalltree_jerdown_njets);
        fChain->SetBranchAddress("smalltree_jet_jerdown_pt",        smalltree_jet_jerdown_pt,           &b_smalltree_jet_jerdown_pt);
        fChain->SetBranchAddress("smalltree_jet_jerdown_eta",       smalltree_jet_jerdown_eta,          &b_smalltree_jet_jerdown_eta);
        fChain->SetBranchAddress("smalltree_jet_jerdown_phi",       smalltree_jet_jerdown_phi,          &b_smalltree_jet_jerdown_phi);
        fChain->SetBranchAddress("smalltree_jet_jerdown_btagdiscri",smalltree_jet_jerdown_btagdiscri,   &b_smalltree_jet_jerdown_btagdiscri);
        fChain->SetBranchAddress("smalltree_jet_jerdown_flav",      smalltree_jet_jerdown_flav,         &b_smalltree_jet_jerdown_flav);
        fChain->SetBranchAddress("smalltree_met_jesup_pt",         &smalltree_met_jesup_pt,             &b_smalltree_met_jesup_pt);
        fChain->SetBranchAddress("smalltree_met_jesup_phi",        &smalltree_met_jesup_phi,            &b_smalltree_met_jesup_phi);
        fChain->SetBranchAddress("smalltree_met_jesdown_pt",       &smalltree_met_jesdown_pt,           &b_smalltree_met_jesdown_pt);
        fChain->SetBranchAddress("smalltree_met_jesdown_phi",      &smalltree_met_jesdown_phi,          &b_smalltree_met_jesdown_phi);
        fChain->SetBranchAddress("smalltree_met_jerup_pt",         &smalltree_met_jerup_pt,             &b_smalltree_met_jerup_pt);
        fChain->SetBranchAddress("smalltree_met_jerup_phi",        &smalltree_met_jerup_phi,            &b_smalltree_met_jerup_phi);
        fChain->SetBranchAddress("smalltree_met_jerdown_pt",       &smalltree_met_jerdown_pt,           &b_smalltree_met_jerdown_pt);
        fChain->SetBranchAddress("smalltree_met_jerdown_phi",      &smalltree_met_jerdown_phi,          &b_smalltree_met_jerdown_phi);
        fChain->SetBranchAddress("smalltree_met_unclsup_pt",       &smalltree_met_unclsup_pt,           &b_smalltree_met_unclsup_pt);
        fChain->SetBranchAddress("smalltree_met_unclsup_phi",      &smalltree_met_unclsup_phi,          &b_smalltree_met_unclsup_phi);
        fChain->SetBranchAddress("smalltree_met_unclsdown_pt",     &smalltree_met_unclsdown_pt,         &b_smalltree_met_unclsdown_pt);
        fChain->SetBranchAddress("smalltree_met_unclsdown_phi",    &smalltree_met_unclsdown_phi,        &b_smalltree_met_unclsdown_phi);
        fChain->SetBranchAddress("smalltree_weight_trigup",        &smalltree_weight_trigup,            &b_smalltree_weight_trigup);
        fChain->SetBranchAddress("smalltree_weight_trigdown",      &smalltree_weight_trigdown,          &b_smalltree_weight_trigdown);
        fChain->SetBranchAddress("smalltree_weight_leptup",        &smalltree_weight_leptup,            &b_smalltree_weight_leptup);
        fChain->SetBranchAddress("smalltree_weight_leptdown",      &smalltree_weight_leptdown,          &b_smalltree_weight_leptdown);
        fChain->SetBranchAddress("smalltree_weight_PDFup",         &smalltree_weight_PDFup,             &b_smalltree_weight_PDFup);
        fChain->SetBranchAddress("smalltree_weight_PDFdown",       &smalltree_weight_PDFdown,           &b_smalltree_weight_PDFdown);
        fChain->SetBranchAddress("smalltree_weight_PUup",          &smalltree_weight_PUup,              &b_smalltree_weight_PUup);
        fChain->SetBranchAddress("smalltree_weight_PUdown",        &smalltree_weight_PUdown,            &b_smalltree_weight_PUdown);
        fChain->SetBranchAddress("smalltree_weight_toppt",         &smalltree_weight_toppt,             &b_smalltree_weight_toppt);
        fChain->SetBranchAddress("smalltree_metnosmear_pt",        &smalltree_metnosmear_pt,            &b_smalltree_metnosmear_pt);
        fChain->SetBranchAddress("smalltree_metnosmear_phi",       &smalltree_metnosmear_phi,           &b_smalltree_metnosmear_phi);
    //   fChain->SetBranchAddress("smalltree_IChannel",            &smalltree_IChannel,                 &b_smalltree_IChannel);
        fChain->SetBranchAddress("smalltree_nvertex",              &smalltree_nvertex,                  &b_smalltree_nvertex);
   }

   if(CorrOption != 0 && (sample == "TTMSDecays_central" || sample == "TTMSDecays_matchingup" || sample == "TTMSDecays_matchingdown" || sample == "TTMSDecays_scaledown" ||  sample == "TTMSDecays_scaleup" ))
   {
        fChain->SetBranchAddress("smalltree_tmeme",                &smalltree_tmeme,                    &b_smalltree_tmeme);
   }

   Notify();
}

Bool_t TreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TreeReader_cxx
