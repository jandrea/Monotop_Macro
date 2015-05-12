#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TFile.h>
#include <TString.h>

#include "TLorentzVector.h"

#include <iostream>
using namespace std;

class TreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

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
   Float_t         smalltree_met_pt;
   Float_t         smalltree_met_phi;
   Float_t         smalltree_evtweight;

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
   TBranch        *b_smalltree_met_pt;   //!
   TBranch        *b_smalltree_met_phi;   //!
   TBranch        *b_smalltree_evtweight;   //!


   TreeReader(vector<TString> samplelist, TTree* tree, TString sample, bool useElectronChannel);
   virtual void     Init(TString sample, TTree* tree);

   //void             initializeHisto(TString sample);
   //void             addHisto( TString var, TString selstep, TString sample, TString syst, int nbins, float min, float max, TString flavtag);
   //void             fillHisto(TString channel, TString var, TString selstep, TString sample, TString syst, float val, float weight, TString flavtag);

};

TreeReader::TreeReader(vector<TString> samplelist, TTree* tree, TString sample, bool useElectronChannel): fChain(0)
{

   TFile *f = 0;
   if (tree == 0)
   {
          if(!useElectronChannel) f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles_IsoSup0p4/proof_IsoSup0p4_merged.root");
          else                    f = (TFile*)gROOT->GetListOfFiles()->FindObject("../InputFiles_IsoSup0p4/proof_IsoSup0p4_monotop_electron.root");

        if (!f || !f->IsOpen())
        {
            if(!useElectronChannel) f = new TFile("../InputFiles_IsoSup0p4/proof_IsoSup0p4_merged.root","READ");
            else                    f = new TFile("../InputFiles_IsoSup0p4/proof_IsoSup0p4_monotop_electron.root","READ");
        }
   }
   f->GetObject(("SmallTree_"+sample).Data()             ,tree);

   Init(sample, tree);
}

void TreeReader::Init(TString sample,  TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("smalltree_evtweight",     &smalltree_evtweight,        &b_smalltree_evtweight);
   fChain->SetBranchAddress("smalltree_nlepton",       &smalltree_nlepton,          &b_smalltree_nlepton);
   fChain->SetBranchAddress("smalltree_lept_pt",       &smalltree_lept_pt,          &b_smalltree_lept_pt);
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

}

void createQCDfromData()
{

    double InvertedIso = 0.6;
    bool useElectronChannel = true;

    std::vector<TString> samplelist;
    if(!useElectronChannel)
    {
        samplelist.push_back("NTuple_53X_SingleMuRun2012A"   );
        samplelist.push_back("NTuple_53X_SingleMuRun2012B"   );
        samplelist.push_back("NTuple_53X_SingleMuRun2012C"   );
        samplelist.push_back("NTuple_53X_SingleMuRun2012D"   );
        samplelist.push_back("NTuple_53X_TTJetsMadgraphZ2"   );
        samplelist.push_back("NTuple_53X_WJetsToLNu"         );
        samplelist.push_back("NTuple_53X_DYJetsToLL_M-10To50");
        samplelist.push_back("NTuple_53X_DYJetsToLL_M-50"    );
        samplelist.push_back("NTuple_53X_T_s-channel"        );
        samplelist.push_back("NTuple_53X_T_t-channel"        );
        samplelist.push_back("NTuple_53X_T_tW-channel"       );
        samplelist.push_back("NTuple_53X_Tbar_t-channel"     );
        samplelist.push_back("NTuple_53X_Tbar_tW-channel"    );
        samplelist.push_back("NTuple_53X_WZJetsIncl"         );
        samplelist.push_back("NTuple_53X_WWJetsIncl"         );
        samplelist.push_back("NTuple_53X_ZZJetsIncl"         );
    }
    else
    {
        samplelist.push_back("SingleElA"            );
        samplelist.push_back("SingleElB"            );
        samplelist.push_back("SingleElC"            );
        samplelist.push_back("SingleElD"            );
        samplelist.push_back("TTMSDecays_central"   );
        samplelist.push_back("W0Jets"               );
        samplelist.push_back("W1Jets"               );
        samplelist.push_back("W2Jets"               );
        samplelist.push_back("W3Jets"               );
        samplelist.push_back("W4Jets"               );
        samplelist.push_back("DYJetsToLL_M-10To50"  );
        samplelist.push_back("DYJetsToLL_M-50"      );
        samplelist.push_back("T_s"                  );
        samplelist.push_back("T_t"                  );
        samplelist.push_back("T_tW"                 );
        samplelist.push_back("Tbar_t"               );
        samplelist.push_back("Tbar_tW"              );
        samplelist.push_back("WZ"                   );
        samplelist.push_back("WW"                   );
        samplelist.push_back("ZZ"                   );
    }



    TFile * theoutputfile;
    if(!useElectronChannel)
    {
        if(InvertedIso == 0.4)      theoutputfile= new TFile( "proof_QCDdatadriven_iso0p4.root", "recreate");
        else if(InvertedIso == 0.5) theoutputfile= new TFile( "proof_QCDdatadriven_iso0p5.root", "recreate");
        else if(InvertedIso == 0.6) theoutputfile= new TFile( "proof_QCDdatadriven_iso0p6.root", "recreate");
    }
    else
    {
        if(InvertedIso == 0.4)      theoutputfile= new TFile( "proof_QCDdatadriven_iso0p4_electron.root", "recreate");
        else if(InvertedIso == 0.5) theoutputfile= new TFile( "proof_QCDdatadriven_iso0p5_electron.root", "recreate");
        else if(InvertedIso == 0.6) theoutputfile= new TFile( "proof_QCDdatadriven_iso0p6_electron.root", "recreate");
    }
    TTree* tree = 0;

    // Declaration of leaf types
    Int_t          tree_nlepton;
    Float_t        tree_lept_pt[3];   //[smalltree_nlepton]
    Float_t        tree_lept_eta[3];   //[smalltree_nlepton]
    Float_t        tree_lept_phi[3];   //[smalltree_nlepton]
    Float_t        tree_lept_iso[3];   //[smalltree_nlepton]
    Int_t          tree_lept_flav[3];   //[smalltree_nlepton]
    Int_t          tree_njets;
    Float_t        tree_jet_pt[100];   //[smalltree_njets]
    Float_t        tree_jet_eta[100];   //[smalltree_njets]
    Float_t        tree_jet_phi[100];   //[smalltree_njets]
    Float_t        tree_jet_btagdiscri[100];   //[smalltree_njets]
    Float_t        tree_jet_btagdiscri_up[100];   //[smalltree_njets]
    Float_t        tree_jet_btagdiscri_down[100];   //[smalltree_njets]
    Int_t          tree_jet_flav[100];   //[smalltree_njets]
    Float_t        tree_met_pt;
    Float_t        tree_met_phi;
    Float_t        tree_evtweight_nominal;
    Float_t        tree_evtweight_minus;
    Float_t        tree_evtweight_plus;

    TTree* outputTree = new TTree("SmallTree_QCDdatadriven", "SmallTree_QCDdatadriven");
    outputTree->Branch("smalltree_evtweight_nominal",   &tree_evtweight_nominal,"smalltree_evtweight_nominal" );
    outputTree->Branch("smalltree_evtweight_minus",     &tree_evtweight_minus,  "smalltree_evtweight_minus"   );
    outputTree->Branch("smalltree_evtweight_plus",      &tree_evtweight_plus,   "smalltree_evtweight_plus"    );
    outputTree->Branch("smalltree_nlepton",             &tree_nlepton,          "smalltree_nlepton/I"         );
    outputTree->Branch("smalltree_lept_pt",             &tree_lept_pt,          "smalltree_lept_pt"       );
    outputTree->Branch("smalltree_lept_eta",            &tree_lept_eta,         "smalltree_lept_eta"      );
    outputTree->Branch("smalltree_lept_phi",            &tree_lept_phi,         "smalltree_lept_phi"      );
    outputTree->Branch("smalltree_lept_iso",            &tree_lept_iso,         "smalltree_lept_iso"      );
    outputTree->Branch("smalltree_lept_flav",           &tree_lept_flav,        "smalltree_lept_flav/I"     );
    outputTree->Branch("smalltree_njets",               &tree_njets,            "smalltree_njets/I"         );
    outputTree->Branch("smalltree_jet_pt",              &tree_jet_pt,           "smalltree_jet_pt"        );
    outputTree->Branch("smalltree_jet_eta",             &tree_jet_eta,          "smalltree_jet_eta"       );
    outputTree->Branch("smalltree_jet_phi",             &tree_jet_phi,          "smalltree_jet_phi"       );
    outputTree->Branch("smalltree_jet_flav",            &tree_jet_flav,         "smalltree_jet_flav/I"      );
    outputTree->Branch("smalltree_jet_btagdiscri",      &tree_jet_btagdiscri,   "smalltree_jet_btagdiscri");
    outputTree->Branch("smalltree_met_pt",              &tree_met_pt,           "smalltree_met_pt"        );
    outputTree->Branch("smalltree_met_phi",             &tree_met_phi,          "smalltree_met_phi"       );

    for (unsigned int i = 0; i < samplelist.size(); i++)
    {
        bool isData = false;
        if(samplelist[i] == "NTuple_53X_SingleMuRun2012A" || samplelist[i] == "NTuple_53X_SingleMuRun2012B" || samplelist[i] == "NTuple_53X_SingleMuRun2012C" || samplelist[i] == "NTuple_53X_SingleMuRun2012D") isData = true;
        if(samplelist[i] == "SingleElA" || samplelist[i] == "SingleElB" || samplelist[i] == "SingleElC" || samplelist[i] == "SingleElD") isData = true;

        TreeReader * tree_ = new TreeReader(samplelist,  tree, samplelist[i], useElectronChannel);

        if (tree_->fChain == 0) return;

        cout << endl;
        cout << "Processing the events... " << endl;
        cout << endl;

        Long64_t nentries = tree_->fChain->GetEntriesFast();
        //Long64_t nentries = tree->GetEntriesFast();

        Long64_t nbytes = 0, nb = 0;
        for (Long64_t jentry=0; jentry<nentries;jentry++)
        {
            Long64_t ientry = tree_->fChain->LoadTree(jentry);
            if (ientry < 0) break;
            nb = tree_->fChain->GetEntry(jentry);   nbytes += nb;
            //nb = tree->GetEntry(jentry);   nbytes += nb;
           // cout << "pT(muon)= " << tree_->smalltree_lept_pt[0] << endl;
           // cout << "MET     = " << tree_->smalltree_met_pt << endl;
           // cout << "Nlepton     = " << tree_->smalltree_nlepton << endl;
            if(isData)
            {
                tree_evtweight_nominal =   tree_->smalltree_evtweight;
                tree_evtweight_minus   =   tree_->smalltree_evtweight;
                tree_evtweight_plus    =   tree_->smalltree_evtweight;
            }
            else
            {
                tree_evtweight_nominal = - tree_->smalltree_evtweight;
                tree_evtweight_minus   = - 2*tree_->smalltree_evtweight;
                tree_evtweight_plus    = 0;
            }
            //if (tree_evtweight < 0) cout << "Evt_weight_bgkd= " << tree_evtweight << endl;
            tree_nlepton        = tree_->smalltree_nlepton;
            //cout << "Evtweight_ex   = " << tree_->smalltree_evtweight << endl;
            //cout << "Evtweight_new  = " << tree_evtweight << endl;
            //cout << "NLepton_ex     = " << tree_->smalltree_nlepton << endl;
            //cout << "NLepton_new    = " << tree_nlepton << endl;
            for(unsigned short int lep = 0; lep < tree_->smalltree_nlepton; lep++)
            {
                tree_lept_pt[lep]   = tree_->smalltree_lept_pt[lep];
                tree_lept_eta[lep]  = tree_->smalltree_lept_eta[lep];
                tree_lept_phi[lep]  = tree_->smalltree_lept_phi[lep];
                tree_lept_iso[lep]  = tree_->smalltree_lept_iso[lep];
                tree_lept_flav[lep] = tree_->smalltree_lept_flav[lep];
            }
            tree_njets          = tree_->smalltree_njets;
            //cout << "NJets_ex   = " << tree_->smalltree_njets << endl;
            //cout << "NJets_new  = " << tree_njets << endl;

            for(unsigned short int jet = 0; jet < tree_->smalltree_njets; jet++)
            {
                tree_jet_pt[jet]         = tree_->smalltree_jet_pt[jet];
                tree_jet_eta[jet]        = tree_->smalltree_jet_eta[jet];
                tree_jet_phi[jet]        = tree_->smalltree_jet_phi[jet];
                tree_jet_flav[jet]       = tree_->smalltree_jet_flav[jet];
                tree_jet_btagdiscri[jet] = tree_->smalltree_jet_btagdiscri[jet];
            }

            tree_met_pt         = tree_->smalltree_met_pt;
            tree_met_phi        = tree_->smalltree_met_phi;

            if(tree_lept_iso[0] < InvertedIso) continue;
            outputTree->Fill();
        }

        delete tree_;
    }

    theoutputfile->cd();
    outputTree->Write();
    theoutputfile->Write();
    theoutputfile->Close();
    delete theoutputfile;

}
