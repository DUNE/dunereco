//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 15 08:40:44 2016 by ROOT version 5.34/32
// from TTree CVNTrainTree/Training records
// found on file: cvn_event_dump_r20000001_s101_hist.root
//////////////////////////////////////////////////////////

#ifndef Analyze_h
#define Analyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Analyze {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //cvn::TrainingData *train;
   Int_t           fInt;
   Float_t         fNuEnergy;
   Float_t         fLepEnergy;
   UInt_t          fPMap_fNWire;
   UInt_t          fPMap_fNTdc;
   vector<float>   fPMap_fPE;
   vector<float>   fPMap_fPEX;
   vector<float>   fPMap_fPEY;
   vector<float>   fPMap_fPEZ;
   vector<double>  fPMap_fPur;
   vector<double>  fPMap_fPurX;
   vector<double>  fPMap_fPurY;
   vector<double>  fPMap_fPurZ;
 //vector<cvn::HType> fPMap_fLab;
 //vector<cvn::HType> fPMap_fLabX;
 //vector<cvn::HType> fPMap_fLabY;
 //vector<cvn::HType> fPMap_fLabZ;
   Int_t           fPMap_fBound_fFirstWire[3];
   Int_t           fPMap_fBound_fLastWire[3];
   Double_t        fPMap_fBound_fFirstTDC[3];
   Double_t        fPMap_fBound_fLastTDC[3];

   // List of branches
   TBranch        *b_train_fInt;   //!
   TBranch        *b_train_fNuEnergy;   //!
   TBranch        *b_train_fLepEnergy;   //!
   TBranch        *b_train_fPMap_fNWire;   //!
   TBranch        *b_train_fPMap_fNTdc;   //!
   TBranch        *b_train_fPMap_fPE;   //!
   TBranch        *b_train_fPMap_fPEX;   //!
   TBranch        *b_train_fPMap_fPEY;   //!
   TBranch        *b_train_fPMap_fPEZ;   //!
   TBranch        *b_train_fPMap_fPur;   //!
   TBranch        *b_train_fPMap_fPurX;   //!
   TBranch        *b_train_fPMap_fPurY;   //!
   TBranch        *b_train_fPMap_fPurZ;   //!
   TBranch        *b_train_fPMap_fBound_fFirstWire;   //!
   TBranch        *b_train_fPMap_fBound_fLastWire;   //!
   TBranch        *b_train_fPMap_fBound_fFirstTDC;   //!
   TBranch        *b_train_fPMap_fBound_fLastTDC;   //!

   Analyze(TTree *tree=0);
   virtual ~Analyze();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyze_cxx
Analyze::Analyze(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("cvn_event_dump_r20000001_s101_hist.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("cvn_event_dump_r20000001_s101_hist.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("cvn_event_dump_r20000001_s101_hist.root:/cvndump");
      dir->GetObject("CVNTrainTree",tree);

   }
   Init(tree);
}

Analyze::~Analyze()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze::LoadTree(Long64_t entry)
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

void Analyze::Init(TTree *tree)
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

   fChain->SetBranchAddress("fInt", &fInt, &b_train_fInt);
   fChain->SetBranchAddress("fNuEnergy", &fNuEnergy, &b_train_fNuEnergy);
   fChain->SetBranchAddress("fLepEnergy", &fLepEnergy, &b_train_fLepEnergy);
   fChain->SetBranchAddress("fPMap.fNWire", &fPMap_fNWire, &b_train_fPMap_fNWire);
   fChain->SetBranchAddress("fPMap.fNTdc", &fPMap_fNTdc, &b_train_fPMap_fNTdc);
   fChain->SetBranchAddress("fPMap.fPE", &fPMap_fPE, &b_train_fPMap_fPE);
   fChain->SetBranchAddress("fPMap.fPEX", &fPMap_fPEX, &b_train_fPMap_fPEX);
   fChain->SetBranchAddress("fPMap.fPEY", &fPMap_fPEY, &b_train_fPMap_fPEY);
   fChain->SetBranchAddress("fPMap.fPEZ", &fPMap_fPEZ, &b_train_fPMap_fPEZ);
   fChain->SetBranchAddress("fPMap.fPur", &fPMap_fPur, &b_train_fPMap_fPur);
   fChain->SetBranchAddress("fPMap.fPurX", &fPMap_fPurX, &b_train_fPMap_fPurX);
   fChain->SetBranchAddress("fPMap.fPurY", &fPMap_fPurY, &b_train_fPMap_fPurY);
   fChain->SetBranchAddress("fPMap.fPurZ", &fPMap_fPurZ, &b_train_fPMap_fPurZ);
   fChain->SetBranchAddress("fPMap.fBound.fFirstWire[3]", fPMap_fBound_fFirstWire, &b_train_fPMap_fBound_fFirstWire);
   fChain->SetBranchAddress("fPMap.fBound.fLastWire[3]", fPMap_fBound_fLastWire, &b_train_fPMap_fBound_fLastWire);
   fChain->SetBranchAddress("fPMap.fBound.fFirstTDC[3]", fPMap_fBound_fFirstTDC, &b_train_fPMap_fBound_fFirstTDC);
   fChain->SetBranchAddress("fPMap.fBound.fLastTDC[3]", fPMap_fBound_fLastTDC, &b_train_fPMap_fBound_fLastTDC);
   Notify();
}

Bool_t Analyze::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyze::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyze_cxx
