#include <iostream>
#include <memory>
#include <random>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TTree.h"

void ReadData() {

  // 4-vectors are stored with the convention
  // p = (px, py, pz, E)
  float photon1[4], photon2[4];

  auto inFile = TFile::Open("data/data_high.root");
  auto tree = (TTree *)inFile->Get("Data");

  tree->SetBranchAddress("photon1", &photon1);
  tree->SetBranchAddress("photon2", &photon2);

  for (Long64_t iEv = 0; iEv < tree->GetEntries(); iEv++) {
    tree->GetEntry(iEv);

    // do stuff here...
  }
}
