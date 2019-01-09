#include <iostream>
#include <memory>
#include <random>

#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"

// clang-format off
#ifdef __CLING__
#pragma link C++ class vector<TLorentzVector>+;
#endif
// clang-format on

double AngRes(double E) {
  // in mrad
  double t1 = 36.5 / sqrt(E);
  double t2 = 4.1;

  // in rad
  return 1e-3 * sqrt(t1 * t1 + t2 * t2);
}

void MCGeneration() {
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

  // gRandom = new TRandom3(42);
  gRandom = new TRandom3(0);

  float HMass = 750.;
  float sqrts = 1300.;
  float pMass = 0.938;
  float phRes = 0.02;
  float jRes = 0.3;
  int nJetMean = 2;

  int nBkg = 1000000;
  // int nSig = gRandom->Poisson(200);
  // int nSig = gRandom->Poisson(10000);

  cout << nBkg << " background events" << endl;
  // cout << nSig << " signal events" << endl;

  TLorentzVector p1, p2;
  p1.SetXYZM(0, 0, sqrts / 2, pMass);
  p2.SetXYZM(0, 0, -sqrts / 2, pMass);

  TLorentzVector W = p1 + p2;

  auto evGen = std::unique_ptr<TGenPhaseSpace>(new TGenPhaseSpace());

  std::vector<TLorentzVector> jets;
  std::vector<TLorentzVector> photons;

  auto outFile = std::unique_ptr<TFile>(new TFile("data.root", "recreate"));
  auto tree = new TTree("temp", "temp");
  tree->SetDirectory(0);
  tree->Branch("photons", &photons);
  tree->Branch("jets", &jets);

  TH1D *hCheck =
      new TH1D("hCheck", ";M_{#gamma#gamma} (GeV);Counts)", 200, 10, 1301);
  TH1D *hCheckSel =
      new TH1D("hCheckSel", ";M_{#gamma#gamma} (GeV);Counts)", 200, 10, 1301);

  TH1D *hAngle =
      new TH1D("hAngle", ";#theta_{#gamma#gamma};Counts)", 200, 0, TMath::Pi());
  TH1D *hAngleCheckBkg = new TH1D(
      "hAngleCheckBkg", ";#theta_{#gamma#gamma};Counts)", 200, 0, TMath::Pi());
  TH1D *hAngleCheckSig = new TH1D(
      "hAngleCheckSig", ";#theta_{#gamma#gamma};Counts)", 200, 0, TMath::Pi());

  TH2D *hAngleMassSig =
      new TH2D("hAngleMassSig", ";#theta_{#gamma#gamma};M_{#gamma#gamma}", 200,
               0, TMath::Pi(), 200, 10, 1301);
  TH2D *hAngleMassBkg =
      new TH2D("hAngleMassBkg", ";#theta_{#gamma#gamma};M_{#gamma#gamma}", 200,
               0, TMath::Pi(), 200, 10, 1301);

  for (int iEv = 0; iEv < nBkg; iEv++) { // event body: bkg events
    jets.clear();
    photons.clear();

    int nJets = gRandom->Poisson(nJetMean);
    if (nJets == 0) {
      iEv--;
      continue;
    }
    jets.resize(nJets);

    const int nP = 2 + nJets;
    double masses[nP];
    masses[0] = masses[1] = 0.;
    for (int ij = 0; ij < nJets; ij++)
      masses[2 + ij] = 3.;

    evGen->SetDecay(W, nP, masses);
    evGen->Generate();

    photons.push_back(*(evGen->GetDecay(0)));
    photons.push_back(*(evGen->GetDecay(1)));

    for (int ij = 0; ij < nJets; ij++)
      jets.push_back(*(evGen->GetDecay(2 + ij)));

    for (auto &photon : photons) {
      photon.SetTheta(photon.Theta() * gRandom->Gaus(1, AngRes(photon.E())));
      photon.SetE(photon.E() * gRandom->Gaus(1, phRes));
    }
    for (auto &jet : jets) {
      jet.SetE(jet.E() * gRandom->Gaus(1, jRes));
    }

    tree->Fill();

    hCheck->Fill((photons.at(0) + photons.at(1)).M());
    hAngle->Fill(photons.at(0).Angle(photons.at(1).Vect()));
    hAngleCheckBkg->Fill(photons.at(0).Angle(photons.at(1).Vect()));

    hAngleMassBkg->Fill(photons.at(0).Angle(photons.at(1).Vect()),
                        (photons.at(0) + photons.at(1)).M());

    if (photons.at(0).Angle(photons.at(1).Vect()) > 2) {
      hCheckSel->Fill((photons.at(0) + photons.at(1)).M());
    }
  }

  // for (int iEv = 0; iEv < nSig; iEv++) { // event body: sig events
  //   jets.clear();
  //   photons.clear();
  //
  //   int nJets = gRandom->Poisson(nJetMean);
  //   if (nJets == 0) {
  //     iEv--;
  //     continue;
  //   }
  //   jets.resize(nJets);
  //
  //   const int nP = 1 + nJets;
  //   double masses[nP];
  //   masses[0] = HMass;
  //   for (int ij = 0; ij < nJets; ij++)
  //     masses[1 + ij] = 3.;
  //
  //   evGen->SetDecay(W, nP, masses);
  //   evGen->Generate();
  //
  //   for (int ij = 0; ij < nJets; ij++)
  //     jets.push_back(*(evGen->GetDecay(1 + ij)));
  //
  //   double dummy_mass[2] = {0., 0.};
  //   evGen->SetDecay(*(evGen->GetDecay(0)), 2, dummy_mass);
  //   evGen->Generate();
  //
  //   photons.push_back(*(evGen->GetDecay(0)));
  //   photons.push_back(*(evGen->GetDecay(1)));
  //
  //   for (auto &photon : photons) {
  //     photon.SetTheta(photon.Theta() * gRandom->Gaus(1, AngRes(photon.E())));
  //     photon.SetE(photon.E() * gRandom->Gaus(1, phRes));
  //   }
  //   for (auto &jet : jets) {
  //     jet.SetE(jet.E() * gRandom->Gaus(1, jRes));
  //   }
  //
  //   tree->Fill();
  //
  //   hCheck->Fill((photons.at(0) + photons.at(1)).M());
  //   hAngle->Fill(photons.at(0).Angle(photons.at(1).Vect()));
  //   hAngleCheckSig->Fill(photons.at(0).Angle(photons.at(1).Vect()));
  //
  //   hAngleMassSig->Fill(photons.at(0).Angle(photons.at(1).Vect()),
  //                       (photons.at(0) + photons.at(1)).M());
  //
  //   if (photons.at(0).Angle(photons.at(1).Vect()) > 2) {
  //     hCheckSel->Fill((photons.at(0) + photons.at(1)).M());
  //   }
  // }

  //---------------------------
  // Plotting mass no selection
  //---------------------------
  hCheck->SetMarkerStyle(20);
  hCheck->SetLineColor(1);
  auto cc = new TCanvas("cc", "", 1200 + 4, 600 + 28);
  hCheck->Draw("E");

  auto fB = new TF1("fB", "[0]*exp([1]*x)", 150, 1250);
  fB->SetParameters(100, -3);
  auto fS = new TF1("fS", "gaus", 150, 1250);
  fS->SetParameters(1, 750, 13);

  auto fTot = new TF1("fTot", "fB+fS", 150, 1250);
  fTot->SetParLimits(3, 200, 1200);
  // fTot->SetParLimits(4, -100, 100);
  fTot->FixParameter(4, 13.4);

  hCheck->Fit(fTot, "LRN");

  fB->SetNpx(1000);
  fS->SetNpx(1000);
  fTot->SetNpx(1000);

  fTot->SetLineColor(kBlue);
  fTot->SetLineWidth(2);

  fB->SetLineColor(kRed);
  fS->SetLineColor(kGreen + 2);

  for (int ip = 0; ip < fB->GetNpar(); ip++) {
    fB->SetParameter(ip, fTot->GetParameter(ip));
  }
  for (int ip = 0; ip < fS->GetNpar(); ip++) {
    fS->SetParameter(ip, fTot->GetParameter(ip + fB->GetNpar()));
  }

  float NSel = fS->GetParameter(0) * fabs(fS->GetParameter(2)) *
               sqrt(2 * TMath::Pi()) / hCheck->GetBinWidth(0);
  cout << "Reconstructed " << NSel << " events" << endl;

  fS->Draw("same");
  fB->Draw("same");
  fTot->Draw("same");
  cc->Print("controlplot1.pdf");

  int startbin = hCheck->FindBin(700);
  int endbin = hCheck->FindBin(800);
  int ntot = 0;
  for(int ibin=startbin; ibin<endbin; ibin++){
    ntot += hCheck->GetBinContent(ibin);
  }
  // std::cout << "Bkg = " << ntot - nSig << std::endl;

  //---------------------------
  // Plotting angle
  //---------------------------
  hAngle->SetMarkerStyle(20);
  hAngle->SetLineColor(1);
  hAngle->SetMinimum(0.1);
  auto c2 = new TCanvas("c2", "", 1200 + 4, 600 + 28);
  hAngleCheckBkg->SetLineColor(2);
  hAngleCheckSig->SetLineColor(4);

  hAngle->Draw("AXIS");
  hAngleCheckSig->Draw("same");
  hAngleCheckBkg->Draw("same");
  hAngle->Draw("E same");
  c2->Print("angle.pdf");

  //---------------------------
  // Plotting mass w selection
  //---------------------------
  hCheckSel->SetMarkerStyle(20);
  hCheckSel->SetLineColor(1);
  auto c3 = new TCanvas("c3", "", 1200 + 4, 600 + 28);
  hCheckSel->Draw("E");

  auto fB2 = new TF1("fB2", "gaus", 150, 1250);
  auto fS2 = new TF1("fS2", "gaus", 150, 1250);
  auto fTot2 = new TF1("fTot2", "gaus(0) + gaus(3)", 150, 1250);
  fTot2->SetParameters(100, 200, 100, 1, 750, 13);
  fTot2->SetParLimits(4, 600, 800);
  // fTot2->SetParLimits(5, -100, 100);
  fTot2->FixParameter(5, 13.4);

  hCheckSel->Fit(fTot2, "LRN");

  fB2->SetNpx(1000);
  fS2->SetNpx(1000);
  fTot2->SetNpx(1000);

  fTot2->SetLineColor(kBlue);
  fTot2->SetLineWidth(2);

  fB2->SetLineColor(kRed);
  fS2->SetLineColor(kGreen + 2);

  for (int ip = 0; ip < fB2->GetNpar(); ip++) {
    fB2->SetParameter(ip, fTot2->GetParameter(ip));
  }
  for (int ip = 0; ip < fS2->GetNpar(); ip++) {
    fS2->SetParameter(ip, fTot2->GetParameter(ip + fB2->GetNpar()));
    fS2->SetParError(ip, fTot2->GetParError(ip + fB2->GetNpar()));
  }

  NSel = fS2->GetParameter(0) * fabs(fS2->GetParameter(2)) *
         sqrt(2 * TMath::Pi()) / hCheckSel->GetBinWidth(0);
  float NErr = pow(fS2->GetParError(0) * fabs(fS2->GetParameter(2)), 2);
  NErr += pow(fabs(fS2->GetParameter(0)) * fS2->GetParError(2), 2);
  NErr = sqrt(NErr);
  cout << "Reconstructed " << NSel << " +- " << NErr << " events" << endl;

  fS2->Draw("same");
  fB2->Draw("same");
  fTot2->Draw("same");
  c3->Print("controlplot2.pdf");

  auto c4 = new TCanvas("c4", "", 1200 + 4, 600 + 28);
  c4->Divide(2, 1);
  c4->cd(1);
  hAngleMassSig->SetStats(0);
  hAngleMassSig->Draw("colz");
  c4->cd(2);
  hAngleMassBkg->SetStats(0);
  hAngleMassBkg->Draw("colz");
  c4->Print("controlplot3.pdf");

  // auto nentries = tree->GetEntries();
  // std::list<int> l(nentries);
  // std::iota(l.begin(), l.end(), 0);
  // std::vector<std::list<int>::iterator> v(l.size());
  // std::iota(v.begin(), v.end(), l.begin());
  // std::shuffle(v.begin(), v.end(), std::mt19937{std::random_device{}()});

  auto shuffled_ch = new TTree("Data", "Data");
  float ph1[4], ph2[4];
  shuffled_ch->Branch("photon1", ph1, "ph1[4]/F");
  shuffled_ch->Branch("photon2", ph2, "ph2[4]/F");
  for (int i=0; i<tree->GetEntries(); i++) {
    tree->GetEntry(i);

    photons.at(0).GetXYZT(ph1);
    photons.at(1).GetXYZT(ph2);

    shuffled_ch->Fill();
  }

  delete tree;

  outFile->WriteTObject(shuffled_ch, "Data");
  outFile->Close();
}
