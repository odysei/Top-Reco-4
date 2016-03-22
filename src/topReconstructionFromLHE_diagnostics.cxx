#ifndef topReconstructionFromLHE_diagnostics_cxx
#define topReconstructionFromLHE_diagnostics_cxx

#include "topReconstructionFromLHE.h"
#include "Math/GenVector/LorentzVector.h"

using namespace std;

void topReconstructionFromLHE::printVector(XYZTLorentzVector &v)
{
    cout << "px = " << v.Px() << endl;
    cout << "py = " << v.Py() << endl;
    cout << "pz = " << v.Pz() << endl;
    cout << "E  = " << v.E() << endl;
}

void topReconstructionFromLHE::PrintPt(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Pt true = " << evh.trueParticlesLH[whichParticle]->Pt() << endl;
    cout << "Pt smeared = " << evh.smearedParticlesLH[whichParticle]->Pt()
         << endl;
    cout << "Pt best = " << evh.bestParticlesLH[whichParticle]->Pt() << endl;
    cout << "Pt smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Pt() -
                evh.trueParticlesLH[whichParticle]->Pt()
         << endl;
    cout << "Pt best - true = "
         << evh.bestParticlesLH[whichParticle]->Pt() -
                evh.trueParticlesLH[whichParticle]->Pt()
         << endl;
}

void topReconstructionFromLHE::PrintPhi(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Phi true = " << evh.trueParticlesLH[whichParticle]->Phi() << endl;
    cout << "Phi smeared = " << evh.smearedParticlesLH[whichParticle]->Phi()
         << endl;
    cout << "Phi best = " << evh.bestParticlesLH[whichParticle]->Phi() << endl;
    cout << "Phi smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Phi() -
                evh.trueParticlesLH[whichParticle]->Phi()
         << endl;
    cout << "Phi best - true = "
         << evh.bestParticlesLH[whichParticle]->Phi() -
                evh.trueParticlesLH[whichParticle]->Phi()
         << endl;
}

void topReconstructionFromLHE::PrintEta(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Eta true = " << evh.trueParticlesLH[whichParticle]->Eta() << endl;
    cout << "Eta smeared = " << evh.smearedParticlesLH[whichParticle]->Eta()
         << endl;
    cout << "Eta best = " << evh.bestParticlesLH[whichParticle]->Eta() << endl;
    cout << "Eta smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->Eta() -
                evh.trueParticlesLH[whichParticle]->Eta()
         << endl;
    cout << "Eta best - true = "
         << evh.bestParticlesLH[whichParticle]->Eta() -
                evh.trueParticlesLH[whichParticle]->Eta()
         << endl;
}

void topReconstructionFromLHE::PrintMass(string whichParticle, handleEvent &evh)
{
    cout << whichParticle << endl;
    cout << "Mass true = " << evh.trueParticlesLH[whichParticle]->M() << endl;
    cout << "Mass smeared = " << evh.smearedParticlesLH[whichParticle]->M()
         << endl;
    cout << "Mass best = " << evh.bestParticlesLH[whichParticle]->M() << endl;
    cout << "Mass smeared - true = "
         << evh.smearedParticlesLH[whichParticle]->M() -
                evh.trueParticlesLH[whichParticle]->M()
         << endl;
    cout << "Mass best - true = "
         << evh.bestParticlesLH[whichParticle]->M() -
                evh.trueParticlesLH[whichParticle]->M()
         << endl;
}

void topReconstructionFromLHE::Print()
{
    inFilePlot = new TFile("output_files/output_0.root");
    inTreePlot = (TTree *)inFilePlot->Get("tree");

    handleEvent evh(particleNames, names, chinames);
    DeclareInBranchesForPlotting(evh);
    DeclareHists();

    int numEvents = inTreePlot->GetEntries();

    for (int i = 0; i < numEvents; i++) {
        inTreePlot->GetEntry(i);
        if ((innerMinStatus == 0 or innerMinStatus == 1) and
            (outerMinStatus == 0 or outerMinStatus == 1)) {
            cout << endl;
            cout << "eventNumber = " << eventNumber << endl;
            PrintMass("Leptonic_W", evh);
            PrintPt("Leptonic_Bottom", evh);
            PrintPhi("Leptonic_Bottom", evh);
            PrintEta("Leptonic_Bottom", evh);
            PrintPt("Neutrino_or_AntiNeutrino", evh);
            PrintPhi("Neutrino_or_AntiNeutrino", evh);
            PrintEta("Neutrino_or_AntiNeutrino", evh);
            PrintMass("Leptonic_Top", evh);
            PrintMass("Hadronic_W", evh);
            PrintPt("Hadronic_Bottom", evh);
            PrintPt("Quark_from_W", evh);
            PrintPt("Antiquark_from_W", evh);

            for (vector<string>::const_iterator t = chinames.begin();
                 t < chinames.end(); t++) {
                cout << *t << " " << evh.chiSquareds[*t] << endl;
            }

            /*     cout<<"Px total = "<<
               evh.bestParticlesLH["Leptonic_Bottom"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Hadronic_Bottom"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Lepton_or_Antilepton"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px()
               <<endl;
                 cout<<evh.bestParticlesLH["Quark_from_W"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Antiquark_from_W"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["B_from_H"]->Px() <<endl;
                 cout<<evh.bestParticlesLH["Bbar_from_H"]->Px() << endl; */

            cout << "Px total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Px() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Px() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Px() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px() +
                        evh.bestParticlesLH["Quark_from_W"]->Px() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Px() +
                        evh.bestParticlesLH["B_from_H"]->Px() +
                        evh.bestParticlesLH["Bbar_from_H"]->Px()
                 << endl;
            cout << "Py total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Py() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Py() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Py() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Py() +
                        evh.bestParticlesLH["Quark_from_W"]->Py() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Py() +
                        evh.bestParticlesLH["B_from_H"]->Py() +
                        evh.bestParticlesLH["Bbar_from_H"]->Py()
                 << endl;
            cout << "Top Px total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Px() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Px() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Px() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Px() +
                        evh.bestParticlesLH["Quark_from_W"]->Px() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Px()
                 << endl;
            cout << "Nontop Px total = "
                 << evh.bestParticlesLH["B_from_H"]->Px() +
                        evh.bestParticlesLH["Bbar_from_H"]->Px()
                 << endl;

            cout << "Top Py total = "
                 << evh.bestParticlesLH["Leptonic_Bottom"]->Py() +
                        evh.bestParticlesLH["Hadronic_Bottom"]->Py() +
                        evh.bestParticlesLH["Lepton_or_AntiLepton"]->Py() +
                        evh.bestParticlesLH["Neutrino_or_AntiNeutrino"]->Py() +
                        evh.bestParticlesLH["Quark_from_W"]->Py() +
                        evh.bestParticlesLH["Antiquark_from_W"]->Py()
                 << endl;
            cout << "Nontop Py total = "
                 << evh.bestParticlesLH["B_from_H"]->Py() +
                        evh.bestParticlesLH["Bbar_from_H"]->Py()
                 << endl;
        }
        // cout<<"after fillhists"<<endl;
    }
}

void topReconstructionFromLHE::Loop_diagnostics(handleEvent &evh)
{
    XYZTLorentzVector gen_all =
        *evh.smearedParticles["bottom"] + *evh.smearedParticles["antiBottom"] +
        *evh.smearedParticles["qFromW"] + *evh.smearedParticles["qbarFromW"] +
        *evh.smearedParticles["bFromH"] + *evh.smearedParticles["bbarFromH"] +
        *evh.smearedParticles["lepton"] + *evh.smearedParticles["antiNeutrino"];
    XYZTLorentzVector best_all =
        *evh.bestParticles["top"] + *evh.bestParticles["antiTop"] +
        *evh.bestParticles["bFromH"] + *evh.bestParticles["bbarFromH"];
    XYZTLorentzVector best_tt =
        *evh.bestParticles["top"] + *evh.bestParticles["antiTop"];

    cout << "Event, pT(gen_all), pT(best_all), pT(tt) = " << eventNumber << ", "
         << gen_all.Pt() << ", " << best_all.Pt() << ", " << best_tt.Pt()
         << endl;

    if (debug_verbosity < 3)
        return;

    const double chi2 =
        evh.chiSquareds["bottom"] + evh.chiSquareds["antiBottom"] +
        evh.chiSquareds["qFromW"] + evh.chiSquareds["qbarFromW"] +
        evh.chiSquareds["lepton_or_antilepton"] +
        evh.chiSquareds["leptonicWMass"];

    cout << "Event, rel_error(ttbb), chi2(bbqqlW) = " << eventNumber << ", "
         << rel_error << ", " << chi2 << endl;
}

void topReconstructionFromLHE::Print_smear_bs_SM(handleEvent &evh)
{
    cout << "Smeared bottom (E, px, py, pz) = ("
         << evh.smearedParticles["bottom"]->E() << ", "
         << evh.smearedParticles["bottom"]->Px() << ", "
         << evh.smearedParticles["bottom"]->Py() << ", "
         << evh.smearedParticles["bottom"]->Pz() << ")\n";

    cout << "Smeared bfromH (E, px, py, pz) = ("
         << evh.smearedParticles["bFromH"]->E() << ", "
         << evh.smearedParticles["bFromH"]->Px() << ", "
         << evh.smearedParticles["bFromH"]->Py() << ", "
         << evh.smearedParticles["bFromH"]->Pz() << ")\n";

    cout << "Smeared bbarFromH (E, px, py, pz) = ("
         << evh.smearedParticles["bbarFromH"]->E() << ", "
         << evh.smearedParticles["bbarFromH"]->Px() << ", "
         << evh.smearedParticles["bbarFromH"]->Py() << ", "
         << evh.smearedParticles["bbarFromH"]->Pz() << ")\n";
}

void topReconstructionFromLHE::Print_smear_tt_SM(handleEvent &evh)
{
    cout << "Smeared top (E, px, py, pz) = ("
         << evh.smearedParticles["top"]->E() << ", "
         << evh.smearedParticles["top"]->Px() << ", "
         << evh.smearedParticles["top"]->Py() << ", "
         << evh.smearedParticles["top"]->Pz() << ")\n";

    cout << "Smeared antiTop (E, px, py, pz) = ("
         << evh.smearedParticles["antiTop"]->E() << ", "
         << evh.smearedParticles["antiTop"]->Px() << ", "
         << evh.smearedParticles["antiTop"]->Py() << ", "
         << evh.smearedParticles["antiTop"]->Pz() << ")\n";
}

#endif
