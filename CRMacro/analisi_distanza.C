
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TTree.h"

void analisi_distanza()
{
    // =========================================
    // PARAMETRI MODIFICABILI
    // =========================================

    const char* nomeFile = "CosmicMC_6boxes.root";
    const char* nomeNtuple = "tracks_box1";

    const int numeroBin = 100;
    const double distanzaMin = 0.0;
    const double distanzaMax = 0.5; // metri

    // =========================================
    // APERTURA FILE
    // =========================================

    TFile *f = TFile::Open(nomeFile, "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "ERRORE: impossibile aprire "
                  << nomeFile << std::endl;
        return;
    }

    std::cout << "\nContenuto del file ROOT:\n";
    f->ls();

    // =========================================
    // LETTURA TNtuple
    // =========================================

    TNtuple *t = dynamic_cast<TNtuple*>(f->Get(nomeNtuple));

    if (!t) {
        std::cerr << "\nERRORE: TNtuple "
                  << nomeNtuple << " non trovata!"
                  << std::endl;
        f->Close();
        return;
    }

    std::cout << "\nTNtuple trovata: " << nomeNtuple << std::endl;
    std::cout << "Numero di eventi: " << t->GetEntries() << std::endl;

    std::cout << "\nStruttura della TNtuple:\n";
    t->Print();

    // =========================================
    // VERIFICA RAMI
    // =========================================

    const char* rami[] = {
        "xIn", "yIn", "zIn",
        "xOut", "yOut", "zOut"
    };

    for (const char* ramo : rami) {
        if (!t->GetBranch(ramo)) {
            std::cerr << "ERRORE: ramo " << ramo
                      << " non presente!" << std::endl;
            f->Close();
            return;
        }
    }

    // =========================================
    // CREAZIONE ISTOGRAMMA
    // =========================================

    TH1F *h = new TH1F(
        "hDistanza",
        "Distanza percorsa nel parallelepipedo;"
        "Distanza [m];Numero di raggi",
        numeroBin,
        distanzaMin,
        distanzaMax
    );

    h->SetDirectory(nullptr);

    // =========================================
    // CALCOLO DISTANZA
    // =========================================

    Long64_t n = t->GetEntries();
    Long64_t validi = 0;

    for (Long64_t i = 0; i < n; i++) {

        t->GetEntry(i);

        double xIn  = t->GetLeaf("xIn")->GetValue();
        double yIn  = t->GetLeaf("yIn")->GetValue();
        double zIn  = t->GetLeaf("zIn")->GetValue();

        double xOut = t->GetLeaf("xOut")->GetValue();
        double yOut = t->GetLeaf("yOut")->GetValue();
        double zOut = t->GetLeaf("zOut")->GetValue();

        double dx = xOut - xIn;
        double dy = yOut - yIn;
        double dz = zOut - zIn;

        double distanza = std::sqrt(
            dx*dx + dy*dy + dz*dz
        );

        h->Fill(distanza);
        validi++;
    }

    std::cout << "\nRaggi analizzati: " << validi << std::endl;
    std::cout << "Distanza media: " << h->GetMean() << " m"
              << std::endl;
    std::cout << "Distanza massima nell'istogramma: "
              << h->GetXaxis()->GetXmax() << " m"
              << std::endl;

    // =========================================
    // DISEGNO
    // =========================================

    TCanvas *c = new TCanvas(
        "cDistanza",
        "Distanza percorsa",
        900,
        600
    );

    h->SetLineWidth(2);
    h->Draw("HIST");

    c->Update();

    f->Close();

    std::cout << "\nAnalisi completata." << std::endl;
}
