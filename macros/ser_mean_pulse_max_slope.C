/*
   save_mean_scintillation_pulse_aligned_posSlope.cpp
   -----------------------------------------------------------------------------
   Construye el pulso promedio (SER o el que toque) alineando **cada** waveform
   al bin donde la pendiente (flanco de subida) es **máxima y positiva**.

   Estrategia de alineado:
     1. Restar el baseline (idéntico a la versión original).
     2. Calcular la pendiente discreta:  slope(k) = P(k+1) - P(k).
     3. Buscar el bin con slope máxima (positiva).
     4. Desplazar la forma de onda para que ese bin coincida con `alignBin`.

   -----------------------------------------------------------------------------
*/

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TMath.h>
#include <iostream>
#include <sstream>
#include <cmath>

//--------------------------------------------------------------------
double bsl_simple(TH1F* hAux, int nPtsBsl, int /*sIni*/, double* rms, int nSamples)
{
    double bslTot = 0.;
    double bslTotF = 0.;
    double rmsTot = 0.;
    double rmsTotF = 0.;

    for (int k = 0; k < nPtsBsl; ++k)
        bslTot += hAux->GetBinContent(k + 1);
    bslTot /= static_cast<double>(nPtsBsl);

    for (int k = 0; k < nPtsBsl; ++k)
        bslTotF += hAux->GetBinContent(nSamples - nPtsBsl + k);
    bslTotF /= static_cast<double>(nPtsBsl);

    for (int k = 0; k < nPtsBsl; ++k)
        rmsTot += TMath::Power(hAux->GetBinContent(k + 1) - bslTot, 2);
    rmsTot = std::sqrt(rmsTot / static_cast<double>(nPtsBsl));

    for (int k = 0; k < nPtsBsl; ++k)
        rmsTotF += TMath::Power(hAux->GetBinContent(nSamples - nPtsBsl + k) - bslTotF, 2);
    rmsTotF = std::sqrt(rmsTotF / static_cast<double>(nPtsBsl));

    if (rmsTot < rmsTotF)
    {
        *rms = rmsTot;
        return bslTot;
    }
    else
    {
        *rms = rmsTotF;
        return bslTotF;
    }
}

//--------------------------------------------------------------------
void save_mean_scintillation_pulse_aligned(int run = 0,
                                           double th = 20., double th2 = 30.,
                                           int nSamples = 1000,
                                           int alignBin = 500)
{
    // ----------------------------------------------------------------
    // 1. Abrir archivos ROOT de entrada
    // ----------------------------------------------------------------
    TFile* fileW = nullptr;
    TFile* fileL1 = nullptr;
    {
        std::ostringstream pathW, pathL1;
        if (run < 10)
        {
            pathW  << "../analyzed/ANAISplusWaveforms/run0000" << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run0000"     << run << ".mid.lz4_analyzed_1.root";
        }
        else if (run < 100)
        {
            pathW  << "../analyzed/ANAISplusWaveforms/run000"  << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run000"      << run << ".mid.lz4_analyzed_1.root";
        }
        else if (run < 1000)
        {
            pathW  << "../analyzed/ANAISplusWaveforms/run00"   << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run00"       << run << ".mid.lz4_analyzed_1.root";
        }
        else
        {
            pathW  << "../analyzed/ANAISplusWaveforms/run0"    << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run0"        << run << ".mid.lz4_analyzed_1.root";
        }

        fileW  = TFile::Open(pathW.str().c_str(),  "READ");
        fileL1 = TFile::Open(pathL1.str().c_str(), "READ");
    }

    if (!fileW || fileW->IsZombie() || !fileL1 || fileL1->IsZombie())
    {
        std::cerr << "Error abriendo archivos de entrada.\n";
        return;
    }

    TTree* ta         = static_cast<TTree*>(fileW ->Get("ta"));
    TTree* ta_control = static_cast<TTree*>(fileL1->Get("ta"));
    if (!ta || !ta_control)
    {
        std::cerr << "No se encontró el TTree 'ta'.\n";
        return;
    }

    // ----------------------------------------------------------------
    // 2. Set-up de ramas de interes
    // ----------------------------------------------------------------
    double area[8] = {0.};
    double t00 = 0.;

    for (int i = 0; i < 8; ++i)
        ta_control->SetBranchAddress(Form("areaFix%d", i), &area[i]);
    ta_control->SetBranchAddress("tmax0", &t00);

    const int nChannels = 1;  // cambiar si necesitas >1 canal
    TH1F** hP = new TH1F*[nChannels];
    for (int j = 0; j < nChannels; ++j)
    {
        hP[j] = new TH1F(Form("hP_%d", j), Form("hP_%d", j), nSamples, 0, nSamples);
        ta->SetBranchAddress(Form("pulse%d", j), &hP[j]);
    }

    TH1F* mean_waveform = new TH1F("mean_pulse", "mean_pulse", nSamples, 0, nSamples);

    // ----------------------------------------------------------------
    // 3. Recorrido de eventos
    // ----------------------------------------------------------------
    const int nPtsBsl = 50;
    double bsl_provisional = 0.;
    double rms_provisional = 0.;
    double selected_entries = 0.;

    for (Long64_t i = 0; i < ta->GetEntries(); ++i)
    {
        ta->GetEntry(i);
        ta_control->GetEntry(i);

        if (i % 10000 == 0)
            std::cout << "Evento " << i << "\n";

        double area_tot = 0.;
        for (int j = 0; j < nChannels; ++j)
            area_tot += area[j];

        if (area_tot < th || area_tot > th2)
            continue;  // fuera de la ventana de energía

        TH1F* pulse = hP[0];

        // --- 3a) Baseline ---------------------------------------------------
        bsl_provisional = bsl_simple(pulse, nPtsBsl, 0, &rms_provisional, nSamples);
        for (int k = 1; k <= nSamples; ++k)
            pulse->SetBinContent(k, pulse->GetBinContent(k) - bsl_provisional);

        // --- 3b) BIN DE MÁXIMA PENDIENTE POSITIVA ---------------------------
        int slopeBin = 1;
        double maxSlope = -1e30;
        for (int k = 1; k <= nSamples - 1; ++k) {   // hasta nSamples-1 por el término k+1
            double slope = pulse->GetBinContent(k + 1) - pulse->GetBinContent(k);
            if (slope > maxSlope) {
                maxSlope = slope;
                slopeBin = k;
            }
        }

        // --- 3c) Desplazamiento y suma en histograma promedio ---------------
        int shift = alignBin - slopeBin;

        for (int k = 1; k <= nSamples; ++k)
        {
            int newBin = k + shift;
            if (newBin < 1 || newBin > nSamples)
                continue;  // evita salir de rango
            mean_waveform->AddBinContent(newBin, pulse->GetBinContent(k));
        }

        ++selected_entries;
    }

    // ----------------------------------------------------------------
    // 4. Normalización y escritura a disco
    // ----------------------------------------------------------------
    if (selected_entries > 0)
        mean_waveform->Scale(1.0 / selected_entries);
    else
        std::cerr << "Advertencia: no se seleccionó ningún pulso.\n";

    TString outdir = TString(gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/"));
    gSystem->mkdir(outdir, true);

    TString foutname = Form("aligned_SER_posSlope_run%05d.root", run);
    TFile fout(outdir + foutname, "RECREATE");
    mean_waveform->Write();
    fout.Close();

    std::cout << "Pulso promedio alineado (por máxima pendiente positiva) guardado en: "
              << (outdir + foutname) << "\n";

    // ----------------------------------------------------------------
    // 5. Limpieza
    // ----------------------------------------------------------------
    delete[] hP;
    delete mean_waveform;
}

