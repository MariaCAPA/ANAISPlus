/*
   save_mean_scintillation_pulse_aligned.cpp
   -----------------------------------------------------------------------------
   Versión modificada para construir el pulso promedio (SER o cualquiera) pero
   alineando previamente cada forma de onda al bin de su máximo.

   Estrategia de alineado:
     1. Se resta el baseline como en la versión original.
     2. Se busca el bin del máximo del pulso (maxBin).
     3. Se define un bin de referencia (alignBin), por ejemplo 500.
     4. Se desplaza la forma de onda para que su máximo coincida con alignBin.

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
    bslTot /= double(nPtsBsl);

    for (int k = 0; k < nPtsBsl; ++k)
        bslTotF += hAux->GetBinContent(nSamples - nPtsBsl + k);
    bslTotF /= double(nPtsBsl);

    for (int k = 0; k < nPtsBsl; ++k)
        rmsTot += TMath::Power(hAux->GetBinContent(k + 1) - bslTot, 2);
    rmsTot = std::sqrt(rmsTot / double(nPtsBsl));

    for (int k = 0; k < nPtsBsl; ++k)
        rmsTotF += TMath::Power(hAux->GetBinContent(nSamples - nPtsBsl + k) - bslTotF, 2);
    rmsTotF = std::sqrt(rmsTotF / double(nPtsBsl));

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
    TFile* file0 = nullptr;
    TFile* file1 = nullptr;
    {
        std::ostringstream pathW;
        std::ostringstream pathL1;
        if (run < 10)
        {
            pathW << "../analyzed/ANAISplusWaveforms/run0000" << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run0000" << run << ".mid.lz4_analyzed_1.root";
        }
        else if (run < 100)
        {
            pathW << "../analyzed/ANAISplusWaveforms/run000" << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run000" << run << ".mid.lz4_analyzed_1.root";
        }
        else if (run < 1000)
        {
            pathW << "../analyzed/ANAISplusWaveforms/run00" << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run00" << run << ".mid.lz4_analyzed_1.root";
        }
        else
        {
            pathW << "../analyzed/ANAISplusWaveforms/run0" << run << ".mid.lz4_waveforms.root";
            pathL1 << "../analyzed/ANAISplusLevel1/run0" << run << ".mid.lz4_analyzed_1.root";
        }
        file0 = TFile::Open(pathW.str().c_str(), "READ");
        file1 = TFile::Open(pathL1.str().c_str(), "READ");
    }

    if (!file0 || file0->IsZombie() || !file1 || file1->IsZombie())
    {
        std::cerr << "Error abriendo archivos de entrada.\n";
        return;
    }

    TTree* ta = static_cast<TTree*>(file0->Get("ta"));
    TTree* ta_control = static_cast<TTree*>(file1->Get("ta"));
    if (!ta || !ta_control)
    {
        std::cerr << "No se encontró el TTree 'ta'.\n";
        return;
    }

    double area[8] = {0.};
    double t00 = 0.;

    for (int i = 0; i < 8; ++i)
        ta_control->SetBranchAddress(Form("areaFix%d", i), &area[i]);
    ta_control->SetBranchAddress("tmax0", &t00);

    const int nChannels = 1;
    TH1F** hP = new TH1F*[nChannels];
    for (int j = 0; j < nChannels; ++j)
    {
        hP[j] = new TH1F(Form("hP_%d", j), Form("hP_%d", j), nSamples, 0, nSamples);
        ta->SetBranchAddress(Form("pulse%d", j), &hP[j]);
    }

    TH1F* mean_waveform = new TH1F("mean_pulse", "mean_pulse", nSamples, 0, nSamples);

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
            continue;

        TH1F* pulse = hP[0];
        bsl_provisional = bsl_simple(pulse, nPtsBsl, 0, &rms_provisional, nSamples);
        for (int k = 1; k <= nSamples; ++k)
            pulse->SetBinContent(k, pulse->GetBinContent(k) - bsl_provisional);

        int maxBin = pulse->GetMaximumBin();
        int shift = alignBin - maxBin;

        for (int k = 1; k <= nSamples; ++k)
        {
            int newBin = k + shift;
            if (newBin < 1 || newBin > nSamples)
                continue;
            mean_waveform->AddBinContent(newBin, pulse->GetBinContent(k));
        }

        ++selected_entries;
    }

    if (selected_entries > 0)
        mean_waveform->Scale(1.0 / selected_entries);
    else
        std::cerr << "Advertencia: no se seleccionó ningún pulso.\n";

    TString outdir = TString(gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/"));
    gSystem->mkdir(outdir, true);

    TString foutname = Form("aligned_SER_max_run%05d.root", run);
    TFile fout(outdir + foutname, "RECREATE");
    mean_waveform->Write();
    fout.Close();

    std::cout << "Pulso promedio alineado (por máximo) guardado en: " << (outdir + foutname) << "\n";

    delete[] hP;
    delete mean_waveform;
}

