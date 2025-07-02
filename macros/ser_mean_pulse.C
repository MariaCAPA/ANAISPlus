/*
   save_mean_scintillation_pulse_aligned.cpp
   -----------------------------------------------------------------------------
   Versión modificada para construir el pulso promedio (SER o cualquiera) pero
   alineando previamente cada forma de onda al cruce del 20 % de su amplitud
   máxima.

   Estrategia de alineado (método robusto contra jitter):
     1. Se resta el baseline como en la versión original.
     2. Se busca el valor máximo del pulso (peakMax).
     3. Se define un umbral: thCross = peakMax * fracThreshold  (fracThreshold = 0.20).
     4. Se localiza el primer bin cuyo contenido ≥ thCross (cruce por subida).
     5. Se desplaza la forma de onda para que ese bin coincida con el bin
        "alignBin" (p.ej. 1000) del histograma promedio.

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
// Función auxiliar: baseline simple (sin cambios respecto a tu versión)
//--------------------------------------------------------------------
double bsl_simple(TH1F* hAux, int nPtsBsl, int /*sIni*/, double* rms, int nSamples)
{
    double bslTot = 0.;
    double bslTotF = 0.;
    double rmsTot = 0.;
    double rmsTotF = 0.;

    // Baseline al principio
    for (int k = 0; k < nPtsBsl; ++k)
        bslTot += hAux->GetBinContent(k + 1);
    bslTot /= double(nPtsBsl);

    // Baseline al final
    for (int k = 0; k < nPtsBsl; ++k)
        bslTotF += hAux->GetBinContent(nSamples - nPtsBsl + k);
    bslTotF /= double(nPtsBsl);

    // RMS inicio
    for (int k = 0; k < nPtsBsl; ++k)
        rmsTot += TMath::Power(hAux->GetBinContent(k + 1) - bslTot, 2);
    rmsTot = std::sqrt(rmsTot / double(nPtsBsl));

    // RMS final
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
// Función principal: genera el pulso promedio alineado
//--------------------------------------------------------------------
void save_mean_scintillation_pulse_aligned(int run = 0,
                                           double th = 20., double th2 = 30.,
                                           int nSamples = 200000,
                                           int alignBin = 1000,
                                           double fracThreshold = 0.20)
{
    // -----------------------------------------------------------------
    // 1. Abrir ficheros de entrada
    // -----------------------------------------------------------------
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

    // -----------------------------------------------------------------
    // 2. Configurar ramas de control (áreas, tmax)
    // -----------------------------------------------------------------
    double area[8] = {0.};
    double t00 = 0.;

    for (int i = 0; i < 8; ++i)
        ta_control->SetBranchAddress(Form("areaFix%d", i), &area[i]);
    ta_control->SetBranchAddress("tmax0", &t00);

    // -----------------------------------------------------------------
    // 3. Configurar lectura de pulsos (usaremos sólo canal 0 para SER)
    // -----------------------------------------------------------------
    const int nChannels = 1;  // solo canal 0 para SER
    TH1F** hP = new TH1F*[nChannels];
    for (int j = 0; j < nChannels; ++j)
    {
        hP[j] = new TH1F(Form("hP_%d", j), Form("hP_%d", j), nSamples, 0, nSamples);
        ta->SetBranchAddress(Form("pulse%d", j), &hP[j]);
    }

    TH1F* mean_waveform = new TH1F("mean_pulse", "mean_pulse", nSamples, 0, nSamples);

    // Auxiliares
    const int nPtsBsl = 500;
    double bsl_provisional = 0.;
    double rms_provisional = 0.;
    double selected_entries = 0.;

    // -----------------------------------------------------------------
    // 4. Recorrer eventos y construir promedio alineado
    // -----------------------------------------------------------------
    for (Long64_t i = 0; i < ta->GetEntries(); ++i)
    {
        ta->GetEntry(i);
        ta_control->GetEntry(i);

        if (i % 10000 == 0)
            std::cout << "Evento " << i << "\n";

        // Área total (suma de canales relevantes)
        double area_tot = 0.;
        for (int j = 0; j < nChannels; ++j)
            area_tot += area[j];

        // Selección 1 phe (rango [th, th2])
        if (area_tot < th || area_tot > th2)
            continue;

        // ----------------------------------------------------------------
        // 4a. Baseline y pulso ya con baseline corregido (solo canal 0)
        // ----------------------------------------------------------------
        TH1F* pulse = hP[0];
        bsl_provisional = bsl_simple(pulse, nPtsBsl, 0, &rms_provisional, nSamples);
        for (int k = 1; k <= nSamples; ++k)
            pulse->SetBinContent(k, pulse->GetBinContent(k) - bsl_provisional);

        // ----------------------------------------------------------------
        // 4b. Alineación por cruce de umbral (20 % del pico)
        // ----------------------------------------------------------------
        int maxBin = pulse->GetMaximumBin();
        double maxVal = pulse->GetBinContent(maxBin);
        if (maxVal <= 0.)
            continue;  // seguridad
        double thresholdVal = maxVal * fracThreshold;

        // Buscar primer bin que cruza threshold
        int crossingBin = -1;
        for (int k = 1; k <= maxBin; ++k)  // buscamos desde el inicio hasta el pico
        {
            if (pulse->GetBinContent(k) >= thresholdVal)
            {
                crossingBin = k;
                break;
            }
        }
        if (crossingBin < 0)
            continue;  // no encontró cruce

        int shift = alignBin - crossingBin;  // desplazamiento

        // ----------------------------------------------------------------
        // 4c. Acumular en histograma promedio con shift
        // ----------------------------------------------------------------
        for (int k = 1; k <= nSamples; ++k)
        {
            int newBin = k + shift;
            if (newBin < 1 || newBin > nSamples)
                continue;  // evitamos overflow
            mean_waveform->AddBinContent(newBin, pulse->GetBinContent(k));
        }

        ++selected_entries;
    }

    // -----------------------------------------------------------------
    // 5. Normalizar promedio final
    // -----------------------------------------------------------------
    if (selected_entries > 0)
        mean_waveform->Scale(1.0 / selected_entries);
    else
        std::cerr << "Advertencia: no se seleccionó ningún pulso.\n";

    // -----------------------------------------------------------------
    // 6. Guardar resultado en archivo ROOT
    // -----------------------------------------------------------------
    TString outdir = TString(gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/"));
    gSystem->mkdir(outdir, true);

    TString foutname = Form("aligned_SER_run%05d.root", run);
    TFile fout(outdir + foutname, "RECREATE");
    mean_waveform->Write();
    fout.Close();

    std::cout << "Pulso promedio alineado guardado en: " << (outdir + foutname) << "\n";

    // Limpieza
    delete[] hP;
    delete mean_waveform;
}

