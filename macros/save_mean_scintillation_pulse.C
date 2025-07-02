#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <cmath>

// DEFINICIÓN DE bsl_simple
double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples)
{
  double bslTot = 0;
  double bslTotF = 0;
  double rmsTot = 0;
  double rmsTotF = 0;
  int k;

  for(k = 0; k < nPtsBsl; k++) bslTot += hAux->GetBinContent(k+1)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) bslTotF += hAux->GetBinContent(nSamples - nPtsBsl + k)/double(nPtsBsl);
  
  for(k = 0; k < nPtsBsl; k++) rmsTot += pow(hAux->GetBinContent(k+1) - bslTot, 2)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) rmsTotF += pow(hAux->GetBinContent(nSamples - nPtsBsl + k) - bslTotF, 2)/double(nPtsBsl);
  
  rmsTot = sqrt(rmsTot);
  rmsTotF = sqrt(rmsTotF);

  if(rmsTot < rmsTotF) { *rms = rmsTot; return bslTot; }
  else { *rms = rmsTotF; return bslTotF; }
}

// FUNCIÓN PRINCIPAL
void save_mean_scintillation_pulse(int run = 0, double th = 20, double th2 = 30, int nSamples = 200000)
{
   TString filename_wf = Form("../analyzed/ANAISplusWaveforms/run%05d.mid.lz4_waveforms.root", run);
   TString filename_ctrl = Form("../analyzed/ANAISplusLevel1/run%05d.mid.lz4_analyzed_1.root", run);

   TFile *file0 = TFile::Open(filename_wf, "READ");
   TFile *file1 = TFile::Open(filename_ctrl, "READ");
   if (!file0 || file0->IsZombie() || !file1 || file1->IsZombie()) {
      std::cerr << "Error al abrir los archivos." << std::endl;
      return;
   }

   TTree* ta = (TTree*)file0->Get("ta");
   TTree* ta_control = (TTree*)file1->Get("ta");

   double area[8], t00;
   for (int i = 0; i < 8; i++)
      ta_control->SetBranchAddress(Form("areaFix%d", i), &area[i]);
   ta_control->SetBranchAddress("tmax0", &t00);

   int nChannels = 1;
   TH1F** hP = new TH1F*[nChannels];
   for (int j = 0; j < nChannels; j++)
      ta->SetBranchAddress(Form("pulse%d", j), &hP[j]);

   TH1F* mean_waveform = new TH1F("mean_pulse", "Mean Scintillation Pulse", nSamples, 0., nSamples * 0.002); // x en µs

   double bsl, rms;
   int nPtsBsl = 500;
   int sIni = 0;
   double area_tot = 0;
   int selected_entries = 0;

   Long64_t nEntries = ta->GetEntries();
   for (Long64_t i = 0; i < nEntries; i++) {
      ta->GetEntry(i);
      ta_control->GetEntry(i);

      if (i % 10000 == 0) std::cout << "Procesando entrada " << i << " / " << nEntries << std::endl;

      area_tot = 0;
      for (int j = 0; j < nChannels; j++)
         area_tot += area[j];

      if (area_tot > th && area_tot < th2) {
         for (int j = 0; j < 1; j++) {
            bsl = bsl_simple(hP[j], nPtsBsl, sIni, &rms, nSamples);
            for (int k = 1; k <= nSamples; k++) {
               double val = hP[j]->GetBinContent(k) - bsl;
               mean_waveform->AddBinContent(k, val);
            }
         }
         selected_entries++;
      }
   }

   if (selected_entries > 0)
      mean_waveform->Scale(1.0 / selected_entries);
   else {
      std::cerr << "No se seleccionaron eventos. Nada que guardar." << std::endl;
      return;
   }

   TString outdir = TString(gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/"));
   gSystem->mkdir(outdir, true); // crea el directorio si no existe

   TString outfile = Form("run%05d_%.0f_%.0f.root", run, th, th2);
   TFile *fout = new TFile(outdir + outfile, "RECREATE");
   mean_waveform->Write();
   fout->Close();

   std::cout << "Pulso promedio guardado en: " << outdir + outfile << std::endl;

   delete[] hP;
}

