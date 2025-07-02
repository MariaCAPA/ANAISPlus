double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples)
{
  double bslTot = 0;
  double bslTotF = 0;
  double rmsTot = 0;
  double rmsTotF = 0;

  for (int k = 0; k < nPtsBsl; k++) bslTot += hAux->GetBinContent(k + 1) / double(nPtsBsl);
  for (int k = 0; k < nPtsBsl; k++) bslTotF += hAux->GetBinContent(nSamples - nPtsBsl + k) / double(nPtsBsl);

  for (int k = 0; k < nPtsBsl; k++) rmsTot += pow(hAux->GetBinContent(k + 1) - bslTot, 2) / double(nPtsBsl); 
  for (int k = 0; k < nPtsBsl; k++) rmsTotF += pow(hAux->GetBinContent(nSamples - nPtsBsl + k) - bslTotF, 2) / double(nPtsBsl);

  rmsTot = sqrt(rmsTot);
  rmsTotF = sqrt(rmsTotF);

  if (rmsTot < rmsTotF) { *rms = rmsTot; return bslTot; }
  else { *rms = rmsTotF; return bslTotF; }
}


void save_mean_scintillation_pulse(int run=0, double th=20, double th2=30, int nSamples=200000)
{

  //---------------------------------------------------
  // 1. Abrir ficheros y árboles
  //---------------------------------------------------
  TString wfname  = Form("../analyzed/ANAISplusWaveforms/run%05d.mid.lz4_waveforms.root",  run);
  TString cfname  = Form("../analyzed/ANAISplusLevel1/run%05d.mid.lz4_analyzed_1.root",    run);

  TFile  fw(wfname, "READ");             if (fw.IsZombie()) { printf("No waveforms\n"); return; }
  TFile  fc(cfname, "READ");             if (fc.IsZombie()) { printf("No control\n");   return; }

  TTree *ta         = (TTree*)fw.Get("ta");
  TTree *ta_control = (TTree*)fc.Get("ta");

  //---------------------------------------------------
  // 2. Crear una entry-list solo con eventos deseados
  //---------------------------------------------------
  // NOTA: areaFixTot está siempre; si un día cambiase el nombre, basta tocar aquí.
  TString cut = Form("areaFixTot > %f && areaFixTot < %f", th, th2);
  ta_control->Draw(">>elist", cut, "entrylist");
  auto elist = (TEntryList*)gDirectory->Get("elist");
  Long64_t nSel = elist->GetN();
  printf("Eventos seleccionados: %lld\n", nSel);
  if (nSel == 0) return;

  //---------------------------------------------------
  // 3. Detectar cuántos canales de pulsos hay realmente
  //---------------------------------------------------
  std::vector<int> chList;
  for (int j = 0; j < 8; ++j)
    if (ta->GetBranch(Form("pulse%d", j))) chList.push_back(j);

  int nChannels = chList.size();
  printf("Canales detectados: %d\n", nChannels);

  //---------------------------------------------------
  // 4. Reservar histogramas y enganchar ramas
  //---------------------------------------------------
  std::vector<TH1F*> hP(nChannels, nullptr);
  for (int idx = 0; idx < nChannels; ++idx) {
      int j = chList[idx];
      hP[idx] = new TH1F(Form("hP_%d", j), "", nSamples, 0, nSamples);
      ta->SetBranchAddress(Form("pulse%d", j), &hP[idx]);
  }

  TH1F *mean_waveform = new TH1F("mean_pulse", "mean_pulse", nSamples, 0, nSamples);

  //---------------------------------------------------
  // 5. Loop sobre SOLO los eventos válidos
  //---------------------------------------------------
  constexpr int nPtsBsl = 500;
  double bsl, rms;
  int    sIni = 0;
  for (Long64_t i = 0; i < nSel; ++i) {
      Long64_t ev = elist->GetEntry(i);
      ta->GetEntry(ev);

      for (auto *hp : hP) {
          bsl = bsl_simple(hp, nPtsBsl, sIni, &rms, nSamples);
          for (int k = 1; k <= nSamples; ++k)
              hp->SetBinContent(k, hp->GetBinContent(k) - bsl);

          mean_waveform->Add(hp);
      }
  }
  mean_waveform->Scale(1.0 / nSel / nChannels);

  //---------------------------------------------------
  // 6. Guardar resultado
  //---------------------------------------------------
  TString outdir = gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/");
  gSystem->mkdir(outdir, /*recursive=*/true);

  TString foutname = Form("fast_fixed_area_range_run%05d_81keV_m2.root", run);
  TFile   fout(outdir + foutname, "RECREATE");
  mean_waveform->Write();
  fout.Close();

  printf("Pulso promedio guardado en: %s\n", (outdir + foutname).Data());
}
