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

void fast_mean_pulse(int run       = 0,
                     double thLow  = 1.6e6,
                     double thHigh = 2.6e6,
                     int nSamples  = 200000,
                     bool saveSub  = true)
{
  // --------------------------------------------
  // 1. Abrir archivos y árboles
  // --------------------------------------------
  TString wfName = Form("../analyzed/ANAISplusWaveforms/run%05d.mid.lz4_waveforms.root",  run);
  TString ctName = Form("../analyzed/ANAISplusLevel1/run%05d.mid.lz4_analyzed_1.root",    run);

  TFile fw(wfName, "READ");   if (fw.IsZombie()) { printf("❌ waveforms\n"); return; }
  TFile fc(ctName, "READ");   if (fc.IsZombie()) { printf("❌ control\n");   return; }

  TTree *ta  = (TTree*)fw.Get("ta");
  TTree *ctl = (TTree*)fc.Get("ta");
  ta->AddFriend(ctl, "ctl");

  // --------------------------------------------
  // 2. Activar solo ramas necesarias
  // --------------------------------------------
  ta->SetBranchStatus("*", 0);  // desactiva todo
  std::vector<int> chList;
  for (int j = 0; j < 8; ++j) {
    if (ta->GetBranch(Form("pulse%d", j))) {
      ta->SetBranchStatus(Form("pulse%d", j), 1);
      chList.push_back(j);
    }
  }
  ta->SetBranchStatus("ctl.areaFixTot", 1);  // ✅ habilitar para el corte

  const int nCh = chList.size();
  printf("Canales activos: %d  ->  ", nCh);
  for (int c : chList) printf("pulse%d ", c);
  printf("\n");

  // --------------------------------------------
  // 3. Crear archivo de salida antes de CopyTree()
  // --------------------------------------------
  TString outdir = gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/");
  gSystem->mkdir(outdir, true);
  TString outName = Form("%sfast_mean_run%05d.root", outdir.Data(), run);
  TFile *fout = new TFile(outName, "RECREATE");  // ✅ archivo antes del árbol

  // --------------------------------------------
  // 4. Copiar subárbol con corte
  // --------------------------------------------
  TString cut = Form("ctl.areaFixTot > %f && ctl.areaFixTot < %f", thLow, thHigh);
  printf("Copiando con corte:  %s\n", cut.Data());
  TTree *sub = ta->CopyTree(cut);  // ✅ ahora sí estará vinculado a `fout`
  Long64_t nSel = sub->GetEntries();
  printf("Eventos seleccionados: %lld\n", nSel);
  if (nSel == 0) { printf("Nada que hacer.\n"); fout->Close(); return; }

  // --------------------------------------------
  // 5. Enganchar las ramas pulseX del sub-árbol
  // --------------------------------------------
  std::vector<TH1F*> pulse(nCh, nullptr);
  for (int k = 0; k < nCh; ++k)
    sub->SetBranchAddress(Form("pulse%d", chList[k]), &pulse[k]);

  // --------------------------------------------
  // 6. Calcular el pulso promedio
  // --------------------------------------------
  TH1F *mean = new TH1F("mean_pulse", "mean_pulse", nSamples, 0, nSamples);
  constexpr int nPtsBsl = 500;
  double bsl, rms;
  int sIni = 0;

  for (Long64_t i = 0; i < nSel; ++i) {
    sub->GetEntry(i);
    for (auto *h : pulse) {
      bsl = bsl_simple(h, nPtsBsl, sIni, &rms, nSamples);
      for (int b = 1; b <= nSamples; ++b)
        h->SetBinContent(b, h->GetBinContent(b) - bsl);
      mean->Add(h);
    }
  }
  mean->Scale(1.0 / nSel / nCh);

  // --------------------------------------------
  // 7. Guardar resultados
  // --------------------------------------------
  mean->Write();
  if (saveSub) {
    sub->Write("ta_reduced");  // ya está vinculado a fout
    printf("Árbol reducido escrito como 'ta_reduced' (%lld evts)\n", nSel);
  }

  fout->Close();
  printf("✅ Pulso promedio guardado en: %s\n", outName.Data());
}

