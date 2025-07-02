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
   TFile *file0;
   if(run<10) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run0000%d.mid.lz4_waveforms.root",run),"READ");
   else if(run<100) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run000%d.mid.lz4_waveforms.root",run),"READ");
   else if(run<1000) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run00%d.mid.lz4_waveforms.root",run),"READ");
   else file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run0%d.mid.lz4_waveforms.root",run),"READ");

   TTree* ta = (TTree*)file0->Get("ta");

   TFile *file1; 
   if(run<10) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
   else if(run<100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
   else if(run<1000) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
   else file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0%d.mid.lz4_analyzed_1.root",run),"READ");

   TTree* ta_control = (TTree*)file1->Get("ta");

   double area[8]; //8 channels digitizer
   double area_tot = 0;
   double t00;
   
   //if some channels are not being used, nothing happens, just don't take them into account
   ta_control->SetBranchAddress("areaFix0",&area[0]);
   ta_control->SetBranchAddress("areaFix1",&area[1]);
   ta_control->SetBranchAddress("areaFix2",&area[2]);
   ta_control->SetBranchAddress("areaFix3",&area[3]);
   ta_control->SetBranchAddress("areaFix4",&area[4]);
   ta_control->SetBranchAddress("areaFix5",&area[5]);
   ta_control->SetBranchAddress("areaFix6",&area[6]);
   ta_control->SetBranchAddress("areaFix7",&area[7]);
   ta_control->SetBranchAddress("tmax0",&t00);

   int nChannels=4;
   TH1F ** hP = new TH1F * [nChannels];

   for(int j = 0; j < nChannels; j++) 
   {
      hP[j] = new TH1F(Form("hP_%d",j),Form("hP_%d",j),nSamples,0,nSamples);
      ta->SetBranchAddress(Form("pulse%d",j),&hP[j]);
   }

   TH1F* hP_aux = new TH1F("hP_aux","hP_aux",nSamples,0,nSamples);

   // X-axis in microseconds
   TH1F * mean_waveform = new TH1F("mean_pulse","mean_pulse",nSamples,0,nSamples); // 0.002 µs per bin

   double bsl_provisional = 0, rms_provisional;
   int nPtsBsl = 500;
   int sIni = 0;
   double selected_entries = 0;

   for(int i = 0; i < ta->GetEntries(); i++)
   { 
      ta->GetEntry(i);
      ta_control->GetEntry(i);

      if(i % 10000 == 0) std::cout << "Entrada " << i << std::endl;

      area_tot = 0;
      for(int j = 0; j < nChannels; j++)
         area_tot += area[j];

      if(area_tot > th && area_tot < th2)
      {
         for(int j = 0; j < 1; j++)
         {
            bsl_provisional = bsl_simple(hP[j],nPtsBsl,sIni,&rms_provisional,nSamples);
            for(int k = 1; k <= nSamples; k++) 
               hP[j]->SetBinContent(k, hP[j]->GetBinContent(k) - bsl_provisional);

            mean_waveform->Add(hP[j],1.0);
            hP_aux->Reset();
         }
         selected_entries++;
      }
   }

   mean_waveform->Scale(1.0 / selected_entries);

   // CREAR DIRECTORIO DE SALIDA
   TString outdir = TString(gSystem->ExpandPathName("~/Escritorio/mean_scintillation_pulse/"));
   gSystem->mkdir(outdir, true);

   // NOMBRE DEL ARCHIVO
   TString foutname = Form("fixed_area_range_run%05d_81keV_m2.root", run);
   TFile *fout = new TFile(outdir + foutname, "RECREATE");
   mean_waveform->Write();
   fout->Close();

   std::cout << "Pulso promedio guardado en: " << (outdir + foutname) << std::endl;

   delete[] hP;
}

