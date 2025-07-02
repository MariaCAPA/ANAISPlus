double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples);

void mean_scintillation_pulse(int run=0, double th=20, double th2=30, int nSamples=4000)
{
   TFile *file0;
   if(run<10) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run0000%d.mid.lz4_waveforms.root",run),"READ");
   if(run>=10 && run<100) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run000%d.mid.lz4_waveforms.root",run),"READ");
   if(run>=100 && run<1000) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run00%d.mid.lz4_waveforms.root",run),"READ");
   if(run>=1000) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run0%d.mid.lz4_waveforms.root",run),"READ");
   TTree* ta = (TTree*)file0->Get("ta");
    
   TFile *file1; 
   if(run<10) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=10 && run<100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=100 && run<1000) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=1000) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0%d.mid.lz4_analyzed_1.root",run),"READ");
   TTree* ta_control = (TTree*)file1->Get("ta");

   double area[8];
   double area_tot=0;
   double t00;
   
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
   
   TH1F ** hP = new TH1F * [nChannels]; //Pulses histograms 
   
   for(int j = 0; j < nChannels; j++) 
   {
    hP[j] = new TH1F(Form("hP_%d",j),Form("hP_%d",j),nSamples,0,nSamples);
    ta->SetBranchAddress(Form("pulse%d",j),&hP[j]);
   }
   
   TH1F* hP_aux = new TH1F(Form("hP_aux"),Form("hP_aux"),nSamples,0,nSamples);
   
   TH1F * mean_waveform = new TH1F(Form("mean_%d",run),Form("mean_%d_%f",run,th),nSamples,0,nSamples);
   
   double bsl_provisional=0, rms_provisional;
   int nPtsBsl=500;
   int sIni=0;
   
   double selected_entries=0;
   
   int init_time;
   
   for(int i=0; i<ta->GetEntries(); i++)
    { 
      ta->GetEntry(i);
      ta_control->GetEntry(i);
      
      if(i%10000==0) cout << i << endl;
      
      for(int j=0; j<nChannels; j++)
      {
        area_tot+=area[j];
      }
      
      //if(area_tot>th && area_tot<th2 && t00>950 && t00<940)
      if(area_tot>th && area_tot<th2)
      //if(area_tot>th)
      {
          for(int j=0; j<1; j++)
          {
           bsl_provisional=bsl_simple(hP[j],nPtsBsl,sIni,&rms_provisional,nSamples);
           
           for(int k=1; k<=nSamples; k++) hP[j]->SetBinContent(k,(hP[j]->GetBinContent(k)-bsl_provisional));
           
           /*for(int k=hP[j]->GetMaximumBin(); k>0; k--) 
            {if(hP[j]->GetBinContent(k)<hP[j]->GetMaximum()*0.2){init_time=k; break;}}
           
           for(int k=1; k<=nSamples; k++) hP_aux->SetBinContent(nSamples*0.2-init_time+k,hP[j]->GetBinContent(k));*/
           
           mean_waveform->Add(hP[j],1.0);
           
           hP_aux->Reset();
          } 
          
          selected_entries++;
      } 
      
      area_tot=0;  
    }
   
   mean_waveform->Scale(1./double(selected_entries));
   
   //bsl_provisional=bsl_simple(mean_waveform,nPtsBsl,sIni,&rms_provisional,nSamples);
   
   //for(int k=1; k<=nSamples; k++) mean_waveform->SetBinContent(k,(mean_waveform->GetBinContent(k)-bsl_provisional));
   
   mean_waveform->Draw("HIST");
   
    TFile * foutroot = new TFile(Form("Scintillation_mean_pulses_OCTOPUS.root"),"UPDATE");
    mean_waveform->Write();
    foutroot->Close();
 }
 
 double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples)
{
  double bslTot=0;
  double bslTotF=0;
  
  double rmsTot=0;
  double rmsTotF=0;
  
  int k;

  for(k = 0; k < nPtsBsl; k++) bslTot += hAux->GetBinContent(k+1)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) bslTotF += hAux->GetBinContent(nSamples-nPtsBsl+k)/double(nPtsBsl);
  
  for(k = 0; k < nPtsBsl; k++) rmsTot += (hAux->GetBinContent(k+1)-bslTot)*(hAux->GetBinContent(k+1)-bslTot)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) rmsTotF += (hAux->GetBinContent(nSamples-nPtsBsl+k)-bslTotF)*(hAux->GetBinContent(nSamples-nPtsBsl+k)-bslTotF)/double(nPtsBsl);
  
  rmsTot = sqrt(rmsTot);
  rmsTotF = sqrt(rmsTotF);

  if(rmsTot < rmsTotF) {*rms=rmsTot; return bslTot;}
  else { *rms=rmsTotF; return bslTotF;}
}
