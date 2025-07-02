double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples);

void mean_phe(int run=0, int nSamples=250)
{
   TFile *file0;
   if(run<10) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run0000%d.mid.lz4_waveforms.root",run),"READ");
   if(run>=10 && run<100) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run000%d.mid.lz4_waveforms.root",run),"READ");
   if(run>=100) file0 = TFile::Open(Form("../analyzed/ANAISplusWaveforms/run00%d.mid.lz4_waveforms.root",run),"READ");
   TTree* ta = (TTree*)file0->Get("ta");
    
   TFile *file1; 
   if(run<10) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=10 && run<100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
   TTree* ta_control = (TTree*)file1->Get("ta");

   double high[8];
   double area_tot=0;
   double t00;
   
   double th[8];
   double th2[8];
   
   TH1F* high_spectrum = new TH1F("high_spectrum","high_spectrum",600,0,50000);
   
   ta_control->SetBranchAddress("areaFix0",&high[0]);
/*   ta_control->SetBranchAddress("high1",&high[1]);
   ta_control->SetBranchAddress("high2",&high[2]);
   ta_control->SetBranchAddress("high3",&high[3]);
   ta_control->SetBranchAddress("high4",&high[4]);
   ta_control->SetBranchAddress("high5",&high[5]);
   ta_control->SetBranchAddress("high6",&high[6]);
   ta_control->SetBranchAddress("high7",&high[7]);
  */ 
   int nChannels=1;
   
   TH1F ** hP = new TH1F * [nChannels]; //Pulses histograms 
   
   for(int j = 0; j < nChannels; j++) 
   {
    hP[j] = new TH1F(Form("hP_%d",j),Form("hP_%d",j),nSamples,0,nSamples);
    ta->SetBranchAddress(Form("pulse%d",j),&hP[j]);
   }
   
   TH1F* hP_aux = new TH1F(Form("hP_aux"),Form("hP_aux"),nSamples,0,nSamples);
   
   TSpectrum *s = new TSpectrum();
   
   double* xpeaks;
   //int index;
   int nfound;
   
   for(int j = 0; j < nChannels; j++) 
   {
    ta_control->Draw(Form("areaFix%d>>high_spectrum",j));
    
    nfound = s->Search(high_spectrum);
    
    xpeaks = s->GetPositionX();
    
    int index[nfound];
    
    TMath::Sort(nfound,xpeaks,index,0);
    
    th[j]=xpeaks[index[1]]-1000;
    th2[j]=xpeaks[index[1]]+1000;
    
    cout << th[j] << endl;
    cout << th2[j] << endl;
    
    cout << "" << endl;
   }
   
   TH1F * mean_waveform = new TH1F(Form("mean_%d",run),Form("mean_%d",run),nSamples,0,nSamples);
   
   double bsl_provisional=0, rms_provisional;
   int nPtsBsl=50;
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
          if(high[j]>th[j] && high[j]<th2[j])
          {
           bsl_provisional=bsl_simple(hP[j],nPtsBsl,sIni,&rms_provisional,nSamples);
           
           for(int k=1; k<=nSamples; k++) hP[j]->SetBinContent(k,(hP[j]->GetBinContent(k)-bsl_provisional));
           
           for(int k=hP[j]->GetMaximumBin(); k>0; k--) 
            if(hP[j]->GetBinContent(k)<hP[j]->GetMaximum()*0.2){init_time=k; break;}
           
           for(int k=1; k<=nSamples; k++) hP_aux->SetBinContent(k,hP[j]->GetBinContent(k)*(-1));
           //for(int k=1; k<=nSamples; k++) hP_aux->SetBinContent(100-init_time+k,hP[j]->GetBinContent(k));
           //for(int k=1; k<=nSamples; k++) hP_aux->SetBinContent(100-hP[j]->GetMaximumBin()+k,hP[j]->GetBinContent(k));
           
           mean_waveform->Add(hP_aux,1.0);
           
           hP_aux->Reset();
           
           selected_entries++;
          } 
      } 

    }
   
   mean_waveform->Scale(1./(double(selected_entries)));
   
   //bsl_provisional=bsl_simple(mean_waveform,nPtsBsl,sIni,&rms_provisional,nSamples);
   
   //for(int k=1; k<=nSamples; k++) mean_waveform->SetBinContent(k,(mean_waveform->GetBinContent(k)-bsl_provisional));
   
   mean_waveform->Draw("HIST");
   
    TFile * foutroot = new TFile(Form("phe_mean_pulses_FIRST_PROTOTYPE.root"),"UPDATE");
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
