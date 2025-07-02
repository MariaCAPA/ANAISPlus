void DCR_time_histogram(int run=0, double th=20)
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

   double height_value;
   ta_control->SetBranchAddress("high0",&height_value);
   
   TH1F* time_hist = new TH1F(Form("time_hist_%d",run),"time_hist",500000,0,2.0);
   TH1F* time_hist2 = new TH1F(Form("time_hist2_%d",run),"time_hist",500000,0,2.0);
   TH1F* time_hist3 = new TH1F(Form("time_hist3_%d",run),"time_hist",500000,0,2.0);
   TH1F* time_hist4 = new TH1F(Form("time_hist4_%d",run),"time_hist",500000,0,2.0);
   
   double time_i, time_var;
   
   ta->SetBranchAddress("timeFile",&time_var);
   
   for(int i=0; i<ta->GetEntries(); i++)
   {
       ta->GetEntry(i);
       ta_control->GetEntry(i);

      if(height_value>th)
       { 
       time_hist->Fill((time_var-time_i)*1.e-9);
       time_i=time_var;
       }
   }
   /*
   for(int i=0; i<ta->GetEntries(); i++)
   {
       ta->GetEntry(i); 
       ta_control->GetEntry(i);
       
       if(height_value>th)
       { 
       time_i=time_var;
       ta->GetEntry(i+2);
       time_hist2->Fill((time_var-time_i)*1.e-9);
       }
   }
   
   for(int i=0; i<ta->GetEntries(); i++)
   {
       ta->GetEntry(i); 
       ta_control->GetEntry(i);
       
       if(height_value>25)
       { 
       time_i=time_var;
       ta->GetEntry(i+3);
       time_hist3->Fill((time_var-time_i)*1.e-9);
       }
   }
   
   for(int i=0; i<ta->GetEntries(); i++)
   {
       ta->GetEntry(i); 
       ta_control->GetEntry(i);
      
       if(height_value>th)
       { 
       time_i=time_var;
       ta->GetEntry(i+4);
       time_hist4->Fill((time_var-time_i)*1.e-9);
       }
   }
   */
   TFile *myfile = TFile::Open("DCR_time_histograms.root","UPDATE");
   time_hist->Write();
   //time_hist2->Write();
   //time_hist3->Write();
   //time_hist4->Write();
   
   myfile->Close();
}
