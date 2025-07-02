void area_tot_hist(int run=0)
{    
   TFile *file1; 
   if(run<10) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=10 && run<100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
   if(run>=100) file1 = TFile::Open(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
   TTree* ta = (TTree*)file1->Get("ta");

   TH1F* area_tot = new TH1F(Form("area_tot_%d",run),Form("area_tot_%d",run),10000,0,10000000);
   
   ta->Draw(Form("areaFix0+areaFix1+areaFix2+areaFix3+areaFix4+areaFix5+areaFix6+areaFix7>>area_tot_%d",run),"","");
   
    TFile * foutroot = new TFile(Form("Ba_spectra_Cycle2.root"),"UPDATE");
    area_tot->Write();
    foutroot->Close();
 }
