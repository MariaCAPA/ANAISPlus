//C,C++ Libraries
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iomanip>

// ROOT libraries
#include <TROOT.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TObject.h>
#include <TStyle.h>
#include <TBrowser.h>
#include <TApplication.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCutG.h>
#include <TString.h>
#include <TObjArray.h>
#include <TPaveLabel.h>

void heigh_cal (TFile *frun, int run, int ch=0, int Height_nbins=300, double* phe_height=0, double* phe_sigma=0, double* mean_number_phe=0, double* mean_number_phe_err=0);

void light_channel_distribution (int run, int run_f, int Height_nbins=300, int area_nbins=350)
{
    gStyle->SetImageScaling(3.);
    gStyle->SetHistLineWidth(3);
    gStyle->SetTitleFont(102,"XYZ");
    gStyle->SetTitleFont(102,"");
    gStyle->SetTitleSize(0.07,"");
    gStyle->SetLegendFont(102);
    gStyle->SetLabelFont(102,"XYZ");
    gStyle->SetTextFont(102);
    
    TFile *frun;
    
    if(run<10) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_anakyzed_1.root",run),"READ");
    if(run>=10 && run<100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
    
    int ch=0;
    
    double phe_height[8], phe_sigma[8];
    double mean_number_phe[8];
    double mean_number_phe_err[8];
    
    TCanvas *c3 = new TCanvas("c3","The Graph example",10,10,1200,700);
    
    gStyle->SetOptStat(0);
    
    c3->Divide(4,2);
    
    if(run==run_f)
    {
    for(ch=0; ch<=7; ch++)
        { 
          c3->cd(ch+1);

          heigh_cal (frun, run, ch, Height_nbins, &phe_height[ch], &phe_sigma[ch], &mean_number_phe[ch], &mean_number_phe_err[ch]);
        }  
    }
    
    else
    {
      for(int j=0; j<=6; j++)
      { 
        ch=j;
        
        if(j>=1) ch++;
        
        c3->cd(ch+1);
        
        if(run+j<10) frun = new TFile(Form("../../../analyzed/run0000%d.mid.lz4.root",run+j),"READ");
        if(run+j>=10 && run<100) frun = new TFile(Form("../../../analyzed/run000%d.mid.lz4.root",run+j),"READ");
        if(run+j>=100) frun = new TFile(Form("../../../analyzed/run00%d.mid.lz4.root",run+j),"READ");
          
        heigh_cal (frun, run+j, 0, Height_nbins, &phe_height[ch], &phe_sigma[ch], &mean_number_phe[ch], &mean_number_phe_err[ch]);
      }  
    }
    
    c3->SaveAs(Form("./Plots/Light_dsitribution/Nphe_channels_spectra_%d.png",run));
    
    //Light distribution plot
    TCanvas *c1 = new TCanvas("c1","The Graph example",200,200,700,700);
    
    c1->SetGrid();
    
    TH2F* light_distribution = new TH2F(Form("Light_distribution_%d",run),"H2",4,0,4,2,0,4);
    light_distribution->SetTitle("Channels light distribution");
    light_distribution->GetXaxis()->SetTitle("4x4 SiPMs");
    light_distribution->GetXaxis()->SetNdivisions(4);
    light_distribution->GetYaxis()->SetNdivisions(2);
    
    TGraphErrors* light_distribution_graph = new TGraphErrors(Form("Light_distribution_graph_%d",run));
    light_distribution_graph->SetMarkerStyle(8);
    light_distribution_graph->SetLineStyle(9);
    light_distribution_graph->SetLineWidth(2);
    
    light_distribution_graph->GetXaxis()->SetRangeUser(0,7);
    light_distribution_graph->GetXaxis()->SetTitle("Channel");
    light_distribution_graph->GetYaxis()->SetTitle("<nphe>");
    light_distribution_graph->SetTitle(Form("Light distribution [Run %d]",run));
    
    
    for(ch=0; ch<=7; ch++)
    { 
      if(ch<4) light_distribution->SetBinContent(ch+1,2,mean_number_phe[ch]);
      if(ch>=4) light_distribution->SetBinContent(8-ch,1,mean_number_phe[ch]);
         
      light_distribution_graph->SetPoint(ch,ch,mean_number_phe[ch]);
      light_distribution_graph->SetPointError(ch,0,mean_number_phe_err[ch]);
    }
    
    light_distribution->Draw("colz");
    
    c1->SaveAs(Form("./Plots/Light_distribution/Nphe_channels_distribution_%d.png",run));
    
    light_distribution_graph->Draw("APL");
    
    c1->SaveAs(Form("./Plots/Light_distribution/Nphe_channels_distribution_graph_%d.png",run));
    
    //Guardamos los datos de distribucion de luz en un fichero
    FILE *f_text;
    
    f_text=fopen(Form("./Plots/Light_distribution/Light_distribution_%d.txt",run),"a");
    
    for(ch=0; ch<=7; ch++)
    {
      fprintf(f_text,"%d %f %f\n",ch,mean_number_phe[ch],mean_number_phe_err[ch]);
    }
    
    fclose(f_text);     
}

void heigh_cal (TFile *frun, int run, int ch=0, int Height_nbins=300, double* phe_height=0, double* phe_sigma=0, double* mean_number_phe=0, double* mean_number_phe_err=0)
{   
    TTree * pulses = (TTree*)frun->Get("ta");
    
    //
    //Definimos propiedades histograma de altura
    //

    TH1F *Height_aux = new TH1F("Height_aux","Height_aux",10000,100,50);

    pulses->Draw(Form("high%d >> Height_aux",ch),Form("high%d>0",ch),"");

    int Height_max = Height_aux->GetXaxis()->GetXmax();
    int Height_min = Height_aux->GetXaxis()->GetXmin();

    double Height_bin_size = (Height_max-Height_min)/10000.;

    int Height_final_bin = Height_aux->FindLastBinAbove(1);
    int Height_initial_bin = Height_aux->FindFirstBinAbove(1);

    Height_max = Height_final_bin*Height_bin_size+Height_min;
    Height_min = Height_initial_bin*Height_bin_size+Height_min;

    TH1F *Height = new TH1F(Form("Height_%d_%d",run,ch), "Height",Height_nbins,0,Height_max);

    pulses->Draw(Form("high%d >> Height_%d_%d",ch,run,ch),Form("high%d>0",ch),"");
    
    //
    //Se buscan los picos
    //
    TSpectrum *s = new TSpectrum(1000);
    int nfound;
    nfound = s->Search(Height,1.5,"",0.01);
    double* x_peaks = s->GetPositionX();
    
    std::sort(x_peaks, x_peaks+nfound);
    
    for(int i=0; i<nfound; i++) cout << x_peaks[i] << endl;
    
    if(nfound>4) nfound=4;
    
    //Ajustes gaussianos 1D
    
    TF1 **gaus = new TF1*[nfound];
    
    int i, j;
     
     for (i=0;i<nfound;i++)
        {
          gaus[i] = new TF1(Form("gaus%d",i+1),"gaus",0,Height_max);
          gaus[i]->SetParameter(0,Height->GetBinContent(x_peaks[i]/double(Height_max)*double(Height_nbins)));
          gaus[i]->SetParameter(1,x_peaks[i]);
          gaus[i]->SetParameter(2,20);
          Height->Fit(gaus[i],"0","",x_peaks[i]-10,x_peaks[i]+10);
        }
     
     TF1* gaus_tot;
     
     if(nfound==1) gaus_tot=new TF1("gaus_tot","gaus(0)",0,Height_max);
     if(nfound==2) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)",0,Height_max);
     if(nfound==3) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)",0,Height_max);
     if(nfound==4) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",0,Height_max);
     
     for(i=0; i<nfound; i++)
       for(j=0; j<3; j++)
         gaus_tot->SetParameter(i*3+j,gaus[i]->GetParameter(j));
     
     Height->Fit(gaus_tot,"0","",0,x_peaks[nfound-1]+30);
     
      for(i=0; i<nfound; i++)
      {
       for(j=0; j<3; j++)
         gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
      }
    
    //Height calibration
    
    TGraphErrors* height_cal = new TGraphErrors();
    TF1 *height_cal_fit = new TF1("height_cal_fit","[0]+[1]*x");
    
    height_cal->SetPoint(0,1.,gaus[1]->GetParameter(1));
    height_cal->SetPointError(0,0.,gaus[1]->GetParameter(2));

    height_cal->SetPoint(1,2.,gaus[2]->GetParameter(1));
    height_cal->SetPointError(1,0.,gaus[2]->GetParameter(2));

    height_cal->SetTitle("Height calibration");

    height_cal->Fit(height_cal_fit);
    
    double cal, cal_err;
    
    cal=height_cal_fit->GetParameter(1);
    cal_err=height_cal_fit->GetParError(1);
    
    //Nphe spectra
    Height->GetXaxis()->SetLimits(0-height_cal_fit->GetParameter(0)/cal, Height_max/cal-height_cal_fit->GetParameter(0)/cal);
    
    //Mean number of phe
    *mean_number_phe=(Height->GetMean()-height_cal_fit->GetParameter(0))/cal;
    *mean_number_phe_err=(Height->GetMeanError())/cal;
    
    //
    //Plot
    //
    Height->SetLineColor(kBlack);
    
    Height->GetXaxis()->SetRangeUser(0,10);
    
    Height->SetTitle(Form("Run %d [Channel %d]",run,ch));
    
    Height->GetXaxis()->SetTitle("Nphe");
    
    Height->Draw("HIST");
    
    TLegend *legend2 = new TLegend(0.70,0.87-nfound*0.05,0.87,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend2->SetTextSize(0.04); //Tamaño del texto
    legend2->SetFillColor(kWhite); //Fondo de la leyenda
    
    /*for(i=0; i<nfound; i++)
      { 
        gaus[i]->SetLineWidth(4.);
        if(i==0) gaus[i]->SetLineColor(kCyan+1);
        if(i==1) gaus[i]->SetLineColor(kRed+1);
        if(i==2) gaus[i]->SetLineColor(kGreen-3);
        if(i==3) gaus[i]->SetLineColor(kMagenta-7);
        
        legend2->AddEntry(gaus[i],Form("%d phe",i),"l");
        
        gaus[i]->Draw("same");
      }
      
    legend2->Draw("same");*/
}
    
