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

void phe_calibration (TFile *frun, int run, int ch=0, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0, int Height_nbins=300, int area_nbins=350);
void areaFix_cal (TFile *frun, int run, int ch=0, int Area_nbins=300, double* phe_area=0, double* phe_sigma=0, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0, int plot_control=0);

void SiPM_calibration (int run, int Height_nbins=300, int area_nbins=350)
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
    
    if(run<10) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=10 && run<100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
    
    int ch=0;
    
    TGraphErrors* channels_cal = new TGraphErrors();
    
    double cal, cal_err, res, res_err;
    double phe_height[8], phe_sigma[8];
    
    //
    //Calibraci髇 en 醨ea vs altura
    //
    
    TCanvas *c1 = new TCanvas("c1","The Graph example",10,10,1200,700);
    
    c1->Divide(4,2);
    
    FILE *f;
    f=fopen(Form("./Plots/Calibration_analysis/Channel_calibration_%d.txt",run),"w");
    
    FILE *f_res;
    f_res=fopen(Form("./Plots/Calibration_analysis/Channel_resolution_%d.txt",run),"w");
    
    for(ch=0; ch<=7; ch++)
    { 
      c1->cd(ch+1);
      
      phe_calibration(frun, run, ch,&cal,&cal_err,&res,&res_err,Height_nbins,area_nbins);
      
      channels_cal->SetPoint(ch,ch,cal);
      channels_cal->SetPointError(ch,0,cal_err);
      
      fprintf(f,"%d %f %f\n",ch, cal, cal_err);
      fprintf(f_res,"%d %f %f\n",ch, res, res_err);
    }  
    
    c1->SaveAs(Form("./Plots/Calibration_analysis/Finger_fit_all_channels_Run_%d.png",run));
     
    fclose(f);
    fclose(f_res);
    
    channels_cal->SetMarkerStyle(8);
    
    channels_cal->SetLineStyle(9);
    
    channels_cal->SetLineWidth(2);
    
    TCanvas *c2 = new TCanvas("c2","The Graph example",10,10,700,500);
    
    channels_cal->SetTitle(Form("Calibration Run %d",run));
    
    channels_cal->GetXaxis()->SetTitle("Channel");
    
    channels_cal->GetYaxis()->SetTitle("Calibration (a.u)");
    
    channels_cal->GetYaxis()->SetRangeUser(0,2000);
    
    channels_cal->Draw("APL");
    
    c2->SetGrid();
    
    c2->SaveAs(Form("./Plots/Calibration_analysis/Calibration_all_channels_Run_%d.png",run));
    
    //
    //Calibraci髇 y ajuste en area
    //
    
    TCanvas *c3 = new TCanvas("c3","The Graph example",10,10,1200,700);
    
    c3->Divide(4,2);
    
    f=fopen(Form("./Plots/Calibration_analysis/Channel_calibration_areaFix_%d.txt",run),"w");
    
    f_res=fopen(Form("./Plots/Calibration_analysis/Channel_resolution_areaFix_%d.txt",run),"w");
    
    for(ch=0; ch<=7; ch++)
    { 
      c3->cd(ch+1)->SetGrid();
      
      areaFix_cal (frun, run, ch, area_nbins, &phe_height[ch], &phe_sigma[ch],&cal,&cal_err,&res,&res_err);
      
      fprintf(f,"%d %f %f\n",ch, cal, cal_err);
      fprintf(f_res,"%d %f %f\n",ch,res, res_err);
    }  
    
    c3->SaveAs(Form("./Plots/Calibration_analysis/AreaFix_channels_fit_%d.png",run));
    
    for(ch=0; ch<=7; ch++)
    { 
      c3->cd(ch+1)->SetGrid();
      
      areaFix_cal (frun, run, ch, area_nbins, &phe_height[ch], &phe_sigma[ch],&cal,&cal_err,&res,&res_err,1);
      
      fprintf(f,"%d %f %f\n",ch, cal, cal_err);
      fprintf(f_res,"%d %f %f\n",ch,res, res_err);
    }  
    
    c3->SaveAs(Form("./Plots/Calibration_analysis/AreaFix_channels_cal_%d.png",run));
    
    cout << phe_height << " " << phe_sigma << endl;
    
    fclose(f);
    fclose(f_res);
    
}

void phe_calibration (TFile *frun, int run, int ch=0, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0, int Height_nbins=300, int area_nbins=350)
{
    gStyle->SetImageScaling(3.);
    gStyle->SetHistLineWidth(3);
    gStyle->SetTitleFont(102,"XYZ");
    gStyle->SetTitleFont(102,"");
    gStyle->SetTitleSize(0.07,"");
    gStyle->SetLegendFont(102);
    gStyle->SetLabelFont(102,"XYZ");
    gStyle->SetTextFont(102);
    

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
    
    Height_max = 700;

    TH1F *Height = new TH1F("Height", "Height",Height_nbins,0,Height_max);

    pulses->Draw(Form("high%d >> Height",ch),Form("high%d>0",ch),"");
    
    //
    //Se buscan los picos
    //
    TSpectrum *s = new TSpectrum(1000);
    int nfound;
    nfound = s->Search(Height,1.7,"",0.01);
    double* x_peaks = s->GetPositionX();
    
    std::sort(x_peaks, x_peaks+nfound);
    
    for(int i=0; i<nfound; i++) cout << x_peaks[i] << endl;
    
    if(nfound>4) nfound=4;   
    
    //Ajustes gaussianos 1D
    TF1 **gaus = new TF1*[nfound];
    
    int i, j;
     
     for (i=0;i<nfound;i++)
        {
          gaus[i] = new TF1(Form("gaus%d",i+1),"gaus",Height_min,Height_max);
          gaus[i]->SetParameter(0,Height->GetBinContent(x_peaks[i]/double(Height_max-Height_min)*double(Height_nbins)));
          gaus[i]->SetParameter(1,x_peaks[i]);
          if(i==0) gaus[i]->SetParameter(2,10);
          else gaus[i]->SetParameter(2,sqrt(x_peaks[i]));
          if(i==0) Height->Fit(gaus[i],"0","",x_peaks[i]-100,x_peaks[i]+10);
          else Height->Fit(gaus[i],"0","",x_peaks[i]-sqrt(x_peaks[i])*2,x_peaks[i]+sqrt(x_peaks[i])*2);
        }
     
     TF1* gaus_tot;
     
     if(nfound==1) gaus_tot=new TF1("gaus_tot","gaus(0)",Height_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==2) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)",Height_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==3) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)",Height_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==4) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",Height_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     
     for(i=0; i<nfound; i++)
       for(j=0; j<3; j++)
         gaus_tot->SetParameter(i*3+j,gaus[i]->GetParameter(j));
     
     Height->Fit(gaus_tot,"0","",Height_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     
      for(i=0; i<nfound; i++)
      {
       for(j=0; j<3; j++)
         gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
         gaus[i]->SetParError(j,gaus_tot->GetParError(j+i*3));
      }

    double par[12];
    double par_err[12];

    gaus_tot->GetParameters(&par[0]);
    
    for(i=0; i<12; i++)
      par_err[i]=gaus_tot->GetParError(i);

    TF1* gaus_0 = new TF1("gaus 0","gaus",0,Height_max);
    TF1* gaus_1 = new TF1("gaus 1","gaus",0,Height_max);
    TF1* gaus_2 = new TF1("gaus 2","gaus",0,Height_max);
    TF1* gaus_3 = new TF1("gaus 3","gaus",0,Height_max);
    
    gaus_0->SetParameters(&par[0]);
    gaus_1->SetParameters(&par[3]);
    gaus_2->SetParameters(&par[6]);
    gaus_3->SetParameters(&par[9]);
    
    gaus_0->SetParErrors(&par_err[0]);
    gaus_1->SetParErrors(&par_err[3]);
    gaus_2->SetParErrors(&par_err[6]);
    gaus_3->SetParErrors(&par_err[9]);
    
    //
    //Histograma de 醨ea
    //
    TH1F *area_aux = new TH1F("area_aux","area_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> area_aux",ch),Form("high%d>0",ch),"");

    gStyle->SetOptStat(0);

    int area_max = area_aux->GetXaxis()->GetXmax();
    int area_min = area_aux->GetXaxis()->GetXmin();
    //int area_nbins = 350;

    double area_bin_size = double(area_max-area_min)/area_nbins;

    int area_final_bin = area_aux->FindLastBinAbove(10);
    int area_initial_bin = area_aux->FindFirstBinAbove(10);

    area_max = area_final_bin*area_bin_size+area_min;
    area_min = area_initial_bin*area_bin_size+area_min;
    
    area_max=8000;
    
    //
    //Ajustes 2D
    // 
    TH2F* finger = new TH2F("finger",Form("Run %d [Channel %d]",run,ch),Height_nbins,0,Height_max,area_nbins,area_min,area_max);
    TF2 * finger_cal_0 = new TF2("finger_cal_0","bigaus",0,Height_max,area_min,area_max);
    TF2 * finger_cal_1 = new TF2("finger_cal_1","bigaus",0,Height_max,area_min,area_max);
    TF2 * finger_cal_2 = new TF2("finger_cal_2","bigaus",0,Height_max,area_min,area_max);
    
    TF2 * finger_cal_tot = new TF2("finger_cal_tot","bigaus(0)+bigaus(6)+bigaus(12)",0,gaus_2->GetParameter(1)+gaus_2->GetParameter(2)*4,area_min,area_max);
    
    //finger->GetXaxis()->SetRangeUser(0,300);
    //finger->GetYaxis()->SetRangeUser(-1000,8000);

    pulses->Draw(Form("areaFix%d:high%d>>finger",ch,ch),"","");
    
    finger->GetXaxis()->SetRangeUser(0,300);
    finger->GetYaxis()->SetRangeUser(-1000,5000);

    finger_cal_0->SetParameter(0,gaus_0->GetParameter(0));
    finger_cal_0->SetParameter(1,gaus_0->GetParameter(1));
    finger_cal_0->SetParameter(2,gaus_0->GetParameter(2));
    //finger_cal_0->SetParLimits(1,0+10,Height_cut[0]-10);
    finger_cal_0->SetParLimits(0,0,100000);
    
    finger_cal_1->SetParameter(0,gaus_1->GetParameter(0));
    finger_cal_1->SetParameter(1,gaus_1->GetParameter(1));
    finger_cal_1->SetParameter(2,gaus_1->GetParameter(2));
    //finger_cal_1->SetParLimits(1,Height_cut[0]+10,Height_cut[1]-10);
    finger_cal_1->SetParLimits(0,0,100000);
    
    finger_cal_2->SetParameter(0,gaus_2->GetParameter(0));
    finger_cal_2->SetParameter(1,gaus_2->GetParameter(1));
    finger_cal_2->SetParameter(2,gaus_2->GetParameter(2));
    //finger_cal_2->SetParLimits(1,Height_cut[1]+10,Height_cut[2]-10);
    finger_cal_2->SetParLimits(0,0,100000);
    
    finger->Print();
    finger->Fit("finger_cal_0","0","",0,gaus_0->GetParameter(1)+gaus_0->GetParameter(2)*2);
    
    finger->Fit("finger_cal_1","0","",gaus_1->GetParameter(1)-gaus_1->GetParameter(2)*2,gaus_1->GetParameter(1)+gaus_1->GetParameter(2)*2);
    
    finger->Fit("finger_cal_2","0","",gaus_2->GetParameter(1)-gaus_2->GetParameter(2)*2,gaus_2->GetParameter(1)+gaus_2->GetParameter(2)*2);

    for(i=0; i<6; i++)
      finger_cal_tot->SetParameter(i,finger_cal_0->GetParameter(i));

    for(i=6; i<12; i++)
      finger_cal_tot->SetParameter(i,finger_cal_1->GetParameter(i-6));
    
    for(i=12; i<18; i++)
      finger_cal_tot->SetParameter(i,finger_cal_2->GetParameter(i-12));

    finger->Fit("finger_cal_tot","0","",0,gaus_2->GetParameter(1)+gaus_2->GetParameter(2)*2);
      
     TF2 * finger_cal_plot_0 = new TF2("finger_cal_plot_0","bigaus",-100,500,-1000,8000);  
     TF2 * finger_cal_plot_1 = new TF2("finger_cal_plot_1","bigaus",0,500,0,8000);
     TF2 * finger_cal_plot_2 = new TF2("finger_cal_plot_2","bigaus",0,500,0,8000);

    for(i=0; i<6; i++)
      finger_cal_plot_0->SetParameter(i,finger_cal_tot->GetParameter(i));

    for(i=0; i<6; i++)
      finger_cal_plot_1->SetParameter(i,finger_cal_tot->GetParameter(i+6));
    
    for(i=0; i<6; i++)
      finger_cal_plot_2->SetParameter(i,finger_cal_tot->GetParameter(i+12));

    finger_cal_plot_1->SetLineColor(2);
    finger_cal_plot_2->SetLineColor(3);
    finger_cal_plot_0->SetLineColor(5);
    finger_cal_plot_1->SetFillColor(2);
    finger_cal_plot_2->SetFillColor(3);
    finger_cal_tot->SetLineColor(1);
    finger_cal_tot->SetLineWidth(3);


    finger->GetXaxis()->SetTitle("Height (a.u)");
    finger->GetYaxis()->SetTitle("Area (a.u)");

    finger->Draw("colz");

    finger_cal_tot->SetNpx(500);
    finger_cal_plot_0->SetNpx(500);
    finger_cal_plot_1->SetNpx(500);
    finger_cal_plot_2->SetNpx(500);
    
    finger_cal_plot_0->SetNpy(500);
    finger_cal_plot_1->SetNpy(500);
    finger_cal_plot_2->SetNpy(500);
    
    finger_cal_plot_0->SetLineWidth(3);
    finger_cal_plot_1->SetLineWidth(3);
    finger_cal_plot_2->SetLineWidth(3);

    finger_cal_plot_0->SetContour(2);
    finger_cal_plot_1->SetContour(2);
    finger_cal_plot_2->SetContour(2);
    finger_cal_tot->SetContour(2);
    
    finger_cal_plot_0->Draw("cont3 && same");
    finger_cal_plot_1->Draw("cont3 && same");
    finger_cal_plot_2->Draw("cont3 && same");

    TLegend *legend4 = new TLegend(0.55,0.15,0.89,0.30); //Tama駉 de la leyenda (x0,y0,x1,y1)
    legend4->SetTextSize(0.04); //Tama駉 del texto
    legend4->SetFillColor(kWhite); //Fondo de la leyenda

    legend4->AddEntry(finger_cal_plot_1,Form("1 phe"),"l");
    legend4->AddEntry(finger_cal_plot_2,Form("2 phe"),"l");
    legend4->AddEntry(finger_cal_tot,Form("#chi^{2}/ndf=%0.2f",finger_cal_tot->GetChisquare()/finger_cal_tot->GetNDF()),"l");

    legend4->Draw("same");
    
    //
    //Calibraci髇 phe
    //
    TGraphErrors* area_cal = new TGraphErrors();
    TF1 *area_cal_fit = new TF1("area_cal_fit","[0]+[1]*x");
    
    area_cal->SetPoint(0,0.,finger_cal_tot->GetParameter(3));
    area_cal->SetPointError(0,0.,finger_cal_tot->GetParError(3));
    
    area_cal->SetPoint(1,1.,finger_cal_tot->GetParameter(9));
    area_cal->SetPointError(1,0.,finger_cal_tot->GetParError(9));

    area_cal->SetPoint(2,2.,finger_cal_tot->GetParameter(15));
    area_cal->SetPointError(2,0.,finger_cal_tot->GetParError(15));

    area_cal->SetTitle("Area 2D calibration");

    area_cal->Fit(area_cal_fit);
    
     *cal=area_cal_fit->GetParameter(1);
     *cal_err=area_cal_fit->GetParError(1);
     
     *res=finger_cal_tot->GetParameter(10)/finger_cal_tot->GetParameter(9)*100;
     *res_err=100*sqrt(finger_cal_tot->GetParError(10)*finger_cal_tot->GetParError(10)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)+finger_cal_tot->GetParError(9)*finger_cal_tot->GetParError(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)*finger_cal_tot->GetParameter(10)*finger_cal_tot->GetParameter(10));
}

void areaFix_cal (TFile *frun, int run, int ch=0, int Area_nbins=300, double* phe_area, double* phe_sigma, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0, int plot_control=0)
{   
    TTree * pulses = (TTree*)frun->Get("ta");
    
    //
    //Definimos propiedades histograma de altura
    //

    TH1F *Area_aux = new TH1F("Area_aux","Area_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> Area_aux",ch),"","");

    int Area_max = Area_aux->GetXaxis()->GetXmax();
    int Area_min = Area_aux->GetXaxis()->GetXmin();

    double Area_bin_size = (Area_max-Area_min)/10000.;

    int Area_final_bin = Area_aux->FindLastBinAbove(1);
    int Area_initial_bin = Area_aux->FindFirstBinAbove(1);

    Area_max = Area_final_bin*Area_bin_size+Area_min;
    Area_min = Area_initial_bin*Area_bin_size+Area_min;

    TH1F *Area = new TH1F("Area", "Area",Area_nbins,Area_min,Area_max);

    pulses->Draw(Form("areaFix%d >> Area",ch),"","");
    
    //
    //Se buscan los picos
    //
    TSpectrum *s = new TSpectrum(1000);
    int nfound;
    nfound = s->Search(Area,1.7,"",0.01);
    double* x_peaks = s->GetPositionX();
    
    std::sort(x_peaks, x_peaks+nfound);
    
    for(int i=0; i<nfound; i++) cout << x_peaks[i] << endl;
    
    if(nfound>4) nfound=4;
    
    //Ajustes gaussianos 1D
    
    TF1 **gaus = new TF1*[nfound];
    
    int i, j;
     
     for (i=0;i<nfound;i++)
        {
          gaus[i] = new TF1(Form("gaus%d",i+1),"gaus",Area_min,Area_max);
          gaus[i]->SetParameter(0,Area->GetBinContent(x_peaks[i]/double(Area_max-Area_min)*double(Area_nbins)));
          gaus[i]->SetParameter(1,x_peaks[i]);
          if(i==0) gaus[i]->SetParameter(2,50);
          else gaus[i]->SetParameter(2,sqrt(x_peaks[i]));
          if(i==0) Area->Fit(gaus[i],"0","",x_peaks[i]-100,x_peaks[i]+100);
          else Area->Fit(gaus[i],"0","",x_peaks[i]-sqrt(x_peaks[i])*2,x_peaks[i]+sqrt(x_peaks[i])*2);
        }
     
     TF1* gaus_tot;
     
     if(nfound==1) gaus_tot=new TF1("gaus_tot","gaus(0)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==2) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==3) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     if(nfound==4) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     
     for(i=0; i<nfound; i++)
       for(j=0; j<3; j++)
         gaus_tot->SetParameter(i*3+j,gaus[i]->GetParameter(j));
     
     Area->Fit(gaus_tot,"0","",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*4);
     
      for(i=0; i<nfound; i++)
      {
       for(j=0; j<3; j++)
         gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
         gaus[i]->SetParError(j,gaus_tot->GetParError(j+i*3));
      }
    
    *phe_area = gaus[1]->GetParameter(1);
    *phe_sigma = gaus[1]->GetParameter(2);
    
    //
    //Calibraci髇 phe
    //
    TGraphErrors* area_cal = new TGraphErrors();
    TF1 *area_cal_fit = new TF1("area_cal_fit","[0]+[1]*x",0,10);
    
    for(i=0; i<nfound; i++)
    {
      area_cal->SetPoint(i,i,gaus[i]->GetParameter(1));
      area_cal->SetPointError(i,0,gaus[i]->GetParError(1));
    }
    
    area_cal->SetTitle("Area 1D calibration");

    area_cal->Fit(area_cal_fit,"0");
    
     *cal=area_cal_fit->GetParameter(1);
     *cal_err=area_cal_fit->GetParError(1);
     
     *res=gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)*100;
     *res_err=100*sqrt(gaus[1]->GetParError(2)*gaus[1]->GetParError(2)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)+gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)*gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)*gaus[1]->GetParError(1)*gaus[1]->GetParError(1));
    
    //
    //Plot
    //
    
    if(plot_control==0)
    {
    Area->SetLineColor(kBlack);
    
    Area->SetTitle(Form("Run %d [Channel %d]",run,ch));
    
    Area->GetXaxis()->SetTitle("Area (ADC)");
    
    Area->Draw("HIST");
    
    TLegend *legend2 = new TLegend(0.70,0.87-nfound*0.05,0.87,0.85); //Tama駉 de la leyenda (x0,y0,x1,y1)
    legend2->SetTextSize(0.04); //Tama駉 del texto
    legend2->SetFillColor(kWhite); //Fondo de la leyenda
    
    for(i=0; i<nfound; i++)
      { 
        gaus[i]->SetLineWidth(4.);
        if(i==0) gaus[i]->SetLineColor(kCyan+1);
        if(i==1) gaus[i]->SetLineColor(kRed+1);
        if(i==2) gaus[i]->SetLineColor(kGreen-3);
        if(i==3) gaus[i]->SetLineColor(kMagenta-7);
        
        legend2->AddEntry(gaus[i],Form("%d phe",i),"l");
        
        gaus[i]->Draw("same");
      }
      
      gaus_tot->SetLineColor(kRed-4);
      gaus_tot->SetLineWidth(4.);
      gaus_tot->Draw("same");
      
    legend2->Draw("same");
    }
    
    else
    {
       area_cal->SetTitle(Form("Run %d [Channel %d]",run,ch));
       area_cal->GetXaxis()->SetTitle("Nphe");
       area_cal->SetMarkerStyle(8);
       area_cal->SetLineWidth(4);
       area_cal->SetMarkerSize(1.0);
       
       area_cal_fit->SetLineStyle(9);
       area_cal_fit->SetLineWidth(4.);
       area_cal_fit->SetLineColor(kRed+1);
       
       area_cal->Draw("AP");
       area_cal_fit->Draw("same");
       
       TLegend *legend2 = new TLegend(0.13,0.80,0.70,0.85); //Tama駉 de la leyenda (x0,y0,x1,y1)
       legend2->SetTextSize(0.04); //Tama駉 del texto
       legend2->SetFillColor(kWhite); //Fondo de la leyenda
       legend2->AddEntry(area_cal_fit,Form("Cal = %.1f #pm %.1f",area_cal_fit->GetParameter(1),area_cal_fit->GetParError(1)),"l");
       
       legend2->Draw();
    }
}

