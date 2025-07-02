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
void areaFix_cal (TFile *frun, int run, int ch=0, int Height_nbins=300, double* phe_height=0, double* phe_sigma=0, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0);
void fwhm_distribution (TFile *frun, int run, int ch=0, int Height_nbins=300, double phe_height=0, double phe_sigma=0, TProfile* one_phe=0);

void calibration_areaFix (int run, int Height_nbins=300, int area_nbins=350)
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
    
    TProfile *one_phe = new TProfile(Form("run%d",run),Form("run%d",run),250,0,250*4);
    
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
    //Calibraci髇 y ajuste en altura
    //
    
    TCanvas *c3 = new TCanvas("c3","The Graph example",10,10,1200,700);
    
    c3->Divide(4,2);
    
    f=fopen(Form("./Plots/Calibration_analysis/Channel_calibration_areaFix_%d.txt",run),"w");
    
    f_res=fopen(Form("./Plots/Calibration_analysis/Channel_resolution_areaFix_%d.txt",run),"w");
    
    for(ch=0; ch<=7; ch++)
    { 
      c3->cd(ch+1);
      
      areaFix_cal (frun, run, ch, area_nbins, &phe_height[ch], &phe_sigma[ch],&cal,&cal_err,&res,&res_err);
      
      fprintf(f,"%d %f %f\n",ch, cal, cal_err);
      fprintf(f_res,"%d %f %f\n",ch,res, res_err);
    }  
    
    c3->SaveAs(Form("./Plots/Calibration_analysis/AreaFix_channels_fit_%d.png",run));
    
    cout << phe_height << " " << phe_sigma << endl;
    
    fclose(f);
    fclose(f_res);
  
    //
    //Estudio FWHM
    //
    
    /*TCanvas *c4 = new TCanvas("c4","The Graph example",10,10,1200,700);
    TCanvas *c5 = new TCanvas("c5","The Graph example",10,10,1200,700);
    
    c4->Divide(4,2);
    c5->Divide(4,2);
    
    for(ch=0; ch<=7; ch++)
    { 
      c4->cd(ch+1);
      
      one_phe = new TProfile(Form("run%d_%d",run,ch),Form("Run %d [Channel %d]",run,ch),250,0,250*4);
      
      one_phe->Reset();
      
      one_phe->SetMarkerStyle(8);
    
      one_phe->SetMarkerSize(0.3);
      
      one_phe->SetMarkerColor(2);
      
      one_phe->SetLineColor(2);
      
      fwhm_distribution (frun, run, ch, Height_nbins, phe_height[ch], phe_sigma[ch], one_phe);
      
      c5->cd(ch+1);
      
      one_phe->Draw();
    } 
    
    c4->SaveAs(Form("./Plots/Calibration_analysis/FWHM_channels_%d.png",run));
    c5->SaveAs(Form("./Plots/Calibration_analysis/Average_phe_%d.png",run));*/

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
    //Se buscan cortes en altura
    //
    int i, j, k=0;
    int Height_cut[3];

    double media=0, control=1, media_ant=0, control_ant=1;

    //Height->GetXaxis()->SetRangeUser(0,500);

    for(i=0; i<3; i++)
    	Height_cut[i]=0;

    for(i=0; i<Height_nbins/5.; i++)
    {
    	media_ant=media;
    	media=0;
    	control_ant=control;

    	for(j=1; j<=5; j++)
    		 media+=Height->GetBinContent(i*5+j);

    	media=media/5.;

    	if(media>media_ant)
    		control=1;

    	if(media<media_ant)
    		control=-1;

    	if(control-control_ant==2 && k<3)
     //if(control!=control_ant && control==1 && k<3)
    		{
    		 Height_cut[k]=((i-1)*5+2.5)*(double(Height_max)/double(Height_nbins)); //Se pasa de bin a high0
    		 k++;
    		}
    }

    if(Height_cut[2]-Height_cut[1]>1.2*(Height_cut[1]-Height_cut[0]) || Height_cut[2]-Height_cut[1]<0.6*(Height_cut[1]-Height_cut[0])) //Ajuste por si se va mucho el limite superior de 2 phe
	      Height_cut[2]=Height_cut[1]+Height_cut[1]-Height_cut[0];

    /*if(Height->GetBinContent(1)==0)
    {
      Height_cut[2]=Height_cut[1];
      Height_cut[1]=Height_cut[0];
      Height_cut[0]=0;
    }*/

	if(run==321)
  {
    Height_cut[2]=Height_cut[1];  
    Height_cut[1]=130;  
    Height_cut[0]=85;
  }
  
  if(run==734)
  {
    Height_cut[2]=180;  
    Height_cut[1]=150;  
    Height_cut[0]=100;
  }

    Height->Draw();

    TLine * line;

    for(k=0; k<3; k++)
    {
    	line = new TLine(Height_cut[k],0,Height_cut[k],Height->GetBinContent(Height->GetMaximumBin())); //Coordenadas de la l韓ea
    	line->SetLineWidth(4.); //Tama駉 de la l韓ea
    	line->SetLineColor(kRed); //Color de la l韓ea
    	line->Draw("same"); //Para dibujar la l韓ea
     cout << Height_cut[k] << endl;
    }
    
    //Ajustes gaussianos 1D
    
    TF1 * gauss_Height_1 = new TF1("gauss_Height_1","gaus",0,500);
    TF1 * gauss_Height_2 = new TF1("gauss_Height_2","gaus",0,500);
    Height->Fit("gauss_Height_1","","",Height_cut[0],Height_cut[1]);
    Height->Fit("gauss_Height_2","","",Height_cut[1],Height_cut[2]);
    
    TF1* gaus_0 = new TF1("gaus 0","gaus",0,Height_max);
    TF1* gaus_1 = new TF1("gaus 1","gaus",0,Height_max);
    TF1* gaus_2 = new TF1("gaus 2","gaus",0,Height_max);
    TF1* gaus_3 = new TF1("gaus 3","gaus",0,Height_max);

    Height->GetXaxis()->SetRangeUser(0,Height_max);
    
    if(run!=741 && run!=735 && run!=734)
    {
      Height->Fit(gaus_0,"","",0,Height_cut[0]);
      Height->Fit(gaus_1,"","",Height_cut[0],Height_cut[1]);
      Height->Fit(gaus_2,"","",Height_cut[1],Height_cut[2]);
      Height->Fit(gaus_3,"","",Height_cut[2],Height_cut[2]+Height_cut[1]-Height_cut[0]);
    }
    
    if(run==741 || run==735 || run==734)
    {
      Height->Fit(gaus_0,"","",0,Height_cut[0]);
      Height->Fit(gaus_1,"","",Height_cut[0]/2.,Height_cut[0]);
      Height->Fit(gaus_2,"","",Height_cut[0],Height_cut[1]);
      Height->Fit(gaus_3,"","",Height_cut[1],Height_cut[2]);
    }

    TF1* gaus_tot = new TF1("gaus tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",0,Height_max);

    double par[12];
    double par_err[12];
    
    gaus_0->GetParameters(&par[0]);
    gaus_1->GetParameters(&par[3]);
    gaus_2->GetParameters(&par[6]);
    gaus_3->GetParameters(&par[9]);

    gaus_tot->SetParameters(par);
    
    if(run==321) gaus_tot->SetParLimits(4,80,100);
    //if(run==734) gaus_tot->SetParLimits(10,120,180);

    if(run!=741 && run!=735 && run!=734) Height->Fit(gaus_tot,"","",0,Height_cut[2]+Height_cut[1]-Height_cut[0]);
    if(run==741 || run==735 || run==734) Height->Fit(gaus_tot,"","",0,Height_cut[2]);

    gaus_tot->GetParameters(&par[0]);
    
    for(i=0; i<12; i++)
      par_err[i]=gaus_tot->GetParError(i);

    gauss_Height_1->SetParameters(&par[3]);
    gauss_Height_2->SetParameters(&par[6]);
    
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
    
    area_min = 0;
    area_max = 12000;
    
    //
    //Ajustes 2D
    //
    cout << area_min << " " << area_max << endl;
    
    TH2F* finger = new TH2F("finger",Form("Run %d [Channel %d]",run,ch),Height_nbins,0,Height_max,area_nbins,area_min,area_max);
    TF2 * finger_cal_0 = new TF2("finger_cal_0","bigaus",0,Height_cut[0],area_min,area_max);
    TF2 * finger_cal_1 = new TF2("finger_cal_1","bigaus",Height_cut[0],Height_cut[1],area_min,area_max);
    TF2 * finger_cal_2 = new TF2("finger_cal_2","bigaus",Height_cut[1],Height_cut[2],area_min,area_max);
    
    TF2 * finger_cal_tot = new TF2("finger_cal_tot","bigaus(0)+bigaus(6)+bigaus(12)",0,Height_cut[2],area_min,area_max);
    
    finger->GetXaxis()->SetRangeUser(0,700);
    finger->GetYaxis()->SetRangeUser(0,8000);
    
    //finger->GetXaxis()->SetRangeUser(0,300);
    //finger->GetYaxis()->SetRangeUser(0,3000);

    pulses->Draw(Form("areaFix%d:high%d>>finger",ch,ch),"","");

    finger_cal_0->SetParameter(0,gaus_0->GetParameter(0));
    finger_cal_0->SetParameter(1,gaus_0->GetParameter(1));
    finger_cal_0->SetParameter(2,gaus_0->GetParameter(2));
    finger_cal_0->SetParLimits(1,0+10,Height_cut[0]-10);
    finger_cal_0->SetParLimits(0,0,100000);
    
    finger_cal_1->SetParameter(0,gaus_1->GetParameter(0));
    finger_cal_1->SetParameter(1,gaus_1->GetParameter(1));
    finger_cal_1->SetParameter(2,gaus_1->GetParameter(2));
    finger_cal_1->SetParLimits(1,Height_cut[0]+10,Height_cut[1]-10);
    finger_cal_1->SetParLimits(0,0,100000);
    
    finger_cal_2->SetParameter(0,gaus_2->GetParameter(0));
    finger_cal_2->SetParameter(1,gaus_2->GetParameter(1));
    finger_cal_2->SetParameter(2,gaus_2->GetParameter(2));
    finger_cal_2->SetParLimits(1,Height_cut[1]+10,Height_cut[2]-10);
    finger_cal_2->SetParLimits(0,0,100000);
    
    finger->Fit("finger_cal_0","0","",0,Height_cut[0]);
    
    finger->Fit("finger_cal_1","0","",Height_cut[0],Height_cut[1]);
      
    finger->Fit("finger_cal_2","0","",Height_cut[1],Height_cut[2]);

    for(i=0; i<6; i++)
      finger_cal_tot->SetParameter(i,finger_cal_0->GetParameter(i));

    for(i=6; i<12; i++)
      finger_cal_tot->SetParameter(i,finger_cal_1->GetParameter(i-6));
    
    for(i=12; i<18; i++)
      finger_cal_tot->SetParameter(i,finger_cal_2->GetParameter(i-12));

    finger->Fit("finger_cal_tot","R && 0");
      
     TF2 * finger_cal_plot_0 = new TF2("finger_cal_plot_0","bigaus",0,500,0,8000);  
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

    /*for(k=0; k<3; k++)
    {
    	line = new TLine(Height_cut[k],0,Height_cut[k],3000); //Coordenadas de la l韓ea
    	line->SetLineWidth(4.); //Tama駉 de la l韓ea
      line->SetLineStyle(9);
    	line->SetLineColor(kBlack); //Color de la l韓ea
    	line->Draw("same"); //Para dibujar la l韓ea
    }*/

    TLegend *legend4 = new TLegend(0.48,0.15,0.89,0.30); //Tama駉 de la leyenda (x0,y0,x1,y1)
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
    
    area_cal->SetPoint(0,1.,finger_cal_tot->GetParameter(9));
    area_cal->SetPointError(0,0.,finger_cal_tot->GetParError(9));

    area_cal->SetPoint(1,2.,finger_cal_tot->GetParameter(15));
    area_cal->SetPointError(1,0.,finger_cal_tot->GetParError(15));

    area_cal->SetTitle("Area 2D calibration");

    area_cal->Fit(area_cal_fit);
    
     *cal=area_cal_fit->GetParameter(1);
     *cal_err=area_cal_fit->GetParError(1);
     
     *res=finger_cal_tot->GetParameter(10)/finger_cal_tot->GetParameter(9)*100;
     *res_err=100*sqrt(finger_cal_tot->GetParError(10)*finger_cal_tot->GetParError(10)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)+finger_cal_tot->GetParError(9)*finger_cal_tot->GetParError(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)*finger_cal_tot->GetParameter(10)*finger_cal_tot->GetParameter(10));
}

void areaFix_cal (TFile *frun, int run, int ch=0, int Height_nbins=300, double* phe_height, double* phe_sigma, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0)
{   
    TTree * pulses = (TTree*)frun->Get("ta");
    
    //
    //Definimos propiedades histograma de altura
    //

    TH1F *Height_aux = new TH1F("Height_aux","Height_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> Height_aux",ch),"","");

    int Height_max = Height_aux->GetXaxis()->GetXmax();
    int Height_min = Height_aux->GetXaxis()->GetXmin();

    double Height_bin_size = (Height_max-Height_min)/10000.;

    int Height_final_bin = Height_aux->FindLastBinAbove(1);
    int Height_initial_bin = Height_aux->FindFirstBinAbove(1);

    Height_max = Height_final_bin*Height_bin_size+Height_min;
    Height_min = Height_initial_bin*Height_bin_size+Height_min;

    TH1F *Height = new TH1F("Height", "Height",Height_nbins,Height_min,Height_max);

    pulses->Draw(Form("areaFix%d >> Height",ch),"","");
    
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
          gaus[i] = new TF1(Form("gaus%d",i+1),"gaus",Height_min,Height_max);
          gaus[i]->SetParameter(0,Height->GetBinContent(x_peaks[i]/double(Height_max-Height_min)*double(Height_nbins)));
          gaus[i]->SetParameter(1,x_peaks[i]);
          gaus[i]->SetParameter(2,150);
          Height->Fit(gaus[i],"0","",x_peaks[i]-100,x_peaks[i]+100);
        }
     
     TF1* gaus_tot;
     
     if(nfound==1) gaus_tot=new TF1("gaus_tot","gaus(0)",Height_min,Height_max);
     if(nfound==2) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)",Height_min,Height_max);
     if(nfound==3) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)",Height_min,Height_max);
     if(nfound==4) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",Height_min,Height_max);
     
     for(i=0; i<nfound; i++)
       for(j=0; j<3; j++)
         gaus_tot->SetParameter(i*3+j,gaus[i]->GetParameter(j));
     
     Height->Fit(gaus_tot,"0","",0,x_peaks[nfound-1]+30);
     
      for(i=0; i<nfound; i++)
      {
       for(j=0; j<3; j++)
         {gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
         gaus[i]->SetParError(j,gaus_tot->GetParError(j+i*3));}
      }
    
    *phe_height = gaus[1]->GetParameter(1);
    *phe_sigma = gaus[1]->GetParameter(2);
    
    //
    //Calibraci髇 phe
    //
    TGraphErrors* area_cal = new TGraphErrors();
    TF1 *area_cal_fit = new TF1("area_cal_fit","[0]+[1]*x");
    
    area_cal->SetPoint(0,1.,gaus[1]->GetParameter(1));
    area_cal->SetPointError(0,0.,gaus[1]->GetParError(1));

    area_cal->SetPoint(1,2.,gaus[2]->GetParameter(1));
    area_cal->SetPointError(1,0.,gaus[2]->GetParError(1));

    area_cal->SetTitle("Area 1D calibration");

    area_cal->Fit(area_cal_fit);
    
     *cal=area_cal_fit->GetParameter(1);
     *cal_err=area_cal_fit->GetParError(1);
     
     *res=gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)*100;
     *res_err=100*sqrt(gaus[1]->GetParError(2)*gaus[1]->GetParError(2)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)+gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)*gaus[1]->GetParameter(2)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)/gaus[1]->GetParameter(1)*gaus[1]->GetParError(1)*gaus[1]->GetParError(1));
    
    //
    //Plot
    //
    Height->SetLineColor(kBlack);
    
    //Height->GetXaxis()->SetRangeUser(0,700);
    
    Height->SetTitle(Form("Run %d [Channel %d]",run,ch));
    
    Height->GetXaxis()->SetTitle("Area (ADC)");
    
    Height->Draw("HIST");
    
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
      
    legend2->Draw("same");
}


void fwhm_distribution (TFile *frun, int run, int ch=0, int Height_nbins=300, double phe_height=0, double phe_sigma=0, TProfile* one_phe=0)
{
  TTree * pulses = (TTree*)frun->Get("ta");
  
  double Height, tFix, tWin;
  
  int bin1, bin2;
  
  TH1F *pulse = new TH1F();
  
  TH1F *fwhm_hist = new TH1F("FWHM",Form("Run %d [Channel %d]",run,ch),50,0,200);
  
  TF1 * gaus_fit = new TF1("gaus","gaus",0,200);
  
  pulses->SetBranchAddress(Form("high%d",ch),&Height);
  pulses->SetBranchAddress(Form("pulse%d",ch),&pulse);
  pulses->SetBranchAddress(Form("tFix_mod%d",ch),&tFix);
  pulses->SetBranchAddress(Form("tWin_mod%d",ch),&tWin);
  
  for(int i=0; i<pulses->GetEntries(); i++)
  {
    pulses->GetEntry(i);
    
    if(Height<=(phe_height+phe_sigma/2.) && Height>=(phe_height-phe_sigma/2.))
      {
        pulse->GetXaxis()->SetRangeUser(tFix,tFix+tWin);
        for(int j=pulse->GetMaximumBin(); j>=1; j--) if(pulse->GetBinContent(j)<pulse->GetMaximum()/2.){bin1=j; break;}
        for(int j=pulse->GetMaximumBin(); j<=tFix+tWin+10; j++) if(pulse->GetBinContent(j)<pulse->GetMaximum()/2.){bin2=j; break;}
        double fwhm = (bin2 - bin1)*4;
        
        for(int k=1;k<=250;k++) one_phe->Fill(k*4-2,pulse->GetBinContent(k+(pulse->GetMaximumBin()-100)));
        
        fwhm_hist->Fill(fwhm);
      }
  }
  
  fwhm_hist->GetXaxis()->SetTitle("FWHM (ns)");
  
  fwhm_hist->SetLineColor(kBlack);
  
  fwhm_hist->Draw();
  
  fwhm_hist->Fit("gaus","0","",0,fwhm_hist->GetMaximumBin()*4+20);
  
  gaus_fit->SetLineWidth(4);
  gaus_fit->Draw("same");
  
  TPaveText *pv = new TPaveText(0.50,0.75,0.83,0.85,"NDC"); //Coordenadas del texto (x0,y0,x1,y1)
    pv->SetBorderSize(0); //Tama駉 del borde
    pv->SetFillColor(0); //Color del fondo
    pv->SetTextSize(.05); //Tama駉 del texto
    pv->SetTextColor(kBlack); //Color del texto
    pv->SetMargin(0.01); //Margen del texto
    pv->AddText(Form("FWHM = %.0f ns",gaus_fit->GetParameter(1)));
    pv->AddText(Form("#sigma FWHM = %.0f ns",gaus_fit->GetParameter(2)));
    pv->SetTextFont(102); //Fuente del texto
    pv->SetTextAlign(22); //Alineaci髇 del texto 22=centrada en horizontal y vertical  
    pv->Draw(); //Para dibujar el cuadro de texto     
  
}
