// Example Program for converting MIDAS format to ROOT format.
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMath.h>
#include <vector>
#include <TH1F.h>
#include <TSpectrum.h>
#include <TStyle.h> //añadido

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h" //añadido

/*
#include "TRootanaEventLoop.hxx"
#include "TDataContainer.hxx"
#include "TDartEvent.hxx"
#include "TV1730RawData.hxx"

#include "TDartAnaManager.hxx"
#include "TEventProcessor.hxx"
#include "TEventMusic.hxx"

#include <TMidasEvent.h>
#include <TDataContainer.hxx>
#include <THistogramArrayBase.h>

#include "midasio.h"
*/

const int VMEBUS_BOARDNO = 0;

void phe_calibration (TFile *frun, int run, int ch=0, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0, int Height_nbins=300, int area_nbins=350);
void phe_calibration_2 (TFile *frun, int run, int ch=0, int Area_nbins=300, double* cal=0, double* cal_err=0, double* res=0, double* res_err=0);

int main(int argc, char *argv[])
{ 
  //
  //Plots configuration
  //
  gStyle->SetImageScaling(3.);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleFont(102,"XYZ");
  gStyle->SetTitleFont(102,"");
  gStyle->SetTitleSize(0.07,"");
  gStyle->SetLegendFont(102);
  gStyle->SetLabelFont(102,"XYZ");
  gStyle->SetTextFont(102);
  gStyle->SetOptStat(0);
    
  //
  //Analysis configuration
  //
  
  //We load the input file (Level 1)
  const char * file = argv[1];
  
  std::string input_file_name = Form("../analyzed/ANAISplusLevel1/%s_analyzed_1.root",file);
  
  TFile * input_file = new TFile(input_file_name.c_str(), "read");
  
  TTree * Level1Information = (TTree*)input_file->Get("ta");
  
  //Configuration file
  int nEvents;
  
  int nSamples;
  int nChannels;
  bool cal;
  double amp;
  int nPtsBsl;
  bool metBsl;
  bool metInt;
  int nProm;
  double thProm;
  int sIni;
  int sWin;
  double threshold;
  
  std::string nameConfig = std::string(argv[2]);
  
  std::string fileConfig = Form("%s.txt",nameConfig.c_str());
  
  std::ifstream fConfig (fileConfig.c_str(), std::ios::in);
  
  //Read configuration
	fConfig >> nSamples;
	fConfig >> nChannels;
  fConfig >> cal;
	fConfig >> amp;
  fConfig >> nPtsBsl;
	fConfig >> metBsl;
  fConfig >> metInt;
	fConfig >> nProm;
  fConfig >> thProm;
	fConfig >> sIni;
  fConfig >> sWin;
	fConfig >> threshold;
 
  nEvents=Level1Information->GetEntries();
 
  //Configuration check
  if(amp==0)
  {
    printf("Amplification must be > 0 \n");
    return -1;
  }
  
  if(nPtsBsl<=10 || nPtsBsl>nSamples)
  {
    printf("Baseline points must be > 10 and < nSamples \n");
    return -1;
  }
  
  if(nProm>nEvents)
  {
    printf("Promedium events higher than recorded events \n");
    return -1;
  }
  
  if(thProm<=0 || thProm>100)
  {
    printf("Promedium threshold must be between 0 and 100 %% \n");
    return -1;
  }
  
  if(sIni<0 || sIni>nSamples-sWin)
  {
    printf("sIni must be > 0 and < nSamples-sWin \n");
    return -1;
  }
  
  if(sWin<0 || sWin>nSamples)
  {
    printf("sWin must be > 0 and < nSamples \n");
    return -1;
  }
  
  if(threshold < 0)
  {
    printf("Threshold must be > 0 \n");
    return -1;
  }
  
  //Parameters used in this level
  double Vov=std::atof(argv[3]);
   
  int run_cal_i=std::stoi(argv[4]); 
  
  int run_cal_f=std::stoi(argv[5]); 
  
  //Input variables
  double bsl[nChannels];
  double rms[nChannels];
  double t0[nChannels];
  double area[nChannels];
  double areaFix[nChannels];
  double tmax[nChannels];
  double high[nChannels];
  double minBsl[nChannels];
  
  for(int j = 0; j < nChannels; j++) 
  {
    Level1Information->SetBranchAddress(Form("bsl%d",j),&bsl[j]);
    Level1Information->SetBranchAddress(Form("rms%d",j),&rms[j]);
    Level1Information->SetBranchAddress(Form("t0%d",j),&t0[j]);
    Level1Information->SetBranchAddress(Form("area%d",j),&area[j]);
    Level1Information->SetBranchAddress(Form("areaFix%d",j),&areaFix[j]);
    Level1Information->SetBranchAddress(Form("tmax%d",j),&tmax[j]);
    Level1Information->SetBranchAddress(Form("high%d",j),&high[j]);
    Level1Information->SetBranchAddress(Form("minBsl%d",j),&minBsl[j]);
  }
 
  //
  //
  //
  
  //
  //Photoelectron calibration
  //
  if(cal==true)
  {
    int ch=0;
   
    double cal_value, cal_err, res, res_err;
    
    TCanvas *c1 = new TCanvas("c1","The Graph example",10,10,1200,700);
    
    if(nChannels<=4 && nChannels>1) {c1 = new TCanvas("c1","The Graph example",10,10,1200,355); c1->Divide(nChannels,1);}
    else {c1->Divide(int(nChannels/2.)+nChannels%2,2);}
    
    FILE *f;
    f=fopen(Form("../analyzed/ANAISplusLevel2/Calibration/Channel_calibration_%d.txt",run_cal_i),"w");
    
    FILE *f_res;
    f_res=fopen(Form("../analyzed/ANAISplusLevel2/Calibration/Channel_resolution_%d.txt",run_cal_i),"w");
    
    for(ch=0; ch<nChannels; ch++)
    { 
      c1->cd(ch+1)->SetGrid();
      
      //phe_calibration(input_file, run_cal_i, ch, &cal_value, &cal_err, &res, &res_err, std::stoi(argv[7]), std::stoi(argv[6]));
      phe_calibration_2 (input_file, run_cal_i, ch, std::stoi(argv[6]), &cal_value, &cal_err, &res, &res_err);
      
      fprintf(f,"%d %f %f %f\n", ch, cal_value, cal_err, Vov);
      fprintf(f_res,"%d %f %f %f\n", ch, res, res_err, Vov);
    }  
    
    c1->SaveAs(Form("../analyzed/ANAISplusLevel2/Calibration/Plots/Channel_calibration_%d.png",run_cal_i));
      
    fclose(f);
    fclose(f_res);
  }
  
  //
  //
  //
  
  //
  //Apply the calibration to extract number of phes
  //
  
  //Output preparation
  
  std::string filenew = Form("../analyzed/ANAISplusLevel2/%s_analyzed_2.root",file);
  
  TFile * fnew = new TFile(filenew.c_str(), "recreate");
  
  TTree * t2 = new TTree("ta","");
  
  //Output variables
  double n[nChannels];
  double nFix[nChannels];
  double nTot=0;
  double nFixTot=0;
  
  for(int j = 0; j < nChannels; j++) 
  {
    t2->Branch(Form("n%d",j),&n[j],Form("n%d/D",j));
    t2->Branch(Form("nFix%d",j),&nFix[j],Form("nFix%d/D",j));
  }
  
  t2->Branch("nTot",&nTot,"nTot/D");
  t2->Branch("nFixTot",&nFixTot,"nFixTot/D");
  
  //We load the calibration
  
  double channel_calibration[nChannels];
  double channel_calibration_err[nChannels];
  double channel_cal_origin[nChannels];
  double channel_cal_origin_err[nChannels];
  
  int channel;
  
  if(run_cal_f==0) //Same Vov in calibration
  {
    std::ifstream calibration_file (Form("../analyzed/ANAISplusLevel2/Calibration/Channel_calibration_%d.txt",run_cal_i), std::ios::in);
    
    for(int j=0; j<nChannels; j++)
    {
      calibration_file >> channel >> channel_calibration[j] >> channel_calibration_err[j] >> Vov;
    }
  }
  
  else //Extrapolation from Vov
  {
    std::ifstream calibration_file (Form("../analyzed/ANAISplusLevel2/Calibration/Channel_calibration_%d_%d.txt",run_cal_i,run_cal_f), std::ios::in);
    
    for(int j=0; j<nChannels; j++)
    {
      calibration_file >> channel >> channel_calibration[j] >> channel_calibration_err[j] >> channel_cal_origin[j] >> channel_cal_origin_err[j];
      
      channel_calibration[j]=channel_calibration[j]*Vov/2.+channel_cal_origin[j];
    }
  }
  
  //We compute the number of phe
  int ev=0;
  
  for(int i=0; i<nEvents; i++)
    { 
      Level1Information->GetEntry(i);
      
      for(int j=0; j<nChannels; j++)
      {
        n[j]=area[j]/channel_calibration[j];
        nFix[j]=areaFix[j]/channel_calibration[j];
        
        nTot+=n[j];
        nFixTot+=nFix[j];
      }
      
      t2->Fill();
      
      nTot=0;
      nFixTot=0;
      
      if(ev%1000 == 0) {t2->AutoSave(); std::cout << ev << " " << file<< std::endl;}
      
      ev++;
    }
  
  t2->Write();

  fnew->Close();
  
  return 0;
}

void phe_calibration (TFile *frun, int run, int ch, double* cal, double* cal_err, double* res, double* res_err, int Height_nbins, int area_nbins)
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
    		{
    		 Height_cut[k]=((i-1)*5+2.5)*(double(Height_max)/double(Height_nbins)); //Se pasa de bin a high0
    		 k++;
    		}
    }

    if(Height_cut[2]-Height_cut[1]>1.2*(Height_cut[1]-Height_cut[0]) || Height_cut[2]-Height_cut[1]<0.6*(Height_cut[1]-Height_cut[0])) //Ajuste por si se va mucho el limite superior de 2 phe
	      Height_cut[2]=Height_cut[1]+Height_cut[1]-Height_cut[0];

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
    //Histograma de área
    //
    TH1F *area_aux = new TH1F("area_aux","area_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> area_aux",ch),Form("high%d>0",ch),"");

    gStyle->SetOptStat(0);

    int area_max = area_aux->GetXaxis()->GetXmax();
    int area_min = area_aux->GetXaxis()->GetXmin();

    double area_bin_size = double(area_max-area_min)/area_nbins;

    int area_final_bin = area_aux->FindLastBinAbove(1);
    int area_initial_bin = area_aux->FindFirstBinAbove(1);

    area_max = area_final_bin*area_bin_size+area_min;
    area_min = area_initial_bin*area_bin_size+area_min;
    
    area_min = 0;
    area_max = 12000;
    
    //
    //Ajustes 2D
    //
    TH2F* finger = new TH2F("finger",Form("Run %d [Channel %d]",run,ch),Height_nbins,0,Height_max,area_nbins,area_min,area_max);
    TF2 * finger_cal_0 = new TF2("finger_cal_0","bigaus",0,Height_cut[0],area_min,area_max);
    TF2 * finger_cal_1 = new TF2("finger_cal_1","bigaus",Height_cut[0],Height_cut[1],area_min,area_max);
    TF2 * finger_cal_2 = new TF2("finger_cal_2","bigaus",Height_cut[1],Height_cut[2],area_min,area_max);
    
    TF2 * finger_cal_tot = new TF2("finger_cal_tot","bigaus(0)+bigaus(6)+bigaus(12)",0,Height_cut[2],area_min,area_max);
    
    finger->GetXaxis()->SetRangeUser(0,700);
    finger->GetYaxis()->SetRangeUser(0,8000);

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


    TLegend *legend4 = new TLegend(0.48,0.15,0.89,0.30); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend4->SetTextSize(0.04); //Tamaño del texto
    legend4->SetFillColor(kWhite); //Fondo de la leyenda

    legend4->AddEntry(finger_cal_plot_1,Form("1 phe"),"l");
    legend4->AddEntry(finger_cal_plot_2,Form("2 phe"),"l");
    legend4->AddEntry(finger_cal_tot,Form("#chi^{2}/ndf=%0.2f",finger_cal_tot->GetChisquare()/finger_cal_tot->GetNDF()),"l");

    legend4->Draw("same");
    
    //
    //Calibración phe
    //
    TGraphErrors* area_cal = new TGraphErrors();
    TF1 *area_cal_fit = new TF1("area_cal_fit","[0]*x");
    
    area_cal->SetPoint(0,1.,finger_cal_tot->GetParameter(9));
    area_cal->SetPointError(0,0.,finger_cal_tot->GetParError(9));

    area_cal->SetPoint(1,2.,finger_cal_tot->GetParameter(15));
    area_cal->SetPointError(1,0.,finger_cal_tot->GetParError(15));

    area_cal->SetTitle("Area 2D calibration");

    area_cal->Fit(area_cal_fit);
    
     *cal=area_cal_fit->GetParameter(0);
     *cal_err=area_cal_fit->GetParError(0);
     
     *res=finger_cal_tot->GetParameter(10)/finger_cal_tot->GetParameter(9)*100;
     *res_err=100*sqrt(finger_cal_tot->GetParError(10)*finger_cal_tot->GetParError(10)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)+finger_cal_tot->GetParError(9)*finger_cal_tot->GetParError(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)/finger_cal_tot->GetParameter(9)*finger_cal_tot->GetParameter(10)*finger_cal_tot->GetParameter(10));
}

void phe_calibration_2 (TFile *frun, int run, int ch, int Area_nbins, double* cal, double* cal_err, double* res, double* res_err)
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

    int Area_final_bin = Area_aux->FindLastBinAbove(10);
    int Area_initial_bin = Area_aux->FindFirstBinAbove(10);

    Area_max = Area_final_bin*Area_bin_size+Area_min;
    Area_min = Area_initial_bin*Area_bin_size+Area_min;

    TH1F *Area = new TH1F("Area", "Area",Area_nbins,Area_min,Area_max);

    pulses->Draw(Form("areaFix%d >> Area",ch),"","");
    
    //
    //Se buscan los picos
    //
    TSpectrum *s = new TSpectrum(1000);
    int nfound;
    nfound = s->Search(Area,1.7,"",0.005);
    double* x_peaks = s->GetPositionX();
    
    std::sort(x_peaks, x_peaks+nfound);
    
    if(nfound>5) nfound=5;
    
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
          if(i==0) Area->Fit(gaus[i],"0","",x_peaks[i]-1000,x_peaks[i]+1000);
          else Area->Fit(gaus[i],"0","",x_peaks[i]-sqrt(x_peaks[i])*15,x_peaks[i]+sqrt(x_peaks[i])*15);
        }
     
     TF1* gaus_tot;
     
     if(nfound==1) gaus_tot=new TF1("gaus_tot","gaus(0)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*15);
     if(nfound==2) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*15);
     if(nfound==3) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*15);
     if(nfound==4) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*15);
     if(nfound==5) gaus_tot=new TF1("gaus_tot","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*15);
     
     for(i=0; i<nfound; i++)
       for(j=0; j<3; j++)
       {
         gaus_tot->SetParameter(i*3+j,gaus[i]->GetParameter(j));
         if(j==1) gaus_tot->SetParLimits(i*3+j,gaus[i]->GetParameter(j)-gaus[i]->GetParameter(j+1),gaus[i]->GetParameter(j)+gaus[i]->GetParameter(j+1));
       }
     
     Area->Fit(gaus_tot,"0","",Area_min,x_peaks[nfound-1]+sqrt(x_peaks[nfound-1])*10);
     
      for(i=0; i<nfound; i++)
      {
       for(j=0; j<3; j++)
         {gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
         gaus[i]->SetParError(j,gaus_tot->GetParError(j+i*3));}
      }
    
    //
    //Calibración phe
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
    int plot_control=0;
    
    if(plot_control==0)
    {
    Area->SetLineColor(kBlack);
    
    Area->SetTitle(Form("Run %d [Channel %d]",run,ch));
    
    Area->GetXaxis()->SetTitle("Area (ADC)");
    
    Area->Draw("HIST");
    
    TLegend *legend2 = new TLegend(0.70,0.87-nfound*0.05,0.87,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend2->SetTextSize(0.04); //Tamaño del texto
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
       
       TLegend *legend2 = new TLegend(0.13,0.80,0.70,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
       legend2->SetTextSize(0.04); //Tamaño del texto
       legend2->SetFillColor(kWhite); //Fondo de la leyenda
       legend2->AddEntry(area_cal_fit,Form("Cal = %.1f #pm %.1f",area_cal_fit->GetParameter(1),area_cal_fit->GetParError(1)),"l");
       
       legend2->Draw();
    }
}
/*{   
    TTree * pulses = (TTree*)frun->Get("ta");
    
    //
    //Definimos propiedades histograma de altura
    //

    TH1F *Height_aux = new TH1F("Height_aux","Height_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> Height_aux",ch));

    int Height_max = Height_aux->GetXaxis()->GetXmax();
    int Height_min = Height_aux->GetXaxis()->GetXmin();

    double Height_bin_size = (Height_max-Height_min)/10000.;

    int Height_final_bin = Height_aux->FindLastBinAbove(1);
    int Height_initial_bin = Height_aux->FindFirstBinAbove(1);

    Height_max = Height_final_bin*Height_bin_size+Height_min;
    Height_min = Height_initial_bin*Height_bin_size+Height_min;

    TH1F *Height = new TH1F("Height", "Height",Height_nbins,Height_min,Height_max);

    pulses->Draw(Form("areaFix%d >> Height",ch));
    
    //
    //Se buscan los picos
    //
    TSpectrum *s = new TSpectrum(1000,1.0);
    int nfound;
    nfound = s->Search(Height,1.5,"",0.01);
    double* x_peaks = s->GetPositionX();
    
    std::sort(x_peaks, x_peaks+nfound);
    
    if(nfound>4) nfound=4;
    
    //if(ch==0) nfound=3;
    
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
         {
             gaus[i]->SetParameter(j,gaus_tot->GetParameter(j+i*3));
             gaus[i]->SetParError(j,gaus_tot->GetParError(j+i*3));
         }
      }
    
    //
    //Calibración phe
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
    
    TLegend *legend2 = new TLegend(0.70,0.87-nfound*0.05,0.87,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend2->SetTextSize(0.04); //Tamaño del texto
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
}*/
