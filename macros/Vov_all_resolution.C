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
#include <TGraphErrors.h>

void Vov_all_resolution (int run_i=532, int run_f=545, const char* dir="LC_cold_2", double Vov_i=4, double delta_Vov=0.5, double Vov_err=0.55)
{
  
  gStyle->SetImageScaling(3.);
  gStyle->SetHistLineWidth(3);
  gStyle->SetTitleFont(102,"XYZ");
  gStyle->SetTitleFont(102,"");
  gStyle->SetTitleSize(0.07,"");
  gStyle->SetLegendFont(102);
  gStyle->SetLabelFont(102,"XYZ");
  gStyle->SetTextFont(102);
  
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit","Simplex");
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
  
  //
  //Graph cal
  // 
  TGraphErrors* cal_ch0 = new TGraphErrors();
  TGraphErrors* cal_ch1 = new TGraphErrors(); 
  TGraphErrors* cal_ch2 = new TGraphErrors();
  TGraphErrors* cal_ch3 = new TGraphErrors();
  TGraphErrors* cal_ch4 = new TGraphErrors();
  TGraphErrors* cal_ch5 = new TGraphErrors();
  TGraphErrors* cal_ch6 = new TGraphErrors();
  TGraphErrors* cal_ch7 = new TGraphErrors(); 
  
    cal_ch0->SetMarkerStyle(8);
    cal_ch1->SetMarkerStyle(8);
    cal_ch2->SetMarkerStyle(8);
    cal_ch3->SetMarkerStyle(8);
    cal_ch4->SetMarkerStyle(8);
    cal_ch5->SetMarkerStyle(8);
    cal_ch6->SetMarkerStyle(8);
    cal_ch7->SetMarkerStyle(8);
    
    cal_ch0->SetMarkerColor(kGray);
    cal_ch1->SetMarkerColor(kRed-4);
    cal_ch2->SetMarkerColor(kGreen-4);
    cal_ch3->SetMarkerColor(kCyan-3);
    cal_ch4->SetMarkerColor(kMagenta-7);
    cal_ch5->SetMarkerColor(kOrange-3);
    cal_ch6->SetMarkerColor(kBlue-3);
    cal_ch7->SetMarkerColor(kViolet-5);
    
    cal_ch0->SetLineColor(kGray);
    cal_ch1->SetLineColor(kRed-4);
    cal_ch2->SetLineColor(kGreen-4);
    cal_ch3->SetLineColor(kCyan-3);
    cal_ch4->SetLineColor(kMagenta-7);
    cal_ch5->SetLineColor(kOrange-3);
    cal_ch6->SetLineColor(kBlue-3);
    cal_ch7->SetLineColor(kViolet-5);
    
    int ch=0;
  
  double cal[8];
  double cal_err[8];
  
  for(ch=0;ch<8;ch++)
  {
    cal[ch]=186;
    cal_err[ch]=0;
  }
  //
  //
  //
  
  //
  //Llenamos el graph de cal
  //
  double energy, Vov=Vov_i;
  int run=run_i, j, k=0;
  
  double area, area_err, sigma, sigma_err;
  double cal_val, cal_val_err;
  
  for(run=run_i; run<=run_f; run++)
  {
    std::ifstream inputfile(Form("./Plots/Calibration_analysis/Channel_resolution_%d.txt",run), std::ios::in );
    
    for(j=0;j<=8;j++)
    {
      inputfile >> ch >> area >> area_err;
      
      if(ch==0) {cal_ch0->SetPoint(k, Vov, area); cal_ch0->SetPointError(k, Vov_err, area_err);}
      if(ch==1) {cal_ch1->SetPoint(k, Vov, area); cal_ch1->SetPointError(k, Vov_err, area_err);}
      if(ch==2) {cal_ch2->SetPoint(k, Vov, area); cal_ch2->SetPointError(k, Vov_err, area_err);}
      if(ch==3) {cal_ch3->SetPoint(k, Vov, area); cal_ch3->SetPointError(k, Vov_err, area_err);}
      if(ch==4) {cal_ch4->SetPoint(k, Vov, area); cal_ch4->SetPointError(k, Vov_err, area_err);}
      if(ch==5) {cal_ch5->SetPoint(k, Vov, area); cal_ch5->SetPointError(k, Vov_err, area_err);}
      if(ch==6) {cal_ch6->SetPoint(k, Vov, area); cal_ch6->SetPointError(k, Vov_err, area_err);}
      if(ch==7) {cal_ch7->SetPoint(k, Vov, area); cal_ch7->SetPointError(k, Vov_err, area_err);}
    }
    
    Vov+=delta_Vov;
    k++;
  }
  
  //
  //Modelamos la resolucion
  //
  //TF1*Fano = new TF1("F","[0]*(1-[1]*x*(1-exp(-x/[2])))^[3]",0,11);
  //TF1*Fano = new TF1("F","[0]*(1+2*[1]*x*(1-exp(-x/[2]))/(1-[1]*x*(1-exp(-x/[2]))))",0,11);
  
  TF1* cal0 = new TF1("cal0","([0]/x)",0,11);
  TF1* cal1 = new TF1("cal1","([0]/x)",0,11);
  TF1* cal2 = new TF1("cal2","([0]/x)",0,11);
  TF1* cal3 = new TF1("cal3","([0]/x)",0,11);
  TF1* cal4 = new TF1("cal4","([0]/x)",0,11);
  TF1* cal5 = new TF1("cal5","([0]/x)",0,11);
  TF1* cal6 = new TF1("cal6","([0]/x)",0,11);
  TF1* cal7 = new TF1("cal7","([0]/x)",0,11);
  
  FILE *f;
  f=fopen(Form("./Plots/%s/Channel_Vov_res_%d_%d.txt",dir,run_i,run_f),"w");
  
  cal_ch0->Fit("cal0","0","",0,6.5);
  cal_ch1->Fit("cal1","0","",0,6.5);
  cal_ch2->Fit("cal2","0","",0,6.5);
  cal_ch3->Fit("cal3","0","",0,6.5);
  cal_ch4->Fit("cal4","0","",0,6.5);
  cal_ch5->Fit("cal5","0","",0,6.5);
  cal_ch6->Fit("cal6","0","",0,6.5);
  cal_ch7->Fit("cal7","0","",0,6.5);
  
    cal0->SetLineColor(kGray);
    cal1->SetLineColor(kRed-4);
    cal2->SetLineColor(kGreen-4);
    cal3->SetLineColor(kCyan-3);
    cal4->SetLineColor(kMagenta-7);
    cal5->SetLineColor(kOrange-3);
    cal6->SetLineColor(kBlue-3);
    cal7->SetLineColor(kViolet-4);
    
    cal0->SetLineStyle(9);
    cal1->SetLineStyle(9);
    cal2->SetLineStyle(9);
    cal3->SetLineStyle(9);
    cal4->SetLineStyle(9);
    cal5->SetLineStyle(9);
    cal6->SetLineStyle(9);
    cal7->SetLineStyle(9);
  
  fprintf(f,"%d %f %f \n",0,cal0->GetParameter(0),cal0->GetParError(0));
  fprintf(f,"%d %f %f \n",1,cal1->GetParameter(0),cal1->GetParError(0));
  fprintf(f,"%d %f %f \n",2,cal2->GetParameter(0),cal2->GetParError(0));
  fprintf(f,"%d %f %f \n",3,cal3->GetParameter(0),cal3->GetParError(0));
  fprintf(f,"%d %f %f \n",4,cal4->GetParameter(0),cal4->GetParError(0));
  fprintf(f,"%d %f %f \n",5,cal5->GetParameter(0),cal5->GetParError(0));
  fprintf(f,"%d %f %f \n",6,cal6->GetParameter(0),cal6->GetParError(0));
  fprintf(f,"%d %f %f \n",7,cal7->GetParameter(0),cal7->GetParError(0));
  
  fclose(f);
  
  TCanvas *c1 = new TCanvas("c1","The Graph example",200,10,700,500);
  cal0->SetTitle("Resolution");
  cal0->GetYaxis()->SetTitle("Phe resolution (%)");
  cal0->GetXaxis()->SetTitle("Vov (V)");
  
  cal0->GetYaxis()->SetRangeUser(30,70);
  
  cal0->Draw("");
  cal1->Draw("same");
  cal2->Draw("same");
  cal3->Draw("same");
  cal4->Draw("same");
  cal5->Draw("same");
  cal6->Draw("same");
  cal7->Draw("same");

  cal_ch0->Draw("P");
  cal_ch1->Draw("P");
  cal_ch2->Draw("P");
  cal_ch3->Draw("P");
  cal_ch4->Draw("P");
  cal_ch5->Draw("P");
  cal_ch6->Draw("P");
  cal_ch7->Draw("P");
  
  TLegend *legend = new TLegend(0.13,0.42,0.23,0.87); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend->SetTextSize(0.04); //Tamaño del texto
    legend->SetFillColor(kWhite); //Fondo de la leyenda

    legend->AddEntry(cal_ch0,"Ch 0","p");
    legend->AddEntry(cal_ch1,"Ch 1","p");
    legend->AddEntry(cal_ch2,"Ch 2","p");
    legend->AddEntry(cal_ch3,"Ch 3","p");
    legend->AddEntry(cal_ch4,"Ch 4","p");
    legend->AddEntry(cal_ch5,"Ch 5","p");
    legend->AddEntry(cal_ch6,"Ch 6","p");
    legend->AddEntry(cal_ch7,"Ch 7","p");
  
  legend->Draw();
  
  c1->SaveAs(Form("./Plots/%s/Vov_resolution_%d_%d.png",dir,run_i,run_f));
  
  //Plot separados
  TCanvas *c2 = new TCanvas("c2","The Graph example",10,10,1200,700);
  c2->Divide(4,2,0,0.01);
  
  cal0->SetTitle("Channel 0");
  cal0->GetYaxis()->SetTitle("Phe resolution (%)");
  cal0->GetXaxis()->SetTitle("Vov (V)");
  cal1->SetTitle("Channel 1");
  cal1->GetYaxis()->SetTitle("Phe resolution (%)");
  cal1->GetXaxis()->SetTitle("Vov (V)");
  cal2->SetTitle("Channel 2");
  cal2->GetYaxis()->SetTitle("Phe resolution (%)");
  cal2->GetXaxis()->SetTitle("Vov (V)");
  cal3->SetTitle("Channel 3");
  cal3->GetYaxis()->SetTitle("Phe resolution (%)");
  cal3->GetXaxis()->SetTitle("Vov (V)");
  cal4->SetTitle("Channel 4");
  cal4->GetYaxis()->SetTitle("Phe resolution (%)");
  cal4->GetXaxis()->SetTitle("Vov (V)");
  cal5->SetTitle("Channel 5");
  cal5->GetYaxis()->SetTitle("Phe resolution (%)");
  cal5->GetXaxis()->SetTitle("Vov (V)");
  cal6->SetTitle("Channel 6");
  cal6->GetYaxis()->SetTitle("Phe resolution (%)");
  cal6->GetXaxis()->SetTitle("Vov (V)");
  cal7->SetTitle("Channel 7");
  cal7->GetYaxis()->SetTitle("Phe resolution (%)");
  cal7->GetXaxis()->SetTitle("Vov (V)");
  
  cal_ch0->SetTitle("Channel 0");
  cal_ch0->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch0->GetXaxis()->SetTitle("Vov (V)");
  cal_ch1->SetTitle("Channel 1");
  cal_ch1->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch1->GetXaxis()->SetTitle("Vov (V)");
  cal_ch2->SetTitle("Channel 2");
  cal_ch2->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch2->GetXaxis()->SetTitle("Vov (V)");
  cal_ch3->SetTitle("Channel 3");
  cal_ch3->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch3->GetXaxis()->SetTitle("Vov (V)");
  cal_ch4->SetTitle("Channel 4");
  cal_ch4->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch4->GetXaxis()->SetTitle("Vov (V)");
  cal_ch5->SetTitle("Channel 5");
  cal_ch5->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch5->GetXaxis()->SetTitle("Vov (V)");
  cal_ch6->SetTitle("Channel 6");
  cal_ch6->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch6->GetXaxis()->SetTitle("Vov (V)");
  cal_ch7->SetTitle("Channel 7");
  cal_ch7->GetYaxis()->SetTitle("Phe resolution (%)");
  cal_ch7->GetXaxis()->SetTitle("Vov (V)");
  
  cal0->SetLineStyle(9);
    cal1->SetLineStyle(9);
    cal2->SetLineStyle(9);
    cal3->SetLineStyle(9);
    cal4->SetLineStyle(9);
    cal5->SetLineStyle(9);
    cal6->SetLineStyle(9);
    cal7->SetLineStyle(9);
    
  cal_ch0->GetYaxis()->SetRangeUser(30,70);
  cal_ch1->GetYaxis()->SetRangeUser(30,70);
  cal_ch2->GetYaxis()->SetRangeUser(30,70);
  cal_ch3->GetYaxis()->SetRangeUser(30,70);
  cal_ch4->GetYaxis()->SetRangeUser(30,70);
  cal_ch5->GetYaxis()->SetRangeUser(30,70);
  cal_ch6->GetYaxis()->SetRangeUser(30,70);
  cal_ch7->GetYaxis()->SetRangeUser(30,70);
  
  TPaveText *info = new TPaveText(0.57,0.05,0.95,0.15,"brNDC"); //Tamaño de la leyenda (x0,y0,x1,y1)
  info->SetTextSize(0.04); //Tamaño del texto
  info->SetFillColor(kWhite); //Fondo de la leyenda
  info->SetBorderSize(0);
  
  TPaveText *info2 = new TPaveText(0.57,0.1,0.95,0.2,"brNDC"); //Tamaño de la leyenda (x0,y0,x1,y1)
  info2->SetTextSize(0.04); //Tamaño del texto
  info2->SetFillColor(kWhite); //Fondo de la leyenda
  info2->SetBorderSize(0);
  
  c2->SetLeftMargin(0.35);
  
  for(ch=0; ch<=7; ch++)
    { 
      c2->cd(ch+1);
      
      info = new TPaveText(0.17,0.80,0.4,0.88,"brNDC");
      info->SetTextSize(0.05); //Tamaño del texto
      info->SetFillColor(kWhite); //Fondo de la leyenda
      info->SetBorderSize(0);
  
      info2 = new TPaveText(0.17,0.80,0.4,0.88,"brNDC");
      info2->SetTextSize(0.05); //Tamaño del texto
      info2->SetFillColor(kWhite); //Fondo de la leyenda
      info2->SetBorderSize(0);
      
      if(ch==0)
      {
        info = new TPaveText(0.27,0.80,0.5,0.88,"brNDC");
        info->SetTextSize(0.05); //Tamaño del texto
        info->SetFillColor(kWhite); //Fondo de la leyenda
        info->SetBorderSize(0);
      }
      
      if(ch==4)
      {
        info2 = new TPaveText(0.27,0.80,0.5,0.88,"brNDC");
        info2->SetTextSize(0.05); //Tamaño del texto
        info2->SetFillColor(kWhite); //Fondo de la leyenda
        info2->SetBorderSize(0);
      }
  
      
      if(ch==0 || ch==4) {gPad->SetLeftMargin(0.15);}
      
      if(ch==0) {cal_ch0->Draw("AP"); cal0->Draw("same");}
      if(ch==1) {cal_ch1->Draw("AP"); cal1->Draw("same");}
      if(ch==2) {cal_ch2->Draw("AP"); cal2->Draw("same");}
      if(ch==3) {cal_ch3->Draw("AP"); cal3->Draw("same");}
      if(ch==4) {cal_ch4->Draw("AP"); cal4->Draw("same");}
      if(ch==5) {cal_ch5->Draw("AP"); cal5->Draw("same");}
      if(ch==6) {cal_ch6->Draw("AP"); cal6->Draw("same");}
      if(ch==7) {cal_ch7->Draw("AP"); cal7->Draw("same");}
      
      if(ch==0) {info->AddText(Form("#chi^{2}/ndf = %.2f",cal0->GetChisquare()/cal0->GetNDF()));}
      if(ch==1) {info->AddText(Form("#chi^{2}/ndf = %.2f",cal1->GetChisquare()/cal1->GetNDF()));}
      if(ch==2) {info->AddText(Form("#chi^{2}/ndf = %.2f",cal2->GetChisquare()/cal2->GetNDF()));}
      if(ch==3) {info->AddText(Form("#chi^{2}/ndf = %.2f",cal3->GetChisquare()/cal3->GetNDF()));}
      
      if(ch==4) {info2->AddText(Form("#chi^{2}/ndf = %.2f",cal4->GetChisquare()/cal4->GetNDF()));}
      if(ch==5) {info2->AddText(Form("#chi^{2}/ndf = %.2f",cal5->GetChisquare()/cal5->GetNDF()));}
      if(ch==6) {info2->AddText(Form("#chi^{2}/ndf = %.2f",cal6->GetChisquare()/cal6->GetNDF()));}
      if(ch==7) {info2->AddText(Form("#chi^{2}/ndf = %.2f",cal7->GetChisquare()/cal7->GetNDF()));}
      
      if(ch<=3) info->Draw("same");
      if(ch>3) info2->Draw("same");
    } 
  
  c2->SaveAs(Form("./Plots/%s/Vov_resolution_individual_channels_%d_%d.png",dir,run_i,run_f));
}
