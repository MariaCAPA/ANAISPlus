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

void crosstalk_stimation (int run, int ch=0, int Height_nbins=300)
{
    gStyle->SetImageScaling(3.);
    gStyle->SetHistLineWidth(3);
    gStyle->SetTitleFont(102,"XYZ");
    gStyle->SetTitleFont(102,"");
    gStyle->SetTitleSize(0.07,"");
    gStyle->SetLegendFont(102);
    gStyle->SetLabelFont(102,"XYZ");
    gStyle->SetTextFont(102);
    gStyle->SetOptStat(0);
    
    TFile *frun;
    
    if(run<10) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=10 && run<100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=100) frun = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");

    TTree * pulses = (TTree*)frun->Get("ta");
    
    //
    //Definimos propiedades histograma de altura
    //

    TH1F *Height_aux = new TH1F("Height_aux","Height_aux",10000,100,50);

    pulses->Draw(Form("areaFix%d >> Height_aux",ch),Form("areaFix%d>0",ch),"");

    int Height_max = Height_aux->GetXaxis()->GetXmax();
    int Height_min = Height_aux->GetXaxis()->GetXmin();

    double Height_bin_size = (Height_max-Height_min)/10000.;

    int Height_final_bin = Height_aux->FindLastBinAbove(1);
    int Height_initial_bin = Height_aux->FindFirstBinAbove(1);

    Height_max = Height_final_bin*Height_bin_size+Height_min;
    Height_min = Height_initial_bin*Height_bin_size+Height_min;
    
    //Height_max = 700;

    TH1F *Height = new TH1F("Height", "Height",Height_nbins,0,Height_max);

    pulses->Draw(Form("areaFix%d >> Height",ch),Form("areaFix%d>0",ch),"");
    
    //
    //Se buscan cortes en altura
    //

    TCanvas *c1 = new TCanvas("c1","The Graph example",10,10,700,500);

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
    	line = new TLine(Height_cut[k],0,Height_cut[k],Height->GetBinContent(Height->GetMaximumBin())); //Coordenadas de la línea
    	line->SetLineWidth(4.); //Tamaño de la línea
    	line->SetLineColor(kRed); //Color de la línea
    	line->Draw("same"); //Para dibujar la línea
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
    TF1* gaus_3 = new TF1("gaus 2","gaus",0,Height_max);

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
    //Análisis
    //
    TH1F *cross_phe = new TH1F(Form("cross_phe%d",run),Form("cross_phe%d",run),Height_nbins,0,Height_max);
    for(i=gaus_1->GetParameter(1)/(double(Height_max)/double(Height_nbins));i<=Height_nbins;i++) cross_phe->SetBinContent(i,Height->GetBinContent(i));
    
    double int_phe, err_phe;
    double int_correct_phe, err_correct_phe;
    double int_1_phe, err_1_phe;
    
    int_phe=cross_phe->IntegralAndError(0,Height_nbins,err_phe);
    int_correct_phe=gaus_1->Integral(gaus_1->GetParameter(1),Height_max);
    err_correct_phe=gaus_1->IntegralError(gaus_1->GetParameter(1),Height_max);
    int_1_phe=gaus_1->Integral(0,Height_max);
    err_1_phe=gaus_1->IntegralError(0,Height_max);
    
    int_phe=int_phe*(double(Height_max)/double(Height_nbins));
    err_phe=err_phe*(double(Height_max)/double(Height_nbins));
    
    double CTprob, CTprob_err;
    
    double err_sustr;
    err_sustr=sqrt(err_phe*err_phe+err_correct_phe*err_correct_phe);
    
    CTprob=(int_phe-int_correct_phe)/(int_phe-int_correct_phe+int_1_phe)*100;
    CTprob_err=100*sqrt(int_1_phe*int_1_phe*err_sustr*err_sustr+(int_phe-int_correct_phe)*(int_phe-int_correct_phe)*err_1_phe*err_1_phe)/(int_phe-int_correct_phe+int_1_phe)/(int_phe-int_correct_phe+int_1_phe);
    
    cout << CTprob << endl;
    
    //
    //Plot
    //
    c1->SetLogy();
    
    gaus_1->SetLineWidth(4.);
    gaus_1->SetNpx(500);
    gaus_1->SetLineColor(kRed+1);
    gaus_1->SetFillColor(kRed-7);
    gaus_1->SetFillStyle(1001);

    gaus_2->SetLineWidth(4.);
    gaus_2->SetLineColor(3);

    gaus_tot->SetLineWidth(4.);
    gaus_tot->SetLineColor(1);
    
    cross_phe->SetLineColor(1);
    cross_phe->SetFillColor(15);

    Height->GetXaxis()->SetRangeUser(0,400);
    Height->GetXaxis()->SetTitle("Height (a.u)");
    //Height->GetYaxis()->SetRangeUser(0,700);
    
    Height->GetXaxis()->SetRangeUser(0,Height_max);

    Height->SetTitle(Form("Height spectrum Run %d",run));
    Height->SetLineColor(1);

    Height->Draw("HIST");
    cross_phe->Draw("same");
    //gaus_tot->Draw("same");
    gaus_1->Draw("same");
    
    TLegend *legend = new TLegend(0.65,0.7,0.87,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend->SetTextSize(0.04); //Tamaño del texto
    legend->SetFillColor(kWhite); //Fondo de la leyenda

    legend->AddEntry(gaus_1,"1 phe","f");
    legend->AddEntry(cross_phe,"Cross talk","f");

    legend->Draw("same");
    
    TPaveText *pv2 = new TPaveText(0.57,0.57,0.87,0.67,"NDC"); //Coordenadas del texto (x0,y0,x1,y1)
    pv2->SetBorderSize(1); //Tamaño del borde
    pv2->SetFillColor(0); //Color del fondo
    pv2->SetTextSize(.04); //Tamaño del texto
    pv2->SetTextColor(kBlack); //Color del texto
    pv2->SetMargin(0.01); //Margen del texto
    pv2->AddText(Form("CT = %.2f #pm %.2f %%",CTprob,CTprob_err)); //Texto a añadir
    pv2->SetTextFont(102); //Fuente del texto
    pv2->SetTextAlign(22); //Alineación del texto 22=centrada en horizontal y vertical
    pv2->Draw(); //Para dibujar el cuadro de texto
    
    c1->SaveAs(Form("./Plots/Dark_current/CT_Stimation_%d.png",run));    
    
    //
    //Más análisis
    //
    
    //Eventos con 0, 1, 2 phe secundarios
    gaus_0->SetLineWidth(4.);
    gaus_0->SetLineColor(kCyan+1);
    
    gaus_1->SetLineWidth(4.);
    gaus_1->SetLineColor(kRed+1);
    gaus_1->SetFillStyle(0);
    
    gaus_2->SetLineWidth(4.);
    gaus_2->SetLineColor(kGreen-3);
    
    gaus_3->SetLineWidth(4.);
    gaus_3->SetLineColor(kMagenta-7);
    
    Height->Draw("HIST");
    gaus_0->Draw("same");
    gaus_1->Draw("same");
    gaus_2->Draw("same");
    gaus_3->Draw("same");
    
    TLegend *legend2 = new TLegend(0.72,0.63,0.87,0.85); //Tamaño de la leyenda (x0,y0,x1,y1)
    legend2->SetTextSize(0.04); //Tamaño del texto
    legend2->SetFillColor(kWhite); //Fondo de la leyenda

    legend2->AddEntry(gaus_0,"Noise","l");
    legend2->AddEntry(gaus_1,"1 phe","l");
    legend2->AddEntry(gaus_2,"2 phe","l");
    legend2->AddEntry(gaus_3,"3 phe","l");

    legend2->Draw("same");
    
    c1->SaveAs(Form("./Plots/Dark_current/CT_phe_fits_%d.png",run));
    
    TGraphErrors* CT_graph = new TGraphErrors();
    
    CT_graph->SetPoint(0,0,gaus_1->Integral(0,Height_max));
    CT_graph->SetPoint(1,1,gaus_2->Integral(0,Height_max));
    CT_graph->SetPoint(2,2,gaus_3->Integral(0,Height_max));
    
    CT_graph->SetPointError(0,0,gaus_1->IntegralError(0,Height_max));
    CT_graph->SetPointError(1,0,gaus_2->IntegralError(0,Height_max));
    CT_graph->SetPointError(2,0,gaus_3->IntegralError(0,Height_max));
    
    CT_graph->SetMarkerStyle(8);
    CT_graph->SetMarkerSize(1.5);
    CT_graph->SetTitle(Form("Cross-talk phe production Run %d",run));
    CT_graph->GetXaxis()->SetTitle("Secondary phe");
    CT_graph->GetYaxis()->SetTitle("Number of events");
    
    CT_graph->GetXaxis()->SetNdivisions(3);
    
    TF1* CT_model_pop = new TF1("CT_model_pop","[0]*([1]^x-[1]^(x+1))",0,4);
    
    TCanvas *c2 = new TCanvas("c2","The Graph example",10,10,700,700);
    
    CT_graph->Fit("CT_model_pop","0");
    
    CT_model_pop->SetLineStyle(9);
    CT_model_pop->SetLineWidth(4);
    CT_model_pop->SetLineColor(kRed-3);
    
    CT_graph->Draw("AP");
    CT_model_pop->Draw("same");
    
    TPaveText *pv3 = new TPaveText(0.45,0.69,0.87,0.85,"NDC"); //Coordenadas del texto (x0,y0,x1,y1)
    pv3->SetBorderSize(1); //Tamaño del borde
    pv3->SetFillColor(0); //Color del fondo
    pv3->SetTextSize(.03); //Tamaño del texto
    pv3->SetTextColor(kBlack); //Color del texto
    pv3->SetMargin(0.01); //Margen del texto
    pv3->AddText(Form("CT [1] = %.2f #pm %.2f %%",CTprob,CTprob_err)); //Texto a añadir
    pv3->AddText(Form("CT [2] = %.2f #pm %.2f %%",CT_model_pop->GetParameter(1)*100,CT_model_pop->GetParError(1)*100));
    pv3->SetTextFont(102); //Fuente del texto
    pv3->SetTextAlign(22); //Alineación del texto 22=centrada en horizontal y vertical
    pv3->Draw(); //Para dibujar el cuadro de texto   
    
    c2->SaveAs(Form("./Plots/Dark_current/CT_phe_populations_%d.png",run));
    /*
    //Poblaciones de cada número de phe
    cout << gaus_2->Integral(0,1000)/gaus_1->Integral(0,1000) << endl;
    cout << gaus_3->Integral(0,1000)/gaus_2->Integral(0,1000) << endl;
    cout << gaus_3->Integral(0,1000) << endl;
    
    //Beta (promedio de secundarios)
    Height->Add(gaus_0,-1);
    
    Height->GetXaxis()->SetLimits(-gaus_0->GetParameter(1), Height->GetXaxis()->GetXmax()-gaus_0->GetParameter(1));
    
    Height->GetXaxis()->SetRangeUser(0,(Height->GetXaxis()->GetXmax()-gaus_0->GetParameter(1)));
    
    Height->GetXaxis()->SetLimits(0,(Height->GetXaxis()->GetXmax()/(gaus_2->GetParameter(1)-gaus_1->GetParameter(1))));
    
    TCanvas *c3 = new TCanvas("c3","The Graph example",10,10,700,500);
    
    Height->SetTitle(Form("Phe spectrum Run %d",run));
    
    Height->GetXaxis()->SetTitle("nphe");
    
    Height->Draw("HIST");
    
    TPaveText *pv4 = new TPaveText(0.5,0.75,0.87,0.85,"NDC"); //Coordenadas del texto (x0,y0,x1,y1)
    pv4->SetBorderSize(1); //Tamaño del borde
    pv4->SetFillColor(0); //Color del fondo
    pv4->SetTextSize(.03); //Tamaño del texto
    pv4->SetTextColor(kBlack); //Color del texto
    pv4->SetMargin(0.01); //Margen del texto
    pv4->AddText(Form("Average nphe = %.2f #pm %.2f phe",Height->GetMean(),Height->GetMeanError())); //Texto a añadir
    pv4->SetTextFont(102); //Fuente del texto
    pv4->SetTextAlign(22); //Alineación del texto 22=centrada en horizontal y vertical
    pv4->Draw(); //Para dibujar el cuadro de texto   
    
    c3->SaveAs(Form("./Plots/Dark_current/CT_average_nphe_%d.png",run));
    
    //Dark current rate stimation
    double time_stamp=0;
    double time_stamp_i=0;
    double time_stamp_f=0;
    
    int time_control=0;
    
    pulses->SetBranchAddress("timeFile",&time_stamp);
    
    pulses->GetEntry(0);
    time_stamp_i=time_stamp*1e-9;
    
    pulses->GetEntry(pulses->GetEntries()-1);
    time_stamp_f=time_stamp*1e-9;
    
    double DCR=0;
    
    DCR=Height->Integral()/(time_stamp_f-time_stamp_i);
    
    cout << "Time = " << time_stamp_f << " " << time_stamp_i << " s" << endl;
    cout << "DCR = " << DCR << " Hz" << endl;
    
    //Dark current rate stimation v.2
    double Height_value;
    double time_value;
    double time_i, delta_t;
    
    pulses->SetBranchAddress(Form("high%d",ch),&Height_value);
    pulses->SetBranchAddress("timeFile",&time_value);
    
    TH1F* DC_time = new TH1F("DC_time",Form("Phe time difference [Run %d]",run),1000,0,0.004);

    for(int i=0; i<pulses->GetEntries(); i++)
    {
      pulses->GetEntry(i);
      
      //if(Height_value>Height_cut[0])
      if(Height_value>gaus_1->GetParameter(1)-gaus_1->GetParameter(2))
      {
       if(time_control==1)
        {
          delta_t=(time_value-time_i)/1.e+9;
          
          if(delta_t>0) DC_time->Fill(delta_t);
          
          time_control=0;
        }
        
        if(time_control==0) {time_i=time_value;  time_control=1;} 
      }
    }
    
    TCanvas *c4 = new TCanvas("c4","The Graph example",10,10,700,500);
    
    c4->SetGrid();
    c4->SetLogy();
    
    DC_time->SetLineColor(kBlack);
    
    DC_time->GetXaxis()->SetTitle("Time (s)");
    
    DC_time->Draw();
    
    TF1* DC_fit = new TF1("DC_fit","[0]*exp(-x*[1])+[2]",0,0.004);
    
    if(run>400) DC_time->Fit(DC_fit,"0","",0.0005,0.004);
    if(run<400)  DC_time->Fit(DC_fit,"0","",0.0001,0.004);
    
    DC_fit->SetLineStyle(9);
    DC_fit->SetLineWidth(4);
    DC_fit->SetLineColor(kRed-3);
    
    DC_fit->Draw("same");
    
    TPaveText *pv5 = new TPaveText(0.57,0.72,0.87,0.87,"NDC"); //Coordenadas del texto (x0,y0,x1,y1)
    pv5->SetBorderSize(1); //Tamaño del borde
    pv5->SetFillColor(0); //Color del fondo
    pv5->SetTextSize(.04); //Tamaño del texto
    pv5->SetTextColor(kBlack); //Color del texto
    pv5->SetMargin(0.01); //Margen del texto
    pv5->AddText(Form("DCR = %.0f #pm %.0f Hz",DC_fit->GetParameter(1),DC_fit->GetParError(1))); //Texto a añadir
    pv5->AddText(Form("#chi^{2}/ndf = %.2f",DC_fit->GetChisquare()/DC_fit->GetNDF())); //Texto a añadir
    pv5->SetTextFont(102); //Fuente del texto
    pv5->SetTextAlign(22); //Alineación del texto 22=centrada en horizontal y vertical
    pv5->Draw(); //Para dibujar el cuadro de texto
    
    c4->SaveAs(Form("./Plots/Dark_current/DCR_estimation_%d.png",run));
    
    FILE *f_text;
    
    f_text=fopen(Form("./Plots/Dark_current/DCR_estimation.txt"),"a");
    
    fprintf(f_text,"%d %f %f\n",run,DC_fit->GetParameter(1),DC_fit->GetParError(1));

    fclose(f_text);  
    */   
}
