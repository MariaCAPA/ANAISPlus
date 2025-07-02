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

void DC_stabilization (int run)
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
    
    TFile *frun1;
    
    if(run==238) 
    {
    if(run<10) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run+1),"READ");
    if(run>=10 && run<100) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run+1),"READ");
    if(run>=100) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run+1),"READ");
    }
    
    else
    {
    if(run<10) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=10 && run<100) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run),"READ");
    if(run>=100) frun1 = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run),"READ");
    }
    
    TTree * pulses1 = (TTree*)frun1->Get("ta");
    
    TFile *frun2;
    
    if(run<10) frun2 = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run+1),"READ");
    if(run>=10 && run<100) frun2 = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run+1),"READ");
    if(run>=100) frun2 = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run+1),"READ");
    
    TTree * pulses2 = (TTree*)frun2->Get("ta");
    
    TFile *frun3;
    
    if(run<10) frun3 = new TFile(Form("../analyzed/ANAISplusLevel1/run0000%d.mid.lz4_analyzed_1.root",run+2),"READ");
    if(run>=10 && run<100) frun3 = new TFile(Form("../analyzed/ANAISplusLevel1/run000%d.mid.lz4_analyzed_1.root",run+2),"READ");
    if(run>=100) frun3 = new TFile(Form("../analyzed/ANAISplusLevel1/run00%d.mid.lz4_analyzed_1.root",run+2),"READ");
    
    TTree * pulses3 = (TTree*)frun3->Get("ta");

    TH1F* h1 = new TH1F("h1","Dark current",400,0,600);
    TH1F* h2 = new TH1F("h2","Dark current",400,0,600);
    TH1F* h3 = new TH1F("h3","Dark current",400,0,600);
    
    pulses1->Draw("high0>>h1");
    pulses2->Draw("high0>>h2");
    pulses3->Draw("high0>>h3");
    
    TCanvas *c1 = new TCanvas("c1","The Graph example",10,10,700,500);
    
    c1->SetGrid();
    c1->SetLogy();
    
    h1->SetLineColor(kRed-3);
    h2->SetLineColor(kMagenta-7);
    h3->SetLineColor(kCyan-3);
    
    h1->GetXaxis()->SetTitle("Height (a.u)");
    
    if(run==100) h1->SetTitle("Dark current [T_{A} = 30 K]");
    if(run==125) h1->SetTitle("Dark current [T_{A} = 80 K]");
    if(run==151) h1->SetTitle("Dark current [T_{A} = 130 K]");
    if(run==180) h1->SetTitle("Dark current [T_{A} = 180 K]");
    if(run==209) h1->SetTitle("Dark current [T_{A} = 230 K]");
    if(run==238) h1->SetTitle("Dark current [T_{A} = 280 K]");
    
    h1->Draw();
    h2->Draw("same");
    h3->Draw("same");
    
    TLegend *legend4 = new TLegend(0.70,0.72,0.87,0.87); //Tama駉 de la leyenda (x0,y0,x1,y1)
    legend4->SetTextSize(0.04); //Tama駉 del texto
    legend4->SetFillColor(kWhite); //Fondo de la leyenda

    legend4->AddEntry(h1,"10 min","l");
    legend4->AddEntry(h2,"60 min","l");
    legend4->AddEntry(h3,"120 min","l");

    legend4->Draw("same");
    
    c1->SaveAs(Form("./Plots/Dark_current/C_stabilization_%d.png",run));
}
