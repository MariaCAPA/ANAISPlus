// Example Program for converting MIDAS format to ROOT format.
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TMath.h>
#include <vector>
#include <TH1F.h>
#include <TStyle.h> //añadido

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
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
*/

//#include "midasio.h"

const int VMEBUS_BOARDNO = 0;

double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples);
double bsl_complex(TH1F* hAux, int nPtsBsl, int nSamples, int sIni, int sWin, double *rms);

int main(int argc, char *argv[])
{
  //We load the input file (waveforms)
  const char * file = argv[1];
  
  std::string input_file_name = Form("../analyzed/ANAISplusWaveforms/%s_waveforms.root",file);
  
  TFile * input_file = new TFile(input_file_name.c_str(), "read");
  
  TTree * Waveforms = (TTree*)input_file->Get("ta");
  
  //
  //Analysis configuration
  //
  
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
 
  nEvents=Waveforms->GetEntries();
 
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
 
  //
  //
  //
  
  //
  //Load input waveforms
  //
  
  TH1F ** hP = new TH1F * [nChannels]; //Pulses histograms

  for(int j = 0; j < nChannels; j++) 
  {
    hP[j] = new TH1F(Form("hP_%d",j),Form("hP_%d",j),nSamples,0,nSamples);
    Waveforms->SetBranchAddress(Form("pulse%d",j),&hP[j]);
  }
  
  //
  //sIni and sWin computation (metInt = true) We use channel 0
  //
  
  double bsl_provisional=0;
  double rms_provisional;
  
  TH1F * mean_waveform = new TH1F("meanhP","meanhP",nSamples,0,nSamples);
  
  int sFin;
  
  if(metInt==true)
  {
    sIni=nSamples;
    sWin=0;
    
    for(int i=0; i<nProm; i++)
    { 
      Waveforms->GetEntry(i);
      
      for(int j=0; j<nChannels; j++)
      {
       if(metBsl==true) bsl_provisional=bsl_simple(hP[j],nPtsBsl,sIni,&rms_provisional,nSamples);
       else bsl_provisional=bsl_complex(hP[j],nPtsBsl,nSamples,sIni,sWin,&rms_provisional);
      
       for(int k=1; k<=nSamples; k++) hP[j]->SetBinContent(k,(hP[j]->GetBinContent(k)-bsl_provisional)/amp);
       
       mean_waveform->Add(hP[j],1.0);
      }    
    }
    
    mean_waveform->Scale(1.0/double(nProm*nChannels));
    
    TFile * foutroot = new TFile(Form("../analyzed/ANAISplusLevel1/%s_average.root",file),"recreate");
    mean_waveform->Write();
    foutroot->Close();
    
    for(int k=1; k<=nSamples; k++)
    {
      if(mean_waveform->GetBinContent(k)>mean_waveform->GetMaximum()*thProm/100.) 
      {
        sIni=k;
        break;
      }
    }
    
    for(int k=nSamples; k>0; k--)
    {
      if(mean_waveform->GetBinContent(k)>mean_waveform->GetMaximum()*thProm/100.) 
      {
        sFin=k;
        break;
      }
    }
    
    sWin=sFin-sIni;
  }
  
  printf("Integration window: \n sIni = %d \n sWin = %d \n",sIni,sWin);
  
  //
  //
  //
  
 
  //
  //Output preparation
  // 
  std::string filenew = Form("../analyzed/ANAISplusLevel1/%s_analyzed_1.root",file);
  
  TFile * fnew = new TFile(filenew.c_str(), "recreate");
  
  TTree * t2 = new TTree("ta","");
  
  //Output variables
  double bsl[nChannels];
  double rms[nChannels];
  double t0[nChannels];
  double area[nChannels];
  double areaFix[nChannels];
  double tmax[nChannels];
  double high[nChannels];
  double minBsl[nChannels];
  double P1[nChannels];
  double areaTot;
  double areaFixTot;
  
  for(int j = 0; j < nChannels; j++) 
  {
    t2->Branch(Form("bsl%d",j),&bsl[j],Form("bsl%d/D",j));
    t2->Branch(Form("rms%d",j),&rms[j],Form("rms%d/D",j));
    t2->Branch(Form("t0%d",j),&t0[j],Form("t0%d/D",j));
    t2->Branch(Form("area%d",j),&area[j],Form("area%d/D",j));
    t2->Branch(Form("areaFix%d",j),&areaFix[j],Form("areaFix%d/D",j));
    t2->Branch(Form("tmax%d",j),&tmax[j],Form("tmax%d/D",j));
    t2->Branch(Form("high%d",j),&high[j],Form("high%d/D",j));
    t2->Branch(Form("minBsl%d",j),&minBsl[j],Form("minBsl%d/D",j));
    t2->Branch(Form("P1%d",j),&P1[j],Form("P1%d/D",j));
    
  }
  
  t2->Branch("sIni",&sIni,Form("sIni/I"));
  t2->Branch("sWin",&sWin,Form("sWin/I"));
  t2->Branch("areaTot",&areaTot,"areaTot/D");
  t2->Branch("areaFixTot",&areaFixTot,"areaFixTot/D");
  
  sIni=sIni;
  sWin=sWin;
  
  //
  //Variables computation
  //
  
  int ev=0;
  
  for(int i=0; i<nEvents; i++)
    { 
      Waveforms->GetEntry(i);
      
      for(int j=0; j<nChannels; j++)
      {
       for(int k = 2; k <= nSamples; k++) hP[j]->SetBinContent(k,hP[j]->GetBinContent(k-1)+0.5*(hP[j]->GetBinContent(k)-hP[j]->GetBinContent(k-1))); 
        
       //bsl & rms
       if(metBsl==true) bsl[j]=bsl_simple(hP[j],nPtsBsl,sIni,&rms[j],nSamples);
       else bsl[j]=bsl_complex(hP[j],nPtsBsl,nSamples,sIni,sWin,&rms[j]);
       
       for(int k=1; k<=nSamples; k++) hP[j]->SetBinContent(k,(hP[j]->GetBinContent(k)-bsl[j])/amp); //substraction of bsl
       
       //t0
       for(int k=sIni; k<=nSamples; k++)
       {
         if(hP[j]->GetBinContent(k)>threshold*rms[j])
         {
           t0[j]=double(k);
           break;
         }
       }
       
       //area
       area[j]=hP[j]->Integral(int(t0[j]),nSamples);
       
       //areaFix
       areaFix[j]=hP[j]->Integral(sIni,sIni+sWin);
       
       //high and tmax
       hP[j]->GetXaxis()->SetRangeUser(sIni,sIni+sWin);
       
       tmax[j]=hP[j]->GetMaximumBin();
       
       high[j]=hP[j]->GetMaximum();
       
       //P1
       P1[j]=hP[j]->Integral(hP[j]->GetMaximumBin(),sIni+sWin)/areaFix[j];
       
       //minBsl
       hP[j]->GetXaxis()->SetRangeUser(0,nPtsBsl);
       
       minBsl[j]=hP[j]->GetMinimum();
       
       //Restore X range
       hP[j]->GetXaxis()->SetRangeUser(0,nSamples); 
      }     
      
      areaTot=0;
      areaFixTot=0;
      
      for(int j = 0; j < nChannels; j++)
      {
        areaTot+=area[j];
        areaFixTot+=areaFix[j];
      }
      
      t2->Fill();
      
      if(ev%1000 == 0) {t2->AutoSave(); std::cout << ev << " " << file<< std::endl;}
      
      ev++;
    }
 
	t2->Write();

  fnew->Close();
  
  return 0;
}

double bsl_simple(TH1F* hAux, int nPtsBsl, int sIni, double *rms, int nSamples)
{
  double bslTot=0;
  double bslTotF=0;
  
  double rmsTot=0;
  double rmsTotF=0;
  
  int k;

  /*for(k = 0; k < nPtsBsl; k++) bslTot += hAux->GetBinContent(k+1)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) bslTotF += hAux->GetBinContent(nSamples-nPtsBsl+k)/double(nPtsBsl);
  
  for(k = 0; k < nPtsBsl; k++) rmsTot += (hAux->GetBinContent(k+1)-bslTot)*(hAux->GetBinContent(k+1)-bslTot)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) rmsTotF += (hAux->GetBinContent(nSamples-nPtsBsl+k)-bslTotF)*(hAux->GetBinContent(nSamples-nPtsBsl+k)-bslTotF)/double(nPtsBsl);
  
  rmsTot = sqrt(rmsTot);
  rmsTotF = sqrt(rmsTotF);

  if(rmsTot < rmsTotF) {*rms=rmsTot; return bslTot;}
  else { *rms=rmsTotF; return bslTotF;}*/
  
  for(k = 0; k < nPtsBsl; k++) bslTot += hAux->GetBinContent(k+1)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) bslTotF += hAux->GetBinContent(nPtsBsl+k+1)/double(nPtsBsl);
  
  for(k = 0; k < nPtsBsl; k++) rmsTot += (hAux->GetBinContent(k+1)-bslTot)*(hAux->GetBinContent(k+1)-bslTot)/double(nPtsBsl);
  for(k = 0; k < nPtsBsl; k++) rmsTotF += (hAux->GetBinContent(nPtsBsl+k+1)-bslTotF)*(hAux->GetBinContent(nPtsBsl+k+1)-bslTotF)/double(nPtsBsl);
  
  rmsTot = sqrt(rmsTot);
  rmsTotF = sqrt(rmsTotF);

  if(bslTot > bslTotF) {*rms=rmsTot; return bslTot;}
  else {*rms=rmsTotF; return bslTotF;}
}

double bsl_complex(TH1F* hAux, int nPtsBsl, int nSamples, int sIni, int sWin, double *rms)
{
  double rmsMin = 1e6;
  double bslRmsMin = 0;
  double bslCambiante = 0;
  double rmsCambiante = 0;

  for(int l = 0; l < nSamples - nPtsBsl; l++)
  {
   if(l < sIni - nPtsBsl || l > sIni + sWin)
   {
     bslCambiante = 0;
     rmsCambiante = 0;

     for(int k = l; k < l + nPtsBsl; k++) bslCambiante += hAux->GetBinContent(k+1)/double(nPtsBsl);
     for(int k = l; k < l + nPtsBsl; k++) rmsCambiante += (hAux->GetBinContent(k+1)-bslCambiante)*(hAux->GetBinContent(k+1)-bslCambiante)/double(nPtsBsl);
     rmsCambiante = sqrt(rmsCambiante);
     if(rmsCambiante < rmsMin)
     {
      rmsMin = rmsCambiante;
      bslRmsMin = bslCambiante;
     }
   }
  }
  
  *rms=rmsMin;

  return bslRmsMin;
}
