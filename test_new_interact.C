#include "TH1.h"
#include <iostream>
#include "TMath.h"
{
  //  TFile* datafile = TFile::Open("test_merged.root");
  TFile* datafile1 = TFile::Open("ZnnG_histos_above0p5_Pt_SF.root");
  TFile *f1 = new TFile ("data.root","recreate");
  f1->cd();
  //  TH1F* h1 =dynamic_cast<TH1F*>(datafile->Get("histo_Zggh3Z0p0h4Z0p0"));  //N(0)  
  /*  TH1F* h2 =dynamic_cast<TH1F*>(datafile->Get("histo_aTGC_ZZg_h30p0_h40p0_interact_JESUp"));  //N(0)  
  TH1F* h3 =dynamic_cast<TH1F*>(datafile->Get("histo_aTGC_ZZg_h30p0_h40p0_interact_JESDown"));  //N(0)  
  TH1F* h4 =dynamic_cast<TH1F*>(datafile->Get("histo_aTGC_ZZg_h30p0_h40p0_interact_PESUp"));  //N(0)  
  TH1F* h5 =dynamic_cast<TH1F*>(datafile->Get("histo_aTGC_ZZg_h30p0_h40p0_interact_PESDown"));  //N(0)  
  */
  //systematic for data and background for second file 
  TH1F* h6 =dynamic_cast<TH1F*>(datafile1->Get("data_obs"));
  /*  TH1F* h7 =dynamic_cast<TH1F*>(datafile->Get("histo_jetfake"));
  TH1F* h8 =dynamic_cast<TH1F*>(datafile->Get("histo_jetfake_errUp"));
  TH1F* h9 =dynamic_cast<TH1F*>(datafile->Get("histo_jetfake_errDown"));
  TH1F* h10 =dynamic_cast<TH1F*>(datafile->Get("histo_elefake"));
  //systematic for mc background for scond file 
  TH1F* h11 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG"));
  TH1F* h12 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG_uncorrected"));
  TH1F* h13 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG_JESUp"));
  TH1F* h14 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG_JESDown"));
  TH1F* h15 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG_PESUp"));
  TH1F* h16 =dynamic_cast<TH1F*>(datafile->Get("histo_ZNuNuG_PESDown"));
  */
  // std::cout<<"SM integral "<<h1->Integral()<<endl;
  
  const int nBins = h6->GetXaxis()->GetNbins();
  Float_t PtBins[6]={225.,250., 300., 400., 600., 1000.0};
      //  TH1F* SM1 = new TH1F("sm","standardmodel",10,0,1200);                                                                                                  
  
  //  TH1F* SM1 = new TH1F("sm","standardmodel",5,PtBins);                                                                                                  
  //TH1F* SMJESUp = new TH1F("sm_JESUp","standardmodel",5,PtBins);                                                                                                        
  //  TH1F* SMJESDown = new TH1F("sm_JESDown","standardmodel",5,PtBins);                                                                                                  
  // TH1F* SMPESUp = new TH1F("sm_PESUp","standardmodel",5,PtBins);
  //TH1F* SMPESDown = new TH1F("sm_PESDown","standardmodel",5,PtBins);                                                                                                    
  TH1F* data_obs = new TH1F("data_obs","data_observa",5,PtBins);








  for(int i=1; i <= nBins; i++)
    {
      //  Float_t a=h1->GetBinContent(i);  //Sm N(0)
      //  Float_t b=h2->GetBinContent(i); //N(h3,0)
      //Float_t c=h3->GetBinContent(i); //N(-h3,0)
      // Float_t d=h4->GetBinContent(i); //N(0,h4)
      //Float_t e=h5->GetBinContent(i); //N(0,-h4)
      Int_t data=h6->GetBinContent(i);

      // float  sm = a; //N(0)
      //float sm_jesup= b;
      // float sm_jesdown= c;
      // float sm_pesup= d;
      // float sm_pesdown= e;

      //SM1->SetBinContent(i,a);   
      //SMJESUp->SetBinContent(i,b); 
      //SMJESDown->SetBinContent(i,c);
      // SMPESUp->SetBinContent(i,d);  
      //SMPESDown->SetBinContent(i,e);
      data_obs->SetBinContent(i,data); 
    }


 
  //  TH1F* SM1 = (TH1F*)SM->Clone();
  
  
  // TFile *f1 = new TFile ("signal_ATGC_locally.root","new");

  //  cout<<SM1->Integral()<<"StandardModel"<<endl;
  //  cout<<SMJESUp->Integral()<<"Jesup"<<endl;
  //cout<<SMJESDown->Integral()<<"jesdown"<<endl;
  // cout<<SMPESUp->Integral()<<"pesup"<<endl;
  // cout<<SMPESDown->Integral()<<"pesdown"<<endl;
  cout<<data_obs->Integral()<<"dataobservat"<<endl;
  data_obs->Write();
  //  SM1->Write();
  // SMJESUp->Write();  
  // SMJESDown->Write();  
  // SMPESUp->Write();  
  // SMPESDown->Write();  

  f1->Close();
  


}
