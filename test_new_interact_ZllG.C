#include "TH1.h"
#include <iostream>
#include "TMath.h"
{
 
  TFile* datafile = TFile::Open("ZnnG_histos_above0p5_Pt_SF.root");
  TFile *f1 = new TFile ("ZllG_combined.root","recreate");
  f1->cd();
  TH1F* h1 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_combined"));  //N(0)  
  TH1F* h2 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_JESUp_combined"));  //N(0)  
  TH1F* h3 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_JESDown_combined"));  //N(0)  
  TH1F* h4 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_PESUp_combined"));  //N(0)  
  TH1F* h5 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_PESDown_combined"));  //N(0)  
  TH1F* h6 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_phoscaleUp_combined"));  //N(0)  
  TH1F* h7 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_phoscaleDown_combined"));  //N(0)  
  TH1F* h8 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_straightUp_combined"));  //N(0)  
  TH1F* h9 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_straightDown_combined"));  //N(0)  
  TH1F* h10 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_twistedUp_combined"));  //N(0)  
  TH1F* h11 =dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_twistedDown_combined"));  //N(0)  
  TH1F* h12=dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_gammaUp_combined"));  //N(0)  
  TH1F* h13=dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_gammaDown_combined"));  //N(0)  
  TH1F* h14=dynamic_cast<TH1F*>(datafile->Get("histo_ZllG_qcdscale_combined"));  //N(0)  
 

  //systematic for data and background for second file 
  //  TH1F* h6 =dynamic_cast<TH1F*>(datafile1->Get("data_obs"));
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
  //  std::cout<<"SM integral "<<h1->Integral()<<endl;
  
  const int nBins = h1->GetXaxis()->GetNbins();
  Float_t PtBins[6]={225.,250., 300., 400., 600., 1000.0};
      //  TH1F* ZllG1 = new TH1F("sm","standardmodel",10,0,1200);                                                                                                  
  
  TH1F* ZllG = new TH1F("ZllG","standardmodel",5,PtBins);                                                                                                  
  TH1F* ZllGJESUp = new TH1F("ZllG_JESUp","standardmodel",5,PtBins);                                                                                                      TH1F* ZllGJESDown = new TH1F("ZllG_JESDown","standardmodel",5,PtBins);                                                                                                  TH1F* ZllGPESUp = new TH1F("ZllG_PESUp","standardmodel",5,PtBins);
  TH1F* ZllGPESDown = new TH1F("ZllG_PESDown","standardmodel",5,PtBins);                                                                                                
  TH1F* ZllGphoscaleUp = new TH1F("ZllG_phoscaleUp","standardmodel",5,PtBins);                                                                                          
  TH1F* ZllGphoscaleDown = new TH1F("ZllG_phoscaleDown","standardmodel",5,PtBins);                                                                                        TH1F* ZllGstraightUp = new TH1F("ZllG_straightUp","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllGstraightDown = new TH1F("ZllG_straightDown","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllG_twistedUp = new TH1F("ZllG_twistedUp","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllG_twistedDown = new TH1F("ZllG_twistedDown","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllG_gammaUp = new TH1F("ZllG_gammaUp","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllG_gammaDown = new TH1F("ZllG_gammaDown","standardmodel",5,PtBins);                                                                                     
  TH1F* ZllG_qcdscale = new TH1F("ZllG_qcdscale","standardmodel",5,PtBins);                                                                                     

  for(int i=1; i <= nBins; i++)
    {
      Float_t a=h1->GetBinContent(i);  //Sm N(0)
      Float_t b=h2->GetBinContent(i); //N(h3,0)
      Float_t c=h3->GetBinContent(i); //N(-h3,0)
      Float_t d=h4->GetBinContent(i); //N(0,h4)
      Float_t e=h5->GetBinContent(i); //N(0,-h4)
      Float_t f=h6->GetBinContent(i); //N(0,-h4)
      Float_t g=h7->GetBinContent(i);
      Float_t h=h8->GetBinContent(i);      
      Float_t I=h9->GetBinContent(i);
      Float_t J=h10->GetBinContent(i);
      Float_t K=h11->GetBinContent(i);
      Float_t l=h12->GetBinContent(i);
      Float_t m=h13->GetBinContent(i);
      Float_t n=h14->GetBinContent(i);

      float  sm = a; //N(0)
      float sm_jesup= b;
      float sm_jesdown= c;
      float sm_pesup= d;
      float sm_pesdown= e;

      ZllG->SetBinContent(i,a);   
      ZllGJESUp->SetBinContent(i,b); 
      ZllGJESDown->SetBinContent(i,c);
      ZllGPESUp->SetBinContent(i,d);  
      ZllGPESDown->SetBinContent(i,e);
      ZllGphoscaleUp->SetBinContent(i,f);
      ZllGphoscaleDown->SetBinContent(i,g);
      ZllGstraightUp->SetBinContent(i,h);
      ZllGstraightDown->SetBinContent(i,I);
      ZllG_twistedUp->SetBinContent(i,J);
      ZllG_twistedDown->SetBinContent(i,K);
      ZllG_gammaUp->SetBinContent(i,l);
      ZllG_gammaDown->SetBinContent(i,m);
      ZllG_qcdscale->SetBinContent(i,n);
    }
  ZllG->Write();
  ZllGJESUp->Write();
  ZllGJESDown->Write();
  ZllGPESUp->Write();
  ZllGPESDown->Write();
  ZllGphoscaleUp->Write();
  ZllGphoscaleDown->Write();
  ZllGstraightUp->Write();
  ZllGstraightDown->Write();
  ZllG_twistedUp->Write();
  ZllG_twistedDown->Write();
  ZllG_gammaUp->Write();
  ZllG_gammaDown->Write();
  ZllG_qcdscale->Write();
  f1->Close();
  


}
