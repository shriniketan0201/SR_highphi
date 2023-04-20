#define postAnalyzer_cxx
#include "postAnalyzer_qcd.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
  {
    std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
    return 1;
  }
  postAnalyzer t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void postAnalyzer::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
  nTotal = 0;
  int nFilters, nHLT, nPhoCand, nWorstChIso, nMET170, nDphiPhoMET, nLeptonVeto, nNoisyCrystalVeto, nDphiJetsMET;
  nFilters = nHLT = nPhoCand = nWorstChIso = nMET170 = nDphiPhoMET = nLeptonVeto = nNoisyCrystalVeto = nDphiJetsMET = 0;

  std::vector<Int_t> runlist;
  runlist.clear();
  std::vector<Long64_t> eventlist;
  eventlist.clear();
  std::vector<Int_t> lumilist;
  lumilist.clear();

  std::vector<int> phoCand;
  phoCand.clear();

  std::vector<int> phoCand1;
  phoCand1.clear();


  std::vector<int> qcdden;
  qcdden.clear();

  std::vector<int> jetveto;
  jetveto.clear();

  int iphi = 41;
  int ieta = 5;

  bool debug=true;
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  //Look at up to maxEvents events, or all if maxEvents == -1.
  Long64_t nentriesToCheck = nentries;
  if (maxEvents != -1LL && nentries > maxEvents)
  nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;

  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  TStopwatch sw;
  sw.Start();
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++)
  {

    event_.clear();
    event_info.clear();

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //=1.0 for real data
    double event_weight=1.0;

    //Sequential cuts

    phoCand1   = getQcdden(225,1.4442,1);
    
    if(phoCand1.size()>0) 
      {
        jetveto = JetVetoDecision(phoCand1[0]);
      }
        
    //    if(metFilters==1536)
    if (!getMetFilter()) continue;  
    {
        // if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1))
    if(HLTPho>>11&1 == 1)
      {
        if(phoCand1.size() >0)
        {
	  //          Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCand1[0]]/TMath::CosH((*phoSCEta)[phoCand1[0]]));
          Float_t uncorrectedPhoEt = ((*phoEt)[phoCand1[0]]);
          event_weight =FakeRatePt(uncorrectedPhoEt,0);
          double event_weight_sidebandUp =FakeRatePt(uncorrectedPhoEt,1);
          double event_weight_sidebandDown =FakeRatePt(uncorrectedPhoEt,2);
          double event_weight_METUp =FakeRatePt(uncorrectedPhoEt,3);
          double event_weight_METDown =FakeRatePt(uncorrectedPhoEt,4);
          double event_weight_binningUp =FakeRatePt(uncorrectedPhoEt,5);
          double event_weight_binningDown =FakeRatePt(uncorrectedPhoEt,6);
          double event_weight_sieieLeft =FakeRatePt(uncorrectedPhoEt,7);
          double event_weight_sieieRight =FakeRatePt(uncorrectedPhoEt,8);
          double event_weight_templateUp =FakeRatePt(uncorrectedPhoEt,9);
          double event_weight_templateDown =FakeRatePt(uncorrectedPhoEt,10);
      //     if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCand1[0]]  - rho*EAchargedworst((*phoSCEta)[phoCand1[0]]) ), 0.0) < 1.37 )
      // {
        fillHistos(12,event_weight,phoCand1[0],jetveto);
        
        float phi_prime = 0.0;
        if (fabs(phoSCPhi->at(phoCand1[0])) < 3.14159/2.0)
          phi_prime = phoSCPhi->at(phoCand1[0]);
        else if (phoSCPhi->at(phoCand1[0]) > 3.14159/2.0 && phoSCPhi->at(phoCand1[0]) < 3.1416)
          phi_prime = 3.14159 - phoSCPhi->at(phoCand1[0]);
        else if (phoSCPhi->at(phoCand1[0]) < -3.14159/2.0 && phoSCPhi->at(phoCand1[0]) > -3.1416)
          phi_prime = -3.14159 - phoSCPhi->at(phoCand1[0]);
        else{
          cout<<"phoSCPhi can't be converted to phi_prime"<<endl;
          continue;
        }
        //if(pfMET>200)  //for both phi region 
	if(pfMET>200 && fabs(phi_prime) > 0.5)  //considering high phi region
          {
            fillHistos(13,event_weight,phoCand1[0],jetveto);
            if(DeltaPhi(phoPhi->at(phoCand1[0]),pfMETPhi)>0.5)
        {
          fillHistos(14,event_weight,phoCand1[0],jetveto);
            if(electron_veto_looseID(phoCand1[0],10) &&  muon_veto_looseID(phoCand1[0],10))
        {
          fillHistos(15,event_weight,phoCand1[0],jetveto);
          
          /*if(!(phoIEta->at(phoCand1[0])==-24 && phoIPhi->at(phoCand1[0])==141) && !(phoIEta->at(phoCand1[0])==4 && phoIPhi->at(phoCand1[0])==41) && !(phoIEta->at(phoCand1[0])==5 && phoIPhi->at(phoCand1[0])==41) && !(phoIEta->at(phoCand1[0])==1 && phoIPhi->at(phoCand1[0])==81) && !(phoIEta->at(phoCand1[0])==4 && phoIPhi->at(phoCand1[0])==21))
            {
              fillHistos(16,event_weight,phoCand1[0],jetveto);
	  */
	      double jetpt =0.0;
              for(int i=0;i<jetveto.size(); i++)
          {
            jetpt+=jetPt->at(jetveto[i]) ;
          }
              
              double minDPhiJetMET = TMath::Pi();
              double minDPhiJetMET_first4 = TMath::Pi();
              for(int j = 0; j < jetveto.size(); j++)
          {
            if(DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi) < minDPhiJetMET)
              {
                minDPhiJetMET = DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi);
                if(j < 4)
            minDPhiJetMET_first4 = DeltaPhi(jetPhi->at(jetveto[j]),pfMETPhi);
              }
          }



              if(dPhiJetMET_veto(jetveto)){
          std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
               
          if(jetveto.size()==2) std::cout<<"for size 2 run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
          if(jetveto.size()>=8) std::cout<<"for size 8 run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
          fillHistos(16,event_weight,phoCand1[0],jetveto);
          fillHistos(29,event_weight_sidebandUp,phoCand1[0],jetveto);
          fillHistos(30,event_weight_sidebandDown,phoCand1[0],jetveto);
          fillHistos(31,event_weight_METUp,phoCand1[0],jetveto);
          fillHistos(32,event_weight_METDown,phoCand1[0],jetveto);
          fillHistos(33,event_weight_binningUp,phoCand1[0],jetveto);
          fillHistos(34,event_weight_binningDown,phoCand1[0],jetveto);
          fillHistos(35,event_weight_sieieLeft,phoCand1[0],jetveto);
          fillHistos(36,event_weight_sieieRight,phoCand1[0],jetveto);
          fillHistos(37,event_weight_templateUp,phoCand1[0],jetveto);
          fillHistos(38,event_weight_templateDown,phoCand1[0],jetveto);
          fillHistos(39,1.0,phoCand1[0],jetveto);
	  //            if((phoSCRawE->at(phoCand1[0])/TMath::CosH(phoSCEta->at(phoCand1[0])))/pfMET < 1.4){
	  if(((*phoEt)[phoCand1[0]])/pfMET < 1.4){
              fillHistos(17,event_weight,phoCand1[0],jetveto);
              fillHistos(18,event_weight_sidebandUp,phoCand1[0],jetveto);
              fillHistos(19,event_weight_sidebandDown,phoCand1[0],jetveto);
              fillHistos(20,event_weight_METUp,phoCand1[0],jetveto);
              fillHistos(21,event_weight_METDown,phoCand1[0],jetveto);
              fillHistos(22,event_weight_binningUp,phoCand1[0],jetveto);
              fillHistos(23,event_weight_binningDown,phoCand1[0],jetveto);
              fillHistos(24,event_weight_sieieLeft,phoCand1[0],jetveto);
              fillHistos(25,event_weight_sieieRight,phoCand1[0],jetveto);
              fillHistos(26,event_weight_templateUp,phoCand1[0],jetveto);
              fillHistos(27,event_weight_templateDown,phoCand1[0],jetveto);
              fillHistos(28,1.0,phoCand1[0],jetveto);
            }
          
          
          
              }
              
            }
        
          }
      }
        // }
    }
      }
      }
    //qcd sequential cuts

  
      



    
    
    tree->Fill();
    
    if (jentry%reportEvery == 0)
      {
        std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
      }
  }

  if((nentriesToCheck-1)%reportEvery != 0)
    std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
  sw.Stop();
  std::cout<<"All events checked."<<std::endl;
  //Report
  std::cout << "RealTime : " << sw.RealTime() / 60.0 << " minutes" << std::endl;
    std::cout << "CPUTime  : " << sw.CpuTime()  / 60.0 << " minutes" << std::endl;
  std::cout << std::endl;
  std::cout<<"Total number of events: "<<nTotal<<std::endl;
  std::cout << "Number of events inspected: " << nTotal << std::endl;
  std::cout<<std::endl;
}

void postAnalyzer::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();

  Float_t PtBins[6]={225.,250., 300., 400., 600., 1000.0};
  Float_t MetBins[6]={200.,250., 300., 400., 600., 1000.0};
  Float_t dPhiJetMETBins[14]={0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.142};
  
  /* Float_t PtBins[7]={175.,200.,250., 300., 400., 600., 1000.0};
  Float_t MetBins[7]={170.,200.,250., 300., 400., 600., 1000.0};
  Float_t MTBins[10]={0.,200.,300.,400.,500.,600.,700.,800.,1000.,1200.};
  Float_t dPhiJetMETBins[14]={0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.1416};
  */
 //h_phoIEtaIPhi = new TH2F("h_phoIEtaIPhi","Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi->Sumw2();

  /*  h_dPhijetmet_min = new TH1F("h_dPhijetmet_min","h_dPhijetmet_min",40,0,3.2);h_dPhijetmet_min->Sumw2();
        h_dPhijetmet_min_first4 = new TH1F("h_dPhijetmet_min_first4","h_dPhijetmet_min_first4",40,0,3.2);h_dPhijetmet_min_first4->Sumw2();
  h_HT = new TH1F("h_HT","h_HT",50,0,1000);h_HT->Sumw2();
  h_HTMET = new TH2F("h_HTMET","h_HTMET",100,0,1000,86,140,1000);h_HTMET->Sumw2();
  h_njetMET = new TH2F("h_njetMET","h_njetMET",20,0,20,86,140,1000);h_njetMET->Sumw2();
  */
  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<45; i++)
  {
    char ptbins[100];
    sprintf(ptbins, "_%d", i);
    std::string histname(ptbins);
    h_nVtx[i] = new TH1F(("nVtx"+histname).c_str(), "nVtx",40,0,40);h_nVtx[i]->Sumw2();
    h_photon_Et[i] = new TH1F(("Photon_Et"+histname).c_str(), "Photon_Et",5,PtBins);h_photon_Et[i]->Sumw2();
    h_photon_Et_range[i] = new TH1F(("Photon_Et_range"+histname).c_str(), "Photon_Et",5,PtBins);h_photon_Et_range[i]->Sumw2();
    h_photon_eta[i] = new TH1F(("Photon_eta"+histname).c_str(), "Photon_eta",20,-1.4442,1.4442);h_photon_eta[i]->Sumw2();
    h_photon_SCEta[i] = new TH1F(("Photon_SCeta"+histname).c_str(), "Photon_SCeta",20,-1.4442,1.4442);h_photon_SCEta[i]->Sumw2();
    h_photon_phi[i] = new TH1F(("Photon_phi"+histname).c_str(), "Photon_phi", 64,0,3.2);h_photon_phi[i]->Sumw2();
    h_photon_SCPhi[i] = new TH1F(("Photon_SCphi"+histname).c_str(), "Photon_SCphi", 20,0,3.1416);h_photon_SCPhi[i]->Sumw2();
    h_photon_IDbit[i] = new TH1F(("Photon_ID_bit"+histname).c_str(), "Photon_ID_bit",8,0,8);h_photon_IDbit[i]->Sumw2();
    h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",5,MetBins);h_pfMET[i]->Sumw2();
    h_pfMET_300[i] = new TH1F(("h_pfMET_300"+histname).c_str(), "pfMET",25,0,300);h_pfMET_300[i]->Sumw2();
    h_dPhi[i] = new TH1F(("h_dPhi"+histname).c_str(),"h_dPhi",40,0,3.2);h_dPhi[i]->Sumw2();
    h_nJet[i] = new TH1F(("nJet"+histname).c_str(), "nJet",20,0,20);h_nJet[i]->Sumw2();
    h_leadingJetPt[i] = new TH1F(("leadingJetPt"+histname).c_str(),"leadingJetPt",30,20,1000);h_leadingJetPt[i]->Sumw2();
    h_leadingJetPt_300[i] = new TH1F(("leadingJetPt_300"+histname).c_str(),"leadingJetPt_300",25,0,300);h_leadingJetPt_300[i]->Sumw2();
    h_leadingJetEta[i] = new TH1F(("h_leadingJetEta"+histname).c_str(),"h_leadingJetEta",40,-1.4442,1.4442);h_leadingJetEta[i]->Sumw2();
    //    h_phoIEtaIPhi[i] = new TH2F(("h_phoIEtaIPhi"+histname).c_str(),"h_phoIEtaIPhi",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi[i]->Sumw2();
    h_PTMET[i] = new TH1F(("PTMET"+histname).c_str(),"P_{T}/Missing E_{T}",50,0,3);h_PTMET[i]->Sumw2();
    //    h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",9,MTBins);h_Mt[i]->Sumw2();
    h_min_dphijetmet[i] = new TH1F(("h_min_dphijetmet"+histname).c_str(),"h_min_dphijetmet",13,dPhiJetMETBins);h_min_dphijetmet[i]->Sumw2();
    h_pfMETsumEt[i] = new TH1F(("pfMETsumEt"+histname).c_str(),"pfMETsumEt",5,MetBins);
    h_METoverSqrtSumEt_extended[i] = new TH1F(("METoverSqrtSumEt_extended"+histname).c_str(),"METoverSqrtSumEt",30,0,30);
    h_METoverSqrtSumEt[i] = new TH1F(("METoverSqrtSumEt"+histname).c_str(),"METoverSqrtSumEt",30,0,10);
    h_r9[i] = new TH1F(("r9"+histname).c_str(),"r9",60,0,1);
  }
}

//Fill the sequential histos at a particular spot in the sequence
void postAnalyzer::fillHistos(int histoNumber, double event_weight,int index,std::vector<int> jets)
{
  Float_t uncorrectedPhoEt = ((*phoEt)[index]);
  h_photon_Et[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
  h_photon_Et_range[histoNumber]->Fill(uncorrectedPhoEt,event_weight);
  h_photon_eta[histoNumber]->Fill(phoEta->at(index),event_weight);
  h_photon_SCEta[histoNumber]->Fill(phoSCEta->at(index),event_weight);
  h_photon_phi[histoNumber]->Fill(phoPhi->at(index),event_weight);
  h_photon_SCPhi[histoNumber]->Fill(phoSCPhi->at(index),event_weight);
  h_pfMET[histoNumber]->Fill(pfMET,event_weight);
  h_pfMET_300[histoNumber]->Fill(pfMET,event_weight);
  double dPhi = DeltaPhi(phoPhi->at(index),pfMETPhi);
  h_dPhi[histoNumber]->Fill(dPhi,event_weight);
  h_nJet[histoNumber]->Fill(jets.size(),event_weight);
  // h_phoIEtaIPhi[histoNumber]->Fill(phoIPhi->at(index),phoIEta->at(index),event_weight);
  //h_PTMET[histoNumber]->Fill(phoEt->at(index)/pfMET,event_weight);
  h_PTMET[histoNumber]->Fill(uncorrectedPhoEt/pfMET,event_weight);
  //  h_Mt[histoNumber]->Fill(sqrt(2*(phoEt->at(index))*pfMET*(1-TMath::Cos(dPhi))),event_weight);
  h_nVtx[histoNumber]->Fill(nVtx);
  if(jets.size()>0){
    h_leadingJetPt[histoNumber]->Fill(jetPt->at(jets[0]),event_weight);
    h_leadingJetEta[histoNumber]->Fill(jetEta->at(jets[0]),event_weight);
    int max_njets = jets.size();
    if(jets.size() > 4)
      max_njets = 4;
    double min_dphijetmet = TMath::Pi();
    for(int i = 0; i < max_njets; i++)
      {
        double dphijetmet = DeltaPhi(jetPhi->at(jets[i]),pfMETPhi);
        if(dphijetmet < min_dphijetmet)
          min_dphijetmet = dphijetmet;
      }
    h_min_dphijetmet[histoNumber]->Fill(min_dphijetmet,event_weight);
  }
  h_pfMETsumEt[histoNumber]->Fill(pfMETsumEt,event_weight);
  h_METoverSqrtSumEt_extended[histoNumber]->Fill(pfMET/sqrt(pfMETsumEt),event_weight);
  h_METoverSqrtSumEt[histoNumber]->Fill(pfMET/sqrt(pfMETsumEt),event_weight);
  h_r9[histoNumber]->Fill(phoR9->at(index),event_weight);
}



std::vector<int> postAnalyzer::getQcdden(double phoPtCut, double phoEtaCut, Short_t isoBit){

  std::vector<int> tmpCand;
  tmpCand.clear();
    
  //Loop over photons
  for(int p=0;p<nPho;p++)
  {
    Float_t uncorrectedPhoEt = ((*phoEt)[p]);
    //Fail loose iso

    bool passChIsoLoose = TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < 8.61;
    bool passNeuIsoLoose = TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < (7.26 + (9.73714e-03* uncorrectedPhoEt) + (1.01605e-05* pow(uncorrectedPhoEt, 2.0)));
    bool passPhoIsoLoose = TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < (0.18 + (2.99336e-03* uncorrectedPhoEt));

    if(!passChIsoLoose || !passNeuIsoLoose || !passPhoIsoLoose)
      {
	bool kinematic = uncorrectedPhoEt > phoPtCut  && fabs((*phoSCEta)[p])<phoEtaCut;
	//"Very loose" ID cuts with inverted shape cut                                                                                                                                                                                                   
	bool photonID = (
			 ((*phoSigmaIEtaIEtaFull5x5)[p]  >  0.0103 ) &&
			 ((*phoHoverE)[p]                <  0.031   ) &&
			 ((*phohasPixelSeed)[p]              ==  0      ) &&
			 ( TMath::Max( ( (*phoPFChWorstIso)[p]  - rho*EAchargedworst((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(8.61) , 0.20*uncorrectedPhoEt) )  &&
			 ( TMath::Max( ( (*phoPFNeuIso)[p] - rho*EAneutral((*phoSCEta)[p]) ), 0.0) < TMath::Min(5.0*(7.26 + (9.73714e-03 * uncorrectedPhoEt) + (1.01605e-05* pow(uncorrectedPhoEt, 2.0))) , 0.20*uncorrectedPhoEt) )  &&
			 ( TMath::Max( ( (*phoPFPhoIso)[p] - rho*EAphoton((*phoSCEta)[p])  ), 0.0) < TMath::Min(5.0*(0.18 + (2.99336e-03 * uncorrectedPhoEt)) , 0.20*uncorrectedPhoEt) )

			 );

	bool noncoll = fabs((*phoSeedTime)[p]) < 3. && (*phoMIPTotEnergy)[p] < 4.9 && (*phoSigmaIEtaIEtaFull5x5)[p] > 0.001 && (*phoSigmaIPhiIPhiFull5x5)[p] > 0.001;
      if(kinematic && photonID && noncoll)
      {
        tmpCand.push_back(p);
      }
    }
  }

  return tmpCand;

}


// Effective area to be needed in PF Iso for photon ID
// https://indico.cern.ch/event/455258/contribution/0/attachments/1173322/1695132/SP15_253rd.pdf -- slide-5


Double_t postAnalyzer::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.06118;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.05545;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.05000;
  return EffectiveArea;
}

Double_t postAnalyzer::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.03279;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.05682;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.07969;
  return EffectiveArea;
}

Double_t postAnalyzer::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.1078;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.1068;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.08661;
  return EffectiveArea;
}


//Returns true if veto passed                                                                                                                                                                                                                                                                                                
//Veto failed if an electron is found that passes Loose Electron ID and elePtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                                 
//Always true if no electrons with |SC eta| < 2.5, since ID always fails for |SC eta| > 2.

bool postAnalyzer::electron_veto_looseID(int pho_index, float elePtCut)
{
  bool veto_passed = true; //pass veto if no good electron found                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      //Make sure these get reset for every electron                                                                                                                        
      /*                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;
      //Find EA for corrected relative iso.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      if(abs(eleSCEta->at(i)) <= 1.0)
  EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
  EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
  EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
  EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
  EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
  EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
  EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

      if(abs(eleSCEta->at(i)) <= 1.479)
  {
    pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
    pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
    pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
    pass_HoverE = eleHoverE->at(i) < 0.104;
    pass_iso = EAcorrIso < 0.0893;
    pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
    pass_d0 = abs(eleD0->at(i)) < 0.0261;
    pass_dz = abs(eleDz->at(i)) < 0.41;
    pass_missingHits = eleMissHits->at(i) <= 2;
    pass_convVeto = eleConvVeto->at(i) == 1;
  }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
  {
    pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
    pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
    pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
    pass_HoverE = eleHoverE->at(i) < 0.0897;
    pass_iso = EAcorrIso < 0.121;
    pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
    pass_d0 = abs(eleD0->at(i)) < 0.118;
    pass_dz = abs(eleDz->at(i)) < 0.822;
    pass_missingHits = eleMissHits->at(i) <= 1;
    pass_convVeto = eleConvVeto->at(i) == 1;
  }
      //Electron passes Loose Electron ID cuts                                                                                                                                                                                                                                                                     
      // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
      */
      
      if(eleIDbit->at(i)>>0&1==1)
  {
    //Electron passes pt cut                                                                                                                                                                                                                                                                             
    if(elePt->at(i) > elePtCut)
      {
        //Electron does not overlap photon                                                                                                                                                                                                                                                           
	//        if(dR(eleSCEta->at(i),eleSCPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
	  if(dR(eleEta->at(i),elePhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
    {
      veto_passed = false;
      break;
    }
      }
  }
    }
  return veto_passed;
}



//Veto failed if a muon is found that passes Loose Muon ID, Loose Muon Isolation, and muPtcut, and does not overlap the candidate photon within dR of 0.5                                                                                                                                                                    
bool postAnalyzer::muon_veto_looseID(int pho_index, float muPtCut)
{
  bool veto_passed = true; //pass veto if no good muon found                                                                                                                                                                                                                                                           
  bool pass_PFMuon = true;
  bool pass_globalMuon = true;
  bool pass_trackerMuon = true;
  bool pass_iso = false;
  //Explicitly stating types to avoid a TMath::Max conversion issue                                                                                                                                                                                                                                                    
  Float_t zero = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      // pass_PFMuon = muIsPFMuon->at(i);
      // pass_globalMuon = muIsGlobalMuon->at(i);
      // pass_trackerMuon = muIsTrackerMuon->at(i);
      /*      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      //Muon passes Loose Muon ID and PF-based combined relative, dBeta-corrected Loose Muon Isolation cuts                                                                  */                                                                                                                                      
      // if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon))
      //      if(muIDbit->at(i)>>0&1==1)
	if(muIDbit->at(i)>>0&1==1 && muIDbit->at(i)>>7&1==1)
  {
    //Muon passes pt cut                                                                                                                                                                                                                                                                                 
    if(muPt->at(i) > muPtCut)
      {
        //Muon does not overlap photon                                                                                                                                                                                                                                                               
	// if(dR(muEta->at(i),muPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index)) > 0.5)
	if(dR(muEta->at(i),muPhi->at(i),phoEta->at(pho_index),phoPhi->at(pho_index)) > 0.5)
	{
      veto_passed = false;
      break;
    }
      }
  }
    }
  return veto_passed;
}


std::vector<int> postAnalyzer::JetVetoDecision(int pho_index) {

  bool jetVeto=true;
  std::vector<int> jetindex;
  float value =0.0;

  for(int i = 0; i < nJet; i++)
    {

      if(0.0 < abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.5)
        value =-0.8;
      else if(2.5 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <2.75)
        value =-0.95;
      else if(2.75 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <3.0)
        value =-0.97;
      else if(3.00 <= abs(jetEta->at(i)) && abs(jetEta->at(i)) <5.0)
        value =-0.99;
      else
        continue;



      //      std::cout<<"Jet size: "<<nJet<<std::endl;
      double deltar = 0.0 ;
      //      std::cout<<"Jet no:"<<i<<"coming here pujetid: "<<pfJet_pt[i]<<std::endl;
      //if(OverlapWithMuon(jetEta->at(i),jetPhi->at(i)))     continue;
      //      std::cout<<"Jet no:"<<i<<"coming here OverlapWithMuon: "<<pfJet_pt[i]<<std::endl;
      //      if(OverlapWithElectron(jetEta->at(i),jetPhi->at(i)))   continue;
      if(pho_index>=0){
        deltar= dR(jetEta->at(i),jetPhi->at(i),phoSCEta->at(pho_index),phoSCPhi->at(pho_index));
        //std::cout<<"Delta R between photon and jet="<<dR<< "jet pt"<<pfJet_pt[i]<<std::endl;                                                                                       
      }
      //      if(deltar>0.4 && jetPt->at(i) >30.0 && jetPFLooseId->at(i)==1 && jetPUID->at(i)>value)
      if(deltar>0.4 && jetPt->at(i) >30.0 && (jetID->at(i) >>1&1==1))  
      {
          jetindex.push_back(i);
        }

      
    }


  //  std::cout<<"Jet size: "<< jetindex.size()<<std::endl;
  //if(jetindex.size()>1)jetVeto = false;                                                                                                                                            
  return jetindex;

}



bool postAnalyzer::OverlapWithMuon(double eta, double phi){
  
  bool overlap = false;
  //  std::cout<<"No of muon:"<<Muon_n<<std::endl;
  bool pass_PFMuon = false;
  bool pass_globalMuon = false;
  bool pass_trackerMuon = false;
  bool pass_iso = false;

  /*  Float_t zero1 = 0.0;
  Float_t muPhoPU = 999.9;
  Float_t tightIso_combinedRelative = 999.9;
  for(int i = 0; i < nMu; i++)
    {
      // pass_PFMuon = muIsPFMuon->at(i);
      // pass_globalMuon = muIsGlobalMuon->at(i);
      // pass_trackerMuon = muIsTrackerMuon->at(i);
      muPhoPU = muPFNeuIso->at(i) + muPFPhoIso->at(i) - 0.5*muPFPUIso->at(i);
      tightIso_combinedRelative = (muPFChIso->at(i) + TMath::Max(zero1,muPhoPU))/(muPt->at(i));
      pass_iso = tightIso_combinedRelative < 0.25;
      if(pass_PFMuon && (pass_globalMuon || pass_trackerMuon) && pass_iso)
        {
          if(muPt->at(i) > 10.)
            {
              if(dR(muEta->at(i),muPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }
  */
  return overlap;


}



bool postAnalyzer::OverlapWithElectron(double eta, double phi){
  bool overlap = false;

  bool pass_SigmaIEtaIEtaFull5x5 = false;
  bool pass_dEtaIn = false;
  bool pass_dPhiIn = false;
  bool pass_HoverE = false;
  bool pass_iso = false;
  bool pass_ooEmooP = false;
  bool pass_d0 = false;
  bool pass_dz = false;
  bool pass_missingHits = false;
  bool pass_convVeto = false;

  /*  Float_t EA = 0.0;
  Float_t zero = 0.0;
  Float_t EAcorrIso = 999.9;
  for(int i = 0; i < nEle; i++)
    {
      pass_SigmaIEtaIEtaFull5x5 = false;
      pass_dEtaIn = false;
      pass_dPhiIn = false;
      pass_HoverE = false;
      pass_iso = false;
      pass_ooEmooP = false;
      pass_d0 = false;
      pass_dz = false;
      pass_missingHits = false;
      pass_convVeto = false;

      if(abs(eleSCEta->at(i)) <= 1.0)
        EA = 0.1752;
      else if(1.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 1.479)
        EA = 0.1862;
      else if(1.479 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.0)
        EA = 0.1411;
      else if(2.0 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.2)
        EA = 0.1534;
      else if(2.2 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.3)
        EA = 0.1903;
      else if(2.3 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) <= 2.4)
        EA = 0.2243;
      else if(2.4 < abs(eleSCEta->at(i)) && abs(eleSCEta->at(i)) < 2.5)
        EA = 0.2687;
      EAcorrIso = (elePFChIso->at(i) + TMath::Max(zero,elePFNeuIso->at(i) + elePFPhoIso->at(i) - rho*EA))/(elePt->at(i));

      if(abs(eleSCEta->at(i)) <= 1.479)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0103;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.0105;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.115;
          pass_HoverE = eleHoverE->at(i) < 0.104;
          pass_iso = EAcorrIso < 0.0893;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.102;
          pass_d0 = abs(eleD0->at(i)) < 0.0261;
          pass_dz = abs(eleDz->at(i)) < 0.41;
          pass_missingHits = eleMissHits->at(i) <= 2;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      else if(1.479 < abs(eleSCEta->at(i)) < 2.5)
        {
          pass_SigmaIEtaIEtaFull5x5 = eleSigmaIEtaIEtaFull5x5->at(i) < 0.0301;
          pass_dEtaIn = abs(eledEtaAtVtx->at(i)) < 0.00814;
          pass_dPhiIn = abs(eledPhiAtVtx->at(i)) < 0.182;
          pass_HoverE = eleHoverE->at(i) < 0.0897;
          pass_iso = EAcorrIso < 0.121;
          pass_ooEmooP = eleEoverPInv->at(i) < 0.126;
          pass_d0 = abs(eleD0->at(i)) < 0.118;
          pass_dz = abs(eleDz->at(i)) < 0.822;
          pass_missingHits = eleMissHits->at(i) <= 1;
          pass_convVeto = eleConvVeto->at(i) == 1;
        }
      //Electron passes Loose Electron ID cuts                                                                                                                                                                                                                                                                              
      // if(pass_SigmaIEtaIEtaFull5x5 && pass_dEtaIn && pass_dPhiIn && pass_HoverE && pass_iso && pass_ooEmooP && pass_d0 && pass_dz && pass_missingHits && pass_convVeto)
      if(eleIDbit->at(i)>>1&1==1)
        {
          //Electron passes pt cut                                                                                                                                                                                                                                                                                          
          if(elePt->at(i) > 10.)
            {
              //Electron does not overlap photon                                                                                                                                                                                                                                                                            
              if(dR(eleSCEta->at(i),eleSCPhi->at(i),eta,phi) < 0.5)
                {
                  overlap = true;
                  break;
                }
            }
        }
    }

  */
  return overlap;



}


double postAnalyzer::dR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaeta = abs(eta1 - eta2);
  double deltaphi = DeltaPhi(phi1, phi2);
  double deltar = sqrt(deltaeta*deltaeta + deltaphi*deltaphi);
  return deltar;
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
double postAnalyzer::DeltaPhi(double phi1, double phi2)
{
  double pi = TMath::Pi();
  double dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}



bool postAnalyzer::dPhiJetMET_veto(std::vector<int> jets)
{
  //Pass veto if no jet is found within DeltaPhi(jet,MET) < 0.5                                                                                         
  bool passes = false;
  int njetsMax = jets.size();
  //Only look at first four jets                                                                                                                        
  if(njetsMax > 4)
    njetsMax = 4;
  int j = 0;
  for(; j < njetsMax; j++)
    {
      //fail veto if a jet is found with DeltaPhi(jet,MET) < 0.5                                                                                    
      if(DeltaPhi(jetPhi->at(jets[j]),pfMETPhi) < 0.5)
  break;
    }
  if(j==njetsMax)
    passes = true;

  return passes;
}

double postAnalyzer::FakeRatePt(Float_t phoPt, int syst_number){

  Float_t bin1[11] = {0.163464,0.160769,0.166192,0.165585,0.165818,0.163856,0.161117,0.163464,0.163464,0.166652,0.160276};
  Float_t bin2[11] = {0.0589625,0.0584719,0.0600986,0.0597755,0.060472,0.0592857,0.0586997,0.0589625,0.0589625,0.0600032,0.0579219};
  Float_t bin3[11] = {0.030392,0.0301976,0.0309801,0.0304578,0.0301257,0.03048,0.0299161,0.030392,0.030392,0.0310653,0.0297186};
  Float_t bin4[11] = {0.0119666,0.0123211,0.0129265,0.0122159,0.011731,0.0119247,0.0118692,0.0119666,0.0119666,0.0123333,0.0115998};
  Float_t bin5[11] = {0.0163508,0.0160008,0.0140943,0.016862,0.0122593,0.0147936,0.0181351,0.0163508,0.0163508,0.0188182,0.0138833};

  double weight = 1.0;
  if(phoPt <= 225.0)
    weight = 0;
  else if(225.0 < phoPt && phoPt <= 250.0)    weight = bin1[syst_number];
  else if(250.0 < phoPt && phoPt <= 300.0)    weight = bin2[syst_number];
  else if(300.0 < phoPt && phoPt <= 400.0)    weight = bin3[syst_number];
  else if(400.0 < phoPt && phoPt <= 600.0)    weight = bin4[syst_number];
  else if(600.0 < phoPt && phoPt <= 1000.0 )    weight = bin5[syst_number];
    //  else if(600.0 < phoPt <=1000.0 )                                                                                                                                   

  return weight;

}

bool postAnalyzer::getMetFilter(){
  bool decision = true;
  if(metFilters>>0 & 1) decision = false;
  if(metFilters>>1 & 1) decision = false;
  if(metFilters>>2 & 1) decision = false;
  if(metFilters>>3 & 1) decision = false;
  if(metFilters>>4 & 1) decision = false;
  if(metFilters>>5 & 1) decision = false;
  if(metFilters>>7 & 1) decision = false;
  if(metFilters>>8 & 1) decision = false;
  return decision;
}
