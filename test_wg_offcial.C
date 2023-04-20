#define postAnalyzer_cxx
#include "test_wg_offcial.h"
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
  int nInspected;
  nInspected = 0;
  double nInspected_genWeighted;
  nInspected_genWeighted = 0.0;
  int nPhoCand, nMET170, nDphiPhoMET, nLeptonVeto, nDphiJetMET;
  nPhoCand = nMET170 = nDphiPhoMET = nLeptonVeto = nDphiJetMET = 0;
  double nPhoCand_weighted, nMET170_weighted, nDphiPhoMET_weighted, nLeptonVeto_weighted, nDphiJetMET_weighted;
  nPhoCand_weighted = nMET170_weighted = nDphiPhoMET_weighted = nLeptonVeto_weighted = nDphiJetMET_weighted = 0.0;

  std::vector<int> phoCand;
  phoCand.clear();
  std::vector<int> phoCandUp;
  phoCandUp.clear();
  std::vector<int> phoCandDown;
  phoCandDown.clear();

   std::vector<int> qcdden;
   qcdden.clear();
   
   std::vector<int> jetveto;
   jetveto.clear();
   
   /*  TFile *file = new TFile("ewk_corr.root");
  TH1D *ewkCorrection = (TH1D*)file->Get("wnlg-130-o_p");
  TH1D *ewkCorrection_straightUp = (TH1D*)file->Get("wnlg-130-o_p_straightUp");
  TH1D *ewkCorrection_straightDown = (TH1D*)file->Get("wnlg-130-o_p_straightDown");
  TH1D *ewkCorrection_twistedUp = (TH1D*)file->Get("wnlg-130-o_p_twistedUp");
  TH1D *ewkCorrection_twistedDown = (TH1D*)file->Get("wnlg-130-o_p_twistedDown");
  TH1D *ewkCorrection_gammaUp = (TH1D*)file->Get("wnlg-130-o_p_gammaUp");
  TH1D *ewkCorrection_gammaDown = (TH1D*)file->Get("wnlg-130-o_p_gammaDown");
  ewkCorrection->Scale(0.62);
  ewkCorrection_straightUp->Scale(0.62);
  ewkCorrection_straightDown->Scale(0.62);
  ewkCorrection_twistedUp->Scale(0.62);
  ewkCorrection_twistedDown->Scale(0.62);
  ewkCorrection_gammaUp->Scale(0.62);
  ewkCorrection_gammaDown->Scale(0.62);
  TH1D *ewkCorrection_m = (TH1D*)file->Get("wnlg-130-o_m");
  TH1D *ewkCorrection_straightUp_m = (TH1D*)file->Get("wnlg-130-o_m_straightUp");
  TH1D *ewkCorrection_straightDown_m = (TH1D*)file->Get("wnlg-130-o_m_straightDown");
  TH1D *ewkCorrection_twistedUp_m = (TH1D*)file->Get("wnlg-130-o_m_twistedUp");
  TH1D *ewkCorrection_twistedDown_m = (TH1D*)file->Get("wnlg-130-o_m_twistedDown");
  TH1D *ewkCorrection_gammaUp_m = (TH1D*)file->Get("wnlg-130-o_m_gammaUp");
  TH1D *ewkCorrection_gammaDown_m = (TH1D*)file->Get("wnlg-130-o_m_gammaDown");
  ewkCorrection->Add(ewkCorrection_m, 0.38);
  ewkCorrection_straightUp->Add(ewkCorrection_straightUp_m, 0.38);
  ewkCorrection_straightDown->Add(ewkCorrection_straightDown_m, 0.38);
  ewkCorrection_twistedUp->Add(ewkCorrection_twistedUp_m, 0.38);
  ewkCorrection_twistedDown->Add(ewkCorrection_twistedDown_m, 0.38);
  ewkCorrection_gammaUp->Add(ewkCorrection_gammaUp_m, 0.38);
  ewkCorrection_gammaDown->Add(ewkCorrection_gammaDown_m, 0.38);

  //Extracting Scalefactors value from the root file                                                                                                                     
  TFile *pfroot = new TFile("egammaEffi.txt_EGM2D_uncer_mySF.root");
  TH2D *phcorrect = (TH2D*)pfroot->Get("EGamma_SF2D");
  TH2D *phcorrect_sys = (TH2D*)pfroot->Get("h2_uncertaintiesEGamma_copy");
   */
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
    
    double inspected_event_weight = 1.0;
    fabs(genWeight) > 0.0 ? inspected_event_weight *= genWeight/fabs(genWeight) : inspected_event_weight = 0.0; //Generator may have given event negative weight
    nInspected_genWeighted += inspected_event_weight;
    nInspected += 1;
    //=1.0 for real data
    //     int ewkCorrectionBin = min(max(ewkCorrection->GetXaxis()->FindBin(LHephoET), 1), ewkCorrection->GetXaxis()->GetNbins());
    double event_weight=1.0;
    double event_weight_straightUp=1.0;
    double event_weight_straightDown=1.0;
    double event_weight_twistedUp=1.0;
    double event_weight_twistedDown=1.0;
    double event_weight_gammaUp=1.0;
    double event_weight_gammaDown=1.0;
    double uncorrected_weight=1.0;
    cout<< nPho << "  " <<  nEle <<  " " << nMu << " " << nJet <<  " " << pfMET << endl;
    
    TH1F *h_nPho = new TH1F("h_nPho", "", 50, -0.5, 49.5);                                                                                                             
    //  TH1F *h_nEle = new TH1F("h_nEle", "", 50, -0.5, 49.5);                                                                                                             
    //  TH1F *h_nMuo = new TH1F("h_nMuo", "", 50, -0.5, 49.5);                                                                                                             
    //  TH1F *h_nJet = new TH1F("h_nJet", "", 50, -0.5, 49.5);                                                                                                             
    //  TH1F *h_MET  = new TH1F("h_MET",  "", 100, -0.5, 999.5);
   
    h_nPho->Fill(nPho);                                                                                                                                                
    //    h_nEle->Fill(nEle);                                                                                                                                            
    //    h_nMuo->Fill(nMu);                                                                                                                                             
    //    h_nJet->Fill(nJet);                                                                                                                                            
    //    h_MET->Fill(pfMET);


    //    phoCand   = getPhoCand(225,1.4442,1);
    //phoCandUp   = getPhoCandUp(225,1.4442,1);
    // phoCandDown   = getPhoCandDown(225,1.4442,1);

    //    if(true || metFilters==1536)
    /*      if (!getMetFilter()) continue;
    {
    //  nMETFiltersPassed++;
    //  if((HLTPho>>7&1 == 1)||(HLTPho>>8&1 == 1)|| (HLTPho>>9&1 == 1) || (HLTPho>>10&1 == 1) || (HLTPho>>11&1 == 1) || (HLTPho>>12&1 == 1) || (HLTPho>>22&1 == 1))
    //  {
    //    nHLTPassed++;
    if(phoCand.size() >0)
    {
      // if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCand[0]] - rho*EAchargedworst((*phoSCEta)[phoCand[0]]) ), 0.0) < 1.37 )
      // {
      //        Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCand[0]]/TMath::CosH((*phoSCEta)[phoCand[0]]));
	Float_t uncorrectedPhoEt = ((*phoEt)[phoCand[0]]);
	event_weight = 1.0;
        Double_t phweight = 1.0;
	Double_t phsys=0.0;
        if( uncorrectedPhoEt >= 1000  || fabs(phoEta->at(phoCand[0]))>1.4442){phweight = 1;}
        else{
          phweight = phcorrect->GetBinContent(phcorrect->GetXaxis()->FindBin(fabs(phoEta->at(phoCand[0]))),phcorrect->GetYaxis()->FindBin(phoEt->at(phoCand[0])));
	  phsys = phcorrect_sys->GetBinContent(phcorrect_sys->GetXaxis()->FindBin(fabs(phoEta->at(phoCand[0]))),phcorrect_sys->GetYaxis()->FindBin(phoEt->at(phoCand[0])));
        }
	Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	event_weight *= EWK_percent_adjustment;
	event_weight *= NNLOCorrection(LHephoET);
        event_weight*=(1.002 - 0.00004395*phoEt->at(phoCand[0])); //Trigger inefficiency correction
	event_weight=event_weight*phweight;
	fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
       

        event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
        event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
        event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
        event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
        event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
        event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
	double systematic_up;
	Float_t systematic_down;
	double event_weight_phoscaleUp;
	double event_weight_phoscaledown;

	systematic_up = phweight + phsys;
	systematic_down = phweight - phsys;
	//              cout<<"systematic_up"<<systematic_up<<endl;                                                                                                        
	//cout<<"systematic_down"<<systematic_down<<endl;                                                                                                                  
	event_weight_phoscaleUp = event_weight*systematic_up;
	event_weight_phoscaledown = event_weight*systematic_down;
	//my changes closed                                                                                                                                                
	cout<<"event_weight_phoscaleUp"<<event_weight_phoscaleUp<<endl;
	cout<<"event_weight_phoscaledown"<<event_weight_phoscaledown<<endl;
	
	uncorrected_weight *= (1.002 - 0.00004395*phoEt->at(phoCand[0]));
        fabs(genWeight) > 0.0 ? uncorrected_weight *= genWeight/fabs(genWeight) : uncorrected_weight = 0;

	nPhoCand++;
        nPhoCand_weighted += event_weight;
        
        jetveto = JetVetoDecision(phoCand[0]);
                
        if(pfMET>200)
        {
          nMET170++;
          nMET170_weighted += event_weight;
          
          Float_t MET_to_use = pfMET;
          Float_t METPhi_to_use = pfMETPhi;
          
          fillHistos(0,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
          if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
          {
            nDphiPhoMET++;
            nDphiPhoMET_weighted += event_weight;
            fillHistos(1,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
            if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
            {
              nLeptonVeto++;
              nLeptonVeto_weighted += event_weight;
              fillHistos(2,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
              if(dPhiJetMET_veto(jetveto,METPhi_to_use))
              {
                nDphiJetMET++;
                nDphiJetMET_weighted += event_weight;
                fillHistos(3,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                fillHistos(30,event_weight*NNLOCorrection_err(LHephoET),phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                if(uncorrectedPhoEt/MET_to_use < 1.4){
                  fillHistos(24,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(21,event_weight_straightUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(22,event_weight_straightDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(31,event_weight_twistedUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(32,event_weight_twistedDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(33,event_weight_gammaUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(34,event_weight_gammaDown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(20,uncorrected_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                  fillHistos(23,event_weight*NNLOCorrection_err(LHephoET),phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  fillHistos(35,event_weight_phoscaleUp,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		  fillHistos(36,event_weight_phoscaledown,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
		}
                
                //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
              }
            }
          }
        }
        if(pfMET_T1JESUp > 200)
        {
          Float_t MET_to_use = pfMET_T1JESUp;
          Float_t METPhi_to_use = pfMETPhi_T1JESUp;
          
          fillHistos(4,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
          if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
          {
            fillHistos(5,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
            if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
            {
              fillHistos(6,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
              if(dPhiJetMET_veto(jetveto,METPhi_to_use))
              {
                fillHistos(7,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                if(uncorrectedPhoEt/MET_to_use < 1.4){
                  fillHistos(25,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                }
                //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
              }
            }
          }
        }
        if(pfMET_T1JESDo > 200)
        {
          Float_t MET_to_use = pfMET_T1JESDo;
          Float_t METPhi_to_use = pfMETPhi_T1JESDo;
          
          fillHistos(8,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
          if(DeltaPhi(phoPhi->at(phoCand[0]),METPhi_to_use)>0.5)
          {
            fillHistos(9,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
            if(electron_veto_looseID(phoCand[0],10) && muon_veto_looseID(phoCand[0],10))
            {
              fillHistos(10,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
              if(dPhiJetMET_veto(jetveto,METPhi_to_use))
              {
                fillHistos(11,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                if(uncorrectedPhoEt/MET_to_use < 1.4){
                  fillHistos(26,event_weight,phoCand[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                }
                //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
              }
            }
          }
        }
      // }
    }
    if(phoCandUp.size() > 0)
    {
      // if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCandUp[0]] - rho*EAchargedworst((*phoSCEta)[phoCandUp[0]]) ), 0.0) < 1.37 )
      // {
        event_weight=1.0;
	//        Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCandUp[0]]/TMath::CosH((*phoSCEta)[phoCandUp[0]]));
	Float_t uncorrectedPhoEt = ((*phoEt)[phoCandUp[0]]);
	uncorrectedPhoEt += 0.015*uncorrectedPhoEt;
	Double_t phweight = 1.0;
        if( uncorrectedPhoEt >= 1000 || fabs(phoEta->at(phoCandUp[0]))>1.4442){phweight = 1;}
        else{
          phweight = phcorrect->GetBinContent(phcorrect->GetXaxis()->FindBin(fabs(phoEta->at(phoCandUp[0]))),phcorrect->GetYaxis()->FindBin(phoEt->at(phoCandUp[0])));
        }

	Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	event_weight *= EWK_percent_adjustment;
	event_weight *= NNLOCorrection(LHephoET);
        event_weight *= (1.002 - 0.00004395*phoEt->at(phoCandUp[0])); //Trigger inefficiency correction
	event_weight=event_weight*phweight;
	fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
        event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
        event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
        event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
        event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
        event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
        event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
        
        jetveto = JetVetoDecision(phoCandUp[0]);
        
                
        if(pfMET_T1PESUp>200)
        {
          Float_t MET_to_use = pfMET_T1PESUp;
          Float_t METPhi_to_use = pfMETPhi_T1PESUp;
          
          fillHistos(12,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
          if(DeltaPhi(phoPhi->at(phoCandUp[0]),METPhi_to_use)>0.5)
          {
            fillHistos(13,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
            if(electron_veto_looseID(phoCandUp[0],10) && muon_veto_looseID(phoCandUp[0],10))
            {
              fillHistos(14,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
              if(dPhiJetMET_veto(jetveto,METPhi_to_use))
              {
                fillHistos(15,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                if(uncorrectedPhoEt/MET_to_use < 1.4){
                  fillHistos(27,event_weight,phoCandUp[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                }
                //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
              }
            }
          }
        }
      // }
    }
    if(phoCandDown.size() > 0)
      {
	// if( TMath::Max( ( (*phoYuPFChWorstIso)[phoCandDown[0]] - rho*EAchargedworst((*phoSCEta)[phoCandDown[0]]) ), 0.0) < 1.37 )
	// {
        event_weight=1.0;
	//        Float_t uncorrectedPhoEt = ((*phoSCRawE)[phoCandDown[0]]/TMath::CosH((*phoSCEta)[phoCandDown[0]]));
	Float_t uncorrectedPhoEt = ((*phoEt)[phoCandDown[0]]);
	uncorrectedPhoEt -= 0.015*uncorrectedPhoEt;
	Double_t phweight = 1.0;
	if( uncorrectedPhoEt >= 1000 || fabs(phoEta->at(phoCandDown[0]))>1.4442){phweight = 1;}
	else{
	  phweight = phcorrect->GetBinContent(phcorrect->GetXaxis()->FindBin(fabs(phoEta->at(phoCandDown[0]))),phcorrect->GetYaxis()->FindBin(phoEt->at(phoCandDown[0])));
	}
	

	Double_t EWK_percent_adjustment = ewkCorrection->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightUp = ewkCorrection_straightUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_straightDown = ewkCorrection_straightDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedUp = ewkCorrection_twistedUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_twistedDown = ewkCorrection_twistedDown->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaUp = ewkCorrection_gammaUp->GetBinContent(ewkCorrectionBin);
	Double_t EWK_percent_adjustment_gammaDown = ewkCorrection_gammaDown->GetBinContent(ewkCorrectionBin);
	event_weight *= EWK_percent_adjustment;
	event_weight *= NNLOCorrection(LHephoET);
        event_weight *= (1.002 - 0.00004395*phoEt->at(phoCandDown[0])); //Trigger inefficiency correction
	event_weight=event_weight*phweight;
	fabs(genWeight) > 0.0 ? event_weight *= genWeight/fabs(genWeight) : event_weight = 0; //Generator may have given event negative weight
        event_weight_straightUp = event_weight*EWK_percent_adjustment_straightUp/EWK_percent_adjustment;
        event_weight_straightDown = event_weight*EWK_percent_adjustment_straightDown/EWK_percent_adjustment;
        event_weight_twistedUp = event_weight*EWK_percent_adjustment_twistedUp/EWK_percent_adjustment;
        event_weight_twistedDown = event_weight*EWK_percent_adjustment_twistedDown/EWK_percent_adjustment;
        event_weight_gammaUp = event_weight*EWK_percent_adjustment_gammaUp/EWK_percent_adjustment;
        event_weight_gammaDown = event_weight*EWK_percent_adjustment_gammaDown/EWK_percent_adjustment;
        
        jetveto = JetVetoDecision(phoCandDown[0]);
        
        
        if(pfMET_T1PESDo>200)
        {
          Float_t MET_to_use = pfMET_T1PESDo;
          Float_t METPhi_to_use = pfMETPhi_T1PESDo;
          
          fillHistos(16,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
          if(DeltaPhi(phoPhi->at(phoCandDown[0]),METPhi_to_use)>0.5)
          {
            fillHistos(17,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
            if(electron_veto_looseID(phoCandDown[0],10) && muon_veto_looseID(phoCandDown[0],10))
            {
              fillHistos(18,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
              if(dPhiJetMET_veto(jetveto,METPhi_to_use))
              {
                fillHistos(19,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                if(uncorrectedPhoEt/MET_to_use < 1.4){
                  fillHistos(28,event_weight,phoCandDown[0],jetveto,MET_to_use,uncorrectedPhoEt,METPhi_to_use);
                }
                //std::cout<<"run:lumi:event:"<<run<<":"<<lumis<<":"<<event<<std::endl;
              }
            }
          }
        }
      // }
    }
    }*/
    
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
  std::cout << "Number of events inspected: " << nInspected << std::endl;
  std::cout << "Number of events inspected (minus negative gen. weights): " << nInspected_genWeighted << std::endl;
  std::cout<<std::endl;
  /*  cout<<"Unweighted"<<endl;
  cout<<"----------"<<endl;
  cout<<"nPhoCand: "<<nPhoCand; if(nInspected>0.){cout<<", "<<nPhoCand/nInspected<<" of previous";} cout<<endl;
  cout<<"nMET170: "<<nMET170; if(nPhoCand>0.){cout<<", "<<nMET170/nPhoCand<<" of previous";} cout<<endl;
  cout<<"nDphiPhoMET: "<<nDphiPhoMET; if(nMET170>0.){cout<<", "<<nDphiPhoMET/nMET170<<" of previous";} cout<<endl;
  cout<<"nLeptonVeto: "<<nLeptonVeto; if(nDphiPhoMET>0.){cout<<", "<<nLeptonVeto/nDphiPhoMET<<" of previous";} cout<<endl;
  cout<<"nDphiJetMET: "<<nDphiJetMET; if(nLeptonVeto>0.){cout<<", "<<nDphiJetMET/nLeptonVeto<<" of previous";} cout<<endl;
  cout<<endl;
  cout<<"Weighted"<<endl;
  cout<<"--------"<<endl;
  cout<<"nPhoCand_weighted: "<<nPhoCand_weighted; if(nInspected_genWeighted>0.){cout<<", "<<nPhoCand_weighted/nInspected_genWeighted<<" of previous";} cout<<endl;
  cout<<"nMET170_weighted: "<<nMET170_weighted; if(nPhoCand_weighted>0.){cout<<", "<<nMET170_weighted/nPhoCand_weighted<<" of previous";} cout<<endl;
  cout<<"nDphiPhoMET_weighted: "<<nDphiPhoMET_weighted; if(nMET170_weighted>0.){cout<<", "<<nDphiPhoMET_weighted/nMET170_weighted<<" of previous";} cout<<endl;
  cout<<"nLeptonVeto_weighted: "<<nLeptonVeto_weighted; if(nDphiPhoMET_weighted>0.){cout<<", "<<nLeptonVeto_weighted/nDphiPhoMET_weighted<<" of previous";} cout<<endl;
  cout<<"nDphiJetMET_weighted: "<<nDphiJetMET_weighted; if(nLeptonVeto_weighted>0.){cout<<", "<<nDphiJetMET_weighted/nLeptonVeto_weighted<<" of previous";} cout<<endl;*/
}

void postAnalyzer::BookHistos(const char* file2)
{
  
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();

  fileName->Close();
  /*  tree = new TTree("ADD","ADD");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);

  
 //h_phoIEtaIPhi = new TH2F("h_phoIEtaIPhi","Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi->Sumw2();
  Float_t PtBins[6]={225.,250., 300., 400., 600., 1000.0};
  Float_t MetBins[6]={200.,250., 300., 400., 600., 1000.0};
  Float_t dPhiJetMETBins[14]={0.0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,2.50,2.75,3.00,3.142};
  //Set up the histos to be filled with method fillHistos
  for(int i=0; i<40; i++)
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
    // h_phoIEtaIPhi[i] = new TH2F(("h_phoIEtaIPhi"+histname).c_str(),"Photon p_{T} > 175 GeV, E^{miss}_{T} > 140 GeV",360,0.5,360.5,186,-85.5,100.5);h_phoIEtaIPhi[i]->Sumw2();
    h_PTMET[i] = new TH1F(("PTMET"+histname).c_str(),"P_{T}/Missing E_{T}",50,0,3);h_PTMET[i]->Sumw2();
    //    h_Mt[i]= new TH1F(("Mt"+histname).c_str(),"MT",9,MTBins);h_Mt[i]->Sumw2();
    h_min_dphijetmet[i] = new TH1F(("h_min_dphijetmet"+histname).c_str(),"h_min_dphijetmet",13,dPhiJetMETBins);h_min_dphijetmet[i]->Sumw2();
    h_pfMETsumEt[i] = new TH1F(("pfMETsumEt"+histname).c_str(),"pfMETsumEt",5,MetBins);
    h_METoverSqrtSumEt_extended[i] = new TH1F(("METoverSqrtSumEt_extended"+histname).c_str(),"METoverSqrtSumEt",30,0,30);
    h_METoverSqrtSumEt[i] = new TH1F(("METoverSqrtSumEt"+histname).c_str(),"METoverSqrtSumEt",30,0,10);
    h_r9[i] = new TH1F(("r9"+histname).c_str(),"r9",60,0,1);
  }
  */}

//Fill the sequential histos at a particular spot in the sequence
