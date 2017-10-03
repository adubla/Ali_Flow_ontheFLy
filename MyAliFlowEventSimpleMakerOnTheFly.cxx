/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/************************************
 * Create an event and perform full *
 * flow analysis 'on the fly'.      *
 *                                  *
 * author: Ante Bilandzic           *
 *         (abilandzic@gmail.com)   *
 ************************************/

#include "Riostream.h"
#include "TMath.h"
#include "TF1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "MyAliFlowEventSimpleMakerOnTheFly.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowVector.h"

class TH2D;

using std::endl;
using std::cout;
ClassImp(MyAliFlowEventSimpleMakerOnTheFly)

//========================================================================================================================================

MyAliFlowEventSimpleMakerOnTheFly::MyAliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed):
fCount(0),
fMinMult(0),
fMaxMult(0),
fPtSpectra(NULL),
fMass(0.13957),
fTemperature(0.44),
fPhiDistribution(NULL),
fV1(0.),
fV2(0.05),
fV3(0.),
fV4(0.),
fV5(0.),
fV6(0.),
fUniformFluctuationsV1(kFALSE),
fMinV1(0.04),
fMaxV1(0.06),
fUniformFluctuationsV2(kFALSE),
fMinV2(0.04),
fMaxV2(0.06),
fGausFluctuationsV1(kFALSE),
fSigV1(0.01),
fMeanV1(0.),
fPtDependentV2(kFALSE),
fV2vsPtCutOff(2.0),
fV2vsPtMax(0.2),
fEtaMinA(-0.8),
fEtaMaxA(-0.5),
fEtaMinB(0.5),
fEtaMaxB(0.8),
fNTimes(1),
fUniformAcceptance(kTRUE),
fPhiMin1(0.),
fPhiMax1(0.),
fProbability1(0.),
fPhiMin2(0.),
fPhiMax2(0.),
fProbability2(0.),
fPi(TMath::Pi()),
fUniformEfficiency(kTRUE),
fNUEPtMin(0.5),
fNUEPtMax(1.0),
fNUEPtProbability(0.75),
fV1B(0.),
fV3B(0.),
fV1EbE(0.),
fCME(0.),
fSimCRC(kFALSE),
fSimCME(kFALSE),
fSimEtaDep(kFALSE),
fSimPtDep(kFALSE),
fUseVZERO(kFALSE),
fChargeDepNUA(kFALSE),
fFuckEP(kFALSE)
{
  // Constructor.

  // Determine seed for gRandom:
  delete gRandom;
  gRandom = new TRandom3(uiSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

} // end of MyAliFlowEventSimpleMakerOnTheFly::MyAliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed):

//====================================================================================================================

MyAliFlowEventSimpleMakerOnTheFly::~MyAliFlowEventSimpleMakerOnTheFly()
{
  // Destructor.

  if(fPtSpectra){delete fPtSpectra;}
  if(fPhiDistribution){delete fPhiDistribution;}

} // end of MyAliFlowEventSimpleMakerOnTheFly::~MyAliFlowEventSimpleMakerOnTheFly()

//====================================================================================================================

void MyAliFlowEventSimpleMakerOnTheFly::Init()
{
  // Book all objects in this method.

  // a) Define the pt spectra;
  // b) Define the phi distribution.

  // a) Define the pt spectra:
  fPtMin = 0.2;
  fPtMax = 8.0;
  fPtSpectra = new TF1("fPtSpectra","x*TMath::Exp(-pow([0]*[0]+x*x,0.5)/[1])",fPtMin,fPtMax); // hardwired is Boltzmann distribution
  fPtSpectra->SetParName(0,"Mass");
  fPtSpectra->SetParameter(0,fMass);
  fPtSpectra->SetParName(1,"Temperature");
  fPtSpectra->SetParameter(1,fTemperature);
  fPtSpectra->SetTitle("Boltzmann Distribution: f(p_{t}) = p_{t}exp[-(m^{2}+p_{t}^{2})^{1/2}/T];p_{t};f(p_{t})");

  // b) Define the phi distribution:
  Double_t dPhiMin = 0.;
  Double_t dPhiMax = TMath::TwoPi();

  // fPhiDistribution = new TF1("fPhiDistribution","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[3]*TMath::Cos(3.*(x-[0]))+2.*[4]*TMath::Cos(4.*(x-[0]))+2.*[5]*TMath::Cos(5.*(x-[0]))+2.*[6]*TMath::Cos(6.*(x-[0]))",dPhiMin,dPhiMax);
  // fPhiDistribution->SetParName(0,"Reaction Plane");
  // fPhiDistribution->SetParameter(0,0.);
  // fPhiDistribution->SetParName(1,"Directed Flow (v1)");
  // fPhiDistribution->SetParameter(1,fV1);
  // fPhiDistribution->SetParName(2,"Elliptic Flow (v2)");
  // fPhiDistribution->SetParameter(2,fV2);
  // fPhiDistribution->SetParName(3,"Triangular Flow (v3)");
  // fPhiDistribution->SetParameter(3,fV3);
  // fPhiDistribution->SetParName(4,"Quadrangular Flow (v4)");
  // fPhiDistribution->SetParameter(4,fV4);
  // fPhiDistribution->SetParName(5,"Pentagonal Flow (v5)");
  // fPhiDistribution->SetParameter(5,fV5);
  // fPhiDistribution->SetParName(6,"Hexagonal Flow (v6)");
  // fPhiDistribution->SetParameter(6,fV6);

  fPhiDistribution = new TF1("fPhiDistribution","1+2.*[1]*TMath::Cos(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[3]*TMath::Cos(3.*(x-[0]))+2.*[4]*TMath::Cos(4.*(x-[0]))",dPhiMin,dPhiMax);
  fPhiDistribution->SetParName(0,"Reaction Plane");
  fPhiDistribution->SetParameter(0,0.);
  fPhiDistribution->SetParName(1,"Directed Flow (v1)");
  fPhiDistribution->SetParameter(1,fV1);
  fPhiDistribution->SetParName(2,"Elliptic Flow (v2)");
  fPhiDistribution->SetParameter(2,fV2);
  fPhiDistribution->SetParName(3,"Triangular Flow (v3)");
  fPhiDistribution->SetParameter(3,fV3);
  fPhiDistribution->SetParName(4,"Quadrangular Flow (v4)");
  fPhiDistribution->SetParameter(4,fV4);

  fPhiDistributionP = new TF1("fPhiDistributionP","1+2.*[1]*TMath::Cos(x-[0])+2.*[3]*TMath::Sin(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[4]*TMath::Cos(3.*(x-[0]))",dPhiMin,dPhiMax);
  fPhiDistributionP->SetParName(0,"Reaction Plane");
  fPhiDistributionP->SetParameter(0,0.);
  fPhiDistributionP->SetParName(1,"Directed Flow (v1)");
  fPhiDistributionP->SetParameter(1,fV1);
  fPhiDistributionP->SetParName(2,"Elliptic Flow (v2)");
  fPhiDistributionP->SetParameter(2,fV2);
  fPhiDistributionP->SetParName(3,"CME");
  fPhiDistributionP->SetParameter(3,fCME);
  fPhiDistributionP->SetParName(4,"Triangular Flow (v3)");
  fPhiDistributionP->SetParameter(4,fV3);
  fPhiDistributionN = new TF1("fPhiDistributionN","1+2.*[1]*TMath::Cos(x-[0])+2.*[3]*TMath::Sin(x-[0])+2.*[2]*TMath::Cos(2.*(x-[0]))+2.*[4]*TMath::Cos(3.*(x-[0]))",dPhiMin,dPhiMax);
  fPhiDistributionN->SetParName(0,"Reaction Plane");
  fPhiDistributionN->SetParameter(0,0.);
  fPhiDistributionN->SetParName(1,"Directed Flow (v1)");
  fPhiDistributionN->SetParameter(1,fV1);
  fPhiDistributionN->SetParName(2,"Elliptic Flow (v2)");
  fPhiDistributionN->SetParameter(2,fV2);
  fPhiDistributionN->SetParName(3,"CME");
  fPhiDistributionN->SetParameter(3,fCME);
  fPhiDistributionN->SetParName(4,"Triangular Flow (v3)");
  fPhiDistributionN->SetParameter(4,fV3);

  // c) Define eta-dependent phi distr.
  fEtaMin = -10.;
  fEtaMax = 10.;
  fEtaRange = fEtaMax-fEtaMin;
  fEtaBinWidth = 1.6/fnEtaBins;

  for (Int_t e=0; e<fnEtaBins; e++) {
    fPhiDistrEtaP[e] = new TF1(*fPhiDistributionP);
    fPhiDistrEtaN[e] = new TF1(*fPhiDistributionN);
  }

  // d) Define pt-dependent phi distr.
  fPtRange = fPtMax-fPtMin;
  fPtBinWidth = fPtRange/fnPtBins;

  for (Int_t p=0; p<fnPtBins; p++) {
    fPhiDistrPtP[p] = new TF1(*fPhiDistribution);
    fPhiDistrPtN[p] = new TF1(*fPhiDistribution);
  }

  fPsiDistribution = new TF1("fPsiDistribution","1+2.*[1]*TMath::Cos(2.*(x-[0]))",dPhiMin,dPhiMax);
  fPsiDistribution->SetParName(0,"Reaction Plane");
  fPsiDistribution->SetParameter(0,0.);
  fPsiDistribution->SetParName(1,"Elliptic Flow (v2)");
  fPsiDistribution->SetParameter(1,fV2);

  for (Int_t p=0; p<fnvonMis; p++) {
    fvonMisesDistrib[p] = new TF1("fvonMisesDistrib","TMath::Exp([0]*TMath::Cos(x))/(TMath::TwoPi()*TMath::BesselI0([0]))",-TMath::Pi(),TMath::Pi());
    fvonMisesDistrib[p]->SetParName(0,"Concentration");
    fvonMisesDistrib[p]->SetParameter(0,10.-0.09*p);
  }

  fPhiDistrEtaVA = new TF1(*fPhiDistributionP);
  fPhiDistrEtaVC = new TF1(*fPhiDistributionP);

  Double_t fV1Odd = -0.2;//0.2;
  Double_t fV1Eve = 0.;//0.1;

  for (Int_t e=0; e<fnEtaBins; e++) {
    Double_t V1Be = fV1B*GetEtafromBin(e);//(GetEtafromBin(e)<0.?fV1B:-fV1B);
    Double_t V1Odde = fV1*GetEtafromBin(e);
    Double_t V1Evee = 0.;//fV1/2.;
    Double_t V1ChInd = V1Evee+V1Odde;
    fPhiDistrEtaP[e]->SetParameter(1,V1ChInd+V1Be);
    fPhiDistrEtaN[e]->SetParameter(1,V1ChInd-V1Be);
  }
  fPhiDistrEtaVA->SetParameter(1,2.*fV1Odd);
  fPhiDistrEtaVC->SetParameter(1,-2.*fV1Odd);

  fZDCMult = new TProfile("fZDCMult","fZDCMult",8,0,8,"s");

  BalanFunc[0] = new TH1D("BalanFunc0","BalanFunc0",100,-TMath::Pi(),TMath::Pi());
  BalanFunc[1] = new TH1D("BalanFunc1","BalanFunc1",100,-TMath::Pi(),TMath::Pi());

  RefEtaHist = new TH1D("RefEtaHist","RefEtaHist",fnEtaBins,-0.8,0.8);

} // end of void MyAliFlowEventSimpleMakerOnTheFly::Init()

//====================================================================================================================

Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPhi(AliFlowTrackSimple *pTrack)
{
  // For the case of non-uniform acceptance determine in this method if particle is accepted or rejected for a given phi.
  Bool_t bAccept = kTRUE;
  if((pTrack->Phi() >= fPhiMin1*fPi/180.) && (pTrack->Phi() < fPhiMax1*fPi/180.) && gRandom->Uniform(0,1) > fProbability1) {
    bAccept = kFALSE; // particle is rejected in the first non-uniform sector
  } else if((pTrack->Phi() >= fPhiMin2*fPi/180.) && (pTrack->Phi() < fPhiMax2*fPi/180.) && gRandom->Uniform(0,1) > fProbability2) {
    bAccept = kFALSE; // particle is rejected in the second non-uniform sector
  }
  return bAccept;
} // end of Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPhi(AliFlowTrackSimple *pTrack);

//====================================================================================================================

Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPhiCharge(AliFlowTrackSimple *pTrack, Int_t dCharge, Double_t dEta)
{
  // For the case of non-uniform acceptance determine in this method if particle is accepted or rejected for a given phi.
  Bool_t bAccept = kTRUE;
  if(dCharge > 0. && dEta > 0.) {
    if((pTrack->Phi() >= fPhiMin1*fPi/180.) && (pTrack->Phi() < fPhiMax1*fPi/180.) && gRandom->Uniform(0,1) > fProbability1) {
      bAccept = kFALSE; // particle is rejected in the first non-uniform sector
    }
  }
  if(dCharge < 0. && dEta < 0.) {
    if((pTrack->Phi() >= fPhiMin2*fPi/180.) && (pTrack->Phi() < fPhiMax2*fPi/180.) && gRandom->Uniform(0,1) > fProbability2) {
      bAccept = kFALSE; // particle is rejected in the second non-uniform sector
    }
  }
  return bAccept;
} // end of Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPhi(AliFlowTrackSimple *pTrack);

//====================================================================================================================

Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPt(AliFlowTrackSimple *pTrack)
{
  // For the case of non-uniform efficiency determine in this method if particle is accepted or rejected for a given pT.
  Bool_t bAccept = kTRUE;
  if((pTrack->Pt() >= fNUEPtMin) && (pTrack->Pt() < fNUEPtMax) && gRandom->Uniform(0,1) > fNUEPtProbability) {
    bAccept = kFALSE; // no mercy!
  }
  return bAccept;
} // end of Bool_t MyAliFlowEventSimpleMakerOnTheFly::AcceptPt(AliFlowTrackSimple *pTrack);

//====================================================================================================================

AliFlowEventSimple* MyAliFlowEventSimpleMakerOnTheFly::CreateEventOnTheFly(AliFlowTrackSimpleCuts const *cutsRP, AliFlowTrackSimpleCuts const *cutsPOI)
{
  // Method to create event 'on the fly'.

  // a) Determine the multiplicity of an event;
  // b) Determine the reaction plane of an event;
  // c) If vn fluctuates uniformly event-by-event, sample its value from [fMinVn,fMaxVn];
  // d) Create event 'on the fly';
  // e) Cosmetics for the printout on the screen.

  // a) Determine the multiplicity of an event:
  Int_t iMult = (Int_t)gRandom->Uniform(fMinMult,fMaxMult);
  Double_t Centrality = (1.-1.*(iMult-fMinMult)/(fMaxMult-fMinMult))*100.;

  // b) Determine the reaction plane of an event and set it in distributions:
  Double_t dReactionPlane = gRandom->Uniform(0.,TMath::TwoPi());

  // d) Create event 'on the fly':
  AliFlowEventSimple *pEvent = new AliFlowEventSimple(iMult);
  pEvent->SetReferenceMultiplicity(iMult);
  pEvent->SetCentrality(Centrality);
  pEvent->SetMCReactionPlaneAngle(dReactionPlane);
  Int_t nRPs = 0; // number of particles tagged RP in this event
  Int_t nPOIs = 0; // number of particles tagged POI in this event
  Double_t dPt, dEta, dPhi, dEta2, dPhi2, dPhiB = 0.;
  Double_t dPi = TMath::Pi();
  Int_t dCharge, dCharge2;
  Bool_t bSkip = kFALSE;
  Bool_t bPass = kTRUE;

  Int_t Gran = 2;
  TH2D* fZDCA = new TH2D("","",Gran,-3.5,3.5,Gran,-3.5,3.5);
  TH2D* fZDCC = new TH2D("","",Gran,-3.5,3.5,Gran,-3.5,3.5);
  Double_t QSpread = 0.5;
  Double_t QACentr[2]={-0.03,-0.02}, QCCentr[2]={0.02,-0.02};
  Double_t fProbZDCA[4]={0.8,1.,1.,0.9}, fProbZDCC[4]={1.,0.9,0.8,1.};

  // Generate particles
  for(Int_t p=0;p<iMult;p++) {

    AliFlowTrackSimple *pTrack = new AliFlowTrackSimple();
    bSkip = kFALSE;
    dPt = fPtSpectra->GetRandom();
    pTrack->SetPt(dPt);
    dCharge = (gRandom->Integer(2)>0.5 ? 1 : -1);
    pTrack->SetCharge(dCharge);

    if (fUseVZERO) {
      if(gRandom->Integer(2)>0.5) {
        dEta = (gRandom->Integer(2)>0.5 ? gRandom->Uniform(-10.,-0.8) : gRandom->Uniform(0.8,10.));
        if(dEta>0.) dPhi = fPhiDistrEtaVA->GetRandom();
        if(dEta<0.) dPhi = fPhiDistrEtaVC->GetRandom();
      } else {
        dEta = gRandom->Uniform(-0.8,0.8);
        dPhi = GetPhi(dEta,dPt,dCharge);
      }
    } else {
      dEta = gRandom->Uniform(-0.8,0.8);
      dPhi = GetPhi(dEta,dPt,dCharge);
    }

    pTrack->SetEta(dEta);

    dPhi += dReactionPlane;
    if(dPhi>TMath::TwoPi()) dPhi -= TMath::TwoPi();

    pTrack->SetPhi(dPhi);

    // Check uniform acceptance:
    if(!fUniformAcceptance && !fChargeDepNUA && !this->AcceptPhi(pTrack)) {continue;}
    if(!fUniformAcceptance && fChargeDepNUA && !this->AcceptPhiCharge(pTrack,dCharge,dEta)) {continue;}
    // Check pT efficiency:
    if(!fUniformEfficiency && !this->AcceptPt(pTrack)){continue;}

    // Checking the RP cuts:
    if(cutsRP->PassesCuts(pTrack)) {
      pTrack->TagRP(kTRUE);
      nRPs++;
    }
    // Checking the POI cuts:
    if(cutsPOI->PassesCuts(pTrack)) {
      pTrack->TagPOI(kTRUE);
      nPOIs++;
    }
    // Assign particles to eta subevents (needed only for Scalar Product method):
    if(pTrack->Eta()>=fEtaMinA && pTrack->Eta()<fEtaMaxA) pTrack->SetForSubevent(0);
    if(pTrack->Eta()>=fEtaMinB && pTrack->Eta()<fEtaMaxB) pTrack->SetForSubevent(1);

    pEvent->AddTrack(pTrack);

    for(Int_t i=1; i<fNTimes; i++) {
      AliFlowTrackSimple *pTrackclone = new AliFlowTrackSimple(*pTrack);
      pTrackclone->SetEta(gRandom->Gaus(dEta,0.6));
      // Checking the RP cuts:
      if(cutsRP->PassesCuts(pTrackclone)) {
        pTrackclone->TagRP(kTRUE);
        nRPs++;
      }
      // Checking the POI cuts:
      if(cutsPOI->PassesCuts(pTrackclone)) {
        pTrackclone->TagPOI(kTRUE);
        nPOIs++;
      }
      pEvent->AddTrack(pTrackclone);
    }

  } // end of for(Int_t p=0;p<iMult;p++)

  pEvent->SetNumberOfRPs(nRPs);
  pEvent->SetNumberOfPOIs(nPOIs);

  // e) Cosmetics for the printout on the screen:
  Int_t cycle = (fPtDependentV2 ? 10 : 1000);
  if((++fCount % cycle) == 0) {
    if(TMath::Abs(dReactionPlane)>1.e-44) {
      cout<<" MC Reaction Plane Angle = "<<dReactionPlane<<endl;
    } else {
      cout<<" MC Reaction Plane Angle is unknown :'( "<< endl;
    }
    cout<<" # of simulated tracks  = "<<iMult<<endl;
    cout<<" # of RP tagged tracks  = "<<nRPs<<endl;
    cout<<" # of POI tagged tracks = "<<nPOIs<<endl;
    cout <<"  .... "<<fCount<< " events processed ...."<<endl;
  } // end of if((++fCount % cycle) == 0)

 Double_t PsiA = dReactionPlane;
 Double_t PsiC = dReactionPlane+TMath::Pi();
 if(PsiC>TMath::TwoPi()) PsiC -= TMath::TwoPi();
 pEvent->SetMCReactionPlaneAngle(dReactionPlane);

 // fuck resolution
 PsiA += gRandom->Gaus(0.,TMath::Pi()/1.5); //1,5
 PsiC += gRandom->Gaus(0.,TMath::Pi()/1.5);

 Double_t QA[2]={0.,0.}, QC[2]={0.,0.};
 QA[0] = TMath::Cos(PsiA);
 QA[1] = TMath::Sin(PsiA);
 QC[0] = TMath::Cos(PsiC);
 QC[1] = TMath::Sin(PsiC);

 pEvent->SetZDC2Qsub(QC,10.,QA,10.);

  // set 2Qsub as ZDC2Qsub
  if(fUseVZERO) {
    Double_t QA[2]={0.,0.}, QC[2]={0.,0.}, j=1.;
    Int_t MA=0, MC=0;
    Double_t p[2] = {-1.75, 1.75};
    Int_t Tile = 1.;
    Double_t AvMultTow[8]={330.855,
      229.08,
      381.796,
      330.195,
      351.635,
      340.737,
      315.291,
      248.372};
    Double_t AvMultA=0.,AvMultC=0.;
    for(Int_t i=0; i<4; i++) {
      AvMultC += AvMultTow[i];
      AvMultA += AvMultTow[i+4];
    }
    AvMultC /= 4.;
    AvMultA /= 4.;

    for (Int_t i=1; i<=Gran; i++) {
      for (Int_t k=1; k<=2; k++) {
        if(fZDCC->GetBinContent(i,k)>0.) {
          // correct !
          //    j = AvMultC/AvMultTow[Tile-1];
          Double_t MTow = fZDCC->GetBinContent(i,k)*j;

          QC[0] += p[i-1]*MTow;
          QC[1] += p[k-1]*MTow;
          MC    += MTow;
          fZDCMult->Fill(Tile-0.5,MTow);
          Tile++;
        }
      }
    }
    if(MC>0.) {
      QC[0] /= MC;
      QC[1] /= MC;
    }
    for (Int_t i=1; i<=2; i++) {
      for (Int_t k=1; k<=2; k++) {
        if(fZDCA->GetBinContent(i,k)>0.) {
          // correct !
          //    j = AvMultA/AvMultTow[Tile-1];
          Double_t MTow = fZDCA->GetBinContent(i,k)*j;

          QA[0] += p[i-1]*MTow;
          QA[1] += p[k-1]*MTow;
          MA    += MTow;
          fZDCMult->Fill(Tile-0.5,MTow);
          Tile++;
        }
      }
    }
    if(MA>0.) {
      QA[0] /= MA;
      QA[1] /= MA;
    }

    pEvent->SetZDC2Qsub(QC,MC,QA,MA);

    //// AliFlowVector* vQarray = new AliFlowVector[2];
    //// Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
    //// pEvent->Get2Qsub(vQarray,1);
    ////  AliFlowVector vQC = vQarray[0];
    ////  AliFlowVector vQA = vQarray[1];
    ////  Double_t dMC = vQC.GetMult();
    ////  Double_t dMA = vQA.GetMult();
    ////  xyZNC[0] = -vQC.X()/dMC;
    ////  xyZNC[1] = -vQC.Y()/dMC;
    ////  xyZNA[0] = -vQA.X()/dMA;
    ////  xyZNA[1] = -vQA.Y()/dMA;
    ////  pEvent->SetZDC2Qsub(xyZNC,dMC,xyZNA,dMA);
  }

  return pEvent;

} // end of CreateEventOnTheFly()

//====================================================================================================================

Double_t MyAliFlowEventSimpleMakerOnTheFly::GetDist(Double_t A, Double_t B)
{
  Double_t Dist = TMath::ATan2(TMath::Sin(A-B),TMath::Cos(A-B));
  if(Dist > TMath::Pi()/2.)  Dist = TMath::Pi() - Dist;
  if(Dist < -TMath::Pi()/2.) Dist = -TMath::Pi() - Dist;
  return Dist;
}

//====================================================================================================================

Double_t MyAliFlowEventSimpleMakerOnTheFly::GetPhi(Double_t dEta, Double_t dPt, Int_t dCharge)
{
  Double_t dPhi = 0.;
  if (!fSimCRC) {
    dPhi = fPhiDistribution->GetRandom();
  }
  else {
    if (!fSimEtaDep && !fSimPtDep) {
      if(dEta>0. && dCharge==1)  dPhi = fPhiDistributionP->GetRandom();
      if(dEta>0. && dCharge==-1) dPhi = fPhiDistributionN->GetRandom();
      if(dEta<0. && dCharge==1)  dPhi = fPhiDistributionN->GetRandom();
      if(dEta<0. && dCharge==-1) dPhi = fPhiDistributionP->GetRandom();
    }
    if(fSimEtaDep && !fSimPtDep) {
      Int_t e = GetBinEta(dEta);
      if(dCharge==1)  dPhi = fPhiDistrEtaP[e]->GetRandom();
      if(dCharge==-1) dPhi = fPhiDistrEtaN[e]->GetRandom();
    }
    if(!fSimEtaDep && fSimPtDep) {
      Int_t p = (Int_t)((dPt-fPtMin)/fPtBinWidth);
      if(dEta>0. && dCharge==1)  dPhi = fPhiDistrPtP[p]->GetRandom();
      if(dEta>0. && dCharge==-1) dPhi = fPhiDistrPtN[p]->GetRandom();
      if(dEta<0. && dCharge==1)  dPhi = fPhiDistrPtN[p]->GetRandom();
      if(dEta<0. && dCharge==-1) dPhi = fPhiDistrPtP[p]->GetRandom();
    }
  }
  if(fSimCME) {
    if(dEta<0. && dCharge==1)  dPhi = fPhiDistrEtaP[0]->GetRandom();
    if(dEta<0. && dCharge==-1) dPhi = fPhiDistrEtaN[0]->GetRandom();
    if(dEta>0. && dCharge==1)  dPhi = fPhiDistrEtaP[1]->GetRandom();
    if(dEta>0. && dCharge==-1) dPhi = fPhiDistrEtaN[1]->GetRandom();
  }
  return dPhi;
}

//====================================================================================================================

Int_t MyAliFlowEventSimpleMakerOnTheFly::GetBinEta(Double_t dEta)
{
  Int_t e = RefEtaHist->FindBin(dEta)-1;
  return e;
}

//====================================================================================================================

Double_t MyAliFlowEventSimpleMakerOnTheFly::GetEtafromBin(Int_t dEta)
{
  Double_t e = (dEta+0.5)*fEtaBinWidth-0.8;
  return e;
}
