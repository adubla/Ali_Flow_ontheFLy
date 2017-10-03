/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

/************************************
 * Create an event and perform full *
 * flow analysis 'on the fly'.      *
 *                                  *
 * author: Ante Bilandzic           *
 *         (abilandzic@gmail.com)   *
 ************************************/

#ifndef MYALIFLOWEVENTSIMPLEMAKERONTHEFLY_H
#define MYALIFLOWEVENTSIMPLEMAKERONTHEFLY_H

class TF1;
class TRandom3;
class TH3F;
class TProfile;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowTrackSimpleCuts;

class MyAliFlowEventSimpleMakerOnTheFly{
public:
 MyAliFlowEventSimpleMakerOnTheFly(UInt_t uiSeed = 0); // constructor
 virtual ~MyAliFlowEventSimpleMakerOnTheFly(); // destructor
 virtual void Init();
 Bool_t AcceptPhi(AliFlowTrackSimple *pTrack);
 Bool_t AcceptPhiCharge(AliFlowTrackSimple *pTrack, Int_t dCharge, Double_t dEta);
 Bool_t AcceptPt(AliFlowTrackSimple *pTrack);
 AliFlowEventSimple* CreateEventOnTheFly(AliFlowTrackSimpleCuts const *cutsRP, AliFlowTrackSimpleCuts const *cutsPOI);
 Double_t GetPhi(Double_t dEta, Double_t dPt, Int_t dCharge);
 Double_t GetDist(Double_t A, Double_t B);
 Int_t GetBinEta(Double_t eta);
 Double_t GetEtafromBin(Int_t dEta);
 // Setters and getters:
 void SetMinMult(Int_t iMinMult) {this->fMinMult = iMinMult;}
 Int_t GetMinMult() const {return this->fMinMult;}
 void SetMaxMult(Int_t iMaxMult) {this->fMaxMult = iMaxMult;}
 Int_t GetMaxMult() const {return this->fMaxMult;}
 void SetMass(Double_t dMass) {this->fMass = dMass;}
 Double_t GetMass() const {return this->fMass;}
 void SetTemperature(Double_t dT) {this->fTemperature = dT;}
 Double_t GetTemperature() const {return this->fTemperature;}
 void SetV1(Double_t dV1) {this->fV1 = dV1;}
 Double_t GetV1() const {return this->fV1;}
 void SetV2(Double_t dV2) {this->fV2 = dV2;}
 Double_t GetV2() const {return this->fV2;}
 void SetV3(Double_t dV3) {this->fV3 = dV3;}
 Double_t GetV3() const {return this->fV3;}
 void SetV4(Double_t dV4) {this->fV4 = dV4;}
 Double_t GetV4() const {return this->fV4;}
 void SetV5(Double_t dV5) {this->fV5 = dV5;}
 Double_t GetV5() const {return this->fV5;}
 void SetV6(Double_t dV6) {this->fV6 = dV6;}
 Double_t GetV6() const {return this->fV6;}
 void SetUniformFluctuationsV1(Bool_t b) {this->fUniformFluctuationsV1 = b;}
 Bool_t GetUniformFluctuationsV1() const {return this->fUniformFluctuationsV1;}
 void SetMinV1(Double_t dMinV1) {this->fMinV1 = dMinV1;}
 Double_t GetMinV1() const {return this->fMinV1;}
 void SetMaxV1(Double_t dMaxV1) {this->fMaxV1 = dMaxV1;}
 Double_t GetMaxV1() const {return this->fMaxV1;}
 void SetUniformFluctuationsV2(Bool_t b) {this->fUniformFluctuationsV2 = b;}
 Bool_t GetUniformFluctuationsV2() const {return this->fUniformFluctuationsV2;}
 void SetMinV2(Double_t dMinV2) {this->fMinV2 = dMinV2;}
 Double_t GetMinV2() const {return this->fMinV2;}
 void SetMaxV2(Double_t dMaxV2) {this->fMaxV2 = dMaxV2;}
 Double_t GetMaxV2() const {return this->fMaxV2;}
 void SetV1B(Double_t V1B) {this->fV1B = V1B;}
 Double_t GetV1B() const {return this->fV1B;}
 void SetV3B(Double_t V3B) {this->fV3B = V3B;}
 Double_t GetV3B() const {return this->fV3B;}
 void SetV1CME(Double_t CME) {this->fCME = CME;}
 Double_t GetV1CME() const {return this->fCME;}
 void SetSimCRC(Bool_t b) {this->fSimCRC = b;}
 Bool_t GetSimCRC() const {return this->fSimCRC;}
 void SetSimCRCv3(Bool_t b) {this->fSimCRCv3 = b;}
 Bool_t GetSimCRCv3() const {return this->fSimCRCv3;}
 void SetSimCME(Bool_t b) {this->fSimCME = b;}
 Bool_t GetSimCME() const {return this->fSimCME;}
 void SetSimEtaDep(Bool_t b) {this->fSimEtaDep = b;}
 Bool_t GetSimEtaDep() const {return this->fSimEtaDep;}
 void SetSimPtDep(Bool_t b) {this->fSimPtDep = b;}
 Bool_t GetSimPtDep() const {return this->fSimPtDep;}
 void SetChargeDepNUA(Bool_t b) {this->fChargeDepNUA = b;}
 Bool_t GetChargeDepNUA() const {return this->fChargeDepNUA;}
 void SetFuckEP(Bool_t b) {this->fFuckEP = b;}
 Bool_t GetFuckEP() const {return this->fFuckEP;}
 TH1D* GetBalanFunc(Int_t i) const {return this->BalanFunc[i];}
 
 void SetGausFluctuationsV1(Bool_t b) {this->fGausFluctuationsV1 = b;}
 Bool_t GetGausFluctuationsV1() const {return this->fGausFluctuationsV1;}
 void SetMeanSigV1(Double_t dMeanV1, Double_t dSigV1) {this->fMeanV1 = dMeanV1; this->fSigV1 = dSigV1;}
 Double_t GetSigV1() const {return this->fSigV1;}
 
 void SetPtDependentV2(Bool_t b) {this->fPtDependentV2 = b;}
 Bool_t GetPtDependentV2() const {return this->fPtDependentV2;}
 void SetV2vsPtCutOff(Double_t dV2vsPtCutOff) {this->fV2vsPtCutOff = dV2vsPtCutOff;}
 Double_t GetV2vsPtCutOff() const {return this->fV2vsPtCutOff;}
 void SetV2vsPtMax(Double_t dV2vsPtMax) {this->fV2vsPtMax = dV2vsPtMax;}
 Double_t GetV2vsPtMax() const {return this->fV2vsPtMax;}
 void SetSubeventEtaRange(Double_t minA, Double_t maxA, Double_t minB, Double_t maxB)
 {this->fEtaMinA = minA;this->fEtaMaxA = maxA;this->fEtaMinB = minB;this->fEtaMaxB = maxB;};
 void SetNTimes(Int_t nt) {this->fNTimes = nt;}
 Int_t GetNTimes() const {return this->fNTimes;}
 void SetUniformAcceptance(Bool_t ua) {this->fUniformAcceptance = ua;}
 Bool_t GetUniformAcceptance() const {return this->fUniformAcceptance;}
 void SetFirstSectorPhiMin(Double_t dPhiMin1) {this->fPhiMin1 = dPhiMin1;}
 Double_t GetFirstSectorPhiMin() const {return this->fPhiMin1;}
 void SetFirstSectorPhiMax(Double_t dPhiMax1) {this->fPhiMax1 = dPhiMax1;}
 Double_t GetFirstSectorPhiMax() const {return this->fPhiMax1;}
 void SetFirstSectorProbability(Double_t dProbability1) {this->fProbability1 = dProbability1;}
 Double_t GetFirstSectorProbability() const {return this->fProbability1;}
 void SetSecondSectorPhiMin(Double_t dPhiMin2) {this->fPhiMin2 = dPhiMin2;}
 Double_t GetSecondSectorPhiMin() const {return this->fPhiMin2;}
 void SetSecondSectorPhiMax(Double_t dPhiMax2) {this->fPhiMax2 = dPhiMax2;}
 Double_t GetSecondSectorPhiMax() const {return this->fPhiMax2;}
 void SetSecondSectorProbability(Double_t dProbability2) {this->fProbability2 = dProbability2;}
 Double_t GetSecondSectorProbability() const {return this->fProbability2;}
 void SetUniformEfficiency(Bool_t ue) {this->fUniformEfficiency = ue;}
 Bool_t GetUniformEfficiency() const {return this->fUniformEfficiency;}
 void SetNUEPtMin(Double_t ptMin) {this->fNUEPtMin = ptMin;}
 Double_t GetNUEPtMin() const {return this->fNUEPtMin;}
 void SetNUEPtMax(Double_t ptMax) {this->fNUEPtMax = ptMax;}
 Double_t GetNUEPtMax() const {return this->fNUEPtMax;}
 void SetNUEPtProbability(Double_t ptp) {this->fNUEPtProbability = ptp;}
 Double_t GetNUEPtProbability() const {return this->fNUEPtProbability;}
 void SetUseVZERO(Bool_t ua) {this->fUseVZERO = ua;}
 Bool_t GetUseVZERO() const {return this->fUseVZERO;}
 TProfile* GetZDCMult() const {return this->fZDCMult;}
 
private:
 MyAliFlowEventSimpleMakerOnTheFly(const MyAliFlowEventSimpleMakerOnTheFly& anAnalysis); // copy constructor
 MyAliFlowEventSimpleMakerOnTheFly& operator=(const MyAliFlowEventSimpleMakerOnTheFly& anAnalysis); // assignment operator
 Int_t fCount; // count number of events
 Int_t fMinMult; // uniformly sampled multiplicity is >= iMinMult
 Int_t fMaxMult; // uniformly sampled multiplicity is < iMaxMult
 TF1 *fPtSpectra; // transverse momentum distribution (pt is sampled from hardwired Boltzmann distribution)
 Double_t fMass; // mass in pt distribution (hardwired is Boltzmann pt distribution)
 Double_t fTemperature; // "temperature" in pt distribution (hardwired is Boltzmann pt distribution)
 TF1 *fPhiDistribution; // azimuthal distribution (phi is sampled from hardwired Fourier-like distribution)
 TF1 *fPhiDistributionP;
 TF1 *fPhiDistributionN;
 TF1 *fPsiDistribution;
 Double_t fV1; // harmonic v1
 Double_t fV2; // harmonic v2
 Double_t fV3; // harmonic v3
 Double_t fV4; // harmonic v4
 Double_t fV5; // harmonic v5
 Double_t fV6; // harmonic v6
 Bool_t fUniformFluctuationsV1; // v1 is sampled uniformly for each event and for all particles from [fMinV1,fMaxV1]
 Double_t fMinV1; // if v1 is sampled uniformly for each event, this is lower boundary on its value
 Double_t fMaxV1; // if v1 is sampled uniformly for each event, this is upper boundary on its value
 Bool_t fUniformFluctuationsV2; // v2 is sampled uniformly for each event and for all particles from [fMinV2,fMaxV2]
 Double_t fMinV2; // if v2 is sampled uniformly for each event, this is lower boundary on its value
 Double_t fMaxV2; // if v2 is sampled uniformly for each event, this is upper boundary on its value
 Bool_t fGausFluctuationsV1; // v1 is sampled uniformly for each event and for all particles from [fMinV1,fMaxV1]
 Double_t fSigV1; //
 Double_t fMeanV1; //
 Bool_t fPtDependentV2; // v2 is pt-dependent
 Double_t fV2vsPtCutOff; // if v2 is pt-dependent: for v2 < fV2vsPtCutOff v2 is growing linearly, otherwise v2 = fV2vsPtMax
 Double_t fV2vsPtMax; // if v2 is pt-dependent: v2 = fV2vsPtMax for v2 >= fV2vsPtCutOff
 Double_t fEtaMinA; // minimum eta of subevent A
 Double_t fEtaMaxA; // maximum eta of subevent A
 Double_t fEtaMinB; // minimum eta of subevent B
 Double_t fEtaMaxB; // maximum eta of subevent B
 Int_t fNTimes; // number of times to use the same particle in the analysis (simulating nonflow)
 Bool_t fUniformAcceptance; // detector has uniform azimuthal acceptance or not
 Double_t fPhiMin1; // first sector with non-uniform acceptance starts at azimuth fPhiMin1
 Double_t fPhiMax1; // first sector with non-uniform acceptance ends at azimuth fPhiMax1
 Double_t fProbability1; // particles emitted in fPhiMin1 < phi < fPhiMax1 are taken with probability fProbability1
 Double_t fPhiMin2; // second sector with non-uniform acceptance starts at azimuth fPhiMin2
 Double_t fPhiMax2; // second sector with non-uniform acceptance ends at azimuth fPhiMax2
 Double_t fProbability2; // particles emitted in fPhiMin2 < phi < fPhiMax2 are taken with probability fProbability2
 Double_t fPi; // pi
 Bool_t fUniformEfficiency; // detector has uniform efficiency vs pT, or perhaps not...
 Double_t fNUEPtMin; // non-uniform efficiency vs pT starts at pT = fPtMin
 Double_t fNUEPtMax; // non-uniform efficiency vs pT ends at pT = fPtMax
 Double_t fNUEPtProbability; // particles emitted in fPtMin <= pT < fPtMax are taken with probability fPtProbability
 Double_t fV1B;
 Double_t fV3B;
 Double_t fV1EbE;
 Double_t fCME;
 Bool_t fSimCRC;
 Bool_t fSimCRCv3;
 Bool_t fSimCME;
 Bool_t fSimEtaDep;
 Bool_t fUseVZERO;
 const static Int_t fnEtaBins = 32;
 Double_t fEtaMin;
 Double_t fEtaMax;
 Double_t fEtaRange;
 Double_t fEtaBinWidth;
 Bool_t fSimPtDep;
 Bool_t fChargeDepNUA;
 Bool_t fFuckEP;
 const static Int_t fnPtBins = 100;
 Double_t fPtMin;
 Double_t fPtMax;
 Double_t fPtRange;
 Double_t fPtBinWidth;
 TF1 *fPhiDistrEtaP[fnEtaBins];
 TF1 *fPhiDistrEtaN[fnEtaBins];
 TF1 *fPhiDistrPtP[fnPtBins];
 TF1 *fPhiDistrPtN[fnPtBins];
 TF1 *fPhiDistrEtaVA;
 TF1 *fPhiDistrEtaVC;
 const static Int_t fnvonMis = 100;
 TF1 *fvonMisesDistrib[fnvonMis];
 TH1D *BalanFunc[2];
 TProfile *fZDCMult;
 TH1D *RefEtaHist;
 
 ClassDef(MyAliFlowEventSimpleMakerOnTheFly,1) // macro for rootcint
};

#endif



