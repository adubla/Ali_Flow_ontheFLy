// Settings for the simulation of events 'on the fly':
//  a) Determine how many events you want to create;
//  b) Set random or same seed for random generator;
//  c) Determine multiplicites of events;
//  d) Parametrize the phi distribution;
//   d1) Enable/disable uniform event-wise fluctuations of v2;
//   d2) Enable/diable pt dependence of v2;
//  e) Parametrize the pt distribution;
//  f) Determine how many times each sampled particle will be taken (simulating nonflow);
//  g) Configure detector's:
//   g1) acceptance;
//   g2) efficiency;
//  h) Decide which flow analysis methods you will use;
//  i) Define simple cuts for Reference Particle (RP) selection;
//  j) Define simple cuts for Particle of Interest (POI) selection;
//  k) Define the ranges for two subevents separated with eta gap (needed only for SP method);
//  l) Enable/disable usage of particle weights.

// folder:
//TString folderName = "output/testEG";

// a) Determine how many events you want to create:
Int_t iNevts = (Int_t)437000*0.6; //120000 total statistics
// b) Set random or same seed for random generator:
Bool_t bSameSeed = kFALSE; // if kTRUE, the created events are the same when re-doing flow analysis 'on the fly'

// c) Determine multiplicites of events:
//    Remark 1: Multiplicity M for each event is sampled uniformly from interval iMinMult <= M < iMaxMult;
//    Remark 2: For constant M of e.g. 500 for each event, set iMinMult = 500 and iMaxMult = 501.
Int_t iMinMult = 1000; // uniformly sampled multiplicity is >= iMinMult
Int_t iMaxMult = 1001; // uniformly sampled multiplicity is < iMaxMult

Int_t nCen = 1;
Int_t wCen = 100.;

// Enable/disable directed flow due to randomly directed magnetic field:
Bool_t bSimulateCRC = kTRUE;
Double_t dV1B = 1.5E-2;
Bool_t bSimulateCRCv3 = kFALSE;
Double_t dV3B = 0.; // NB: only bSimEtaDep = kFALSE;

// Enable/disable eta or pt dependency:
Bool_t bSimEtaDep = kTRUE;
Bool_t bSimPtDep = kFALSE;
Bool_t bRunDepNUA = kFALSE;
Bool_t bChargeDepNUA = kFALSE;
Bool_t bFuckEP = kFALSE;

// Set what to calculate:
Bool_t bUseVZERO = kFALSE;
Bool_t bUseZDC = kTRUE;
Bool_t bCalculateCRC2 = kFALSE;
Bool_t bCalculateCRCVZ = kFALSE;
Int_t  nCRC2EtaBins=10;
Bool_t bCalculateCRCPt = kFALSE;
Bool_t bUseCRCRecentering = kFALSE;
Bool_t bRecenderZDC = kFALSE;
Bool_t bNUACorr = kTRUE; // correct for NUA
TString QVecWeightsFileName = "QVecWeights.root";

// Enable/disable usage of particle weights: PhiEtaWeights.root
Bool_t usePhiWeights = kFALSE; // phi weights
Bool_t usePtWeights  = kFALSE; // pt weights
Bool_t useEtaWeights = kFALSE; // eta weights
Bool_t usePhiEtaWeights = kFALSE; // eta weights
TString PhiEtaWeightsFileName = "PhiEtaWeights.root";

// Simulate CME:
Bool_t bCalculateCME = kFALSE;
Bool_t bSimCME = kFALSE;
Double_t dCME = 0.;

// Set Harmonic
Int_t iHarmonic = 1;

// Flow SP
Bool_t bCalculateFlow = kTRUE;

// d) Parametrize the phi distribution:
//    Remark 1: Hardwired is Fourier-like distribution f(phi) = (1/2pi)(1+sum_{n=1}^{4} 2v_n cos[n(phi-rp)]),
//              where reaction plane (rp) is sampled uniformly for each event from interval [0,2pi]
Double_t dV1 = 0.;//-3.2E-1; // constant harmonic v1
Double_t dV2 = 1.E-1; // constant harmonic v2
Double_t dV3 = 1.E-2; // constant harmonic v3
Double_t dV4 = 0.0; // constant harmonic v4
Double_t dV5 = 0.0; // constant harmonic v5
Double_t dV6 = 0.0; // constant harmonic v6
//    Remark 2: By default all harmonics are constant for each event and for each particle. However, for v2
//              the uniform event-wise fluctuations or pt dependence can be enabled:
//  d1) Enable/disable uniform event-wise fluctuations of v2:
Bool_t bUniformFluctuationsV1 = kFALSE; // enable uniform event-wise flow fluctuations (set than also dMinV2 and dMaxV2 bellow)
Double_t dMinV1 = 0.15; // lower boundary on v2, when bUniformFluctuationsV2 = kTRUE
Double_t dMaxV1 = 0.25; // upper boundary on v2, when bUniformFluctuationsV2 = kTRUE
Bool_t bUniformFluctuationsV2 = kFALSE; // enable uniform event-wise flow fluctuations (set than also dMinV2 and dMaxV2 bellow)
Double_t dMinV2 = 0.04; // lower boundary on v2, when bUniformFluctuationsV2 = kTRUE
Double_t dMaxV2 = 0.06; // upper boundary on v2, when bUniformFluctuationsV2 = kTRUE
Bool_t bGausFluctuationsV1 = kFALSE; // enable uniform event-wise flow fluctuations (set than also dMinV2 and dMaxV2 bellow)
Double_t dSigV1 = 0.05; // lower boundary on v2, when bUniformFluctuationsV2 = kTRUE
//  d2) Enable/disable pt dependence of v2:
Bool_t bPtDependentV2 = kFALSE; // enable pt dependence of v2 (set then also dV2vsPtMax and dV2vsPtCutOff bellow)
Double_t dV2vsPtCutOff = 2.0; // up to pt = dV2vsPtCutOff v2 is growing linearly as a function of pt
Double_t dV2vsPtMax = 0.20; // for pt >= dV2vsPtCutOff, v2(pt) = dV2vsPtMax

// e) Parametrize the pt distribution:
//    Remark: Hardwired is Boltzmann distribution f(pt) = pt*exp[-sqrt(dMass^2+pt^2)/dT]
Double_t dMass = 0.13957; // mass in GeV/c^2 (e.g. m_{pions} = 0.13957)
Double_t dTemperature = 0.44; // "temperature" in GeV/c (increase this parameter to get more high pt particles)

// f) Determine how many times each sampled particle will be taken in the analysis (simulating nonflow):
Int_t nTimes = 1; // e.g. for nTimes = 2, strong 2-particle nonflow correlations are introduced

// g1) Configure detector's acceptance:
Bool_t uniformAcceptance = kTRUE; // if kTRUE: detectors has uniform azimuthal acceptance.
// if kFALSE: you will simulate detector with non-uniform acceptance in one or
// two sectors. For each sector you specify phiMin, phiMax and probability p.
// Then all particles emitted in direction phiMin < phi < phiMax will be taken
// with probability p. If p = 0, that sector is completely blocked. Set bellow
// phiMin1, phiMax1, p1 for the first sector and phiMin2, phiMax2, p2 for the second
// sector. If you set phiMin2 = phiMax2 = p2 = 0, only first non-uniform sector is
// simulated.
// 1st non-uniform sector:
Double_t phiMin1 = 0; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax1 = 90; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p1 = 0.5; // probability that particles emitted in [phiMin1,phiMax1] are taken
// 2nd non-uniform sector:
Double_t phiMin2 = 180.; // first non-uniform sector starts at this azimuth (in degrees)
Double_t phiMax2 = 315.; // first non-uniform sector ends at this azimuth (in degrees)
Double_t p2 = 0.8; // probablitity that particles emitted in [phiMin2,phiMax2] are taken

// g2) Configure detector's efficiency:
Bool_t uniformEfficiency = kTRUE; // if kTRUE: detectors has uniform pT efficiency
// if kFALSE: you will simulate detector with non-uniform pT efficiency.
// Then all particles emitted in ptMin <= pt < ptMax will be taken
// with probability p, to be specified in lines just below.
Double_t NUEptMin = 0.8; // non-uniform efficiency vs pT starts at pT = fPtMin
Double_t NUEptMax = 1.2; // non-uniform efficiency vs pT ends at pT = fPtMax
Double_t NUEp = 0.5; // probablitity that particles emitted in [ptMin,ptMax> are taken

// h) Decide which flow analysis methods you will use:
Bool_t MCEP    = kTRUE; // Monte Carlo Event Plane
Bool_t QCwPOI  = kTRUE; // Q-cumulants POI-POI

// i) Define simple cuts for Reference Particle (RP) selection:
Double_t ptMinRP = 0.2; // in GeV
Double_t ptMaxRP = 8.0; // in GeV
Double_t etaMinRP = -10.;
Double_t etaMaxRP = 10.;
Double_t phiMinRP = 0.0; // in degrees
Double_t phiMaxRP = 360.0; // in degrees
Bool_t bUseChargeRP = kFALSE; // if kFALSE, RPs with both sign of charges are taken
Int_t chargeRP = 1; // +1 or -1

// j) Define simple cuts for Particle of Interest (POI) selection:
Double_t ptMinPOI = 0.2; // in GeV
Double_t ptMaxPOI = 8.0; // in GeV
Double_t etaMinPOI = -0.8; //
Double_t etaMaxPOI = 0.8;
Double_t phiMinPOI = 0.0; // in degrees
Double_t phiMaxPOI = 360.0; // in degrees
Bool_t bUseChargePOI = kFALSE; // if kFALSE, POIs with both sign of charges are taken
Int_t chargePOI = -1; // +1 or -1

// k) Define the ranges for two subevents separated with eta gap (needed only for SP method):
Double_t etaMinA = -10.; // minimum eta of subevent A
Double_t etaMaxA = -0.8; // maximum eta of subevent A
Double_t etaMinB = 0.8; // minimum eta of subevent B
Double_t etaMaxB = 10.; // maximum eta of subevent B

enum anaModes {mLocal,mLocalSource,mLocalPAR};
// mLocal: Analyze data on your computer using aliroot
// mLocalPAR: Analyze data on your computer using root + PAR files
// mLocalSource: Analyze data on your computer using root + source files

#include "TStopwatch.h"
#include "TRandom3.h"
#include "TObjArray"
#include "Riostream.h"
#include "TFile.h"

int runFlowAnalysisOnTheFly(Int_t mode=mLocal)
{
 // Beging analysis 'on the fly'.

 // a) Formal necessities....;
 // b) Initialize the flow event maker 'on the fly';
 // c) If enabled, access particle weights from external file;
 // d) Configure the flow analysis methods;
 // e) Simple cuts for RPs;
 // f) Simple cuts for POIs;
 // g) Create and analyse events 'on the fly';
 // h) Create the output file and directory structure for the final results of all methods;
 // i) Calculate and store the final results of all methods.

 // a) Formal necessities....:
 CheckUserSettings();
 LoadLibraries(mode);
 TStopwatch timer;
 timer.Start();

 // b) Initialize the flow event maker 'on the fly':
 UInt_t uiSeed = 0; // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID
 if(bSameSeed){uiSeed = 44;}
 MyAliFlowEventSimpleMakerOnTheFly* eventMakerOnTheFly = new MyAliFlowEventSimpleMakerOnTheFly(uiSeed);
 eventMakerOnTheFly->SetMinMult(iMinMult);
 eventMakerOnTheFly->SetMaxMult(iMaxMult);
 eventMakerOnTheFly->SetMass(dMass);
 eventMakerOnTheFly->SetTemperature(dTemperature);
 eventMakerOnTheFly->SetV1(dV1);
 eventMakerOnTheFly->SetV2(dV2);
 eventMakerOnTheFly->SetV3(dV3);
 eventMakerOnTheFly->SetV4(dV4);
 eventMakerOnTheFly->SetV5(dV5);
 eventMakerOnTheFly->SetV6(dV6);
 if(bUniformFluctuationsV1) {
  eventMakerOnTheFly->SetUniformFluctuationsV1(bUniformFluctuationsV1);
  eventMakerOnTheFly->SetMinV1(dMinV1);
  eventMakerOnTheFly->SetMaxV1(dMaxV1);
 }
 if(bUniformFluctuationsV2) {
  eventMakerOnTheFly->SetUniformFluctuationsV2(bUniformFluctuationsV2);
  eventMakerOnTheFly->SetMinV2(dMinV2);
  eventMakerOnTheFly->SetMaxV2(dMaxV2);
 }
 if(bGausFluctuationsV1) {
  eventMakerOnTheFly->SetGausFluctuationsV1(bGausFluctuationsV1);
  eventMakerOnTheFly->SetMeanSigV1(dV1,dSigV1);
 }
 if(bPtDependentV2) {
  eventMakerOnTheFly->SetPtDependentV2(bPtDependentV2);
  eventMakerOnTheFly->SetV2vsPtCutOff(dV2vsPtCutOff);
  eventMakerOnTheFly->SetV2vsPtMax(dV2vsPtMax);
 }
 eventMakerOnTheFly->SetSubeventEtaRange(etaMinA,etaMaxA,etaMinB,etaMaxB);
 eventMakerOnTheFly->SetNTimes(nTimes);
 if(!uniformAcceptance) {
  eventMakerOnTheFly->SetUniformAcceptance(kFALSE);
  eventMakerOnTheFly->SetFirstSectorPhiMin(phiMin1);
  eventMakerOnTheFly->SetFirstSectorPhiMax(phiMax1);
  eventMakerOnTheFly->SetFirstSectorProbability(p1);
  eventMakerOnTheFly->SetSecondSectorPhiMin(phiMin2);
  eventMakerOnTheFly->SetSecondSectorPhiMax(phiMax2);
  eventMakerOnTheFly->SetSecondSectorProbability(p2);
 }
 if(!uniformEfficiency) {
  eventMakerOnTheFly->SetUniformEfficiency(kFALSE);
  eventMakerOnTheFly->SetNUEptMinPtMin(NUEptMin);
  eventMakerOnTheFly->SetNUEptMinPtMax(NUEptMax);
  eventMakerOnTheFly->SetNUEptMinPtProbability(NUEp);
 }
 if(bSimulateCRC) {
  eventMakerOnTheFly->SetSimCRC(bSimulateCRC);
  eventMakerOnTheFly->SetV1B(dV1B);
  eventMakerOnTheFly->SetSimCRCv3(bSimulateCRCv3);
  eventMakerOnTheFly->SetV3B(dV3B);
  eventMakerOnTheFly->SetSimEtaDep(bSimEtaDep);
  eventMakerOnTheFly->SetSimPtDep(bSimPtDep);
  eventMakerOnTheFly->SetSimCME(bSimCME);
  eventMakerOnTheFly->SetV1CME(dCME);
 }
 if(bUseVZERO) {
  eventMakerOnTheFly->SetUseVZERO(kTRUE);
 }
 if(bChargeDepNUA) {
  eventMakerOnTheFly->SetChargeDepNUA(bChargeDepNUA);
 }
 eventMakerOnTheFly->SetFuckEP(bFuckEP);
 eventMakerOnTheFly->Init();

 // welcome message with settings
 cout<<endl;
 cout<<"      ---- ARE YOU READY TO FLY ? ----      "<<endl;
 cout<<endl;
 cout<<" ---- BEGIN FLOW ANALYSIS 'ON THE FLY' ---- "<<endl;
 cout<<endl;
 cout<<" ---- settings: ---- "<<endl;
 cout<<endl;
 cout<<" v_1: "<<eventMakerOnTheFly->GetV1()<<" (";
 if (!eventMakerOnTheFly->GetUniformFluctuationsV1()) {cout<<"costant)"<<endl;}
 else {cout<<"fluctuating)"<<endl;}
 cout<<" non-flow: ";
 if (eventMakerOnTheFly->GetNTimes()>1) {cout<<"ON"<<endl;}
 else {cout<<"OFF"<<endl;}
 cout<<" azimuthal acceptance: ";
 if (eventMakerOnTheFly->GetUniformAcceptance() && !bRunDepNUA) {cout<<"uniform"<<endl;}
 else {cout<<"non-uniform"<<endl;}
 cout<<" magnetic field: ";
 if (!eventMakerOnTheFly->GetSimCRC()) {cout<<"OFF"<<endl;}
 else {
  cout<<"ON"<<endl;
  cout<<"  v_1B: "<<eventMakerOnTheFly->GetV1B()<<endl;
  cout<<"  eta-dep: ";
  if (!eventMakerOnTheFly->GetSimEtaDep()) {cout<<"OFF"<<endl;}
  else {cout<<"ON"<<endl;}
  cout<<"  pt-dep: ";
  if (!eventMakerOnTheFly->GetSimPtDep()) {cout<<"OFF"<<endl;}
  else {cout<<"ON"<<endl;}
 }
 cout<<endl;

 // c) If enabled, access particle weights from external file:
 TFile *fileWithWeights = NULL;
 TList *listWithWeights = NULL;
 if(usePhiWeights||usePtWeights||useEtaWeights)
 {
  fileWithWeights = TFile::Open("weights.root","READ");
  if(fileWithWeights)
  {
   listWithWeights = (TList*)fileWithWeights->Get("weights");
  }
  else
  {
   cout << " WARNING: the file <weights.root> with weights from the previous run was not found."<<endl;
   break;
  }
 } // end of if(usePhiWeights||usePtWeights||useEtaWeights)

 // d) Configure the flow analysis methods:
 gROOT->LoadMacro("AliFlowAnalysisSPZDC.cxx++");
 AliFlowAnalysisSPZDC *qcPOI = NULL;
 // AliFlowAnalysisCRC *qcPOI = NULL;
 // gROOT->LoadMacro("MyAliFlowAnalysisWithMCEventPlane.cxx+");
 AliFlowAnalysisWithMCEventPlane *mcep = NULL;
 // MCEP = monte carlo event plane
 if(MCEP)
 {
  mcep = new AliFlowAnalysisWithMCEventPlane();
  mcep->SetHarmonic(iHarmonic); // default is v2
  mcep->Init();
 } // end of if(MCEP)
 // QC = Q-cumulants
 if(QCwPOI)
 {
   qcPOI = new AliFlowAnalysisSPZDC("AliFlowAnalysisSPZDC",nCen,wCen);
  //  qcPOI = new AliFlowAnalysisCRC("AliFlowAnalysisCRC",nCen,wCen);
   if(listWithWeights){qcPOI->SetWeightsList(listWithWeights);}
   if(usePhiWeights){qcPOI->SetUsePhiWeights(usePhiWeights);}
   if(usePtWeights){qcPOI->SetUsePtWeights(usePtWeights);}
   if(useEtaWeights){qcPOI->SetUseEtaWeights(useEtaWeights);}
   qcPOI->SetHarmonic(iHarmonic);
   qcPOI->SetCalculateDiffFlow(kFALSE);
   qcPOI->SetCalculate2DDiffFlow(kFALSE); // vs (pt,eta)
   qcPOI->SetApplyCorrectionForNUA(bNUACorr);
   qcPOI->SetFillMultipleControlHistograms(kFALSE);
   qcPOI->SetMultiplicityWeight("combinations"); // default (other supported options are "unit" and "multiplicity")
   qcPOI->SetCalculateCumulantsVsM(kFALSE);
   qcPOI->SetCalculateAllCorrelationsVsM(kFALSE); // calculate all correlations in mixed harmonics "vs M"
   qcPOI->SetnBinsMult(10000);
   qcPOI->SetMinMult(0);
   qcPOI->SetMaxMult(10000);
   qcPOI->SetBookOnlyBasicCCH(kTRUE); // book only basic common control histograms
   qcPOI->SetCalculateDiffFlowVsEta(kFALSE); // if you set kFALSE only differential flow vs pt is calculated
   qcPOI->SetCalculateMixedHarmonics(kFALSE); // calculate all multi-partice mixed-harmonics correlators
   qcPOI->SetStoreVarious(kTRUE);
   qcPOI->SetDataSet(AliFlowAnalysisSPZDC::k2010);
   qcPOI->SetCalculateCRC(kTRUE);
   qcPOI->SetCalculateCRCPt(kFALSE);
   qcPOI->SetCalculateCME(kFALSE);
   qcPOI->SetCalculateCRC2(bCalculateCRC2);
   qcPOI->SetCalculateFlowQC(bCalculateFlow);
   qcPOI->SetCalculateFlowZDC(bCalculateFlow);
   qcPOI->SetCalculateFlowVZ(kFALSE);
   qcPOI->SetUseVZERO(bUseVZERO);
   qcPOI->SetCRCEtaRange(-0.8,0.8);
   qcPOI->SetCRC2nEtaBins(nCRC2EtaBins);
   qcPOI->SetNUAforCRC(bNUACorr);
   qcPOI->SetUseZDC(bUseZDC);
   qcPOI->SetDivSigma(kFALSE);
   qcPOI->SetRecenterZDC(bRecenderZDC);
   qcPOI->SetQAZDCCuts(kTRUE);
   qcPOI->SetMinMulZN(0);
   qcPOI->SetMaxDevZN(100.);
   qcPOI->SetCalculateFlowQC(kTRUE);
   qcPOI->SetFlowQCCenBin(100);
   qcPOI->SetUseCRCRecenter(bUseCRCRecentering);
   qcPOI->SetCorrWeightTPC(AliFlowAnalysisSPZDC::kMultiplicity);
   qcPOI->SetCorrWeightZDC(AliFlowAnalysisSPZDC::kMultiplicity);
   qcPOI->SetCorrWeightVZ(AliFlowAnalysisSPZDC::kMultiplicity);
   qcPOI->SetInvertZDC(kFALSE);
  if(bUseCRCRecentering || bRecenderZDC) {
   TFile* QVecWeightsFile = TFile::Open(QVecWeightsFileName,"READ");
   if(!QVecWeightsFile) {
    cout << "ERROR: QVecWeightsFile not found!" << endl;
    exit(1);
   }
   TDirectory* QCFlowDir = dynamic_cast<TDirectory*>(QVecWeightsFile->FindObjectAny("output_QC"));
   TList* QCFlowList = dynamic_cast<TList*>(QCFlowDir->FindObjectAny("cobjQC"));
   TList* QVecWeightsList = dynamic_cast<TList*>(QCFlowList->FindObject("Q Vectors"));
   if(QVecWeightsList) {
    qcPOI->SetCRCQVecWeightsList(QVecWeightsList);
    cout << "Q Vector weights set (from " <<  QVecWeightsFileName.Data() << ")" << endl;
   }
   else{
    cout << "ERROR: QVecWeightsList not found!" << endl;
    exit(1);
   }
  }
  qcPOI->Init();
 } // end of if(qcPOI)

 // e) Simple cuts for RPs:
 AliFlowTrackSimpleCuts *cutsRP = new AliFlowTrackSimpleCuts();
 cutsRP->SetPtMax(ptMaxRP);
 cutsRP->SetPtMin(ptMinRP);
 cutsRP->SetEtaMax(etaMaxRP);
 cutsRP->SetEtaMin(etaMinRP);
 cutsRP->SetPhiMax(phiMaxRP*TMath::Pi()/180.);
 cutsRP->SetPhiMin(phiMinRP*TMath::Pi()/180.);
 if(bUseChargeRP){cutsRP->SetCharge(chargeRP);}

 // f) Simple cuts for POIs:
 AliFlowTrackSimpleCuts *cutsPOI = new AliFlowTrackSimpleCuts();
 cutsPOI->SetPtMax(ptMaxPOI);
 cutsPOI->SetPtMin(ptMinPOI);
 cutsPOI->SetEtaMax(etaMaxPOI);
 cutsPOI->SetEtaMin(etaMinPOI);
 cutsPOI->SetPhiMax(phiMaxPOI*TMath::Pi()/180.);
 cutsPOI->SetPhiMin(phiMinPOI*TMath::Pi()/180.);
 if(bUseChargePOI){cutsPOI->SetCharge(chargePOI);}

 Int_t runlist[] = {139510,139507};
 gRandom = new TRandom3(uiSeed); // if uiSeed is 0, the seed is determined uniquely in space and time via TUUID

 // g) Create and analyse events 'on the fly':
 for(Int_t i=0;i<iNevts;i++)
 {
  Int_t r = 0;
  if(bRunDepNUA) {
   r = gRandom->Integer(2);
   if(r==1) {
    eventMakerOnTheFly->SetUniformAcceptance(kFALSE);
    eventMakerOnTheFly->SetFirstSectorPhiMin(phiMin1);
    eventMakerOnTheFly->SetFirstSectorPhiMax(phiMax1);
    eventMakerOnTheFly->SetFirstSectorProbability(p1);
    eventMakerOnTheFly->SetSecondSectorPhiMin(phiMin2);
    eventMakerOnTheFly->SetSecondSectorPhiMax(phiMax2);
    eventMakerOnTheFly->SetSecondSectorProbability(p2);
   } else {
    eventMakerOnTheFly->SetUniformAcceptance(kTRUE);
   }
  } else {
   eventMakerOnTheFly->SetUniformAcceptance(kTRUE);
  }
  // Creating the event 'on the fly':
  AliFlowEventSimple *event = eventMakerOnTheFly->CreateEventOnTheFly(cutsRP,cutsPOI);
  // Passing the created event to flow analysis methods:
  if(MCEP){mcep->Make(event);}
  if(QCwPOI) {
   qcPOI->SetRunNumber(runlist[r]);
   qcPOI->Make(event);
  }
  delete event;
 } // end of for(Int_t i=0;i<iNevts;i++)

 // h) Create the output file and directory structure for the final results of all methods:
 //TString outputFileName = folderName+"/AnalysisResults.root";
 TString outputFileName = "AnalysisResults.root";
 TFile *outputFile = new TFile(outputFileName.Data(),"RECREATE");
 TDirectoryFile *dirFileFinalMC = new TDirectoryFile("output_MCEP","output_MCEP");
 TDirectoryFile *dirFileFinalQC = new TDirectoryFile("output_QC","output_QC");

 // i) Calculate and store the final results of all methods:
 if(MCEP) {
  mcep->Finish();
  mcep->WriteHistograms(dirFileFinalMC);
 }
 if(QCwPOI) {
  qcPOI->Finish();
  qcPOI->WriteHistograms(dirFileFinalQC);
 }

 outputFile->Close();
 delete outputFile;

 cout<<endl;
 cout<<endl;
 cout<<" ---- LANDED SUCCESSFULLY ---- "<<endl;
 cout<<endl;

 if(bFuckEP) {
  cout << "ZDC towers multiplicity map:" << endl;
  TProfile* ZDCMult = dynamic_cast<TProfile*>(eventMakerOnTheFly->GetZDCMult());
  for(Int_t i=0; i<8; i++) {
   cout << ZDCMult->GetBinContent(i+1) << endl;
  }
  cout<<endl;
 }

 timer.Stop();
 cout << endl;
 timer.Print();

} // end of int runFlowAnalysisOnTheFly(Int_t mode=mLocal)

void SetupPar(char* pararchivename)
{
 //Load par files, create analysis libraries
 //For testing, if par file already decompressed and modified
 //classes then do not decompress.

 TString cdir(Form("%s", gSystem->WorkingDirectory() )) ;
 TString parpar(Form("%s.par", pararchivename)) ;
 if ( gSystem->AccessPathName(parpar.Data()) ) {
  gSystem->ChangeDirectory(gSystem->Getenv("ALICE_ROOT")) ;
  TString processline(Form(".! make %s", parpar.Data())) ;
  gROOT->ProcessLine(processline.Data()) ;
  gSystem->ChangeDirectory(cdir) ;
  processline = Form(".! mv /tmp/%s .", parpar.Data()) ;
  gROOT->ProcessLine(processline.Data()) ;
 }
 if ( gSystem->AccessPathName(pararchivename) ) {
  TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
  gROOT->ProcessLine(processline.Data());
 }

 TString ocwd = gSystem->WorkingDirectory();
 gSystem->ChangeDirectory(pararchivename);

 // check for BUILD.sh and execute
 if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
  printf("*******************************\n");
  printf("*** Building PAR archive    ***\n");
  cout<<pararchivename<<endl;
  printf("*******************************\n");

  if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
   Error("runProcess","Cannot Build the PAR Archive! - Abort!");
   return -1;
  }
 }
 // check for SETUP.C and execute
 if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
  printf("*******************************\n");
  printf("*** Setup PAR archive       ***\n");
  cout<<pararchivename<<endl;
  printf("*******************************\n");
  gROOT->Macro("PROOF-INF/SETUP.C");
 }

 gSystem->ChangeDirectory(ocwd.Data());
 printf("Current dir: %s\n", ocwd.Data());
}

void CheckUserSettings()
{
 // Check if user settings make sense before taking off.

 if(iNevts <= 0)
 {
  printf("\n WARNING: nEvts <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult < 0.)
 {
  printf("\n WARNING: iMinMult < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMaxMult <= 0.)
 {
  printf("\n WARNING: iMaxMult <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(iMinMult >= iMaxMult)
 {
  printf("\n WARNING: iMinMult >= iMaxMult !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dMass < 0.)
 {
  printf("\n WARNING: dMass < 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(dTemperature <= 1e-44)
 {
  printf("\n WARNING: dTemperature <= 0 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV1) > 0.5)
 {
  printf("\n WARNING: |dV1| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV2) > 0.5)
 {
  printf("\n WARNING: |dV2| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV3) > 0.5)
 {
  printf("\n WARNING: |dV3| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV4) > 0.5)
 {
  printf("\n WARNING: |dV4| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV5) > 0.5)
 {
  printf("\n WARNING: |dV5| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }
 if(TMath::Abs(dV6) > 0.5)
 {
  printf("\n WARNING: |dV6| > 0.5 !!!! Please check your settings before taking off.\n\n");
  exit(0);
 }

 if(!uniformAcceptance && phiMin1 > phiMax1)
 {
  cout<<" WARNING: You must have phiMin1 < phiMax1 !!!!"<<endl;
  exit(0);
 }
 if(!uniformAcceptance && !((TMath::Abs(phiMin2) < 1.e-44) && (TMath::Abs(phiMax2) < 1.e-44) && (TMath::Abs(p2) < 1.e-44))
    && (phiMin2 < phiMax1 || phiMin2 > phiMax2))
 {
  cout<<" WARNING: You must have phiMin2 > phiMax1 and phiMin2 < phiMax2 !!!!"<<endl;
  exit(0);
 }
 if((phiMin1 < 0 || phiMin1 > 360) || (phiMax1 < 0 || phiMax1 > 360) ||
    (phiMin2 < 0 || phiMin2 > 360) || (phiMax2 < 0 || phiMax2 > 360) )
 {
  cout<<" WARNING: You must take azimuthal angles from interval [0,360] !!!!"<<endl;
  exit(0);
 }
 if((p1 < 0 || p1 > 1) || (p2 < 0 || p2 > 1))
 {
  cout<<" WARNING: you must take p1 and p2 from interval [0,1] !!!!"<<endl;
  exit(0);
 }
 if(bPtDependentV2 && bUniformFluctuationsV2)
 {
  cout<<" WARNING: Uniform fluctuations not supported for pt denependent v2 !!!!"<<endl;
  exit(0);
 }

} // end of void CheckUserSettings()

void LoadLibraries(const anaModes mode) {

 //--------------------------------------
 // Load the needed libraries most of them already loaded by aliroot
 //--------------------------------------
 //gSystem->Load("libTree");
 gSystem->Load("libGeom");
 gSystem->Load("libVMC");
 gSystem->Load("libXMLIO");
 gSystem->Load("libPhysics");
 gSystem->Load("libCore.so");
 gSystem->Load("libTree.so");

 //----------------------------------------------------------
 // >>>>>>>>>>> Local mode <<<<<<<<<<<<<<
 //----------------------------------------------------------
 if (mode==mLocal) {
  //--------------------------------------------------------
  // If you want to use already compiled libraries
  // in the aliroot distribution
  //--------------------------------------------------------
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libOADB.so");
  gSystem->Load("libCORRFW");
  cerr<<"libCORRFW loaded..."<<endl;
  gSystem->Load("libPWGflowBase");
  cerr<<"libPWGflowBase loaded..."<<endl;
  gSystem->Load("libPWGflowTasks");
  cerr<<"libPWGflowTasks loaded..."<<endl;
 }

 else if (mode == mLocalPAR) {
  //--------------------------------------------------------
  //If you want to use root and par files from aliroot
  //--------------------------------------------------------
  //If you want to use root and par files from aliroot
  //--------------------------------------------------------
  SetupPar("STEERBase");
  SetupPar("ESD");
  SetupPar("AOD");
  SetupPar("ANALYSIS");
  SetupPar("ANALYSISalice");
  SetupPar("PWG2AOD");
  SetupPar("CORRFW");
  SetupPar("PWGflowBase");
  cerr<<"PWGflowBase.par loaded..."<<endl;
  SetupPar("PWGflowTasks");
  cerr<<"PWGflowTasks.par loaded..."<<endl;
 }

 //---------------------------------------------------------
 // <<<<<<<<<< Source mode >>>>>>>>>>>>
 //---------------------------------------------------------
 else if (mode==mLocalSource) {

  // In root inline compile

  // Constants
  gROOT->LoadMacro("Base/AliFlowCommonConstants.cxx+");
  gROOT->LoadMacro("Base/AliFlowLYZConstants.cxx+");

  // Flow event
  gROOT->LoadMacro("Base/AliFlowVector.cxx+");
  gROOT->LoadMacro("Base/AliFlowTrackSimple.cxx+");
  gROOT->LoadMacro("Base/AliFlowTrackSimpleCuts.cxx+");
  gROOT->LoadMacro("Base/AliFlowEventSimple.cxx+");

  // Output histosgrams
  gROOT->LoadMacro("Base/AliFlowCommonHist.cxx+");
  gROOT->LoadMacro("Base/AliFlowCommonHistResults.cxx+");
  gROOT->LoadMacro("Base/AliFlowLYZHist1.cxx+");
  gROOT->LoadMacro("Base/AliFlowLYZHist2.cxx+");

  // Functions needed for various methods
  gROOT->LoadMacro("Base/AliCumulantsFunctions.cxx+");
  gROOT->LoadMacro("Base/AliFlowLYZEventPlane.cxx+");

  // Flow Analysis code for various methods
  gROOT->LoadMacro("Base/AliFlowAnalysisWithMCEventPlane.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithScalarProduct.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithLYZEventPlane.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithLeeYangZeros.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithCumulants.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithQCumulants.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithFittingQDistribution.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithMixedHarmonics.cxx+");
  gROOT->LoadMacro("Base/AliFlowAnalysisWithNestedLoops.cxx+");

  // Class to fill the FlowEvent on the fly (generate Monte Carlo events)
  //    gROOT->LoadMacro("Base/AliFlowEventSimpleMakerOnTheFly.cxx+");

  cout << "finished loading macros!" << endl;

 }

 // Use AliRoot includes to compile our task
 gROOT->ProcessLine(".include $ALICE_ROOT/include");
 gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
 gROOT->LoadMacro("MyAliFlowEventSimpleMakerOnTheFly.cxx++g");

}
