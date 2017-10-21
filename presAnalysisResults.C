
// Set number of tasks:
const Int_t nTasks = 1;
const Int_t nQCum = 2;
// Set canvas division:
const Int_t nRows = nQCum;
const Int_t nColumns = nTasks;
const Int_t nPtbins = 6;
// Set name of the output file of flow analysis to be accessed:
TString folder = "output/testBASEwEG";
TString outputFileName = folder+"/correctedAnalysisResults.root"; //correctedAnalysisResults
TString taskSuffix[nTasks] = {"QCCEA_POIPOI"};//,AODFB), Form("QCCEA_NUA",AODFB)}; //,"QCCEA_PW_MC"};
const Int_t nTasksBegin = 0;
const Int_t nTasksEnd = 1;

const Int_t nPtIntervals = 3;
Double_t ptInterval[nPtIntervals+1] = {0.,2.,5.,10.}; // in GeV
Int_t nMergedBins[nPtIntervals] = {1,2,5}; // for instance in 2nd pt interval (2-5 GeV) 5 pt bins will be merged into 1
// Set here if you want to use default values for harmonic, pt and eta binning:
Bool_t useDefaultValues = kTRUE;
TString methodForSettings = "QC"; // alternatively set here method from whose output files harmonic, pt and eta binning will be accessed
// Set here if you want to show error mesh:
Bool_t showErrorMesh = kTRUE;
Bool_t showErrorMeshDiffFlow = kTRUE;
// Set here if both the error mesh and markers will be shown for specified method in the plots for differential flow:
Bool_t showBothErrorMeshAndMarkers = kFALSE;

// Do not touch this unless you are looking for a trouble:
const Int_t nMethods = 12;
TString method[nMethods] = {"MCEP","SP","GFC","QC","FQD","LYZ1SUM","LYZ1PROD","LYZ2SUM","LYZ2PROD","LYZEP","MH","NL"};
TList *list[nMethods] = {NULL}; // lists holding histograms for each flow analysis method

enum libModes{mLocal,mLocalSource};
AliFlowAnalysisCRC* qc;

void presAnalysisResults(TString analysisType="AOD",Int_t analysisMode=mLocal)
{
    // Load needed libraries:
    LoadLibrariesCFR();
    GlobalSettings();
    
    // plotZDCv1D0(5,kFALSE);
    // plotZDCv1D0(6,kTRUE);
    plotZDCv1D0(8,kTRUE);
    
    // plotZDCv1D0Compare(6,kTRUE);
    // plotZDCv1D0CompareRes(5,kTRUE);
    
    // plotZDCv1Ch(5);
    // plotZDCv1Ch(6);
    // plotZDCv1Ch(7);
    
} // end of void plotAnalysisResults(TString analysisType="AOD",Int_t analysisMode=mLocal)

// ===========================================================================================

void plotZDCv1D0(Int_t ord=6, Bool_t bSave = kTRUE, Bool_t bSetErrorFromSignificance = kFALSE)
{
    TFile *mergedFile = TFile::Open(Form("AnalysisResultsD0_%d.root",ord),"READ");
    // TFile *mergedFile = TFile::Open("AnalysisResults.root","READ");
    TList *FlowList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Flow SP ZDC"));
    
    TList *VariousList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Various"));
    TH1D* EventCounter = (TH1D*)VariousList->FindObject("EventCounter");
    Double_t NDs = EventCounter->GetBinContent(1)*1000.;
    
    TProfile* flowpro[4] = {NULL};
    TH1D* flowhist[4] = {NULL};
    for(Int_t i=0; i<4; i++) {
        flowpro[i] = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][0][%d]",3+i));
        flowpro[i]->GetXaxis()->SetRange(1,flowpro[i]->GetNbinsX());
        flowhist[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),flowpro[i]->GetNbinsX(),-0.8,0.8);
        // TCanvas* canvas = new TCanvas(Form("Cum_%d",i),Form("Cum_%d",i),1000,1000);
        // flowpro[i]->Draw();
    }
    
    TProfile* resA = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][0]"));
    TProfile* resC = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][1]"));
    TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
    Double_t QAQB = resA->GetBinContent(3)*resC->GetBinContent(3);
    Double_t QAQBer = res->GetBinError(3);
    Double_t den = sqrt(fabs(QAQB));
    printf("resolution: %e %e --> %e \n",resA->GetBinContent(3),resC->GetBinContent(3),den);
    
    TH1D* flowaver = (TH1D*)flowhist[0]->Clone("aver");
    flowaver->Rebin(5);
    Double_t averval = flowaver->GetBinContent(1);
    
    //positive particles
    for(Int_t eb=0; eb<flowpro[0]->GetNbinsX(); eb++) {
        Double_t v1P = flowpro[0]->GetBinContent(eb+1);
        Double_t v1Per = flowpro[0]->GetBinError(eb+1);
        Double_t v1T = flowpro[1]->GetBinContent(eb+1);
        Double_t v1Ter = flowpro[1]->GetBinError(eb+1);
        
        if(bSetErrorFromSignificance) {
            // set stat errors equal to 1/significance
            // significance/sqrt(events) = 5E-3 (from AliceITSupgrade) --> significance = 5E-3 * sqrt(events) = 5E-3 * sqrt(NDs / (NDs/event))
            // how did we compute NDs/event?
            // 2.3 D0 produced in |y|<0.5 per MB event (from AliceITSupgrade)
            // reconstruction efficiency 0.05 (at <pt>) (from AliceITSupgrade)
            // BR = 3.8%--> 2.3*0.05*0.038 = 0.00437   D0 / MB event
            sign = 5.E-3 * sqrt(1.E9);
            v1Per = v1P*sqrt(pow(v1Per/v1P,2.) + pow(1./sign,2.));
            v1Ter = v1T*sqrt(pow(v1Ter/v1T,2.) + pow(1./sign,2.));
        }
        
        Double_t v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
        Double_t v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
        Double_t v1odder = 0.5*sqrt(v1oddSq);
        flowhist[0]->SetBinContent(eb+1,v1odd);
        flowhist[0]->SetBinError(eb+1,v1odder);
        
        //negative particles
        v1P = flowpro[2]->GetBinContent(eb+1);
        v1Per = flowpro[2]->GetBinError(eb+1);
        v1T = flowpro[3]->GetBinContent(eb+1);
        v1Ter = flowpro[3]->GetBinError(eb+1);
        
        if(bSetErrorFromSignificance) {
            sign = 5.E-3 * sqrt(1.E9);
            v1Per = v1P*sqrt(pow(v1Per/v1P,2.) + pow(1./sign,2.));
            v1Ter = v1T*sqrt(pow(v1Ter/v1T,2.) + pow(1./sign,2.));
        }
        
        v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
        v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
        v1odder = 0.5*sqrt(v1oddSq);
        flowhist[1]->SetBinContent(eb+1,v1odd);
        flowhist[1]->SetBinError(eb+1,v1odder);
    }
    
    
    
    TH1F *systDo = (TH1F*)flowhist[0]->Clone();
    TH1F *systDobar = (TH1F*)flowhist[0]->Clone();
    
    TH1F *relsysterrMinMax = (TH1F*)flowhist[0]->Clone();
    relsysterrMinMax->Reset();
    TGraphErrors* systerrorInclv1 = new TGraphErrors(5);
    systerrorInclv1->SetFillStyle(0);
    systerrorInclv1->SetLineColor(kBlue);
    TGraphErrors* systerrorInclv2 = new TGraphErrors(5);
    systerrorInclv2->SetFillStyle(0);
    systerrorInclv2->SetLineColor(kRed);
    Double_t centrality[20];
    Double_t NUA[20];
    Double_t ZDC[20];
    Double_t magnetpol[20];
    Double_t pileup[20];
    Double_t signal[20];
    Double_t Bfeed[20];
    Double_t gey[20];
    
    TFile *systPhotonic = new TFile("AbsoluteSystv1odd.root");
    TH1D *syst2 = (TH1D*)systPhotonic->Get("v1+");
    
    
    
    
    for(int j=0;j<5;j++){
        //        centrality[j] = (2*flowhist[0]->GetBinContent(j+1))/100;
        //        NUA[j] = (1*flowhist[0]->GetBinContent(j+1))/100;
        //        ZDC[j] = (7*flowhist[0]->GetBinContent(j+1))/100;
        //        magnetpol[j] = (4*flowhist[0]->GetBinContent(j+1))/100;
        //        pileup[j] = (2*flowhist[0]->GetBinContent(j+1))/100;
        centrality[j] = syst2->GetBinContent(j+1)*25;
        signal[j] = (10*flowhist[0]->GetBinContent(j+1))/100;
        Bfeed[j] = (48*flowhist[0]->GetBinContent(j+1))/100;
        
        gey[j] = TMath::Sqrt((centrality[j]*centrality[j])+(signal[j]*signal[j]));
        systerrorInclv2->SetPoint(j,flowhist[0]->GetBinCenter(j+1),flowhist[0]->GetBinContent(j+1));
        systerrorInclv2->SetPointError(j,flowhist[0]->GetBinWidth(j+1)/4,gey[j]);
        systDo->SetBinContent(j+1,flowhist[0]->GetBinContent(j+1));
        systDo->SetBinError(j+1,gey[j]);
        
    }
    
    Double_t centrality1[20];
    Double_t NUA1[20];
    Double_t ZDC1[20];
    Double_t magnetpol1[20];
    Double_t pileup1[20];
    Double_t signal1[20];
    Double_t Bfeed1[20];
    
    
    Double_t gey1[20];
    for(int j=0;j<5;j++){
        //        centrality[j] = (2*flowhist[0]->GetBinContent(j+1))/100;
        //        NUA[j] = (1*flowhist[0]->GetBinContent(j+1))/100;
        //        ZDC[j] = (7*flowhist[0]->GetBinContent(j+1))/100;
        //        magnetpol[j] = (4*flowhist[0]->GetBinContent(j+1))/100;
        //        pileup[j] = (2*flowhist[0]->GetBinContent(j+1))/100;
        centrality1[j] = syst2->GetBinContent(j+1)*25;
        signal1[j] = (10*flowhist[0]->GetBinContent(j+1))/100;
        Bfeed1[j] = (48*flowhist[0]->GetBinContent(j+1))/100;
        
        gey1[j] = TMath::Sqrt((centrality1[j]*centrality1[j])+(signal1[j]*signal1[j]));
        systerrorInclv1->SetPoint(j,flowhist[1]->GetBinCenter(j+1),flowhist[1]->GetBinContent(j+1));
        systerrorInclv1->SetPointError(j,flowhist[1]->GetBinWidth(j+1)/4,gey1[j]);
        systDobar->SetBinContent(j+1,flowhist[1]->GetBinContent(j+1));
        systDobar->SetBinError(j+1,gey1[j]);
    }
    
    
    TCanvas *u = new TCanvas();
    relsysterrMinMax->SetTitle();
    relsysterrMinMax->SetLineColor(kBlack);
    relsysterrMinMax->GetXaxis()->SetTitle("p_{T}");
    relsysterrMinMax->GetYaxis()->SetTitle("relative error (%)");
    relsysterrMinMax->Draw();
    
    
    
    TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
    // res->Draw();
    
    TCanvas* canvas = new TCanvas(Form("v1odd_%.2e",NDs),Form("v1odd_%.2e",NDs),1000,1000);
    SetColorHisto(flowhist[0],kRed);
    flowhist[0]->GetYaxis()->SetRangeUser(-8.E-2,10.E-2);
    ALICESettings(flowhist[0],"#eta","v_{1}^{odd}");
 //   flowhist[0]->Scale(-1.);
    flowhist[0]->Draw();
    SetColorHisto(flowhist[1],kBlue);
    flowhist[1]->Draw("same");
    systerrorInclv2->Draw("2same");
    systerrorInclv1->Draw("2same");
    
    TLine *line3 = new TLine(-0.8,0.,0.8,0.);
    line3->SetLineStyle(2);
    line3->SetLineColor(1);
    line3->Draw("same");
    
    
    
    
    TLegend *lg = new TLegend(0.2,0.65,0.5,0.85);
    lg->SetBorderSize(0.1); lg->SetFillColor(0);
    lg->SetHeader(Form("%.2e D^{0}, #bar{D}^{0}",NDs));
    lg->AddEntry(flowhist[0],"D^{0}","lp");
    lg->AddEntry(flowhist[1],"#bar{D}^{0}","lp");
    lg->Draw();
    TPaveText *pt = new TPaveText(0.55,0.81,0.85,0.85,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("#Deltav_{1}^{odd} = 3e-02");
    pt->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/v1Ds_stat%d.png",ord));
    
    TCanvas* canvas = new TCanvas(Form("Deltav1odd_%.2e",NDs),Form("Deltav1odd_%.2e",NDs),1000,1000);
    TH1D* flowdiff = (TH1D*)flowhist[0]->Clone("diff");
    SetColorHisto(flowdiff,kBlack);
    flowdiff->Add(flowhist[1],-1.);
    
    TH1F* flowdiffsyst = (TH1F*)systDo->Clone("diffsyst");
    flowdiffsyst->Add(systDobar,-1.);
    
    
        TGraphErrors* systerrorInclv1delta = new TGraphErrors(5);
        systerrorInclv1delta->SetFillStyle(0);
        systerrorInclv1delta->SetLineColor(kBlack);
    //    TH1D *syst2 = (TH1D*)systPhotonic->Get("#Deltav1");
    //    Double_t chpart[20];
        for(int j=0;j<5;j++){
//            chpart[j] = (syst2->GetBinContent(j+1)*flowdiff->GetBinContent(j+1))/100;
//            gey[j] = TMath::Sqrt((chpart[j]*chpart[j]));
            systerrorInclv1delta->SetPoint(j,flowdiffsyst->GetBinCenter(j+1),flowdiffsyst->GetBinContent(j+1));
            systerrorInclv1delta->SetPointError(j,flowdiffsyst->GetBinWidth(j+1)/4,flowdiffsyst->GetBinError(j+1));
    
        }
    
    
    
    flowdiff->GetYaxis()->SetRangeUser(-8.E-2,10.E-2);
    flowdiff->Draw();
    systerrorInclv1delta->Draw("2same");
    TF1* fitpol = new TF1("fitpol","[0]*x",-0.8,0.8);
    flowdiff->Fit(fitpol,"QN");
    Double_t kappa = fitpol->GetParameter(0);
    Double_t kappaer = fitpol->GetParError(0);
    fitpol->Draw("same");
    
    TLegend *lg = new TLegend(0.2,0.7,0.5,0.85);
    lg->SetBorderSize(0.1); lg->SetFillColor(0);
    lg->SetHeader(Form("%.2e D^{0}, #bar{D}^{0}",NDs));
    lg->AddEntry(flowdiff,"v_{1}^{odd}[D^{0}] - v_{1}^{odd}[#bar{D}^{0}]","lp");
    lg->Draw();
    TPaveText *pt = new TPaveText(0.55,0.81,0.85,0.85,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("#Deltav_{1}^{odd} = 3e-02");
    pt->Draw("same");
    
    TPaveText *pt = new TPaveText(0.55,0.67,0.85,0.79,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("fit function: k #times #eta");
    pt->AddText(Form("k = %.1e #pm %.1e",kappa,kappaer));
    pt->Draw("same");
    
    line3->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/Deltav1Ds_stat%d.png",ord));
    
    // TFile *outputFile = new TFile("v1D0_Dubla.root","RECREATE");
    // flowdiff->SetName("Deltav1_100x");
    // outputFile->Add(flowdiff);
    // flowhist[0]->SetName("v1plus_100x");
    // outputFile->Add(flowhist[0]);
    // flowhist[1]->SetName("v1minus_100x");
    // outputFile->Add(flowhist[1]);
    // outputFile->Write();
    // outputFile->Close();
    
}
// ===========================================================================================

void plotZDCv1D0Compare(Int_t ord=6, Bool_t bSave = kTRUE)
{
    TFile *mergedFile = TFile::Open(Form("AnalysisResultsD0_%d.root",ord),"READ");
    // TFile *mergedFile = TFile::Open("AnalysisResults.root","READ");
    TList *FlowList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Flow SP ZDC"));
    
    TList *VariousList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Various"));
    TH1D* EventCounter = (TH1D*)VariousList->FindObject("EventCounter");
    Double_t NDs = EventCounter->GetBinContent(1)*1000.;
    
    TProfile* flowpro[4] = {NULL};
    TH1D* flowhist[4] = {NULL};
    for(Int_t i=0; i<4; i++) {
        flowpro[i] = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][0][%d]",3+i));
        flowpro[i]->GetXaxis()->SetRange(1,flowpro[i]->GetNbinsX());
        flowhist[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),flowpro[i]->GetNbinsX(),-0.8,0.8);
        // TCanvas* canvas = new TCanvas(Form("Cum_%d",i),Form("Cum_%d",i),1000,1000);
        // flowpro[i]->Draw();
    }
    
    TProfile* resA = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][0]"));
    TProfile* resC = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][1]"));
    TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
    Double_t QAQB = resA->GetBinContent(3)*resC->GetBinContent(3);
    Double_t QAQBer = res->GetBinError(3);
    Double_t den = sqrt(fabs(QAQB));
    printf("resolution: %e %e --> %e \n",resA->GetBinContent(3),resC->GetBinContent(3),den);
    
    TH1D* flowaver = (TH1D*)flowhist[0]->Clone("aver");
    flowaver->Rebin(5);
    Double_t averval = flowaver->GetBinContent(1);
    TCanvas* canvas = new TCanvas(Form("Deltav1odd_%.2e",NDs),Form("Deltav1odd_%.2e",NDs),1000,1000);
    TLegend *lg = new TLegend(0.2,0.7,0.5,0.85);
    lg->SetBorderSize(0.1); lg->SetFillColor(0);
    
    //positive particles
    for(Int_t k=0; k<2; k++) {
        for(Int_t eb=0; eb<flowpro[0]->GetNbinsX(); eb++) {
            Double_t v1P = flowpro[0]->GetBinContent(eb+1);
            Double_t v1Per = flowpro[0]->GetBinError(eb+1);
            Double_t v1T = flowpro[1]->GetBinContent(eb+1);
            Double_t v1Ter = flowpro[1]->GetBinError(eb+1);
            
            if(k==1) {
                // set stat errors equal to 1/significance
                // significance/sqrt(events) = 5E-3 (from AliceITSupgrade) --> significance = 5E-3 * sqrt(events) = 5E-3 * sqrt(NDs / (NDs/event))
                // how did we compute NDs/event?
                // 2.3 D0 produced in |y|<0.5 per MB event (from AliceITSupgrade)
                // reconstruction efficiency 0.05 (at <pt>) (from AliceITSupgrade)
                // BR = 3.8%--> 2.3*0.05*0.038 = 0.00437   D0 / MB event
                sign = 5.E-3 * sqrt(1.E9);
                v1Per = v1P*sqrt(pow(v1Per/v1P,2.) + pow(1./sign,2.));
                v1Ter = v1T*sqrt(pow(v1Ter/v1T,2.) + pow(1./sign,2.));
            }
            
            Double_t v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
            Double_t v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
            Double_t v1odder = 0.5*sqrt(v1oddSq);
            flowhist[0]->SetBinContent(eb+1,v1odd);
            flowhist[0]->SetBinError(eb+1,v1odder);
            
            //negative particles
            v1P = flowpro[2]->GetBinContent(eb+1);
            v1Per = flowpro[2]->GetBinError(eb+1);
            v1T = flowpro[3]->GetBinContent(eb+1);
            v1Ter = flowpro[3]->GetBinError(eb+1);
            
            if(k==1) {
                sign = 5.E-3 * sqrt(1.E9);
                v1Per = v1P*sqrt(pow(v1Per/v1P,2.) + pow(1./sign,2.));
                v1Ter = v1T*sqrt(pow(v1Ter/v1T,2.) + pow(1./sign,2.));
            }
            
            v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
            v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
            v1odder = 0.5*sqrt(v1oddSq);
            flowhist[1]->SetBinContent(eb+1,v1odd);
            flowhist[1]->SetBinError(eb+1,v1odder);
        }
        
        TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
        // res->Draw();
        
        TH1D* flowdiff = (TH1D*)flowhist[0]->Clone("diff");
        SetColorHisto(flowdiff,kBlack);
        flowdiff->Add(flowhist[1],-1.);
        
        
        
        
        
        flowdiff->GetYaxis()->SetRangeUser(-8.E-2,10.E-2);
        if(k==0) flowdiff->SetLineColor(kRed);
        if(k==0) flowdiff->SetLineWidth(10);
        if(k==1) flowdiff->SetLineWidth(5);
        flowdiff->Draw("same");
        
        if(k==0) lg->AddEntry(flowdiff,"#sigma_{v1}","lp");
        if(k==1) lg->AddEntry(flowdiff,"v_{1} #sqrt{(#sigma_{v1}/v_{1})^{2}+(1/s)^{2}}","lp");
    }
    
    lg->SetHeader(Form("%.2e D^{0}, #bar{D}^{0}",NDs));
    lg->Draw();
    TPaveText *pt = new TPaveText(0.55,0.81,0.85,0.85,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("#Deltav_{1}^{odd} = 3e-02");
    pt->Draw("same");
    
    TLine *line3 = new TLine(-0.8,0.,0.8,0.);
    line3->SetLineStyle(2);
    line3->SetLineColor(1);
    line3->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/Deltav1Ds_err_stat%d.png",ord));
    
}
// ===========================================================================================

void plotZDCv1D0CompareRes(Int_t ord=6, Bool_t bSave = kTRUE)
{
    TCanvas* canvas = new TCanvas(Form("Deltav1odd"),Form("Deltav1odd"),1000,1000);
    TLegend *lg = new TLegend(0.2,0.6,0.65,0.85);
    lg->SetBorderSize(0.1); lg->SetFillColor(0); lg->SetFillStyle(0);
    
    Float_t resk[] = {0.1,0.15,1.};
    
    //positive particles
    for(Int_t k=0; k<3; k++) {
        
        TFile *mergedFile = TFile::Open(Form("AnalysisResultsD0_%d_res%d.root",ord,k),"READ");
        // TFile *mergedFile = TFile::Open("AnalysisResults.root","READ");
        TList *FlowList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Flow SP ZDC"));
        
        TList *VariousList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Various"));
        TH1D* EventCounter = (TH1D*)VariousList->FindObject("EventCounter");
        Double_t NDs = EventCounter->GetBinContent(1)*1000.;
        
        TProfile* flowpro[4] = {NULL};
        TH1D* flowhist[4] = {NULL};
        for(Int_t i=0; i<4; i++) {
            flowpro[i] = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][0][%d]",3+i));
            flowpro[i]->GetXaxis()->SetRange(1,flowpro[i]->GetNbinsX());
            flowhist[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),flowpro[i]->GetNbinsX(),-0.8,0.8);
            // TCanvas* canvas = new TCanvas(Form("Cum_%d",i),Form("Cum_%d",i),1000,1000);
            // flowpro[i]->Draw();
        }
        
        TProfile* resA = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][0]"));
        TProfile* resC = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][1]"));
        TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
        Double_t QAQB = resA->GetBinContent(3)*resC->GetBinContent(3);
        Double_t QAQBer = res->GetBinError(3);
        Double_t den = sqrt(fabs(QAQB));
        printf("resolution: %e %e --> %e \n",resA->GetBinContent(3),resC->GetBinContent(3),den);
        
        TH1D* flowaver = (TH1D*)flowhist[0]->Clone("aver");
        flowaver->Rebin(5);
        Double_t averval = flowaver->GetBinContent(1);
        
        for(Int_t eb=0; eb<flowpro[0]->GetNbinsX(); eb++) {
            Double_t v1P = flowpro[0]->GetBinContent(eb+1);
            Double_t v1Per = flowpro[0]->GetBinError(eb+1);
            Double_t v1T = flowpro[1]->GetBinContent(eb+1);
            Double_t v1Ter = flowpro[1]->GetBinError(eb+1);
            
            Double_t v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
            Double_t v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
            Double_t v1odder = 0.5*sqrt(v1oddSq);
            flowhist[0]->SetBinContent(eb+1,v1odd);
            flowhist[0]->SetBinError(eb+1,v1odder);
            
            //negative particles
            v1P = flowpro[2]->GetBinContent(eb+1);
            v1Per = flowpro[2]->GetBinError(eb+1);
            v1T = flowpro[3]->GetBinContent(eb+1);
            v1Ter = flowpro[3]->GetBinError(eb+1);
            
            v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
            v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
            v1odder = 0.5*sqrt(v1oddSq);
            flowhist[1]->SetBinContent(eb+1,v1odd);
            flowhist[1]->SetBinError(eb+1,v1odder);
        }
        
        TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
        // res->Draw();
        
        TH1D* flowdiff = (TH1D*)flowhist[0]->Clone("diff");
        SetColorHisto(flowdiff,k+1);
        flowdiff->Add(flowhist[1],-1.);
        flowdiff->GetYaxis()->SetRangeUser(-8.E-2,10.E-2);
        flowdiff->Draw("same");
        
        TF1* fitpol = new TF1("fitpol","[0]*x",-0.8,0.8);
        flowdiff->Fit(fitpol,"QN");
        Double_t kappa = fitpol->GetParameter(0);
        Double_t kappaer = fitpol->GetParError(0);
        fitpol->SetLineColor(k+1);
        fitpol->Draw("same");
        
        lg->SetHeader(Form("%.2e D^{0}, #bar{D}^{0}, fit function: k #times #eta",NDs));
        lg->AddEntry(flowdiff,Form("R_{ZDC} = %.2f, k = %.1e #pm %.1e",resk[k],kappa,kappaer),"lp");
    }
    
    lg->Draw();
    
    TLine *line3 = new TLine(-0.8,0.,0.8,0.);
    line3->SetLineStyle(2);
    line3->SetLineColor(1);
    line3->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/Deltav1Ds_res_stat%d.png",ord));
    
}
// ===========================================================================================

void plotZDCv1Ch(Int_t ord=6)
{
    TFile *mergedFile = TFile::Open(Form("AnalysisResultsCh_%d.root",ord),"READ");
    // TFile *mergedFile = TFile::Open("AnalysisResults.root","READ");
    TList *FlowList = dynamic_cast<TList*>(mergedFile->FindObjectAny("cobjQC")->FindObject("Flow SP ZDC"));
    Bool_t bSave = kTRUE;
    Double_t FakeScale = 0.0025;
    
    TProfile* flowpro[4] = {NULL};
    TH1D* flowhist[4] = {NULL};
    for(Int_t i=0; i<4; i++) {
        flowpro[i] = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][0][%d]",3+i));
        flowpro[i]->GetXaxis()->SetRange(1,flowpro[i]->GetNbinsX());
        flowpro[i]->Scale(FakeScale);
        flowhist[i] = new TH1D(Form("hist_%d",i),Form("hist_%d",i),flowpro[i]->GetNbinsX(),-0.8,0.8);
        // TCanvas* canvas = new TCanvas(Form("Cum_%d",i),Form("Cum_%d",i),1000,1000);
        // flowpro[i]->Draw();
    }
    
    TProfile* resA = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][0]"));
    TProfile* resC = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][1]"));
    TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
    Double_t QAQB = resA->GetBinContent(3)*resC->GetBinContent(3);
    Double_t QAQBer = res->GetBinError(3);
    Double_t den = sqrt(fabs(QAQB));
    printf("resolution: %e %e --> %e \n",resA->GetBinContent(3),resC->GetBinContent(3),den);
    
    TH1D* flowaver = (TH1D*)flowhist[0]->Clone("aver");
    flowaver->Rebin(5);
    Double_t averval = flowaver->GetBinContent(1);
    
    //positive particles
    for(Int_t eb=0; eb<flowpro[0]->GetNbinsX(); eb++) {
        Double_t v1P = flowpro[0]->GetBinContent(eb+1);
        Double_t v1Per = flowpro[0]->GetBinError(eb+1);
        Double_t v1T = flowpro[1]->GetBinContent(eb+1);
        Double_t v1Ter = flowpro[1]->GetBinError(eb+1);
        
        Double_t v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
        Double_t v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
        Double_t v1odder = 0.5*sqrt(v1oddSq);
        flowhist[0]->SetBinContent(eb+1,v1odd);
        flowhist[0]->SetBinError(eb+1,v1odder);
        
        //negative particles
        v1P = flowpro[2]->GetBinContent(eb+1);
        v1Per = flowpro[2]->GetBinError(eb+1);
        v1T = flowpro[3]->GetBinContent(eb+1);
        v1Ter = flowpro[3]->GetBinError(eb+1);
        
        v1odd = (v1P-v1T)/(sqrt(fabs(QAQB))*2.);
        v1oddSq = pow(v1Per/den,2.) + pow(v1Ter/den,2.) + pow(QAQBer*0.5*(averval)/pow(fabs(QAQB),1.5),2.);
        v1odder = 0.5*sqrt(v1oddSq);
        flowhist[1]->SetBinContent(eb+1,v1odd);
        flowhist[1]->SetBinError(eb+1,v1odder);
    }
    
    TProfile* res = (TProfile*)FlowList->FindObject(Form("fFlowSPZDCv1etaPro[0][2][2]"));
    res->Draw();
    
    TCanvas* canvas = new TCanvas("Cum","Cum",1000,1000);
    SetColorHisto(flowhist[0],kRed);
    flowhist[0]->GetYaxis()->SetRangeUser(-8.E-4,8.E-4);
    ALICESettings(flowhist[0],"#eta","v_{1}^{odd}");
    flowhist[0]->Draw();
    SetColorHisto(flowhist[1],kBlue);
    flowhist[1]->Draw("same");
    
    TLine *line3 = new TLine(-0.8,0.,0.8,0.);
    line3->SetLineStyle(2);
    line3->SetLineColor(1);
    line3->Draw("same");
    
    TLegend *lg = new TLegend(0.6,0.65,0.85,0.78);
    lg->SetBorderSize(0.1); lg->SetFillColor(0);
    lg->SetHeader(Form("4.8e%02d h^{#pm}",ord+5));
    lg->AddEntry(flowhist[0],"h^{+}","lp");
    lg->AddEntry(flowhist[1],"h^{-}","lp");
    lg->Draw();
    TPaveText *pt = new TPaveText(0.55,0.81,0.85,0.85,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("#Deltav_{1}^{odd} = 5e-05");
    pt->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/v1Ch_stat%d.png",ord));
    
    TCanvas* canvas = new TCanvas("CumDiff","CumDiff",1000,1000);
    TH1D* flowdiff = (TH1D*)flowhist[0]->Clone("diff");
    SetColorHisto(flowdiff,kBlack);
    flowdiff->Add(flowhist[1],-1.);
    flowdiff->GetYaxis()->SetRangeUser(-0.4E-3,0.4E-3);
    flowdiff->Draw();
    
    TF1* fitpol = new TF1("fitpol","[0]*x",-0.8,0.8);
    flowdiff->Fit(fitpol,"QN");
    Double_t kappa = fitpol->GetParameter(0);
    Double_t kappaer = fitpol->GetParError(0);
    fitpol->Draw("same");
    
    TLegend *lg = new TLegend(0.2,0.7,0.5,0.85);
    lg->SetBorderSize(0.1); lg->SetFillColor(0);
    lg->SetHeader(Form("4.8e%02d  h^{#pm}",ord+5));
    lg->AddEntry(flowdiff,"v_{1}^{odd}[h^{+}] - v_{1}^{odd}[h^{-}]","lp");
    lg->Draw();
    TPaveText *pt = new TPaveText(0.55,0.81,0.85,0.85,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("#Deltav_{1}^{odd} = 5e-05");
    pt->Draw("same");
    
    TPaveText *pt = new TPaveText(0.55,0.67,0.85,0.79,"brNDC");
    pt->SetBorderSize(0.1); pt->SetFillColor(0);
    pt->AddText("fit function: k #times #eta");
    pt->AddText(Form("k = %.1e #pm %.1e",kappa,kappaer));
    pt->Draw("same");
    
    line3->Draw("same");
    if(bSave) canvas->SaveAs(Form("RUN3Studyv1odd/Deltav1Ch_stat%d.png",ord));
    
}
// ===========================================================================================

void plotCME()
{
    TCanvas* canvas = new TCanvas("Cum","Cum",1000,800);
    Int_t const size = 4;
    canvas->Divide(2,2);
    
    canvas->cd(1);
    TH1D* CMEB = qc->GetCMEZDCCorrPro(2);
    CMEB->SetLineColor(kBlack);
    CMEB->Draw("E1same");
    TH1D* CMEP = qc->GetCMEZDCCorrPro(0);
    CMEP->SetLineColor(kRed);
    CMEP->Draw("E1same");
    TH1D* CMEN = qc->GetCMEZDCCorrPro(1);
    CMEN->SetLineColor(kBlue);
    CMEN->Draw("E1same");
    
    canvas->cd(2);
    TH1D* CMER = qc->GetCMEZDCCorrPro(4);
    CMER->SetLineColor(kRed);
    CMER->Draw("E1same");
    TH1D* CMED = qc->GetCMEZDCCorrPro(3);
    CMED->SetLineColor(kBlue);
    CMED->Draw("E1same");
    
    TLegend* legend = new TLegend(0.1,0.7,0.75,0.9);
    legend->AddEntry(CMER, Form("%f (%f)",CMER->GetMean(),CMER->GetMeanError()),"lp");
    legend->AddEntry(CMED, Form("%f (%f)",CMED->GetMean(),CMED->GetMeanError()),"lp");
    legend->Draw();
    
    canvas->cd(3);
    TH1D* CMERS = qc->GetCMEZDCCorrPro(6);
    CMERS->SetLineColor(kBlue);
    CMERS->Draw("E1same");
    TH1D* CMEDS = qc->GetCMEZDCCorrPro(5);
    CMEDS->SetLineColor(kRed);
    CMEDS->Draw("E1same");
    
    TLegend* legend = new TLegend(0.1,0.7,0.75,0.9);
    legend->AddEntry(CMERS, Form("%f (%f)",CMERS->GetMean(),CMERS->GetMeanError()),"lp");
    legend->AddEntry(CMEDS, Form("%f (%f)",CMEDS->GetMean(),CMEDS->GetMeanError()),"lp");
    legend->Draw();
    
    canvas->cd(4);
    TH1D* CMERS = qc->GetCMEZDCCorrPro(8);
    CMERS->SetLineColor(kBlack);
    CMERS->Draw("E1same");
    TH1D* CMEDS = qc->GetCMEZDCCorrPro(7);
    CMEDS->SetLineColor(kRed);
    CMEDS->Draw("E1same");
    
    TLegend* legend = new TLegend(0.1,0.7,0.75,0.9);
    legend->AddEntry(CMERS, Form("%f (%f)",CMERS->GetMean(),CMERS->GetMeanError()),"lp");
    legend->AddEntry(CMEDS, Form("%f (%f)",CMEDS->GetMean(),CMEDS->GetMeanError()),"lp");
    legend->Draw();
    
} // end of void plotCum()

// ===========================================================================================

void plotIntFlow(Double_t size, Bool_t bPlotLeg = kFALSE) {
    
    // Access common output file:
    TFile* outputFile = AccessOutputFile(outputFileName);
    TList* outputFileKeys = outputFile->GetListOfKeys();
    TDirectory* directory = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(1)->GetName()));
    TString taskName = "cobjQC";
    TList* QCFlowList = dynamic_cast<TList*>(directory->FindObjectAny(taskName));
    
    TList* IntFlowList = dynamic_cast<TList*>(QCFlowList->FindObject("Integrated Flow"));
    TList* IntFlowResultsList = dynamic_cast<TList*>(IntFlowList->FindObject("Results"));
    TH1D *fIntFlowCorrelationsHist = new TH1D("Average correlations for all events","Average correlations for all events",4,0,4);
    fIntFlowCorrelationsHist = dynamic_cast<TH1D*>(IntFlowResultsList->FindObject("fIntFlowCorrelationsHist"));
    Double_t Int2p = fIntFlowCorrelationsHist->GetBinContent(1);
    Double_t Int2pErr = fIntFlowCorrelationsHist->GetBinError(1);
    Double_t Int4p = fIntFlowCorrelationsHist->GetBinContent(2);
    Double_t IntCum4p = pow(2.*pow(Int2p,2.)-Int4p,0.5);
    TH1D* hInt2p = new TH1D("hInt2p","hInt2p",size,0.,size*1.);
    TH1D* hInt4p = new TH1D("hInt4p","hInt4p",size,0.,size*1.);
    
    TDirectory* MCDir = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(0)->GetName()));
    TList* MCFlowList = dynamic_cast<TList*>(MCDir->FindObjectAny("cobjMCEP"));
    TProfile *MCFlowCumHist = new TProfile("MC Cum","MC Cum",1,0,1);
    MCFlowCumHist = dynamic_cast<TProfile*>(MCFlowList->FindObject("FlowPro_V_MCEP"));
    Double_t MCCum = MCFlowCumHist->GetBinContent(1);
    
    TH1D* hIntMC = new TH1D("hIntMC","hIntMC",size,0.,size*1.);
    for(Int_t bin=0; bin<16; bin++)
    {
        hInt2p->SetBinContent(bin+1,Int2p);
        hInt2p->SetBinError(bin+1,0.);
        hInt4p->SetBinContent(bin+1,IntCum4p);
        hInt4p->SetBinError(bin+1,0.);
        hIntMC->SetBinContent(bin+1,pow(MCCum,2.));
        hIntMC->SetBinError(bin+1,0.);
    }
    hInt2p->SetStats(0);
    hInt2p->SetLineWidth(2);
    hInt2p->SetLineStyle(2);
    hInt2p->SetLineColor(kBlue);
    hInt2p->Draw("E1same");
    // hInt4p->SetStats(0);
    // hInt4p->SetLineWidth(2);
    // hInt4p->SetLineStyle(2);
    // hInt4p->SetLineColor(kRed);
    // hInt4p->Draw("E1same");
    // hIntMC->SetStats(0);
    // hIntMC->SetLineWidth(2);
    // hIntMC->SetLineStyle(2);
    // hIntMC->SetLineColor(kBlack);
    // hIntMC->Draw("E1same");
    
    if(bPlotLeg) {
        TLegend* legend = new TLegend(0.1,0.7,0.35,0.9);
        legend->SetHeader("FlowPackage: ");
        legend->AddEntry(hInt2p, "v_{1}{2}^{2} ","l");
        //  legend->AddEntry(hInt4p, "v_{1}{4}^{2} ","l");
        //  legend->AddEntry(hIntMC, "v_{1}{MC}^{2} ","l");
        legend->Draw();
    }
}

// ===========================================================================================

void plotCum2()
{
    TCanvas* canvas = new TCanvas("Cum2","Cum2",1200,600);
    Int_t const size = 4;
    
    // Access common output file:
    TFile* outputFile = AccessOutputFile(outputFileName);
    TList* outputFileKeys = outputFile->GetListOfKeys();
    TDirectory* directory = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(1)->GetName()));
    
    TString taskName = "cobjQC";
    TList* QCFlowList = dynamic_cast<TList*>(directory->FindObjectAny(taskName));
    TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
    
    TH1D *CEACumHist[2][2][2][2]; // cumulants
    TString HistName[2][2][2][2];
    
    TString plotname = "Something";
    if(TString(directory->GetName()).Contains("POI")) plotname = "Products of flow coefficients";
    else if(TString(directory->GetName()).Contains("RP")) plotname = "Flow products with cumulants (POI-RP)";
    
    TH1D* hFinal = new TH1D("hFinal",plotname,size,0.,size*1.);
    TH1D* hFinal2 = new TH1D("hFinal2",plotname,size,0.,size*1.);
    TH1D* hFinal3 = new TH1D("hFinal3",plotname,size,0.,size*1.);
    TH1D* hFinal4 = new TH1D("hFinal4",plotname,size,0.,size*1.);
    
    Int_t j = 1;
    for(Int_t c=0;c<2;c++)
    {
        for(Int_t y=0;y<2;y++)
        {
            for(Int_t c2=0;c2<2;c2++)
            {
                for(Int_t y2=0;y2<2;y2++)
                {
                    if(j<size+1)
                    {
                        HistName[c][y][c2][y2] = Form("fCEACumulants[%d][%d][%d][%d]",c,y,c2,y2);
                        CEACumHist[c][y][c2][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(HistName[c][y][c2][y2]));
                        hFinal->SetBinContent(j,CEACumHist[c][y][c2][y2]->GetBinContent(1));
                        hFinal->SetBinError(j,CEACumHist[c][y][c2][y2]->GetBinError(1));
                        hFinal2->SetBinContent(j,CEACumHist[c][y][c2][y2]->GetBinContent(2));
                        hFinal2->SetBinError(j,CEACumHist[c][y][c2][y2]->GetBinError(2));
                        hFinal3->SetBinContent(j,CEACumHist[c][y][c2][y2]->GetBinContent(3));
                        hFinal3->SetBinError(j,CEACumHist[c][y][c2][y2]->GetBinError(3));
                        TString charge, charge2;
                        if(c!=c2 && y==y2)
                            hFinal->GetXaxis()->SetBinLabel(j,"#color[4]{#LTv_{1}(c,#eta) v_{1}(-c,#eta)#GT_{c,#eta}}");
                        else if(c==c2 && y==y2)
                            hFinal->GetXaxis()->SetBinLabel(j,"#color[2]{#LTv_{1}(c,#eta) v_{1}(c,#eta)#GT_{c,#eta}}");
                        else if(c==c2 && y!=y2)
                            hFinal->GetXaxis()->SetBinLabel(j,"#color[4]{#LTv_{1}(c,#eta) v_{1}(c,-#eta)#GT_{c,#eta}}");
                        else
                            hFinal->GetXaxis()->SetBinLabel(j,"#color[2]{#LTv_{1}(c,#eta) v_{1}(-c,-#eta)#GT_{c,#eta}}");
                        hFinal->SetStats(0);
                        hFinal->GetXaxis()->SetLabelSize(0.08);
                        hFinal->GetYaxis()->SetLabelSize(0.05);
                        j++;
                    }
                }
            }
        }
    }
    
    Double_t max = hFinal->GetBinContent(hFinal->GetMaximumBin());
    if (hFinal2->GetBinContent(hFinal2->GetMaximumBin()) > max) max = 2.*hFinal2->GetBinContent(hFinal2->GetMaximumBin());
    if (hFinal3->GetBinContent(hFinal3->GetMaximumBin()) > max) max = 2.*hFinal3->GetBinContent(hFinal3->GetMaximumBin());
    if (hFinal4->GetBinContent(hFinal4->GetMaximumBin()) > max) max = 2.*hFinal4->GetBinContent(hFinal4->GetMaximumBin());
    Double_t min = hFinal->GetBinContent(hFinal->GetMinimumBin());
    if (hFinal2->GetBinContent(hFinal2->GetMinimumBin()) < min) min = 0.5*hFinal2->GetBinContent(hFinal2->GetMinimumBin());
    if (hFinal3->GetBinContent(hFinal3->GetMinimumBin()) < min) min = 0.5*hFinal3->GetBinContent(hFinal3->GetMinimumBin());
    if (hFinal4->GetBinContent(hFinal4->GetMinimumBin()) < min) min = 0.5*hFinal4->GetBinContent(hFinal4->GetMinimumBin());
    hFinal->SetAxisRange(0.052,0.07,"Y");
    hFinal->SetLineWidth(2);
    hFinal->SetMarkerColor(kBlue);
    hFinal->SetMarkerStyle(20);
    hFinal->SetLineColor(kBlue);
    hFinal->Draw("E1same");
    hFinal2->SetStats(0);
    hFinal2->SetLineWidth(2);
    hFinal2->SetMarkerColor(kBlack);
    hFinal2->SetMarkerStyle(20);
    hFinal2->SetLineColor(kBlack);
    hFinal2->Draw("E1same");
    hFinal3->SetStats(0);
    hFinal3->SetLineWidth(2);
    hFinal3->SetMarkerColor(kViolet);
    hFinal3->SetMarkerStyle(20);
    hFinal3->SetLineColor(kViolet);
    hFinal3->Draw("E1same");
    
    plotIntFlow(size, kTRUE);
    
    TLegend* legend = new TLegend(0.35,0.5,0.9,0.9);
    legend->SetHeader("v_{1B} OFF, v_{1} = cos(#eta #pi/4) :");
    legend->AddEntry(hFinal, "v_{1}(c,#eta)v_{1}(c',#eta') {2},  |#Delta#eta|>0.4 ","pl");
    legend->AddEntry(hFinal2,"v_{1}(#eta)v_{1}(#eta') {2},  |#Delta#eta|>0.4 ","pl");
    legend->AddEntry(hFinal3,"v_{1}(c,#eta)v_{1}(c',#eta') {2}, |#Delta#eta|>0.4,  #Deltav_{#eta} corr ","pl");
    legend->Draw();
    
    line1 = new TLine(0.,0.,size*1.,0.);
    line1->SetLineStyle(2);
    line1->SetLineColor(1);
    line1->Draw();
    
}

// ===========================================================================================

void plotCEA()
{
    TCanvas* canvas = new TCanvas("CEA","CEA",700,600);
    Int_t const size = 4;
    
    // Access common output file:
    TFile* outputFile = AccessOutputFile(outputFileName);
    TList* outputFileKeys = outputFile->GetListOfKeys();
    TDirectory* directory = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(1)->GetName()));
    
    TString taskName = "cobjQC";
    TList* QCFlowList = dynamic_cast<TList*>(directory->FindObjectAny(taskName));
    TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
    
    TH1D *CEACumHist[2][2]; // cumulants
    TString HistName[2][2];
    
    TString plotname = "Correlation function";
    
    TH1D* hFinal = new TH1D("hFinal",plotname,size,0.,size*1.);
    Double_t AvCEA=0.,AvCEAErr=0.;
    
    Int_t j = 1;
    for(Int_t c=0;c<2;c++) {
        for(Int_t y=0;y<2;y++) {
            if(j<size+1) {
                HistName[c][y] = Form("fCEAFunctions[%d][%d][%d][%d]",c,y,c,y);
                CEACumHist[c][y] = dynamic_cast<TH1D*>(CEAList->FindObject(HistName[c][y]));
                hFinal->SetBinContent(j,CEACumHist[c][y]->GetBinContent(1));
                hFinal->SetBinError(j,CEACumHist[c][y]->GetBinError(1));
                AvCEA += CEACumHist[c][y]->GetBinContent(1);
                AvCEAErr += pow(CEACumHist[c][y]->GetBinError(1),2.);
                TString charge;
                if (c==0) { charge = "+";}
                else { charge = "-";}
                TString eta;
                if (y==0) {eta = "-#eta";}
                else {eta = "#eta";}
                hFinal->GetXaxis()->SetBinLabel(j,Form("C(%s,%s)",charge.Data(),eta.Data()));
                hFinal->SetStats(0);
                hFinal->GetXaxis()->SetLabelSize(0.09);
                hFinal->GetYaxis()->SetLabelSize(0.06);
                j++;
            }
        }
    }
    
    AvCEAErr = pow(AvCEAErr,0.5);
    cout << endl;
    cout << AvCEA << " pm " << AvCEAErr;
    cout << endl;
    AvCEA /= 16.;
    AvCEAErr /= 16;
    cout << endl;
    cout << AvCEA << " pm " << AvCEAErr;
    cout << endl;
    
    Double_t max = hFinal->GetBinContent(hFinal->GetMaximumBin());
    Double_t min = hFinal->GetBinContent(hFinal->GetMinimumBin());
    hFinal->SetAxisRange(0.15,0.45,"Y");
    hFinal->SetLineWidth(2);
    hFinal->SetLineColor(kBlue);
    hFinal->SetMarkerColor(kBlue);
    hFinal->SetMarkerStyle(20);
    hFinal->Draw("E1same");
    
    line2 = new TLine(0.,0.25,size,0.25);
    line2->SetLineStyle(2);
    line2->SetLineColor(1);
    line2->Draw();
    
    TLegend* legend = new TLegend(0.45,0.7,0.9,0.9);
    legend->AddEntry(hFinal, "C(c,#eta) {2},  |#Delta#eta|>0.4 ","pl");
    legend->AddEntry(line2, "4 v_{1B}^{2} ","l");
    legend->Draw();
    
} // end of void plotCEA()


// ===========================================================================================

void plotCorrProd()
{
    TCanvas* canvas = new TCanvas("CorProd","CorProd",1200,600);
    canvas->Divide(1,nTasks);
    Int_t const size = 16;
    
    // #eta
    for(Int_t i=nTasksBegin; i<nTasksEnd+1; i++)
    {
        canvas->cd(i+1);
        // Access common output file:
        TFile* outputFile = AccessOutputFile(outputFileName);
        TList* outputFileKeys = outputFile->GetListOfKeys();
        TDirectory* directory = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(i+1)->GetName()));
        
        TString taskName = "cobjQC";
        TList* QCFlowList = dynamic_cast<TList*>(directory->FindObjectAny(taskName));
        TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
        
        TH1D *CEACorrCovHist[2][2][2][2]; // cumulants
        
        TString plotname = "Something";
        if(TString(directory->GetName()).Contains("POI")) plotname = "Products of flow coefficients (POI-POI)";
        else if(TString(directory->GetName()).Contains("RP")) plotname = "Flow products with cumulants (POI-RP)";
        
        TH1D* hFinal = new TH1D("hFinal",plotname,size,0.,size*1.);
        TH1D* hFinal2 = new TH1D("hFinal2",plotname,size,0.,size*1.);
        TH1D* hFinal3 = new TH1D("hFinal3",plotname,size,0.,size*1.);
        TH1D* hFinal4 = new TH1D("hFinal4",plotname,size,0.,size*1.);
        
        TString HistName[2][2][2][2];
        TString CanvTitle;
        
        Int_t j = 1;
        for(Int_t c=0;c<2;c++)
        {
            for(Int_t y=0;y<2;y++)
            {
                for(Int_t c2=0;c2<2;c2++)
                {
                    for(Int_t y2=0;y2<2;y2++)
                    {
                        if(j<size+1)
                        {
                            HistName[c][y][c2][y2] = Form("fCEACorrProdHist[%d][%d][%d][%d]",c,y,c2,y2);
                            CEACorrCovHist[c][y][c2][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(HistName[c][y][c2][y2]));
                            hFinal->SetBinContent(j,CEACorrCovHist[c][y][c2][y2]->GetBinContent(1));
                            hFinal->SetBinError(j,CEACorrCovHist[c][y][c2][y2]->GetBinError(1));
                            hFinal2->SetBinContent(j,CEACorrCovHist[c][y][c2][y2]->GetBinContent(2));
                            hFinal2->SetBinError(j,CEACorrCovHist[c][y][c2][y2]->GetBinError(2));
                            hFinal3->SetBinContent(j,CEACorrCovHist[c][y][c2][y2]->GetBinContent(3));
                            hFinal3->SetBinError(j,CEACorrCovHist[c][y][c2][y2]->GetBinError(3));
                            TString charge, charge2;
                            if (c==0) { charge = "+";}
                            else { charge = "-";}
                            if (c2==0) { charge2 = "+";}
                            else { charge2 = "-";}
                            TString eta, eta2;
                            if (y==0) {eta = "-y";}
                            else {eta = "y";}
                            if (y2==0) {eta2 = "-y";}
                            else {eta2 = "y";}
                            if(c!=c2 && y==y2)
                                hFinal->GetXaxis()->SetBinLabel(j,Form("#color[4]{(%s,%s)(%s,%s)}",charge.Data(),eta.Data(),charge2.Data(),eta2.Data()));
                            else if(c==c2 && y==y2)
                                hFinal->GetXaxis()->SetBinLabel(j,Form("#color[2]{(%s,%s)(%s,%s)}",charge.Data(),eta.Data(),charge2.Data(),eta2.Data()));
                            else if(c==c2 && y!=y2)
                                hFinal->GetXaxis()->SetBinLabel(j,Form("#color[4]{(%s,%s)(%s,%s)}",charge.Data(),eta.Data(),charge2.Data(),eta2.Data()));
                            else
                                hFinal->GetXaxis()->SetBinLabel(j,Form("#color[2]{(%s,%s)(%s,%s)}",charge.Data(),eta.Data(),charge2.Data(),eta2.Data()));
                            hFinal->SetStats(0);
                            hFinal->GetXaxis()->SetLabelSize(0.09);
                            hFinal->GetYaxis()->SetLabelSize(0.06);
                            j++;
                        }
                    }
                }
            }
        }
        
        Double_t max = hFinal->GetBinContent(hFinal->GetMaximumBin());
        if (hFinal2->GetBinContent(hFinal2->GetMaximumBin()) > max) max = 2.*hFinal2->GetBinContent(hFinal2->GetMaximumBin());
        if (hFinal3->GetBinContent(hFinal3->GetMaximumBin()) > max) max = 2.*hFinal3->GetBinContent(hFinal3->GetMaximumBin());
        if (hFinal4->GetBinContent(hFinal4->GetMaximumBin()) > max) max = 2.*hFinal4->GetBinContent(hFinal4->GetMaximumBin());
        Double_t min = hFinal->GetBinContent(hFinal->GetMinimumBin());
        if (hFinal2->GetBinContent(hFinal2->GetMinimumBin()) < min) min = 0.5*hFinal2->GetBinContent(hFinal2->GetMinimumBin());
        if (hFinal3->GetBinContent(hFinal3->GetMinimumBin()) < min) min = 0.5*hFinal3->GetBinContent(hFinal3->GetMinimumBin());
        if (hFinal4->GetBinContent(hFinal4->GetMinimumBin()) < min) min = 0.5*hFinal4->GetBinContent(hFinal4->GetMinimumBin());
        hFinal->SetAxisRange(min,max,"Y");
        hFinal->SetLineWidth(2);
        hFinal->SetLineColor(kBlue);
        hFinal->Draw("E1same");
        hFinal2->SetStats(0);
        hFinal2->SetLineWidth(2);
        hFinal2->SetLineColor(kViolet);
        hFinal2->Draw("E1same");
        hFinal3->SetStats(0);
        hFinal3->SetLineWidth(2);
        hFinal3->SetLineColor(kRed);
        hFinal3->Draw("E1same");
        //  hFinal4->SetStats(0);
        //  hFinal4->SetLineWidth(2);
        //  hFinal4->SetLineColor(kOrange);
        //  hFinal4->Draw("E1same");
        
        
        TList* IntFlowList = dynamic_cast<TList*>(QCFlowList->FindObject("Integrated Flow"));
        TList* IntFlowResultsList = dynamic_cast<TList*>(IntFlowList->FindObject("Results"));
        TH1D *fIntFlowCorrelationsHist = new TH1D("Average correlations for all events","Average correlations for all events",4,0,4);
        fIntFlowCorrelationsHist = dynamic_cast<TH1D*>(IntFlowResultsList->FindObject("fIntFlowCorrelationsHist"));
        Double_t Int2p = fIntFlowCorrelationsHist->GetBinContent(1);
        Double_t Int2pErr = fIntFlowCorrelationsHist->GetBinError(1);
        Double_t Int4p = fIntFlowCorrelationsHist->GetBinContent(2);
        Double_t IntCum4p = pow(2.*pow(Int2p,2.)-Int4p,0.5);
        TH1D* hInt2p = new TH1D("hInt2p","hInt2p",size,0.,size*1.);
        TH1D* hInt4p = new TH1D("hInt4p","hInt4p",size,0.,size*1.);
        for(Int_t bin=0; bin<16; bin++)
        {
            hInt2p->SetBinContent(bin+1,Int2p);
            hInt2p->SetBinError(bin+1,0.);
            hInt4p->SetBinContent(bin+1,IntCum4p);
            hInt4p->SetBinError(bin+1,0.);
        }
        hInt2p->SetStats(0);
        hInt2p->SetLineWidth(2);
        hInt2p->SetLineStyle(2);
        hInt2p->SetLineColor(kViolet);
        hInt2p->Draw("E1same");
        hInt4p->SetStats(0);
        hInt4p->SetLineWidth(2);
        hInt4p->SetLineStyle(2);
        hInt4p->SetLineColor(kOrange);
        hInt4p->Draw("E1same");
        
        TLegend* legend = new TLegend(0.1,0.5,0.45,0.9);
        legend->AddEntry(hFinal,"<v_{1}{2}(c_{1},y_{1}) v_{1}{2}(c_{2},y_{2})>, #Delta#eta>0 ","l");
        legend->AddEntry(hInt2p,"<v_{1}{2}>^{2}","l");
        legend->AddEntry(hFinal2,"<v_{1}{4}(c_{1},y_{1}) v_{1}{4}(c_{2},y_{2})>, #Delta#eta>0 ","l");
        legend->AddEntry(hInt4p,"<v_{1}{4}>^{2}","l");
        //     legend->Draw();
        
        line1 = new TLine(0.,0.,size*1.,0.);
        line1->SetLineStyle(2);
        line1->SetLineColor(1);
        line1->Draw();
    } // end of for(Int_t i=0; i<nTasksPlot; i++)
    
} // end of void plotCorrProd()

// ===========================================================================================

void plotCEAcompare(Int_t np)
{
    TCanvas* canvas = new TCanvas("CEAsimple","CEAsimple",800,600);
    canvas->Divide(1,2);
    TH1D *hFinal[nTasks];
    TH1D *CEACFunctionsHist[2][2][2][2]; // correlation functions, [c2][c3][y][y2], c=pos,c4=neg
    TLegend* legend = new TLegend(0.7,0.6,0.9,0.9);
    
    Double_t Ymin, Ymax;
    if(np==2)
    {
        Ymin = -3.E-7;
        Ymax = 3.E-7;
    }
    else if (np==4)
    {
        Ymin = -3.E-4;
        Ymax = 3.E-4;
    }
    
    // #eta
    for(Int_t i=0; i<nTasks; i++)
    {
        canvas->cd(1);
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TString taskName = "AnalysisResults.root:";
        taskName += taskSuffix[i];
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(taskName));
        TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
        
        for(Int_t c=0;c<2;c++)
        {
            for(Int_t y=0;y<2;y++)
            {
                for(Int_t c2=0;c2<2;c2++)
                {
                    for(Int_t y2=0;y2<2;y2++)
                    {
                        CEACFunctionsHist[c][c2][y][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(Form("fCEAFunctions[%d][%d][%d][%d]",c,c2,y,y2)));
                    }
                }
            }
        }
        TString plotname;
        if(taskSuffix[i].EqualTo("QCCEA")) plotname = "Correlation Functions";
        else if(taskSuffix[i].EqualTo("QCCEA_PW")) plotname = "Correlation Functions (using particle weights)";
        else if(taskSuffix[i].EqualTo("QCCEA_NUA")) plotname = "Correlation Functions (using NUA corrections)";
        
        hFinal[i] = new TH1D("hFinal",plotname,3,0.,3.);
        hFinal[i]->Sumw2();
        
        if(np==2) Int_t bin = 1;
        else if(np==4) Int_t bin = 2;
        
        hFinal[i]->SetBinContent(1,CEACFunctionsHist[1][0][1][1]->GetBinContent(bin));
        hFinal[i]->SetBinError(1,CEACFunctionsHist[1][0][1][1]->GetBinError(bin));
        hFinal[i]->GetXaxis()->SetBinLabel(1,"C^{+-,+-}(Y,Y)");
        
        hFinal[i]->SetBinContent(2,CEACFunctionsHist[1][0][1][0]->GetBinContent(bin));
        hFinal[i]->SetBinError(2,CEACFunctionsHist[1][0][1][0]->GetBinError(bin));
        hFinal[i]->GetXaxis()->SetBinLabel(2,"C^{+-,+-}(Y,-Y)");
        
        hFinal[i]->SetBinContent(3,CEACFunctionsHist[0][1][1][0]->GetBinContent(bin));
        hFinal[i]->SetBinError(3,CEACFunctionsHist[0][1][1][0]->GetBinError(bin));
        hFinal[i]->GetXaxis()->SetBinLabel(3,"C^{++,--}(Y,-Y)");
        
        hFinal[i]->GetXaxis()->SetLabelSize(0.09);
        hFinal[i]->GetYaxis()->SetLabelSize(0.035);
        hFinal[i]->SetLineWidth(2);
        hFinal[i]->SetLineColor(i+1);
        hFinal[i]->SetStats(0);
        hFinal[i]->SetMarkerStyle(20);
        hFinal[i]->SetMarkerColor(i+1);
        hFinal[i]->SetAxisRange(Ymin,Ymax,"Y");
        hFinal[i]->Draw("E1same");
        
        TString legname;
        if(taskSuffix[i].EqualTo("QCCEA")) legname = Form("CF{%d,QC}",np);
        else if(taskSuffix[i].EqualTo("QCCEA_NUA")) legname = Form("CF{%d,QC} (NUA)",np);
        
        legend->AddEntry(hFinal[i],legname,"l");
        
        for(Int_t c=0;c<2;c++)
        {
            for(Int_t y=0;y<2;y++)
            {
                for(Int_t c2=0;c2<2;c2++)
                {
                    for(Int_t y2=0;y2<2;y2++)
                    {
                        CEACFunctionsHist[c][c2][y][y2] = NULL;
                    }
                }
            }
        }
        
    } // end of for(Int_t i=0; i<nTasks; i++)
    
    legend->Draw();
    
    line1 = new TLine(0.,0.,3.,0.);
    line1->SetLineStyle(2);
    line1->SetLineColor(1);
    line1->Draw();
    
    TH1D* hTh = new TH1D("hTh","hTh",3,0.,3.);
    hTh->SetFillColor(kRed);
    hTh->SetFillStyle(3003);
    hTh->SetMarkerStyle(20);
    hTh->SetBinContent(1,Ymax);
    hTh->SetBinContent(2,0.);
    hTh->SetBinContent(3,Ymin);
    hTh->Draw("bsame");
    
    canvas->cd(2);
    TH1D *hRatio = new TH1D(*hFinal[0]);
    hRatio->SetTitle("Ratio (uncorrelated)");
    
    for(Int_t bin = 1; bin<4; bin++)
    {
        Double_t x = hFinal[0]->GetBinContent(bin);
        Double_t xErr = hFinal[0]->GetBinError(bin);
        Double_t b = hFinal[1]->GetBinContent(bin);
        Double_t bErr = hFinal[1]->GetBinError(bin);
        Double_t f = x/b;
        Double_t fErr = f*pow(pow(xErr/x,2.)+pow(bErr/b,2.),0.5);
        
        hRatio->SetBinContent(bin,f);
        hRatio->SetBinError(bin,fErr);
    }
    
    Ymin = 0.;
    Ymax = 2.;
    hRatio->SetLineWidth(2);
    hRatio->SetLineColor(kBlue);
    hRatio->SetStats(0);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerColor(kBlue);
    hRatio->SetAxisRange(Ymin,Ymax,"Y");
    hRatio->Draw("E1same");
    
    line2 = new TLine(0.,1.,3.,1.);
    line2->SetLineStyle(2);
    line2->SetLineColor(1);
    line2->Draw();
    
} // end of void plotCEAsimple()

// ===========================================================================================

void plotCEAsimple()
{
    TCanvas* canvas = new TCanvas("CEAsimple","CEAsimple",1920,600);
    canvas->Divide(2,1);
    
    // #eta
    for(Int_t i=0; i<nTasks; i++)
    {
        canvas->cd(i+1);
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TString taskName = "AnalysisResults.root:";
        taskName += taskSuffix[i];
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(taskName));
        TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
        
        TH1D *CEACFunctionsHist[2][2][2][2]; // correlation functions, [c2][c3][y][y2], c=pos,c4=neg
        TH1D *CEACorrelationsHist[2][2][2][2]; // mixed <<2'>>, [0=pos,1=neg][0=back,1=forw][0=pos,1=neg][0=back,1=forw]
        
        Int_t j = 1;
        for(Int_t c=0;c<2;c++)
        {
            for(Int_t y=0;y<2;y++)
            {
                for(Int_t c2=0;c2<2;c2++)
                {
                    for(Int_t y2=0;y2<2;y2++)
                    {
                        CEACFunctionsHist[c][c2][y][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(Form("fCEAFunctions[%d][%d][%d][%d]",c,c2,y,y2)));
                        CEACorrelationsHist[c][y][c2][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(Form("fCEACumulants[%d][%d][%d][%d]",c,y,c2,y2)));
                    }
                }
            }
        }
        TString plotname;
        if(taskSuffix[i].EqualTo("QCCEA")) plotname = "Correlation Functions";
        else if(taskSuffix[i].EqualTo("QCCEA_PW")) plotname = "Correlation Functions (using particle weights)";
        else if(taskSuffix[i].EqualTo("QCCEA_NUA")) plotname = "Correlation Functions (using NUA corrections)";
        
        Double_t Ymin = -0.02;
        Double_t Ymax = 0.02;
        
        TH1D* hFinal = new TH1D("hFinal",plotname,3,0.,3.);
        
        hFinal->SetBinContent(1,CEACFunctionsHist[1][0][1][1]->GetBinContent(1));
        hFinal->SetBinError(1,CEACFunctionsHist[1][0][1][1]->GetBinError(1));
        hFinal->GetXaxis()->SetBinLabel(1,"C^{+-,+-}(Y,Y)");
        
        hFinal->SetBinContent(2,CEACFunctionsHist[1][0][1][0]->GetBinContent(1));
        hFinal->SetBinError(2,CEACFunctionsHist[1][0][1][0]->GetBinError(1));
        hFinal->GetXaxis()->SetBinLabel(2,"C^{+-,+-}(Y,-Y)");
        
        hFinal->SetBinContent(3,CEACFunctionsHist[0][1][1][0]->GetBinContent(1));
        hFinal->SetBinError(3,CEACFunctionsHist[0][1][1][0]->GetBinError(1));
        hFinal->GetXaxis()->SetBinLabel(3,"C^{++,--}(Y,-Y)");
        
        hFinal->GetXaxis()->SetLabelSize(0.09);
        hFinal->GetYaxis()->SetLabelSize(0.035);
        hFinal->SetLineWidth(2);
        hFinal->SetLineColor(kBlue);
        hFinal->SetStats(0);
        hFinal->SetMarkerStyle(20);
        hFinal->SetMarkerColor(kBlue);
        hFinal->SetAxisRange(Ymin,Ymax,"Y");
        hFinal->Draw("E1same");
        
        TH1D* hFinal4p = new TH1D("hFinal4p",plotname,3,0.,3.);
        
        hFinal4p->SetBinContent(1,CEACFunctionsHist[1][0][1][1]->GetBinContent(2));
        hFinal4p->SetBinError(1,CEACFunctionsHist[1][0][1][1]->GetBinError(2));
        hFinal4p->GetXaxis()->SetBinLabel(1,"C^{+-,+-}(Y,Y)");
        
        hFinal4p->SetBinContent(2,CEACFunctionsHist[1][0][1][0]->GetBinContent(2));
        hFinal4p->SetBinError(2,CEACFunctionsHist[1][0][1][0]->GetBinError(2));
        hFinal4p->GetXaxis()->SetBinLabel(2,"C^{+-,+-}(Y,-Y)");
        
        hFinal4p->SetBinContent(3,CEACFunctionsHist[0][1][1][0]->GetBinContent(2));
        hFinal4p->SetBinError(3,CEACFunctionsHist[0][1][1][0]->GetBinError(2));
        hFinal4p->GetXaxis()->SetBinLabel(3,"C^{++,--}(Y,-Y)");
        
        hFinal4p->SetLineColor(2);
        hFinal4p->SetLineWidth(2);
        hFinal4p->SetStats(0);
        hFinal4p->SetMarkerStyle(20);
        hFinal4p->SetMarkerColor(2);
        hFinal4p->Draw("E1same");
        
        TLegend* legend = new TLegend(0.7,0.1,0.9,0.2);
        legend->AddEntry(hFinal,"CF{2,QC}","l");
        legend->AddEntry(hFinal4p,"CF{4,QC}","l");
        legend->Draw();
        
        line1 = new TLine(0.,0.,3.,0.);
        line1->SetLineStyle(2);
        line1->SetLineColor(1);
        line1->Draw();
        
        TH1D* hTh = new TH1D("hTh","hTh",3,0.,3.);
        hTh->SetFillColor(kRed);
        hTh->SetFillStyle(3003);
        hTh->SetMarkerStyle(20);
        hTh->SetBinContent(1,Ymax);
        hTh->SetBinContent(2,0.);
        hTh->SetBinContent(3,Ymin);
        hTh->Draw("bsame");
        
    } // end of for(Int_t i=0; i<nTasks; i++)
    
} // end of void plotCEAsimple()

// ===========================================================================================

void plotCEAsimpleDiff()
{
    TCanvas* canvas = new TCanvas("CEAsimple","CEAsimple",1920,600);
    canvas->Divide(nRows,nColumns);
    
    // #eta
    for(Int_t i=0; i<nTasks; i++)
    {
        canvas->cd(i+1);
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TString taskName = "AnalysisResults.root:";
        taskName += taskSuffix[i];
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(taskName));
        TList* CEAList = dynamic_cast<TList*>(QCFlowList->FindObject("Charge-Eta Asymmetry"));
        
        TH1D *CEACFun2pHistDiff[2][2][2][2]; // correlation functions, [c2][c3][y][y2], c=pos,c4=neg
        TH1D *CEACFun4pHistDiff[2][2][2][2]; // mixed <<2'>>, [0=pos,1=neg][0=back,1=forw][0=pos,1=neg][0=back,1=forw]
        TH1D* hFinal2p[nPtbins];
        TH1D* hFinal4p[nPtbins];
        
        Int_t j = 1;
        for(Int_t c=0;c<2;c++)
        {
            for(Int_t y=0;y<2;y++)
            {
                for(Int_t c2=0;c2<2;c2++)
                {
                    for(Int_t y2=0;y2<2;y2++)
                    {
                        CEACFun2pHistDiff[c][c2][y][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(Form("fCEACFunctions2pDiffHist[%d][%d][%d][%d]",c,c2,y,y2)));
                        CEACFun4pHistDiff[c][c2][y][y2] = dynamic_cast<TH1D*>(CEAList->FindObject(Form("fCEACFunctions4pDiffHist[%d][%d][%d][%d]",c,c2,y,y2)));
                    }
                }
            }
        }
        
        TLegend* legend = new TLegend(0.7,0.1,0.9,0.3);
        
        for(Int_t b=1; b<=nPtbins; b++)
        {
            hFinal2p[b] = new TH1D(Form("hFinal2p%d",b),Form("hFinal2p%d",b),3,0.,3.);
            hFinal2p[b]->GetXaxis()->SetBinLabel(1,"C^{+-,+-}(Y,Y)");
            hFinal2p[b]->GetXaxis()->SetBinLabel(2,"C^{+-,+-}(Y,-Y)");
            hFinal2p[b]->GetXaxis()->SetBinLabel(3,"C^{++,--}(Y,-Y)");
            
            hFinal2p[b]->SetBinContent(1,CEACFun2pHistDiff[1][0][1][1]->GetBinContent(b));
            hFinal2p[b]->SetBinError(1,CEACFun2pHistDiff[1][0][1][1]->GetBinError(b));
            
            hFinal2p[b]->SetBinContent(2,CEACFun2pHistDiff[1][0][1][0]->GetBinContent(b));
            hFinal2p[b]->SetBinError(2,CEACFun2pHistDiff[1][0][1][0]->GetBinError(b));
            
            hFinal2p[b]->SetBinContent(3,CEACFun2pHistDiff[0][1][1][0]->GetBinContent(b));
            hFinal2p[b]->SetBinError(3,CEACFun2pHistDiff[0][1][1][0]->GetBinError(b));
            
            hFinal2p[b]->GetXaxis()->SetLabelSize(0.09);
            hFinal2p[b]->GetYaxis()->SetLabelSize(0.035);
            hFinal2p[b]->SetLineWidth(2);
            hFinal2p[b]->SetLineColor(b);
            hFinal2p[b]->SetStats(0);
            hFinal2p[b]->SetMarkerStyle(20);
            hFinal2p[b]->SetMarkerColor(b);
            hFinal2p[b]->SetAxisRange(-0.02,0.02,"Y");
            hFinal2p[b]->Draw("E1same");
            
            //      hFinal4p[b] = new TH1D(Form("hFinal4p%d",b),Form("hFinal4p%d",b),3,0.,3.);
            //      hFinal4p[b]->GetXaxis()->SetBinLabel(1,"C^{+-,+-}(Y,Y)");
            //      hFinal4p[b]->GetXaxis()->SetBinLabel(2,"C^{+-,+-}(Y,-Y)");
            //      hFinal4p[b]->GetXaxis()->SetBinLabel(3,"C^{++,--}(Y,-Y)");
            //
            //      hFinal4p[b]->SetBinContent(1,CEACFun4pHistDiff[1][0][1][1]->GetBinContent(b));
            //      hFinal4p[b]->SetBinError(1,CEACFun4pHistDiff[1][0][1][1]->GetBinError(b));
            //
            //      hFinal4p[b]->SetBinContent(2,CEACFun4pHistDiff[1][0][1][0]->GetBinContent(b));
            //      hFinal4p[b]->SetBinError(2,CEACFun4pHistDiff[1][0][1][0]->GetBinError(b));
            //
            //      hFinal4p[b]->SetBinContent(3,CEACFun4pHistDiff[0][1][1][0]->GetBinContent(b));
            //      hFinal4p[b]->SetBinError(3,CEACFun4pHistDiff[0][1][1][0]->GetBinError(b));
            //
            //      hFinal4p[b]->SetAxisRange(-0.01,0.01,"Y");
            //      hFinal4p[b]->SetLineColor(b+8);
            //      hFinal4p[b]->SetLineWidth(2);
            //      hFinal4p[b]->SetStats(0);
            //      hFinal4p[b]->SetMarkerStyle(20);
            //      hFinal4p[b]->SetMarkerColor(b+8);
            //      hFinal4p[b]->Draw("E1same");
            
            legend->AddEntry(hFinal2p[b],Form("CF{2,QC} Ptbin%d",b),"l");
            //     legend->AddEntry(hFinal4p[b],Form("CF{4,QC} Ptbin%d",b),"l");
        }
        
        legend->Draw();
        
        line1 = new TLine(0.,0.,3.,0.);
        line1->SetLineStyle(2);
        line1->SetLineColor(1);
        line1->Draw();
        
    } // end of for(Int_t i=0; i<nTasks; i++)
    
} // end of void plotCEAsimple()

// ===========================================================================================

void plotDiffFlow(TString n)
{
    TCanvas* canvas = new TCanvas("POI#eta","POI#eta",1000,800);
    canvas->Divide(2,2);
    TLegend* legend = new TLegend(0.6,0.8,1.0,1.0);
    
    // #eta
    for(Int_t i=0; i<2; i++)
    {
        canvas->cd(i+1);
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(Form("AnalysisResults.root:QCFlow_v%s",n.Data())));
        TList* DiffFlowList = dynamic_cast<TList*>(QCFlowList->FindObject("Differential Flow"));
        TList* ResultsList = dynamic_cast<TList*>(DiffFlowList->FindObject("Results"));
        TList* POIetaList = dynamic_cast<TList*>(ResultsList->FindObject("Differential flow (POI, #eta)"));
        
        TH1F* hVarPos = dynamic_cast<TH1F*>(POIetaList->At(i));
        //if(i=1) {canvas->SetLogy();}
        hVarPos->SetAxisRange(-0.79,0.79,"X");
        hVarPos->SetAxisRange(-0.1,0.1,"Y");
        //hVarPos->Rebin(8);
        hVarPos->SetLineColor(2);
        hVarPos->Draw("Esame");
        //hVarNeg->Rebin(8);
        if(i==1) {legend->AddEntry(hVarPos,outputFileName,"l");}
        TLine *line = new TLine(-0.8,0,0.8,0);
        line->SetLineStyle(2);
        line->Draw();
    } // end of for(Int_t i=0; i<2; i++)
    
    // p_{T}
    for(Int_t i=2; i<4; i++)
    {
        canvas->cd(i+1);
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(Form("AnalysisResults.root:QCFlow_v%s",n.Data())));
        TList* DiffFlowList = dynamic_cast<TList*>(QCFlowList->FindObject("Differential Flow"));
        TList* ResultsList = dynamic_cast<TList*>(DiffFlowList->FindObject("Results"));
        TList* POIpTList = dynamic_cast<TList*>(ResultsList->FindObject("Differential flow (POI, p_{T})"));
        
        TH1F* hVarPos = dynamic_cast<TH1F*>(POIpTList->At(i-2));
        for (Int_t j=1; j<=hVarPos->GetNbinsX(); j++)
        {
            Double_t inv = hVarPos->GetBinContent(j);
            hVarPos->SetBinContent(j,inv);
        }
        //if(i=1) {canvas->SetLogy();}
        hVarPos->SetAxisRange(0.,5.,"X");
        hVarPos->SetAxisRange(-0.2,0.2,"Y");
        //hVarPos->Rebin(8);
        hVarPos->SetLineColor(2);
        hVarPos->Draw("Esame");
        //hVarNeg->Rebin(8);
        TLine *line = new TLine(-0.1,0,5.1,0);
        line->SetLineStyle(2);
        line->Draw();
    } // end of for(Int_t i=0; i<2; i++)
    
    legend->Draw();
    
    
} // end of void plotDiffFlow()

// ===========================================================================================

void plotVariousCompare()
{
    // Access common output file:
    TFile *outputFile = AccessOutputFile("correctedAnalysisResults.root");
    TString taskName = "AnalysisResults.root:";
    taskName += taskSuffix[0];
    TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(taskName));
    TList* VariousList = dynamic_cast<TList*>(QCFlowList->FindObject("Various"));
    
    TCanvas* canvas = new TCanvas("VariousCompare_PW","VariousCompare_PW",1000,800);
    canvas->Divide(nColumns,2);
    
    for(Int_t k=0; k<nTasks; k++)
    {
        // Access common output file:
        TFile *outputFile = AccessOutputFile(outputFileName);
        TString taskName = "AnalysisResults.root:";
        taskName += taskSuffix[k];
        TList* QCFlowList = dynamic_cast<TList*>(outputFile->FindObjectAny(taskName));
        TList* VariousList = dynamic_cast<TList*>(QCFlowList->FindObject("Various"));
        
        canvas->cd(k+1);
        
        TH1F* hVarP = dynamic_cast<TH1F*>(VariousList->FindObject("phi_weights+"));
        TH1F* hVarN = dynamic_cast<TH1F*>(VariousList->FindObject("phi_weights-"));
        //if(i=1) {canvas->SetLogy();}
        //hVarPos->Rebin(8);
        hVarP->SetLabelSize(0.05,"X");
        hVarP->GetXaxis()->SetTitleSize(0.05);
        hVarP->SetLabelSize(0.05,"Y");
        hVarP->GetYaxis()->SetTitleSize(0.05);
        hVarP->GetYaxis()->SetTitle("dN/d#phi nBins/Ntot (1/rad)");
        hVarP->SetLineColor(2);
        hVarP->SetTitle("#phi distribution");
        hVarP->SetStats(0);
        hVarP->SetAxisRange(0.,1.,"Y");
        hVarP->Draw("Esame");
        hVarN->SetLineColor(4);
        hVarN->Draw("Esame");
        
        TLegend* legend = new TLegend(0.8,0.7,0.9,0.9);
        legend->AddEntry(hVarP," #pi^{+}","l");
        legend->AddEntry(hVarN," #pi^{-}","l");
        legend->Draw();
        
        canvas->cd(k+nTasks+1);
        
        TH1F* hVarP = dynamic_cast<TH1F*>(VariousList->FindObject("eta_weights+"));
        TH1F* hVarN = dynamic_cast<TH1F*>(VariousList->FindObject("eta_weights-"));
        //if(i=1) {canvas->SetLogy();}
        //hVarPos->Rebin(8);
        hVarP->SetLabelSize(0.05,"X");
        hVarP->GetXaxis()->SetTitleSize(0.05);
        hVarP->SetLabelSize(0.05,"Y");
        hVarP->GetYaxis()->SetTitleSize(0.05);
        hVarP->GetYaxis()->SetTitle("dN/d#eta nBins/Ntot");
        hVarP->SetLineColor(2);
        hVarP->SetTitle("#eta distribution");
        hVarP->SetStats(0);
        hVarP->SetAxisRange(0.,1.,"Y");
        hVarP->Draw("Esame");
        hVarN->SetLineColor(4);
        hVarN->Draw("Esame");
        
        TLegend* legend = new TLegend(0.8,0.7,0.9,0.9);
        legend->AddEntry(hVarP," #pi^{+}","l");
        legend->AddEntry(hVarN," #pi^{-}","l");
        legend->Draw();
        
        TLegend* legend = new TLegend(0.8,0.7,0.9,0.9);
        legend->AddEntry(hVarP," #pi^{+}","l");
        legend->AddEntry(hVarN," #pi^{-}","l");
        legend->Draw();
    }
    
} // end of void plotVariousCompare()

// ===========================================================================================

void plotEventCutQA()
{
    // Access common output file:
    TFile *outputFile = AccessOutputFile(outputFileName);
    
    TList* CutsQAList = dynamic_cast<TList*>(outputFile->FindObjectAny("AnalysisResults.root:CutsQA+"));
    if(!CutsQAList) {Warning("AnalysisResults.root:CutsQA+");}
    TList* EventCutsQAList = dynamic_cast<TList*>(CutsQAList->FindObject("EventCuts QA"));
    if(!EventCutsQAList) {Warning("EventCuts QA");}
    
    TList* beforeList = dynamic_cast<TList*>(EventCutsQAList->FindObject("before"));
    TList* afterList = dynamic_cast<TList*>(EventCutsQAList->FindObject("after"));
    
    TCanvas* canvas = new TCanvas("EventCutsQA","EventCutsQA",750,600);
    canvas->Divide(beforeList->GetEntries(),2);
    for(Int_t i=0; i<beforeList->GetEntries(); i++)
    {
        TH1F* hBefore = dynamic_cast<TH1F*>(beforeList->At(i));
        TH1F* hAfter = dynamic_cast<TH1F*>(afterList->At(i));
        //canvas->SetLogy();
        //hbefore->SetAxisRange(0.1,hbefore->GetBinContent(hbefore->GetMaximumBin()),"Y");
        hBefore->Rebin(4);
        hBefore->SetLineColor(kBlue);
        canvas->cd(i+1);
        hBefore->DrawCopy("E");
        hAfter->Rebin(4);
        hAfter->SetLineColor(kRed);
        canvas->cd(i+1+beforeList->GetEntries());
        hAfter->DrawCopy("E");
    } // end of for(Int_t i=0; i<VariousPList->GetEntries(); i++)
    TLegend* legend = new TLegend(0.8,0.8,1.0,1.0);
    legend->AddEntry(hBefore,"before","l");
    legend->AddEntry(hAfter,"after","l");
    legend->Draw();
    
} // end of void plotEventCutQA()

// ===========================================================================================

void plotVarious()
{
    // Access common output file:
    TFile* outputFile = AccessOutputFile(outputFileName);
    TList* outputFileKeys = outputFile->GetListOfKeys();
    TDirectory* directory = dynamic_cast<TDirectory*>(outputFile->Get(outputFileKeys->At(1)->GetName()));
    
    TString taskName = "cobjQC";
    TList* QCFlowList = dynamic_cast<TList*>(directory->FindObjectAny(taskName));
    TList* VariousList = dynamic_cast<TList*>(QCFlowList->FindObject("Various"));
    
    if(VariousList->GetEntries())
    {
        TCanvas* canvas = new TCanvas("Various","Various",1000,800);
        canvas->Divide(2,2);
        
        canvas->cd(1);
        TH1D* fPhiDistrRPs = dynamic_cast<TH1D*>(VariousList->FindObject("fPhiDistrRPs"));
        fPhiDistrRPs->Draw("E");
        
        canvas->cd(2);
        TH1D* fPhiDistrPOIsP = dynamic_cast<TH1D*>(VariousList->FindObject("fPhiDistrPOIs+"));
        fPhiDistrPOIsP->SetLineColor(3);
        fPhiDistrPOIsP->Draw("Esame");
        TH1D* fPhiDistrPOIsN = dynamic_cast<TH1D*>(VariousList->FindObject("fPhiDistrPOIs-"));
        fPhiDistrPOIsN->SetLineColor(2);
        fPhiDistrPOIsN->Draw("Esame");
        
        canvas->cd(3);
        TH1D* fEtaDistrRPs = dynamic_cast<TH1D*>(VariousList->FindObject("fEtaDistrRPs"));
        fEtaDistrRPs->Draw("E");
        
        canvas->cd(4);
        TH1D* fEtaDistrPOIsP = dynamic_cast<TH1D*>(VariousList->FindObject("fEtaDistrPOIs+"));
        fEtaDistrPOIsP->SetLineColor(3);
        fEtaDistrPOIsP->Draw("Esame");
        TH1D* fEtaDistrPOIsN = dynamic_cast<TH1D*>(VariousList->FindObject("fEtaDistrPOIs-"));
        fEtaDistrPOIsN->SetLineColor(2);
        fEtaDistrPOIsN->Draw("Esame");
        
    } // end of if(VariousPList->GetEntries() == VariousMList->GetEntries())
    
} // end of void plotVarious()

// ===========================================================================================

TFile* AccessOutputFile(TString outputFileName)
{
    // Access the common output file.
    
    TFile *outputFile = NULL;
    if(!(gSystem->AccessPathName(Form("%s%s%s",gSystem->pwd(),"/",outputFileName.Data()),kFileExists)))
    {
        outputFile = TFile::Open(outputFileName.Data(),"READ");
    } else
    {
        cout<<endl;
        cout<<"WARNING: Couldn't find the file "<<outputFileName.Data()<<" in "<<endl;
        cout<<"         directory "<<gSystem->pwd()<<" !!!!"<<endl;
        cout<<endl;
        exit(0);
    }
    
    return outputFile;
    
} // end of TFile* AccessOutputFile(TString outputFileName)

// ===========================================================================================

void GetListsWithHistograms(TFile *outputFile, TString analysisType)
{
    // Access from common output file the TDirectoryFile's for each flow analysis method
    // and from them the TList's holding histograms with final results:
    
    const Int_t n = 2;
    TString fileName[n];
    TDirectoryFile *dirFile[n] = {NULL};
    TString listName[n];
    Int_t failureCounter = 0;
    for(Int_t i=0;i<n;i++)
    {
        // Form a file name for each method:
        fileName[i]+="QC_output_for_n=";
        fileName[i]+=n;
        // Access this file:
        dirFile[i] = (TDirectoryFile*)outputFile->FindObjectAny(fileName[i]);
        // Form a list name for each method:
        if(dirFile[i])
        {
            TList* listTemp = dirFile[i]->GetListOfKeys();
            if(listTemp && listTemp->GetEntries() == 1)
            {
                listName[i] = listTemp->At(0)->GetName(); // to be improved - implemented better (use dynamic_cast instead)
                dirFile[i]->GetObject(listName[i].Data(),list[i]);
            } else
            {
                cout<<" WARNING: Accessing TList from TDirectoryFile failed for method "<<method[i].Data()<<" !!!!"<<endl;
                cout<<"          Did you actually use "<<method[i].Data()<<" in the analysis?"<<endl;
                cout<<endl;
            }
        } else
        {
            cout<<" WARNING: Couldn't find a TDirectoryFile "<<fileName[i].Data()<<".root !!!!"<<endl;
            failureCounter++;
        }
    } // end of for(Int_t i=0;i<nMethods;i++)
    
    // If no files were found most probably the 'TString analysisType' was specified wrongly:
    if(failureCounter == nMethods)
    {
        cout<<endl;
        cout<<"Did you specify 'TString analysisType' correctly? Can be \"\", \"ESD\", \"AOD\", etc."<<endl;
        cout<<endl;
        exit(0);
    }
    
} // end of void GetListsWithHistograms(TFile *outputFile, TString analysisType)

// ===========================================================================================

void ALICESettings(TH1D* histo, TString xAxisTitle="", TString yAxisTitle="")
{
    histo->SetMarkerStyle(20);
    histo->GetXaxis()->SetLabelSize(0.045);
    histo->GetYaxis()->SetLabelSize(0.045);
    histo->GetXaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetLabelOffset(0.01);
    histo->GetYaxis()->SetTitleSize(0.05);
    histo->GetXaxis()->SetTitleSize(0.05);
    histo->GetYaxis()->SetTitleOffset(1.25);
    histo->GetXaxis()->SetTitleOffset(1.1);
    histo->GetXaxis()->SetTitle(xAxisTitle);
    histo->GetYaxis()->SetTitle(yAxisTitle);
}

void SetColorHisto(TH1D* v24standard, TColor any)
{
    v24standard->SetLineColor(any);
    v24standard->SetMarkerStyle(20);
    v24standard->SetMarkerColor(any);
}

void SetColorHisto(TH1D* v24standard, Int_t any)
{
    v24standard->SetLineColor(any);
    v24standard->SetMarkerStyle(20);
    v24standard->SetMarkerColor(any);
}

void GlobalSettings()
{
    // Settings which will affect all plots.
    
    gROOT->SetStyle("Plain"); // default color is white instead of gray
    gROOT->ForceStyle();
    gStyle->SetOptStat(0); // remove stat. box from all histos
    gStyle->SetTitleX(0.2);
    gStyle->SetLabelSize(0.05,"X");
    gStyle->SetLineWidth(2);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerStyle(20);
    //gStyle->SetErrorX(0);
    //  gStyle->SetNdivisions(10,"Y");
    TGaxis::SetMaxDigits(3); // prefer exp notation for 2 and more significant digits
    gStyle->SetErrorX(0.);
    gStyle->SetMarkerStyle(20);
    gStyle->SetErrorX(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.25,"y");
    gStyle->SetTitleOffset(1.1,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
    
} // end of void GlobalSettings()

// ===========================================================================================

// ===========================================================================================

void LoadLibrariesCFR() {
    
    //--------------------------------------
    // Load the needed libraries most of them already loaded by aliroot
    //--------------------------------------
    //gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libVMC");
    gSystem->Load("libXMLIO");
    gSystem->Load("libPhysics");
    //--------------------------------------------------------
    // If you want to use already compiled libraries
    // in the aliroot distribution
    //--------------------------------------------------------
    
    //==================================================================================
    //load needed libraries:
    gSystem->AddIncludePath("-I$ROOTSYS/include");
    //gSystem->Load("libTree");
    
    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWGflowBase");
    //cerr<<"libPWGflowBase loaded ..."<<endl;
    
} // end of void LoadLibrariesCFR(const libModes analysisMode)

void Warning(TString name)
{
    cout<<endl;
    cout<<" WARNING: Couldn't find "<<name<<" !!!"<<endl;
    cout<<endl;
    exit(0);
}
