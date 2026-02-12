// #include <iostream>
// using std::cout;
// using std::endl;

// #include "TRandom.h"
// #include "TH1.h"
// #include "TH2.h" 
// #include "TCanvas.h"
// #include "TFile.h"

// #include "RooUnfoldResponse.h"
// #include "RooUnfoldBayes.h"

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif



void Unfolding(){

    /*read the histogram from a file.root*/
    TFile *f = TFile::Open("Unfolding_output/SH_pre_unfolding.root");

    /*MC signal*/
    TH2F *res = (TH2F*)f->Get("h_matrix");
    TH1F *tru = (TH1F*)f->Get("h_true");
    TH1F *reco = (TH1F*)f->Get("h_reco");

    /*C3H6 measured*/
    TH1F *CH = (TH1F*)f->Get("CH_reco");

    /*C background*/
    TH1F *C = (TH1F*)f->Get("C_reco");

    /*H from substraction*/
    //TH1F *H_sub= (TH1F*)(*CH - *C);
    TH1F *H_sub = (TH1F*)CH->Clone("H_sub");
    H_sub->Add(C, -5.2);

    /*acceptance and efficiency*/
    // TH1F *acc = (TH1F*)f->Get("h_accep");
    // TH1F *eff = (TH1F*)f->Get("h_eff");
    TH1D* hRecoProj  = res->ProjectionX("hRecoProj");
    TH1D* hTruthProj = res->ProjectionY("hTruthProj");
    RooUnfoldResponse response(hRecoProj, hTruthProj, res);

    TH1D* acc = dynamic_cast<TH1D*>(hRecoProj->Clone("acc"));
    acc->Divide(hRecoProj, reco, 1.0, 1.0, "B");

    TH1D* eff = dynamic_cast<TH1D*>(hTruthProj->Clone("eff"));
    eff->Divide(hTruthProj, tru, 1.0, 1.0, "B");

    /*correction to h_sub, multiply for the acceptance*/
    TH1F *H_corr = (TH1F*)H_sub->Clone("H_corr");
    H_corr->Multiply(acc);


    /*Unfolded distribution*/
    RooUnfoldBayes    unfold (&response, H_corr, 4);
    //RooUnfoldBayes    unfold (&response, reco, 4);
    TH1D *hUnfold = (TH1D *)unfold.Hunfold();

    /*Correction to the unfolded distribution*/
    TH1D *hUnfCorr = (TH1D*) hUnfold->Clone("hUnfCorr");
    hUnfCorr->Divide(eff); 
    
    
    
    
    
    TFile *fout = new TFile("Unfolding_output/Unfolding_results.root", "RECREATE");
    
    H_sub->Write();
    H_corr->Write();
    hUnfold->Write();
    hUnfCorr->Write();
    acc->Write();
    eff->Write();
    res->Write();
    tru->Write();

    TCanvas *c = new TCanvas("c", "Canvas", 1600, 600);
    tru->Draw();
    hUnfCorr->Draw("SAME");

    c->Write();


    fout->Close();

}