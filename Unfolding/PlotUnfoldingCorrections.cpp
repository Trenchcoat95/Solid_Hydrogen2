//==============================================================
// PlotUnfoldingCorrections.C (robust version + raw response view)
// - Acceptance, Efficiency, Response (normalized + raw)
// - tolerant to different histogram names and missing objects
//==============================================================

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TColor.h"
#include "TAxis.h"
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

const char* tryNamesAcceptance[] = {"hAcceptance", "Acceptance", "acc", nullptr};
const char* tryNamesEfficiency[] = {"hEfficiency", "Efficiency", "eff", nullptr};
const char* tryNamesResponse[]   = {"hResponse", "hResponse2D", "h_matrix", "response", nullptr};

TH1D* Find1D(TFile* f, const char** names) {
    for (int i=0; names[i]; ++i) {
        TH1D* h = dynamic_cast<TH1D*>(f->Get(names[i]));
        if (h) {
            cout << "[OK] Found 1D: " << names[i] << endl;
            return h;
        }
    }
    cout << "[WARN] No 1D found among provided names." << endl;
    return nullptr;
}

TH2D* Find2D(TFile* f, const char** names) {
    for (int i=0; names[i]; ++i) {
        TH2D* h = dynamic_cast<TH2D*>(f->Get(names[i]));
        if (h) {
            cout << "[OK] Found 2D: " << names[i] << endl;
            return h;
        }
    }
    cout << "[WARN] No 2D response found among provided names." << endl;
    return nullptr;
}

void SaveCanvasIfValid(TCanvas* c, const TString& base, const char* suffix) {
    if (!c) return;
    //TString png = TString::Format("%s_%s.png", base.Data(), suffix);
    TString pdf = TString::Format("%s_%s.pdf", base.Data(), suffix);
    //c->SaveAs(png);
    c->SaveAs(pdf);
    //cout << "[SAVED] " << png.Data() << " , " << pdf.Data() << endl;
    cout << "[SAVED] " << pdf.Data() << endl;
}

//==============================================================
// Main macro
//==============================================================
void PlotUnfoldingCorrections(const char* filename = "Unfolding_output/Unfolding_results.root") {
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("5.1f");
    gStyle->SetTitleFontSize(0.045);

    TFile* f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        ::Error("PlotUnfoldingCorrections", "Cannot open file %s", filename);
        return;
    }

    // find histos (tolerant names)
    TH1D* hAcc = Find1D(f, tryNamesAcceptance);
    TH1D* hEff = Find1D(f, tryNamesEfficiency);
    TH2D* hResp = Find2D(f, tryNamesResponse);

    TString base = filename;
    base.ReplaceAll(".root", "");

    // --- 1) Acceptance canvas ---
    TCanvas* cAcc = new TCanvas("cAcceptance", "Acceptance", 800, 600);
    if (hAcc) {
        cAcc->SetGrid();
        hAcc->SetLineColor(kBlue+1);
        hAcc->SetMarkerColor(kBlue+1);
        hAcc->SetMarkerStyle(21);
        hAcc->SetTitle(";Reconstructed neutrino energy (MeV);Acceptance ");
        hAcc->SetMinimum(0);
        double ymaxAcc = hAcc->GetMaximum();
        if (ymaxAcc > 0) hAcc->SetMaximum(1.2 * ymaxAcc);
        hAcc->Draw("E");
    } else {
        cAcc->cd();
        TH1D* hframe = new TH1D("frameAcc","Acceptance (not found)", 1, 0, 1);
        hframe->SetTitle("Acceptance not found in file");
        hframe->Draw();
    }
    SaveCanvasIfValid(cAcc, base, "Acceptance");

    // --- 2) Efficiency canvas ---
    TCanvas* cEff = new TCanvas("cEfficiency", "Efficiency", 800, 600);
    if (hEff) {
        cEff->SetGrid();
        hEff->SetLineColor(kGreen+2);
        hEff->SetMarkerColor(kGreen+2);
        hEff->SetMarkerStyle(20);
        hEff->SetTitle(";True neutrino energy (MeV);Efficiency ");
        hEff->SetMinimum(0);
        double ymaxEff = hEff->GetMaximum();
        if (ymaxEff > 0) hEff->SetMaximum(1.2 * ymaxEff);
        hEff->Draw("E");
    } else {
        cEff->cd();
        TH1D* hframe2 = new TH1D("frameEff","Efficiency (not found)", 1, 0, 1);
        hframe2->SetTitle("Efficiency not found in file");
        hframe2->Draw();
    }
    SaveCanvasIfValid(cEff, base, "Efficiency");

    // --- 3) Response (raw, not normalized) ---
    TCanvas* cRespRaw = new TCanvas("cResponseRaw", "Response Matrix (Raw)", 800, 600);
    if (hResp) {
        cRespRaw->cd();
        cRespRaw->SetRightMargin(0.15);
        cRespRaw->SetGrid();
        gStyle->SetPalette(kBird); // default nice ROOT palette
	hResp->SetTitle("");
        hResp->GetXaxis()->SetTitle("Reconstructed neutrino energy (MeV)");
        hResp->GetYaxis()->SetTitle("True neutrino energy (MeV)");
        hResp->Draw("COLZ TEXT ");
    } else {
        cRespRaw->cd();
        TH1D* hframe3 = new TH1D("frameRespRaw","Response (not found)", 1, 0, 1);
        hframe3->SetTitle("Response matrix not found in file");
        hframe3->Draw();
    }
    SaveCanvasIfValid(cRespRaw, base, "ResponseRaw");

    // --- 4) Response normalized matrix ---
    TCanvas* cResp = new TCanvas("cResponseNorm", "Normalized Response Matrix", 800, 600);
    if (!hResp) {
        cResp->cd();
        TH1D* hframe4 = new TH1D("frameResp","Response (not found)", 1, 0, 1);
        hframe4->SetTitle("Response matrix not found in file");
        hframe4->Draw();
        SaveCanvasIfValid(cResp, base, "ResponseNorm");
        f->Close();
        return;
    }

    int nBinsX = hResp->GetNbinsX();
    int nBinsY = hResp->GetNbinsY();
    std::vector<double> binEdgesX(nBinsX + 1);
    std::vector<double> binEdgesY(nBinsY + 1);

    for (int i = 1; i <= nBinsX; ++i)
        binEdgesX[i-1] = hResp->GetXaxis()->GetBinLowEdge(i);
    binEdgesX[nBinsX] = hResp->GetXaxis()->GetBinUpEdge(nBinsX);

    for (int j = 1; j <= nBinsY; ++j)
        binEdgesY[j-1] = hResp->GetYaxis()->GetBinLowEdge(j);
    binEdgesY[nBinsY] = hResp->GetYaxis()->GetBinUpEdge(nBinsY);

    // normalize raws
    std::vector<double> Psum(nBinsY, 0.0);
    for (int j = 1; j <= nBinsY; ++j) {
        double sum = 0.0;
        for (int i = 1; i <= nBinsX; ++i) sum += hResp->GetBinContent(i,j);
        Psum[j-1] = sum;
    }

    TH2D* hRespNorm = new TH2D("hRespNorm",
			       " ",
                               nBinsX, binEdgesX.data(),
                               nBinsY, binEdgesY.data());

    for (int i = 1; i <= nBinsX; ++i) {
        for (int j = 1; j <= nBinsY; ++j) {
            double val = hResp->GetBinContent(i,j);
            double norm = (Psum[j-1] != 0.0) ? (val / Psum[j-1] * 100.0) : 0.0;
            hRespNorm->SetBinContent(i, j, norm);
        }
    }

    // orange pastel palette
    Double_t r[]    = {1.00, 1.00, 0.95, 0.85};
    Double_t g[]    = {1.00, 0.85, 0.60, 0.40};
    Double_t b[]    = {1.00, 0.60, 0.20, 0.00};
    Double_t stop[] = {0.00, 0.33, 0.66, 1.00};
    TColor::CreateGradientColorTable(4, stop, r, g, b, 100);
    gStyle->SetPaintTextFormat("5.1f");

    cResp->cd();
    cResp->SetRightMargin(0.15);
    cResp->SetGrid();
    //    hRespNorm->SetTitle("Response Matrix (normalized to 100% per truth bin)");
    hRespNorm->GetXaxis()->SetTitle("Reconstructed neutrino energy (MeV)");
    hRespNorm->GetYaxis()->SetTitle("True neutrino energy (MeV)");
    hRespNorm->SetMinimum(1.0);
    hRespNorm->SetMarkerSize(1.5); //in origine era 1.0
    hRespNorm->Draw("COLZ TEXT");

    SaveCanvasIfValid(cResp, base, "ResponseNorm");

    f->Close();
}





