/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : l0ex1_sb_single.cc
 * @created     : Monday Dec 02, 2024 22:10:49 CET
 */

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"

/**
 * @struct ParameterDef_t
 * @brief Structure to store parameter attributes
 */
struct ParameterDef_t {
  double fMean; 
  double fStdDev; 
  double fVal;
  double fError;
  double fTruth;
  TString fName; 
  TString fTitle;

  ParameterDef_t() 
    : fMean(0), fStdDev(0), fVal(0), fError(0), fTruth(0), fName("no_name"), fTitle("no_title") {}

  ParameterDef_t(const ParameterDef_t& p) {
    fMean = p.fMean; 
    fStdDev = p.fStdDev; 
    fVal = p.fVal; 
    fError = p.fError;
    fTruth = p.fTruth;
    fName = p.fName; 
    fTitle = p.fTitle; 
  }

  ParameterDef_t(const TString name, const TString title, const double true_val = 0.0) 
    : fMean(0), fStdDev(0), fVal(0), fError(0), fTruth(true_val), fName(name), fTitle(title) {}
}; 


void produce_dataset(TH1D* h_data, const double exposure, 
    const double sig_rate, const double sig_mu, const double sig_sigma, 
    const double bkg_rate, const double bkg_tau);

TH1D* make_residuals(const TH1D* h_data, const TF1* fModel, const double range_min, const double range_max);

TH1D* make_pulls(const TH1D* h_data, const TF1* fModel, const double range_min, const double range_max);

TH1D* get_par_values(TTree* t, const ParameterDef_t par);

int l0ex1_sb_multi(const int n_experiments = 1)
{
  gStyle->SetOptFit(1111); 
  // --------------------------------------------------------------------
  // Define data PDFs
  // --------------------------------------------------------------------
  // define experimental parameters
  const double exposure = 1e6; // exposure time in seconds
  const double signal_rate = 0.005; // signal rate in Hz
  const double background_rate = 0.05; // background rate in Hz
  // Signal parameters
  const double sig_mean = 511.0; 
  const double sig_sigma = 45.0; 
  // Background parameters
  const double bkg_slope = -1.5e-3; 

  printf("********************************************************************\n");
  printf("Expected signal counts: %g\n", exposure * signal_rate);
  printf("Expected background counts: %g\n", exposure * background_rate);
  printf("********************************************************************\n\n");

  // --------------------------------------------------------------------
  // Create data histogram 
  // --------------------------------------------------------------------
  TH1D* h_data = new TH1D("h_data", "data;Channel;Counts", 256, 0, 2048);
  const double bin_width = h_data->GetBinWidth(10);

  // --------------------------------------------------------------------
  // Define model for the signal, background and data
  // --------------------------------------------------------------------
  TF1* fSignal = new TF1("signal", 
      [bin_width](double* x, double* p) {return bin_width*p[0]*TMath::Gaus(x[0], p[1], p[2], true);},
      0, 2048, 3);
  fSignal->SetParNames("sig N", "sig #mu", "sig #sigma");
  fSignal->SetLineColor(kRed+1); fSignal->SetLineWidth(2);

  TF1* fBackground = new TF1("background", 
      [bin_width](double* x, double* p) {return -bin_width*p[0]*p[1]*exp(x[0]*p[1]);}, 
      0, 2048, 2);
  fBackground->SetParNames("bkg N", "bkg slope");
  fBackground->SetLineColor(kBlue+1); fBackground->SetLineWidth(2); fBackground->SetLineStyle(7);

  TF1* fModel = new TF1("fModel",
      [bin_width, fBackground, fSignal](double* x, double* p) {
        fBackground->SetParameters(p[0], p[1]);
        fSignal->SetParameters(p[2], p[3], p[4]);
        double y = fBackground->Eval(x[0]) + fSignal->Eval(x[0]);
        return y;
      },
      0, 2048, 5); 
  fModel->SetLineColor(kBlack); fModel->SetLineWidth(2);

  std::vector<ParameterDef_t> fit_pars = {
    ParameterDef_t("bkg_norm", "Bkg norm", exposure*background_rate),
    ParameterDef_t("bkg_slope", "Bkg slope", bkg_slope),
    ParameterDef_t("signal_norm", "Sig norm", exposure*signal_rate),
    ParameterDef_t("signal_mean", "Sig #mu", sig_mean),
    ParameterDef_t("signal_sigma", "Sig #sigma", sig_sigma)
  };
  for (int ipar = 0; ipar < fit_pars.size(); ipar++) {
    fModel->SetParName(ipar, fit_pars[ipar].fName);
  }

  const double range_min = sig_mean - 5*sig_sigma;
  const double range_max = sig_mean + 5*sig_sigma;
  const int ibin_min = h_data->FindBin(range_min);
  const int ibin_max = h_data->FindBin(range_max);
  const int ndf = ibin_max - ibin_min - fModel->GetNpar();
  double minFCN = 0.0;
  double pval = 0.0;
  printf("Fitting range: [%g, %g] - %i degrees of freedom\n", range_min, range_max, ndf);

  TFile* file_output = nullptr;
  TTree* t_results = nullptr;
  if (n_experiments > 1) {
    file_output = new TFile("pseudo_data.root", "recreate");
    t_results = new TTree("t_results", "t_results");
    t_results->Branch("minFCN", &minFCN, "minFCN/D");
    t_results->Branch("pval", &pval, "pval/D");
    for (auto& par : fit_pars) {
      t_results->Branch(par.fName, &par.fVal, Form("%s/D", par.fName.Data()));
      t_results->Branch(par.fName + "_err", &par.fError, Form("%s_err/D", par.fName.Data()));
    }
  }

  for (int iexp = 0; iexp < n_experiments; iexp++) {
    h_data->Reset();
    produce_dataset(h_data, exposure, signal_rate, sig_mean, sig_sigma, background_rate, -1.0/bkg_slope);
    
    // initialize parameters
    double entries = h_data->GetEntries();
    fModel->SetParameters( entries , -0.001, exposure*signal_rate, sig_mean, sig_sigma);

    TFitResultPtr fit_results = h_data->Fit(fModel, "SQ", "", range_min, range_max);

    minFCN = fit_results->MinFcnValue();
    pval = fit_results->Prob();
    for (int ipar = 0; ipar < fit_pars.size(); ipar++) {
      fit_pars[ipar].fVal = fModel->GetParameter(ipar);
      fit_pars[ipar].fError = fModel->GetParError(ipar);
    }

    if (n_experiments > 1) {
      if (fit_results->Status() == 0) t_results->Fill();
    }
    else {
      fit_results->Print("V");
    }
  }

  TCanvas* c_data = new TCanvas("c_data", "c_data", 800, 600);
  if (n_experiments == 1) {
    TVirtualPad* pad_data = c_data->cd(0);
    pad_data->Divide(1, 2);
    pad_data->cd(1); gPad->SetGrid(); gPad->SetTicks();
    h_data->GetXaxis()->SetRangeUser(range_min, range_max);
    h_data->GetYaxis()->SetRangeUser(0, 1.2*h_data->GetMaximum());
    h_data->Draw();
    fSignal->SetParameters(fModel->GetParameter(2), fModel->GetParameter(3), fModel->GetParameter(4));
    fBackground->SetParameters(fModel->GetParameter(0), fModel->GetParameter(1));
    std::vector<TF1*> func_list = {fSignal, fBackground, fModel};
    for (auto& f : func_list) {
      f->SetRange(range_min, range_max);
      f->Draw("same");
    }
    pad_data->cd(2); gPad->SetGrid(); gPad->SetTicks();
    TH1D* h_res = make_residuals(h_data, fModel, range_min, range_max);
    h_res->Draw("p");
  }
  else {
    c_data->DivideSquare(fModel->GetNpar());
    int ipar = 0; 
    for (const auto& par : fit_pars) {
      TH1D* h_par = get_par_values(t_results, par);
      c_data->cd(ipar+1);
      h_par->Draw();
      h_par->Fit("gaus", "Q");
      TF1* fgaus = h_par->GetFunction("gaus");
      const double bias = (par.fTruth - fgaus->GetParameter(1)) / fgaus->GetParError(1);
      printf("%s: bias = %g [%g%%]\n", par.fName.Data(), bias, (par.fTruth - fgaus->GetParameter(1)) / par.fTruth * 100);
      ipar++;
    }

    TCanvas* c_gof = new TCanvas("c_gof", "c_gof", 800, 600);
    c_gof->Divide(2, 1);
    c_gof->cd(1);
    TH1D* h_gof = new TH1D("h_gof", "h_gof", 100, ndf - 4*sqrt(2*ndf), ndf + 4*sqrt(2*ndf));
    t_results->Draw("minFCN>>h_gof", "", "goff");
    h_gof->DrawClone();
    TF1* fchi2 = new TF1("fchi2", Form("[0]*ROOT::Math::chisquared_pdf(x, %i)", ndf), ndf - 3*sqrt(2*ndf), 4*sqrt(2*ndf));
    fchi2->SetParameter(0, n_experiments);
    h_gof->Fit(fchi2, "QL");

    c_gof->cd(2);
    TH1D* h_pval = new TH1D("h_pval", "h_pval", 100, 0, 1);
    t_results->Draw("pval>>h_pval", "", "goff");
    h_pval->DrawClone();
    h_pval->Fit("pol0", "QL");

    //Write the results to the output file
    t_results->Write();
  }

  return 0;
}

void produce_dataset(TH1D* h_data, const double exposure, 
    const double sig_rate, const double sig_mu, const double sig_sigma,
    const double bkg_rate, const double bkg_tau) 
{
  const int n_ev_expc_signal = exposure * sig_rate;
  const int n_ev_expc_background = exposure * bkg_rate;

  for (int i = 0; i < n_ev_expc_signal; i++) {
    h_data->Fill( gRandom->Gaus(sig_mu, sig_sigma) );
  }
  for (int i = 0; i < n_ev_expc_background; i++) {
    h_data->Fill( gRandom->Exp( bkg_tau ) );
  }

  return; 
}

TH1D* make_residuals(const TH1D* h_data, const TF1* fModel, const double xmin, const double xmax) 
{
  TH1D* h_res = static_cast<TH1D*>(h_data->Clone("h_res"));
  h_res->Reset();

  double res_max = 0;

  for (int ibin = 1; ibin <= h_data->GetNbinsX(); ibin++) {
    double x = h_data->GetBinCenter(ibin);
    if (x < xmin || x > xmax) continue;
    double y_data = h_data->GetBinContent(ibin);
    double y_model = fModel->Eval(x);
    double y_res = (y_data - y_model) / TMath::Sqrt(y_model);
    if (TMath::Abs(y_res) > res_max) res_max = TMath::Abs(y_res);
    h_res->SetBinContent(ibin, y_res);
  }

  h_res->SetYTitle("(data - model) / #sqrt{model}");
  h_res->GetYaxis()->SetRangeUser(-1.1*res_max, 1.1*res_max);
  h_res->SetMarkerStyle(20);

  return h_res;
}

TH1D* make_pulls(const TH1D* h_data, const TF1* fModel, const double xmin, const double xmax) 
{
  TH1D* h_pull = new TH1D("h_pull", "pull distribution", 50, -5, 5);

  for (int ibin = 1; ibin <= h_data->GetNbinsX(); ibin++) {
    double x = h_data->GetBinCenter(ibin);
    if (x < xmin || x > xmax) continue;
    double y_data = h_data->GetBinContent(ibin);
    double y_model = fModel->Eval(x);
    double y_res = (y_data - y_model) / TMath::Sqrt(y_model);
    h_pull->Fill(y_res);
  }
  return h_pull;
}

TH1D* get_par_values(TTree* t, const ParameterDef_t par) {
  t->Draw(Form("%s>>htmp", par.fName.Data()), "", "goff"); 
  TH1D* htmp = (TH1D*)gDirectory->Get("htmp");
  double mean = htmp->GetMean();
  double rms = htmp->GetRMS();
  delete htmp;
  TString hist_name = Form("h%s", par.fName.Data());
  TH1D* h = new TH1D(hist_name, "h", 50, mean - 6*rms, mean + 6*rms);
  t->Draw(Form("%s>>%s", par.fName.Data(), hist_name.Data()), "", "goff");
  return h; 
}



