/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : l0ex1_sb_single.cc
 * @created     : Monday Dec 02, 2024 22:10:49 CET
 */

#include <iostream>
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TRandom3.h"

void produce_dataset(TH1D* h_data, const double exposure, 
    const double sig_rate, const double sig_mu, const double sig_sigma, 
    const double bkg_rate, const double bkg_tau);

TH1D* make_residuals(const TH1D* h_data, const TF1* fModel, const double range_min, const double range_max);

TH1D* make_pulls(const TH1D* h_data, const TF1* fModel, const double range_min, const double range_max);

int l0ex1_sb_single()
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

  produce_dataset(h_data, exposure, signal_rate, sig_mean, sig_sigma, background_rate, -1.0 / bkg_slope);

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
  fModel->SetParNames("bkg N", "bkg slope", "sig N", "sig #mu", "sig #sigma");
  fModel->SetLineColor(kBlack); fModel->SetLineWidth(2);

  const double range_min = sig_mean - 5*sig_sigma;
  const double range_max = sig_mean + 5*sig_sigma;
  printf("Fitting range: [%g, %g]\n", range_min, range_max);

  // initialize parameters
  double entries = h_data->GetEntries();
  fModel->SetParameters( entries , -0.001, exposure*signal_rate, sig_mean, sig_sigma);

  TFitResultPtr fit_results = h_data->Fit(fModel, "SQ", "", range_min, range_max);

  // --------------------------------------------------------------------
  // Print results
  // --------------------------------------------------------------------
  fit_results->Print("V");

  // --------------------------------------------------------------------
  // Plot results
  // --------------------------------------------------------------------
  TCanvas* c_data = new TCanvas("c_data", "c_data", 800, 600);
  //c_data->Divide(2, 1);
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
  //c_data->cd(2); gPad->SetGrid(); gPad->SetTicks();
  //TH1D* h_pull = make_pulls(h_data, fModel, range_min, range_max);
  //h_pull->Draw("hist");

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
