/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : p_value_test.C
 * @created     : Sunday Dec 01, 2024 09:48:23 CET
 */

#include <iostream>
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "Math/PdfFuncMathCore.h"
#include "Math/PdfFuncMathMore.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"


/**
 * @brief A simple macro to produce the distribution of a test-statistics 
 * and of fit p-values
 *
 * @param n_experiments number of pseudo-experiments to perform
 * @return 0 if successful
 */
int p_value_test(const int n_experiments = 1e4)
{
  // define the gaussian parameters for producing the data
  const double gauss_mu = 150; 
  const double gauss_sigma = 20;

  // Create a histogram to store the synthetic data
  TH1D* h = new TH1D("h", "My first histogram;Channel;Counts", 2048, 0, 2048);

  // Define the model function
  TF1* fGaus = new TF1("fGaus", "gausn", 0, 2048);

  // Set the number of entries in a single pseudo-experiment
  const int n_events = 200000;

  // Set the fit range
  const double xmin = TMath::Max(h->GetXaxis()->GetXmin(), gauss_mu - 3*gauss_sigma);
  const double xmax = TMath::Min(h->GetXaxis()->GetXmax(), gauss_mu + 3*gauss_sigma);
  
  // Estimate the ndf for this fit
  const int ibin0 = h->FindBin(xmin);
  const int ibin1 = h->FindBin(xmax);
  const int expc_ndf = ibin1 - ibin0 + 1 - 3;
  printf("Expected number of degrees of freedom: %i\n", expc_ndf);

  // Create histograms to store the fit parameters
  TH1D* h_par_val[3];
  TH1D* h_par_err[3];
  h_par_val[0] = new TH1D("h_const", "Constant distribution;Constant;Counts", 100, 
      n_events - sqrt(n_events), n_events + sqrt(n_events));
  h_par_val[1] = new TH1D("h_mean", "Mean distribution;Mean;Counts", 100, 
      gauss_mu - gauss_mu/sqrt(n_events), gauss_mu + gauss_mu/sqrt(n_events));
  h_par_val[2] = new TH1D("h_sigma", "Sigma distribution;Sigma;Counts", 100, 
      gauss_sigma - 3.0*gauss_sigma/sqrt(n_events), gauss_sigma + 3.0*gauss_sigma/sqrt(n_events));
  h_par_err[0] = new TH1D("h_const_err", "Constant error distribution;Constant error;Counts",
      100, 0.001*n_events, 0.01*n_events);
  h_par_err[1] = new TH1D("h_mean_err", "Mean error distribution;Mean error;Counts", 
      100, 0.001*gauss_mu, 0.01*gauss_mu);
  h_par_err[2] = new TH1D("h_sigma_err", "Sigma error distribution;Sigma error;Counts", 
      100, 0.001*gauss_sigma, 0.01*gauss_sigma);

  TH1D* h_chi2 = new TH1D("h_chi2", "#chi^{2} distribution;#chi^{2};Counts", 
      100, expc_ndf-3*sqrt(2*expc_ndf), expc_ndf+3*sqrt(2*expc_ndf));
  TH1D* h_pval = new TH1D("h_pval", "p-value distribution;p-value;Counts", 100, 0, 1);

  // Loop over the pseudo-experiments
  for (int iexp = 0; iexp < n_experiments; iexp++) {
    h->Reset();
    // Fill the histogram with synthetic data
    for (int i = 0; i < n_events; i++) {
      h->Fill(gRandom->Gaus(gauss_mu, gauss_sigma));
    }
    // Fit the histogram
    fGaus->SetParameters(200, 150, 20);
    TFitResultPtr result = h->Fit(fGaus, "S0Q", "", xmin, xmax);
    if (result->Status() != 0) {
      std::cout << "Fit failed with status " << result->Status() << std::endl;
      continue;
    }

    // Collect the fit results
    for (int i=0; i<3; i++) {
      h_par_val[i]->Fill(result->Parameter(i));
      h_par_err[i]->Fill(result->ParError(i));
    }

    const double chi2 = result->Chi2();
    const double ndf = result->Ndf();
    const double pval = TMath::Prob(chi2, ndf);
    h_chi2->Fill(chi2);
    h_pval->Fill(pval);
  }

  // -----------------------------------------------------------
  // Plot the distributions
  // -----------------------------------------------------------

  // Activate the display of fit results on the stat box
  gStyle->SetOptFit(1111);

  // Draw the fit best-estimates of the parameters and their errors
  TCanvas* c_par = new TCanvas("c_par", "Parameter canvas", 800, 600);
  c_par->Divide(3, 2);
  for (int i=0; i<3; i++) {
    c_par->cd(i+1); h_par_val[i]->Draw();
    h_par_val[i]->Fit("gaus");
    c_par->cd(i+4); h_par_err[i]->Draw();
  }

  // define chi2 pdf
  TF1* fChi2 = new TF1("fChi2", "[0]*ROOT::Math::chisquared_pdf(x, [1])", 60, 260);
  fChi2->SetParameters(n_experiments, 160);
  fChi2->SetParName(0, "Constant");
  fChi2->SetParName(1, "ndf");
  h_chi2->Fit(fChi2, "0");

  // Draw the chi2 distribution and fit it with its PDF
  TCanvas* c_chi2 = new TCanvas("c_chi2", "chi2 canvas", 800, 600);
  h_chi2->Draw();
  fChi2->Draw("same");

  // Draw the p-value distribution and fit it with a constant
  TCanvas* c_pval = new TCanvas("c_pval", "p-value canvas", 800, 600);
  h_pval->GetYaxis()->SetRangeUser(0., h_pval->GetMaximum()*1.1);
  h_pval->Draw();
  h_pval->Fit("pol0");


  return 0;
}

