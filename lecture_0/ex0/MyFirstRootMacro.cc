#include <iostream>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"


/**
 * @brief A function that checks the bin errors of a histogram
 *
 * @param h input histogram
 * @return true if all bin errors are correct, false otherwise
 */
bool check_errors(const TH1D* h)
{
  bool bin_check = true;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    const double err = h->GetBinError(i);
    const double content = h->GetBinContent(i);
    
    if ( err != sqrt(content) ) {
      std::cerr << "Error: bin " << i << " has an error of " << err << " and a content of " << content << std::endl;
      bin_check = false;
    }
  }
  return bin_check;
}

/**
 * @brief A simple implementation of a Gaussian function
 *
 * @param x input variable
 * @param par Gaussian parameters (constant, mean, sigma)
 * @return the value of the Gaussian function at x[0]
 */
double my_gaus(double* x, double* par) 
{
  double constant = par[0]; 
  double mean = par[1];
  double sigma = par[2];
  
  double delta = (x[0] - mean) / sigma;
  double y = constant * exp(-0.5 * delta * delta);
  
  return y;
}

/**
 * @brief Fill a histogram with data from a text file
 *
 * @param h target histogram
 * @param filename txt file path
 * @return 0 if successful
 */
int fill_hist_from_txt(TH1D* h, const TString filename);

/**
 * @brief A simple ROOT macro that reads data from a file and fits it with a Gaussian
 *
 * @param filename input file name
 * @return 0 if successful
 */
int MyFirstRootMacro(const TString filename = "test_data.txt") 
{
  // --------------------------------------------------------------------
  // Build the data histogram
  // --------------------------------------------------------------------

  // Create a histogram to store the data
  TH1D* h = new TH1D("h", "My first histogram;Channel;Counts", 2048, 0, 2048);

  // Fill the histogram with data from the file
  int fill_status = fill_hist_from_txt(h, filename);
  if (fill_status != 0) {
    std::cout << "Error: fill_hist_from_txt failed with status " << fill_status << std::endl;
    return 1;
  }

  // Check bin errors
  bool bin_check = check_errors(h); 

  // Draw the histogram
  TCanvas* c = new TCanvas("c", "My first canvas", 800, 600);
  h->Draw();

  // --------------------------------------------------------------------
  // Define model 
  // --------------------------------------------------------------------
  // Create a fit model
  TF1* f = nullptr;

  // A quick (and partial) look at your options for definining a function
  f = new TF1("f", "gaus", 0, 2048); // Built-in Gaussian function
  //f = new TF1("f", "[0]*exp(-0.5*TMath::Sq((x-[1])/[2]))", 0, 2048); // Inline function definition
  //f = new TF1("f", my_gaus, 0, 2048, 3); // External function definition
  //f = new TF1("f", [](double* x, double* par) { // Lambda function definition
      //double constant = par[0]; 
      //double mean = par[1];
      //double sigma = par[2];
      //double y = constant * exp(-0.5 * TMath::Sq((x[0] - mean) / sigma));
      //return y;
      //}, 0, 2048, 3);

  f->SetParName(0, "Constant");
  f->SetParName(1, "#mu");
  f->SetParName(2, "#sigma");

  // Initialize the parameters
  f->SetParameter(0, h->GetMaximum());
  f->SetParameter(1, h->GetMean());
  f->SetParameter(2, h->GetRMS());
  
  // --------------------------------------------------------------------
  // Fit the data
  // --------------------------------------------------------------------
  // Set the fit range
  double xmin = h->GetMean() - 3*h->GetRMS();
  double xmax = h->GetMean() + 3*h->GetRMS();

  // Fit the histogram
  TFitResultPtr result = h->Fit(f, "S", "", xmin, xmax);

  // retrieve the χ² value and the ndf
  const double chi2 = result->Chi2();
  const double ndf = result->Ndf();

  // compute the p-value
  double p_value = TMath::Prob(chi2, ndf);
  std::cout << "p-value: " << p_value << std::endl;

  if (p_value < 0.05) {
    std::cout << "Warning: the p-value is less than 0.05" << std::endl;
    return 1;
  }

  return 0; 
}

int fill_hist_from_txt(TH1D* h, const TString filename)
{
  // Open and read the input file in a stream
  std::ifstream inputfile;
  inputfile.open(filename);
  if (!inputfile.is_open()) {
    std::cerr << "Error: file " << filename << " not found." << std::endl;
    return 1;
  }
  std::string line;
  std::getline(inputfile, line); // Skip the header
                                 
  // Read the data from the file
  int channel = 0;
  double counts = 666;

  while ( std::getline(inputfile, line) ) {
    // Read the data from the file
    std::istringstream sstream(line);
    sstream >> channel >> counts;
    h->SetBinContent(channel, counts);
  }

  // Close the input file
  inputfile.close();

  return 0;
}

