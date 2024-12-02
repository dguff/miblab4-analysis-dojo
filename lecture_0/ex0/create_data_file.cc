#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TRandom3.h"

void create_data_file(const TString outputname) 
{
  // Open the output file in a stream
  std::ofstream outputfile;
  outputfile.open(outputname);
  if (!outputfile.is_open()) {
    std::cerr << "Error: file " << outputname << " not found." << std::endl;
    return;
  }

  // Create a histogram to store the syntentic data
  TH1D* h = new TH1D("h", "My first histogram;Channel;Counts", 2048, 0, 2048);

  // Fill the histogram with synthetic data
  const int nentries = 10000;

  for (int i = 0; i < nentries; i++) {
    h->Fill(gRandom->Gaus(150, 20));
  }

  // Write the channel and counts to the output file
  outputfile << "Channel\tCounts\n";

  for (int i=1; i<=h->GetNbinsX(); i++) {
    outputfile << i << "\t" << h->GetBinContent(i) << std::endl;
  }

  outputfile.close();
  
  return;
}
