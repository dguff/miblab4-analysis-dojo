/**
 * @author      : Daniele Guffanti (daniele.guffanti@mib.infn.it)
 * @file        : root_io.C
 * @created     : Wednesday Nov 27, 2024 19:09:55 CET
 */

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"

#include "TRandom3.h"

int root_io()
{
  TH1D* hh = new TH1D("my_hist", "histogram", 100, -5, 5);
  for (int i = 0; i < 10000; i++) {
    hh->Fill(gRandom->Gaus(0, 1));
  }

  TFile* f = new TFile("output.root", "RECREATE");
  hh->Write();
  f->Close();
    
  return 0;
}

