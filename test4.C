#include "GaussianKernelSmoother.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
//#include "TF1.h"
#include "TH1D.h"
#include <stdlib.h>


void test4(){
  gStyle->SetOptStat(0);

  //  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gks;

  //  int ret=gks.setInputHisto( "FF_corr_TT_MC_noGen_nonclosure.root" , "nonclosure_fit" );
  //  int ret=gks.setInputHisto( "FF_corr_QCD_MCsum_noGen_nonclosure.root" , "nonclosure_fit" );
  //  int ret=gks.setInputHisto( "FF_corr_Wjets_MCsum_noGen_nonclosure.root" , "nonclosure_fit" );
  int ret=gks.setInputHisto( "FF_corr_QCD_MCsum_noGen_OSSScorr.root" , "OSSS_corr" );
  TH1D *h=gks.returnInputHisto();

  if ( ret != 0 ) return;

  gks.set_doWeights(1);
  gks.set_doIgnoreZeroBins(0);

  //  gks.set_kernelDistance( "err" );
  //  gks.setWidth(0.01); //fitWidth
  //  gks.set_doWidthInBins(1);

  gks.setWidth(2.0); //fitWidth
  gks.set_doWidthInBins(1);
  gks.set_widthInBins_sf(1.115);
  //  gks.setWidth(20); //fitWidth
  //  gks.set_doWidthInBins(0);

  gks.set_doErrors(1);

  //if(mode & _QCD) 
  //  gks.set_lastBinFrom(150);

  gks.getSmoothHisto();
  TH1D *h2=gks.returnSmoothedHisto();
  //  h2->Write();
  //  gks.set_doErrors(1);
  gks.getContSmoothHisto();
  TGraphAsymmErrors *g=   gks.returnSmoothedGraph();
  //  g->SetTitle("nonclosure"+sample);
  //  g->SetName("nonclosure"+sample);
  //  g->Write();


  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  h->GetYaxis()->SetRangeUser(0.4,2.0);
  g->Draw("same LP");
  h->Draw("E same");

  //  h2->Draw();
  gPad->SaveAs("test4.png");

}

#ifndef __CINT__
int main(int argc, char* argv[]) {
  //  int channel = argc > 1 ? atoi(argv[1]) : 0;
  test4();
}
#endif
