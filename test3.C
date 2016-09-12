
void test3(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gks;

  int ret=gks.setInputHisto( "FF_corr_QCD_MCsum_noGen_nonclosure.root" , "nonclosure" );
  TH1D *h=gks.returnInputHisto();

  if ( ret != 0 ) return;
  gks.set_doWeights(1);
  gks.set_doIgnoreZeroBins(1);
  //  gks.set_kernelDistance( "err" );
  gks.set_doWidthInBins(1);
  gks.set_doErrors(1);
  gks.set_lastBinFrom( 120 );
  Double_t fitWidth=1.5;
  gks.setWidth(fitWidth);

  gks.getSmoothHisto();
  TH1D *h2=gks.returnSmoothedHisto();
  //  h2->Write();
  gks.set_doErrors(1);
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

}
