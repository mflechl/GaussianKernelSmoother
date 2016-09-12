
void test3(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gsk;

  int ret=gsk.setInputHisto( "FF_corr_QCD_MCsum_noGen_nonclosure.root" , "nonclosure" );
  TH1D *h=gsk.returnInputHisto();

  if ( ret != 0 ) return;
  gsk.set_doWeights(1);
  gsk.set_doIgnoreZeroBins(1);
  //  gsk.set_kernelDistance( "err" );
  gsk.set_doWidthInBins(1);
  gsk.set_doErrors(1);
  gsk.set_lastBinFrom( 120 );
  Double_t fitWidth=1.5;
  gsk.setWidth(fitWidth);

  gsk.getSmoothHisto();
  TH1D *h2=gsk.returnSmoothedHisto();
  //  h2->Write();
  gsk.set_doErrors(1);
  gsk.getContSmoothHisto();
  TGraphAsymmErrors *g=   gsk.returnSmoothedGraph();
  //  g->SetTitle("nonclosure"+sample);
  //  g->SetName("nonclosure"+sample);
  //  g->Write();


  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  h->GetYaxis()->SetRangeUser(0.4,2.0);
  g->Draw("same LP");
  h->Draw("E same");

}
