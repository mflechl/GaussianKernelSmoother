
void test2(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gks;

  //  int ret=gks.createTestHisto();
  //  int ret=gks.setInputHisto( "input_nonclosure_2bins.root" , "nonclosure" );
  int ret=gks.setInputHisto( "FF_corr_Wjets_MCsum_noGen_nonclosure.root" , "nonclosure" );
  //  int ret=gks.setInputHisto( "FF_corr_QCD_MCsum_noGen_nonclosure.root" , "nonclosure" );
  if ( ret != 0 ) return;
  TH1D *h=gks.returnInputHisto();

  gks.set_doWeights(1);
  gks.set_doIgnoreZeroBins(1);
  gks.set_kernelDistance( "lin" );

  //option 1: give smoothing scale parameter as multiples of bins (->i.e. non-const for non-equidistant histos!)
  //  gks.set_doWidthInBins(1);
  //  gks.setWidth( 2 );

  //option 2: give smoothing scale parameter directly
  //  gks.set_doWidthInBins(0);
  //  gks.setWidth( 2* h->GetBinWidth(1) );

  //option 3: log scaling
  //  gks.set_kernelDistance( "log" );
  //  gks.setWidth( 0.2 );

  //option 4: error scaling
  gks.setWidth( 0.08 );
  //gks.setWidth( 0.22 );
  gks.set_kernelDistance( "err" );
  gks.set_doWidthInBins(0);
  gks.set_doErrors( 1 );
  //gks.setWidth( 0.03 );

  gks.getSmoothHisto();
  TH1D *h2=gks.returnSmoothedHisto();

  gks.set_doErrors(1);
  gks.getContSmoothHisto();
  TGraphAsymmErrors *g=   gks.returnSmoothedGraph();
  //X  TGraphAsymmErrors *gerr=gks.returnSmoothedGraphErr();


  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  g->Draw("same LP");
  h->Draw("E same");

}
