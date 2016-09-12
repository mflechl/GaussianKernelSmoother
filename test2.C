
void test2(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gsk;

  //  int ret=gsk.createTestHisto();
  //  int ret=gsk.setInputHisto( "input_nonclosure_2bins.root" , "nonclosure" );
  int ret=gsk.setInputHisto( "FF_corr_Wjets_MCsum_noGen_nonclosure.root" , "nonclosure" );
  //  int ret=gsk.setInputHisto( "FF_corr_QCD_MCsum_noGen_nonclosure.root" , "nonclosure" );
  if ( ret != 0 ) return;
  TH1D *h=gsk.returnInputHisto();

  gsk.set_doWeights(1);
  gsk.set_doIgnoreZeroBins(1);
  gsk.set_kernelDistance( "lin" );

  //option 1: give smoothing scale parameter as multiples of bins (->i.e. non-const for non-equidistant histos!)
  //  gsk.set_doWidthInBins(1);
  //  gsk.setWidth( 2 );

  //option 2: give smoothing scale parameter directly
  //  gsk.set_doWidthInBins(0);
  //  gsk.setWidth( 2* h->GetBinWidth(1) );

  //option 3: log scaling
  //  gsk.set_kernelDistance( "log" );
  //  gsk.setWidth( 0.2 );

  //option 4: error scaling
  gsk.setWidth( 0.08 );
  //gsk.setWidth( 0.22 );
  gsk.set_kernelDistance( "err" );
  gsk.set_doWidthInBins(0);
  gsk.set_doErrors( 1 );
  //gsk.setWidth( 0.03 );

  gsk.getSmoothHisto();
  TH1D *h2=gsk.returnSmoothedHisto();

  gsk.set_doErrors(1);
  gsk.getContSmoothHisto();
  TGraphAsymmErrors *g=   gsk.returnSmoothedGraph();
  //X  TGraphAsymmErrors *gerr=gsk.returnSmoothedGraphErr();


  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  g->Draw("same LP");
  h->Draw("E same");

}
