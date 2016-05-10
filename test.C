
void test(){
  gStyle->SetOptStat(0);

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gsk;

  int ret=gsk.createTestHisto();
  //  int ret=gsk.setInputHisto( "/afs/hephy.at/work/m/mflechl/git/HephyHiggs/FakeFactor2015/fakefactor/data_mutau/sys_TT_J_nonclosure_mvis.root" , "sys" );
  if ( ret != 0 ) return;

  gsk.set_doWeights(1);

  TH1D *h=gsk.returnInputHisto();
  gsk.setWidth( 5* h->GetBinWidth(1) );

  gsk.getSmoothHisto();
  TH1D *h2=gsk.returnSmoothedHisto();

  //  gsk.set_kernelDistance( "log" );
  gsk.set_kernelDistance( "lin" );
  gsk.set_doErrors(1);
  gsk.getContSmoothHisto();
  TGraphAsymmErrors *g=   gsk.returnSmoothedGraph();
  //X  TGraphAsymmErrors *gerr=gsk.returnSmoothedGraphErr();


  TCanvas *c2=new TCanvas();
  h->Draw("E");     //to set the axis
  g->Draw("same LP");
  h->Draw("E same");

}
