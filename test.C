
void test(){

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gsk;

  int ret=gsk.createTestHisto();
  //  int ret=gsk.setInputHisto( "/afs/hephy.at/work/m/mflechl/git/HephyHiggs/FakeFactor2015/fakefactor/data_mutau/sys_TT_J_nonclosure_mvis.root" , "sys" );
  if ( ret != 0 ) return;

  TH1D *h=gsk.returnInputHisto();
  gsk.setWidth( 5* h->GetBinWidth(1) );
  
  gsk.getSmoothHisto();
  TH1D *h2=gsk.returnSmoothedHisto();

  gsk.set_doErrors( 1 );
  //  gsk.set_kernelDistance( "log" );
  gsk.set_kernelDistance( "lin" );
  gsk.getContSmoothHisto();
  TGraphAsymmErrors *g=   gsk.returnSmoothedGraph();
  TGraphAsymmErrors *gerr=gsk.returnSmoothedGraphErr();

  TCanvas *c1=new TCanvas();
  h->Draw("E");
  h2->SetLineColor(kRed);
  h2->Draw("HIST C same");
  gPad->Update();

  TCanvas *c2=new TCanvas();
  c2->Divide(2);
  c2->cd(1);
  g->Draw("AP");
  c2->cd(2);
  gerr->Draw("AP");

}
