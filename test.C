
void test(){

  gROOT->ProcessLine(".L GaussianKernelSmoother.C+");
  GaussianKernelSmoother gsk;

  //  gsk.createTestHisto();
  int ret=gsk.setInputHisto( "../fakefactor/data_mutau/sys_TT_J_nonclosure_mvis.root" , "sys" );
  if ( ret != 0 ) return;

  TH1D *h=gsk.returnInputHisto();
  gsk.setWidth( 2* h->GetBinWidth(1) );
  
  gsk.computeSmoothHisto();

  TH1D *h2=gsk.returnSmoothedHisto();

  h->Draw("E");
  h2->SetLineColor(kRed);
  h2->Draw("HIST C same");
  gPad->Update();


}
