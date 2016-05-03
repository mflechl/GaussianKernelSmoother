#include "GaussianKernelSmoother.h"

ClassImp(GaussianKernelSmoother)

GaussianKernelSmoother::GaussianKernelSmoother(){
  this->width=1.;
  this->doErrors=0;
}

/*
TH1D* GaussianKernelSmoother::fluctuateHisto(){
  TH1D *h = this->h_in->Clone(this->h_in->GetName()+'_toy');
  Int_t nbins = this->h_in->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++){
    double val = this->h_in->GetBinContent(ib);
    double err = this->h_in->GetBinError(ib);
    double newval = this->random.Gaus(val,err);
    h->SetBinContent(b,newval);
  }
  return h;
}
*/

void GaussianKernelSmoother::createTestHisto(){

  TF1 *p1 = new TF1("p1","[0]+x*[1]+x*x*[2]");
  p1->SetParameters(50., 10., 0.);
  //p1->SetParameters(0, 5.);
  //p1->SetParameters(1,10.);
  //p1->SetParameters(2,-1.);

  std::cout << p1->Eval(-5) << "\t" << p1->Eval(0) << "\t" << p1->Eval(5) << std::endl;

  this->h_in = new TH1D("","h_in",100,-5,5);
  this->h_in->FillRandom("p1",1000);

}

double GaussianKernelSmoother::getSmoothedValue(double x){

  double sumw = 0.;
  double sumwy = 0. ;
  Int_t nbins = this->h_in->GetNbinsX();

  int firstBinWithContent=0;
  for (int ib=1; ib<=nbins; ib++){
    double yi = this->h_in->GetBinContent(ib);
    if ( !firstBinWithContent && fabs(yi)<1e-8 ) continue; else{ firstBinWithContent=ib; } //skip all bins until the first with non-zero content
    double xi = this->h_in->GetXaxis()->GetBinCenter(ib);

    double dx = ( x - xi ) / this->width;

    double wi = TMath::Gaus(dx);
    //    if (this->weights) wi *= self.weights.GetBinContent(b);
    sumw += wi;
    sumwy += wi*yi;
  }
   
  double val = 0.;
  if (sumw>0.) val = sumwy/sumw;
    
  return val;
}

void GaussianKernelSmoother::computeSmoothHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: GaussianKernelSmoother::computeSmoothHisto: input histo does not exist." << std::endl;

  this->h_out = (TH1D*)this->h_in->Clone(this->h_in->GetName()+(TString)"_smoothed");
  this->h_out->SetDirectory(0);

  Int_t nbins = this->h_out->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++){
    double x = this->h_out->GetBinCenter(ib);
    double smooth_val = this->getSmoothedValue(x);
    this->h_out->SetBinContent(ib,smooth_val);
  }
}
