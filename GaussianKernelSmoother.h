#ifndef __GaussianKernelSmoother__
#define __GaussianKernelSmoother__

//#include <algorithm>
#include <stdlib.h>
#include<math.h>
#include<iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TString.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;

class GaussianKernelSmoother
{

  public:
  GaussianKernelSmoother();
  void  computeSmoothHisto();
  void  createTestHisto();
  void setWidth(const double w){ this->width = w; }
  TH1D* returnSmoothedHisto(){ return this->h_out; }
  TH1D* returnInputHisto(){ return this->h_in; }
  int setInputHisto(TString fname , TString hname ){ 
    TFile *f=new TFile( fname );
    this->h_in=(TH1D*) f->Get( hname );
    if ( this->h_in ) return 0;
    else{ cerr << "ERROR: " << hname << " not found in file " << fname << endl; return 1; }
  }

  private:
  double getSmoothedValue(double x);
  //  TH1D* fluctuateHisto();

  TH1D *h_in;
  TH1D *h_out;
  TH1D *weights;

  double width;
  int doErrors;
  TRandom3 rand;
  //        self.weights = None
  //	  self.rescaling = lambda x:x
  //        self.invertRescaling = self.rescaling

  public:
  //    virtual ~GaussianKernelSmoother();

  public:
    ClassDef(GaussianKernelSmoother,0)

};

#endif
