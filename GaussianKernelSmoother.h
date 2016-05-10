#ifndef __GaussianKernelSmoother__
#define __GaussianKernelSmoother__

//#include <algorithm>
#include <stdlib.h>
#include<math.h>
#include<iostream>

#include "TMath.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;

class GaussianKernelSmoother
{

  public:
  GaussianKernelSmoother();
  void set_doErrors( int m_doErrors ){ this->doErrors = m_doErrors; }
  void set_doWeights( int m_doWeights ){ this->doWeights = m_doWeights; if (m_doWeights) h_w=this->makeWeights( this->h_in ); }
  void set_kernelDistance( TString m_kernelDistance ){ this->kernelDistance = m_kernelDistance; }
  void getSmoothHisto();
  void getContSmoothHisto();
  void createTestHisto();
  void setWidth(const double w){ this->width = w; }
  TH1D* returnSmoothedHisto(){ return this->h_out; }
  TH1D* returnInputHisto(){ return this->h_in; }
  TGraphAsymmErrors* returnSmoothedGraph(){ return this->g_out; }
  TGraph* returnSmoothedGraphErr(){ return this->g_out_err; }
  int setInputHisto(TString fname , TString hname ){ 
    TFile *f=new TFile( fname );
    this->h_in=(TH1D*) f->Get( hname );
    if ( this->h_in ) return 0;
    else{ cerr << "ERROR: " << hname << " not found in file " << fname << endl; return 1; }
  }

  private:
  double getSmoothedValue(TH1D* m_h , const double x);
  TH1D* makeWeights( TH1D* h);
  double rescaling( double val );
  double invertRescaling( double val );
  TH1D* fluctuateHisto();
  double std_dev(std::vector<double> v );

  TH1D *h_in;
  TH1D *h_out;
  TH1D *h_w;
  TGraphAsymmErrors *g_out;
  TGraph *g_out_err;
  TH1D *weights;

  double width;
  int doErrors;
  int doWeights;
  TString kernelDistance;

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
