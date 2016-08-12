#include "GaussianKernelSmoother.h"

ClassImp(GaussianKernelSmoother)

GaussianKernelSmoother::GaussianKernelSmoother(){
  this->doIgnoreZeroBins=1.;
  this->doWidthInBins=0.;
  this->width=1.;
  this->doErrors=0;
  this->doWeights=0;
  this->kernelDistance="lin";
}

TH1D* GaussianKernelSmoother::fluctuateHisto(){
  TH1D *h = (TH1D*)this->h_in->Clone((TString)this->h_in->GetName()+"_toy");
  Int_t nbins = this->h_in->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++){
    double val = this->h_in->GetBinContent(ib);
    double err = this->h_in->GetBinError(ib);
    double newval = this->rand.Gaus(val,err);
    h->SetBinContent(ib,newval);
  }
  return h;
}

double GaussianKernelSmoother::getSmoothedValue(TH1D* m_h , const double x){

  double sumw = 0.;
  double sumwy = 0. ;
  Int_t nbins = m_h->GetNbinsX();

  double m_width=this->width;

  if ( this->doWidthInBins ){
    m_width=0;
    double wsum=0;

    double bin_width=0;
    //    bin_width=( m_h->GetBinLowEdge(nbins+1)-m_h->GetBinLowEdge(1) ) / nbins; //option 1: average bin width
    for (int ib=1; ib<=nbins; ib++){
      if ( bin_width<m_h->GetBinWidth(ib) ) bin_width=m_h->GetBinWidth(ib);  //option 2: max bin width
      bin_width/=2;                                                          //option 2b: half max bin width
    }

    //    int m_bin=m_h->FindBin(x);

    for (int ib=1; ib<=nbins; ib++){
      double w=TMath::Gaus( fabs( x-m_h->GetBinCenter(ib) )/(bin_width) );
      m_width+=m_h->GetBinWidth(ib)* w;
      wsum+=w;
    }
    m_width/=wsum;
    m_width*=this->width;
    //    cout << x << " " << m_width << " " << wsum << " " << bin_width*nbins << " " <<bin_width << endl;
  }


  //int firstBinWithContent=0;
  for (int ib=1; ib<=nbins; ib++){
    double yi = m_h->GetBinContent(ib);
    double ei = m_h->GetBinError(ib);
    //  if ( !firstBinWithContent && fabs(yi)<1e-8 ) continue; else{ firstBinWithContent=ib; } //skip all bins until the first with non-zero content
    if ( fabs(yi)<1e-8 && ei<1e-8 ) continue;
    double xi = m_h->GetXaxis()->GetBinCenter(ib);

    //    double dx = ( x - xi ) / this->width;
    double dx = ( this->rescaling(x) - this->rescaling(xi) ) / m_width;

    double wi = TMath::Gaus(dx);
    if (this->doWeights) wi *= this->h_w->GetBinContent(ib);
    sumw += wi;
    sumwy += wi*yi;
  }
   
  double val = 0.;
  if (sumw>0.) val = sumwy/sumw;
    
  return val;
}

void GaussianKernelSmoother::getSmoothHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: GaussianKernelSmoother::getSmoothHisto: input histo does not exist." << std::endl;

  this->h_out = (TH1D*)this->h_in->Clone(this->h_in->GetName()+(TString)"_smoothed");
  this->h_out->SetDirectory(0);

  Int_t nbins = this->h_out->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++){
    double x = this->h_out->GetBinCenter(ib);
    double smooth_val = this->getSmoothedValue(this->h_in , x);
    this->h_out->SetBinContent(ib,smooth_val);
  }
}

double GaussianKernelSmoother::rescaling( double val ){
  if ( kernelDistance == "lin" ) return val;
  if ( kernelDistance == "log" ) return log(val);
  else return val;
}

double GaussianKernelSmoother::invertRescaling( double val ){
  if ( kernelDistance == "lin" ) return val;
  if ( kernelDistance == "log" ) return exp(val);
  else return val;
}

void GaussianKernelSmoother::getContSmoothHisto(){

  if ( ! this->h_in ) std::cerr << "ERROR: GaussianKernelSmoother::getSmoothHisto: input histo does not exist." << std::endl;

  double mini = this->h_in->GetXaxis()->GetBinLowEdge(1);
  double maxi = this->h_in->GetXaxis()->GetBinUpEdge(this->h_in->GetNbinsX());

  if ( fabs(mini)<1e-8 ) mini = this->h_in->GetXaxis()->GetBinUpEdge(1)/10.;

  const int nbins = 500;
  double bins[nbins+1]={0};
  double xs[nbins]={0};
  double ys[nbins]={0};

  double dx = (this->rescaling(maxi) - this->rescaling(mini))/nbins;

  for (int i=0; i<=nbins; i++){
    bins[i] = this->invertRescaling(  this->rescaling(mini) + i*dx  );
  }
  
  for (int i=0; i<nbins; i++){
    xs[i] = ( bins[i]+bins[i+1] )/2;
    ys[i] = this->getSmoothedValue( this->h_in , xs[i] );
  } 

  if ( this->doErrors ){
    const int NTOYS=100;
    double ys_rnd[nbins]={0};
    double ys_err[nbins]={0};
    double ys_low[nbins]={0};
    double ys_high[nbins]={0};

    std::vector<double> toys[nbins];
    for (int it=0; it<NTOYS; it++){
      TH1D *h_rnd=this->fluctuateHisto();
      for (int ib=0; ib<nbins; ib++){
	ys_rnd[ib]=this->getSmoothedValue( h_rnd , xs[ib] );
	toys[ib].push_back( ys_rnd[ib] );
      }
      delete h_rnd;
    }
    for (int i=0; i<nbins; i++){
      ys_err[i]  = this->std_dev( toys[i] );
      ys_low[i]  = ys[i]-ys_err[i];
      ys_high[i] = ys[i]+ys_err[i];
    }

    g_out=new TGraphAsymmErrors( nbins , xs , ys , 0 , 0 , ys_err , ys_err );
    g_out->SetLineColor(38);
    g_out->SetMarkerColor(kRed);
    g_out->SetMarkerStyle(20);
    g_out->SetMarkerSize(0.5);
    g_out_err=new TGraphAsymmErrors( nbins , xs , ys_err );
  } else{
    g_out = new TGraphAsymmErrors( nbins , xs , ys , 0 , 0 , 0 , 0 );
    g_out_err=new TGraphAsymmErrors( nbins , xs , 0 );
  }

}

TH1D* GaussianKernelSmoother::makeWeights( TH1D* h ){

  TH1D *h_w = (TH1D*)h->Clone( (TString)h->GetName()+"_weights" );
  h_w->SetDirectory(0);

  double sum_weights=0;

  for (int ib=1; ib<=h->GetNbinsX(); ib++){
    double weight = 0.;
    if ( fabs( h->GetBinError(ib) ) > 1e-8 ){
      weight = 1./( pow(h->GetBinError(ib),2) );
    }
    h_w->SetBinContent( ib , weight );
    h_w->SetBinError(   ib , 0      );
    sum_weights+=weight;
  }

  h_w->Scale( 1/sum_weights );
  //  std::cout << "XXX " << h_w->Integral(-1,-1) << std::endl;

  return h_w;
}


double GaussianKernelSmoother::std_dev( std::vector<double> v ){
  unsigned size=v.size();

  double mean=0;
  for (unsigned i=0; i<size; i++) mean+=v.at(i)/size;

  double std_dev=0;
  for (unsigned i=0; i<size; i++){
    double diff = v.at(i) - mean;
    std_dev+=diff*diff;
  }
  std_dev/=size;
  std_dev=sqrt(std_dev);
 
  return std_dev;
}

void GaussianKernelSmoother::createTestHisto(){

  TF1 *p1 = new TF1("p1","[0]+x*[1]+x*x*[2]");
  p1->SetParameters(5., 10., 0.);

  //  std::cout << p1->Eval(-5) << "\t" << p1->Eval(0) << "\t" << p1->Eval(5) << std::endl;

  this->h_in = new TH1D("h_in","",100,1,11);
  this->h_in->FillRandom("p1",1000);
  this->h_in->SetMarkerColor(kBlack);
  this->h_in->SetLineColor(kBlack);

}

