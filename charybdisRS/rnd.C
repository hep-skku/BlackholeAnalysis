#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TF1.h>
#include <TLegend.h>
#include <cmath>
 
using namespace std;
 
const char* prefix = "20130417";
 
//const double conv = 3.894e5; // nb*GeV2
const double conv = 3.894e-1; // nb*TeV2
const int dim = 11;
const double mD = 1.;
const double sqrts = 14;
const double pi = TMath::Pi();
const double fiop = 4.*pow(1.+(dim-2.)*(dim-2.)/4., -2./(dim-3.));
const double kn = pow(pow(2.,dim-4.)*pow(pi,(dim-7.)/2.)*TMath::Gamma((dim-1.)/2)/(dim-2.),1./(dim-3));
 
const int n = 1000000;
const double mMin = 5., mMax = sqrts;
const int nBinsX = 100;
const double binW = (mMax-mMin)/nBinsX;
 
double ffAN(double* xx, double*)
{
  const double x = xx[0];
  const double rs2 = kn*kn*pow(x/mD, 2./(dim-3))/mD/mD;
  const double y = -2*conv*(fiop*pi*rs2)*x/sqrts/sqrts*2*log(x/sqrts);
  if ( y < 0 ) return 0;
 
  return y;
}
 
void rnd()
{
  gROOT->ProcessLine(".x ~/work/rootlogon.C");
  TStopwatch time;
 
  TLegend* leg = new TLegend(0.70, 0.70, 0.92, 0.92);
 
  TH1F* hFrame = new TH1F("hFrame", "hFrame;Mass (TeV/c^{2});Parton level d#sigma/dm (nb/TeVc^{2})", nBinsX, mMin, mMax);
 
  TH1F* hCh = new TH1F("hCh", "hCh;Mass (TeV/c^{2})", nBinsX, mMin, mMax);
  time.Reset();
  time.Start();
  for ( int i=0; i<n; ++i )
  {
    const double bPow = 2./(dim-3) - 7;
    const double m = pow(gRandom->Uniform(pow(mMax, bPow), pow(mMin, bPow)), 1./bPow);
    if ( m < mMin || m > mMax || m > sqrts )
    {
      --i;
      continue;
    }
 
    const double u = m*m/sqrts/sqrts;
    const double fact = 2.*kn*kn/bPow*pi/pow(mD, 2./(dim-3)+2);
    //const double w = fact*pow(m, 7)*log(u)*(pow(mMin, bPow) - pow(mMax, bPow));
    const double w = fiop*fact*pow(m, 7)*u*log(u)*(pow(mMin, bPow) - pow(mMax, bPow));
 
    hCh->Fill(m, w);
    //hCh->Fill(m);
  }
  time.Stop();
  time.Print();
 
  TH1F* hBM = new TH1F("hBM", "hBM;Mass (TeV/c^{2});Normalized", nBinsX, mMin, mMax);
  time.Reset();
  time.Start();
  for ( int i=0; i<n; ++i )
  {
    const double m = 1./gRandom->Uniform(1./mMax, 1./mMin);
    if ( m < mMin || m > mMax || m > sqrts )
    {
      --i;
      continue;
    }
 
    const double u = m*m/sqrts/sqrts;
    //const double v = pow(u, gRandom->Uniform());
    const double rs2 = kn*kn*pow(m/mD, 2./(dim-3))/mD/mD;
    const double w = 2*(fiop*pi*rs2)*(1./mMax-1./mMin)*m*u*log(u);
 
    hBM->Fill(m, w);
    //hBM->Fill(m);
  }
  time.Stop();
  time.Print();
 
  const double uu = mMin*mMin/sqrts/sqrts;
  TH1F* hXY = new TH1F("hXY", "hXY;Mass (TeV/c^{2});Normalized", nBinsX, mMin, mMax);
  time.Reset();
  time.Start();
  for ( int i=0; i<n; ++i )
  {
    const double x1 = gRandom->Uniform(uu, 1);
    const double x2 = gRandom->Uniform(uu, 1);
    const double m = sqrt(x1*x2)*sqrts;
 
    if ( m < mMin || m > mMax )
    {
      --i;
      continue;
    }
 
    const double rs2 = kn*kn*pow(m, 2./(dim-3));
    const double w = fiop*pi*rs2*(1-uu*(1-log(uu)));
 
    hXY->Fill(m, w);
  }
  time.Stop();
  time.Print();
 
  TH1F* hAN = new TH1F("hAN", "hAN;Mass (TeV/c^{2})", nBinsX, mMin, mMax);
  TF1* fAN = new TF1("fAN", ffAN, mMin, mMax);
  time.Reset();
  time.Start();
  for ( int i=0; i<n; ++i )
  {
    const double m = fAN->GetRandom(mMin, mMax);
    const double f = -kn*kn*fiop*pi/sqrts/sqrts*(dim-3)/(dim-2);
    const double w1 = pow(mMax/mD, 2*(1./(dim-3)+1))*(2*log(mMax/sqrts)-(dim-3)/(dim-2.));
    const double w2 = pow(mMin/mD, 2*(1./(dim-3)+1))*(2*log(mMin/sqrts)-(dim-3)/(dim-2.));
    hAN->Fill(m, f*(w1-w2));
  }
  time.Stop();
  time.Print();
 
  TCanvas* c = new TCanvas("c", "c", 500, 500);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  hAN->Scale(conv/binW/n);
  hXY->Scale(conv/binW/n);
  hBM->Scale(conv/binW/n);
  hCh->Scale(conv/binW/n);
 
  cout << "fAN : " << fAN->Integral(mMin, mMax) << endl;
  cout << "hAN : " << hAN->Integral()*binW << endl;
  cout << "hXY : " << hXY->Integral()*binW << endl;
  cout << "hBM : " << hBM->Integral()*binW << endl;
  cout << "hCh : " << hCh->Integral()*binW << endl;
 
  const double yMax = fAN->GetMaximum();
 
  hFrame->SetMaximum(1.2*yMax);
  hFrame->SetMinimum(0);
 
  fAN->SetLineWidth(10);
  fAN->SetLineColor(kGreen);
  hAN->SetLineColor(kBlack);
  hXY->SetLineColor(kBlue);
  hBM->SetLineColor(kMagenta);
  hCh->SetLineColor(kRed);
 
  hFrame->Draw();
  fAN->Draw("same");
  hBM->Draw("same");
  hCh->Draw("same");
  hAN->Draw("same");
  hXY->Draw("same");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(fAN, "Formula", "l");
  leg->AddEntry(hAN, "Analytic", "l");
  leg->AddEntry(hBM, "BlackMax", "l");
  leg->AddEntry(hCh, "Charybdis", "l");
  leg->AddEntry(hXY, "ToyMC", "l");
  leg->Draw();
 
  c->Print(Form("%s_Mass_D%d_%.0fTeV.png", prefix, dim, sqrts));
}
