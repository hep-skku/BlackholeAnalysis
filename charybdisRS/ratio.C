void ratio()
{
  const double pi = TMath::Pi();
  const double kl = 11.3*pi;
  double k = 10;
  const double mBH = 10000;
  const double mPL = 1000;
  const double nu[12] = {-0.579,-0.517,-0.473,-0.579,-0.517,-0.473,0.742,0.558,
  -0.339,0.711,0.666,0.533};
  const double dHiggs = 1;
 
  double dofTop = 2*3; //(ptl+anti-ptl)*3color
  double dofHiggs = 1; // scalar
  double dofGauge = 8*2+2+2*(1+2)+(1+2); //(8gluon*2polar)+(2polar photon)+ 
  //(ptl+anti-ptl W)*(1scalar+2vector)+Z(1scalar+2vector)
  double dofAllQuark = 2*3; //(ptl+anti-ptl)*3color 

  double wFctrQuark = 3./4.;
  double wFctrHiggs = 1;
  double wFctrGauge = 1;
 
  double dFermionL[6];
  double dFermionR[6];

  const int n = 100;
  double step = pow(100,1./n);
  TGraph *g = new TGraph();
  TGraph *g1 = new TGraph();
  TGraph *g2 = new TGraph();
  TGraph *g3 = new TGraph();
  TGraph *g4 = new TGraph();

  for(int i=1; i<=n+1; ++i )
  { 
    double epsln = sqrt(2.*mBH/3./pi/mPL)/kl*k/mPL;
    double dGauge = epsln;
    double dFermionSum = 0;    
    for(int j=0; j<=5; ++j)
    {
      dFermionL[j] = (1.-exp(-1.*(1.+2*nu[j])*kl*epsln))/(1.-exp(-1.*kl*(1.+2.*nu[j])));
      dFermionR[j] = (1.-exp(-1.*(1.-2*nu[j+6])*kl*epsln))/(1.-exp(-1.*kl*(1.-2.*nu[j+6])));
      dFermionSum = dFermionSum + dFermionL[j] + dFermionR[j];
    }       
    double decayAll = dofAllQuark*wFctrQuark*dFermionSum+dofGauge*wFctrGauge*dGauge+
                      dofHiggs*wFctrHiggs*dHiggs;
    double decayTop = dofTop*wFctrQuark*(dFermionL[2]+dFermionR[2]);
    double decayHiggs = dofHiggs*wFctrHiggs*dHiggs;

    g->SetPoint(i-1,k/mPL,dFermionL[2]+dFermionR[2]);
    g1->SetPoint(i-1,k/mPL,decayTop/decayHiggs);
    
    g2->SetPoint(i-1,k/mPL,(decayTop+decayHiggs)/decayAll);
    g3->SetPoint(i-1,k/mPL,decayTop/decayAll);
    g4->SetPoint(i-1,k/mPL,decayHiggs/decayAll);
    k = k*step;
  }
  TCanvas *c = new TCanvas("c","c",500);
  c->SetLogx();
  g->GetXaxis()->SetTitle("k/M_{pl}");
  g->GetYaxis()->SetTitle("D_{top}/D_{higgs}");
  g->SetMarkerSize(0.5);
  g->Draw("al");

  TCanvas *cc = new TCanvas("cc","cc",500);
  cc->SetLogx();
  g1->GetXaxis()->SetTitle("k/M_{pl}");
  g1->GetYaxis()->SetTitle("Br(BH->BH+Top)/Br(BH->BH+Higgs)");
  g1->SetMarkerSize(0.5);
  g1->Draw("al");

  TCanvas *ccc = new TCanvas("ccc","ccc",500);
  ccc->SetLogx();
  g3->SetLineColor(kRed);
  g4->SetLineColor(kBlue);
  g2->SetMarkerSize(0.5);
  g3->SetMarkerSize(0.5);
  g4->SetMarkerSize(0.5);
  g3->GetXaxis()->SetTitle("k/M_{pl}");
  g3->GetYaxis()->SetTitle("Branching Ratio");
  g2->GetXaxis()->SetTitle("k/M_{pl}");
  g2->GetYaxis()->SetTitle("Branching Ratio");
  g4->GetXaxis()->SetTitle("k/M_{pl}");
  g4->GetYaxis()->SetTitle("Branching Ratio");
  g2->SetMinimum(0);

  g2->Draw("AL");
  g3->Draw("L");
  g4->Draw("L");
  TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->AddEntry(g2,"Top+Higgs","l");
  l->AddEntry(g3,"Top","l");
  l->AddEntry(g4,"Higgs","l");
  l->Draw();

}
