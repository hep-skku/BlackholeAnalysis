void data()
{
  ifstream flux("data_files/flux_n2.dat");
  TGraph* g1 = new TGraph();
  TGraph* g2 = new TGraph();
  TGraph* g3 = new TGraph();
  
  double f1,f2,f3;
  double cfluxke,dummy;
  int i=0;
  int MK0 = 3;
  int MK1 = 3;
  int MK2 = 3;
  int ME = 4;
  while(!flux.eof())
  {
    flux >> f1 >> f2 >> f3;
    g1->SetPoint(i,i+1,f1);
    g2->SetPoint(i,i+1,f2);
    g3->SetPoint(i,i+1,f3);
    ++i;
  }
  flux.close();

  ifstream fluxke("data_files/cfluxke_n2.dat");
  TGraph* CFLUX2S0K = new TGraph(); 
  TGraph* CFLUX2S1K = new TGraph(); 
  TGraph* CFLUX2S2K = new TGraph(); 
  for(int ma=1; ma <= 26; ++ma)
  { 
    for(int k=1; k <=170;++k)
    {
      
      fluxke >> dummy;
      if(k==MK0) cfluxke = dummy;
    }
    CFLUX2S0K->SetPoint(ma-1,ma,cfluxke);
  }
  for(int ma=1; ma <= 26; ++ma)
  { 
    for(int k=1; k <=157;++k)
    {
      fluxke >> dummy;
      if(k==MK1) cfluxke = dummy;
    }
    CFLUX2S1K->SetPoint(ma-1,ma,cfluxke);
  }
  for(int ma=1; ma <= 26; ++ma)
  { 
    for(int k=1; k <=121;++k)
    {
      fluxke >> dummy;
      if(k==MK2) cfluxke = dummy;
    }
    CFLUX2S2K->SetPoint(ma-1,ma,cfluxke);
  }
  TGraph *CFLUX2S0KE = new TGraph();
  TGraph *CFLUX2S1KE = new TGraph();
  TGraph *CFLUX2S2KE = new TGraph();
  
  for(int ma=1; ma <= 26; ++ma)
  {
    for(int k=1; k <= 169; ++k)
    {
      for(int e=1; e <= 101; ++e)
      {
        fluxke >> dummy;
        if(k==MK0 && e==ME) cfluxke = dummy;
      }
    }
    CFLUX2S0KE->SetPoint(ma-1,ma,cfluxke);
  }
  for(int ma=1; ma <= 26; ++ma)
  {
    for(int k=1; k <= 156; ++k)
    {
      for(int e=1; e <= 101; ++e)
      {
        fluxke >> dummy;
        if(k==MK1 && e==ME) cfluxke = dummy;
      }
    }
    CFLUX2S1KE->SetPoint(ma-1,ma,cfluxke);
  }
  for(int ma=1; ma <= 26; ++ma)
  {
    for(int k=1; k <= 120; ++k)
    {
      for(int e=1; e <= 101; ++e)
      {
        fluxke >> dummy;
        if(k==MK2 && e==ME) cfluxke = dummy;
      }
    }
    CFLUX2S2KE->SetPoint(ma-1,ma,cfluxke);
  }
  fluxke.close();
  
  g1->SetLineColor(kRed);
  g2->SetLineColor(kBlue);
  g3->SetLineColor(kGreen);
  
  TCanvas *c = new TCanvas("c","c",500,500);
  //c->SetLogy();
  g3->Draw("al");
  g2->Draw("same");
  g1->Draw("same");
  TLegend *l = new TLegend(0.6, 0.6, 0.9, 0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->AddEntry(g1,"spin 0","l"); 
  l->AddEntry(g2,"spin 1/2","l"); 
  l->AddEntry(g3,"spin 1","l");
  l->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",500,500);
  CFLUX2S0K->SetLineColor(kRed);
  CFLUX2S1K->SetLineColor(kBlue);
  CFLUX2S2K->SetLineColor(kGreen);
  CFLUX2S0K->Draw("ALP"); 
  TCanvas *c3 = new TCanvas("c3","c3",500,500);
  CFLUX2S1K->Draw("ALP"); 
  TCanvas *c4 = new TCanvas("c4","c4",500,500);
  CFLUX2S2K->Draw("ALP"); 
  TCanvas *c5 = new TCanvas("c5","c5",500,500);
  CFLUX2S0KE->Draw("ALP");
  TCanvas *c6 = new TCanvas("c6","c6",500,500);
  CFLUX2S1KE->Draw("ALP");
  TCanvas *c7 = new TCanvas("c7","c7",500,500);
  CFLUX2S2KE->Draw("ALP");
  
  
}
