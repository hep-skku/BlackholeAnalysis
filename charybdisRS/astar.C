void astar()
{
  double dummy;
  ifstream flux("data_files/cfluxke_n2.dat");
  TGraph *g = new TGraph();
  for(int k=1; k <= 26*(170+157+121); k++)
  {
    flux >> dummy;
    //g->SetPoint(k-1,k,dummy);
  }
  for(int k=1; k <=101*5; ++k)
  {
    flux >> dummy;
    g->SetPoint(k-1,k,dummy);
  }
  flux.close();
  TCanvas *c = new TCanvas("c","c",500,500);
  g->Draw("A*");
}

