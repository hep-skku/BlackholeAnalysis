void hist()
{
  TH1F* h1 = new TH1F("h1","h1",9,0,9);
  TH1F* h2 = new TH1F("h2","h2",9,0,9);
  TLegend* leg = new TLegend(0.50, 0.50, 0.92, 0.92);
  ifstream bt("bt.txt");
  string dummy;
  int id;
  //while(test.eof()==false)
  //{
  while ( bt >> dummy >> id )
  {
    if(id==5 || id==-5)
      h1->AddBinContent(3);
    if(id==6 || id==-6)
      h1->AddBinContent(4);
    if(id==21)
      h1->AddBinContent(5);
    if(id==22)
      h1->AddBinContent(6);
    if(id==23)
      h1->AddBinContent(7);
    if(id==24 || id==-24)
      h1->AddBinContent(8);
    if(id==25)
      h1->AddBinContent(9);
  } 
  bt.close();
  ifstream all("all.txt");
  while ( all >> dummy >> id )
  {
    if(id==5 || id==-5)
      h2->AddBinContent(3);
    else if(id==6 || id==-6)
      h2->AddBinContent(4);
    else if(id==21)
      h2->AddBinContent(5);
    else if(id==22)
      h2->AddBinContent(6);
    else if(id==23)
      h2->AddBinContent(7);
    else if(id==24 || id==-24)
      h2->AddBinContent(8);
    else if(id==25)
      h2->AddBinContent(9);
    else if(id==12 || id==14 || id==16)
      h2->AddBinContent(2);
    else
      h2->AddBinContent(1);
  }
  all.close();
  h2->GetXaxis()->SetBinLabel(1,"light");
  h2->GetXaxis()->SetBinLabel(2,"#nu");
  h2->GetXaxis()->SetBinLabel(3,"b");
  h2->GetXaxis()->SetBinLabel(4,"t");
  h2->GetXaxis()->SetBinLabel(5,"g");
  h2->GetXaxis()->SetBinLabel(6,"#gamma");
  h2->GetXaxis()->SetBinLabel(7,"Z");
  h2->GetXaxis()->SetBinLabel(8,"W");
  h2->GetXaxis()->SetBinLabel(9,"H");
  

  TCanvas* c = new TCanvas("c","c",500,500);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h2->Draw();
  h1->Draw("same");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"b&t only","l");
  leg->AddEntry(h2,"all","l");
  leg->Draw();
  cout << h1->Integral() << endl;
  cout << h2->Integral() << endl;
  
}
