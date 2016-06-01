void test()
{
  ifstream fin;
  fin.open("test.txt");
  //string dummy1, dummy2;
  int id,q1,q2,q3,q4,q5;
  double pX,pY,topPt,btmPt,energy,mass;
  double q6,q7;
  gStyle->SetPalette(1);

  double topSt = 0;
  double btmSt = 0;
  int nTop = 0;
  int nBtm = 0;
  //int nBtm = 0;

  //TH1F* hNTop = new TH1F("hNTop", "hNTop", 7, 0, 7);
  //TH1F* hNBtm = new TH1F("hNBtm", "hNBtm", 7, 0, 7);
  TH2F* hNHvy = new TH2F("hNHvy", "hNHvy;Top;Bottom", 6, 0, 6, 6, 0, 6);
  TH1F* hHvyPt = new TH1F("hHvyPt","hHvyPt;sT(GeV/c);Number of Events/20GeV/c",100,0,2000);
  TH1F* hTopPt = new TH1F("hTopPt","hTopPt;sT(GeV/c);Number of Events/20GeV/c",100,0,2000);
  TH1F* hBtmPt = new TH1F("hBtmPt","hBtmPt;sT(GeV/c);Number of Events/20GeV/c",100,0,2000);

  int nEvent = 0;
  string line;
  while( getline(fin, line) && !fin.eof() )
  {
    const bool isEventEnd = line[7] == '<';
    if ( isEventEnd )
    {
      ++nEvent;
      // Fill histograms
      //hNTop->Fill(nTop);
      hNHvy->Fill(nTop,nBtm);
      // Reset multiplicity
      nTop = 0;
      nBtm = 0;

      //nBtm = 0;

      cout << nEvent << endl;
    }
    else
    {
      stringstream ss;
      ss = line;
      string dummy;
      ss >> dummy;
      ss >> id;

      //if ( abs(id) == 6 ) ++nTop;
      if ( abs(id) == 6 ) 
        {  
          ++nTop;
          ss >> q1 >> q2 >> q3 >> q4 >> q5;
          ss >> pX >> pY;
          topPt = sqrt(pX*pX+pY*pY);
          hTopPt->Fill(topPt);
          hHvyPt->Fill(topPt);
        }

      if ( abs(id) == 5 ) 
        {  
          ++nBtm;
          ss >> q1 >> q2 >> q3 >> q4 >> q5;
          ss >> pX >> pY;
          btmPt = sqrt(pX*pX+pY*pY);
          hBtmPt->Fill(btmPt);
          hHvyPt->Fill(btmPt);
        }
        
    }
  }
  fin.close();
  TList* list = new TList;
  //list->Add(hNTop);
  //list->Add(hNBtm);
  //hNHvy.Merge(list);

  hNHvy->SetLineColor(kBlack);
  //hNTop->SetLineColor(kBlue);
  //hNBtm->SetLineColor(kRed);
  hNHvy->Draw("colztext");
  TCanvas* c = new TCanvas("c","c", 500, 500);
  hTopPt->Draw();
  TCanvas* cc = new TCanvas("cc","cc", 500, 500);
  hBtmPt->Draw();
  TCanvas* ccc = new TCanvas("ccc","ccc", 500, 500);
  hHvyPt->Draw(); 
  //hNTop->Draw("sames");
  //hNBtm->Draw("sames");
}
