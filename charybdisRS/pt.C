void pt()
{
  ifstream fin;
  fin.open("test.txt");
  //string dummy1, dummy2;
  int id,q1,q2,q3,q4,q5;
  double pX,pY,pT,pTSum,PtSum,energy,mass;
  double q6,q7;

  //int nBtm = 0;

  //TH1F* hNTop = new TH1F("hNTop", "hNTop", 7, 0, 7);
  //TH1F* hNBtm = new TH1F("hNBtm", "hNBtm", 7, 0, 7);
  TH1F* hPtSum = new TH1F("hPtSum", "hPtSum;sT(GeV/c);Number of Events/100GeV/c", 100, 0,10000);

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
      hPtSum->Fill(pTSum);
      // Reset multiplicity
      pTSum = 0;
      //nBtm = 0;

      cout << nEvent << endl;
    }
    else
    {
      stringstream ss;
      ss = line;
      string dummy;
      ss >> dummy >> id  >> q1 >> q2 >> q3 >> q4 >> q5;
      ss >> pX >> pY;
      pT = sqrt(pX*pX+pY*pY);
      pTSum = pTSum + pT;

      //if ( abs(id) == 6 ) ++nTop;
    }

  }
  fin.close();
  TList* list = new TList;
  //list->Add(hNTop);
  //list->Add(hNBtm);
  //hNHvy.Merge(list);

  hPtSum->SetLineColor(kBlack);
  //hNTop->SetLineColor(kBlue);
  //hNBtm->SetLineColor(kRed);
  hPtSum->Draw();
  //hNTop->Draw("sames");
  //hNBtm->Draw("sames");
}
