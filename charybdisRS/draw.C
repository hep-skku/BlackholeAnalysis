void draw()
{
  TFile* f1 = TFile::Open("add.root");
  TH1D* hStFrame = new TH1D("hStFrame", "St;S_{T} (GeV/c);Events per 10GeV/c", 50, 0, 500);
  //TH1D* hMet = new TH1D("hMet", "Missing Et;Missing transverse momentum (GeV/c);Events per 10GeV/c", 50, 0, 500);
  

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->SetLogy();
  hStFrame->SetMinimum(0.1);
  hStFrame->SetMaximum(1.2*(hMET->GetMaximum()));
  hStFrame->Draw();
  hMET->SetLineColor(kBlue);
  hMET->Draw("same");
  TLegend* legN = new TLegend(0.65, 0.65, 0.92, 0.92);
  legN->SetFillStyle(0);
  legN->SetBorderSize(0);
  legN->AddEntry(hSt,"add","l");
  
  TFile* f2 = TFile::Open("rs.root");
  hMET->SetLineColor(kRed); 
  hMET->Draw("same");
  legN->AddEntry(hSt,"rs","l");
  legN->Draw();
  //hST->SetMinimum(0.1);
  //TCanvas* c1 = new TCanvas("c", "c", 500, 500);
  //hSt->Draw();
  
}
