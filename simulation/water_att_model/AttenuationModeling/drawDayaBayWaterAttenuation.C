{
  fstream fin("water_data_DayaBay_eV_vs_cm.txt",ios::in);
  double eV, cm, wavelength;
  // h_bar*c = 197.32698 MeV*fm; MeV*fm = eV*nm
  double constant = 2.0*TMath::Pi()*197.32698;

  TGraph* g = new TGraph();
  int n=0;
  while(fin>>eV>>cm){
    wavelength = constant/eV;
    g->SetPoint(n,wavelength, cm*10.0);
    n++;
  }
  fin.close();

/*

  TFile* rtfile = new TFile("myWaterData.root","read");
  TCanvas* c = (TCanvas*)rtfile->Get("c1_n2");
  TGraph* gMyPureWater = (TGraph*)c->GetListOfPrimitives()->FindObject("gPureWater");
  TGraph* gMyOnetonWater = (TGraph*)c->GetListOfPrimitives()->FindObject("gOnetonWater");
*/

  TFile* rtfile = new TFile("RatPac_WaterAttenuation.root","read");
  TCanvas* c = (TCanvas*)rtfile->Get("c1");
  TGraph* gMyPureWater = (TGraph*)c->GetListOfPrimitives()->FindObject("gWaterRatPac");

  TCanvas* c0 = new TCanvas();

  g->Draw("AP");
  g->SetMarkerStyle(8);
  g->SetMarkerColor(1);
  g->GetXaxis()->SetTitle("Wavelength (nm)");
  g->GetYaxis()->SetTitle("Attenuation length (mm)");

  gMyPureWater->Draw("P");
  gMyPureWater->SetMarkerStyle(8);
  gMyPureWater->SetMarkerColor(2);

/*
  gMyOnetonWater->Draw("P");
  gMyOnetonWater->SetMarkerStyle(8);
  gMyOnetonWater->SetMarkerColor(3);
*/
  TLegend* leg = new TLegend(0.7,0.6,0.9,0.9);
  leg->AddEntry(gMyPureWater,"Data used in RatPac (from my measurements","p");
  //leg->AddEntry(gMyOnetonWater,"Measured one ton water","p");
  leg->AddEntry(g,"Data in DayaBay simulation","p");
  leg->Draw();
}
