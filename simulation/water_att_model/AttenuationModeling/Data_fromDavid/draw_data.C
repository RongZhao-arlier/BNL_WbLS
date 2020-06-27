#include <iostream>     // std::cout
#include <algorithm>    // std::remove

void draw_data(){

  fstream fin2("../water_data_DayaBay_eV_vs_cm.txt",ios::in);
  double eV, cm, wavelength;
  // h_bar*c = 197.32698 MeV*fm; MeV*fm = eV*nm
  double constant = 2.0*TMath::Pi()*197.32698;

  double a_at_200 = 0.001;
  double defa_at_200 = 0.007;

  TGraph* gDayaBay = new TGraph();
  TGraph* gDayaBayModify = new TGraph();

  int n=0;
  while(fin2>>eV>>cm){
    wavelength = constant/eV;
    gDayaBay->SetPoint(n,wavelength, 1./cm );

    double adef = defa_at_200/TMath::Power(wavelength/200.0,4);
    double b = a_at_200/TMath::Power(wavelength/200.0,4);
    double anew = 1./cm-adef+b;
    double ltot = 1./anew;

    gDayaBayModify->SetPoint(n,wavelength,1./ltot);

    n++;
  }
  fin2.close();

  const int N=13;
  
  string aline;

  TGraph* g[N];
  char temp[100];
  for(int i=0;i<N;i++){
    sprintf(temp,"abs%d.vec",i+1);
    g[i] = new TGraph();
    g[i]->SetName(temp);

    fstream fin(temp,ios::in);
    cout<<temp<<" is read."<<endl;
    int cnt=0;
    double wavelength=0; // in nm
    double absorption=0; // in 1/cm
    string delimiter = "\t";
    while(std::getline(fin,aline)){
      if(aline.find("*")==0){
        //cout<<aline<<endl;
        continue;
      }
      stringstream stream(aline);
      stream>>wavelength>>absorption;
      //cout<<wavelength<<"\t"<<absorption<<endl;
      g[i]->SetPoint(cnt, wavelength, absorption);
      cnt++;
    }
    fin.close();
  }

  TCanvas* c = new TCanvas();
  TGraph* g0 = new TGraph();
  g0->SetPoint(0,200,0);
  g0->SetPoint(1,800,0.1);
  g0->Draw("AP");
  g0->GetXaxis()->SetRangeUser(200,800);
  g0->GetYaxis()->SetRangeUser(0.00001,0.1);
  TLegend* leg = new TLegend(0.1,0.6,0.6,0.9);
  for(int i=0; i<N;i++){
    if(i==0){
      g[i]->Draw("P");
    }
    else{g[i]->Draw("P");}
    g[i]->SetMarkerColor(i%3+1);
    g[i]->SetMarkerStyle(21+i/3);
    leg->AddEntry(g[i],g[i]->GetName(),"lp");
  }
  leg->Draw();
  /*
  TCanvas* ctest[N];
  for(int i=0; i<N; i++){
    ctest[i] = new TCanvas();
    g[i]->Draw("AP");
  }
  */
//  c->cd();
  gDayaBay->Draw("P");
  gDayaBay->SetMarkerStyle(24);
  gDayaBay->SetMarkerColor(6);

  
 leg->AddEntry(gDayaBay,"dayabay","lp");

gDayaBayModify->Draw("P");
  gDayaBayModify->SetMarkerStyle(25);
  gDayaBayModify->SetMarkerColor(3);

 leg->AddEntry(gDayaBayModify,"dayabay_modify","lp");

}


