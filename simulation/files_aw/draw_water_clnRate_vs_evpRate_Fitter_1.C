double mypol1func(double* x, double * par){
  return par[0]+par[1]*x[0];
}

double mypol1func0(double* x, double * par){
  return par[0]*x[0];
}


double simple_line(double* x, double* par){
  return par[0]+par[1]*x[0];
}

double cln1[4]={};
double clne1[4]={};
double evp1[4]={};
double evpe1[4]={};
double cln2[4]={};
double clne2[4]={};
double evp2[4]={};
double evpe2[4]={};

const int N_pars = 4;
double par_slope = 0;
double par_c1 = 0;
double par_c2 = 0;
double par_c = 0;

void fcn(Int_t &npar, Double_t *deriv, Double_t &f, Double_t *par, Int_t iflag)
{
  double chi2_1 = 0, chi2_2=0;

  par_slope = par[0];
  par_c1 = par[1];
  par_c2 = par[2];
  par_c = par[3];

  for(int i=0;i<4;i++){
    chi2_1 += TMath::Power( cln1[i]- (par_slope*(evp1[i]+par_c1)+par_c), 2)/TMath::Power(clne1[i], 2) ;
  }
  for(int i=0;i<4;i++){
    chi2_2 += TMath::Power( cln2[i]- (par_slope*(evp2[i]+par_c2)+par_c), 2)/TMath::Power(clne2[i], 2) ;
  }

  f = chi2_1 + chi2_2
      + TMath::Power( par_c1/evpe1[0] , 2)
      + TMath::Power( par_c2/evpe2[0] , 2)
       ;
}

void draw_water_clnRate_vs_evpRate_Fitter_1(){
  fstream fin("new_water_clnRate_vs_evpRate_combineTwoSets_v2.txt",ios::in);
  double evp, evpe, cln, clne;
  int n=0;
  while(fin>>evp>>evpe>>cln>>clne){
    if(n<=3){
        cln1[n] = cln;
        clne1[n] = clne;
        evp1[n] = evp;
        evpe1[n] = evpe;
    }
    if(n>=4 && n<=7){
        cln2[n-4] = cln;
        clne2[n-4] = clne;
        evp2[n-4] = evp;
        evpe2[n-4] = evpe;
    }
    n++;
  }
  fin.close();

  // minimization - 
  TMinuit *gMinuit = new TMinuit(N_pars);  //initialize TMinuit with a maximum of N_pars params
  gMinuit->SetFCN(fcn);

  //Double_t arglist[N_pars];
  //Int_t ierflg = 0;

  //arglist[0] = 1;
  //gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  Double_t vstart[N_pars] = {0.85,0.,0., 0};
  Double_t step[N_pars] = {0.005, 0.0001, 0.0001, 0.0000001};
  Double_t minVal[N_pars] = {0.75, 0, 0,  -0.001};
  Double_t maxVal[N_pars] = {1.2, 0.1, 0.1,  0.001};
  string parName[N_pars]={"Henry", "n1", "n2", "c"};

  for(int i=0; i<N_pars; i++){
      gMinuit->DefineParameter(i, parName[i].c_str(), vstart[i],
                                  step[i], minVal[i], maxVal[i]);
  }
  gMinuit->FixParameter(3);
  gMinuit->Migrad();

  double amin,edm,errdef;
  double nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //cout<<"begin tag ..."<<endl;
  // gMinuit->mnprin(3,amin);
  //cout<<"end tag ..."<<endl;
  cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;
  //double v=0, ve=0;
  //cout<<gMinuit->GetParameter(2, v, ve);
  //cout<<v<<"\t"<<ve<<endl;
  

  double outpar[N_pars], errors[N_pars];
  cout<<"Output fit parameters: "<<endl;
  cout<<"FCN "<<amin<<endl;
  for(int i=0; i<N_pars; i++){
      gMinuit->GetParameter(i, outpar[i], errors[i]);
      cout<<i<<"\t"<<parName[i]<<"\t"<<outpar[i]<<"\t"<<errors[i]<<endl;
  }

  //for(int i=0; i<3; i++) { evp1[i] += outpar[1]; evpe1[i] = errors[1]; }
  //for(int i=0; i<3; i++) { evp2[i] += outpar[2]; evpe2[i] = errors[2]; }
  //for(int i=0; i<2; i++) { evp3[i] += outpar[3]; evpe3[i] = errors[3]; }

  TGraphErrors* g1 = new TGraphErrors();//(3, cln1, evp1, clne1, evpe1);
  for(int i=0; i<4; i++){g1->SetPoint(i, evp1[i], cln1[i]); g1->SetPointError(i, evpe1[i], clne1[i]);}
  TGraphErrors* g2 = new TGraphErrors();//3, cln2, evp2, clne2, evpe2);
  for(int i=0; i<4; i++){g2->SetPoint(i, evp2[i], cln2[i]); g2->SetPointError(i, evpe2[i], clne2[i]);}
  
  TGraph* gFrame = new TGraph();
  gFrame->SetPoint(0, 0, 0.);
  gFrame->SetPoint(1, 0.025, 0.055);
  TCanvas* c = new TCanvas();
  gFrame->Draw("AP");
  gFrame->GetXaxis()->SetTitle("r_{evp} (mole/s)");
  gFrame->GetYaxis()->SetTitle("r_{cln} (mole/s)");
  g1->Draw("P");
  g1->SetMarkerStyle(24);
  g2->Draw("P");
  g2->SetMarkerStyle(8);
  TLegend* leg0 = new TLegend(0.15,0.7,0.45,0.9);
  leg0->SetHeader("Data sets");
  leg0->AddEntry(g1,"Feburary 2016","lp");
  leg0->AddEntry(g2,"December 2018","lp");
  leg0->Draw("same");
  
  TF1* thefit = new TF1("thefit",simple_line,0,0.022,2);
  thefit->SetParameter(0, outpar[3]);
  thefit->SetParameter(1, outpar[0]);
  
  thefit->Draw("same");

  
  TF1* thefit1 = new TF1("thefit1",simple_line,0,0.022,2);
  thefit1->SetParameter(0, 0);
  thefit1->SetParameter(1, 0.94);
  
  thefit1->Draw("same"); thefit1->SetLineColor(4); thefit1->SetLineStyle(2);
  
}


