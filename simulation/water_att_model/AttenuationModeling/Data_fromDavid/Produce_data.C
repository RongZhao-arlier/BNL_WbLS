#include <iostream>     // std::cout
#include <algorithm>    // std::remove
#include <string.h>
#include <sstream>

void Produce_data(){

  string otherdata1="{\nname: \"OPTICS\",\nindex: \"water\",\nvalid_begin : [0, 0],\nvalid_end : [0, 0],\n\/\/LIGHT_YIELD: 0.0\nNEUTRON_CAPTURE_TIME_value1: [0.0, 1.0, ],\nNEUTRON_CAPTURE_TIME_value2: [163000.0, 163000.0, ],\nNEUTRON_SLOW_DIFFUSION_CONST_value1: [0.0, 1.0, ],\nNEUTRON_SLOW_DIFFUSION_CONST_value2: [0.03, 0.03, ],\nNEUTRON_FAST_DIFFUSION_RMS_value1: [0.0, 1.0, ],\nNEUTRON_FAST_DIFFUSION_RMS_value2: [50.0, 50.0, ],\nRINDEX_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\nRINDEX_value1: [60.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, ],\nRINDEX_value2: [1.42516, 1.396, 1.37761, 1.35942, 1.34978, 1.34378, 1.3397, 1.33676, 1.33458, 1.33293, 1.33165, 1.33065, 1.32986, 1.3292, ],\nABSLENGTH_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\n" ;

string otherdata2 = "\nRSLENGTH_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\nRSLENGTH_value1: [60.0, 200.0, 800.0, ],\n\/\/RSLENGTH_value2: [10.0e3, 10.0e3, 10.0e3, ],\nRSLENGTH_value2: [100.0e3, 100.0e3, 100.0e3, ],\nPROPERTY_LIST: [\"NEUTRON_CAPTURE_TIME\", \"NEUTRON_SLOW_DIFFUSION_CONST\", \"NEUTRON_FAST_DIFFUSION_RMS\", \"RINDEX\", \"ABSLENGTH\", \"RSLENGTH\", ]\n}" ;

  //std::string otherdata1 = "abc";
  //std::string otherdata2 = "def";

  fstream fin2("water_data_DayaBay_eV_vs_cm.txt",ios::in);
  double eV, cm, wavelength;
  // h_bar*c = 197.32698 MeV*fm; MeV*fm = eV*nm
  double constant = 2.0*TMath::Pi()*197.32698;
  //double a_at_200 = 0.001;
  //double defa_at_200 = 0.007;
  const int N = 20; // number of curves, each curve use one parameter

  fstream fout_attL_data[N];
  double myparameter[N];
  TGraph* myAbsorptionCurve[N];
  TGraph* myAttenuationCurve[N];
  //string myCurveName[N];
  string wavelengthStr = "ABSLENGTH_value1: [";
  string waveNumStr="";
  string dataNumStr[N];//="";
  string dataStrHead = "ABSLENGTH_value2: [";
  string dataStr[N];
  char temp[100];
  for(int i=0; i<N; i++){
    myparameter[i] = -0.007 + 0.0005*i;
    myAbsorptionCurve[i] = new TGraph();
    sprintf(temp,"par=%.3f",myparameter[i]);
    myAbsorptionCurve[i]->SetName(temp);

    myAttenuationCurve[i] = new TGraph();
    sprintf(temp,"p=%.3f",myparameter[i]);
    myAttenuationCurve[i]->SetName(temp);

    sprintf(temp,"OPTICS_water_par_%.4f.ratdb",myparameter[i]);
    fout_attL_data[i].open(temp,ios::out);
    dataStr[i] = dataStrHead;
	dataNumStr[i]="";
  }

  TGraph* gDayaBay = new TGraph();
  TGraph* gDayaBay_attL = new TGraph();

  int n=0;
  double myshift;
  double anew;
  double ltot;
    
  while(fin2>>eV>>cm){
    wavelength = constant/eV;
    gDayaBay->SetPoint(n,wavelength, 1./cm );
    gDayaBay_attL->SetPoint(n,wavelength, cm*0.01 );

    ostringstream s;
    s << wavelength;
	string mystring = s.str();
	if(mystring.find(".")!=std::string::npos)
      waveNumStr.insert(0,s.str()+", ");
    else waveNumStr.insert(0,s.str()+".0, ");
    //wavelengthStr = wavelengthStr + s.str() + ", "; // std::to_string(wavelength);//s.c_str();
    //wavelengthStr = wavelengthStr + ", ";

    for(int i=0; i<N; i++){
      myshift = myparameter[i]/TMath::Power(wavelength/200.0,4);
      anew = 1./cm + myshift;
      ltot = 1./anew;
      double attenuation = ltot*10.0;// convert cm to mm, rat-pac needs mm for attenuation length
      myAbsorptionCurve[i]->SetPoint(n,wavelength,1./ltot);
      myAttenuationCurve[i]->SetPoint(n,wavelength,ltot*0.01); // use m for attenuation length to plot
      ostringstream ss;
      ss << attenuation;
	  mystring = ss.str();
	  if(mystring.find(".")!=std::string::npos)
        dataNumStr[i].insert(0,ss.str()+", ");
	  else dataNumStr[i].insert(0,ss.str()+".0, ");
      //dataStr[i] = dataStr[i] + ss.str() + ", "; //std::to_string(ltot*0.01);//s.c_str();
      //dataStr[i] = dataStr[i] + ", ";
    }
    n++;
  }
  fin2.close();

  for(int i=0; i<N; i++){
    fout_attL_data[i] << otherdata1 << endl;
    fout_attL_data[i] << wavelengthStr;
	fout_attL_data[i] << waveNumStr ;
    fout_attL_data[i] << "], " << endl;
    fout_attL_data[i] << dataStr[i];
	fout_attL_data[i] << dataNumStr[i] ;
    fout_attL_data[i] << "], " << endl;
    fout_attL_data[i] << otherdata2 << endl;
    fout_attL_data[i].close();
  }


  // absorption curves
  TCanvas* c = new TCanvas();
  TGraph* g0 = new TGraph();
  g0->SetPoint(0,200,0);
  g0->SetPoint(1,800,0.1);
  g0->Draw("AP");
  g0->GetXaxis()->SetRangeUser(200,800);
  g0->GetYaxis()->SetRangeUser(0.00001,0.1);
  g0->GetXaxis()->SetTitle("Wavelength (nm)");
  g0->GetYaxis()->SetTitle("Absorption length (1/cm)");
  TLegend* leg = new TLegend(0.1,0.6,0.6,0.9);

  gDayaBay->Draw("P");
  gDayaBay->SetMarkerStyle(24);
  gDayaBay->SetMarkerColor(1);
  leg->AddEntry(gDayaBay,"DayaBay_ref (par=0)","lp");

  for(int i=0; i<N; i=i+2){
    myAbsorptionCurve[i]->Draw("P");
    myAbsorptionCurve[i]->SetMarkerStyle(21+i/3);
    myAbsorptionCurve[i]->SetMarkerColor(1+i%3);
    leg->AddEntry(myAbsorptionCurve[i],myAbsorptionCurve[i]->GetName(),"lp");
  }
  leg->SetNColumns(2);
  leg->Draw();

  // attenuation length curves
  TCanvas* c2 = new TCanvas();
  TGraph* g2 = new TGraph();
  g2->SetPoint(0,200,0);
  g2->SetPoint(1,800,500);
  g2->Draw("AP");
  g2->GetXaxis()->SetRangeUser(200,800);
  //g2->GetYaxis()->SetRangeUser(0.00001,0.1);
  g2->GetXaxis()->SetTitle("Wavelength (nm)");
  g2->GetYaxis()->SetTitle("Attenuation length (m)");
  TLegend* leg2 = new TLegend(0.1,0.6,0.6,0.9);

  gDayaBay_attL->Draw("P");
  gDayaBay_attL->SetMarkerStyle(24);
  gDayaBay_attL->SetMarkerColor(1);
  leg2->AddEntry(gDayaBay_attL,"DayaBay_ref (p=0)","lp");

  for(int i=0; i<N; i=i+2){
    myAttenuationCurve[i]->Draw("P");
    myAttenuationCurve[i]->SetMarkerStyle(21+i/3);
    myAttenuationCurve[i]->SetMarkerColor(1+i%3);
    leg2->AddEntry(myAttenuationCurve[i],myAttenuationCurve[i]->GetName(),"lp");
  }
  leg2->SetNColumns(2);
  leg2->Draw();

  exit(0);
}


