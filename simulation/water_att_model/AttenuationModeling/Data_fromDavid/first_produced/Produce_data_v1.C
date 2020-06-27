#include <iostream>     // std::cout
#include <algorithm>    // std::remove
#include <string.h>
#include <sstream>

void Produce_data_v1(){

  string otherdata1="{\nname: \"OPTICS\",\nindex: \"water\",\nvalid_begin : [0, 0],\nvalid_end : [0, 0],\n\/\/LIGHT_YIELD: 0.0\nNEUTRON_CAPTURE_TIME_value1: [0.0, 1.0, ],\nNEUTRON_CAPTURE_TIME_value2: [163000.0, 163000.0, ],\nNEUTRON_SLOW_DIFFUSION_CONST_value1: [0.0, 1.0, ],\nNEUTRON_SLOW_DIFFUSION_CONST_value2: [0.03, 0.03, ],\nNEUTRON_FAST_DIFFUSION_RMS_value1: [0.0, 1.0, ],\nNEUTRON_FAST_DIFFUSION_RMS_value2: [50.0, 50.0, ],\nRINDEX_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\nRINDEX_value1: [60.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, ],\nRINDEX_value2: [1.42516, 1.396, 1.37761, 1.35942, 1.34978, 1.34378, 1.3397, 1.33676, 1.33458, 1.33293, 1.33165, 1.33065, 1.32986, 1.3292, ],\nABSLENGTH_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\n" ;

  string otherdata2 = "\nRSLENGTH_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water\nRSLENGTH_value1: [60.0, 200.0, 800.0, ],\n\/\/RSLENGTH_value2: [10.0e3, 10.0e3, 10.0e3, ],\nRSLENGTH_value2: [100.0e3, 100.0e3, 100.0e3, ],\nPROPERTY_LIST: [\"NEUTRON_CAPTURE_TIME\", \"NEUTRON_SLOW_DIFFUSION_CONST\", \"NEUTRON_FAST_DIFFUSION_RMS\", \"RINDEX\", \"ABSLENGTH\", \"RSLENGTH\", ]\n}" ;

  //std::string otherdata1 = "abc";
  //std::string otherdata2 = "def";

  fstream fin2("abs_combin.vec",ios::in);
  string tempstr;
  for(int i=0; i<7; i++) getline(fin2, tempstr);
  double eV, cm, wavelength;
  // h_bar*c = 197.32698 MeV*fm; MeV*fm = eV*nm
  double constant = 2.0*TMath::Pi()*197.32698;
  //double a_at_200 = 0.001;
  //double defa_at_200 = 0.007;
  const int N = 300; // number of curves, each curve use one parameter

  fstream fout_attL_data[N];
  double myparameter[N];
  double myMaxAttL[N];
  double myWaveAtMaxAttL[N];
  TGraph* myAbsorptionCurve[N];
  TGraph* myAttenuationCurve[N];
  TGraph* gMaxAttL_vs_parameter = new TGraph();
  TGraph* gMaxAttL_vs_waveLength = new TGraph();
  //string myCurveName[N];
  string wavelengthStr = "ABSLENGTH_value1: [";
  string waveNumStr="";
  string dataNumStr[N];//="";
  string dataStrHead = "ABSLENGTH_value2: [";
  string dataStr[N];
  char temp[100];
  for(int i=0; i<N; i++){
	myMaxAttL[i]=-1.0;
	myWaveAtMaxAttL[i]=0;
    myparameter[i] = 0.0001 + 0.0002*i;
    myAbsorptionCurve[i] = new TGraph();
    sprintf(temp,"par=%.5f",myparameter[i]);
    myAbsorptionCurve[i]->SetName(temp);

    myAttenuationCurve[i] = new TGraph();
    sprintf(temp,"p=%.5f",myparameter[i]);
    myAttenuationCurve[i]->SetName(temp);

    sprintf(temp,"OPTICS_water_par_%.0f.ratdb",myparameter[i]*10000);
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
  double absorption;
  fstream fout2("parameter_maxL_maxLWave.txt",ios::out);
  while(fin2>>wavelength>>absorption){
	  cout<<wavelength<<"\t"<<absorption<<endl;
    //wavelength = constant/eV;
    gDayaBay->SetPoint(n,wavelength, absorption ); // absorption in 1/cm
	if(absorption!=0)
      gDayaBay_attL->SetPoint(n,wavelength, 1.0/absorption*0.01 ); // convert to attenuation length in m

    ostringstream s;
    s << wavelength;
	string mystring = s.str();
	if(mystring.find(".")!=std::string::npos)
      waveNumStr.insert(waveNumStr.length(),s.str()+", ");
    else waveNumStr.insert(waveNumStr.length(),s.str()+".0, ");
    //wavelengthStr = wavelengthStr + s.str() + ", "; // std::to_string(wavelength);//s.c_str();
    //wavelengthStr = wavelengthStr + ", ";

    for(int i=0; i<N; i++){
      myshift = myparameter[i]/TMath::Power(wavelength/200.0,4);
      anew = absorption + myshift;
	  //anew = 4.425e-5 + myshift;
      ltot = 1./anew;
	  if(ltot>myMaxAttL[i]){ myMaxAttL[i] = ltot;myWaveAtMaxAttL[i]=wavelength;}
      double attenuation = ltot*10.0;// convert cm to mm, rat-pac needs mm for attenuation length
      myAbsorptionCurve[i]->SetPoint(n,wavelength,1./ltot);
      myAttenuationCurve[i]->SetPoint(n,wavelength,ltot*0.01); // use m for attenuation length to plot
      ostringstream ss;
      ss << attenuation;
	  mystring = ss.str();
	  if(mystring.find(".")!=std::string::npos)
        dataNumStr[i].insert(dataNumStr[i].length(),ss.str()+", ");
	  else dataNumStr[i].insert(dataNumStr[i].length(),ss.str()+".0, ");
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
	
	gMaxAttL_vs_parameter->SetPoint(i,myMaxAttL[i]/100.0,myparameter[i]);
	gMaxAttL_vs_waveLength->SetPoint(i,myMaxAttL[i]/100.0,myWaveAtMaxAttL[i]);
	fout2<<myparameter[i]<<"\t"<<myMaxAttL[i]/100.0<<"\t"<<myWaveAtMaxAttL[i]<<endl;
  }
/*
  //int n=N;
  for(double wave=200; wave<380; wave+=0.5){
	  for(int i=0; i<N; i++){
	    double a_shift = myparameter[i]/TMath::Power(wave/200.0,4);
		ltot = 1./a_shift;
		//if(ltot>myMaxAttL[i]) myMaxAttL[i] = ltot;
        double attenuation = ltot*10.0;// convert cm to mm, rat-pac needs mm for attenuation length
        myAbsorptionCurve[i]->SetPoint(n,wave,1./ltot);
        myAttenuationCurve[i]->SetPoint(n,wave,ltot*0.01); // use m for attenuation length to plot
        ostringstream ss;
        ss << attenuation;
	    mystring = ss.str();
	    if(mystring.find(".")!=std::string::npos)
          dataNumStr[i].insert(0,ss.str()+", ");
	    else dataNumStr[i].insert(0,ss.str()+".0, ");
		//gMaxAttL_vs_parameter->SetPoint(n, myparameter[i],my
	  }
	  n++;
  }
*/
  // absorption curves
  TCanvas* c = new TCanvas();
  gPad->SetLogy();
  TGraph* g0 = new TGraph();
  g0->SetPoint(0,200,0);
  g0->SetPoint(1,800,0.1);
  g0->Draw("AP");
  g0->GetXaxis()->SetRangeUser(200,800);
  g0->GetYaxis()->SetRangeUser(0.00001,0.1);
  g0->GetXaxis()->SetTitle("Wavelength (nm)");
  g0->GetYaxis()->SetTitle("Absorption length (1/cm)");
  TLegend* leg = new TLegend(0.1,0.6,0.6,0.9);

  //gDayaBay->Draw("P");
  //gDayaBay->SetMarkerStyle(24);
  //gDayaBay->SetMarkerColor(1);
  //leg->AddEntry(gDayaBay,"Ref (par=0)","lp");
  leg->SetHeader("maxAttL, at wavelength");
  for(int i=0; i<N; i=i+8){
    myAbsorptionCurve[i]->Draw("P");
    myAbsorptionCurve[i]->SetMarkerStyle(21+i/3);
    myAbsorptionCurve[i]->SetMarkerColor(1+i%3);
	sprintf(temp,"%.2f m, %.1f nm",myMaxAttL[i]/100.0, myWaveAtMaxAttL[i]);
    leg->AddEntry(myAbsorptionCurve[i],temp,"lp");
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

  //gDayaBay_attL->Draw("P");
  //gDayaBay_attL->SetMarkerStyle(24);
  //gDayaBay_attL->SetMarkerColor(1);
  //leg2->AddEntry(gDayaBay_attL,"Ref (p=0)","lp");
  leg2->SetHeader("maxAttL, at wavelength");
  for(int i=0; i<N; i=i+15){
	 //if(myAttenuationCurve[i]->GetName()=="0.00170" || myAttenuationCurve[i]->GetName()=="0.01450"){
    myAttenuationCurve[i]->Draw("P");
    myAttenuationCurve[i]->SetMarkerStyle(21+i/15);
    myAttenuationCurve[i]->SetMarkerColor(1+i/15);
	sprintf(temp,"%.2f m, %.1f nm",myMaxAttL[i]/100.0, myWaveAtMaxAttL[i]);
    leg2->AddEntry(myAttenuationCurve[i],temp,"lp");
	 //}
  }
  leg2->SetNColumns(2);
  leg2->Draw();

  TCanvas* c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLogx();gPad->SetLogy();
  gPad->SetGridx();gPad->SetGridy();
  gMaxAttL_vs_parameter->Draw("AP");
  gMaxAttL_vs_parameter->SetMarkerStyle(8);
  gMaxAttL_vs_parameter->GetYaxis()->SetTitle("parameter");
  gMaxAttL_vs_parameter->GetXaxis()->SetTitle("Max. Att Length (m)");
  c3->cd(2);
  gPad->SetLogx();//gPad->SetLogy();
  gPad->SetGridx();gPad->SetGridy();
  gMaxAttL_vs_waveLength->Draw("AP");
  gMaxAttL_vs_waveLength->SetMarkerStyle(8);
  gMaxAttL_vs_waveLength->GetYaxis()->SetTitle("wavelength (nm)");
  gMaxAttL_vs_waveLength->GetXaxis()->SetTitle("Max. Att Length (m)");
  //exit(0);
}


