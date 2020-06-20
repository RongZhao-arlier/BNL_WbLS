double PI = TMath::Pi();

double IdealResponse(double *x,double *par){
    double mu = par[0];
    double q = par[1];
    double sigma = par[2];
    double amplitude = par[3];
    double sum=0;
    for(Int_t n=1; n<50; n++){
        sum += TMath::Power(mu,n)*TMath::Exp(-1.0*mu)/TMath::Factorial(n)*TMath::Exp(-1.0*(x[0]-q*n)*(x[0]-q*n)/(2.0*n*sigma*sigma))/(sigma*TMath::Sqrt(2.0*PI*n));
    }
    return sum*amplitude;
}

void drawIdeal(){
	TF1* Fideal = new TF1("Fideal",IdealResponse,0,500,4);
    Fideal->SetParNames("Npe","Peak","Width","Amplitude");
    Fideal->SetLineColor(2); Fideal->SetLineStyle(1);
    Fideal->SetParameter(0,-0.0005349);
    Fideal->SetParameter(1,166.6);
	Fideal->SetParameter(2,-47.77);
	Fideal->SetParameter(3,2.329*1E7);
	Fideal->Draw();
}

void produceOneSpeSpectrum(fstream &fout, double * par){ // produce PMT spe spectrum given fitted parameters
    TF1* Fideal = new TF1("Fideal",IdealResponse,0,500,4);
    Fideal->SetParNames("Npe","Peak","Width","Amplitude");
	Fideal->SetParameter(0,par[0]);
    Fideal->SetParameter(1,par[1]);
	Fideal->SetParameter(2,par[2]);
	Fideal->SetParameter(3,par[3]);
	//cout<<par[0]<<"\t"<<par[1]<<endl;
	for(double adc = 1; adc<500.0; adc+=1){
		fout<<fixed<<setprecision(6)<<Fideal->Eval(adc)/5000.0<<", ";
	}
	fout<<endl;
}

void produceAllSpeSpectrumForMC(){
	fstream fin("./PMTSpeSpectraCalibration/Oneton_PMT_spe_water.txt",ios::in);
	string strPMT;
	double par[4], parerr[4], chi2, ndf, prob;
	char name[100];
	for(int i=0; i<8; i++){
		sprintf(name,"GeneratedSPE_PMTS%d.txt",i);
		fstream fout(name,ios::out);
		for(double adc=1.0; adc<500.0; adc+=1.0){
			fout<<fixed<<adc<<", ";
		}
		fout<<endl;
		fin>>strPMT;
		for(int j=0;j<4;j++) {fin>>par[j]>>parerr[j]; }
		fin>>chi2>>ndf>>prob;
		produceOneSpeSpectrum(fout,par);
		fout.close();
	}
	fin.close();

	TGraph* gSpe[8];
	double a, b;
	for(int i=0; i<8; i++){
		sprintf(name,"GeneratedSPE_PMTS%d.txt",i);
		fstream fin2(name,ios::in);
		string aline1, aline2;
		getline(fin2,aline1);
		getline(fin2,aline2);
		fin2.close();
		
		vector<double> val1;
		vector<double> val2;

        std::string delimiter = ", ";

        size_t pos = 0;
        std::string token;
        while ((pos = aline1.find(delimiter)) != std::string::npos) {
            token = aline1.substr(0, pos);
            aline1.erase(0, pos + delimiter.length());
	        val1.push_back(atof(token.c_str()));
        }
		pos=0;
		while ((pos = aline2.find(delimiter)) != std::string::npos) {
            token = aline2.substr(0, pos);
            aline2.erase(0, pos + delimiter.length());
	        val2.push_back(atof(token.c_str()));
        }
		
		gSpe[i] = new TGraph();
		for(int j=0; j<val1.size(); j++){
			gSpe[i]->SetPoint(j,val1[j],val2[j]);
		}
		
	}
	TLegend* leg = new TLegend(0.8,0.6,0.9,0.9);
	leg->SetNColumns(2);
	for(int i=0; i<8; i++){
		if(i==0) gSpe[i]->Draw("AL");
		else gSpe[i]->Draw("L");
		gSpe[i]->SetLineColor(i+1);
		gSpe[i]->SetLineWidth(2);
		sprintf(name,"PMTS%d",i);
		leg->AddEntry(gSpe[i], name,"l");
	}
	leg->Draw();
}

void fitIdeal(fstream &fout, TH1F *adc, int iii){
      gStyle->SetOptFit(1111);
      TF1* Fideal = new TF1("Fideal",IdealResponse,0,500,4);
      Fideal->SetParNames("Npe","Peak","Width","Amplitude");
      Fideal->SetLineColor(2); Fideal->SetLineStyle(1);

      Fideal->SetParameter(0,0.1);
	  Fideal->SetParLimits(0,0,0.5);
	  Fideal->SetParLimits(1,100,210);
	  Fideal->SetParLimits(2,10,150);
	  
      if(iii==0){
          Fideal->SetParameter(1,143);
          Fideal->SetParameter(2,51);
      }
      if(iii==1){
          Fideal->SetParameter(1,171);
          Fideal->SetParameter(2,57);
      }
      if(iii==2){
          Fideal->SetParameter(1,178);
          Fideal->SetParameter(2,82);
      }
      if(iii==3){
          Fideal->SetParameter(1,177);
          Fideal->SetParameter(2,54);
      }
      if(iii==4){
          Fideal->SetParameter(1,166);
          Fideal->SetParameter(2,52);
      }
      if(iii==5){
          //Fideal->SetParameter(0,0.8);
          Fideal->SetParameter(1,166);
          Fideal->SetParameter(2,50);
      }
	  if(iii==6){
          Fideal->SetParameter(1,140);
          Fideal->SetParameter(2,70);
      }
	  if(iii==7){
          Fideal->SetParameter(1,135);
          Fideal->SetParameter(2,47);
      }
 
      double amp = 0;
      double par[4];
      for(Int_t i=1;i<adc->GetNbinsX();i++)
            amp += adc->GetBinContent(i+1);
      Fideal->SetParameter(3,amp);
      //adc->Fit("Fideal","RQ");

      for(int i=0;i<3;i++){
          //if(iii !=5) 
          adc->Fit("Fideal","RQ","",par[1]-1.1*par[2],500); //par[1]-1.2*par[2]
          //else adc->Fit("Fideal","R","",60,1000); //par[1]-1.2*par[2]
          Fideal->GetParameters(par);
          Fideal->SetParameters(par);
      }
      Fideal->GetParameters(par);
	  /*
	  if(par[1]<50){
		  cout<<"Problem in fit, let me refit..."<<endl;
		  Fideal->SetParameter(1,150);
		  Fideal->SetParameter(2,50);
		  Fideal->SetParameter(0,0.1);
		  adc->Fit("Fideal","RQ","",100,500);
	  }
	  Fideal->GetParameters(par);
	  */
      fout<<adc->GetName()<<"\t";
      for(Int_t j=0;j<4;j++){
          fout<<par[j]<<"\t"<<Fideal->GetParError(j)<<"\t";
      }
	  fout<<Fideal->GetChisquare()<<"\t"<<Fideal->GetNDF()<<"\t"<<Fideal->GetProb();
      fout<<endl;
      adc->GetXaxis()->SetRangeUser(0,500);
      adc->Draw();
}

