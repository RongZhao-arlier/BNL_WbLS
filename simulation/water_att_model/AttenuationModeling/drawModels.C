double attenuationLengthModel(double* x, double* par){
  return par[0]/TMath::Power(x[0]/200.0,4);
}

void drawModels(){
	TF1* fAttL0 = new TF1("AttLength0",attenuationLengthModel,200,800,1);
	fAttL0->SetParameter(0,0.0078);
	fAttL0->Draw();
}