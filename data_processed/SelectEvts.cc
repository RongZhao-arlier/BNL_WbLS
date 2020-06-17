#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TLine.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TTimeStamp.h>
#include <TDatime.h>
#include <TMath.h>
#include <TPaveStats.h>
#include <TMinuit.h>
#include <stdlib.h>  
using namespace std;
//double getchi2(double factor);
double getchi2(Int_t currentnb, Double_t pmtcf);
//TH1F *hChargeMCtemp[8];
char tempnamemc[100];
double minchi2=1000000;
int nscan=10000;
int pointscnt=0;
int nevent=0;
// sprintf(tempnamemc,"hCharge_50ns_S%d_allTrigMC",0);
TH1F *hChargeMCtemp = new TH1F("tempnamemc","",50,0,50);
TGraphErrors *fitgraph=new TGraphErrors();
void fcn(Int_t &npar, Double_t *deriv, Double_t &f, Double_t *par, Int_t iflag);
void fcnall(Int_t &npar, Double_t *deriv, Double_t &f, Double_t *par, Int_t iflag);
int currentpmt=0;
const int pmtnb = 8;
//double pmtcf[pmtnb]={0.9622,0.8,0.386,0.669,1.578,1.63,.966,.68};
double pmtcfnew[pmtnb]={0.926,0.774,0.375,0.683,1.584,1.74,.919,.697};
//double pmtcf[pmtnb]={.896,0.689,0.457,0.818,1.57,1.69,0.909,.642};
//double pmtcf[pmtnb]={1.01,0.85,0.432,0.741,.82,1.52,1.01,.808};
double pmtcf[pmtnb]={2.199,2.06,1.204,1.979,2.029,2.199,2.192,2.004};
//double pmtcf[pmtnb]={1.98,1.497,1.051,0.869,1.584,1.74,1.056,.836};
//double pmtcfnew[pmtnb]={0.9622,0.8,0.386,0.669,1.578,1.63,.966,.68};
double pmtchi2new[pmtnb]={0};
TH1F *hChargeData[pmtnb];
double allpmtnpe[pmtnb][2000000]={0};
double allpmtnpedata[pmtnb][2000000]={0};
int main(){
//std::string rootfilename="./ResultsV16_allTrigCombinations_run19885_run20387_20190605.root";
std::string rootfilename="./ResultsV16_allTrigCombinations_run19885_run28655_20190605.root";  //the  data file
//std::string rootfilenamenpe="/gpfs01/lbne/users/zhaor/rat-pac/data/testana/PMTnpe26.root";
std::string rootfilenamenpe="/gpfs01/lbne/users/zhaor/rat-pac/data/testana/PMTnpe20m.root";   //the MC file
//std::string rootfilenamenpe="/gpfs01/lbne/users/zhaor/rat-pac/data/testana/PMTnpe10m.root";
const int nbTrig=27;
char dirnames[100];
char tempname[100];
char pngname[100];
string sTriggerType[nbTrig]={"101000", "101010", "101001",   //1010
                           "100100", "100110", "100101",    //1001
                           "011000", "011010", "011001",    //0110
                           "010100", "010110", "010101",    //0101
                           "111000", "111010", "111001",    //1110
                           "110100", "110110", "110101",    //1101
                           "011100", "011110", "011101",    //0111
                           "101100", "101110", "101101",    //1011
                           "111100", "111110", "111101"    //1111
                          };
  TFile* rtfile = new TFile(rootfilename.c_str(),"read");
  TFile* npefile = new TFile(rootfilenamenpe.c_str(),"read");
  TTree* npeTree = nullptr;
  npefile->GetObject("Treenpe",npeTree);
  TTree* datanpeTree = nullptr;
//the data tree and braches
//double datapmtnpe[8][20];
//int trigtype,trigtype_multi,trigtype_hodo,trigtype_led;
//  rtfile->GetObject("OneTonEvent",datanpeTree);
//  TBranch *datanpebranch  = datanpeTree->GetBranch("DigitizerPulseCharge");
//  TBranch *tgbranch  = npeTree->GetBranch("TrigType");
//  TBranch *tg_mbranch  = npeTree->GetBranch("TrigTypeFlag_Multi");
//  TBranch *tg_hbranch  = npeTree->GetBranch("TrigTypeFlag_Hodo");
//  TBranch *tg_lbranch  = npeTree->GetBranch("TrigTypeFlag_Led");
//cout<<"testtttttttttttttt"<<endl;
//  datanpebranch->SetAddress(datapmtnpe);
//  tgbranch->SetAddress(&trigtype);
//  tg_mbranch->SetAddress(&trigtype_multi);
//  tg_hbranch->SetAddress(&trigtype_hodo);
//  tg_lbranch->SetAddress(&trigtype_led);
//  int data_nevent = datanpeTree->GetEntries();
//    for (Int_t i=0;i<data_nevent;i++) {
//      //     bntrac//k->GetEvent(i);
//      datanpeTree->GetEntry(i);
//	if(trigtype_hodo!=1)continue;  //select the hodo trig
//	for(int j=0;j<pmtnb;j++)allpmtnpedata[j][i]=datapmtnpe[j][0];	
//	}

//the MC 
  TBranch *npebranch  = npeTree->GetBranch("pmtnpes");
  double pmtnpe[pmtnb]={0};
  npebranch->SetAddress(pmtnpe);
  nevent = npeTree->GetEntries();
//nevent=100000;
  cout<<"nevent in MC is:"<<nevent<<endl;
    for (Int_t i=0;i<nevent;i++) {
      //     bntrac//k->GetEvent(i);
      npeTree->GetEntry(i);
      for(int j=0;j<pmtnb;j++)
	{allpmtnpe[j][i]=pmtnpe[j];}
//	if(pmtnpe[0]==0){cout<<i<<": ";for(int mm=0;mm<8;mm++){cout<<pmtnpe[mm]<<" ";}cout<<endl;}
    }
  TDirectory* dir_TrigType[nbTrig];
  TH1F *hCharge[pmtnb];
  TH1F *hChargeMC[pmtnb];
  TH1F *H_npe[pmtnb][nbTrig]; 
  for(int j=0;j<pmtnb;j++){
             sprintf(tempname,"hCharge_50ns_S%d_Trig",j);
     hCharge[j] = new TH1F(tempname,"",50,0,50);
             sprintf(tempname,"hCharge_50ns_S%d_allTrig",j);
     hChargeData[j] = new TH1F(tempname,"",50,0,50);
             sprintf(tempname,"hCharge_50ns_S%d_allTrigMC",j);
     hChargeMC[j] = new TH1F(tempname,"",50,0,50);
        for(int i=0; i<nbTrig; i++){
            sprintf(dirnames,"dir_%s/hCharge_50ns_S%d_%sTrig",sTriggerType[i].c_str(),j,sTriggerType[i].c_str());
    	//cout<<dirnames<<endl;
    	//     dir_TrigType[i] = rtfile->GetObject(dirnames);
     rtfile->GetObject(dirnames,H_npe[j][i]);
	if(i<28)  //event slection by type tag
	hCharge[j]->Add(H_npe[j][i]);
//	if(i>24)cout<<"npecounts is:"<<H_npe[j][i]->GetEntries();
	}
    cout<<"entries of hcharge is:"<<hCharge[j]->GetEntries()<<endl;;
    const int nbevt=100000;
   // for(int ii=0;ii<nbevt;ii++){
   int mc_cnt=0,read_cnt=0;
   bool multag=1;
    while(mc_cnt<nbevt&&read_cnt<1e7){
	multag=1;
   for(int tc=0;tc<pmtnb;tc++){if(allpmtnpe[tc][read_cnt]*pmtcf[tc]<1)multag=0;}
//	cout<<multag<<endl;
    //cout<<hCharge[j]->GetRandom()<<endl;
//    if(allpmtnpe[j][read_cnt]*pmtcf[j]>0){     //the threshold of 0-1 p.e   //rong zhao, the cut of 0 p.e 
    if(allpmtnpe[j][read_cnt]*pmtcf[j]>.5&&multag!=0){     //the threshold of 0-1 p.e   //rong zhao, the cut of 0 p.e 
	hChargeMC[j]->Fill(allpmtnpe[j][read_cnt]*pmtcf[j]);
	hChargeData[j]->Fill(hCharge[j]->GetRandom());
		mc_cnt++;
    		}
	read_cnt++;
	}
	cout<<"rong zhao test:"<<mc_cnt++<<endl;
    for(int k=5;k<50;k++){
	if(hChargeData[j]->GetBinContent(k)<0*nbevt*1./50000){  //for the last bin
	int resbinct_data=0;
	int resbinct_mc=0;
	for(int m=k;m<51;m++){resbinct_data+=hChargeData[j]->GetBinContent(m);hChargeData[j]->SetBinContent(m,0);}
	for(int m=k;m<51;m++){resbinct_mc+=hChargeMC[j]->GetBinContent(m);hChargeMC[j]->SetBinContent(m,0);}
	hChargeData[j]->SetBinContent(k,resbinct_data*1);
	hChargeMC[j]->SetBinContent(k,resbinct_mc*1);
	break;}
	}
  }
  //sampling the histogram
  //        tree = (TTree*)rootfile->Get("OneTonEvent");
  //        bEvtNumber = tree->GetBranch("EvtNumber");
  TString outpdf=("./pmtnpe.pdf");
  TString outpdf_start=("./pmtnpe.pdf[");
  TString outpdf_end=("./pmtnpe.pdf]");
  TCanvas *can1=new TCanvas("can1","canvas1",1200,900);
  can1->cd();
  can1->Print(outpdf_start);
  //hChargeMC[0]->Draw();
  for(int hc=0;hc<pmtnb;hc++){
  can1->cd();
TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
pad1->SetBottomMargin(0.00001);
pad1->SetBorderMode(0);
//pad1->SetLogy();
pad2->SetTopMargin(0.00001);
pad2->SetBottomMargin(0.3);
//pad2->SetLeftMargin(0.0001);
//pad2->SetBorderMode(0);
pad1->Draw();
pad2->Draw();
pad1->cd();
  hChargeData[hc]->Draw();
  //hCharge[hc]->Draw("sames");
  hChargeData[hc]->SetLineWidth(2);
  hChargeData[hc]->SetLineColor(2);
  hChargeData[hc]->GetXaxis()->SetTitle("n.p.e");
  hChargeData[hc]->GetYaxis()->SetRangeUser(0,hChargeData[hc]->GetMaximum()*1.5);
//  gPad->Update();   //update to access stat box
//  TPaveStats *ps2 = (TPaveStats*) gPad->GetPrimitive("stats");
  hChargeMC[hc]->Draw("sames");
  hChargeMC[hc]->SetLineColor(kTeal+3);
  hChargeMC[hc]->SetLineWidth(2);
can1->Update();
  gPad->Update();   //update to access stat box
//  TPaveStats *ps = (TPaveStats*) gPad->GetPrimitive("stats");
TPaveStats* ps2= (TPaveStats*)(hChargeData[hc]->FindObject("stats"));
  ps2->SetX1NDC(0.7);
  ps2->SetX2NDC(0.9);
  ps2->SetY1NDC(0.3);
  ps2->SetY2NDC(0.5);
TPaveStats* ps= (TPaveStats*)(hChargeMC[hc]->FindObject("stats"));
  ps->SetX1NDC(0.7);
  ps->SetX2NDC(0.9);
  ps->SetY1NDC(0.7);
  ps->SetY2NDC(0.9);
    TLegend *leg = new TLegend(0.7,0.6,0.9,0.7);
    leg->AddEntry(hChargeData[hc],"Data","l");
    leg->AddEntry(hChargeMC[hc],"MC by rat","l");
    leg->Draw();
/*************the second pad for difference of data and MC****************/
pad2->cd();
	TGraphErrors *datamcgraph=new TGraphErrors();
	for(int binnum=1;binnum<51;binnum++){
	if(hChargeData[hc]->GetBinContent(binnum)<1)continue;
//	datamcgraph->SetPoint(binnum,binnum ,(hChargeMC[hc]->GetBinContent(binnum)-hChargeData[hc]->GetBinContent(binnum)));
	datamcgraph->SetPoint(binnum,binnum ,TMath::Power((hChargeData[hc]->GetBinContent(binnum)-hChargeMC[hc]->GetBinContent(binnum)),1)*1./TMath::Sqrt((hChargeData[hc]->GetBinContent(binnum)+1e-7)+hChargeMC[hc]->GetBinContent(binnum)));
double ddata=(hChargeData[hc]->GetBinContent(binnum)+1e-15);
double dmc=hChargeMC[hc]->GetBinContent(binnum);
//	datamcgraph->SetPointError(binnum,.5,TMath::Sqrt(dmc*1./ddata*(ddata+dmc)/(ddata*ddata)));
	datamcgraph->SetPointError(binnum,.5,(TMath::Sqrt(ddata)+TMath::Sqrt(dmc))*1./TMath::Sqrt(ddata+dmc));
	  }
gPad->Update();
datamcgraph->Draw("AP");
gPad->Update();
 datamcgraph->GetXaxis()->SetTitle("number of p.e");
// datamcgraph->GetXaxis()->SetRangeUser(0,50);
 datamcgraph->GetXaxis()->SetLimits(0,50);
//TStyle* mcStyle = new TStyle("mcStyle","Manuel's RootStyles"); 
//mcStyle->SetLabelSize(0.03,"xyz");
 datamcgraph->GetXaxis()->SetLabelSize(.07);
 datamcgraph->GetYaxis()->SetLabelSize(.07);
 datamcgraph->GetXaxis()->SetTitleSize(.0850);
 datamcgraph->GetYaxis()->SetTitleSize(.0850);
 datamcgraph->GetYaxis()->SetTitleOffset(.320);
can1->Update();
// datamcgraph->GetYaxis()->SetTitleSize(20);
 datamcgraph->GetYaxis()->SetTitle( "(MC-DATA)/#sigma");
  can1->SaveAs(outpdf);
  sprintf(pngname,"figures/PMTnperes%d.png",hc);
  can1->SaveAs(pngname);
delete datamcgraph;
  }
  can1->Print(outpdf_end);
  //hChargeMC[0]->GetYaxis()->SetRangeUser(0,2000);
  //can1->SaveAs("test.png");
//  const int npar = 1;              // the number of parameters
  const int npar = 6;              // the number of parameters
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);
//  double par[npar];               // the start values
//  double stepSize[npar];          // step sizes 
//  double minVal[npar];            // minimum bound on parameter 
//  double maxVal[npar];            // maximum bound on parameter
//  string parName[npar];
//double   par[npar] ={.9,.8,.4,.6,.9,.7};            // a guess
double   par[npar] ={.7,.7,.7,.7,.7,.7};            // a guess
double   stepSize[npar]={0.02,0.02,0.02,0.02,0.02,0.02};       // take e.g. 0.1 of start value
//double   minVal[npar] = {0.2,0.2,0.02,0.2,0.2,0.2};   // if min and max values = 0, parameter is unbounded.
double   minVal[npar] = {0.,0.,0.0,0.,0.,0.};   // if min and max values = 0, parameter is unbounded.
//double   maxVal[npar] = {1.2,1.2,1.,1.,1.2,1.1};
//double   maxVal[npar] = {12000,12000,12000,12000.,12000,11000};
double   maxVal[npar] = {1e7,1e7,1e7,1e7,1e7,1e7};
string   parName[npar]= {"s0","s1","s2","s3","s6","s7"};
  double arglist[10];
  int ierflg = 1;
  	for (int i=0; i<npar; i++){
  	minuit.DefineParameter(i, parName[i].c_str(), 
  	par[i], stepSize[i], minVal[i], maxVal[i]);
  	}
  arglist[0] = 0;
//  arglist[0] = .5;
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
//for(int ii=0;ii<pmtnb;ii++)
for(int ii=0;ii<8;ii++)
    {  //loop of fitting
	currentpmt=ii;
	minchi2=1000000;
	//currentpmt=1;
	// minuit.mnexcm("SIMplex ", arglist ,1,ierflg);
	//minuit.SetMaxIterations(5000);
	arglist[0] = 0.5;
  minuit.SetFCN(fcnall);
  	for (int i=0; i<npar; i++){
  	minuit.DefineParameter(i, parName[i].c_str(), 
  	par[i], stepSize[i], minVal[i], maxVal[i]);
  	}
//	minuit.Migrad();       // Minuit's best minimization algorithm
//	minuit.mnexcm("SCAN", arglist,1,ierflg);
	for(int k=0;k<nscan;k++){getchi2(currentpmt,.6+k*0.0002);}
//minuit.mnsimp();
	double outpar[npar], err[npar];
	for(int i=0; i<npar; i++){
    	minuit.GetParameter(i,outpar[i],err[i]);
	}
	//const int N_pars=1;
	//TMinuit *gMinuit = new TMinuit(N_pars);  //initialize TMinuit with a maximum of N_pars params
	//gMinuit->SetFCN(fcn);
	//gMinuit->DefineParameter(i, parName[i].c_str(), vstart[i],step[i], minVal[i], maxVal[i]);
	//Double_t arglist[10];
	//Int_t ierflg = 0;
	//arglist[0] = 1;
	//gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	//gMinuit->DefineParameter(0, "cf0", 1,0.001, 0.01, 10);
	//arglist[0] = 0;
	//arglist[1] = 1.;
	//gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	//gMinuit->SetMaxIteration(500);
	//gMinuit->Migrad();
	//gMinuit->Simplex();
	//arglist[0] = 0;
	//gMinuit->mnexcm("SCAN", arglist,1,ierflg);
  	minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  	cout<<amin<<"\t"<<edm<<"\t"<<errdef<<"\t"<<nvpar<<"\t"<<nparx<<"\t"<<icstat<<endl;
  	double outpari[npar], errors[npar];
  	cout<<"Output fit parameters: "<<endl;
  	cout<<"FCN "<<amin<<endl;
  	for(int i=0; i<npar; i++){
  	    minuit.GetParameter(i, outpar[i], errors[i]);
  	    cout<<i<<"\t"<<"factor0"<<"\t"<<outpar[i]<<"\t"<<errors[i]<<endl;
  	}
	cout<<"pmtcfnew:"<<currentpmt<<":"<<pmtcfnew[currentpmt]<<endl;
	TCanvas *can2=new TCanvas("can2","canvas2",1200,900);
	can2->cd();
	fitgraph->Draw("ap*");
	fitgraph->GetXaxis()->SetRangeUser(pmtcfnew[ii]-0.05,pmtcfnew[ii]+0.05);
	fitgraph->GetYaxis()->SetRangeUser(pmtchi2new[ii]-5,pmtchi2new[ii]+100);
	TLine *sl=new TLine(pmtcfnew[ii]-0.05,pmtchi2new[ii]+1,pmtcfnew[ii]+0.05,pmtchi2new[ii]+1);
	 sl->SetLineColor(kBlue);
	 sl->Draw();
	TLine *slcf=new TLine(pmtcfnew[ii]-0.,pmtchi2new[ii]-5,pmtcfnew[ii]+0.,pmtchi2new[ii]+10);
	 slcf->SetLineColor(kRed);
	 slcf->Draw();
//	fitgraph->GetXaxis()->SetRangeUser(pmtcf[ii]-0.2,pmtcf[ii]+0.2);
//	fitgraph->GetYaxis()->SetRangeUser(0,pmtchi2new[ii]+15);
	fitgraph->GetYaxis()->SetTitle("#chi ^{2}");
	fitgraph->GetXaxis()->SetTitle("calibration factor");
	//can2->SetLogy();
	gPad->SetGrid();
         TLegend *leg = new TLegend(0.1,0.1,0.3,0.2);
         leg->AddEntry(fitgraph,Form("#p.e of PMTs%d",ii),"l");
         leg->Draw();
	can2->SaveAs(Form("figures/chi2_%d.eps",ii));
	delete can2;
    }//loop of fitting

return 0;
}

void fcn(Int_t &npar, Double_t *deriv, Double_t &f, Double_t *par, Int_t iflag)
{
cout<<"par0:"<<par[0]<<endl;
//f=getchi2(pmtcf[0]);
//}
//double getchi2(double factor){
double chi2=0;
int ndmcnt=0;
for(int i=0;i<nevent;i++){
hChargeMCtemp->Fill(allpmtnpe[currentpmt][i]*par[0]);
}
double scalefactor=(hChargeData[currentpmt]->GetEntries()*1./nevent);
hChargeMCtemp->Scale(scalefactor);
for(int j=0;j<50;j++){
//for(int j=0;j<(hChargeData[0]->GetSize()-2);j++){
//if((hChargeData[0]->GetBinContent(j)+hChargeMC[0]->GetBinContent(j))!=0)chi2+=TMath::Power(hChargeData[0]->GetBinContent(j)-hChargeMC[0]->GetBinContent(j),2)/(1.);
if((hChargeData[currentpmt]->GetBinContent(j+1)+0*hChargeMCtemp->GetBinContent(j+1))!=0)
 {chi2+=TMath::Power(hChargeData[currentpmt]->GetBinContent(j+1)-hChargeMCtemp->GetBinContent(j+1),2)/(hChargeData[currentpmt]->GetBinContent(j+1)+hChargeMCtemp->GetBinContent(j+1)*scalefactor);
//cout<<j<<" "<<hChargeData[0]->GetBinContent(j)<<" "<<hChargeMCtemp->GetBinContent(j)<<endl;
  ndmcnt++;
 }
}
//cout<<"ndmcntis:::::::"<<ndmcnt<<endl;
//f=chi2*1./ndmcnt;
f=chi2*1.;
//cout.precision(12);
//cout<<"f is:"<<f*1<<endl;
//int nscount=pointscnt%nscan;
//fitgraph->SetPoint(nscount, par[0],f);
//fitgraph->SetPointError(nscount,0.002,0);
//return  chi2;
hChargeMCtemp->Reset("ICESM");
//pointscnt++;
//delete hChargeMCtemp;
}
double getchi2(Int_t currentnb, Double_t pmtcf)
{
	double chi2=0;
	int ndmcnt=0;
	//for(int i=0;i<hChargeData[currentnb]->GetEntries();i++){
	for(int i=0;i<nevent;i++){
hChargeMCtemp->Fill(allpmtnpe[currentnb][i]*pmtcf);
	}
double scalefactor=(hChargeData[currentpmt]->GetEntries()*1./nevent);
hChargeMCtemp->Scale(scalefactor);
	for(int j=0;j<50;j++){
//		if((hChargeData[currentpmt]->GetBinContent(j+1))<30) break;
		//for(int j=0;j<(hChargeData[0]->GetSize()-2);j++){
		//if((hChargeData[0]->GetBinContent(j)+hChargeMC[0]->GetBinContent(j))!=0)chi2+=TMath::Power(hChargeData[0]->GetBinContent(j)-hChargeMC[0]->GetBinContent(j),2)/(1.);
		if((hChargeData[currentnb]->GetBinContent(j+1)+0*hChargeMCtemp->GetBinContent(j+1))>=1){
	chi2+=TMath::Power(hChargeData[currentnb]->GetBinContent(j+1)-hChargeMCtemp->GetBinContent(j+1),2)*1./(hChargeData[currentnb]->GetBinContent(j+1)+hChargeMCtemp->GetBinContent(j+1)*scalefactor);
		//cout<<j<<" "<<hChargeData[0]->GetBinContent(j)<<" "<<hChargeMCtemp->GetBinContent(j)<<endl;
		ndmcnt++;
		}
//		else {
//			while(j<50){ chi2+=TMath::Power(hChargeData[currentnb]->GetBinContent(j+1)-hChargeMCtemp->GetBinContent(j+1),2)*1./(hChargeData[currentnb]->GetBinContent(j+1)+hChargeMCtemp->GetBinContent(j+1)*scalefactor);
//   j++;}
//			break;}
	    }
//	chi2=chi2*1./ndmcnt;
	if(chi2<minchi2){minchi2=chi2;pmtcfnew[currentpmt]=pmtcf;pmtchi2new[currentpmt]=minchi2;}
//	cout<<"cfnew: "<<pmtcf<<"chi2 "<<chi2<<"minchi2:"<<minchi2<<endl;}
	cout.precision(12);
	int nscount=pointscnt%nscan;
	fitgraph->SetPoint(nscount,pmtcf ,chi2);
	fitgraph->SetPointError(nscount,0.0,0);
	hChargeMCtemp->Reset("ICESM");
	pointscnt++;
	return  chi2;
	//delete hChargeMCtemp;
}
//fcn all, fitting all pmt togather
void fcnall(Int_t &npar, Double_t *deriv, Double_t &f, Double_t *par, Int_t iflag)
{
//cout<<"par0:"<<par[0]<<endl;
double chi2=0;
double scalefactor;
int ndmcnt=0;
int pmts[6]={0,1,2,3,6,7};
for(int ii=0;ii<6;ii++){
	for(int i=0;i<nevent;i++){
	hChargeMCtemp->Fill(allpmtnpe[pmts[ii]][i]*par[ii]);
	}	
	scalefactor=(hChargeData[pmts[ii]]->GetEntries()*1./nevent);
	hChargeMCtemp->Scale(scalefactor);
	for(int j=0;j<50;j++){
	if((hChargeData[pmts[ii]]->GetBinContent(j+1)+0*hChargeMCtemp->GetBinContent(j+1))!=0)
//	if((hChargeData[pmts[ii]]->GetBinContent(j+1)<10)
	 {chi2+=TMath::Power(hChargeData[pmts[ii]]->GetBinContent(j+1)-hChargeMCtemp->GetBinContent(j+1),2)/(hChargeData[pmts[ii]]->GetBinContent(j+1)+hChargeMCtemp->GetBinContent(j+1)*scalefactor);
  	ndmcnt++;
	  }
	 }
	hChargeMCtemp->Reset("ICESM");
	}

f=chi2*1.;
}

