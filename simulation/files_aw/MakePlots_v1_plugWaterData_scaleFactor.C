#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TTimeStamp.h>
#include <TDatime.h>
#include <TMath.h>
#include <TFractionFitter.h>
#include <TObjArray.h>

#include <bits/stdc++.h> 
#include <stdlib.h>  

using namespace std;

vector<double> fitParabolicAnalyticly(TGraph* gGraph);

void AlignChi2Curves(string name);

// Chi-2 calculation
int main(){
    const int pmtnb = 8;
    const int nbTrigs = 27;
    const int nbOfMCDataPoints = 28+1; // in my example, I have 17 data points (17 runs under different attenuation parameters)
    const int file_index_as_data = nbOfMCDataPoints-1; // maximmum is nbOfMCDataPoints-1
    const int select_for_compare = 4; // select one to compare histograms of npe distribution
    string strig[nbTrigs]={"101000", "101010", "101001",
                         "100100", "100110", "100101",
                         "011000", "011010", "011001",
                         "010100", "010110", "010101",
                         "111000", "111010", "111001",
                         "110100", "110110", "110101",
                         "011100", "011110", "011101",
                         "101100", "101110", "101101",
                         "111100", "111110", "111101"
                        };
    int UseTrigIndex[10]={0, 3, 6, 9, 12, 15, 16, 17, 21, 24};
    string sMCFileNames[nbOfMCDataPoints]; 
    sMCFileNames[file_index_as_data] = "../ResultsV16_allTrigCombinations_run19885_run20890_20190605.root";
    double attenuationPar[nbOfMCDataPoints];
    attenuationPar[file_index_as_data] = 0.0;

    char make_namestr[150];
    for(int i=0; i<file_index_as_data; i++){
      sprintf(make_namestr,"Results_v12_allTrigCombinations_scanFiles4.txt_scaleFactor_%.2f.root",0.1+0.05*(i+0));
      sMCFileNames[i] = make_namestr;
      //attenuationPar[i] = 40.792 * 0.1*(i+1);
      attenuationPar[i] = 0.1+0.05*(i+0);
    }

    for(int i=0; i<nbOfMCDataPoints; i++){
      cout<<i<<" attenuation after scaled efficiency: "<<attenuationPar[i]<<endl;
    }

    // this run corresponds to 19.6837 m max. attenuation length
    // I'll see if this can be fit by scanning chi2
    if(file_index_as_data>=nbOfMCDataPoints) {cout<<"Check the file index as data."<<endl; exit(0);}
    string sDataFileName = sMCFileNames[file_index_as_data];

    // get the histograms first, for both MC and data
    TFile* MCFile[nbOfMCDataPoints];
    TFile* DataFile = new TFile(sDataFileName.c_str(),"read");
    TH1F* hMCnpe[nbOfMCDataPoints][pmtnb][nbTrigs];
    TH1F* hDatanpe[pmtnb][nbTrigs];

    TH1F* hCloneMCnpe[pmtnb][nbTrigs];

    char name[150];
    // get histograms from data 
    double Nevt_data[pmtnb][nbTrigs];
    for(int j=0; j<pmtnb; j++){
        for(int k=0; k<nbTrigs; k++){
            sprintf(name,"dir_%s/hCharge_50ns_S%d_%sTrig",strig[k].c_str(),j,strig[k].c_str());
            hDatanpe[j][k] = (TH1F*)DataFile->Get(name);
            Nevt_data[j][k] = hDatanpe[j][k]->GetEntries();
        }
    }
    // get histograms from MC, and calculate chi2
    double Nevt_MC[nbOfMCDataPoints][pmtnb][nbTrigs];
    for(int i = 0; i<nbOfMCDataPoints; i++){
        if(i==file_index_as_data) continue; // just skip the one that I use as data
        cout<<"MC point "<<sMCFileNames[i]<<endl;
        
        MCFile[i] = new TFile(sMCFileNames[i].c_str(),"read");
        TCanvas* c_temp = (TCanvas*)MCFile[i]->Get("cNumberOfEvents");

        for(int j=0; j<pmtnb; j++){
            for(int k=0; k<nbTrigs; k++){
                sprintf(name,"dir_%s/hCharge_50ns_S%d_%sTrig",strig[k].c_str(),j,strig[k].c_str());
                hMCnpe[i][j][k] = (TH1F*)MCFile[i]->Get(name);
                //cout<<hMCnpe[i][j][k]->GetEntries()<<"\t"<<hMCnpe[i][j][k]->GetNbinsX()<<endl;
                //cout<<hDatanpe[j][k]->GetEntries()<<"\t"<<hMCnpe[i][j][k]->GetEntries()<<endl;
                Nevt_MC[i][j][k] = hMCnpe[i][j][k]->GetEntries();
                hMCnpe[i][j][k]->Scale(1.0*hDatanpe[j][k]->GetEntries()/hMCnpe[i][j][k]->GetEntries());
                
                //if(i==select_for_compare){
                //    hCloneMCnpe[j][k] = (TH1F*)hMCnpe[i][j][k]->Clone();
                //}
            }
        }
        //MCFile[i]->Close();
    }

    /*
    //try the TFractionFitter method
    TH1F* hDataTotalPMTS0 = (TH1F*)hDatanpe[0][0]->Clone();
    for(int k=1;k<nbTrigs;k++) hDataTotalPMTS0->Add(hDatanpe[0][k]);
    TObjArray* oMC = new TObjArray(nbTrigs);
    //oMC->Add(hMCnpe[0][0][0]);
    for(int j=0;j<nbTrigs;j++){
        //if(j==file_index_as_data) continue; // just skip the one that I use as data
        oMC->Add(hMCnpe[6][0][j]);
    }
    TFractionFitter* fit = new TFractionFitter(hDataTotalPMTS0, oMC);
    fit->SetRangeX(1,30);
    int status_fit = fit->Fit();
    cout<<"\tfit status: "<<status_fit<<endl;
    if (status_fit == 0) {                       // check on fit status
        TH1F* result = (TH1F*) fit->GetPlot();
        TCanvas* cResult = new TCanvas();
        hDataTotalPMTS0->Draw("Ep");
        result->Draw("same");
        sprintf(name, "FractionFitResultTest/S%d.png",0);
        cResult->SaveAs(name);
        cResult->Close();
    }
    */

    
    // draw compared npe distributions
    TFile* fcompare = new TFile("npe_dis_compare_with_data_30m_scaleFactor.root","recreate");
    gStyle->SetOptStat(0);
    for(int m=0; m<nbOfMCDataPoints; m++){
        sprintf(name, "dir_fId_%d_atte_%f",m, attenuationPar[m]);
        cout<<name<<endl;
        if(m==file_index_as_data){cout<<"--that is the data, skip."<<endl;continue;}
        TDirectory* adir = fcompare->mkdir(name);
        adir->cd();
        for(int i=0; i<pmtnb; i++){
            for(int j=0;j<nbTrigs;j++){
                TCanvas* c_com = new TCanvas("name","title",1200,800);
                sprintf(name,"com_S%d_trig_%s",i,strig[j].c_str());
                //cout<<"canvas get name "<<name<<endl;
                c_com->SetName(name);
                c_com->SetTitle(name);
                c_com->Divide(1,2);
                c_com->cd(1);
                //gPad->SetPad(0,0.76,1,1);
                //gPad->SetLogy();
                hDatanpe[i][j]->Draw("ep"); hDatanpe[i][j]->SetMarkerStyle(8); hDatanpe[i][j]->SetMarkerColor(1);
                hMCnpe[m][i][j]->Draw("same ep"); hMCnpe[m][i][j]->SetMarkerStyle(8); hMCnpe[m][i][j]->SetMarkerColor(2);
                hDatanpe[i][j]->GetXaxis()->SetTitleSize(0.05);
                hDatanpe[i][j]->GetXaxis()->SetLabelSize(0.05);
                hDatanpe[i][j]->GetYaxis()->SetTitleSize(0.1);
                hDatanpe[i][j]->GetYaxis()->SetLabelSize(0.1);
                hDatanpe[i][j]->GetYaxis()->SetTitleOffset(0.3);
                sprintf(name,"S%d_trig%s, Nevt: Data %.0f / MC %.0f, scaleFactor in MC %f", i, strig[j].c_str(), Nevt_data[i][j], Nevt_MC[m][i][j],  attenuationPar[m] );
                hDatanpe[i][j]->SetTitle(name);
                TH1F* hRes = new TH1F("residual","",50,0,50);
                sprintf(name,"res_fId_%d_S%d_trigId%d",m,i,j);
                hRes->SetName(name);
                TH1F* hChiSqDevi = new TH1F("Chi2Abs","Absolute #chi^{2} value in bins",50,0,50);
                sprintf(name,"chi2abs_%d_S%d_trigId%d",m,i,j);
                hChiSqDevi->SetName(name);
                //TGraphErrors* hRes = new TGraphErrors();
                int nbin = hDatanpe[i][j]->GetNbinsX();
                for(int k=0;k<nbin;k++){
                    double a = hDatanpe[i][j]->GetBinContent(k+1);
                    double b = hMCnpe[m][i][j]->GetBinContent(k+1);
                    //cout<<k<<"\t"<<b<<"\t"<<a<<endl;
                    if(a>=10){
                        //hRes->SetPoint(k, k, b/a);
                        //hRes->SetPointError(k,0,b/a*TMath::Sqrt(1./b+1./a));
                        hRes->SetBinContent(k+1,b/a);
                        hRes->SetBinError(k+1,b/a*TMath::Sqrt(1./b+1./a));
                        hChiSqDevi->SetBinContent(k+1, (a-b)/TMath::Sqrt(a+b));
                    }
                    else{
                        a = hDatanpe[i][j]->Integral(k+1,nbin);
                        b = hMCnpe[m][i][j]->Integral(k+1,nbin);
                        if(a!=0 && b!=0){
                            hRes->SetBinContent( k+1+(nbin-(k+1))/2, a/b);
                            hRes->SetBinError(k+1+(nbin-(k+1))/2,a/b*TMath::Sqrt(1./b+1./a));
                            hChiSqDevi->SetBinContent(k+1+(nbin-(k+1))/2, (a-b)/TMath::Sqrt(a+b));
                            //hChiSqDevi->SetBinError(k+1+(nbin-(k+1))/2,0);
                        }
                        break;
                    }
                }
                /*
                c_com->cd(2);
                gPad->SetPad(0,0.34,1,0.75);
                hRes->Draw("ep"); hRes->SetMarkerStyle(8);
                hRes->GetYaxis()->SetTitle("Ratio: Data/MC");
                hRes->GetYaxis()->SetRangeUser(0,3);
                hRes->GetXaxis()->SetTitleSize(0.05);
                hRes->GetXaxis()->SetLabelSize(0.05);
                hRes->GetYaxis()->SetTitleSize(0.05);
                hRes->GetYaxis()->SetLabelSize(0.05);
                hRes->GetYaxis()->SetTitleOffset(0.3);
                */
                c_com->cd(2);
                //gPad->SetPad(0,0.0,1,0.33);
                hChiSqDevi->Draw("P"); hChiSqDevi->SetMarkerStyle(8);
                hChiSqDevi->GetYaxis()->SetTitle("#Chi^{2} / bin (D - M)");
                hChiSqDevi->GetXaxis()->SetTitleSize(0.05);
                hChiSqDevi->GetXaxis()->SetLabelSize(0.05);
                hChiSqDevi->GetYaxis()->SetTitleSize(0.1);
                hChiSqDevi->GetYaxis()->SetLabelSize(0.1);
                hChiSqDevi->GetYaxis()->SetTitleOffset(0.3);
                hChiSqDevi->GetYaxis()->SetRangeUser(-5,5);

                //if(attenuationPar[m]>=3 && attenuationPar[m]<=15){
                    //sprintf(name, "MC_Data_comp/dir_fId_%d_scaleF_%f_S%d_TrigId_%s.png",m, 0.1*(m+1),i,strig[j].c_str());
                    //c_com->SaveAs(name);
                //}
                c_com->Write();
                c_com->Close();
            }
        }
    }
    fcompare->Close();
    //fstream fout("mytestout.txt",ios::out);

    // now calculate the chi2
    TGraphErrors* gChi2 = new TGraphErrors(); // total chi2
    TGraph* gChi2_TrigPMTSelected=new TGraph(); // total chi2 in those selected trigger types and PMTs
    TGraphErrors* gChi2_indiv[pmtnb][nbTrigs]; // chi2 for each PMT for each trigger type
    TGraph* gChi2_TrigTypes[nbTrigs]; // calculate chi2 for every trig type for all PMTs
    TGraph* gChi2_PMTs[pmtnb]; // calculate chi2 for every PMT for all trig types
    double totalchi2 = 0;
    double totalchi2_selected = 0;
    double totalchi2_indiv[pmtnb][nbTrigs];
    double totalchi2_TrigTypes[nbTrigs];
    double totalchi2_PMTs[pmtnb];
    // initialize chi2
    double chi2[nbOfMCDataPoints][pmtnb][nbTrigs];
    for(int i=0; i<nbOfMCDataPoints; i++){
        for(int j=0; j<pmtnb; j++){
            for(int k=0; k<nbTrigs; k++){
                chi2[i][j][k]=0;
            }
        }
    }
    for(int j=0; j<pmtnb; j++){
        totalchi2_PMTs[j]=0.0;
        gChi2_PMTs[j] = new TGraph();
        sprintf(name,"gChi2_PMT_S%d",j);
        gChi2_PMTs[j]->SetName(name);
        for(int k=0; k<nbTrigs; k++){
            gChi2_indiv[j][k] = new TGraphErrors();
            totalchi2_indiv[j][k]=0;
        }
    }
    for(int k=0; k<nbTrigs;k++){
        totalchi2_TrigTypes[k]=0;
        gChi2_TrigTypes[k]=new TGraph();
        sprintf(name,"gChi2_TrigTypes_%s",strig[k].c_str());
        gChi2_TrigTypes[k]->SetName(name);
    }

    // the core part calculating chi2
    for(int i=0; i<nbOfMCDataPoints; i++){ // looping through MC files
        if(i==file_index_as_data){ // ignore the one that is treated as data, if the index is found,
          totalchi2 = 0;           // reset the values and then move to the next
          totalchi2_selected = 0;
          for(int j=0; j<pmtnb; j++){
              totalchi2_PMTs[j]=0;
              for(int k=0;k<nbTrigs;k++){
                  totalchi2_indiv[j][k]=0;
                  //chi2[i][j][k]=0;
              }
          }
          for(int k=0; k<nbTrigs; k++){
              totalchi2_TrigTypes[k]=0;
          }
          continue;
        }
        
        for(int j=0; j<pmtnb; j++){
            for(int k=0; k<nbTrigs; k++){
                int nbin = hMCnpe[i][j][k]->GetNbinsX();
                double ui =0, ni = 0;
                bool isUseTrig = false;
                for(int m=0;m<10;m++){
                    if(k==UseTrigIndex[m]) {isUseTrig=true;break;}
                }
                bool needs_mask = (j==4 || j==5) ;//|| (isUseTrig==false);
                //cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<needs_mask<<endl;
                int nb_fit=0;
                if(j==4 || j==5) nb_fit = 5;
                else nb_fit = 50;
                double overbinchi2 = 0;
                for(int bin=0; bin<nbin; bin++){ // count only 30 bins, the bins above will be in one bin
                    ui = hMCnpe[i][j][k]->GetBinContent(bin+1);
                    ni = hDatanpe[j][k]->GetBinContent(bin+1);
                    // Poisson chi-square
                    //if(ui!=0 && ni!=0) chi2[i][j][k] += 2.0*( ui - ni + 1.0*ni*TMath::Log(ni/ui) );
                    //if(ui!=0 && ni==0) chi2[i][j][k] += 2.0*ui;
                    //Gaussian chi-square
                    //if(ni==0 && ui==0) cout<<"ni=0 for "<<i<<"\t"<<j<<"\t"<<k<<"\t"<<bin+1<<endl;
                    if(ni>=10){
                        chi2[i][j][k] = 1.0*(ui-ni)*(ui-ni)/(ni+ui);
                        totalchi2 += chi2[i][j][k];
                        totalchi2_indiv[j][k] += chi2[i][j][k];
                        totalchi2_PMTs[j] += chi2[i][j][k];
                        if( !needs_mask ){
                            totalchi2_selected  += chi2[i][j][k];
                        }
                    }
                    else{
                        ui=hMCnpe[i][j][k]->Integral(bin+1,nbin);
                        ni=hDatanpe[j][k]->Integral(bin+1,nbin);
                        if(!(ni==0 && ui==0)){
                            chi2[i][j][k] = 1.0*(ui-ni)*(ui-ni)/(ni+ui);
                            totalchi2 += chi2[i][j][k];
                            totalchi2_indiv[j][k] += chi2[i][j][k];
                            totalchi2_PMTs[j] += chi2[i][j][k];
                            if( !needs_mask ){
                                totalchi2_selected  += chi2[i][j][k];
                            }
                        }
                        else{cout<<"get 0 .... get 0 ... bin "<<bin+1<<ui<<", "<<ni<<", "<<i<<", "<<j<<", "<<k<<endl;}
                        break;
                    }
                    /*
                    totalchi2 += chi2[i][j][k];
                    totalchi2_indiv[j][k] += chi2[i][j][k];
                    totalchi2_PMTs[j] += chi2[i][j][k];
                    if( !needs_mask ){
                       totalchi2_selected  += chi2[i][j][k];
                    }
                    */
                }
                // overflow bins to be counted
                //ui = hMCnpe[i][j][k]->Integral(nb_fit,50);
                //ni = hDatanpe[j][k]->Integral(nb_fit,50);
                //totalchi2  +=  1.0*(ui-ni)*(ui-ni)/(ui+ni);
                //totalchi2_indiv[j][k] += 1.0*(ui-ni)*(ui-ni)/(ui+ni);
                //totalchi2_PMTs[j] += 1.0*(ui-ni)*(ui-ni)/(ui+ni);
                //if( !needs_mask ){
                //   totalchi2_selected  += 1.0*(ui-ni)*(ui-ni)/(ui+ni);
                //}
            }
        }
        MCFile[i]->Close();
        //cout<<"calculation for mc point "<<i+1<<" is done."<<endl;
        int graph_index = 0;
        if(i<file_index_as_data) graph_index = i;
        if(i>file_index_as_data) graph_index = i-1;
        gChi2->SetPoint(graph_index,attenuationPar[i],totalchi2);
        gChi2_TrigPMTSelected->SetPoint(graph_index,attenuationPar[i],totalchi2_selected);
        for(int j=0;j<pmtnb;j++){
            gChi2_PMTs[j]->SetPoint(graph_index,attenuationPar[i],totalchi2_PMTs[j]);
            totalchi2_PMTs[j] = 0;
        }
        for(int k=0; k<nbTrigs; k++){
            for(int j=0; j<pmtnb; j++){
                totalchi2_TrigTypes[k] += totalchi2_indiv[j][k];
            }
            gChi2_TrigTypes[k]->SetPoint(graph_index,0.1*(i+1),totalchi2_TrigTypes[k]);
            totalchi2_TrigTypes[k]=0;
        }
        //fout<<i<<"\t"<<graph_index<<"\t"<<attenuationPar[i]<<"\t"<<totalchi2<<endl;
        totalchi2 = 0;
        totalchi2_selected=0;
        for(int j=0; j<pmtnb; j++){
            for(int k=0; k<nbTrigs; k++){
                gChi2_indiv[j][k]->SetPoint(graph_index,attenuationPar[i],totalchi2_indiv[j][k]);
                if(attenuationPar[i]==0)
                    cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<attenuationPar[i]<<"\t"<<totalchi2_indiv[j][k]<<endl;
                totalchi2_indiv[j][k]=0;
                //chi2[i][j][k] = 0;
            }
        }

    }

    string outfilename = "file_chi2_results_fitData_30m_scaleFactor.root";
    TFile* outfile = new TFile(outfilename.c_str(),"recreate");
    char cName[100];

    TCanvas* cChi2 = new TCanvas();
    cChi2->SetName("cChi2");
    cChi2->SetTitle("cChi2");
    gChi2->Draw("AP");
    gChi2->SetMarkerStyle(8);
    gChi2_TrigPMTSelected->Draw("P");
    gChi2_TrigPMTSelected->SetMarkerStyle(21);
    gChi2_TrigPMTSelected->SetMarkerColor(2);
    gChi2->GetXaxis()->SetTitle("Scale factor");
    gChi2->GetYaxis()->SetTitle("Chi2 value");
    cChi2->Write();
    cChi2->Close();

    TCanvas* cChi2_TrigTypes[nbTrigs];
    TCanvas* cChi2_TrigTypes_together = new TCanvas("cChi2TrigTypesTogether","cChi2TrigTypesTogether",800,600);
    TLegend* leg_TrigTypes_together = new TLegend(0.7,0.6,0.9,0.9);
    for(int k=0; k<nbTrigs; k++){
        cChi2_TrigTypes[k] = new TCanvas();
        sprintf(cName,"chi2_trig_%s",strig[k].c_str());
        cChi2_TrigTypes[k]->SetName(cName);
        cChi2_TrigTypes[k]->SetTitle(cName);
        gChi2_TrigTypes[k]->Draw("AP");
        gChi2_TrigTypes[k]->GetXaxis()->SetTitle("Scale factor");
        gChi2_TrigTypes[k]->GetYaxis()->SetTitle("chi2");
        gChi2_TrigTypes[k]->SetTitle(cName);
        gChi2_TrigTypes[k]->SetMarkerStyle(8);
        cChi2_TrigTypes[k]->Write();
        cChi2_TrigTypes[k]->Close();
        cChi2_TrigTypes_together->cd();
        if(k==0) gChi2_TrigTypes[k]->Draw("AP");
        else gChi2_TrigTypes[k]->Draw("P");
        gChi2_TrigTypes[k]->SetMarkerColor(1+k/3);
        gChi2_TrigTypes[k]->SetMarkerStyle(20+k%9);
        leg_TrigTypes_together->AddEntry(gChi2_TrigTypes[k],cName,"lp");
    }
    cChi2_TrigTypes_together->cd();
    leg_TrigTypes_together->Draw();
    cChi2_TrigTypes_together->Write();
    cChi2_TrigTypes_together->Close();

    TCanvas* cChi2_PMTs[pmtnb];
    TCanvas* cChi2_PMTs_together = new TCanvas("cChi2PMTsTogether","cChi2PMTsTogether",800,600);
    TLegend* leg_PMTs_together = new TLegend(0.7,0.6,0.9,0.9);
    for(int k=0; k<pmtnb; k++){
        cChi2_PMTs[k] = new TCanvas();
        sprintf(cName,"canvas_chi2_PMTs_S%d",k);
        cChi2_PMTs[k]->SetName(cName);
        cChi2_PMTs[k]->SetTitle(cName);
        //cout<<"test out "<<k<<endl;
        gChi2_PMTs[k]->Draw("AP");
        gChi2_PMTs[k]->GetXaxis()->SetTitle("Scale factor");
        gChi2_PMTs[k]->GetYaxis()->SetTitle("chi2");
        gChi2_PMTs[k]->SetTitle(cName);
        gChi2_PMTs[k]->SetMarkerStyle(8);
        cChi2_PMTs[k]->Write();
        cChi2_PMTs[k]->Close();
        cChi2_PMTs_together->cd();
        if(k==0) gChi2_PMTs[k]->Draw("AP");
        else gChi2_PMTs[k]->Draw("P");
        gChi2_PMTs[k]->SetMarkerColor(k+1);
        gChi2_PMTs[k]->SetMarkerStyle(20+k/2);
        leg_PMTs_together->AddEntry(gChi2_PMTs[k],cName,"lp");
    }
    cChi2_PMTs_together->cd();
    leg_PMTs_together->Draw();
    cChi2_PMTs_together->Write();
    cChi2_PMTs_together->Close();

    
    TCanvas* cIndiv[pmtnb][nbTrigs];

    TCanvas* cIndiv_1[nbTrigs];
    for(int i=0;i<nbTrigs;i++){
        sprintf(cName,"trig_%s",strig[i].c_str());
        cIndiv_1[i] = new TCanvas(cName,cName,2000,1000);
        //cIndiv_1[i]->SetName(cName);
        //cIndiv_1[i]->SetTitle(cName);
        cIndiv_1[i]->Divide(4,2);
    }
    
    
    for(int j=0; j<pmtnb; j++){
        for(int k=0; k<nbTrigs; k++){
            sprintf(cName,"pmt_S%d_trig_%s",j,strig[k].c_str());
            cIndiv[j][k] = new TCanvas();
            cIndiv[j][k]->SetName(cName);
            cIndiv[j][k]->SetTitle(cName);
            gChi2_indiv[j][k]->Draw("AP");
            gChi2_indiv[j][k]->SetTitle(cName);
            gChi2_indiv[j][k]->SetName(cName);
            gChi2_indiv[j][k]->SetMarkerStyle(8);
            gChi2_indiv[j][k]->GetXaxis()->SetTitle("Scale factor");
            gChi2_indiv[j][k]->GetYaxis()->SetTitle("Chi2 value");
            cIndiv[j][k]->Write();
            cIndiv[j][k]->Close();
            cIndiv_1[k]->cd(j+1);
            gChi2_indiv[j][k]->Draw("AP");
//            gChi2_indiv[j][k]->SetMarkerStyle(20+k);
//            gChi2_indiv[j][k]->SetMarkerColor(j+1);
        }
    }
    
    
    for(int i=0; i<nbTrigs; i++){
        cIndiv_1[i]->Write();
        sprintf(cName,"fit_data/30m_chi2_results_trigType_%s.png",strig[i].c_str());
        cIndiv_1[i]->SaveAs(cName);
        cIndiv_1[i]->Close();
    }
    

    TGraph* gChi2_TrigTypes_shift[nbTrigs];
    TGraph* gChi2_PMTs_shift[pmtnb];
    TCanvas* c_TrigTypes_shift = new TCanvas();
    c_TrigTypes_shift->SetName("c_TrigTypes_shift");
    c_TrigTypes_shift->SetTitle("c_TrigTypes_shift");
    TLegend* leg_trigTypes_shift = new TLegend(0.4,0.7,0.9,0.9);
    leg_trigTypes_shift->SetNColumns(2);
    for(int i=0; i<nbTrigs; i++){
        gChi2_TrigTypes_shift[i] = new TGraph();
        int N = gChi2_TrigTypes[i]->GetN();
        double x,y;
        double setx,sety;
        gChi2_TrigTypes[i]->GetPoint(17,setx,sety);
        //cout<<x<<"\t"<<y<<endl;
        for(int n=0;n<N;n++){
            gChi2_TrigTypes[i]->GetPoint(n,x,y);
            gChi2_TrigTypes_shift[i]->SetPoint(n,x,y-sety+100);
        }

        if(i==0) gChi2_TrigTypes_shift[i]->Draw("AP");
        else gChi2_TrigTypes_shift[i]->Draw("P");
        gChi2_TrigTypes_shift[i]->SetMarkerStyle(20+i%9);
        gChi2_TrigTypes_shift[i]->SetMarkerColor(1+i/3);
        sprintf(name,"Trig_%s",strig[i].c_str());
        leg_trigTypes_shift->AddEntry(gChi2_TrigTypes_shift[i],name,"p");
    }
    leg_trigTypes_shift->Draw();
    c_TrigTypes_shift->Write();

    TCanvas* c_PMTs_shift = new TCanvas();
    c_PMTs_shift->SetName("c_PMTs_shift");
    c_PMTs_shift->SetTitle("c_PMTs_shift");
    TLegend* leg_PMTs_shift = new TLegend(0.4,0.7,0.9,0.9);
    leg_PMTs_shift->SetNColumns(2);
    for(int i=0; i<pmtnb; i++){
        gChi2_PMTs_shift[i] = new TGraph();
        int N = gChi2_PMTs[i]->GetN();
        double x,y;
        double setx,sety;
        gChi2_PMTs[i]->GetPoint(17,setx,sety);
        //cout<<x<<"\t"<<y<<endl;
        for(int n=0;n<N;n++){
            gChi2_PMTs[i]->GetPoint(n,x,y);
            gChi2_PMTs_shift[i]->SetPoint(n,x,y-sety+100);
        }

        if(i==0) gChi2_PMTs_shift[i]->Draw("AP");
        else gChi2_PMTs_shift[i]->Draw("P");
        gChi2_PMTs_shift[i]->SetMarkerStyle(20+i%4);
        gChi2_PMTs_shift[i]->SetMarkerColor(1+i/2);
        sprintf(name,"PMT_S%d",i);
        leg_PMTs_shift->AddEntry(gChi2_PMTs_shift[i],name,"p");
    }
    leg_PMTs_shift->Draw();
    c_PMTs_shift->Write();
    
    outfile->Close();

    vector<double> fitPara;
    fitPara.clear();
    fstream fout_fit("test_fitData_parameters_20190725_30m_scaleFactor.txt",ios::out);
    fout_fit<<"True attenuation length is unknown."<<endl;
    //fit parabolics
    cout<<"\nfit total chi2 "<<endl;
    fitPara = fitParabolicAnalyticly(gChi2);
    fout_fit<<"use all PMTs all TrigTypes \t"<<fitPara[1]<<"\t"<<1.0/TMath::Sqrt(fitPara[0])<<"\t(minChi2 ="<<fitPara[2]<<")"<<endl;
    cout<<"\nfit PMTs (trigTypes all combined)"<<endl;
    fout_fit<<"\nuse all TrigTypes "<<endl;
    for(int i=0; i<pmtnb; i++){
        cout<<gChi2_PMTs[i]->GetName()<<endl;
        fitPara = fitParabolicAnalyticly(gChi2_PMTs[i]);
        fout_fit<<gChi2_PMTs[i]->GetName()<<"\t"<<fitPara[1]<<"\t"<<1.0/TMath::Sqrt(fitPara[0])<<"\t(minChi2 ="<<fitPara[2]<<")"<<endl;
    }
    cout<<"\nfit trigTypes (PMTs all combined)"<<endl;
    fout_fit<<"\nfit trigTypes (PMTs all combined)"<<endl;
    for(int i=0;i<nbTrigs;i++){
        cout<<gChi2_TrigTypes[i]->GetName()<<endl;
        fitPara = fitParabolicAnalyticly(gChi2_TrigTypes[i]);
        fout_fit<<gChi2_TrigTypes[i]->GetName()<<"\t"<<fitPara[1]<<"\t"<<1.0/TMath::Sqrt(fitPara[0])<<"\t(minChi2 ="<<fitPara[2]<<")"<<endl;
    }
    fout_fit.close();

//    TGraphErrors* gChi2_indiv[pmtnb][nbTrigs]; // chi2 for each PMT for each trigger type

    //AlignChi2Curves(outfilename);
    
}

//void AlignChi2Curves(string aname){
//    TFile* f = TFile(aname.c_str(),"update");
//    
//}

vector<double> fitParabolicAnalyticly(TGraph* gGraph){ // fit parabolic analyticly using three data points
  vector<double> ans; // this vector will hold the three parameters a, b, c in equation a(x-b)^2+c
  for(int i=0; i<3; i++){
    ans.push_back(0.0);
  }
  //get the smallest data points in the graph
  int N = gGraph->GetN();
  double x, y;
  std::vector< std::pair<double, double> > data;
  std::vector< std::pair<double, double> > data1;
  for(int i=0;i<N;i++){
      gGraph->GetPoint(i,x,y);
      //cout<<i<<"\t"<<apair.first<<"\t"<<apair.second<<endl;
      data.push_back( std::make_pair(y,x) );
      data1.push_back( std::make_pair(x,y) );
  }
  sort(data.begin(),data.end());
  sort(data1.begin(),data1.end());
  //take the first three data points for the parabolic
  double x1=data[0].second; // x2=data[1].second, x3=data[2].second;
  double y1=data[0].first; //  y2=data[1].first,  y3=data[2].first;
  double x2, y2, x3, y3;
  for(int i=0; i<N; i++){
      if(data1[i].first == x1){
          x2 = data1[i-1].first;
          y2 = data1[i-1].second;
          x3 = data1[i+1].first;
          y3 = data1[i+1].second;
          break;
      }
  }
  cout<<"fit to ("<<x1<<", "<<y1<<"), ("<<x2<<", "<<y2<<"), ("<<x3<<", "<<y3<<")\t";
  double b = ( (y1-y3)*(x1*x1-x2*x2)-(y1-y2)*(x1*x1-x3*x3) ) / ( (y1-y3)*(x1-x2)-(y1-y2)*(x1-x3) ) /2.0;
  double a = (y1-y2)/( (x1-b)*(x1-b) - (x2-b)*(x2-b) );
  double c = y1 - a*(x1-b)*(x1-b);
  cout<<"get a, b, c: "<<a<<"\t"<<b<<"\t"<<c<<endl;
  ans[0] = a;
  ans[1] = b;
  ans[2] = c;
  cout<<"fitted result: minimum at "<<b<<" m, (chi2="<<c<<")"<<endl;
  cout<<"               uncertainty "<< 1.0/TMath::Sqrt(a)<<endl;  
  return ans;
}
