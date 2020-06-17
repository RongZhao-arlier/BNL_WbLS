
#include <string>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TBranch.h>
#include <stdlib.h>  
using namespace std;
int main(){
std::string rootfilename="../../aw325/OneTon_postAnalysis/Water/processed_runs_20180416_20180612/run28638_Tree_v8.root";
//std::string rootfilename="ResultsV16_allTrigCombinations_run19885_run20085_20190605.root";
  TFile* rtfile = new TFile(rootfilename.c_str(),"read");
  TTree* npeTree = nullptr;
double pmtnpe[8][20];
int trigtype,trigtype_multi,trigtype_hodo,trigtype_led;
  rtfile->GetObject("OneTonEvent",npeTree);
  TBranch *npebranch  = npeTree->GetBranch("DigitizerPulseCharge");
  TBranch *tgbranch  = npeTree->GetBranch("TrigType");
  TBranch *tg_mbranch  = npeTree->GetBranch("TrigTypeFlag_Multi");
  TBranch *tg_hbranch  = npeTree->GetBranch("TrigTypeFlag_Hodo");
  TBranch *tg_lbranch  = npeTree->GetBranch("TrigTypeFlag_Led");
  npebranch->SetAddress(pmtnpe);
  tgbranch->SetAddress(&trigtype);
  tg_mbranch->SetAddress(&trigtype_multi);
  tg_hbranch->SetAddress(&trigtype_hodo);
  tg_lbranch->SetAddress(&trigtype_led);
  int nevent = npeTree->GetEntries();
    for (Int_t i=0;i<nevent;i++) {
      //     bntrac//k->GetEvent(i);
      npeTree->GetEntry(i);
      for(int j=0;j<8;j++)
	{ 
        if(pmtnpe[0][0]==0){for(int k=0;k<1;k++)cout<<pmtnpe[j][k]<<" ";}
          }
	if(pmtnpe[0][0]==0)cout<<"tgigtype:"<<trigtype<<"mlti:"<<trigtype_multi<<"hodo:"<<trigtype_hodo<<"led:"<<trigtype_led<<endl;
	if(i>100000)break;
	}
return 0;
}
