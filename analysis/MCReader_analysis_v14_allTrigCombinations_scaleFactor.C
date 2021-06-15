#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

#include <RAT/DS/Run.hh>
#include <RAT/DS/PMTInfo.hh>
#include <RAT/DS/Root.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TMath.h>
#include <TApplication.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

//this version is used to test my reconstruction idea with single (known) event.
double getratio(double thick,double att,double thetai,double ni,double nr);
double getangle(double thetai,double ni,double nr);

int main(int argc, char **argv) {

    const int numberofinputfiles = atoi(argv[2]);//150; //148 number of input files

    double npescalefactor[8]; // scale the nep, just like there is an overall efficiency on photon detection 
                              // every PMT may have different scale factor
 int hodoselect=atoi(argv[11]);
int d4x=atoi(argv[12]);
int d4y=atoi(argv[13]);
int d4z=atoi(argv[14]);
//	int hodoselect=11;  //0, for Hodo+H6, 1 for Hodo+H5, 2 for Hodo only, 3, for Hodo26, 4 for Hodo25
// 5 for Hodo 24,6for Hodo 21, 7for Hodo 18
       // 11 for 1010; 12 for 1001; 13 for 0110
       // 14 for 0101; 15 for 1110; 16 for 1101
       // 17 for 1101; 18for 1011; 19 for 1111
	TString hsstring[21]={"PMTnpeHodoh5","PMTnpeHodoh4","PMTnpeHodoonly","PMTnpeHodo26","PMTnpeHodo25","PMTnpeHodo24","PMTnpeHodo21","PMTnpeHodo18","08","09","10","PMTnpeHodotrig1010","PMTnpeHodotrig1001","PMTnpeHodotrig0110","PMTnpeHodotrig0101","PMTnpeHodotrig1110","PMTnpeHodotrig1101","PMTnpeHodotrig0111","PMTnpeHodotrigcheckh4","PMTnpeHodotrigcheckh5","PMTnpeHodotrig111101"};
    for(int i=0;i<8;i++){
        npescalefactor[i] = atof(argv[3+i]);
    }
    int nbadruns = 0;
    //this file is without time information
    string headdir = "";
//    string headdir = "./";
    string inputrootfile[numberofinputfiles];

    fstream finnames(argv[1],ios::in);
	
    string onename;
    int cntname = 0;
    while(finnames>>onename){
        inputrootfile[cntname] = onename;
        cntname++;
	if(cntname>=numberofinputfiles) break;
    }
    finnames.close();
//rong zhao define the pdf output
    char tempoutfilechar[100]; 
    sprintf(tempoutfilechar,"Results_%sh%s.txt", argv[1],argv[11]);
  fstream txtout(tempoutfilechar,ios::out);
  TString outpdf=("./pmtznpe.pdf");
  TString outpdf_start=("./pmtznpe.pdf[");
  TString outpdf_end=("./pmtznpe.pdf]");
  TCanvas *can1=new TCanvas("can1","canvas1",1200,900);
  can1->cd();
  can1->Print(outpdf_start);

    for(int i=0;i<numberofinputfiles;i++) inputrootfile[i] = headdir + inputrootfile[i];
	
    std::string sDigiChMap[8]={"Digi1_CH0_S0","Digi1_CH1_S1","Digi1_CH2_S2","Digi1_CH3_S3",
                               "Digi2_CH0_S4","Digi2_CH1_S5","Digi2_CH2_S6","Digi2_CH3_S7"};
    std::string sDigiCh[8] = {"S0","S1","S2","S3","S4","S5","S6","S7"};
    std::string sTDC_ch_Map[16]={"H0","H1","H2","H3","H4","H5",
                                 "HodoTrig","MultTrig","LedTrig","AllTrig",
                                 "S0","S1","S2","S3","S4","S5"}; // here the "AllTrig" means multi_or_hodo_or_led_

    const int nbTrigs = 27;
    // trigger types:
    // I have 27 types of combinations based on the normal hodotrig defined by (H2 or H0) and (H3 or H1)
    // Using the order H2, H0, H3, H1, H4, H5
    // and 0, 1 for hit or no hit, I list the combinations as following:
double   acrylicatt[90]={1 ,0.999924 ,0.999696 ,0.999316 ,0.998783 ,0.998097 ,0.997256 ,0.99626 ,0.995107 ,0.993796 , 0.992325 , 0.990692 , 0.988896 , 0.986932 , 0.9848 , 0.982495 , 0.980015 , 0.977356 , 0.974515 , 0.971487 , 0.968267 , 0.964852 , 0.961235 , 0.957411 , 0.953373 , 0.949116 , 0.944632 , 0.939914 , 0.934952 , 0.929738 , 0.924262 , 0.918514 , 0.912481 , 0.906152 , 0.899511 , 0.892545 , 0.885236 , 0.877566 , 0.869516 , 0.861061 , 0.852178 , 0.842839 , 0.833011 , 0.822659 , 0.811743 , 0.800215 , 0.788023 , 0.775103 , 0.761382 , 0.746773 , 0.73117 , 0.714446 , 0.696443 , 0.676966 , 0.655761 , 0.632502 , 0.606749 , 0.577901 , 0.545095 , 0.50705 , 0.461738 , 0.405708 , 0.33242 , 0.22741 , 0.0545496 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0};
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
     int indexMap_forEvtStartTime[nbTrigs]={2,2,2,2,2,2,  0,0,0,0,0,0,   2,2,2,2,2,2, 0,0,0,  2,2,2,2,2,2};
    // the above scheme, will select exclusive triggers. 
    // I ignore the case where both H4 and H5 are hit (because they are mostly caused by one muon and one scattered electron or gamma)
    // For the case where H4/H5 information should not be considered, the union of H4/H5 "01", "10", "00" will make it ("11" is ignored).

    char temprtfilechar[100]; sprintf(temprtfilechar,"Results_v13_allTrigCombinations_%s_scaleFactor_%s.root", argv[1],argv[3]);
    TFile* foutroot = new TFile(temprtfilechar,"recreate");

    // cut on plastic scintillators, 0.5 MeV
    double edepcut = 0.5;
    // PMT spe calibration
    double converterfactor[8]={146.278, 169.495, 180.849, 173.302, 165.905, 165.331, 140.525, 138.591};
    // the mean values of the pulse start time for muons, I will update the values using the arrival time of first photons in PMTs.    
    float meanPulseStartBin[8]={7.198, 7.098, 7.396, 7.904, 5.509, 4.151, 7.179, 7.026};
    double coincidenceTimeMeanMatrix[8][8]={
              {-100,       -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {-0.194713,   -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {-0.758341,  -0.547047,   -100,      -100,      -100,     -100,      -100,      -100},
              {-1.36467,   -1.17166,  -0.647578,  -100,      -100,     -100,      -100,      -100},
              {0.804999,   1.01294,   1.56821 ,  2.18819,     -100,     -100,      -100,      -100},
              {1.71103,   2.02743,    2.57424, 3.23104,   1.01547, -100,      -100,      -100},
              {3.09198,   -0.981527,   -0.581544, 0.166766,   -2.12371,  -3.31489, -100,      -100},
              {1.36372	, 0.18576,   -0.961526,  -0.460571,   -2.55673, -3.96895, -0.205749, -100}    
                                       };
    double coincidenceTimeWidthMatrix[8][8]={
              {-100,      -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.14204,   -100,       -100,      -100,      -100,     -100,      -100,      -100},
              {1.14403,  0.918625,     -100,      -100,      -100,     -100,      -100,      -100},
              {1.14654,   0.882617,    0.962183,   -100,      -100,     -100,      -100,      -100},
              {1.13599,  0.906206,     1.10838,   1.14001,   -100,     -100,      -100,      -100},
              {1.37561,  1.36813,     1.49745,    1.6342,  1.78241,  -100,      -100,      -100},
              {1.02948,  0.985831,      1.00617,   0.900906,   0.923113,  0.902838,   -100,      -100},
              {2.94202,  3.40961,      1.38917,   1.47374,   1.51825,   1.60992,   0.843143,   -100}    
                                       };
    // directories - by trigger type
    foutroot->cd();
    TDirectory* dir[nbTrigs];
    char dirnames[100];
    for(int i=0; i<nbTrigs; i++){
        sprintf(dirnames,"dir_%s",strig[i].c_str());
        dir[i] = foutroot->mkdir(dirnames);
        dir[i]->cd();
    }
    TDirectory* dir_mcwaveforms = foutroot->mkdir("MCWaveforms");
    const int pmtnb = 8;
    char hName[200];
    char tempname[150];
//rong zhao. save the npe info of PMTs into trees, for further Tminuit analysis 
 //  TFile *outrootnpe =new TFile(hsstring[hodoselect]+".root","recreate");
   TFile *outrootnpe =new TFile(hsstring[hodoselect]+argv[1]+".root","recreate");
//   TFile *outrootnpe =new TFile("PMTnpe.root","recreate");
outrootnpe->cd();
  TTree *npetree=new TTree("Treenpe","PMT npe tree");
  TTree *oltree=new TTree("Treeol","PMT ol tree");
  TTree *npedecaytree=new TTree("Treenpedecay","PMT npe decay tree");
int thistriggerflag=0;
double pmtnpe[pmtnb];
double pmtdis[pmtnb]={0};
double pmtol[pmtnb]={0};
double pmtia[pmtnb]={0};   //incident angle to pmt
double sinphi=0;   // muon track angle
double pmtatt[pmtnb]={0};   //att
double pmtattw[pmtnb]={0};   //att  water
double pmtattao[pmtnb]={0};   //att acrylic and cookie
double pmtattacry[pmtnb]={0};   //att acrylic and cookie
double pmtattcookie[pmtnb]={0};   //att acrylic and cookie
double pmtattsac[pmtnb]={0};   //att acrylic and cookie and solid angle
double pmtatts[pmtnb]={0};   //att  solid angle
double  crlx[8]={0};
double  crly[8]={0};
double  crlz[8]={0};
	    float zratio1[8]={0};
	    float hx1[8]={0};
	    float hy1[8]={0};
	    float zratio2[8]={0}; //hodo 1,3
	    float hx2[8]={0};
	    float hy2[8]={0};
	    float mustopz=0;
	    float mumumentum=0;
	    float muangle=0;
double pmtdecaynpe[pmtnb];
double pmtnpe_all;//rong zhao, to count the total npe for each event
double pmtnpe_alldecay;//rong zhao, to count the total npe for each event
double pmtposx[8]={381.35,131.35,-118.65,-368.65,-165.025,-165.025,381.35,131.35};
double pmtposy[8]={-30.79,-30.79,-30.79,-30.79,280.375,-128.425,169.21,169.21};
double pmtposz[8]={-501.4,-501.4,-501.4,-501.4,801.4,801.4,-501.4,-501.4};
float TrackStartX=0, TrackStartY=0, TrackStartZ=0;
npetree->Branch("pmtnpes", &pmtnpe, "pmtnpe[8]/D");
npetree->Branch("thistriggerglag", &thistriggerflag, "thistriggerflag/I");
npetree->Branch("brsmutopz", &mustopz, "mustopz/F");
npetree->Branch("brmumumentum", &mumumentum, "mumumentum/F");
npetree->Branch("brmuangle", &muangle, "muangle/F");
oltree->Branch("brpmtol", &pmtol, "pmtol[8]/D");
oltree->Branch("brpmtdis", &pmtdis, "pmtdis[8]/D");
oltree->Branch("brpmtia", &pmtia, "pmtia[8]/D");
oltree->Branch("brpmtatt", &pmtatt, "pmtatt[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtattw", &pmtattw, "pmtattw[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtattao", &pmtattao, "pmtattao[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtattacry", &pmtattacry, "pmtattacry[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtattcookie", &pmtattcookie, "pmtattcookie[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtattsac", &pmtattsac, "pmtattsac[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtatts", &pmtatts, "pmtatts[8]/D"); //the total attenuation eefect estimation for understanding MC, rongz hao
oltree->Branch("brpmtolx", &crlx, "crlx[8]/D");
oltree->Branch("brpmtoly", &crly, "crly[8]/D");
oltree->Branch("brpmtolz", &crlz, "crlz[8]/D");
oltree->Branch("brpmtstartx", &TrackStartX, "TrackStartX/F");
oltree->Branch("brpmtstarty", &TrackStartY, "TrackStartY/F");
oltree->Branch("brpmtstartz", &TrackStartZ, "TrackStartZ/F");
oltree->Branch("brhx1", &hx1, "hx1[8]/F");
oltree->Branch("brhy1", &hy1, "hy1[8]/F");
oltree->Branch("brhx2", &hx2, "hx2[8]/F");
oltree->Branch("brhy2", &hy2, "hy2[8]/F");
npedecaytree->Branch("pmtdecaynpes", &pmtdecaynpe, "pmtdecaynpe[8]/D");
    foutroot->cd();
    //input muons's momentum and angles
    TH1F* hMuStopz_all;
    TH1F* hMuStopz_coincidence;
    TH1F* hMuStopz_topcoincidence;
    TH1F* hMuStopz_botcoincidence;
    TH1F* hMuStopz_tbcoincidence;
    TH1F* hMuMomentum[nbTrigs];
    TH1F* hMuMomentum_all;
    TH1F* hMuMomentum_decay;
    TH1F* hMuMomentum_hodotrig;
    TH1F* hMuMomentum_coincidence;
    TH1F* hMuPAngle[nbTrigs];
    TH1F* hMuAAngle[nbTrigs];
    TH1F* hMuAAngletrigger;
    TH1F* hMuPAngletrigger;
    TH2F* hhodoXYtrigger[6]; // end x and y coordinates at  6 PS
    TH1F* hDecayZ[nbTrigs];
    TH2F* hEndXY[nbTrigs]; // end x and y coordinates at Z=H4/H5 
    TH2F* hEndXY0[nbTrigs];// end x and y coordinates when entering liquid (under the top acrylic lid)
    TH2F* hz_npe[nbTrigs];// rong  zhao 
    TH2F* hz_npe_all;// 
    TH2F* hz_mup_all;//
    TH2F* hz_npe_alldecay;// 
    TH2F* hmu_npe_alldecay;// 
    TH1D* hmu_ndecay;// 
    TH1D* hz_npe_all1d;// 
    TH2F* hz_npe_allemu;// 
    TH2F* hrp_npe[nbTrigs];//  the r prime is the distance to hodo_xy center
    TH2F* hrp_npe_all;// 
    TGraph* gNumberOfEvents = new TGraph();
    //TH2F* hEndXY1[nbTrigs];// end x and y coordinates when leaving liquid (under the bottom acrylic lid)
    int hindex[6]={2,0,3,1,4,5};
    for(int j=0;j<6;j++){
        sprintf(hName,"hXY_hodoH%d",hindex[j]); // at z corrdinate where H4 and H5 sit
	if(j<4){
        hhodoXYtrigger[j] = new TH2F(hName,"",200, 100, 500, 200, -200, 200);}
	else{
        hhodoXYtrigger[j] = new TH2F(hName,"",200, -500, 500, 200, -600, 350);}
        hhodoXYtrigger[j]->SetXTitle("X (mm)");
        hhodoXYtrigger[j]->SetYTitle("Y (mm)");
	}
        hMuStopz_all = new TH1F("hMuStopz_all","",60,-600,1000);
        hMuStopz_all->SetXTitle("Stopz position (mm)");
        hMuStopz_all->SetYTitle("Counts");
        hMuStopz_coincidence = new TH1F("hMuStopz_coincidence","",60,-600,1000);
        hMuStopz_coincidence->SetXTitle("Stopz position (mm)");
        hMuStopz_coincidence->SetYTitle("Counts");
        hMuStopz_topcoincidence = new TH1F("hMuStopz_topcoincidence","",60,-600,1000);
        hMuStopz_topcoincidence->SetXTitle("Stopz position (mm)");
        hMuStopz_topcoincidence->SetYTitle("Counts");
        hMuStopz_botcoincidence = new TH1F("hMuStopz_botcoincidence","",60,-600,1000);
        hMuStopz_botcoincidence->SetXTitle("Stopz position (mm)");
        hMuStopz_botcoincidence->SetYTitle("Counts");
        hMuStopz_tbcoincidence = new TH1F("hMuStopz_tbcoincidence","",60,-600,1000);
        hMuStopz_tbcoincidence->SetXTitle("Stopz position (mm)");
        hMuStopz_tbcoincidence->SetYTitle("Counts");

        hMuMomentum_all = new TH1F("hMuMomentum_all","",100,0,5000);
        hMuMomentum_all->SetXTitle("Momentum (MeV/c)");
        hMuMomentum_all->SetYTitle("Counts");
        hMuMomentum_decay = new TH1F("hMuMomentum_decay","",100,0,1000);
        hMuMomentum_decay->SetXTitle("Momentum (MeV/c)");
        hMuMomentum_decay->SetYTitle("Counts");
        hMuMomentum_hodotrig = new TH1F("hMuMomentum_hodotrig","",100,0,1000);
        hMuMomentum_hodotrig->SetXTitle("Momentum (MeV/c)");
        hMuMomentum_hodotrig->SetYTitle("Counts");
        hMuMomentum_coincidence = new TH1F("hMuMomentum_coincidence","",100,0,5000);
        hMuMomentum_coincidence->SetXTitle("Momentum (MeV/c)");
        hMuMomentum_coincidence->SetYTitle("Counts");
    for(int j=0;j<nbTrigs;j++){
        sprintf(hName,"hMuMomentum_%sTrig",strig[j].c_str());
        hMuMomentum[j] = new TH1F(hName,"",100,0,5000);
        hMuMomentum[j]->SetXTitle("Momentum (MeV/c)");
        hMuMomentum[j]->SetYTitle("Counts");
        sprintf(hName,"hMuPolarAngle_%sTrig",strig[j].c_str());
        hMuPAngle[j] = new TH1F(hName,"",170,0.5*TMath::Pi(),1.2*TMath::Pi());
        hMuPAngle[j]->SetXTitle("Polar angle (rad.)");
        hMuPAngle[j]->SetYTitle("Counts");
        sprintf(hName,"hMuAzimuthalAngle_%sTrig",strig[j].c_str());
        hMuAAngle[j] = new TH1F(hName,"",100, -1.0*TMath::Pi(), 1.0*TMath::Pi());
        hMuAAngle[j]->SetXTitle("Azimuthal angle (rad.)");
        hMuAAngle[j]->SetYTitle("Counts");
        hMuAAngletrigger = new TH1F("muon_angle","",100, 0, 360);
//        hMuAAngletrigger = new TH1F("muon_angle","",100, -1.0*TMath::Pi(), 1.0*TMath::Pi());
        hMuAAngletrigger->SetXTitle("Azimuthal angle (degree)");
        hMuAAngletrigger->SetYTitle("Counts");
        hMuPAngletrigger = new TH1F("muon_polar_angle","",100, 0, 180);
        hMuPAngletrigger->SetXTitle("zenith  angle (degree)");
        hMuPAngletrigger->SetYTitle("Counts");

        sprintf(hName,"hDecayZ_%s",strig[j].c_str());
        hDecayZ[j] = new TH1F(hName,"",240, -1200, 1200);
        hDecayZ[j]->SetXTitle("Z coordinate (mm)");
        hDecayZ[j]->SetYTitle("Counts");
        sprintf(hName,"hEndXY_%s_zAtH4H5",strig[j].c_str()); // at z corrdinate where H4 and H5 sit
        hEndXY[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        hEndXY[j]->SetXTitle("X (mm)");
        hEndXY[j]->SetYTitle("Y (mm)");
        sprintf(hName,"hEndXY_%s_z-500",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hEndXY0[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        hEndXY0[j]->SetXTitle("X (mm)");
        hEndXY0[j]->SetYTitle("Y (mm)");
        sprintf(hName,"hz_npe%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_npe[j] = new TH2F(hName,"",50, 0, 50, 100, -1000, 1000);
        hz_npe[j]->SetXTitle("npe ");
        hz_npe[j]->SetYTitle("Z (mm)");
        sprintf(hName,"hz_npe_all%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_npe_all = new TH2F(hName,"",50, -600, 1000, 100, 0, 50);
        hz_npe_all->SetYTitle("npe ");
        hz_npe_all->SetXTitle("Z (mm)");
        sprintf(hName,"hz_mup_all%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_mup_all = new TH2F(hName,"",50,-600, 1000, 50,0,1000);
        hz_mup_all->SetXTitle("stop Z (mm) ");
        hz_mup_all->SetYTitle("momentum (MeV)");
        sprintf(hName,"hz_npe_all1d%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_npe_all1d = new TH1D(hName,"",50, -600, 1000 );
        hz_npe_all1d->SetYTitle("npe ");
        hz_npe_all1d->SetXTitle("Z (mm)");
        sprintf(hName,"hz_npe_alldecay%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_npe_alldecay = new TH2F(hName,"",100, 0, 50, 100, -1000, 1000);
        hz_npe_alldecay->SetXTitle("npe decay");
        hz_npe_alldecay->SetYTitle("Z (mm)");
        sprintf(hName,"hmu_ndecay%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hmu_ndecay = new TH1D(hName,"",100, 0, 50);
        hmu_ndecay->SetYTitle("n decay");
        hmu_ndecay->SetXTitle("npe mu ");
        sprintf(hName,"hmu_npe_alldecay%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hmu_npe_alldecay = new TH2F(hName,"",100, 0, 50, 50, 0, 20);
        hmu_npe_alldecay->SetYTitle("npe decay");
        hmu_npe_alldecay->SetXTitle("npe mu ");
        sprintf(hName,"hz_npe_allemu%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hz_npe_allemu = new TH2F(hName,"",100, 0, 50, 50, 0, 20);
        hz_npe_allemu->SetYTitle("npe decay");
        hz_npe_allemu->SetXTitle("npe mu ");
        sprintf(hName,"hrp_npe%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hrp_npe[j] = new TH2F(hName,"",50, 0, 50, 100, 0, 1500);
        hrp_npe[j]->SetXTitle("npe ");
        hrp_npe[j]->SetYTitle("rp (mm)");
        sprintf(hName,"hrp_npe_all%s",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        hrp_npe_all = new TH2F(hName,"",50, 0, 50, 100,0, 1500);
        hrp_npe_all->SetXTitle("npe ");
        hrp_npe_all->SetYTitle("rp (mm)");
        //sprintf(hName,"hEndXY_%s_z-500",strig[j].c_str()); // at z corrdinate where muon enters the liquid (z=800.4)
        //hEndXY1[j] = new TH2F(hName,"",100, -1000, 1000, 100, -1000, 1000);
        //hEndXY1[j]->SetXTitle("X (mm)");
        //hEndXY1[j]->SetYTitle("Y (mm)");
    }

    // edep and hit time in the 6 plastic scintillators
    TH1F* hEdepInPS[6][nbTrigs];
    TH1F* hHitTimeInPS[6][nbTrigs];
    for(int i=0; i<6; i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(hName,"hEdepInPs_%d_%sTrig",i,strig[j].c_str());
            hEdepInPS[i][j]=new TH1F(hName,"",500,0,10);
            hEdepInPS[i][j]->SetXTitle("Energy (MeV)");
            hEdepInPS[i][j]->SetYTitle("Counts");
            sprintf(hName,"hHitTimeInPs_%d_%sTrig",i,strig[j].c_str());
            hHitTimeInPS[i][j]=new TH1F(hName,"",500,0,50);
            hHitTimeInPS[i][j]->SetXTitle("Time (ns)");
            hHitTimeInPS[i][j]->SetYTitle("Counts");
        }
    }
    // the event time using the hodoscope detector's time.
    TH1F* hEvtTime[nbTrigs]; 
    TH1F* hTrigTime[nbTrigs];
    TH1F* hDecayTime[nbTrigs];    //rong zhao 
    TH1F* hDecayTime_all;    //rong zhao 
        sprintf(tempname,"hEvtDecayTimeall");
        hDecayTime_all = new TH1F(tempname,"",100,0,3000);
        hDecayTime_all->SetXTitle("Event Decay time (ns)");
        hDecayTime_all->SetYTitle("Number of events");
    TH1F* hcoincidence_time;    //rong zhao 
        sprintf(tempname,"hcoincidence_time");
        hcoincidence_time = new TH1F(tempname,"",100,0,3000);
        hcoincidence_time->SetXTitle("coincidence_time (ns)");
        hcoincidence_time->SetYTitle("Number of events");
    TH1F* hcoincidence_npe[pmtnb];    //rong zhao 
    for(int j=0;j<pmtnb;j++){
        sprintf(tempname,"hcoincidence_npe_PMTS%d",j);
        hcoincidence_npe[j] = new TH1F(tempname,"",80,0,20);
        hcoincidence_npe[j]->SetXTitle("NPE ofcoincidence_events");
        hcoincidence_npe[j]->SetYTitle("Number of events");
    }
    for(int j=0;j<nbTrigs;j++){
        sprintf(tempname,"hEvtTime_PS_%s",strig[j].c_str());
        hEvtTime[j] = new TH1F(tempname,"",500,0,10);
        hEvtTime[j]->SetXTitle("Event begin time (ns)");
        hEvtTime[j]->SetYTitle("Number of events");
        sprintf(tempname,"hEvtTriggerTime_%s",strig[j].c_str());
        hTrigTime[j] = new TH1F(tempname,"",500,0,10);
        hTrigTime[j]->SetXTitle("Event begin time (ns)");
        hTrigTime[j]->SetYTitle("Number of events");
        sprintf(tempname,"hEvtDecayTime_%s",strig[j].c_str());
        hDecayTime[j] = new TH1F(tempname,"",500,0,2500);
        hDecayTime[j]->SetXTitle("Event Decay time (ns)");
        hDecayTime[j]->SetYTitle("Number of events");
    }
	
    //hit time differences in hodoscope detectors
    TH1F* hTDC_H0H2 = new TH1F("hTDC_H0_H2","TDC_H0 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H1H2 = new TH1F("hTDC_H1_H2","TDC_H1 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H3H2 = new TH1F("hTDC_H3_H2","TDC_H3 - TDC_H2", 1000,-10,10);
    TH1F* hTDC_H1H0 = new TH1F("hTDC_H1_H0","TDC_H1 - TDC_H0", 1000,-10,10);
    TH1F* hTDC_H3H0 = new TH1F("hTDC_H3_H0","TDC_H3 - TDC_H0", 1000,-10,10);
    TH1F* hTDC_H5H2 = new TH1F("hTDC_H5_H2","TDC_H5 - TDC_H2", 1000,00,20);
    TH1F* hTDC_H4H2 = new TH1F("hTDC_H4_H2","TDC_H4 - TDC_H2", 1000,00,20);
    TH1F* hTDC_H5H0 = new TH1F("hTDC_H5_H0","TDC_H5 - TDC_H0", 1000,00,20);
    TH1F* hTDC_H4H0 = new TH1F("hTDC_H4_H0","TDC_H4 - TDC_H0", 1000,00,20);
    
    TH1F* hTDC_H4H5 = new TH1F("hTDC_H4_H5","TDC_H4 - TDC_H5", 400,-20,20);
    TH2F* hEdep_vs_TDC_H1H2 = new TH2F("hEdep_vs_TDC_H1H2","",1000,-10,10,500,0,10);
    TH2F* hEdepH4_vs_EdepH5 = new TH2F("hEdepH4_vs_EdepH5","",50,0,10,50,0,10);
	
    //charge, time in PMTs one by one
    TH1F* hCharge[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 2.5us (from event start time, 0 ns, in the MC: muon starts at Z=1200 mm).
    TH1F* hCharge1[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 50ns. This might be useful
    TH1F* hCharge2[pmtnb][nbTrigs]; // the total charge in every pmt within a time window of 40ns. This might be useful
    TH1F* hTime[pmtnb][nbTrigs];  // time of arrival of the first photon at a pmt
    TH2F* hPMTNpeVsTime[pmtnb][nbTrigs]; // 2D plot of time of photon vs. charge
    TH2F* hPMTPhotonTimeSpannpe[pmtnb][nbTrigs]; //  time of all photons in each pmt
    TH1F* hPMTPhotonTimeSpan[pmtnb][nbTrigs]; //  time of all photons in each pmt
    for(int i=0;i<pmtnb;i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(hName,"hCharge_2500ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge[i][j] = new TH1F(hName,"",50,0,50);
            hCharge[i][j]->SetXTitle("Number of detected photons");
            hCharge[i][j]->SetYTitle("Counts");
            sprintf(hName,"hCharge_50ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge1[i][j] = new TH1F(hName,"",50,0,50);
            hCharge1[i][j]->SetXTitle("npe in 50 ns window after the first photon");
            hCharge1[i][j]->SetYTitle("Counts");
            sprintf(hName,"hCharge_40ns_S%d_%sTrig",i,strig[j].c_str());
            hCharge2[i][j] = new TH1F(hName,"",50,0,50);
            hCharge2[i][j]->SetXTitle("npe in 40 ns window after the first photon");
            hCharge2[i][j]->SetYTitle("Counts");
            sprintf(hName,"hTime_S%d_%sTrig",i,strig[j].c_str());
            hTime[i][j]=new TH1F(hName,"",500,0,50); // here I only lock myself with a 2560 ns window
            hTime[i][j]->SetXTitle("Time of the first photon arriving at PMT (ns)");
            hTime[i][j]->SetYTitle("Counts");
            sprintf(hName,"hPMTNpeVsFirstPhotonTime_S%d_%sTrig",i,strig[j].c_str());
            hPMTNpeVsTime[i][j]=new TH2F(hName,"",520,-2,50,2560,-0.5,2560-0.5);
            hPMTNpeVsTime[i][j]->SetXTitle("npe");
            hPMTNpeVsTime[i][j]->SetYTitle("Time of photons arriving at PMT (ns)");	
            sprintf(hName,"hPMTAllPhotonTime_S%d_%sTrig",i,strig[j].c_str());
            hPMTPhotonTimeSpan[i][j]=new TH1F(hName,"",500,0,5000);
            hPMTPhotonTimeSpan[i][j]->SetXTitle("Time difference between last and the first photons (ns)");
            hPMTPhotonTimeSpan[i][j]->SetYTitle("Counts");			
            sprintf(hName,"hPMTAllPhotonTimenpe_S%d_%sTrig",i,strig[j].c_str());
            hPMTPhotonTimeSpannpe[i][j]=new TH2F(hName,"",500,0,2500,50,0,50);
            hPMTPhotonTimeSpannpe[i][j]->SetXTitle("time of all  photons (ns)");
            hPMTPhotonTimeSpannpe[i][j]->SetYTitle("NPE");			
	}
    }
    // total charge from all PMTs
    TH1F* hTotalCharge[nbTrigs]; // total charge of all pmts in npe
    TH1F* hTotalChargeTop[nbTrigs]; // total charge of top pmts in npe
    TH2F* hTotalChargeratio_stopz[nbTrigs]; 
    TH1F* hTotalChargeBottom[nbTrigs]; // total charge of bottom pmts in npe
    TH1F* hNpeRatio[nbTrigs]; // charge ratio: total charge in top PMTs over total charge in ALL PMTs, 
	                          // weighted by number of PMTs.
    for(int j=0;j<nbTrigs;j++){
        sprintf(hName,"hTotalCharge_%sTrig",strig[j].c_str());
        hTotalCharge[j] = new TH1F(hName,"",1000,0,1000);
        hTotalCharge[j]->SetXTitle("Number of detected photons");
        hTotalCharge[j]->SetYTitle("Counts");
        sprintf(hName,"hTotalChargeratio_stopz_%sTrig",strig[j].c_str());
        hTotalChargeratio_stopz[j] = new TH2F(hName,"",100,0,1,160,-600,1000);
        hTotalChargeratio_stopz[j]->SetXTitle("NPE top/bot");
        hTotalChargeratio_stopz[j]->SetYTitle("Stopz[mm]");
        sprintf(hName,"hTotalChargeTop_%sTrig",strig[j].c_str());
        hTotalChargeTop[j] = new TH1F(hName,"",100,0,40);
        hTotalChargeTop[j]->SetXTitle("Number of detected photons");
        hTotalChargeTop[j]->SetYTitle("Counts");
        sprintf(hName,"hTotalChargeBottom_%sTrig",strig[j].c_str());
        hTotalChargeBottom[j] = new TH1F(hName,"",100,0,40);
        hTotalChargeBottom[j]->SetXTitle("Number of detected photons");
        hTotalChargeBottom[j]->SetYTitle("Counts");
	sprintf(hName,"hNpeRatio_TopBotPMTs_%sTrig",strig[j].c_str());
        hNpeRatio[j] = new TH1F(hName,"",100,0,1);
        hNpeRatio[j]->SetXTitle("Charge ratio: charge in top PMTs over charge in all PMTs");
        hNpeRatio[j]->SetTitle("Weighted by number of PMTs");
    }
	
    TH2F* hPulseTimeToS0[7][nbTrigs]; // the difference in pulse times between other PMTs and S0
    TH2F* hPulseTimeToS1[6][nbTrigs]; // the difference in pulse times between other PMTs and S1
    TH2F* hPulseTimeToS2[5][nbTrigs]; // the difference in pulse times between other PMTs and S2
    TH2F* hPulseTimeToS3[4][nbTrigs]; // the difference in pulse times between other PMTs and S3
    TH2F* hPulseTimeToS4[3][nbTrigs]; // the difference in pulse times between other PMTs and S4
    TH2F* hPulseTimeToS5[2][nbTrigs]; // the difference in pulse times between other PMTs and S5
    TH2F* hPulseTimeToS6[1][nbTrigs]; // the difference in pulse times between other PMTs and S6
    TH1F* hPulseTimeDiffToS0[7][nbTrigs]; // the difference in pulse times between other PMTs and S0
    TH1F* hPulseTimeDiffToS1[6][nbTrigs]; // the difference in pulse times between other PMTs and S1
    TH1F* hPulseTimeDiffToS2[5][nbTrigs]; // the difference in pulse times between other PMTs and S2
    TH1F* hPulseTimeDiffToS3[4][nbTrigs]; // the difference in pulse times between other PMTs and S3
    TH1F* hPulseTimeDiffToS4[3][nbTrigs]; // the difference in pulse times between other PMTs and S4
    TH1F* hPulseTimeDiffToS5[2][nbTrigs]; // the difference in pulse times between other PMTs and S5
    TH1F* hPulseTimeDiffToS6[1][nbTrigs]; // the difference in pulse times between other PMTs and S6
    for(int i=0; i<pmtnb; i++){
        for(int j=0;j<nbTrigs;j++){
            sprintf(tempname,"hPulseTimeDiffToS0_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
            hPulseTimeDiffToS0[i-1][j] = new TH1F(tempname,"",2560*4,-2560,2560);
            sprintf(tempname,"Pulse time differencce between %s and S0", sDigiChMap[i].c_str());
            hPulseTimeDiffToS0[i-1][j]->SetXTitle(tempname);
            hPulseTimeDiffToS0[i-1][j]->SetYTitle("Counts / 1 ns");
            sprintf(tempname,"hPulseTimeToS0_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
            hPulseTimeToS0[i-1][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
            hPulseTimeToS0[i-1][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
            hPulseTimeToS0[i-1][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",0));
            if(i>1){
              sprintf(tempname,"hPulseTimeDiffToS1_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS1[i-2][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S1", sDigiChMap[i].c_str());
              hPulseTimeDiffToS1[i-2][j]->SetXTitle(tempname);
              hPulseTimeDiffToS1[i-2][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS1_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS1[i-2][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS1[i-2][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS1[i-2][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",1));
            }
            if(i>2){
              sprintf(tempname,"hPulseTimeDiffToS2_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS2[i-3][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S2", sDigiChMap[i].c_str());
              hPulseTimeDiffToS2[i-3][j]->SetXTitle(tempname);
              hPulseTimeDiffToS2[i-3][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS2_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS2[i-3][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS2[i-3][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS2[i-3][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",2));
            }
            if(i>3){
              sprintf(tempname,"hPulseTimeDiffToS3_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS3[i-4][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S3", sDigiChMap[i].c_str());
              hPulseTimeDiffToS3[i-4][j]->SetXTitle(tempname);
              hPulseTimeDiffToS3[i-4][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS3_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS3[i-4][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS3[i-4][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS3[i-4][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",3));
            }
            if(i>4){
              sprintf(tempname,"hPulseTimeDiffToS4_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS4[i-5][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S4", sDigiChMap[i].c_str());
              hPulseTimeDiffToS4[i-5][j]->SetXTitle(tempname);
              hPulseTimeDiffToS4[i-5][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS4_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS4[i-5][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS4[i-5][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS4[i-5][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",4));
            }
            if(i>5){
              sprintf(tempname,"hPulseTimeDiffToS5_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS5[i-6][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S5", sDigiChMap[i].c_str());
              hPulseTimeDiffToS5[i-6][j]->SetXTitle(tempname);
              hPulseTimeDiffToS5[i-6][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS5_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS5[i-6][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS5[i-6][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS5[i-6][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",5));
            }
            if(i>6){
              sprintf(tempname,"hPulseTimeDiffToS6_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeDiffToS6[i-7][j] = new TH1F(tempname,"",2560*4,-2560,2560);
              sprintf(tempname,"Pulse time differencce between %s and S6", sDigiChMap[i].c_str());
              hPulseTimeDiffToS6[i-7][j]->SetXTitle(tempname);
              hPulseTimeDiffToS6[i-7][j]->SetYTitle("Counts / 1 ns");
              sprintf(tempname,"hPulseTimeToS6_%s_%s",sDigiChMap[i].c_str(),strig[j].c_str());
              hPulseTimeToS6[i-7][j] = new TH2F(tempname,"",256,0,2560,256,0,2560);
              hPulseTimeToS6[i-7][j]->SetXTitle(Form("Pulse time of PMTS%d [ns]",i));
              hPulseTimeToS6[i-7][j]->SetYTitle(Form("Pulse time of PMTS%d [ns]",6));
            }
        }
    } // end preparing histograms

    //some variables
    int test_cnt1 = 0, test_cnt2 = 0;
    bool istrigger[nbTrigs];
    int  trigcnts[nbTrigs];
    int trigcnts_hodo = 0;
    int trigcnts_hodo_H5_notH4=0;
    int trigcnts_hodo_notH5_notH4=0;
    int trigcnts_hodo_notH5_H4=0;
    int trigcnts_hodo_H5_H4=0;
    int trigcnts_multi=0;
    int trigcnts_decays = 0;
    int total_number_of_events = 0;
    double eventStartTime[nbTrigs];
double hz1=1133.5; //the z position of H2,H4 
double hz2=896.5; //the z position of H2,H4 
    for(int index=0;index<nbTrigs;index++){
      istrigger[index]=false;
      eventStartTime[index]=0;
      trigcnts[index]=0;
    }
    int testcntTtrigtype1[8]; for(int t=0;t<8;t++) testcntTtrigtype1[t]=0;

    //Getting ready to look through ROOT Tree files
    TFile* f;
    TTree* T;
    TTree* runT;
    for(int ffff=0; ffff<numberofinputfiles; ffff++){ // start loop: reading the root tree files
        //open the root file
        f = new TFile(inputrootfile[ffff].c_str(),"read");
        std::cout<<"Reading "<<inputrootfile[ffff]<<std::endl;
        T = (TTree*) f->Get("T");
        if(T==NULL){
            cout<<"bad tree "<<inputrootfile[ffff]<<endl;
            continue;
        }
        runT = (TTree*) f->Get("runT");
        RAT::DS::Run *run = new RAT::DS::Run();
        runT->SetBranchAddress("run", &run);
        if (runT->GetEntries() != 1) {
          cout << "Funny run tree, aborting" << endl;
          return -1;
        }
        runT->GetEntry(0);
        RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
        RAT::DS::Root *ds = new RAT::DS::Root();
        T->SetBranchAddress("ds", &ds);
        int nEvents = T->GetEntries();
        cout << "Reading in " << nEvents << " events" << endl;

        //reading in events
        for (int i = 0; i < nEvents; i++) {  
            total_number_of_events++;
            //if(i%100==0) cout<<"Event "<<i<<endl;
            T->GetEntry(i);
            if (ds->GetEVCount() != 1) {
                cout << "EV " << i << " is multi-trigger, ignoring" << endl;
                continue;
            }
        
            RAT::DS::MC *mc = ds->GetMC();
            RAT::DS::EV *ev = ds->GetEV(0);

            std::vector<int> myTriggerFlag = ev->GetTriggerFlag();
            std::vector<double> myTriggerEdep = ev->GetTriggerEdep();
            std::vector<double> myTriggerHitTime = ev->GetTriggerHitTime();

            double t_hodo_1 = 0, t_hodo_2=0;
            // determine trigger time in H0 and H2 
            if(myTriggerHitTime[0]>0 && myTriggerHitTime[2]>0){t_hodo_1 = (myTriggerHitTime[0]+myTriggerHitTime[2])/2.0; }
            else if(myTriggerHitTime[0]>0 && myTriggerHitTime[2]==0){t_hodo_1 = myTriggerHitTime[0]; }
            else if(myTriggerHitTime[0]==0 && myTriggerHitTime[2]>0){t_hodo_1 = myTriggerHitTime[2]; }
            // determine trigger time in H1 and H3
            if(myTriggerHitTime[1]>0 && myTriggerHitTime[3]>0){t_hodo_2 = (myTriggerHitTime[1]+myTriggerHitTime[3])/2.0; }
            else if(myTriggerHitTime[1]>0 && myTriggerHitTime[3]==0){t_hodo_2 = myTriggerHitTime[1]; }
            else if(myTriggerHitTime[1]==0 && myTriggerHitTime[3]>0){t_hodo_2 = myTriggerHitTime[3]; }
            // take the earlier time in t_hodo_1 and t_hodo_2 as the trigger time
            double eventTime = (t_hodo_1<t_hodo_2)?t_hodo_1:t_hodo_2;
            double delta_t_hodoDet = t_hodo_2 - t_hodo_1;
            // time of the first photon from muon, time of the first photon from decays        
            double t_top_first=0, t_bot_first=0;

            int nPrim = mc->GetMCParticleCount();
            int npe = mc->GetNumPE();
            //MCParticle
            TVector3 evpos = mc->GetMCParticle(0)->GetPosition(); // position in millimeter
            int pdgcode = mc->GetMCParticle(0)->GetPDGCode();     // particle pdg code
            string particlename = mc->GetMCParticle(0)->GetParticleName(); // particle name
            float evtime = mc->GetMCParticle(0)->GetTime(); // particle start time, it is zero for every event!!!

            float ke = mc->GetMCParticle(0)->GetKE(); // particle ke in MeV
            TVector3 momentum = mc->GetMCParticle(0)->GetMomentum(); // particle momentum in MeV/c
            TVector3 polarization = mc->GetMCParticle(0)->GetPolarization(); // particle polarization, not used in my study
//            float TrackStartX=evpos.x(), TrackStartY=evpos.y(), TrackStartZ=evpos.z();
             TrackStartX=evpos.x(); TrackStartY=evpos.y(); TrackStartZ=evpos.z();
            float TrackPx = momentum.x(), TrackPy=momentum.y(), TrackPz=momentum.z();
	    //cout<<"tx,y,z:"<<TrackPx<<" "<<TrackPy<<" "<<TrackPz<<endl;
            float MomentumSquare = TrackPx*TrackPx+TrackPy*TrackPy+TrackPz*TrackPz;
            float Mass = (MomentumSquare-ke*ke)/(2.0*ke);

            float TrackSpeed=TMath::Sqrt(1.0/(1.0+Mass*Mass/MomentumSquare)); // speed of the input particle, I know it since I know its momentum
	    double costheta=1./1.34/TrackSpeed;   //rong zhao to calculate the track length of optical photons 
		if(costheta>0.999)costheta=0.999;
	    double sintheta=TMath::Sqrt(1-costheta*costheta);
	    //cout<<"the theta angle:"<<costheta<<endl; 
            float X1=evpos.x(), Y1=evpos.y(), Z1=evpos.z();
            float Xprime=TrackPx/TrackPz, Yprime=TrackPy/TrackPz;
            float Xzero=X1-Xprime*Z1, Yzero=Y1-Yprime*Z1;
            float polarAngle = TMath::Pi() - TMath::ATan(TMath::Sqrt(Xprime*Xprime+Yprime*Yprime));
//            float polarAngle = TMath::Pi() - TMath::ATan(TMath::Sqrt(Xprime*Xprime+Yprime*Yprime));
            float azimuthalAngle = TMath::ATan(Yprime/Xprime)+.5*TMath::Pi();
		if(Xprime<0)azimuthalAngle+=TMath::Pi();
//            float azimuthalAngle = TMath::ATan(TMath::Sqrt(Xprime*Xprime+Yprime*Yprime));
	    TVector3 pmtline[8];
	    TVector3 pmtdisv[8];
	    TVector3 pmtiav[8];
	    TVector3 pmtmtv[8]; //the muon track and pmt norm
	    TVector3 pmtolv[8];  //the vector from cherencov light point t pmt,
	    TVector3 pmtnormv;  //the normal vector of PMT cathod
	    double pmtdm[8]={0}; //the distance from muon start to cerencov light generation
	    float dratio[8]={0};
	    for(int ii=0;ii<8;ii++){
		//pmtline[ii]=(pmtposx[ii]-X1,pmtposy[ii]-Y1,pmtposz[ii]-Z1);
	    //cout<<"the xyz11:"<<X1<<" "<<Y1<<" "<<Z1<<endl; 
		pmtline[ii][0]=pmtposx[ii]-X1;
		pmtline[ii][1]=pmtposy[ii]-Y1;
		pmtline[ii][2]=pmtposz[ii]-Z1;
		pmtdisv[ii]=pmtline[ii].Cross(momentum);
		pmtdis[ii]=pmtdisv[ii].Mag()/momentum.Mag();
		pmtol[ii]=pmtdis[ii]/sintheta;
		pmtdm[ii]=TMath::Sqrt(fabs(pmtline[ii].Mag()*pmtline[ii].Mag()-pmtdis[ii]*pmtdis[ii]))-TMath::Sqrt(fabs(pmtol[ii]*pmtol[ii]-pmtdis[ii]*pmtdis[ii]));
		//pmtdm[ii]=fabs(pmtdm[ii]);
		if(pmtdm[ii]<30)pmtdm[ii]=0;
		//if(pmtdm[ii]>8000)cout<<TMath::Sqrt(fabs(pmtline[ii].Mag()*pmtline[ii].Mag()-pmtdis[ii]*pmtdis[ii]))<<" "<<sintheta<<endl;;
		dratio[ii]=pmtdm[ii]/momentum.Mag();
		zratio1[ii]=(hz1-Z1)/TrackPz;
		hx1[ii]=X1+TrackPx*zratio1[ii];
		hy1[ii]=Y1+TrackPy*zratio1[ii];
		zratio2[ii]=(hz2-Z1)/TrackPz;
		hx2[ii]=X1+TrackPx*zratio2[ii];
		hy2[ii]=Y1+TrackPy*zratio2[ii];
		//if(pmtdm[ii]>3000)
		//cout<<pmtdm[ii]<<"ffffffff"<<dratio[ii]<<"ol:"<<pmtol[ii]<<"pmtm0.k"<<TMath::Sqrt(fabs(pmtline[ii].Mag()*pmtline[ii].Mag()-pmtdis[ii]*pmtdis[ii]))<<"pmtm1:"<<TMath::Sqrt(fabs(pmtol[ii]*pmtol[ii]-pmtdis[ii]*pmtdis[ii]))<<endl;
		crlx[ii]=X1+TrackPx*dratio[ii];  //cherencov light generated position x.
		crly[ii]=Y1+TrackPy*dratio[ii];
		crlz[ii]=Z1+TrackPz*dratio[ii];
		pmtolv[ii][0]=pmtposx[ii]-crlx[ii];
		pmtolv[ii][1]=pmtposy[ii]-crly[ii];
		pmtolv[ii][2]=pmtposz[ii]-crlz[ii];
		pmtnormv[0]=0;
		pmtnormv[1]=0;
		pmtnormv[2]=1;
		pmtiav[ii]=pmtolv[ii].Cross(pmtnormv);
		pmtmtv[ii]=momentum.Cross(pmtnormv);
		//sinphi=pmtmtv[ii].Mag()/momentum.Mag();
//cout<<"sssssssss:"<<sinphi<<endl;
		pmtia[ii]=TMath::ASin(pmtiav[ii].Mag()/pmtolv[ii].Mag());
		int iaindex=ceil(pmtia[ii]/TMath::Pi()*180.);
//		pmtatt[ii]=TMath::Exp(-pmtol[ii]/20000.)*acrylicatt[iaindex]*230000./(pmtol[ii]*pmtol[ii]);
		//pmtatt[ii]=TMath::Exp(-pmtol[ii]/20000.)*acrylicatt[iaindex]*23./(4*pmtdis[ii]+4*23);
		//pmtatt[ii]=TMath::Exp(-pmtol[ii]/20000.)*acrylicatt[iaindex]*23./(4*pmtdis[ii]+4*23);
		pmtatt[ii]=TMath::Exp(-pmtol[ii]/20000.)*acrylicatt[iaindex]*pmtiav[ii].Mag()/pmtolv[ii].Mag()*23./(4*pmtdis[ii]+4*23);
		pmtattw[ii]=TMath::Exp(-pmtol[ii]/20000.);
		pmtattao[ii]=acrylicatt[iaindex];
		pmtatts[ii]=pmtiav[ii].Mag()/pmtolv[ii].Mag()*23./(4*pmtdis[ii]+4*23);
// radius of PMT =23mm
pmtattacry[ii]=exp(-(25.4)/103.74)*getratio(25.4,103.74,pmtia[ii],1.3333,1.506);//acrylicthick=25.4,acrylicatt=103.74cookierindex=1.43
double inangle2=getangle(pmtia[ii],1.506,1.43); //acrylicrindex=1.506
//cout<<"the angle:"<<pmtia[ii]<<" ia222:"<<inangle2<<endl;
pmtattcookie[ii]=exp(-(2)/8.612)*getratio(2,8.612,inangle2,1.506,1.43);//cookiethick=2,cookieatt=8.612cookierindex=1.43
pmtattsac[ii]=pmtatts[ii]*pmtattacry[ii]*pmtattcookie[ii];
//ra2[j]=getratio(cookiethick,cookieatt[lambdaindex],inangle2,acrylicrindex[lambdaindex],cookierindex[lambdaindex]);

		
//cout<<"x:yz:"<<crlx[ii]<<" "<<crly[ii]<<" "<<crlz[ii]<<"dm:"<<pmtdm[ii]<<endl;
//		cout<<"ddd"<<ii<<" "<<pmtdis[ii]<<endl;
//		cout<<"crl:x:y:z"<<crlx[0]<<" "<<crly[0]<<" "<<crlz[0]<<endl;
		if(crlz[ii]<-550||crlz[ii]>900){pmtol[ii]=0;pmtia[ii]=0;pmtatt[ii]=0;
		pmtattw[ii]=0;
		pmtattao[ii]=0;
		pmtattacry[ii]=0;
		pmtattcookie[ii]=0;
		pmtattsac[ii]=0;
		pmtatts[ii]=0;
//hx1[ii]=-1000;
//hy1[ii]=-1000;
//hx2[ii]=-1000;
//hy2[ii]=-1000;
}
		if((crlx[ii]*crlx[ii]+crly[ii]*crly[ii])>482*482){pmtol[ii]=0;pmtia[ii]=0;pmtatt[ii]=0;
		pmtattw[ii]=0;
		pmtattao[ii]=0;
		pmtattacry[ii]=0;
		pmtattcookie[ii]=0;
		pmtattsac[ii]=0;
		pmtatts[ii]=0;
//hx1[ii]=-1000;
//hy1[ii]=-1000;
//hx2[ii]=-1000;
//hy2[ii]=-1000;
}
//		if(pmtol[ii]==0&&ii==7)cout<<"cccccc7777777"<<pmtdis[ii]<<"z: "<<crlz[ii]<<"xy:"<<crlx[ii]<<" "<<crly[ii]<<"sintheta:"<<sintheta<<"ol:"<<pmtol[ii]<<endl;
	    }
	    //cout<<"tx,y,z:"<<TrackPx<<" "<<TrackPy<<" "<<TrackPz<<endl;

            // for track reconstruction later
            float Za0 = -500.4+25.4; //800.4-25.4; // Point A at the water surface, where track goes through A
            float Xa0 = Xzero + Xprime*Za0, Ya0 = Yzero + Yprime*Za0;            
            float Za = -1118.8075; // Point B at the surface of H4 and/or H5
            float Xa = Xzero + Xprime*Za, Ya = Yzero + Yprime*Za;
            float distance1 = TMath::Sqrt((TrackStartX-Xa)*(TrackStartX-Xa)+(TrackStartY-Ya)*(TrackStartY-Ya)+(TrackStartZ-Za)*(TrackStartZ-Za));
            float time_stage1 = distance1/(TrackSpeed*300.0); // in ns
		float hz[6]={1144.7,1133.50,896.5,889.7,-1123.9+d4z,-1123.9};
		float hx[6]={0,0,0,0,0,0};
		float hy[6]={0,0,0,0,0,0};
		for(int ii=0;ii<6;ii++){
             	hx[ii] = Xzero + Xprime*hz[ii],
		 hy[ii] = Yzero + Yprime*hz[ii];
		}

            /*
            cout<<"Particle\t"<<particlename<<"\tmass\t"<<Mass<<endl;
            cout<<"Position (mm)\t"<<evpos.x()<<"\t"<<evpos.y()<<"\t"<<evpos.z()<<endl;
            cout<<"Ke (MeV)\t"<<ke<<endl;
            cout<<"Momentum (MeV)\t"<<momentum.x()<<"\t"<<momentum.y()<<"\t"<<momentum.z()<<endl;
            cout<<"Speed (c)\t"<<TrackSpeed<<endl;
            */ 

            //TVector3 energyCentroid = mc->GetMCSummary()->GetEnergyCentroid(); 
            /**
            * Centroid of energy loss.
            *
            * This is the average position of all steps in this event, weighted by
            * the energy lost in that step. Optical photons, and the rest mass of
            * particles when a track terminates are not included.
            */
            //TVector3 energyRMS = mc->GetMCSummary()->GetEnergyRMS();
            //float totalScintEdep = mc->GetMCSummary()->GetTotalScintEdep(); // in MeV
			
            //MCTrack
            int mctrackcnt = mc->GetMCTrackCount();
            double decayElectronEnergy=Mass;
            int mydecayflag = 0;
            double stopx=-1000, stopy=-1000, stopz=-10000, gTime=-1,decaytime=-1.;
	    //cout<<"mctrackcnt="<<mctrackcnt<<endl;
            if(mctrackcnt!=0){ // this means a decay happend

                for(int imctrackcnt=0;imctrackcnt<mctrackcnt;imctrackcnt++){
                    int mctrackstepcnt = mc->GetMCTrack(imctrackcnt)->GetMCTrackStepCount();
					//fstream fout("track_output_test.txt",ios::out | ios::app);
                    //fout<<imctrackcnt<<"\t"<<mc->GetMCTrack(imctrackcnt)->GetID()<<"\t"<<mc->GetMCTrack(imctrackcnt)->GetParentID()<<"\t"<<mctrackstepcnt<<endl;
                    //fout<< "---->Track " << imctrackcnt << " id " << mc->GetMCTrack(imctrackcnt)->GetID() 
                    //    << " pName " << mc->GetMCTrack(imctrackcnt)->GetParticleName() << " parent " << mc->GetMCTrack(imctrackcnt)->GetParentID()<< " stepcnt "<< mctrackstepcnt <<endl;
                    //fout<<"------->stepcnt "<<mctrackstepcnt<<endl;
                    decayElectronEnergy -= mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(0)->GetKE();
                    if(imctrackcnt==0 && mydecayflag==0){
                        stopx = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().x();
                        stopy = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().y();
                        stopz = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetEndpoint().z();
			mustopz=stopz;
                        gTime = mc->GetMCTrack(0)->GetMCTrackStep(0)->GetGlobalTime();
			decaytime=gTime;
                        mydecayflag=1;
                    }

                    //for(int jmcstepcnt=0; jmcstepcnt<mctrackstepcnt; jmcstepcnt++){
                    //    TVector3 stepEndPoint = mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetEndpoint();
                    //    fout<< "----------> steps " << jmcstepcnt << " ke "<< mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetKE()
                    //        << " end_points (" <<stepEndPoint[0]<<", "<<stepEndPoint[1]<<", "<<stepEndPoint[2]<<")" << " time " << mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetGlobalTime()<<" process "<< mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetProcess() << " volume " << mc->GetMCTrack(imctrackcnt)->GetMCTrackStep(jmcstepcnt)->GetVolume() <<endl;
					//}
                }
            } //
            else mustopz = -1500;
	mumumentum=momentum.Mag();

            // define trigger types
            int trigflag = 0;
            int trigflag_sub = 0;
            int nHits = ev->GetPMTCount();
            //if(nHits>6){ // multiplicity triggers
                //  trigflag = 1; istrigger[trigflag]=true; trigcnts[trigflag]++;
                //hEvtTime_multiTrig->Fill(eventTime);
            //}
            // through going triggers
            hMuMomentum_all->Fill(TMath::Sqrt(MomentumSquare));
            if( myTriggerEdep[2]==0 && myTriggerEdep[0]>0 && myTriggerEdep[3]>0 && myTriggerEdep[1]>0 ){
                trigcnts_hodo++;
                if(myTriggerEdep[4]>0. && myTriggerEdep[5]>0.){
                    trigcnts_hodo_H5_H4++;
                    if(mydecayflag==1) trigcnts_decays++;
                    hEdepH4_vs_EdepH5->Fill(myTriggerEdep[4],myTriggerEdep[5]);
                    //fstream testout("testout.txt",ios::out|ios::app);
                    //testout<<ffff<<"\t"<<i<<"\t"<<pdgcode<<"\t"<<X1<<"\t"<<Y1<<"\t"<<TrackPx<<"\t"<<TrackPy<<"\t"<<TrackPz<<endl;
                    //testout<<"\t\t Edep in H4, H5: "<<myTriggerEdep[4]<<"\t"<<myTriggerEdep[5]<<endl;
                    //testout.close();
                    hTDC_H4H5->Fill(myTriggerHitTime[4]-myTriggerHitTime[5]);
                }
            }

            //now make my trigger conbinations based on the normal hodo trigger
            int PsTrigFlag[6]={0,0,0,0,0,0};
            for(int pscnt=0; pscnt<6; pscnt++){
              if(myTriggerEdep[pscnt]>0){
                PsTrigFlag[pscnt]=1;
              }
              else PsTrigFlag[pscnt]=0;
            }

            thistriggerflag =  PsTrigFlag[2]*TMath::Power(2,5)
                                  +PsTrigFlag[0]*TMath::Power(2,4)
                                  +PsTrigFlag[3]*TMath::Power(2,3)
                                  +PsTrigFlag[1]*TMath::Power(2,2)
                                  +PsTrigFlag[4]*TMath::Power(2,1)
                                  +PsTrigFlag[5]*TMath::Power(2,0);
            //cout<<ffff<<"\t"<<i<<"\t"<<PsTrigFlag[2]<<PsTrigFlag[0]<<PsTrigFlag[3]<<PsTrigFlag[1]<<PsTrigFlag[4]<<PsTrigFlag[5]<<endl;
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                int testtrignum = std::strtol(strig[trigindex].c_str(), 0, 2);//astring.Atoi();
                //cout<<strig[trigindex]<<"\t"<<testtrignum<<"\t"<<thistriggerflag<<endl;
	//		cout<<"arlier test"<<thistriggerflag<<"edep"<<myTriggerEdep[1]<<endl;
                if(thistriggerflag == testtrignum){
                    hEndXY[trigindex]->Fill(Xa, Ya); hEndXY0[trigindex]->Fill(Xa0, Ya0);
                    istrigger[trigindex]=true;
                    trigcnts[trigindex]++;
	//		cout<<"arlier test"<<endl;
                    eventStartTime[trigindex] = myTriggerHitTime[indexMap_forEvtStartTime[trigindex]];
                    hTrigTime[trigindex]->Fill(eventStartTime[trigindex]);
                    if(mydecayflag==1)hDecayTime[trigindex]->Fill(decaytime);  //rong zhao
                    if(mydecayflag==1)hDecayTime_all->Fill(decaytime);  //rong zhao

                    break;
                }
            }

	    bool testTrig = false;			
	    //define triggers 
            if( (myTriggerEdep[2]>0 || myTriggerEdep[0]>0) && (myTriggerEdep[3]>0 || myTriggerEdep[1]>0)  // the normal hodoTrig
              ){
            //select in-time events only, time coincidence is made by constraining time difference between hit times in plastic scintillators
                if(myTriggerEdep[4]==0 && myTriggerEdep[5]>0.0) trigcnts_hodo_H5_notH4++;
                if(myTriggerEdep[4]>0. && myTriggerEdep[5]==0.) trigcnts_hodo_notH5_H4++;
                if(myTriggerEdep[4]==0. && myTriggerEdep[5]==0.) trigcnts_hodo_notH5_notH4++;
                if(myTriggerEdep[4]>0. && myTriggerEdep[5]>0.){
                    trigcnts_hodo_H5_H4++;}
                int index1 = -1, index2=-1;
                if(myTriggerEdep[0]>0 && myTriggerEdep[2]<=0) index1 = 0;
                else index1=2;
                if(myTriggerEdep[3]>0 && myTriggerEdep[1]<=0) index2 = 3;
                else index2=1;
                //if(  myTriggerHitTime[index2] - myTriggerHitTime[index1] > 0 
                //  && myTriggerHitTime[index2] - myTriggerHitTime[index1] < 10
                //  //&& myTriggerHitTime[5]-myTriggerHitTime[index1] > 0
                //  )
                testTrig = true;
                				
                trigflag = 0;
                if(testTrig == true){
                  //istrigger[trigflag]=true; trigcnts[trigflag]++;
                  hEvtTime[trigflag]->Fill(eventTime);
		  hDecayZ[trigflag]->Fill(stopz);
                }

		if( 1 ) {
		    // get time difference between plastic scintillators
   			if(myTriggerEdep[2]>0 && myTriggerEdep[0]>0){
	    			if(myTriggerHitTime[2]<=0 || myTriggerHitTime[0]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H0H2->Fill(myTriggerHitTime[0] - myTriggerHitTime[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[3]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[3]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H3H2->Fill(myTriggerHitTime[3] - myTriggerHitTime[2]);
       			}
		    	if(myTriggerEdep[2]>0 && myTriggerEdep[1]>0){
			    	if(myTriggerHitTime[2]<=0 || myTriggerHitTime[1]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
				    hTDC_H1H2->Fill(myTriggerHitTime[1] - myTriggerHitTime[2]);
					hEdep_vs_TDC_H1H2->Fill(myTriggerHitTime[1] - myTriggerHitTime[2], myTriggerEdep[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[4]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[4]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H4H2->Fill(myTriggerHitTime[4] - myTriggerHitTime[2]);
    			}
	    		if(myTriggerEdep[2]>0 && myTriggerEdep[5]>0){
		    		if(myTriggerHitTime[2]<=0 || myTriggerHitTime[5]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H5H2->Fill(myTriggerHitTime[5] - myTriggerHitTime[2]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[3]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[3]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H3H0->Fill(myTriggerHitTime[3] - myTriggerHitTime[0]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[1]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[1]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H1H0->Fill(myTriggerHitTime[1] - myTriggerHitTime[0]);
			    }
    			if(myTriggerEdep[0]>0 && myTriggerEdep[4]>0){
	    			if(myTriggerHitTime[0]<=0 || myTriggerHitTime[4]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
		    		hTDC_H4H0->Fill(myTriggerHitTime[5] - myTriggerHitTime[0]);
    			}
	    		if(myTriggerEdep[0]>0 && myTriggerEdep[5]>0){
		    		if(myTriggerHitTime[0]<=0 || myTriggerHitTime[5]<=0) cout<<"Interesting mismatch edep>0 and hitTime=0"<<endl;
			    	hTDC_H5H0->Fill(myTriggerHitTime[5] - myTriggerHitTime[0]);
			    }
    		}				
            }
	    testTrig=false;
			
            // fill histograms of (1) edep and hittime in plastic scintillators for the different trigger types
            //                    (2) number of hitting PMTs, for the different trigger types
            //                    (3) the initial muons momentum and angles.
            for(int index=0; index<nbTrigs; index++){
                if(istrigger[index]==true){
                    //hMuMomentum[index]->Fill(TMath::Sqrt(MomentumSquare));

                    hMuPAngle[index]->Fill(polarAngle);
                    hMuAAngle[index]->Fill(azimuthalAngle);
                    for(int psid=0; psid<6; psid++){
                        hEdepInPS[psid][index]->Fill(myTriggerEdep[psid]);
                        hHitTimeInPS[psid][index]->Fill(myTriggerHitTime[psid]);
                    }
                }
            }

            //Dealing with photons in PMTs. Getting event info from MCphoton
            double totQ1 = 0; // total charge of all pmts in npe, cut photon arrival time to <2500 ns
            double totQ_top1 = 0, totQ_bot1 = 0;
            // use vectors to hold pulses (charge and time) in the PMTs:
	    	// pmtnb=8 columns for the 8 PMTs, each contains 0 photons when initialized
            std::vector< std::vector<double> > pulseTimeInPMTs1; // pulse time
            pulseTimeInPMTs1.resize(pmtnb,std::vector<double>(0,0)); 
            std::vector< std::vector<double> > pulseChargeInPMTs1; // pulse charge
            pulseChargeInPMTs1.resize(pmtnb,std::vector<double>(0,0)); 

            double topcharge_mu=0, topcharge_e=0;
            double botcharge_mu=0, botcharge_e=0;
            double totalcharge_mu=0, totalcharge_e=0;
            int nhit_topPMT=0, nhit_botPMT=0;
            double top_bot_ratio_eSig=0, top_bot_ratio_muSig=0;
            bool hitflag_topPMT=false, hitflag_botPMT=false;

            //cout<<ffff<<"\t"<<i<<"\t"<<mc->GetMCPMTCount()<<endl;
            double thispmt_charge[8], thispmt_charge_50ns[8], thispmt_charge_40ns[8],thispmt_charge_decay[8];
            for(int t=0; t<8; t++){ thispmt_charge[t]=0; thispmt_charge_50ns[t]=0; thispmt_charge_40ns[t]=0;thispmt_charge_decay[t]=0;}//rong zhao add the decay npe
            //in simulation, if a pmt has no hit it will not record, so mc->GetMCPMTCount()<=8.
            for (int imcpmt=0; imcpmt < mc->GetMCPMTCount(); imcpmt++) { // for each pmt
                RAT::DS::MCPMT *mcpmt = mc->GetMCPMT(imcpmt);
                int pmtID = mcpmt->GetID(); // pmt id
                int nPhotons = mcpmt->GetMCPhotonCount(); // number of photons in this pmt
                int firstphotonflag = 0;
		double firstphotontime = 0;
		//See if there's any no-photon-hit PMT in this event.
                //if(nPhotons<1){cout<<"Evt "<<i<<", pmt "<<pmtID<<", no photons."<<endl; continue;}
                double thispmt_charge_afterdecay=0;
                double iTime_firstPhoton = 0;
                double ichargesum = 0;
		double lasttime=0;
                for (int iphoton=0; iphoton < nPhotons; iphoton++)  { // for all the photons produced in this pmt
                    double iTime = mcpmt->GetMCPhoton(iphoton)->GetHitTime()-evtime; // global (absolute) time of the photon
                    if(iphoton==0) iTime_firstPhoton = iTime;
                    double iCharge = mcpmt->GetMCPhoton(iphoton)->GetCharge()/converterfactor[pmtID]; // charge of this photon (applied my PMT SPE spectrum)
                    //don't consider photons that arrive later than 2500 ns, 
		    //also I set a 0.5 pe threshold, even if it is a photon, but if its charge is <0.5 pe, don't consider
                    //the charge could be <1 pe because the PMT SPE response is used
			if(iTime-lasttime<10){
				ichargesum+=iCharge;
				}
			else{
			if(iphoton>0){hPMTPhotonTimeSpannpe[pmtID][0]->Fill(lasttime,ichargesum);
					}
			else{hPMTPhotonTimeSpannpe[pmtID][0]->Fill(iTime,iCharge);}
			ichargesum=iCharge;
			}
			lasttime=iTime;
			//hPMTPhotonTimeSpannpe[pmtID][0]->Fill(iTime,iCharge);
                    if(iTime>2500 || (iCharge<0.5 && iphoton==0)) {continue;} //rong modified,  for fitting MC, the cut should not be applied here. 
                    pulseTimeInPMTs1[pmtID].push_back(iTime);
                    pulseChargeInPMTs1[pmtID].push_back(iCharge);
		    totQ1 += iCharge;
    	            thispmt_charge[pmtID] += iCharge;
    	            if(iTime<=(50+iTime_firstPhoton)&&iTime<300){// rong modified, since some photon from decay event could be the first photon of one waveform.
                    //if(iTime<=(50+iTime_firstPhoton)){
    	                thispmt_charge_50ns[pmtID]+=iCharge;
    	            }
		          else if(iphoton<=3){
    	                thispmt_charge_decay[pmtID]+=iCharge;
			         }
    	            if(iTime<=(40+iTime_firstPhoton)){
    	                thispmt_charge_40ns[pmtID]+=iCharge;
    	            }
	            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
			    hPMTPhotonTimeSpan[pmtID][trigindex]->Fill(iTime);
			 //	hPMTPhotonTimeSpannpe[pmtID][trigindex]->Fill(iTime,iCharge);
			//if(istrigger[trigindex]==true){
			//    hPMTPhotonTimeSpan[pmtID][trigindex]->Fill(iTime);
			//}
		    }
                    if(pmtID==4||pmtID==5) { // in top PMTs
			             totQ_top1 += iCharge;
                    }
                    else{ // in bottom PMTs
                        totQ_bot1 += iCharge;
                    }
                    firstphotonflag++;
                    if(firstphotonflag==1) {
			firstphotontime = iTime;
                        for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                            if(istrigger[trigindex]==true ){ // first photon in this PMT in this event 
                                hTime[pmtID][trigindex]->Fill(iTime);
                            }
                        }
                    }
                }//end of getting photons in one pmt
		for(int trigindex=0; trigindex<nbTrigs; trigindex++){
		    if(istrigger[trigindex]==true){
            hMuMomentum_hodotrig->Fill(TMath::Sqrt(MomentumSquare));
            if(mydecayflag==1)hMuMomentum_decay->Fill(TMath::Sqrt(MomentumSquare));
			if(firstphotonflag==1){
			    hPMTNpeVsTime[pmtID][trigindex]->Fill(thispmt_charge[pmtID], firstphotontime);							
			}
		    }
		}
            }//end of getting photons from all PMTs
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
		pmtnpe_all=0;
		pmtnpe_alldecay=0;
                for(int pmtID=0; pmtID<8; pmtID++){
			pmtnpe[pmtID]=-1; //rong zhao 
			//pmtol[pmtID]=-1; //rong zhao 
			pmtdecaynpe[pmtID]=-1; //rong zhao 
                    if(istrigger[trigindex]==true){
//                        hCharge[pmtID][trigindex]->Fill(thispmt_charge[pmtID]);
//                        calculate decay charge instead of muon charge 
                        hCharge[pmtID][trigindex]->Fill(thispmt_charge[pmtID]);
                        hCharge1[pmtID][trigindex]->Fill(thispmt_charge_50ns[pmtID]*npescalefactor[pmtID]);
                        hCharge2[pmtID][trigindex]->Fill(thispmt_charge_40ns[pmtID]);
			pmtnpe[pmtID]=thispmt_charge_50ns[pmtID];// rong zhao
			//pmtol[pmtID]=
			pmtdecaynpe[pmtID]=thispmt_charge_decay[pmtID];// rong zhao
                    }
		    pmtnpe_all+=thispmt_charge_50ns[pmtID];//rong zhao
		    pmtnpe_alldecay+=thispmt_charge_decay[pmtID];
                	int testtrignum = std::strtol(strig[trigindex].c_str(), 0, 2);//astring.Atoi();
                	if(pmtID==7&&-600<stopz&&stopz<1000&&pmtnpe_alldecay>0){
                    //if(pmtID==7&&-600<stopz&&stopz<1000&&pmtnpe_alldecay>0){
                	//if(pmtID==7&&-600<stopz&&stopz<600&&pmtnpe_alldecay>5&&thispmt_charge_decay[4]>0){
                	//if(pmtID==7&&-600<stopz&&stopz<600&&pmtnpe_alldecay>0&&thispmt_charge_decay[4]>0&&thispmt_charge_decay[5]>0){
                	//if(thistriggerflag == testtrignum&&pmtID==7&&pmtnpe_alldecay>0&&-600<stopz&&stopz<600){
                	//if(thistriggerflag == testtrignum&&trigindex==nbTrigs-1){
		    	  hz_npe[trigindex]->Fill(pmtnpe_all,stopz);
		    	  hz_npe_all->Fill(stopz,pmtnpe_all);
                  hz_mup_all->Fill(stopz,TMath::Sqrt(MomentumSquare));
		    	  //if(pmtnpe_alldecay>0.0)hz_npe_alldecay->Fill(pmtnpe_alldecay,stopz);
		    	  //if(pmtnpe_alldecay>0.0)hmu_npe_alldecay->Fill(pmtnpe_all,stopz);
		    	  hz_npe_alldecay->Fill(pmtnpe_alldecay,stopz);
		    	  if(pmtnpe_alldecay>0)hmu_npe_alldecay->Fill(pmtnpe_all,pmtnpe_alldecay);
		    	  hz_npe_allemu->Fill(pmtnpe_all,pmtnpe_alldecay);}
		    	  hrp_npe[trigindex]->Fill(pmtnpe_all,TMath::Sqrt(TMath::Power(stopx-280,2)+TMath::Power(stopy-55,2)));
		    	  if(stopz>-400)hrp_npe_all->Fill(pmtnpe_all,TMath::Sqrt(TMath::Power(stopx-280,2)+TMath::Power(stopy-55,2)));}
		//if(((trigindex)%3==0)&&(istrigger[trigindex]==true)){    //rong zhao fill the tree with selected trig type decay events
		//if(((trigindex+1)%3==2)&&(istrigger[trigindex]==true)){    //rong zhao fill the tree with selected trig typeH5
		//if((trigindex==25)&&(istrigger[trigindex]==true)){    //rong zhao fill the tree with selected trig type
	//int hodoselect=0;  //0, for Hodo+H6, 1 for Hodo+H5, 2 for Hodo only, 3, for Hodo26, 4 for Hodo25
	//if((pmtnpe_alldecay>5)&&(istrigger[trigindex]==true)){    //rong zhao for the decay events
//if(pmtdm[4]>0)cout<<"x:yz:"<<crlx[4]<<" "<<crly[4]<<" "<<crlz[4]<<"dm:"<<pmtdm[4]<<endl;
		bool hitinh4=0;
		bool hitinh5=0;
		if(hx[4]<-174.6+d4x+153.99&&hx[4]>-174.6+d4x-153.99&&hy[4]<-141.7+d4y+358.4&&hy[4]>-141.7+d4y-358.4){
			hitinh4=1;
			}
		if(hx[5]<179.4+152.4&&hx[5]>179.4-152.4&&hy[5]<-131.7+357.9&&hy[5]>-131.7-357.9){
			hitinh5=1;
			}
		if((istrigger[trigindex]==true)){    //rong zhao fill the tree with selected trig type
			int saveflag=1;
		switch(hodoselect){
		case 0:
			if((trigindex+1)%3!=0){saveflag=0;} //H5
//			if((trigindex+1)%3!=0&&trigindex%3!=0){saveflag=0;}
			break;
		case 1:
			if((trigindex+1)%3!=2){saveflag=0;}//H4
//			if((trigindex+1)%3!=2&&trigindex%3!=0){saveflag=0;}
			break;
		case 2:
			if((trigindex+1)%3!=1){saveflag=0;} //Hodo only
//			if((trigindex+1)%3!=1){saveflag=0;} //Hodo only
//			if((trigindex+1%3)!=1){saveflag=0;} //Hodo only
//			if((trigindex+1)%3==1&&trigindex%3!=0){saveflag=0;}
			break;
		case 3:
			if(trigindex!=19){saveflag=0;}
			break;
		case 4:
			if(trigindex!=20){saveflag=0;}
			break;
		case 5:
			if(trigindex!=24&&trigindex!=25&&trigindex!=26){saveflag=0;}
			break;
		case 6:
			if(trigindex!=21&&trigindex!=22&&trigindex!=23){saveflag=0;}
			break;
		case 7:
			if(trigindex!=18&&trigindex!=19&&trigindex!=20){saveflag=0;}
			break;
		case 11:
			if(trigindex!=1&&trigindex!=2){saveflag=0;}
			break;
		case 12:
			if(trigindex!=4&&trigindex!=5){saveflag=0;}
			break;
		case 13:
			if(trigindex!=7&&trigindex!=8){saveflag=0;}
			break;
		case 14:
			if(trigindex!=10&&trigindex!=11){saveflag=0;}
			break;
		case 15:
			if(trigindex!=13&&trigindex!=14){saveflag=0;}
			break;
		case 16:
			if(trigindex!=16&&trigindex!=17){saveflag=0;}
			break;
		case 17:
			if(trigindex!=19&&trigindex!=20){saveflag=0;}
			break;
		case 18:
//			if(trigindex!=22&&trigindex!=23){saveflag=0;}
//			H4 check
			if(hitinh4&&testTrig){saveflag=1;}
			break;
		case 19:
//			H5 check
			if(hitinh5&&testTrig){saveflag=1;}
//			if(trigindex!=25&&trigindex!=26){saveflag=0;}
			
			break;
		case 20:
			saveflag=1;break;
		default:
			if(trigindex!=26){saveflag=0;}
		}
		if(saveflag==0)break;	
            for(int index=0; index<6; index++){
		hhodoXYtrigger[index]->Fill(hx[index],hy[index]);
            }
    		outrootnpe->cd();
//cout<<"pmtol"<<pmtol[0]<<" "<<pmtol[1]<<endl;
                    hMuAAngletrigger->Fill(azimuthalAngle*180./TMath::Pi());
                    hMuPAngletrigger->Fill(polarAngle*180./TMath::Pi());
//            for(int index=0; index<6; index++){
//		hhodoXYtrigger[index]->Fill(hx[index],hy[index]);
 //           }
		muangle=azimuthalAngle*180./TMath::Pi();
		npetree->Fill();	
		oltree->Fill();	
		//npedecaytree->Fill();	
		npedecaytree->Fill();	
    		foutroot->cd();
			}
	    }
		
	    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
		if(istrigger[trigindex]==true){
		    hTotalCharge[trigindex]->Fill(totQ1);
		    //hTotalChargeTop[trigindex]->Fill(totQ_top1);
	        //hTotalChargeBottom[trigindex]->Fill(totQ_bot1);
		    if(totQ1>0){
			hNpeRatio[trigindex]->Fill(totQ_top1/totQ1*(2/8)); // charge ratio, weighted by number of PMTs.
	            }
		}
	    }
            // get the time difference between pmts, for all photons
		int coincidence_flag=0;
		float coincidence_time=0;
	    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
            hMuStopz_all->Fill(stopz);
                for(int ngroup = 0; ngroup<7; ngroup++){
                    for(int nbpulseS0=0; nbpulseS0<pulseTimeInPMTs1[ngroup].size(); nbpulseS0++){
                        for(int pmtid=ngroup+1;pmtid<8;pmtid++){
                            for(int nbpulseS1=0; nbpulseS1<pulseTimeInPMTs1[pmtid].size(); nbpulseS1++){
                                double dT = (pulseTimeInPMTs1[pmtid][nbpulseS1]-meanPulseStartBin[pmtid])
            								- (pulseTimeInPMTs1[ngroup][nbpulseS0]-meanPulseStartBin[ngroup]);
				double corrt1=pulseTimeInPMTs1[pmtid][nbpulseS1]-meanPulseStartBin[pmtid];
				double corrt2=pulseTimeInPMTs1[ngroup][nbpulseS0]-meanPulseStartBin[ngroup];
                                if(ngroup==0&&corrt1>100&&corrt2>100) hPulseTimeToS0[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==1&&corrt1>100&&corrt2>100) hPulseTimeToS1[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==2&&corrt1>100&&corrt2>100) hPulseTimeToS2[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==3&&corrt1>100&&corrt2>100) hPulseTimeToS3[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==4&&corrt1>100&&corrt2>100) hPulseTimeToS4[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==5&&corrt1>100&&corrt2>100) hPulseTimeToS5[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==6&&corrt1>100&&corrt2>100) hPulseTimeToS6[pmtid-1-ngroup][trigindex]->Fill( corrt1,corrt2 );
                                if(ngroup==0) hPulseTimeDiffToS0[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==1) hPulseTimeDiffToS1[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==2) hPulseTimeDiffToS2[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==3) hPulseTimeDiffToS3[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==4) hPulseTimeDiffToS4[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==5) hPulseTimeDiffToS5[pmtid-1-ngroup][trigindex]->Fill( dT );
                                if(ngroup==6) hPulseTimeDiffToS6[pmtid-1-ngroup][trigindex]->Fill( dT );
                                int coincidence_tag1=-1;
                                int coincidence_tag2=-1;
                                //hMuStopz_all->Fill(stopz);
				if(fabs(dT)<50&&corrt1>100&&corrt2>100){
                                    coincidence_tag1=(ngroup%4)*(ngroup%5);
                                    coincidence_tag2=(pmtid%4)*(pmtid%5);
                                    coincidence_flag=1;
                                    hMuStopz_coincidence->Fill(stopz);
                                    if(coincidence_tag1==0&&coincidence_tag2==0)hMuStopz_topcoincidence->Fill(stopz);
                                    else if(coincidence_tag1*coincidence_tag2==0)hMuStopz_tbcoincidence->Fill(stopz);
                                    if(coincidence_tag1>0&&coincidence_tag2>0)hMuStopz_botcoincidence->Fill(stopz);
                                    coincidence_time=(corrt1+corrt2)*.5;
                                    //hcoincidence_npe[pmtid]->Fill(pulseChargeInPMTs1[pmtid][nbpulseS1]);
                                }
                            }
                        }
                    }
                }
	    }
		//fill the reconstruncted time histograms.
		if(coincidence_time>0)hcoincidence_time->Fill(coincidence_time);
		if(coincidence_time>100){   //fill the nPE histogram of slected decay events.
            hMuMomentum_coincidence->Fill(TMath::Sqrt(MomentumSquare));
            for(int pmtid=0;pmtid<8;pmtid++){
                //for(int nbpulseS1=0; nbpulseS1<pulseTimeInPMTs1[pmtid].size(); nbpulseS1++){
                //    hcoincidence_npe[pmtid]->Fill(pulseChargeInPMTs1[pmtid][nbpulseS1]); 
                //}
                //if(pulseTimeInPMTs1[pmtid].size()==0)hcoincidence_npe[pmtid]->Fill(0);
               hcoincidence_npe[pmtid]->Fill(thispmt_charge_decay[pmtid]); 
                // cout<<"test coincidence npe"<<thispmt_charge_50ns[pmtid]<<endl;
            }
 
        }
            //reset the flags
            for(int trigindex=0; trigindex<nbTrigs; trigindex++){
                if(coincidence_time>100){
                    hMuMomentum[trigindex]->Fill(TMath::Sqrt(MomentumSquare));
                    
                    hTotalChargeratio_stopz[trigindex]->Fill(totQ_top1*1./(totQ_bot1+totQ_top1),stopz);
                    hTotalChargeTop[trigindex]->Fill(totQ_top1);
                     hTotalChargeBottom[trigindex]->Fill(totQ_bot1);}
                istrigger[trigindex]=false;
                eventStartTime[trigindex]=0;
	    }				
  //  foutroot->cd();
  //  npetree->Write();
  //  foutroot->cd();
        }  // end of for loop reading ONE root file
        f->Close(); // close current root tree file
    }// end of loop for reading ALL root files

	//Now beging writing histograms in the output root file
    foutroot->cd();
    hMuStopz_all->Write();
    hMuStopz_tbcoincidence->Write();
    hMuStopz_topcoincidence->Write();
    hMuStopz_coincidence->Write();
    hMuStopz_botcoincidence->Write();
    hMuMomentum_all->Write();
    hMuMomentum_decay->Write();
    hMuMomentum_hodotrig->Write();
    hMuMomentum_coincidence->Write();
    std::cout<<"Done reading files, now writting to ROOT file ..."<<std::endl;
	dir[0]->cd();
	//plastic scintillators (hodo detectors)
	TCanvas* cEdepInPS = new TCanvas("cEdepInPS","",1200,900);
	TLegend* legEdepInPS = new TLegend(0.6,0.3,0.9,0.5);
	TCanvas* cHitTimeInPS = new TCanvas("cHitTimeInPS","",1200,900);
	TLegend* legHitTimeInPS = new TLegend(0.6,0.3,0.9,0.5);
    for(int i=0; i<6; i++){
        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
			cEdepInPS->cd(); 
			hEdepInPS[i][j]->Draw("sames");
			hEdepInPS[i][j]->SetLineColor(i+1);
			legEdepInPS->AddEntry(hEdepInPS[i][j],hEdepInPS[i][j]->GetName(),"l");
			
			cHitTimeInPS->cd(); 
			hHitTimeInPS[i][j]->Draw("sames");
			hHitTimeInPS[i][j]->SetLineColor(i+1);
			legHitTimeInPS->AddEntry(hHitTimeInPS[i][j],hHitTimeInPS[i][j]->GetName(),"l");
			
            hEdepInPS[i][j]->Write();
            hHitTimeInPS[i][j]->Write();
        }
    }
	cEdepInPS->cd(); legEdepInPS->Draw();  cEdepInPS->Write(); cEdepInPS->Close();
	cHitTimeInPS->cd(); legHitTimeInPS->Draw(); cHitTimeInPS->Write(); cHitTimeInPS->Close();
	
	hTDC_H0H2->Write();
	hTDC_H3H2->Write();
	hTDC_H1H2->Write();
	hTDC_H4H2->Write();
	hTDC_H5H2->Write();
	hTDC_H3H0->Write();
	hTDC_H1H0->Write();
	hTDC_H4H0->Write();
	hTDC_H5H0->Write();
	hTDC_H4H5->Write();
	hEdep_vs_TDC_H1H2->Write();
	hEdepH4_vs_EdepH5->Write();

	//input muon information
	TCanvas* cMuonInfo = new TCanvas("cMuonInfo","",1200,900);
	cMuonInfo->Divide(2,1);
	cMuonInfo->cd(2);
	TLegend* legMuonAngles = new TLegend(0.6,0.3,0.9,0.5);

	for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
	    cMuonInfo->cd(1); 
	    hMuMomentum[j]->Draw();
	    hMuMomentum[j]->Write();
	    cMuonInfo->cd(2);
	    hMuPAngle[j]->Draw("sames");
	    hMuPAngle[j]->SetLineColor(1);
	    hMuAAngle[j]->Draw("same");
	    hMuAAngle[j]->SetLineColor(2);
	    legMuonAngles->AddEntry(hMuPAngle[j],hMuPAngle[j]->GetName(),"l");
	    legMuonAngles->AddEntry(hMuAAngle[j],hMuAAngle[j]->GetName(),"l");
            hMuPAngle[j]->Write();
            hMuAAngle[j]->Write();
	    hEvtTime[j]->Write(); // the time defined by hodoDetectors
	    hTrigTime[j]->Write();
	    hDecayTime[j]->Write();
	    hDecayTime_all->Write();
	    hcoincidence_time->Write();
	}
	cMuonInfo->cd(2); legMuonAngles->Draw(); cMuonInfo->Write();cMuonInfo->Close();
	

    for(int i=0; i<pmtnb; i++){

        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
             hcoincidence_npe[i]->SetLineWidth(2);
             hcoincidence_npe[i]->Write();
	    if(i>0) hPulseTimeDiffToS0[i-1][j]->Write();
            if(i>1) hPulseTimeDiffToS1[i-2][j]->Write();
            if(i>2) hPulseTimeDiffToS2[i-3][j]->Write();
            if(i>3) hPulseTimeDiffToS3[i-4][j]->Write();
            if(i>4) hPulseTimeDiffToS4[i-5][j]->Write();
            if(i>5) hPulseTimeDiffToS5[i-6][j]->Write();
            if(i>6) hPulseTimeDiffToS6[i-7][j]->Write();
	    if(i>0) hPulseTimeToS0[i-1][j]->Write();
            if(i>1) hPulseTimeToS1[i-2][j]->Write();
            if(i>2) hPulseTimeToS2[i-3][j]->Write();
            if(i>3) hPulseTimeToS3[i-4][j]->Write();
            if(i>4) hPulseTimeToS4[i-5][j]->Write();
            if(i>5) hPulseTimeToS5[i-6][j]->Write();
            if(i>6) hPulseTimeToS6[i-7][j]->Write();
        }
    }

    for(int i=0; i<pmtnb; i++){
        for(int j=0; j<nbTrigs; j++){
            dir[j]->cd();
            hCharge[i][j]->Write();
            hCharge1[i][j]->Write();
            hCharge2[i][j]->Write();
	    hTime[i][j]->Write();
            hPMTNpeVsTime[i][j]->Write();
            hPMTPhotonTimeSpan[i][j]->Write();
            hPMTPhotonTimeSpannpe[i][j]->Write();
        }
    }
	
	for(int j=0; j<nbTrigs; j++){
	    dir[j]->cd();
	    hTotalCharge[j]->Write();
	    hTotalChargeTop[j]->Write();
        hTotalChargeratio_stopz[j]->Write();
	    hTotalChargeBottom[j]->Write();
	    hNpeRatio[j]->Write();	
            hDecayZ[j]->Write();
            hEndXY[j]->Write();
            hEndXY0[j]->Write();
            hz_npe[j]->Write();
            hz_npe_all->Write();
            hz_mup_all->Write();
            hz_npe_alldecay->Write();
            hmu_npe_alldecay->Write();
            hz_npe_allemu->Write();
            hrp_npe[j]->Write();
            hrp_npe_all->Write();
		can1->cd();
		hz_npe[j]->Draw("colz");
		can1->SetLogz();
		can1->SaveAs(outpdf);
		hrp_npe[j]->Draw("colz");
		can1->SetLogz();
		can1->SaveAs(outpdf);

	}
  		can1->Print(outpdf_end);

		TCanvas* cz_npe = new TCanvas();
		cz_npe->Divide(2,1);
		cz_npe->cd(1);
	    	cz_npe->SetName("chz_npe_all");
    		hz_npe_all->Draw("colz"); hz_npe_all->SetMarkerStyle(8);
		cz_npe->SetLogz();
		cz_npe->cd(2);
		hz_npe_all1d=hz_npe_all->ProjectionX();
    		hz_npe_all1d->Draw(); 
    		foutroot->cd();
           	hz_npe_all1d->Write();
    		cz_npe->Write();
    		cz_npe->SaveAs("znpeall.png");
    		cz_npe->Close();
		TCanvas* cz_npedecay = new TCanvas();
		cz_npedecay->cd();
	    	cz_npedecay->SetName("chz_npe_alldecay");
    		hz_npe_alldecay->Draw("colz"); hz_npe_alldecay->SetMarkerStyle(8);
		cz_npedecay->SetLogz();
		TCanvas* cmu_npedecay = new TCanvas();
		cmu_npedecay->Divide(2,1);
		cmu_npedecay->cd(1);
	    	cmu_npedecay->SetName("chmu_npe_alldecay");
    		hmu_npe_alldecay->Draw("colz"); hmu_npe_alldecay->SetMarkerStyle(8);
		cmu_npedecay->cd(2);
		hmu_ndecay=hmu_npe_alldecay->ProjectionX();
		hmu_ndecay->Draw();
    		foutroot->cd();
    		cz_npedecay->Write();
    		cz_npedecay->SaveAs("znpealldecay.png");
    		cmu_npedecay->Write();
    		cmu_npedecay->SaveAs("munpealldecay.png");
    		cz_npe->Close();
		TCanvas* cz_npeemu = new TCanvas();
		cz_npeemu->cd();
	    	cz_npeemu->SetName("chz_npe_allemu");
    		hz_npe_allemu->Draw("colz"); hz_npe_allemu->SetMarkerStyle(8);
		cz_npeemu->SetLogz();
    		foutroot->cd();
    		cz_npeemu->Write();
    		cz_npeemu->SaveAs("znpeallemu.png");
    		cz_npe->Close();

		TCanvas* crp_npe = new TCanvas();
		crp_npe->cd();
	    	crp_npe->SetName("chrp_npe_all");
    		hrp_npe_all->Draw("colz"); hrp_npe_all->SetMarkerStyle(8);
		crp_npe->SetLogz();
    		foutroot->cd();
    		crp_npe->Write();
    		crp_npe->SaveAs("rpnpeall.png");
    		crp_npe->Close();
            hMuAAngletrigger->Write();
            hMuPAngletrigger->Write();
	for(int jj=0;jj<6;jj++){
	    hhodoXYtrigger[jj]->Write();
	}

    std::cout<<"Done reading files, now writting to ROOT file ..."<<std::endl;

    cout<<"Trigger type: \tnumber of events\n";
	int tmpcount=0;
	int tmpcounth4=0;
	int tmpcounth5=0;
    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
      cout<<strig[trigindex]<<"\t"<<trigcnts[trigindex]<<endl;
      gNumberOfEvents->SetPoint(trigindex,std::strtol(strig[trigindex].c_str(), 0, 2),trigcnts[trigindex]);
    }
	txtout<<argv[1]<<" ";
    for(int trigindex=0; trigindex<nbTrigs; trigindex++){
	tmpcount+=trigcnts[trigindex];
	if(trigindex%3==1)tmpcounth4+=trigcnts[trigindex];
	if(trigindex%3==2)tmpcounth5+=trigcnts[trigindex];
	//if(trigindex%3==2) {std::cout<<tmpcount<<" ";tmpcount=0;}
	std::cout<<tmpcount<<" ";
	txtout<<tmpcount<<" ";
	tmpcount=0;
    }
	txtout<<endl;
	txtout.close();
    std::cout<<tmpcounth4<<" "<<tmpcounth5<<endl;
    gNumberOfEvents->GetXaxis()->SetTitle("TrigType (converted from binary values)");
    gNumberOfEvents->GetYaxis()->SetTitle("Number Of Events");
    gNumberOfEvents->SetName("gNumberOfEvents_allTrigCombinations");
    TCanvas* cNumberOfEvents = new TCanvas();
    cNumberOfEvents->SetName("cNumberOfEvents");
    gNumberOfEvents->Draw("AP"); gNumberOfEvents->SetMarkerStyle(8);
    foutroot->cd();
    cNumberOfEvents->Write();
    cNumberOfEvents->Close();
	
    cout<<"HodoTrig "<<trigcnts_hodo<<endl;
    cout<<"Hodo + H5 + notH4 "<<trigcnts_hodo_H5_notH4<<endl;
    cout<<"Hodo + notH5 + notH4 "<<trigcnts_hodo_notH5_notH4<<endl;
    cout<<"Hodo + H5 + H4 "<<trigcnts_hodo_H5_H4<<endl;
    cout<<"Hodo + notH5 + H4 "<<trigcnts_hodo_notH5_H4<<endl;
    cout<<" - in which number of decays "<<trigcnts_decays<<endl;


    cout<<"total_number_of_events = "<<total_number_of_events<<endl;

    cout<<"testCntType1 ---- "<<endl;
    for(int t=0;t<8;t++) cout<<testcntTtrigtype1[t]<<"\t";
    cout<<endl;
    outrootnpe->cd();
    npetree->Write();
    oltree->Write();
    npedecaytree->Write();
//    foutroot->Write();
    foutroot->Close();
    outrootnpe->Close();
//rong zhao plot
    
    return 0; 

}


double getratio(double thick,double att,double thetai,double ni,double nr){
//cout<<TMath::Sin(thetai)<<endl;
double sinthetar=ni*TMath::Sin(thetai)/nr;
if(sinthetar>=1) return 0;
double costhetar=TMath::Sqrt((1-sinthetar*sinthetar));
double lprime=thick/costhetar;
double ratio=0;
//cout<<"thick diff:"<<lprime<<" "<<thick<<endl;
ratio=exp(-(lprime-thick)/att);
return ratio;
}
double getangle(double thetai,double ni,double nr){

double sinthetar=ni*TMath::Sin(thetai)/nr;
double angle=TMath::ASin(sinthetar);
//cout<<"thetarp:"<<angle<<endl;
return angle;
}

