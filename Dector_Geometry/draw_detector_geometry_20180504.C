{
  TCanvas* c0 = new TCanvas("1ton_geo_xy","1ton_geo_xy",600,600);
  c0->Draw();
// TFrame *frame1 = (TFrame)c0->GetFrame();
//frame1->SetLineWidth(2);
  gPad->SetTopMargin(0.15); gPad->SetLeftMargin(0.15);
  gPad->SetLineWidth(2);
  c0->SetFrameFillColor(kYellow-9);
  c0->SetFrameLineWidth(2);
  c0->SetFrameFillStyle(3003);
  TGraph* gXYregion = new TGraph();
  gXYregion->SetPoint(0,-500,-500);
  gXYregion->SetPoint(1,500,500);
  gXYregion->Draw("AP");
  gXYregion->GetXaxis()->SetTitle("X (mm)");
  gXYregion->GetYaxis()->SetTitle("Y (mm)");
  //draw the CD: its radius is 497.5
  TGraph* gCD = new TGraph();
  double r=497.5;//mm, radius
  int n000=0;
  for(double a=0; a<=TMath::Pi()*2.0; a=a+0.01){
    double thisx=r*TMath::Cos(a);
    double thisy=r*TMath::Sin(a);
    gCD->SetPoint(n000,thisx,thisy);
    n000++;
  }
  cout<<"n000="<<n000<<endl;
  gCD->Draw("L"); gCD->SetLineColor(2); gCD->SetLineWidth(2);
  //draw the teflon wall, its radius is 481.825 to 485
  TGraph* gteflon1 = new TGraph();
  r=481.825;//mm, radius
  n000=0;
  for(double a=0; a<=TMath::Pi()*2.0; a=a+0.01){
    double thisx=r*TMath::Cos(a);
    double thisy=r*TMath::Sin(a);
    gteflon1->SetPoint(n000,thisx,thisy);
    n000++;
  }
  cout<<"n000="<<n000<<endl;
  gteflon1->Draw("L"); gteflon1->SetLineColor(1); gteflon1->SetLineWidth(2);
  //
  TGraph* gteflon2 = new TGraph();
  r=485;//mm, radius
  n000=0;
  for(double a=0; a<=TMath::Pi()*2.0; a=a+0.01){
    double thisx=r*TMath::Cos(a);
    double thisy=r*TMath::Sin(a);
    gteflon2->SetPoint(n000,thisx,thisy);
    n000++;
  }
  cout<<"n000="<<n000<<endl;
  gteflon2->Draw("L"); gteflon2->SetLineColor(1); gteflon2->SetLineWidth(2);

  //pmt positions
  double x[8] = {
                381.35, 131.35, -118.65, -368.5, -165.025, -165.025, 381.35, 131.35
               };
  double y[8] = {
                -30.79375, -30.79375, -30.79375, -30.79375, 280.375, -128.425, 169.20625, 169.20625
               };
  double pmtz_bottom = -500.4;
  double pmtz_top = 800.4;
  TEllipse* pmts[8];
  for(int i=0;i<8;i++){
    pmts[i] = new TEllipse(x[i],y[i],25.4,25.4);
    pmts[i]->Draw("same"); pmts[i]->SetLineWidth(3);
  }

  //hodo counters x, y and z
  double hodox[6]={266.7, 267.05, 319.0875, 273.05, -174.625, 179.3875};
  double hodoy[6]={24.41575, 24.47925, 26.19375, 94.26575, -141.7125, -131.7625};
  double hodoz[6]={1133.5375, 889.7006, 1144.65, 896.5435, -1123.8875, -1123.8875};
  double hodosizex[6]={50.8, 50.8, 50.8, 50.8, 153.9875, 152.4};
  double hodosizey[6]={44.45, 57.15, 50.8, 57.15, 358.7125, 357.98125};
  double hodosizez[6]={3.175, 4.7244, 3.175, 4.3815, 5.08, 5.08};
  //H0
  double x000=266.7, y000=24.41575, sizex=50.8, sizey=44.45;
  TLine* h01=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h01->Draw("same"); h01->SetLineWidth(3);
  TLine* h02=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h02->Draw("same"); h02->SetLineWidth(3);
  TLine* h03=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h03->Draw("same"); h03->SetLineWidth(3);
  TLine* h04=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h04->Draw("same"); h04->SetLineWidth(3);
  //H1
  x000=267.05, y000=24.47925, sizex=50.8, sizey=57.15;
  TLine* h11=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h11->Draw("same"); h11->SetLineColor(2); h11->SetLineWidth(3);
  TLine* h12=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h12->Draw("same"); h12->SetLineColor(2); h12->SetLineWidth(3);
  TLine* h13=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h13->Draw("same"); h13->SetLineColor(2); h13->SetLineWidth(3);
  TLine* h14=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h14->Draw("same"); h14->SetLineColor(2); h14->SetLineWidth(3);
  //H2
  x000=319.0875, y000=26.19375, sizex=50.8, sizey=50.8;
  TLine* h21=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h21->Draw("same"); h21->SetLineColor(3); h21->SetLineWidth(3);
  TLine* h22=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h22->Draw("same"); h22->SetLineColor(3); h22->SetLineWidth(3);
  TLine* h23=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h23->Draw("same"); h23->SetLineColor(3); h23->SetLineWidth(3);
  TLine* h24=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h24->Draw("same"); h24->SetLineColor(3); h24->SetLineWidth(3);
  //H3
  x000=273.05, y000=94.26575, sizex=50.8, sizey=57.15;
  TLine* h31=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h31->Draw("same"); h31->SetLineColor(4); h31->SetLineWidth(3);
  TLine* h32=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h32->Draw("same"); h32->SetLineColor(4); h32->SetLineWidth(3);
  TLine* h33=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h33->Draw("same"); h33->SetLineColor(4); h33->SetLineWidth(3);
  TLine* h34=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h34->Draw("same"); h34->SetLineColor(4); h34->SetLineWidth(3);
  //H4
  x000=-174.625, y000=-141.7125, sizex=153.9875, sizey=358.40625;
  TLine* h41=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h41->Draw("same"); h41->SetLineColor(7); h41->SetLineWidth(3);
  TLine* h42=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h42->Draw("same"); h42->SetLineColor(7); h42->SetLineWidth(3);
  TLine* h43=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h43->Draw("same"); h43->SetLineColor(7); h43->SetLineWidth(3);
  TLine* h44=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h44->Draw("same"); h44->SetLineColor(7); h44->SetLineWidth(3);
  //H5
  x000=179.3875, y000=-131.7625, sizex=152.4, sizey=357.98125;
  TLine* h51=new TLine(x000-sizex,y000-sizey,x000+sizex,y000-sizey); h51->Draw("same"); h51->SetLineColor(6); h51->SetLineWidth(3);
  TLine* h52=new TLine(x000+sizex,y000-sizey,x000+sizex,y000+sizey); h52->Draw("same"); h52->SetLineColor(6); h52->SetLineWidth(3);
  TLine* h53=new TLine(x000-sizex,y000+sizey,x000+sizex,y000+sizey); h53->Draw("same"); h53->SetLineColor(6); h53->SetLineWidth(3);
  TLine* h54=new TLine(x000-sizex,y000-sizey,x000-sizex,y000+sizey); h54->Draw("same"); h54->SetLineColor(6); h54->SetLineWidth(3);

  //put legend
  TLegend* leg=new TLegend(0.15,0.85,0.90,0.95);
  leg->SetNColumns(3);
  leg->AddEntry(h01,"Hodoscope0","lf");
  leg->AddEntry(h11,"Hodoscope1","lf");
  leg->AddEntry(h21,"Hodoscope2","lf");
  leg->AddEntry(h31,"Hodoscope3","lf");
  leg->AddEntry(h41,"Hodoscope4","lf");
  leg->AddEntry(h51,"Hodoscope5","lf");
//  leg->SetTextColor(kBlue);
  leg->Draw();
  leg->SetLineWidth(2);

   TLatex latex;
   latex.SetTextSize(0.05);
   latex.SetTextAlign(13);  //align at top
   latex.DrawLatex(380,-70,"S_{0}");
   latex.DrawLatex(130,-70,"S_{1}");
   latex.DrawLatex(-110,-70,"S_{2}");
   latex.DrawLatex(-410,-70,"S_{3}");
   latex.DrawLatex(-200,-180,"S_{5}");
   latex.DrawLatex(-130,310,"S_{4}");
   latex.DrawLatex(130,320,"S_{7}");
   latex.DrawLatex(300,310,"S_{6}");

  c0->Print("geoxy.eps");

  // now draw the geometry in x-z plane
  TCanvas* c1 = new TCanvas("1ton_geo_xz","",600,600);
  c1->cd();
  gPad->SetTopMargin(0.15); gPad->SetLeftMargin(0.15);
  gPad->SetLineWidth(2);
  c1->SetFrameFillColor(kYellow-9);
  c1->SetFrameLineWidth(2);
  c1->SetFrameFillStyle(3003);
  // acrylic wall - outer
  TGraph* gAcrylicWall1 = new TGraph();
  gAcrylicWall1->SetPoint(0,-522.9,150-625); 
  gAcrylicWall1->SetPoint(1,522.9, 150-625);
  gAcrylicWall1->SetPoint(2,522.9, 625.0+150);
  gAcrylicWall1->SetPoint(3,-522.9,625.0+150);
  gAcrylicWall1->SetPoint(4,-522.9,150-625);
  gAcrylicWall1->Draw("AL");
  gAcrylicWall1->GetXaxis()->SetTitle("X (mm)");
  gAcrylicWall1->GetYaxis()->SetTitle("Z (mm)");
  gAcrylicWall1->SetLineWidth(2);
//  gAcrylicWall1->GetXaxis()->SetTitle("X (mm)");
  gAcrylicWall1->GetYaxis()->SetRangeUser(-1250,1250);
  // acrylic wall - inner, in fact it is the CD
  TGraph* gAcrylicWall2 = new TGraph();
  gAcrylicWall2->SetPoint(0,-497.5,150-625); 
  gAcrylicWall2->SetPoint(1,497.5, 150-625);
  gAcrylicWall2->SetPoint(2,497.5, 625.0+150);
  gAcrylicWall2->SetPoint(3,-497.5,625.0+150);
  gAcrylicWall2->SetPoint(4,-497.5,150-625);
  gAcrylicWall2->Draw("L");
  gAcrylicWall2->SetLineColor(2);
  gAcrylicWall2->SetLineWidth(2);
  //teflon wall outer
  TGraph* gTeflonWall1 = new TGraph();
  gTeflonWall1->SetPoint(0,-485,150-625); 
  gTeflonWall1->SetPoint(1,485, 150-625);
  gTeflonWall1->SetPoint(2,485, 625.0+150);
  gTeflonWall1->SetPoint(3,-485,625.0+150);
  gTeflonWall1->SetPoint(4,-485,150-625);
  gTeflonWall1->Draw("L");
  gTeflonWall1->SetLineWidth(2);
  //teflon wall inner
  TGraph* gTeflonWall2 = new TGraph();
  gTeflonWall2->SetPoint(0,-481.825,150-625); 
  gTeflonWall2->SetPoint(1,481.825, 150-625);
  gTeflonWall2->SetPoint(2,481.825, 625.0+150);
  gTeflonWall2->SetPoint(3,-481.825,625.0+150);
  gTeflonWall2->SetPoint(4,-481.825,150-625);
  gTeflonWall2->Draw("L");
  gTeflonWall2->SetLineWidth(2);
  //acrylic bottom
  TGraph* gAcrylicBottom = new TGraph();
  gAcrylicBottom->SetPoint(0,-522.9, -487.7-12.7);
  gAcrylicBottom->SetPoint(1,522.9, -487.7-12.7);
  gAcrylicBottom->SetPoint(2,522.9, -487.7+12.7);
  gAcrylicBottom->SetPoint(3,-522.9, -487.7+12.7);
  gAcrylicBottom->SetPoint(4,-522.9, -487.7-12.7);
  gAcrylicBottom->Draw("L");
  gAcrylicBottom->SetLineWidth(2);
  //acrylic lid
  TGraph* gAcrylicLid = new TGraph();
  gAcrylicLid->SetPoint(0,-575, 787.7-12.7);
  gAcrylicLid->SetPoint(1,575, 787.7-12.7);
  gAcrylicLid->SetPoint(2,575, 787.7+12.7);
  gAcrylicLid->SetPoint(3,-575, 787.7+12.7);
  gAcrylicLid->SetPoint(4,-575, 787.7-12.7);
  gAcrylicLid->Draw("L");
  gAcrylicLid->SetLineWidth(2);
  //pmts
  TGraph* gPMT_XZ[8];
  for(int i=0;i<8;i++){
    gPMT_XZ[i] = new TGraph();
    if(i!=4 && i!=5){//bottom pmts
      gPMT_XZ[i]->SetPoint(0,x[i]-25.4,pmtz_bottom-112);
      gPMT_XZ[i]->SetPoint(1,x[i]+25.4,pmtz_bottom-112);
      gPMT_XZ[i]->SetPoint(2,x[i]+25.4,pmtz_bottom);
      gPMT_XZ[i]->SetPoint(3,x[i]-25.4,pmtz_bottom);
      gPMT_XZ[i]->SetPoint(4,x[i]-25.4,pmtz_bottom-112);
      gPMT_XZ[i]->Draw("L");
    }
    else{
      gPMT_XZ[i]->SetPoint(0,x[i]-25.4,pmtz_top);
      gPMT_XZ[i]->SetPoint(1,x[i]+25.4,pmtz_top);
      gPMT_XZ[i]->SetPoint(2,x[i]+25.4,pmtz_top+112);
      gPMT_XZ[i]->SetPoint(3,x[i]-25.4,pmtz_top+112);
      gPMT_XZ[i]->SetPoint(4,x[i]-25.4,pmtz_top);
      gPMT_XZ[i]->Draw("L");
    }
    gPMT_XZ[i]->SetLineColor(1);
    gPMT_XZ[i]->SetLineWidth(2);
  }
  //hodo detectors
  TGraph* gHodo_XZ[6];
  TLegend* legHodoXZ = new TLegend(0.15,0.85,0.90,0.95);
  legHodoXZ->SetLineWidth(2);
  legHodoXZ->SetNColumns(3);
  string strhodoname[6]={"Hodoscope0","Hodoscope1","Hodoscope2","Hodoscope3","Hodoscope4","Hodoscope5"};
  int hodocolor[6]={1,2,3,4,7,6};
  for(int i=0;i<6;i++){
    gHodo_XZ[i] = new TGraph();
    gHodo_XZ[i]->SetPoint(0,hodox[i]-hodosizex[i],hodoz[i]-hodosizez[i]);
    gHodo_XZ[i]->SetPoint(1,hodox[i]+hodosizex[i],hodoz[i]-hodosizez[i]);
    gHodo_XZ[i]->SetPoint(2,hodox[i]+hodosizex[i],hodoz[i]+hodosizez[i]);
    gHodo_XZ[i]->SetPoint(3,hodox[i]-hodosizex[i],hodoz[i]+hodosizez[i]);
    gHodo_XZ[i]->SetPoint(4,hodox[i]-hodosizex[i],hodoz[i]-hodosizez[i]);
    gHodo_XZ[i]->Draw("L");
    gHodo_XZ[i]->SetLineColor(hodocolor[i]);
    gHodo_XZ[i]->SetLineWidth(2);
    legHodoXZ->AddEntry(gHodo_XZ[i],strhodoname[i].c_str(),"l");
  }
  legHodoXZ->Draw();
  c1->Print("geoxz.eps");

  // now draw the geometry in y-z plane
  TCanvas* c2 = new TCanvas("1ton_geo_yz","",600,600);
  c2->cd();
  gPad->SetTopMargin(0.15); gPad->SetLeftMargin(0.15);
  gPad->SetLineWidth(2);
  c2->SetFrameFillColor(kYellow-9);
  c2->SetFrameLineWidth(2);
  c2->SetFrameFillStyle(3003);
  gAcrylicWall1->Draw("AL");
  gAcrylicWall1->GetXaxis()->SetTitle("Y (mm)");
  gAcrylicWall1->GetYaxis()->SetTitle("Z (mm)");

  gAcrylicWall2->Draw("L");
  gTeflonWall1->Draw("L");
  gTeflonWall2->Draw("L");
  gAcrylicBottom->Draw("L");
  gAcrylicLid->Draw("L");

  //pmts
  TGraph* gPMT_YZ[8];
  for(int i=0;i<8;i++){
    gPMT_YZ[i] = new TGraph();
    if(i!=4 && i!=5){//bottom pmts
      gPMT_YZ[i]->SetPoint(0,y[i]-25.4,pmtz_bottom-112);
      gPMT_YZ[i]->SetPoint(1,y[i]+25.4,pmtz_bottom-112);
      gPMT_YZ[i]->SetPoint(2,y[i]+25.4,pmtz_bottom);
      gPMT_YZ[i]->SetPoint(3,y[i]-25.4,pmtz_bottom);
      gPMT_YZ[i]->SetPoint(4,y[i]-25.4,pmtz_bottom-112);
      gPMT_YZ[i]->Draw("L");
    }
    else{
      gPMT_YZ[i]->SetPoint(0,y[i]-25.4,pmtz_top);
      gPMT_YZ[i]->SetPoint(1,y[i]+25.4,pmtz_top);
      gPMT_YZ[i]->SetPoint(2,y[i]+25.4,pmtz_top+112);
      gPMT_YZ[i]->SetPoint(3,y[i]-25.4,pmtz_top+112);
      gPMT_YZ[i]->SetPoint(4,y[i]-25.4,pmtz_top);
      gPMT_YZ[i]->Draw("L");
    }
    gPMT_YZ[i]->SetLineColor(1);
    gPMT_YZ[i]->SetLineWidth(2);
  }
  //hodo detectors
  TGraph* gHodo_YZ[6];
  TLegend* legHodoYZ = new TLegend(0.15,0.85,0.9,0.95);
  legHodoYZ->SetNColumns(3);
  legHodoYZ->SetLineWidth(2);
  for(int i=0;i<6;i++){
    gHodo_YZ[i] = new TGraph();
    gHodo_YZ[i]->SetPoint(0,hodoy[i]-hodosizey[i],hodoz[i]-hodosizez[i]);
    gHodo_YZ[i]->SetPoint(1,hodoy[i]+hodosizey[i],hodoz[i]-hodosizez[i]);
    gHodo_YZ[i]->SetPoint(2,hodoy[i]+hodosizey[i],hodoz[i]+hodosizez[i]);
    gHodo_YZ[i]->SetPoint(3,hodoy[i]-hodosizey[i],hodoz[i]+hodosizez[i]);
    gHodo_YZ[i]->SetPoint(4,hodoy[i]-hodosizey[i],hodoz[i]-hodosizez[i]);
    gHodo_YZ[i]->Draw("L");
    gHodo_YZ[i]->SetLineColor(hodocolor[i]);
    gHodo_YZ[i]->SetLineWidth(2);
    legHodoYZ->AddEntry(gHodo_YZ[i],strhodoname[i].c_str(),"l");
  }
  legHodoYZ->Draw();
  c2->Print("geoyz.eps");

//  TCanvas* c3 = new TCanvas("1ton_geo_XY","",600,600);  
  
  
}
