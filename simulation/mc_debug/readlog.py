import ROOT
import sys

h_posz = ROOT.TH1D("end_z","final z of photon" ,200 ,-650.,801 )
h_zout = ROOT.TH1D("z_out","final z of photon out detector" ,200 ,-1500.,1500 )
h_zlid = ROOT.TH1D("z_lid","final z in acrylic lid" ,25 ,775.,801 )
h_zcd = ROOT.TH1D("z_cd","final z in cd" ,200 ,-475.,775 )
h_zacryb = ROOT.TH1D("z_acry","final z in bottom acrylic" ,25 ,-500.4,-475 )
h_zoc = ROOT.TH1D("z_oc","final z in optical cookie" ,10 ,-502.4,-500.4 )
h_zpmt = ROOT.TH1D("z_pmt","final z in pmt" ,50 ,-614.4,-502.4 )
with open("tmp.txt","r") as f:
	for line in f.readlines():
		line=line.strip('\n')
		linelist=line.split()
		#print(linelist[3])
		h_posz.Fill(float(linelist[3]))
		z=float(linelist[3])
		if z>=775 and z<800.4: # the acrylic_lid
			h_zlid.Fill(float(linelist[3]))
		if z>=-475.5 and z<775: # the cd
			h_zcd.Fill(float(linelist[3]))
		if z>=-500.4 and z<-475.5: # the acryb
			h_zacryb.Fill(float(linelist[3]))
		if z>=-502.4 and z<-500.4: # the cptical cookie
			h_zoc.Fill(float(linelist[3]))
		if z>=-614.4 and z<-502.4: # the pmt
			h_zpmt.Fill(float(linelist[3]))
		if z>=800.4 or z<-614.4: # the outside
			h_zout.Fill(float(linelist[3]))
h_posz.SetDirectory(0)
h_zlid.SetDirectory(0)
h_zcd.SetDirectory(0)
h_zacryb.SetDirectory(0)
h_zoc.SetDirectory(0)
h_zpmt.SetDirectory(0)
h_zout.SetDirectory(0)
outHistFile=ROOT.TFile.Open("mc_op_test.root","RECREATE")
outHistFile.cd()
h_posz.Write()
h_zlid.Write()
h_zcd.Write()
h_zacryb.Write()
h_zoc.Write()
h_zpmt.Write()
h_zout.Write()
outHistFile.Close()
