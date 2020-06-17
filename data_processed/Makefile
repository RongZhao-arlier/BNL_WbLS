ROOTSYS = /gpfs01/lbne/users/aw325/sw/root5.34.36.install
#CXX = $(shell $(ROOTSYS)/bin/root-config --cxx)
#CXX = /opt/rh/devtoolset-3/root/usr/bin/c++
CXX =/usr/bin/g++

all:
	
%:%.cc
	$(CXX) -o $@ $< -g `root-config --cflags` -I$(ROOTSYS)/include -L$(ROOTSYS)/lib -Llib -D__ROOTSHAREDLIBRARY__  -fPIC -lMinuit -lTMVA -lXMLIO -lMLP -lTreePlayer $(shell $(ROOTSYS)/bin/root-config --libs) 
