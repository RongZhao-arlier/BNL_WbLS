{
  fstream fin("testdata.txt",ios::in);
  string s;
  getline(fin,s);
  fin.close();
  double val;

  std::string delimiter = ", ";

  size_t pos = 0;
  std::string token;
  fstream fout("testdata_out.txt",ios::out);
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
	fout<<0.01*atof(token.c_str())<<", ";
  }
  fout<<endl;
  fout.close();
  //std::cout << s << std::endl;  
  exit(0);
}