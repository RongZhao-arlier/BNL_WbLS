{
  fstream fout("temp_insert_val.txt",ios::out);
  double w = 200;
  while(w<=380){
    fout<<w<<"\t"<<"1.0e-20"<<endl;
	w+=1.0;
  }
  fout.close();
  exit(0);
}