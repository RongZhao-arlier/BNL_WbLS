{
  double x1=323.850006, y1=0.000212;
  double x2=380, y2=0.0001137;
  
  double slope = (y2-y1)/(x2-x1);
  double inter = y2-slope*x2;
  
  fstream fout("interpolation_data.txt",ios::out);
  double x= 324;
  double y = 0;
  while(x<380){
	  y = slope*x + inter;
	  fout<<x<<"\t"<<y<<endl;
	  x += 1.0; 
  }
  fout.close();
}
