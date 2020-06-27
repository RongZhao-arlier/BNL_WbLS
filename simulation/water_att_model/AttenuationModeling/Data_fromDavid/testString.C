{
  string a = "
{
name: \"OPTICS\",
index: \"water\",
valid_begin : [0, 0],
valid_end : [0, 0],
\/\/LIGHT_YIELD: 0.0
NEUTRON_CAPTURE_TIME_value1: [0.0, 1.0, ],
NEUTRON_CAPTURE_TIME_value2: [163000.0, 163000.0, ],
NEUTRON_SLOW_DIFFUSION_CONST_value1: [0.0, 1.0, ],
NEUTRON_SLOW_DIFFUSION_CONST_value2: [0.03, 0.03, ],
NEUTRON_FAST_DIFFUSION_RMS_value1: [0.0, 1.0, ], NEUTRON_FAST_DIFFUSION_RMS_value2: [50.0, 50.0, ],
RINDEX_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water 
RINDEX_value1: [60.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, ],
RINDEX_value2: [1.42516, 1.396, 1.37761, 1.35942, 1.34978, 1.34378, 1.3397, 1.33676, 1.33458, 1.33293, 1.33165, 1.33065, 1.32986, 1.3292, ],
ABSLENGTH_option: \"wavelength\", \/\/60 nm is bogus to prevent G4 from complaining when VUVs hit the water 
";

  cout<<a<<endl;

  double b = 100.01;
  ostringstream s;
  s << b;
  string c = "ABSLENGTH_value1: [";
  c = c + s.str();
  cout<<c<<endl;
}
