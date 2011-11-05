double RDIFF_CUT = 27.;//cm
double THETADIFF_CUT = 0.4;//rad
double FIDX_CUT = 3.;//cm
double FIDY_CUT = 4.;//cm
double FIDZup_CUT = 6.;//cm
double FIDZdown_CUT = 4.;//cm

// double RDIFF_CUT = 27.;//cm
// double THETADIFF_CUT = 0.4;//rad
// double FIDX_CUT = 4.;//cm
// double FIDY_CUT = 5.;//cm
// double FIDZup_CUT = 7.;//cm
// double FIDZdown_CUT = 5.;//cm


double delta_X = .021;   //fractional //1.56+/-.03 mm/us
double delta_Y = .0025;   //fractional corresponding to 5 mm
double delta_Z = .0011 ; //fractional corresponding to 5 mm
double NTARGETS = ((47.-(2.*FIDX_CUT))*(40.-(2.*FIDY_CUT))*(90.-FIDZup_CUT-FIDZdown_CUT))*1.399*(1./39.948)*(6.0221415*pow(10,23));//number of argon atoms (not number of neutrons+protons)

double delta_NTARGETS = sqrt(pow(delta_X,2)+pow(delta_Y,2)+pow(delta_Z,2));//fractional

// double INTEGRATEDFLUX_3_50 = 133010./10000.; //particles/gev/cm^2/10^9 POT 
// double INTEGRATEDFLUX_0_3 = 42145./10000.; //particles/gev/cm^2/10^9 POT

double INTEGRATEDFLUX_3_50 = 164588./10000.; //particles/cm^2/10^9 POT 
double INTEGRATEDFLUX_0_3 = 126435./10000.; //particles/cm^2/10^9 POT

double INTEGRATEDFLUX_0_50 = INTEGRATEDFLUX_0_3 + INTEGRATEDFLUX_3_50;

double delta_INTEGRATEDFLUX_3_50 = .067; //fractional //8948.5/133010//correlated errors
double delta_INTEGRATEDFLUX_0_3 = .3500; //fractional  //14751/42145

double delta_INTEGRATEDFLUX_0_50 = .157;//0_3 and 3_50 are considered uncorrelated

double denom=INTEGRATEDFLUX_0_50*(8.518840*pow(10,9))*NTARGETS;//POT

double delta_pot=.01;

double MOM_OFFSET_UNC=0.;

double MOM_OFFSET;