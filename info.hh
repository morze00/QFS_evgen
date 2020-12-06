#ifndef _INFO_HH_
#define _INFO_HH_

const int MAX_STORY=1000000;		// Total number of generated events
const double ENERGY=500.00;		// 12C Beam energy (MeV/u) from Atima
const double AMU=931.4940038;		// 12C Beam energy (MeV/u) from Atima
const int A=17;			// Mass number of the incoming nucleus

const double MASS_A = 17.0177140*AMU; 	// Nuclear mass of 17Ne (MeV) from AME2016 
const double MASS_B = 16.011466 * AMU; 	// Nuclear mass of 16F from AME2016
//const double MASS_A= 11177.900; 	// Nuclear mass of 12C (MeV) from http://t2.lanl.gov/data/astro/molnix96/massd.html
//const double MASS_B= 10255.100; 	// Nuclear mass of 11B (MeV) from http://t2.lanl.gov/data/astro/molnix96/massd.html
const double NEUTRON_MASS= 939.565;	// MeV
const double PROTON_MASS= 938.279;	// MeV	 
const double ALPHA_MASS= 3728.401;	// MeV	 
const double EXE= 0.0;			// Residual excitation energy (MeV) in 10B (change it e.g. for deeply bound states)
const double MOM_SIGMA= 106.5;	 	// Internal momentum spread (Gauss distribution)
//const double MOM_SIGMA= 0.001;	 	// Internal momentum spread (Gauss distribution)

const bool DOPPLER= true;

//Constants
const double PI = 3.141592653589793238;
const double R2D = 57.29577951;
const double D2R = 0.017453292;
const double HBAR = 197.326960;
const double ALPHA = 0.007297352;
const double LN2 = 0.6931471805599452862;// ln2
const double CC = 299792458.;		// Speed of light
const double UNIT = 931.494013; 	// Atomic mass unit  (MeV/c^2)

#endif

