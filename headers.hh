#ifndef __HEADERS_HH__
#define __HEADERS_HH__

//A structure to store momentum and angles of two
//scattered particles in CM system
struct cm_values
{
	//Internal cluster
	double	e_clust;
	double	p_clust;
	double	theta_clust;
	//Scattered particle
	double	e_scat;
	double	p_scat;
	double	theta_scat;
	//indicates satisfactory kinematics (i.e. energy & momentum conservation)
	bool 	good;
	double T;
};

void run(Char_t * root_filename,
        Char_t * ascii_filename,
        Char_t * input_filename, 
        char   * theory_mom_filename,
        int nevents,
        Bool_t is_ppn,
        Bool_t is_ppa, 
        Bool_t is_isotropic,
        Bool_t is_direct,
        Bool_t is_nnp
        );
TVector3 DREHUNG(TVector3 v1,TVector3 v2);
double CINEMA(double,double,double);
cm_values CENMASS(double,double,double,double);
double momentum_CM(double,double,double);
std::pair<double, double> Lorentz(double,double,double,double);
double get_T(double,double);
TH1F* Theoretical_Momentum(char* filename);

#endif //__HEADERS_HH__
