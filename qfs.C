/* Original code was writen in FORTRAN by Leonid Chulkov l.chulkov@gsi.de
 *
 * The program is simulating inverse kinematics of a quasi-elastic scattering 
 * of a proton "i" (being at rest) on a cluster/nucleon "a" inside an incident
 * nucleus A. The kinematics of the residual system "B" and the two scattered particles
 * 1 and 2 is stored in the output tree.
 *
 *  A --> i  ==> B + (i -> a) =>  B + 1 + 2
 *
 * Valerii Panin October 2012 (v.panin@gsi.de)
 *
 * 27/10/2020 -> Modified as following for direct kinematics (n,np) reactions for Vadim
 *
 * With the direct kinematics option the program is simulating the following:
 * The beam particle "i" scatters off a cluster/nucleon "a" inside a target
 * nucleus A. The kinematical variables of the residual fragment "B" and the
 * two scattered particles "1" (scattered) and "2" (knocked-out) are stored 
 * in the output ROOT tree.
 *
 *  i --> A  ==> (i -> a) + B  =>  1 + 2 + B
 *
 * Valerii Panin October 2020 (v.panin@gsi.de)
 */


#include "libs.hh"
#include "headers.hh"
#include "info.hh"

using namespace std;


Double_t *a;
Double_t *b;
Double_t *c;
#define NP 61

Double_t slp[NP]={0};
Double_t offs[NP]={0};

Double_t funct(Double_t *x, Double_t *par){

    Double_t dxs = -1;
    Double_t dxssin = -1;

    if(x[0]>a[NP])
        dxs = (c[NP]+c[NP-1])/2.;
    else if(x[0]<a[0])
        dxs = c[0];

    for(Int_t i=0; i<NP; i++){
        if(x[0]>a[i] && x[0]<a[i+1]){
            dxs = slp[i]*x[0]+offs[i];
            //dxs = (c[i]+c[i+1])/2.;    
            break;
        }
    }

    dxssin = (dxs*2*PI*sin(x[0]/180.*PI));
    //if(x[0]<16.||x[0]>164.) { dxssin = 0.; }
    //return dxs;
    return dxssin;
}


TF1* f1;
bool isotropic = false;

//Main function
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
        )
{
    TTree *tr = new TTree("tr","data");

    cout << "\n-- Starting main function!";

    //tr->ReadFile("data_filtered_calc_385MeV.txt","dxs/F:erdxs/F:energy/F:angle/F"); // Choose Parametrization File
    tr->ReadFile("data_filtered_calc_400MeV.txt","dxs/F:erdxs/F:energy/F:angle/F"); // Choose Parametrization File

    tr->Draw("dxs:erdxs:angle","","goff");
    a = tr->GetV3();
    b = tr->GetV2();
    c = tr->GetV1();
    Double_t aerr[NP]={0};

    for(Int_t i=0; i<tr->GetSelectedRows(); i++){
        //cout << a[i] << " " << b[i] << " " << c[i] << endl;
        slp[i]=(c[i+1]-c[i])/(a[i+1]-a[i]);
        offs[i]=c[i]-slp[i]*a[i];
        aerr[i]=0;
    }
    //Cross section function
    f1 = new TF1("f1",funct,0.,180.,0);

    //Theoretical momentum disribution
    TH1F * h_momentum_theory;
    if(theory_mom_filename)
        h_momentum_theory = Theoretical_Momentum(theory_mom_filename);


    //------- Output tree ----------
    TFile file(root_filename,"RECREATE");
    TTree *tree = new TTree("tree","Tree with simulated QFS kinematics");
#include "tree.hh"

    // Total kinetic energy (MeV) of the projectile
    double Tkin = ENERGY * A;
    if(is_direct) Tkin = ENERGY;

    const double MA   = MASS_A;  	
    const double MB   = MASS_B + EXE; 

    double Ma = PROTON_MASS;

    std::cout << "\n-- Max events = " << nevents;

    if(is_isotropic){
        isotropic = true;
        cout << "\n-- Using isotropic nucleon distribution!" << endl;
    }
    else{
        isotropic = false ;
        cout <<  "\n-- Using parametrized anisotropic nucleon distribution!";
    }

    if(!is_direct){
        cout << "\n-- Assuming inverse kinematics";
    }

    if(theory_mom_filename){
        cout << "\n-- Using theoretical momentum distribution of nucleons inside A";
    }

    if(is_ppn){
        Ma   = NEUTRON_MASS;  
        // 		isotropic = false;
        cout << "\n-- Generating (p,pn) events!" << endl;
        cout << "Neutron Mass: " << Ma << " MeV";
    }
    else if(is_ppa){
        Ma   = ALPHA_MASS;  
        cout << "\n-- Generating (p,pa) events!" << endl;
        cout << "Alpha Mass: " << Ma << " MeV";
    }
    else if(is_nnp){
        cout << "\n-- Generating (n,np) events!" << endl;
        cout << "-- Proton Mass: " << Ma << " MeV";
    }

    else{
        cout << "\n-- Generating (p,2p) events!" << endl;
        cout << "Proton Mass: " << Ma << " MeV";
    }

    double Mi   = PROTON_MASS; 		
    if(is_nnp) Mi = NEUTRON_MASS;

    double PA, EA, bA, gA, Pi, Ei, bi, gi, EaL, EBL;
    double S_first, sigma; 

    //--------------------- Beam parameters -----------------------------

    if(!is_direct)//inverse kinematics
    {
        PA = sqrt(Tkin*(Tkin + 2*MA)); 	// Total 3-momentum of the beam
        EA = sqrt(MA*MA + PA*PA);      	// Total energy of the beam
        bA = -PA/EA;		      	// Beta of the beam
        gA = 1/sqrt(1-bA*bA);	      	// Gamma of the beam
        S_first = (EA+Mi)*(EA+Mi) - PA*PA;	// Invariant mass (Mandelstam S-variable)
        sigma = MOM_SIGMA; 		// Internal momentum spread of the cluster "a" inside "A"

        cout << "\n****** Beam parameters ********";
        cout << "\nMA:\t" << MA << " MeV";
        cout << "\nTotal momentum:\t" << PA << " MeV";
        cout << "\nTotal energy:\t"   << EA << " MeV";
        cout << "\nBeta (beam):\t" << (-bA) << "\nGamma (beam):\t" << gA << "\n\n";
        cout << "\nProcessing " << nevents << " events........\n\n";
    }

    else//direct kinematics
    {
        Pi = sqrt(Tkin*(Tkin + 2*Mi)); 	// Total 3-momentum of the projectile
        Ei = sqrt(Mi*Mi + Pi*Pi);      	// Total energy of the projectile
        bi = -Pi/Ei;		      	// Beta of the projectile
        gi = 1/sqrt(1-bi*bi);	      	// Gamma of the projectile
        S_first = (Ei+MA)*(Ei+MA) - Pi*Pi;	// Invariant mass (Mandelstam S-variable)
        sigma = MOM_SIGMA; 		// Internal momentum spread of the cluster "a" inside "A"

        cout << "\n\n****** Beam parameters (direct kinematics) ********";
        cout << "\nMi:\t" << Mi << " MeV";
        cout << "\nTotal momentum:\t" << Pi << " MeV";
        cout << "\nTotal energy:\t"   << Ei << " MeV";
        cout << "\nBeta (beam):\t" << (-bi) << "\nGamma (beam):\t" << gi;
        cout << "\nMandelstam S:\t" << S_first << " MeV\n\n";
        cout << "\nProcessing " << nevents << " events........\n\n";
    }

    //Random number generators in ROOT 
    // 	gRandom = new TRandom3(); //default 
    // 	gRandom->SetSeed(0); //using computer time


    TRandom3 r1;
    r1.SetSeed(0);

    TH1F* input_xy;
    TH1F* input_z;

    //if(strcmp(input_filename,"gauss")==1)
    //{
    //    TFile* input_p = TFile::Open(input_filename,"READ");

    //    input_xy  = (TH1F*) input_p->Get("p_y")->Clone("input_xy");
    //    input_z   = (TH1F*) input_p->Get("p_z")->Clone("input_z");
    //}

    //Output text file for R3BROOT event generator
    ofstream outfile;
    outfile.open(ascii_filename);

    TLorentzVector LVa; //Lorentz vector of the cluster 
    TLorentzVector LVi; //Lorentz vector of the incident particle

    if(is_direct) LVi.SetPxPyPzE(0.0,0.0,Pi,Ei);
    else LVi.SetPxPyPzE(0.0,0.0,0.0,Mi);

    TVector3 Pa(1e-9,0.0,0.0);
    TVector3 P1cm(1e-9,0.0,0.0);
    TVector3 P2cm(1e-9,0.0,0.0);

    int events = 0;//generated event counter
    while(events<nevents) //eventloop
    {
        //------------ Internal momentum of a cluster -------------------
        //if(strcmp(input_filename,"gauss")==1)
        //{	
        //    Pa.SetX(input_xy->GetRandom());
        //    Pa.SetY(input_xy->GetRandom());
        //    Pa.SetZ(input_z->GetRandom());
        //}
        
        //Random sampling from theoretical momentum if filename is not NULL(e.g. main optin is given)
        if(theory_mom_filename)
        {
            Pa.SetMag(h_momentum_theory->GetRandom());
            Pa.SetTheta( TMath::ACos(2.* (r1.Uniform(0.,1.)) - 1.));
            Pa.SetPhi(2.*TMath::Pi()  * (r1.Uniform(0.,1.)) );

        }
        else
        {
            Pa.SetX(r1.Gaus(0,sigma));
            Pa.SetY(r1.Gaus(0,sigma));
            Pa.SetZ(r1.Gaus(0,sigma));
        }
        //------------ Internal momentum of the residual-----------------
        TVector3 PB;
        PB.SetX(-Pa.X());
        PB.SetY(-Pa.Y());
        PB.SetZ(-Pa.Z());

        //Filling tree variables for the fragmnent
        PBx    = PB.X();
        PBy    = PB.Y();
        PBz_rf = PB.Z();//Longitudinal momentum in the restframe

        //From the energy conservation in the virtual dissociation A->B+a
        double rrtt = MA*MA + MB*MB - 2*MA * sqrt(MB*MB + Pa.Mag2());
        if(rrtt<=0){
            cout<<"\nERROR off-shell mass!!";//non-zero and real off-shell mass
            cout<<"\nP:"<< PBx << "\t" << PBy << "\t" << PBz_rf << "\n";//non-zero and real off-shell mass
            continue;	
        }
        //Off-shell mass of the bound cluster

        double Ma_off = sqrt(rrtt);
        //Total energies of "a" and "B" in the restframe of "A"
        double EaL  = sqrt(Ma_off*Ma_off + Pa.Mag2()); 
        double EBL  = sqrt(MB*MB + PB.Mag2());

        if(!is_direct)
        {
            //------- Lorentz transformations into laboratory system ---------
            std::pair<double, double> lora = Lorentz(gA,bA,EaL,Pa.Z());
            EaL = lora.first; //cluster energy in lab
            Pa.SetZ(lora.second); 	 //cluster Pz in lab

            std::pair<double, double> lorB = Lorentz(gA,bA,EBL,PB.Z());
            EBL = lorB.first; //energy of the residual B in lab
            PB.SetZ(lorB.second);    //Pz of the residual B in lab
            //---------- Generating CM scattering process ----------
            //S = Ma_off*Ma_off + Mi*Mi + 2*Mi*EaL; //Mandelstam invariant
        }

        Ea = EaL - Ma_off;//filling tree variable

        LVa.SetPxPyPzE(Pa.X(), Pa.Y(), Pa.Z(), EaL);
        TLorentzVector LVstart = LVa + LVi;
        double S = LVstart.Mag2(); //Mandelstam invariant
        //cout << "\nS=" <<  LVstart.Mag2();

        Mandelstam_S = S;//filling tree variable
        //Now generate CM scattering kinematics
        cm_values CM = CENMASS(S,Ma_off,Mi,Ma);
        if(!CM.good) continue;//non-physical output

        double phi_rand = r1.Uniform(-PI,PI);

        P2cm.SetMag(CM.p_clust);
        P2cm.SetTheta(CM.theta_clust);
        P2cm.SetPhi(phi_rand);

        P1cm.SetX(-P2cm.X());
        P1cm.SetY(-P2cm.Y());
        P1cm.SetZ(-P2cm.Z());

        //------- Calculate realtive to the direction of the quasi-particle (cluster) --------
        double beta_cm;
        //if(is_direct){ beta_cm = 0.-LVstart. Beta(); } 
        //else{ beta_cm = -Pa.Mag()/(EaL+Mi); }
        beta_cm = 0.-LVstart.Beta();

        double gamma_cm = 1/sqrt(1-beta_cm*beta_cm);

        std::pair<double, double> lora1 = Lorentz(gamma_cm,beta_cm,CM.e_scat,P1cm.Z());	
        std::pair<double, double> lora2 = Lorentz(gamma_cm,beta_cm,CM.e_clust,P2cm.Z());

        P1cm.SetZ(lora1.second);
        P2cm.SetZ(lora2.second);
        //-------- Rotating back to the beam direction -----------
        TVector3 direction = LVstart.Vect();
        TVector3 P1L;
        TVector3 P2L;


        //if(is_direct)
        //{
        direction = direction.Unit();
        P1cm.RotateUz(direction);
        P2cm.RotateUz(direction);
        P1L = P1cm;
        P2L = P2cm;
        //}

        //else
        //{
        //    P1L = DREHUNG(P1cm,Pa);
        //    P2L = DREHUNG(P2cm,Pa);
        //}

        //Double_t coin = r1.Unifrom(-1.,1.);

        //--------- Filling in the ROOTTree variables------------
        theta_1 = P1L.Theta();		
        theta_2 = P2L.Theta();		
        theta_B = PB.Theta();		

        phi_1 = P1L.Phi();		
        phi_2 = P2L.Phi();		
        phi_B = PB.Phi();		

        P1x = P1L.X();
        P1y = P1L.Y();
        P1z = P1L.Z();

        P2x = P2L.X();
        P2y = P2L.Y();
        P2z = P2L.Z();

        PBz_lab = PB.Z();

        E1 = sqrt(Mi*Mi + P1L.Mag2()) - Mi;
        E2 = sqrt(Ma*Ma + P2L.Mag2()) - Ma;
        EB = sqrt(MB*MB + PB.Mag2()) - MB;

        th1_cm  = CM.theta_scat;
        th2_cm  = CM.theta_clust; 
        P1_cm 	= CM.p_scat;
        P2_cm	= CM.p_clust;
        Moff 	= Ma_off;
        Mandelstam_T = CM.T;
        Opang=acos(sin(P1L.Theta())*sin(P2L.Theta())*cos(P1L.Phi()-P2L.Phi())+cos(P1L.Theta())*cos(P2L.Theta()));

        double df = fabs(phi_1-phi_2);

        if(df > 0 && df <= PI)   Dif_phi = df;
        else  Dif_phi = 2*PI - df;

        if(events%10000==0) cout<< "\r" << events <<" of "<<nevents<<" ("<<(float)events/nevents*100<<"%)"<<flush;	

        outfile << PBx/1000 << "\t" << PBy/1000 << "\t" << PBz_lab/1000 << "\t" <<P1x/1000 << "\t"<< P1y/1000 << "\t" << P1z/1000 << "\t"<< P2x/1000 << "\t" << P2y/1000 <<  "\t" << P2z/1000 << "\n";

        tree->Fill();
        events++;

    }

    outfile.close();

    //tree->Print();
    tree->AutoSave();
    if(theory_mom_filename) h_momentum_theory->Write();
    file.Close();

    return;
}

//---- Two consecutive rotations 
//first around Z on <phi>, then around new X' on <theta> (1=Pcm, 2=Pa in lab)
TVector3 DREHUNG(TVector3 v1,TVector3 v2) 
{
    double CT = v2.Z()/v2.Mag(); // cos(theta) of v2 wrt. Z-axis
    double ST = sqrt(1-CT*CT);   // sin(theta)
    double CF = v2.X()/v2.Mag()/ST;
    double SF = v2.Y()/v2.Mag()/ST;

    TVector3 v3;
    double _v3x =  v1.X()*CT*CF - v1.Y()*SF + v1.Z()*ST*CF;
    double _v3y =  v1.X()*CT*SF + v1.Y()*CF + v1.Z()*ST*SF;
    double _v3z = -v1.X()*ST   +  v1.Z()*CT;
    v3.SetXYZ(_v3x,_v3y,_v3z);
    return v3;
}

//Kinematical function
double CINEMA(double x,double y,double z)
{	
    double lambda = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
    return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
cm_values CENMASS(double s,double m2off,double m1,double m2)
{
    cm_values output;
    output.good=false;
    double X = s;
    double Y = m2off*m2off;
    double Z = m1*m1;
    double sqrs = sqrt(s);

    // Kinematics before the scattering process
    // (with one off-shell mass)
    double p2_off = sqrt(CINEMA(X,Y,Z))/2/sqrs;
    double p1_off = p2_off;
    // CM energies
    double e1_off = (s+Z-Y)/2/sqrs;
    double e2_off = (s+Y-Z)/2/sqrs;

    // Now take the real masses (after scattering)
    Y = m2*m2;  Z = m1*m1;
    //And check whether the kinematical function is ok
    //for this specific kinematical case
    double ERROR_CI = CINEMA(X,Y,Z);
    if(ERROR_CI <= 0.){
        return output;
    }

    // Kinematics after the scattering process
    // (with all real masses)
    double p2 = sqrt(CINEMA(X,Y,Z))/2/sqrs;
    double p1 = p2;
    double e1 = (s+Z-Y)/2/sqrs;
    double e2 = (s+Y-Z)/2/sqrs;

    // Let's consider momentum transfer <t> from the
    // target particle 1 to the cluster 2
    double tmax = 2*(m1*m1 - e1_off*e1 - p1_off*p1);//COSINE=(-1)
    double tmin = 2*(m1*m1 - e1_off*e1 + p1_off*p1);//COSINE=(1)

    double cmth;
    double t;

    // Generate random momentum transfer for this kinematical case
    //	if(!isotropic){ t = get_T(s,tmax); }			//Using parameterized cross sections
    if(!isotropic){
        cmth = f1->GetRandom(0.,180.);
        t= m2off*m2off + m2*m2 - 2*e2_off*e2 + 2*p2_off*p2*cos(cmth*PI/180.);
        //	double t_max = -2*p1_off*p1_off*(1-cos(PI));
        //	Double_t rr = gRandom->Uniform(-1.,1.);// to randomize wrt 90 degrees
        //	if(rr>0) t = t_max - t;
    }			//Using parameterized cross sections

    else{ t = gRandom->Uniform(0.,3000000.) * (-1); } //Isotropic scattering


    double COSINE = (t - 2*m1*m1 + 2*e1_off*e1)/(2*p1_off*p1);
    if(fabs(COSINE)>=1){//momentum transfer out of range
        return output;
    }

    //CM scattering angles
    double theta1 = acos(COSINE);
    double theta2 = PI - theta1;

    output.e_clust = e2;
    output.p_clust = p2;
    output.theta_clust = theta2;

    output.e_scat = e1;
    output.p_scat = p1;
    output.theta_scat = theta1;

    output.T = t;
    output.good = true;

    return output;
}

// Calculate 3-momentum in CM system of two particles M1 and M2
// when M1 has kinetic energy TLAB and M2 is at rest
double momentum_CM(double TLAB,double M1,double M2)
{
    //Particle M2 is assumed to be in rest
    double PLAB = sqrt(TLAB*(TLAB + 2*M1));//  Total 3-momentum of an incident particle in Lab
    double ELAB = sqrt(PLAB*PLAB + M1*M1); //  Total energy of an incident particle in lab
    double SLAB = M1*M1 + M2*M2 + 2*M2*ELAB;// Mandelstam invariant S in lab
    double PCM = PLAB * M2 / sqrt(SLAB); // Momentum of both particles in CM frame
    return PCM;
}

std::pair<double, double> Lorentz(double g,double b,double e,double p)
{
    double eL = g*e - g*b*p;
    double pL = g*p - g*b*e;
    return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c)² units)
// distributed according to the parameterized proton-proton
// invarant cross section. Pass as a parameter "sm" the 
// Mandelstam variable S (in MeV²)
// and the maximum possible momentum transfer  
double get_T(double sm, double max)
{
    double Tmax = max*0.000001;//convert to GeV² units

    Double_t rr = gRandom->Uniform(-1.,1.);// to randomize wrt 90 degrees
    double mandels = sm * 0.000001; // in GeV²

    //Probability function from the parameterization
    TF1 * foo = new TF1("foo","[0]*exp(x*[1])*(1+0.02*exp((-6)*x))",Tmax/2,0);
    double c = 0.; 
    if(mandels<=4.79) c = -3283.75 + 3064.11*mandels - 1068.44 *mandels*mandels + 164.844*pow(mandels,3) - 9.48152*pow(mandels,4);
    else if(mandels>4.79) c = -776.822 + 586.016*mandels - 175.347 *mandels*mandels + 26.1823*pow(mandels,3) - 1.94889*pow(mandels,4) + 0.0578352*pow(mandels,5);

    foo->FixParameter(0, 25.); //normalization constant (could be anything)
    foo->FixParameter(1, c);	

    double Trand = foo->GetRandom(Tmax/2,0.); //from 90 to 0 degrees
    if(rr>0) Trand = Tmax - Trand; // symmetrization relative to 90 degrees

    delete foo;
    return (Trand*1000000); // returning value in MeV²
}

//This is a function to read theoretical momentum distribution (e.g. C.Bertulani)
//Input ASCII file has only two columns: 1: Total Momentum 2: Cross section(Probability)
//The returned value is a histogram which will be randomly sampled inside the main function run()
TH1F* Theoretical_Momentum(char* filename)
{
    //===== Reading data from the input file =====
    int NPoints=0;//number of line in the input file
    int i=0;//iterator
    ifstream infile(filename);
    std::string line;
    while ( std::getline(infile, line) ) NPoints++;
    infile.clear();
    infile.seekg( 0, std::ios::beg );
    float Q[NPoints], dS_dQ[NPoints];//data from input file
    while(1)
    {
        infile >> Q[i] >> dS_dQ[i];
        if(!infile.good()) break;
        i++;
    }
    infile.close();

    //===== Now making histogram from Graph and compute Histogram bins
    TGraph * gQ = new TGraph(NPoints,Q,dS_dQ);
    Double_t BinLimits[NPoints+1];
    gQ->Sort();
    // determine lower limit of histogram: half the distance to next point
    Double_t x0,x1,y;
    gQ->GetPoint(0,x0,y);
    gQ->GetPoint(1,x1,y);
    Double_t Distance = TMath::Abs(x0-x1);
    BinLimits[0] = x0 - Distance/2.;
    // now set upper limits for all the other points
    for (Int_t k = 0 ; k<NPoints-1;k++)
    {
        gQ->GetPoint(k,x0,y);
        gQ->GetPoint(k+1,x1,y);
        Distance = TMath::Abs(x0-x1);
        BinLimits[k+1] = x0 + Distance/2.;
    }
    // for the last point set upper limit similar to first point:
    gQ->GetPoint(NPoints-2,x0,y);
    gQ->GetPoint(NPoints-1,x1,y);
    Distance = TMath::Abs(x0-x1);
    BinLimits[NPoints] = x1 + Distance/2.;
    // now we know the binning and can create the histogram:
    TString Name = "Theoretical_Momentum"; 
    // make name unique 
    //Name+= rand();
    TH1F *ThisHist = new TH1F(Name,Name,NPoints,BinLimits);
    // now fill the histogram
    for (Int_t i = 0; i<gQ->GetN();i++)
    {
        Double_t x,y;
        gQ->GetPoint(i,x,y);
        ThisHist->SetBinContent(i+1,y);
        ThisHist->SetBinError(i+1,1e-14);
    }
    return ThisHist;
}


int main(Int_t argc, Char_t* argv[])
{
    Char_t *root_filename=0;
    Char_t *ascii_filename=0;
    Char_t *input_filename=0;

    char* theory_mom_filename=NULL;

    int nevents = MAX_STORY;

    Bool_t NeedHelp = kTRUE;
    Bool_t is_ppn = kFALSE;
    Bool_t is_nnp = kFALSE;
    Bool_t is_ppa = kFALSE;
    Bool_t is_isotropic = kTRUE;
    Bool_t is_direct = kFALSE;

    if (argc > 1)
    {
        for (Int_t i = 0; i < argc; i++)
        {
            if (strncmp(argv[i],"--root=",7) == 0)
            {
                root_filename = argv[i]+7;
                NeedHelp = kFALSE;
            }

            else if (strncmp(argv[i],"--ascii=",8) == 0)
            {
                ascii_filename = argv[i]+8;
            }

            else if (strncmp(argv[i],"--input=",8) == 0)
            {
                input_filename = argv[i]+8;
            }

            else if (strncmp(argv[i],"--theory_mom=",13) == 0)
            {
                theory_mom_filename = argv[i]+13;
            }

            else if (strncmp(argv[i],"--max-events=",13) == 0)
            {
                nevents = atoi(argv[i]+13);
            }

            else if (strncmp(argv[i],"--ppn",5) == 0)
            {
                is_ppn = kTRUE;
                std::cout << "\n-- Generate  (p,pn)-data";
            }

            else if (strncmp(argv[i],"--nnp",5) == 0)
            {
                is_nnp = kTRUE;
                std::cout << "\n-- Generate  (n,np)-data";
            }


            else if (strncmp(argv[i],"--ppa",5) == 0)
            {
                is_ppa = kTRUE;
                std::cout << "\n-- Generate (p,pa)-data";
            }

            else if (strncmp(argv[i],"--direct-kinematics",19) == 0)
            {
                is_direct = kTRUE;
                std::cout << "\n-- Assuming direct kinematics";
            }

            else if (strncmp(argv[i],"--not-iso",9) == 0)
            {
                is_isotropic = kFALSE;
                std::cout << "\n-- Generate anisotropic nucleon-distributions";
            }
        }
    }

    if (NeedHelp)
    {
        std::cout << "\nOptions:\n\n";
        std::cout << "  --direct-kinematics -> Assuming direct reaction kinematics\n\n";
        std::cout << "  --ppn -> Generate (p,pn)-data\n\n";
        std::cout << "  --nnp -> Generate (n,np)-data -> only with is_direct option\n\n";
        std::cout << "  --ppa -> Generate (p,pa)-data.  If no option is specified events will be (p,2p).\n\n";
        std::cout << "  --max-events -> Number of generated events\n\n";
        std::cout << "  --input=/path/to/your/file/filename.root -> root file with momentum histograms, type 'gauss' to use Gaussian distribution with width specified in info.hh instead.\n\n";
        std::cout << "  --ascii=/path/to/your/file/filename.dat -> ascii file name\n\n";
        std::cout << "  --theory_mom=/path/to/your/file/filename.dat -> input ASCII file with Theoretical_Momentum distribution inside A nucleus. Expecting two-column format (1st: Total momentum, 2nd: Cross section(Probability) )\n\n"; 
        std::cout << "  --root=/path/to/your/file/filename.root -> output root file name\n\n";
        std::cout << "  --not_iso -> Generate anisotropic nucleon-distributions.  If no option is specified events will be homegeneous.\n\n";
        return 0;
    }


    //Main generator function
    if(root_filename!=NULL && ascii_filename!=NULL) 
    {
        run(root_filename,
                ascii_filename,
                input_filename,
                theory_mom_filename,
                nevents,
                is_ppn,
                is_ppa,
                is_isotropic,
                is_direct,
                is_nnp);
    }

    else if(ascii_filename==NULL)
    {
        std::cout << "\nASCII FILE IS NOT SPECIFIED!!!\nType: ./qfs \n\n" << std::endl;
        return 0;
    }
    //else if (input_filename==NULL)
    //{
    //    std::cout << "\nMOMENTUM INPUT FILE IS NOT SPECIFIED!!!\nType: ./qfs \n\n" << std::endl;
    //    return 0;
    //}
    else
    {
        std::cout << "\nROOT FILE IS NOT SPECIFIED!!!\nType: ./qfs \n\n" << std::endl;
        return 0;
    }

    return 0;
}
