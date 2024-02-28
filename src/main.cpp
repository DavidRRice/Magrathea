#include "define.h"
#include "EOS.h"
#include "EOSmodify.h"
#include "phase.h"
#include "hydro.h"
#include "compfind.h"

struct timeval start_time, end_time;

const double rho_eps_rel = 1E-11;	// relative error tolerance of the density solver
const double T_eps_rel = 1E-11;	// relative error tolerance of the temperature solver for the adiabatic profile
const double ode_eps_rel0 = 1E-7; // relative error tolerance for first round ode integrator (1E-7) used in mode 0
const double ode_eps_rel1 = 1E-10; // relative error tolerance for second round ode integrator (1E-10) used in mode 0
int fit_iter;
const double R_eps_rel = 2E-5; // relative error tolerance in mode 0 first round radius determination (5E-5).  Should be around sqrt(ode_eps_rel0).
const double ode_eps_rel2 = 1E-10; // relative error tolerance for ode integrator (1E-10) used in mode 1
const double P_eps_rel = 1E-10;	// relative error tolerance in mode 1 central pressure determination (1E-10).  Should not be more restrict than ode_eps_rel.
const double fit_eps_rel = 1E-4; // the relative error tolerance at the fitting point in mode 0 round 2 (1E-4). Should be four orders of magnitudes larger than ode_eps_rel1.
vector<double> ave_rho = {15, 5, 2, 1E-3};// Assuming the density of the core is 15, mantle is 5, water is 2, and gas is 1E-3.
const bool verbose = false;		  // Whether print warnings.
const double P_surface = 1E5;		  // The pressure level that the broad band optical transit radius probes. (in microbar)
int count_shoot = 0;			  // used to count the total number of shootings per each solution
int count_step = 0;			  // used to count the sum of integral steps in all shooting results.

int main()
{
  gsl_set_error_handler_off();  //Dismiss the abortion from execution, which is designed for code testing.
  hydro* planet;
  ifstream fin;
  int input_mode=0;
  // Choose between the 8 available input_mode values:
  // 0: regular solver, 1: temperature-free solver, 2: two-layer solver, 
  // 3: modify a built-in EOS on they fly, 
  // 4: iterate over EOS modifications with two-layer solver, 5: iterate over EOS with regular solver
  // 6: bulk input mode with regular solver
  // 7: composition finder, secant method to find third layer mass to match a mass and radius measurement

  if (input_mode == 0)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0, 0, 0, 300};
    // The temperature of the outer boundary of the inner component minus the inner boundary of the outer component.  A positive number indicates temperature increases inward.  0 indicates the temperature is continuous at the boundary of components.  The last number is the planetary surface temperature.
    vector<double> Mcomp =  {1.0,0.5,0.1,0.00001}; // Mass in Earth Masses of Core, Mantle, Hydrosphere, Atmosphere
    planet=fitting_method(Comp, Mcomp, Tgap, ave_rho, P_surface, false);
    cout<<count_shoot<<' '<<count_step<<endl;
    if (!planet)
    {
      for (unsigned int i=0; i < Mcomp.size(); i++)
	cout<<Mcomp[i]<<", ";
      cout<<"\t No solution found."<<endl;
    }
    else
      planet->print("./result/Structure.txt", true); // Save the result in an asc file with this name.

    delete planet;
  }  

  else if(input_mode == 1)
  {
    planet=getmass(3,3,3,P_surface);
    // Mass in Earth Masses of Core, Mantle, Hydrosphere
    if (planet)
      planet->print("./result/Structure.txt");
    // Save the result in an asc file with this name.
    delete planet;
  }
  
  else if(input_mode == 2)
  {
    vector<double> Mp,Rp;
    double deltat;
    gettimeofday(&start_time,NULL);

    twolayer(0,0,Mp,Rp,P_surface,true);

    for(int i=0; i < int(Mp.size()); i++)
      cout<<Mp[i]<<" \t"<<Rp[i]<<endl;
    gettimeofday(&end_time, NULL);
  
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

    cout<<"running time "<<deltat<<'s'<<endl;
  }

  else if(input_mode == 3)
  {
    vector<double> Mp,Rp;
    double fraction = 0;
    
    int nPhases=3;
    EOS **highEOSs = new EOS*[nPhases];
    for(int i=0; i<nPhases; i++)
      highEOSs[i] = new EOS();
    highEOSs[0] -> setphasename("Ice VII");
    highEOSs[1] -> setphasename("Ice VII'");
    highEOSs[2] -> setphasename("Ice X");

    double IceVII_Grande[5][2]={{0.,0.},{1.,13.62},{2.,10.48},{3.,4.},{5.,18.01528}};
    double IceVIIp_Grande[5][2]={{0.,0.},{1.,12.32},{2.,22.36},{3.,4.},{5.,18.01528}};
    double IceX_Grande[5][2]={{0.,0.},{1.,11.0},{2.,38.6},{3.,4.},{5.,18.01528}};
    double start_P[2]={4.76,34.14};
    highEOSs[0] -> modifyEOS(IceVII_Grande,5);
    highEOSs[1] -> modifyEOS(IceVIIp_Grande,5);
    highEOSs[2] -> modifyEOS(IceX_Grande,5);

    water.set_phase_highP(nPhases, start_P, highEOSs);
    planet = getmass(0,0,0.0999968,P_surface);
  
    if (planet)
      planet->print("./result/Structure.txt");
  
    twolayer(0,fraction,Mp,Rp,P_surface);

    for(int i=0; i < int(Mp.size()); i++)
      cout<<Mp[i]<<" \t"<<Rp[i]<<endl;

    for(int i=0; i<nPhases; i++)
      delete highEOSs[i];
    delete[] highEOSs;
  }  

  else if(input_mode == 4)
  {
    double deltat;
    double fraction = 1;
    int index = 0;
    
    gettimeofday(&start_time,NULL);

    twolayer(index,fraction,P_surface,1,"./run/PosteriorsPPv.txt", "./result/PPvMCMC.txt");
  
    gettimeofday(&end_time, NULL);
  
    deltat = ((end_time.tv_sec  - start_time.tv_sec) * 1000000u + end_time.tv_usec - start_time.tv_usec) / 1.e6;

    cout<<"running time "<<deltat<<'s'<<endl;
  }
  
  else if(input_mode == 5)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0, 0, 0, 300};
    vector<double> Mcomp =  {0.29,0.69,1.02,0.001};

    fullmodel(Comp,Mcomp,Tgap,ave_rho,P_surface,false,2,"./run/PosteriorsVinet.txt", "./result/IceMCMCnew.txt");
  }

  else if(input_mode == 6)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0,0,0,300}; // Using full temperature solver, Temperature gap between each layer and surface temperature.

    string infilename="./run/inputcore.txt";
    string outfilename="./result/coreplanets.txt";

    int solver=1; //Calculation mode. Mode 2 is twice faster, but only applicable for planet with no gas layer and whose temperature effect is not important

    
    multiplanet(Comp, Tgap, solver, ave_rho, P_surface, false, infilename, outfilename);
    

  }

  else if(input_mode == 7)
  {
    vector<PhaseDgm> Comp = {Fe, Si, water, atm};
    vector<double> Tgap = {0,0,0,300}; // Using full temperature solver, Temperature gap between each layer and surface temperature.

    string infilename="./run/T1hpost.txt";
    string outfilename="./result/T1houtput.txt";

    //Solver will hold 2 layers in constant Partial Mass Ratio (PMF) and find the mass of 3rd layer: 
    //PMR(%) is OMF/(IMF+OMF)*100, where OMG is outer-mass fraction, IMF is inner-mass fraction
    //Solver can be looped through steps in PMR
    int findlayer=3; //1 to find core, 2 for mantle, 3 for water, 4 for atmosphere mass fraction
    vector<int> layers={1,1,0,0}; //mark which layers {C,M,W,A} to hold in constant total ratio
    double minPMR=67.0; //Minimum PMR(%) for iteration, must be multiple of 0.1
    double maxPMR=67.0; //Maximum PMR(%) for iteration, must be multiple of 0.1
    float step=1.0; // Step size for partial mass ratio, must be multiple of 0.1

    double rerr=0.001; // Error in the simulated radius to target radius

    //Multi-threading is commented out by default for easy install
    //Must uncomment Line 2 in Makefile and all occurences of "pragma..." in comfind.cpp 
    int num_threads=0; 
    
    compfinder(Comp,findlayer,layers,minPMR,maxPMR,step,rerr,num_threads,Tgap,ave_rho,P_surface,false,infilename, outfilename);

  }

  // ============================
  
  delete Fe_liquid;
  delete Fe_liquid2;
  delete Fe_hcp;
  delete Fe_hcp2;
  delete Fe_hcp3;
  delete Fe_Dummy;
  delete Fe_7Si;
  delete Fe_15Si;
  delete Fe_Seager;
  delete Si_Pv_Shim;
  delete Si_Pv;
  delete Si_PPv;
  delete Si_PPv_Sakai;
  delete Si_PREM;
  delete Si_BM2fit;
  delete Si_Seager;
  delete Si_liquid;
  delete Si_Liquid_Wolf;
  delete Si_Dummy;
  delete Ice_Seager;
  delete Water_ExoPlex;
  delete Water;
  delete IceIh_ExoPlex;
  delete IceVI_ExoPlex;
  delete IceVII_ExoPlex;
  delete IceVII;
  delete IceVIIp;
  delete IceVII_FFH2004;
  delete IceVII_FFH2004fit;
  delete IceVII_FFH2004BM;
  delete IceX_HS;
  delete IceX;
  delete IceZeng2013FFH;
  delete IceZeng2013FMNR;
  delete Ice_Dummy;
  delete Gas;

  return 0;
}
