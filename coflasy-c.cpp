// COFLASY-C: COsmological evolution of FLavour ASYmmetries - C++
// by Valerie Domcke (valerie.domcke@cern.ch), Miguel Escudero (miguel.escudero@cern.ch), Mario Fernandez Navarro (mario.fernandeznavarro@physik.uzh.ch), Stefan Sandner (stefan.sandner@lanl.gov)

// If you use this code for scientific publications, please cite these papers:

// "Lepton Flavour Asymmetries: from the early Universe to BBN."  Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner,
// [JHEP 06, 137] ( https://link.springer.com/article/10.1007/JHEP06(2025)137 ), [arXiv:2502.14960] ( https://arxiv.org/abs/2502.14960 ), [INSPIRE] ( https://inspirehep.net/literature/2893306 ).
// "A Limit on the Total Lepton Number in the Universe from BBN and the CMB." Valerie Domcke, Miguel Escudero, Mario Fernandez Navarro and Stefan Sandner, [arxiv:2510.02438] ( https://arxiv.org/abs/2510.02438 ).

//Last major update: 03-10-2025
//v1.0


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <chrono>
#include <thread>
#include <ctime>  
#include <array>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multiroots.h>

#include "include/common.h"

#include "source/read_ini.cpp"
#include "source/Collision_interp1D.cpp"
#include "source/Collision_interp2D.cpp"
#include "source/Collision_interp2D_me.cpp"
#include "source/Transfer_interp1D.cpp"
#include "source/Transfer_interp2D.cpp"
#include "source/Density_interp.cpp"
#include "source/MatOp.cpp"
#include "source/ParamTranslate.cpp"
#include "source/Xn_interp.cpp"
#include "source/WeakRates_interp1D.cpp"
#include "source/WeakRates_interp2D.cpp"
#include "source/WeakRatesCorrection_interp1D.cpp"


using namespace std;
typedef complex<double> complx;


// Global flag to indicate an error has occurred
bool gsl_error_flag = false;


// Custom GSL error handler
void custom_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
    // std::cerr << "GSL Error: " << reason << " in " << file << " at line " << line << std::endl;
    gsl_error_flag = true;  // Set the error flag
    // Do not call abort(), allowing the program to continue
}



// function for the 8 GellMann matrices l1,l2,..,l8 and the unit matrix l0
void GellMannMatrices(complx (&lvec) [9][3][3]){
      matop.mat33_free(lvec[0]); matop.mat33_free(lvec[1]); matop.mat33_free(lvec[2]); 
      matop.mat33_free(lvec[3]); matop.mat33_free(lvec[4]); matop.mat33_free(lvec[5]); 
      matop.mat33_free(lvec[6]); matop.mat33_free(lvec[7]); matop.mat33_free(lvec[8]);
      
      lvec[0][0][0]=1.;            lvec[0][1][1]=1.;            lvec[0][2][2]=1.;
      lvec[1][0][1]=1.;            lvec[1][1][0]=1.;
      lvec[2][0][1]=-1i;           lvec[2][1][0]=1i;
      lvec[3][0][0]=1;             lvec[3][1][1]=-1;
      lvec[4][0][2]=1;             lvec[4][2][0]=1;
      lvec[5][0][2]=-1i;           lvec[5][2][0]=1i;
      lvec[6][1][2]=1;             lvec[6][2][1]=1;
      lvec[7][1][2]=-1i;           lvec[7][2][1]=1i;
      lvec[8][0][0]=1./sqrt(3);    lvec[8][1][1]= 1./sqrt(3);   lvec[8][2][2]=-2./sqrt(3);
}


// function returning the real part of the Trace of a 3x3 matrix
double trace(complx M[3][3]){
      complx trace = M[0][0] + M[1][1] + M[2][2];
      return trace.real();
}


// initial condition in the flavor basis for r. rbar is then just obtained by xi->-xi
void ini_flavor(double xi1, double xi2, double xi3, complx (&r_flavor_ini) [3][3]){
      double prefac = 2./(3.*_Zeta3_) * gsl_sf_gamma(3.);
      matop.mat33_free(r_flavor_ini);
      r_flavor_ini[0][0]      = (prefac * gsl_sf_fermi_dirac_int(2, xi1));
      r_flavor_ini[1][1]      = (prefac * gsl_sf_fermi_dirac_int(2, xi2));
      r_flavor_ini[2][2]      = (prefac * gsl_sf_fermi_dirac_int(2, xi3));

}


/*
calculate initial epsilon detailed balance parameter, znu, zg such that diagonal elements of collision term
as well as transfer rate vanish, see 2502.14960 Sec. 2B
*/

// Struct for passing multiple parameters
struct Params {
    double p[3];   // Array to store parameters
    double x_min, x_max;
    double y_min, y_max;
    double penalty_factor; // Large penalty factor
};

// Penalty function to add constraints
double penalty(double val, double min_val, double max_val, double k) {
    if (val < min_val) return k * (min_val - val) * (min_val - val);
    if (val > max_val) return k * (val - max_val) * (val - max_val);
    return 0.0;
}


// Define the system of nonlinear equations
int ini_equations(const gsl_vector *variables, void *params, gsl_vector *f){
      Params *p = (Params *)params;
      double xi1  = p->p[0]; 
      double xi2  = p->p[1]; 
      double xi3  = p->p[2]; 

      double penalty_factor = p->penalty_factor;

      double tini= 1;

      double prefac = 2./(3.*_Zeta3_) * gsl_sf_gamma(3.);
      
      
      double rfd1, rfd2, rfd3;
      double rbfd1, rbfd2, rbfd3;
      rfd1  = prefac * gsl_sf_fermi_dirac_int(2, xi1); 
      rfd2  = prefac * gsl_sf_fermi_dirac_int(2, xi2); 
      rfd3  = prefac * gsl_sf_fermi_dirac_int(2, xi3);
      rbfd1 = prefac * gsl_sf_fermi_dirac_int(2, -xi1); 
      rbfd2 = prefac * gsl_sf_fermi_dirac_int(2, -xi2); 
      rbfd3 = prefac * gsl_sf_fermi_dirac_int(2, -xi3);


      double znu  = gsl_vector_get(variables, 0);
      double zg   = gsl_vector_get(variables, 1);



      /* ------- set the rates ------- */
      // collision term 0. moment
      // 1D in znu
      double A12nunuee;        
      // 2D in zg, x
      double A34nunuee;     
      // 2D in znu, zg
      double A123nunuee, A134nunuee;    
      /* ------ */
      // transfer rates 1. moment
      double A12nuenueI1, A34nuenueI1;
      double A12nunueeI1, A34nunueeI1;

      // collision term 0. moment

      // 1D rates in znu
      A12nunuee   = collision1D.A12_nunu_ee_I0(znu);

      // 2D rates in zg, x
      A34nunuee   = collision2D_me.A34_nunu_ee_I0_me(zg, _Tref_/tini);

      // 2D rates in znu, zg
      A123nunuee  = collision2D.A123_nunu_ee_I0(znu, zg);
      A134nunuee  = collision2D.A134_nunu_ee_I0(znu, zg);

      /* --------- */

      // transfer rates 1. moment
      
      A12nuenueI1       = transfer2D.A12_nue_nue_I1(znu, zg);
      A34nuenueI1       = transfer2D.A34_nue_nue_I1(znu,zg);

      A12nunueeI1       = transfer1D.Aij_nunu_ee_I1(znu);
      A34nunueeI1       = transfer1D.Aij_nunu_ee_I1(zg);

      /* ----------------------------- */


      double epsdb1, epsdb2, epsdb3;

      epsdb1      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd1 + rfd1) + A12nunuee*(rbfd1 + rfd1) + 
                        2*A134nunuee*pow(znu,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd1 + rfd1,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znu,6),0.5));

      epsdb2      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd2 + rfd2) + A12nunuee*(rbfd2 + rfd2) + 
                        2*A134nunuee*pow(znu,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd2 + rfd2,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znu,6),0.5));

      epsdb3      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd3 + rfd3) + A12nunuee*(rbfd3 + rfd3) + 
                        2*A134nunuee*pow(znu,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd3 + rfd3,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znu,6),0.5));



      double f1, f2;
      /* deps/dt = 0, eq. (2.16) in 2502.14960 */
      f1 = -0.03125*(pow(_PI_,-5)*pow(znu,-6)*(0. + 4.5804263131994185*A12nunueeI1*(epsdb1 + rbfd1)*(epsdb1 + rfd1) + 
            1.0116205281900887*A12nunueeI1*(epsdb2 + rbfd2)*(epsdb2 + rfd2) + 
            1.0116205281900887*A12nunueeI1*(epsdb3 + rbfd3)*(epsdb3 + rfd3) + 
            2.2902131565997093*A12nuenueI1*(epsdb1 + rbfd1)*pow(znu,3) - 2.2902131565997093*A34nuenueI1*(epsdb1 + rbfd1)*pow(znu,3) - 
            0.5058102640950444*(-A12nuenueI1 + A34nuenueI1)*(epsdb2 + rbfd2)*pow(znu,3) - 
            0.5058102640950444*(-A12nuenueI1 + A34nuenueI1)*(epsdb3 + rbfd3)*pow(znu,3) + 
            2.2902131565997093*A12nuenueI1*(epsdb1 + rfd1)*pow(znu,3) - 2.2902131565997093*A34nuenueI1*(epsdb1 + rfd1)*pow(znu,3) - 
            0.5058102640950444*(-A12nuenueI1 + A34nuenueI1)*(epsdb2 + rfd2)*pow(znu,3) - 
            0.5058102640950444*(-A12nuenueI1 + A34nuenueI1)*(epsdb3 + rfd3)*pow(znu,3) - 6.603667369579596*A34nunueeI1*pow(znu,6)));
      /* Hubble(1,1) = Hubble(znu,zg), see 2502.14960 below eq. (2.31) */
      f2 = (44 + 7*(2*epsdb1 + 2*epsdb2 + 2*epsdb3 + rbfd1 + rbfd2 + rbfd3 + rfd1 + rfd2 + 
            rfd3) - 44*pow(zg,4) - 
            7*(2*epsdb1 + 2*epsdb2 + 2*epsdb3 + rbfd1 + rbfd2 + rbfd3 + rfd1 + rfd2 + 
            rfd3)*pow(znu,4))/8.;




      // Apply penalty function to "discourage" values outside the range
      f1 += penalty(znu, p->x_min, p->x_max, penalty_factor);
      f1 += penalty(zg, p->y_min, p->y_max, penalty_factor);
      f2 += penalty(znu, p->x_min, p->x_max, penalty_factor);
      f2 += penalty(zg, p->y_min, p->y_max, penalty_factor);

      gsl_vector_set(f, 0, f1);
      gsl_vector_set(f, 1, f2);

      return GSL_SUCCESS;
}


// solve the system
void ini_correction(double xi_zon_ini[4], 
                        complx (&r_flavor_ini) [3][3], complx (&rb_flavor_ini) [3][3], 
                        double (&znu_ini), double (&zg_ini) ){
      

      int zon = static_cast<int>(xi_zon_ini[3]);

      double znuSol, zgSol;
      if (zon == 1){
            // Define the dimension (two variables x, y)
            const size_t dim = 2;

            // Set initial guess for znu and zg
            gsl_vector *variables = gsl_vector_alloc(dim);
            gsl_vector_set(variables, 0, 1.01);       // Initial guess for znu
            gsl_vector_set(variables, 1, 0.99);       // Initial guess for zg

            // Define constraints and parameters
            Params params = {{xi_zon_ini[0], xi_zon_ini[1], xi_zon_ini[2]}, 0.8, 1.2, 0.8, 1.2, 1e3};



            gsl_multiroot_function f = {&ini_equations, dim, &params};

            // Choose a solver type
            const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
            gsl_multiroot_fsolver *solver = gsl_multiroot_fsolver_alloc(T, dim);

            // Initialize the solver
            gsl_multiroot_fsolver_set(solver, &f, variables);

            // Iterate until convergence
            int status;
            size_t iter = 0;
            const size_t max_iter = 1000;
            const double tol = 1e-12;
            
            cout << "Find self consistent initial conditions... ";
            do{
                  iter++;
                  status = gsl_multiroot_fsolver_iterate(solver);
                  if (status){
                        cout << "\n\033[31mGSL multiroot failure at iteration i = " << iter << "\033[39;49m" << endl;
                        break;
                  }
                  status = gsl_multiroot_test_residual(solver->f, tol);
            } while (status == GSL_CONTINUE && iter < max_iter);


            // Output the solution
            znuSol      = gsl_vector_get(solver->x, 0);
            zgSol       = gsl_vector_get(solver->x, 1);

            // Ensure the solution is within the desired range
            if (status == GSL_SUCCESS && znuSol > params.x_min && znuSol < params.x_max && zgSol > params.y_min && zgSol < params.y_max){
                  cout << "\033[32m\n...Solution found after " << iter << " iterations\033[0m" << endl;
                  cout << std::fixed << std::setprecision(4) << "znu = " << znuSol << ", " << "zg = " << zgSol;
            }
            // Free multiroot resources
            gsl_multiroot_fsolver_free(solver);
            gsl_vector_free(variables);
      }
      else{
            znuSol= 1; zgSol= 1;
            cout << std::fixed << std::setprecision(4) << "znu = " << znuSol << ", " << "zg = " << zgSol;
      }


      matop.mat33_free(r_flavor_ini);

      double prefac = 2./(3.*_Zeta3_) * gsl_sf_gamma(3.);
      double xi1,xi2,xi3;
      xi1=xi_zon_ini[0]; xi2=xi_zon_ini[1]; xi3=xi_zon_ini[2];

      double rfd1, rfd2, rfd3;
      double rbfd1, rbfd2, rbfd3;
      rfd1  = prefac * gsl_sf_fermi_dirac_int(2, xi1); 
      rfd2  = prefac * gsl_sf_fermi_dirac_int(2, xi2); 
      rfd3  = prefac * gsl_sf_fermi_dirac_int(2, xi3);
      rbfd1 = prefac * gsl_sf_fermi_dirac_int(2, -xi1); 
      rbfd2 = prefac * gsl_sf_fermi_dirac_int(2, -xi2); 
      rbfd3 = prefac * gsl_sf_fermi_dirac_int(2, -xi3);

      double epsdb1, epsdb2, epsdb3;
      double tini= 1;

      /* ------- set the rates ------- */
      // collision term 0. moment
      // 1D in znu
      double A12nunuee;       
      // 2D in zg, x
      double A34nunuee;     
      // 2D in znu, zg
      double A123nunuee, A134nunuee;    

      // collision term 0. moment
      // 1D rates in znu
      A12nunuee   = collision1D.A12_nunu_ee_I0(znuSol);
      // 2D rates in zg, x
      A34nunuee   = collision2D_me.A34_nunu_ee_I0_me(zgSol, _Tref_/tini);
      // 2D rates in znu, zg
      A123nunuee  = collision2D.A123_nunu_ee_I0(znuSol, zgSol);
      A134nunuee  = collision2D.A134_nunu_ee_I0(znuSol, zgSol);

      /* ----------------------------- */

      epsdb1      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd1 + rfd1) + A12nunuee*(rbfd1 + rfd1) + 
                        2*A134nunuee*pow(znuSol,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd1 + rfd1,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znuSol,6),0.5));

      epsdb2      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd2 + rfd2) + A12nunuee*(rbfd2 + rfd2) + 
                        2*A134nunuee*pow(znuSol,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd2 + rfd2,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znuSol,6),0.5));

      epsdb3      = pow(4*A123nunuee - 2*A12nunuee,-1)*
                        (-2*A123nunuee*(rbfd3 + rfd3) + A12nunuee*(rbfd3 + rfd3) + 
                        2*A134nunuee*pow(znuSol,3) - pow(pow(-2*A123nunuee + A12nunuee,2)*
                        pow(-rbfd3 + rfd3,2) + 4*
                        (-2*A123nunuee*A34nunuee + A12nunuee*A34nunuee + pow(A134nunuee,2))*
                        pow(znuSol,6),0.5));




      /* ------ print results ------- */
      cout << ", " << "epsDB1 = " << epsdb1 << ", " << "epsDB2 = " << epsdb2 << ", " << "epsDB3 = " << epsdb3 << endl;



      r_flavor_ini[0][0]      = (rfd1 + epsdb1)/pow(znuSol,3.);
      r_flavor_ini[1][1]      = (rfd2 + epsdb2)/pow(znuSol,3.);
      r_flavor_ini[2][2]      = (rfd3 + epsdb3)/pow(znuSol,3.);

      rb_flavor_ini[0][0]     = (rbfd1 + epsdb1)/pow(znuSol,3.);
      rb_flavor_ini[1][1]     = (rbfd2 + epsdb2)/pow(znuSol,3.);
      rb_flavor_ini[2][2]     = (rbfd3 + epsdb3)/pow(znuSol,3.);

      znu_ini= znuSol; zg_ini= zgSol;

}



// translate general matrix to the Gell-Mann matrix coefficients
void get_GellMann_coefficients_from_mat(complx mat[3][3], complx lvec [9][3][3], double (&coeff_GellMann) [9]){
      // initialize dummy 3x3 matrix to staore matrix multiplication results
      complx matmulaux[3][3];
      matop.mat33_free(matmulaux);

      coeff_GellMann[0]    = 1./3. * trace(mat);
      matop.matmul_33_33(mat,lvec[1],matmulaux); coeff_GellMann[1]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[2],matmulaux); coeff_GellMann[2]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[3],matmulaux); coeff_GellMann[3]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[4],matmulaux); coeff_GellMann[4]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[5],matmulaux); coeff_GellMann[5]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[6],matmulaux); coeff_GellMann[6]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[7],matmulaux); coeff_GellMann[7]    = 1./2. * trace(matmulaux);
      matop.matmul_33_33(mat,lvec[8],matmulaux); coeff_GellMann[8]    = 1./2. * trace(matmulaux);
}


/*
this function translates a set of Gell-Mann coefficients to the matrix representation
*/
void get_mat_from_GellMann_coeff(double coeffvec[9], complx lvec [9][3][3], complx (&mat) [3][3]){
      matop.mat33_free(mat);

      double c0, c1,c2,c3,c4,c5,c6,c7,c8;
      c0= coeffvec[0]; c1= coeffvec[1]; c2= coeffvec[2]; c3= coeffvec[3]; c4= coeffvec[4]; 
      c5= coeffvec[5]; c6= coeffvec[6]; c7= coeffvec[7]; c8= coeffvec[8];

      /* auxilary matrices needed for the multiplication */
      complx mat0_aux[3][3], mat1_aux[3][3], mat2_aux[3][3], mat3_aux[3][3], mat4_aux[3][3], mat5_aux[3][3], mat6_aux[3][3], mat7_aux[3][3], mat8_aux[3][3];
      complx matmul_aux[3][3];

      matop.mat33_free(mat0_aux); matop.mat33_free(mat1_aux); matop.mat33_free(mat2_aux); matop.mat33_free(mat3_aux); matop.mat33_free(mat4_aux);
      matop.mat33_free(mat5_aux); matop.mat33_free(mat6_aux); matop.mat33_free(mat7_aux); matop.mat33_free(mat8_aux); matop.mat33_free(matmul_aux);


      /* calculate r density matrix in IP representation */
      // first multiply all scalar solutions to the GellMann matrices
      matop.matmul_double_33(c0, lvec[0], mat0_aux);
      matop.matmul_double_33(c1, lvec[1], mat1_aux);
      matop.matmul_double_33(c2, lvec[2], mat2_aux);
      matop.matmul_double_33(c3, lvec[3], mat3_aux);
      matop.matmul_double_33(c4, lvec[4], mat4_aux);
      matop.matmul_double_33(c5, lvec[5], mat5_aux);
      matop.matmul_double_33(c6, lvec[6], mat6_aux);
      matop.matmul_double_33(c7, lvec[7], mat7_aux);
      matop.matmul_double_33(c8, lvec[8], mat8_aux);
      // now add all of the results together
      matop.matadd_33_33(mat0_aux, mat1_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat2_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat3_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat4_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat5_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat6_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat7_aux, matmul_aux);
      matop.matadd_33_33(matmul_aux, mat8_aux, mat);


}



/*
function translating the result of the ODE into (n-\bar{n})_\alpha * normalization
the normalization factor encodes if you wish to have (n-\bar{n})/T^3, (n-\bar{n})/n_\gamma, etc..
INPUT: r and rbar matrix in the FLAVOR BASIS.
*/
array<double, 3> get_dn_alpha_from_r_rbar_fb(double normalize, complx r_fb[3][3], complx rbar_fb[3][3]){

      double nnbar11, nnbar22, nnbar33;

      nnbar11 = r_fb[0][0].real() - rbar_fb[0][0].real();
      nnbar22 = r_fb[1][1].real() - rbar_fb[1][1].real();
      nnbar33 = r_fb[2][2].real() - rbar_fb[2][2].real();

      return {normalize*nnbar11, normalize*nnbar22, normalize*nnbar33};
}


/*
function translating the result of the ODE into (n+\bar{n})_\alpha * normalization
the normalization factor encodes if you wish to have (n+\bar{n})/T^3, (n+\bar{n})/n_\gamma, etc..
INPUT: r and rbar matrix in the FLAVOR BASIS.
*/
array<double, 3> get_np_alpha_from_r_rbar_fb(double normalize, complx r_fb[3][3], complx rbar_fb[3][3]){

      double np11, np22, np33;

      np11 = r_fb[0][0].real() + rbar_fb[0][0].real();
      np22 = r_fb[1][1].real() + rbar_fb[1][1].real();
      np33 = r_fb[2][2].real() + rbar_fb[2][2].real();

      return {normalize*np11, normalize*np22, normalize*np33};
}


/*
***************************************
---------------------------------------
Define full system of kinetic equation
---------------------------------------
***************************************
*/
// right hand side of the full kinetic eq for the GellMann coefficients in the flavor basis
int ode_fun(double t, const double y[], double f[], void *params){
      std::vector<double> params_var= *(std::vector<double> *)params;
      double dm21, dm31, theta12, theta13, theta23, zon, adiabatic; 
      dm21= params_var[0]; dm31= params_var[1];
      theta12= params_var[2]; theta13= params_var[3]; theta23= params_var[4]; 
      zon= params_var[5]; adiabatic= params_var[6];


      // define the sin(Weinberg angle)
      double sw;
      sw   = sin(_ThetaW_);

      // set neutrino and photon temperature
      double znu,zg, dznu;
      znu= y[18]; zg= y[19]; dznu= f[18];
      double Tg, Tnu;
      Tg          = _Tref_*zg/t;
      Tnu         = _Tref_*znu/t;


      double znu_lhs;
      // need to correct for initial numerical instability (fluctuation in the system due to dynamical fix point convergence)
      // this hard coded number can vary with the GSL ODE solver options like hstart, length_vec_ode_sol etc.
      // but with the initial condition procedure described in Sec 2B of 2502.14960 this instability should not be present
      // regard it as a backup just in case..
      znu_lhs    = zon * 3.*dznu/znu;
      if ( abs(1 - znu_lhs) > 25. ){
            znu_lhs     = 0;
      }


      /* ------- set the rates ------- */
      // collision term 0. moment
      // 1D in znu
      double A12nunuee;        
      double A12nuan, A123nuan;      
      // 2D in zg, x
      double A34nunuee;     
      // 2D in znu, zg
      double A12nuenue, A124nuenue;       
      double A123nunuee, A134nunuee;    
      /* ------ */
      // transfer rates 1. moment
      double A12nuenueI1, A34nuenueI1;
      double A12nunueeI1, A34nunueeI1;


      // collision term 0. moment

      // 1D rates in znu
      A12nunuee   = collision1D.A12_nunu_ee_I0(znu);

      A12nuan     = collision1D.A12_nuan_I0(znu);
      A123nuan    = collision1D.A123_nuan_I0(znu);

      // 2D rates in zg, x
      //stef, for massive rate use _Tref_/t, for massless use just _Tref_ --> but then the Neff normalization (NeffNorm) has to be adjusted
      A34nunuee   = collision2D_me.A34_nunu_ee_I0_me(zg, _Tref_/t);

      // 2D rates in znu, zg
      A12nuenue   = collision2D.A12_nue_nue_I0(znu, zg);
      A124nuenue  = collision2D.A124_nue_nue_I0(znu, zg);
      A123nunuee  = collision2D.A123_nunu_ee_I0(znu, zg);
      A134nunuee  = collision2D.A134_nunu_ee_I0(znu, zg);


      /* --------- */

      // transfer rates 1. moment; do not include NLO FD correction !!
      
      A12nuenueI1       = transfer2D.A12_nue_nue_I1(znu, zg);
      A34nuenueI1       = transfer2D.A34_nue_nue_I1(znu,zg);

      A12nunueeI1       = transfer1D.Aij_nunu_ee_I1(znu);
      A34nunueeI1       = transfer1D.Aij_nunu_ee_I1(zg);

      
      /* ----------------------------- */



      // translate the Gell-Mann coefficients to the sum of the diagonal elements of the r & rb matrix (Trace)
      double Trace_r_rb;
      Trace_r_rb  = 3. * (y[0] + y[9]);

      // define the "reduced" planck mass, this assumes massless electrons
      double mpl;
      mpl         = _MPL_/sqrt(11./2. * pow(zg,4.) + 7./8. * Trace_r_rb * pow(znu,4.));


      // define rescaled mass splitting entering the kinetic eq.
      double d21, d31;
      d21         = 1./pow(znu,1.) * pow(t/_Tref_,2.) * mpl/_Tref_ * pow(_PI_,2.)/(36.*_Zeta3_) * dm21;
      d31         = 1./pow(znu,1.) * pow(t/_Tref_,2.) * mpl/_Tref_ * pow(_PI_,2.)/(36.*_Zeta3_) * dm31; 


      // V_th terms; note that the interpolation is a function of T instead of x directly
      double e1, e2;
      e1          = mpl/_Tref_ * (2.*sqrt(2)*_GF_/pow(_MW_,2.)) * ( (7.*pow(_PI_,4.))/(180.*_Zeta3_) * pow(znu,1.) ) * ( density.rhoe(Tg) + density.pe(Tg) );
      e2          = mpl/_Tref_ * (2.*sqrt(2)*_GF_/pow(_MW_,2.)) * ( (7.*pow(_PI_,4.))/(180.*_Zeta3_) * pow(znu,1.) ) * ( density.rhomu(Tg) + density.pmu(Tg) );


      // self scattering potential
      double vs;
      vs          = adiabatic * pow(znu,3.) * mpl/_Tref_ * 3.*_Zeta3_/(2. * pow(_PI_,2.) * sqrt(2)) * _GF_ * pow(_Tref_/t,2.);



      // overall pre factor (normalization) for collision term
      double c;
      c           = 1./pow(znu,3.) * 1./(3.*_Zeta3_/2.) * mpl/_Tref_ * pow(_GF_,2.) * pow(_Tref_/t,4.);



      // dummy variables for the collision term
      // c1, c2, c3 isolate nuenue, nunuee, nuan; qubicterms can set to zero epsilon^3 terms
      int c1, c2, c3, qubicTerms;
      c1= 1; c2= 1; c3= 1; qubicTerms= 1;



      /* ---------------- 0. moment evolution ------------------ */

      /* --------- r ---------- */
      f[0]  =     (c*c1*pow(_PI_,-3)*(3*A34nunuee*(3 - 4*pow(sw,2) + 24*pow(sw,4)) - 
                  A134nunuee*(9*y[0] + 9*y[9] + 72*pow(sw,4)*(y[0] + y[9]) + 
                  4*pow(sw,2)*(-3*y[0] + 6*y[3] + 2*pow(3,0.5)*y[8] - 3*y[9] + 6*y[12] + 2*pow(3,0.5)*y[17])) + 
                  (2*A123nunuee - A12nunuee)*(9*y[0]*y[9] - 
                  6*(y[1]*y[10] + y[2]*y[11] - y[3]*y[12] + y[4]*y[13] + y[5]*y[14] - y[6]*y[15] - y[7]*y[16] - 
                  y[8]*y[17]) + 4*pow(sw,2)*(y[0]*(-3*y[9] + 6*y[12] + 2*pow(3,0.5)*y[17]) + 
                  2*(3*y[3]*y[9] + pow(3,0.5)*y[8]*y[9] + pow(3,0.5)*y[8]*y[12] - 3*y[6]*y[15] - 
                  3*y[7]*y[16] + pow(3,0.5)*y[3]*y[17] - 2*y[8]*y[17])) + 
                  24*pow(sw,4)*(3*y[0]*y[9] + 2*(y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16] + y[8]*y[17])))))/144.
                  -
                  znu_lhs*y[0];

      f[1]  =     (e1 - e2)*y[2] - ((-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*(-3 + cos(2*theta23))*y[2])/8. + 
                  ((d21 - 2*d31 + 3*d21*cos(2*theta12))*pow(cos(theta23),2)*y[2])/4. + 
                  cos(theta23)*(sin(theta23)*(-2*d21*cos(theta12)*sin(theta12)*sin(theta13)*y[2] - 
                  d21*pow(cos(theta12),2)*y[5] + d31*pow(cos(theta13),2)*y[5] + 
                  d21*pow(sin(theta12),2)*pow(sin(theta13),2)*y[5]) + 
                  cos(theta13)*(d31 - d21*pow(sin(theta12),2))*sin(theta13)*y[7]) - 
                  d21*cos(theta12)*sin(theta12)*(cos(2*theta23)*sin(theta13)*y[5] + cos(theta13)*sin(theta23)*y[7]) + 
                  vs*(-2*y[3]*y[11] + 2*y[2]*y[12] - y[7]*y[13] + y[6]*y[14] - y[5]*y[15] + y[4]*y[16])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(6*A134nunuee*c1*y[1] - 12*A124nuenue*c2*y[1] + 
                  12*A12nuenue*c2*y[1] + 48*A134nunuee*c1*pow(sw,4)*y[1] - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[1] + 48*A123nuan*c3*qubicTerms*pow(y[10],2)*y[1] - 
                  12*A123nunuee*c1*y[1]*y[9] + 6*A12nunuee*c1*y[1]*y[9] + 72*A12nuan*c3*y[1]*y[9] - 
                  96*A123nunuee*c1*pow(sw,4)*y[1]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[1]*y[9] + 
                  2*c1*(-1 + 8*pow(sw,4))*(3*A134nunuee - (2*A123nunuee - A12nunuee)*(3*y[0] + pow(3,0.5)*y[8]))*
                  y[10] + 48*A123nuan*c3*qubicTerms*y[1]*y[2]*y[11] - 48*A123nunuee*c1*pow(sw,2)*y[1]*y[12] + 
                  24*A12nunuee*c1*pow(sw,2)*y[1]*y[12] + 48*A123nuan*c3*qubicTerms*y[1]*y[3]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[4]*y[13] + 6*A123nunuee*c1*y[6]*y[13] - 
                  3*A12nunuee*c1*y[6]*y[13] + 36*A12nuan*c3*y[6]*y[13] - 48*A123nunuee*c1*pow(sw,4)*y[6]*y[13] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[13] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[13] + 48*A123nuan*c3*qubicTerms*y[1]*y[5]*y[14] + 
                  6*A123nunuee*c1*y[7]*y[14] - 3*A12nunuee*c1*y[7]*y[14] + 36*A12nuan*c3*y[7]*y[14] - 
                  48*A123nunuee*c1*pow(sw,4)*y[7]*y[14] + 24*A12nunuee*c1*pow(sw,4)*y[7]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[14] - 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[14] - 
                  6*A123nunuee*c1*y[4]*y[15] + 3*A12nunuee*c1*y[4]*y[15] + 36*A12nuan*c3*y[4]*y[15] + 
                  24*A123nunuee*c1*pow(sw,2)*y[4]*y[15] - 12*A12nunuee*c1*pow(sw,2)*y[4]*y[15] - 
                  48*A123nunuee*c1*pow(sw,4)*y[4]*y[15] + 24*A12nunuee*c1*pow(sw,4)*y[4]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[15] + 48*A123nuan*c3*qubicTerms*y[1]*y[6]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[15] - 6*A123nunuee*c1*y[5]*y[16] + 
                  3*A12nunuee*c1*y[5]*y[16] + 36*A12nuan*c3*y[5]*y[16] + 24*A123nunuee*c1*pow(sw,2)*y[5]*y[16] - 
                  12*A12nunuee*c1*pow(sw,2)*y[5]*y[16] - 48*A123nunuee*c1*pow(sw,4)*y[5]*y[16] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[16] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[16] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[7]*y[16] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[16] - 
                  4*A123nunuee*c1*pow(3,0.5)*y[1]*y[17] + 2*A12nunuee*c1*pow(3,0.5)*y[1]*y[17] + 
                  24*A12nuan*c3*pow(3,0.5)*y[1]*y[17] - 32*A123nunuee*c1*pow(3,0.5)*pow(sw,4)*y[1]*y[17] + 
                  16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[1]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[1]*y[17] + 48*A123nuan*c3*qubicTerms*y[1]*y[8]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[9]*y[17] + 
                  24*c3*y[10]*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[1],2) + y[2]*y[11] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17]))))
                  -
                  znu_lhs*y[1];

      f[2]  =     (-8*e1*y[1] + 8*e2*y[1] + 4*d21*sin(2*theta12)*sin(theta13)*(sin(2*theta23)*y[1] + cos(2*theta23)*y[4]) - 
                  2*(d21 - 2*d31 + 3*d21*cos(2*theta12))*cos(theta23)*(cos(theta23)*y[1] - sin(theta23)*y[4]) + 
                  (-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*
                  ((-3 + cos(2*theta23))*y[1] - sin(2*theta23)*y[4]) + 
                  2*(-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(-2*sin(theta23)*y[3] + cos(theta23)*y[6]) - 
                  4*d21*cos(theta13)*sin(2*theta12)*(2*cos(theta23)*y[3] + sin(theta23)*y[6]) + 
                  8*vs*(2*y[3]*y[10] - 2*y[1]*y[12] - y[6]*y[13] - y[7]*y[14] + y[4]*y[15] + y[5]*y[16]))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(6*A134nunuee*c1*y[2] - 12*A124nuenue*c2*y[2] + 
                  12*A12nuenue*c2*y[2] + 48*A134nunuee*c1*pow(sw,4)*y[2] - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[2] + 48*A123nuan*c3*qubicTerms*pow(y[11],2)*y[2] - 
                  12*A123nunuee*c1*y[2]*y[9] + 6*A12nunuee*c1*y[2]*y[9] + 72*A12nuan*c3*y[2]*y[9] - 
                  96*A123nunuee*c1*pow(sw,4)*y[2]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[2]*y[9] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[2]*y[10] + 
                  2*c1*(-1 + 8*pow(sw,4))*(3*A134nunuee - (2*A123nunuee - A12nunuee)*(3*y[0] + pow(3,0.5)*y[8]))*
                  y[11] - 48*A123nunuee*c1*pow(sw,2)*y[2]*y[12] + 24*A12nunuee*c1*pow(sw,2)*y[2]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[3]*y[12] + 48*A123nuan*c3*qubicTerms*y[2]*y[4]*y[13] - 
                  6*A123nunuee*c1*y[7]*y[13] + 3*A12nunuee*c1*y[7]*y[13] - 36*A12nuan*c3*y[7]*y[13] + 
                  48*A123nunuee*c1*pow(sw,4)*y[7]*y[13] - 24*A12nunuee*c1*pow(sw,4)*y[7]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[13] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[13] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[5]*y[14] + 6*A123nunuee*c1*y[6]*y[14] - 
                  3*A12nunuee*c1*y[6]*y[14] + 36*A12nuan*c3*y[6]*y[14] - 48*A123nunuee*c1*pow(sw,4)*y[6]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[14] - 6*A123nunuee*c1*y[5]*y[15] + 
                  3*A12nunuee*c1*y[5]*y[15] + 36*A12nuan*c3*y[5]*y[15] + 24*A123nunuee*c1*pow(sw,2)*y[5]*y[15] - 
                  12*A12nunuee*c1*pow(sw,2)*y[5]*y[15] - 48*A123nunuee*c1*pow(sw,4)*y[5]*y[15] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[15] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[15] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[6]*y[15] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[15] + 
                  6*A123nunuee*c1*y[4]*y[16] - 3*A12nunuee*c1*y[4]*y[16] - 36*A12nuan*c3*y[4]*y[16] - 
                  24*A123nunuee*c1*pow(sw,2)*y[4]*y[16] + 12*A12nunuee*c1*pow(sw,2)*y[4]*y[16] + 
                  48*A123nunuee*c1*pow(sw,4)*y[4]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[4]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[16] + 48*A123nuan*c3*qubicTerms*y[2]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[16] - 4*A123nunuee*c1*pow(3,0.5)*y[2]*y[17] + 
                  2*A12nunuee*c1*pow(3,0.5)*y[2]*y[17] + 24*A12nuan*c3*pow(3,0.5)*y[2]*y[17] - 
                  32*A123nunuee*c1*pow(3,0.5)*pow(sw,4)*y[2]*y[17] + 
                  16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[2]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[2]*y[17] + 48*A123nuan*c3*qubicTerms*y[2]*y[8]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[9]*y[17] + 
                  24*c3*y[11]*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[2],2) + y[1]*y[10] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17]))))
                  -
                  znu_lhs*y[2];

      f[3]  =     cos(theta13)*(-(d21*cos(theta12)*sin(theta12)*sin(theta23)*y[5]) + 
                  cos(theta23)*(d21*sin(2*theta12)*y[2] + (d31 - d21*pow(sin(theta12),2))*sin(theta13)*y[5])) + 
                  d21*cos(theta12)*cos(2*theta23)*sin(theta12)*sin(theta13)*y[7] - 
                  d31*cos(theta23)*pow(cos(theta13),2)*sin(theta23)*y[7] + 
                  sin(theta23)*((d31 - d21*pow(sin(theta12),2))*sin(2*theta13)*y[2] + 
                  d21*cos(theta23)*(pow(cos(theta12),2) - pow(sin(theta12),2)*pow(sin(theta13),2))*y[7]) + 
                  vs*(-2*y[2]*y[10] + 2*y[1]*y[11] - y[5]*y[13] + y[4]*y[14] + y[7]*y[15] - y[6]*y[16])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-24*A34nunuee*c1*pow(sw,2) - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[3] + 48*A123nuan*c3*qubicTerms*pow(y[12],2)*y[3] + 
                  24*A12nunuee*c1*pow(sw,2)*y[0]*y[9] + 6*A12nunuee*c1*y[3]*y[9] + 72*A12nuan*c3*y[3]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[3]*y[9] + 8*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[8]*y[9] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[3]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[3]*y[11] + 
                  2*A12nunuee*c1*(3*(1 + 8*pow(sw,4))*y[0] + 12*pow(sw,2)*y[3] + pow(3,0.5)*(1 + 8*pow(sw,4))*y[8])*
                  y[12] - 3*A12nunuee*c1*y[4]*y[13] + 36*A12nuan*c3*y[4]*y[13] + 
                  24*A12nunuee*c1*pow(sw,4)*y[4]*y[13] - 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[13] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[4]*y[13] - 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[13] - 
                  3*A12nunuee*c1*y[5]*y[14] + 36*A12nuan*c3*y[5]*y[14] + 24*A12nunuee*c1*pow(sw,4)*y[5]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[14] + 48*A123nuan*c3*qubicTerms*y[3]*y[5]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[14] - 3*A12nunuee*c1*y[6]*y[15] - 36*A12nuan*c3*y[6]*y[15] + 
                  12*A12nunuee*c1*pow(sw,2)*y[6]*y[15] - 24*A12nunuee*c1*pow(sw,4)*y[6]*y[15] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[15] + 48*A123nuan*c3*qubicTerms*y[3]*y[6]*y[15] + 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[15] - 3*A12nunuee*c1*y[7]*y[16] - 36*A12nuan*c3*y[7]*y[16] + 
                  12*A12nunuee*c1*pow(sw,2)*y[7]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[16] + 48*A123nuan*c3*qubicTerms*y[3]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[16] + 8*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[0]*y[17] + 
                  2*A12nunuee*c1*pow(3,0.5)*y[3]*y[17] + 24*A12nuan*c3*pow(3,0.5)*y[3]*y[17] + 
                  16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[3]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[3]*y[17] + 8*A12nunuee*c1*pow(sw,2)*y[8]*y[17] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[8]*y[17] - 24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[9]*y[17] + 
                  2*A134nunuee*c1*(3*(y[3] + y[12]) + 24*pow(sw,4)*(y[3] + y[12]) + 
                  4*pow(sw,2)*(3*y[0] + pow(3,0.5)*y[8] + 3*y[9] + pow(3,0.5)*y[17])) + 
                  24*c3*y[12]*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[3],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17])) - 
                  2*A123nunuee*c1*(2*(1 + 8*pow(sw,4))*(3*y[0] + pow(3,0.5)*y[8])*y[12] - 
                  3*(y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + y[7]*y[16]) + 
                  2*y[3]*(3*(1 + 8*pow(sw,4))*y[9] + 12*pow(sw,2)*y[12] + pow(3,0.5)*(1 + 8*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(3*y[6]*y[15] + 3*y[7]*y[16] + 
                  6*pow(sw,2)*(y[4]*y[13] + y[5]*y[14] - y[6]*y[15] - y[7]*y[16]) + 
                  2*y[8]*(pow(3,0.5)*y[9] + y[17]) + 2*y[0]*(3*y[9] + pow(3,0.5)*y[17])))))
                  -
                  znu_lhs*y[3];

      f[4]  =     (-8*d21*cos(2*theta23)*sin(2*theta12)*sin(theta13)*y[2] + (d21 - 2*d31 + 16*e1)*y[5] - 
                  3*d21*cos(2*theta23)*y[5] + 6*d31*cos(2*theta23)*y[5] - 
                  2*(d21 - 2*d31)*cos(2*theta13)*(3 + cos(2*theta23))*y[5] + 
                  2*d21*cos(2*theta12)*((-3 + cos(2*theta13))*cos(2*theta23) + 6*pow(cos(theta13),2))*y[5] + 
                  2*d21*pow(cos(theta23),2)*y[5] - 4*d31*pow(cos(theta23),2)*y[5] + 
                  8*sin(2*theta23)*(-(d21*pow(cos(theta12),2)*y[2]) + d31*pow(cos(theta13),2)*y[2] + 
                  d21*sin(theta13)*(pow(sin(theta12),2)*sin(theta13)*y[2] + sin(2*theta12)*y[5])) - 
                  8*d21*cos(theta13)*cos(theta23)*sin(2*theta12)*y[7] - 8*d31*sin(2*theta13)*sin(theta23)*y[7] + 
                  8*d21*pow(sin(theta12),2)*sin(2*theta13)*sin(theta23)*y[7] + 
                  16*vs*(y[7]*y[10] + y[6]*y[11] + y[5]*y[12] - (y[3] + pow(3,0.5)*y[8])*y[14] - y[2]*y[15] - 
                  y[1]*y[16] + pow(3,0.5)*y[5]*y[17]))/16.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-12*A124nuenue*c2*y[4] + 12*A12nuenue*c2*y[4] - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[4] + 48*A123nuan*c3*qubicTerms*pow(y[13],2)*y[4] - 
                  12*A123nunuee*c1*y[4]*y[9] + 6*A12nunuee*c1*y[4]*y[9] + 72*A12nuan*c3*y[4]*y[9] - 
                  96*A123nunuee*c1*pow(sw,4)*y[4]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[4]*y[9] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[4]*y[10] + 6*A123nunuee*c1*y[6]*y[10] - 
                  3*A12nunuee*c1*y[6]*y[10] + 36*A12nuan*c3*y[6]*y[10] - 48*A123nunuee*c1*pow(sw,4)*y[6]*y[10] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[10] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[4]*y[11] - 
                  6*A123nunuee*c1*y[7]*y[11] + 3*A12nunuee*c1*y[7]*y[11] - 36*A12nuan*c3*y[7]*y[11] + 
                  48*A123nunuee*c1*pow(sw,4)*y[7]*y[11] - 24*A12nunuee*c1*pow(sw,4)*y[7]*y[11] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[11] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[11] - 
                  6*A123nunuee*c1*y[4]*y[12] + 3*A12nunuee*c1*y[4]*y[12] + 36*A12nuan*c3*y[4]*y[12] - 
                  24*A123nunuee*c1*pow(sw,2)*y[4]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[4]*y[12] - 
                  48*A123nunuee*c1*pow(sw,4)*y[4]*y[12] + 24*A12nunuee*c1*pow(sw,4)*y[4]*y[12] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[12] + 48*A123nuan*c3*qubicTerms*y[3]*y[4]*y[12] - 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[12] + 
                  6*A134nunuee*c1*(y[4] - y[13] + 8*pow(sw,4)*(y[4] + y[13])) + 
                  48*A123nuan*c3*qubicTerms*y[4]*y[5]*y[14] - 6*A123nunuee*c1*y[1]*y[15] + 
                  3*A12nunuee*c1*y[1]*y[15] + 36*A12nuan*c3*y[1]*y[15] + 24*A123nunuee*c1*pow(sw,2)*y[1]*y[15] - 
                  12*A12nunuee*c1*pow(sw,2)*y[1]*y[15] - 48*A123nunuee*c1*pow(sw,4)*y[1]*y[15] + 
                  24*A12nunuee*c1*pow(sw,4)*y[1]*y[15] - 36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[15] + 
                  48*A123nuan*c3*qubicTerms*y[4]*y[6]*y[15] - 36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[15] + 
                  6*A123nunuee*c1*y[2]*y[16] - 3*A12nunuee*c1*y[2]*y[16] - 36*A12nuan*c3*y[2]*y[16] - 
                  24*A123nunuee*c1*pow(sw,2)*y[2]*y[16] + 12*A12nunuee*c1*pow(sw,2)*y[2]*y[16] + 
                  48*A123nunuee*c1*pow(sw,4)*y[2]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[2]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[16] + 48*A123nuan*c3*qubicTerms*y[4]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[16] + 2*A123nunuee*c1*pow(3,0.5)*y[4]*y[17] - 
                  A12nunuee*c1*pow(3,0.5)*y[4]*y[17] - 12*A12nuan*c3*pow(3,0.5)*y[4]*y[17] - 
                  24*A123nunuee*c1*pow(3,0.5)*pow(sw,2)*y[4]*y[17] + 
                  12*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[4]*y[17] + 
                  16*A123nunuee*c1*pow(3,0.5)*pow(sw,4)*y[4]*y[17] - 
                  8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[4]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[4]*y[17] + 48*A123nuan*c3*qubicTerms*y[4]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[9]*y[17] + 
                  y[13]*(-2*A123nunuee*c1*(-1 + 8*pow(sw,4))*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A12nunuee*c1*(-1 + 8*pow(sw,4))*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + y[0]*(-3*y[3] + pow(3,0.5)*y[8]) - 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*(pow(y[4],2) + y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + 
                  y[5]*y[14] + y[6]*y[15] + y[7]*y[16] + y[8]*y[17]))))))
                  -
                  znu_lhs*y[4];

      f[5]  =     (-((d21 - 2*d31 + 8*e1)*y[4]) + (d21 - 2*d31 + 3*d21*cos(2*theta12))*
                  (sin(2*theta23)*y[1] + cos(2*theta23)*y[4]) + 
                  (d21 - 2*d31)*cos(2*theta13)*(sin(2*theta23)*y[1] + (3 + cos(2*theta23))*y[4]) + 
                  4*d21*sin(2*theta12)*sin(theta13)*(cos(2*theta23)*y[1] - sin(2*theta23)*y[4]) + 
                  2*(d21 - 2*d31)*sin(2*theta13)*(-(sin(theta23)*y[6]) + cos(theta23)*(y[3] + pow(3,0.5)*y[8])) + 
                  4*d21*cos(theta13)*sin(2*theta12)*(cos(theta23)*y[6] + sin(theta23)*(y[3] + pow(3,0.5)*y[8])) - 
                  d21*cos(2*theta12)*(6*pow(cos(theta13),2)*y[4] + 
                  cos(2*theta13)*(sin(2*theta23)*y[1] + cos(2*theta23)*y[4]) + 
                  2*sin(2*theta13)*(-(sin(theta23)*y[6]) + cos(theta23)*(y[3] + pow(3,0.5)*y[8]))) + 
                  8*vs*(-(y[6]*y[10]) + y[7]*y[11] + y[3]*y[13] + pow(3,0.5)*y[8]*y[13] + y[1]*y[15] - y[2]*y[16] - 
                  y[4]*(y[12] + pow(3,0.5)*y[17])))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-12*A124nuenue*c2*y[5] + 12*A12nuenue*c2*y[5] - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[5] + 48*A123nuan*c3*qubicTerms*pow(y[14],2)*y[5] - 
                  12*A123nunuee*c1*y[5]*y[9] + 6*A12nunuee*c1*y[5]*y[9] + 72*A12nuan*c3*y[5]*y[9] - 
                  96*A123nunuee*c1*pow(sw,4)*y[5]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[5]*y[9] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[5]*y[10] + 6*A123nunuee*c1*y[7]*y[10] - 
                  3*A12nunuee*c1*y[7]*y[10] + 36*A12nuan*c3*y[7]*y[10] - 48*A123nunuee*c1*pow(sw,4)*y[7]*y[10] + 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[10] - 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[5]*y[11] + 
                  6*A123nunuee*c1*y[6]*y[11] - 3*A12nunuee*c1*y[6]*y[11] + 36*A12nuan*c3*y[6]*y[11] - 
                  48*A123nunuee*c1*pow(sw,4)*y[6]*y[11] + 24*A12nunuee*c1*pow(sw,4)*y[6]*y[11] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[11] - 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[11] - 
                  6*A123nunuee*c1*y[5]*y[12] + 3*A12nunuee*c1*y[5]*y[12] + 36*A12nuan*c3*y[5]*y[12] - 
                  24*A123nunuee*c1*pow(sw,2)*y[5]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[5]*y[12] - 
                  48*A123nunuee*c1*pow(sw,4)*y[5]*y[12] + 24*A12nunuee*c1*pow(sw,4)*y[5]*y[12] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[12] + 48*A123nuan*c3*qubicTerms*y[3]*y[5]*y[12] - 
                  36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[12] + 48*A123nuan*c3*qubicTerms*y[4]*y[5]*y[13] + 
                  6*A134nunuee*c1*(y[5] - y[14] + 8*pow(sw,4)*(y[5] + y[14])) - 6*A123nunuee*c1*y[2]*y[15] + 
                  3*A12nunuee*c1*y[2]*y[15] + 36*A12nuan*c3*y[2]*y[15] + 24*A123nunuee*c1*pow(sw,2)*y[2]*y[15] - 
                  12*A12nunuee*c1*pow(sw,2)*y[2]*y[15] - 48*A123nunuee*c1*pow(sw,4)*y[2]*y[15] + 
                  24*A12nunuee*c1*pow(sw,4)*y[2]*y[15] - 36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[15] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[6]*y[15] - 36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[15] - 
                  6*A123nunuee*c1*y[1]*y[16] + 3*A12nunuee*c1*y[1]*y[16] + 36*A12nuan*c3*y[1]*y[16] + 
                  24*A123nunuee*c1*pow(sw,2)*y[1]*y[16] - 12*A12nunuee*c1*pow(sw,2)*y[1]*y[16] - 
                  48*A123nunuee*c1*pow(sw,4)*y[1]*y[16] + 24*A12nunuee*c1*pow(sw,4)*y[1]*y[16] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[16] + 48*A123nuan*c3*qubicTerms*y[5]*y[7]*y[16] - 
                  36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[16] + 2*A123nunuee*c1*pow(3,0.5)*y[5]*y[17] - 
                  A12nunuee*c1*pow(3,0.5)*y[5]*y[17] - 12*A12nuan*c3*pow(3,0.5)*y[5]*y[17] - 
                  24*A123nunuee*c1*pow(3,0.5)*pow(sw,2)*y[5]*y[17] + 
                  12*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[5]*y[17] + 
                  16*A123nunuee*c1*pow(3,0.5)*pow(sw,4)*y[5]*y[17] - 
                  8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[5]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[5]*y[17] + 48*A123nuan*c3*qubicTerms*y[5]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[9]*y[17] + 
                  y[14]*(-2*A123nunuee*c1*(-1 + 8*pow(sw,4))*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A12nunuee*c1*(-1 + 8*pow(sw,4))*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + y[0]*(-3*y[3] + pow(3,0.5)*y[8]) - 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*(pow(y[5],2) + y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + 
                  y[4]*y[13] + y[6]*y[15] + y[7]*y[16] + y[8]*y[17]))))))
                  -
                  znu_lhs*y[5];

      f[6]  =     -0.25*((-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(cos(theta23)*y[2] + sin(theta23)*y[5])) + 
                  e2*y[7] + (cos(2*theta23)*(d21*cos(2*theta12)*(-3 + cos(2*theta13)) - 
                  2*(d21 - 2*d31)*pow(cos(theta13),2))*y[7])/4. + 
                  d21*cos(theta12)*sin(theta12)*(cos(theta13)*(sin(theta23)*y[2] - cos(theta23)*y[5]) + 
                  4*cos(theta23)*sin(theta13)*sin(theta23)*y[7]) + 
                  vs*(y[5]*y[10] - y[4]*y[11] - y[7]*y[12] + y[2]*y[13] - y[1]*y[14] + y[3]*y[16] - 
                  pow(3,0.5)*y[8]*y[16] + pow(3,0.5)*y[7]*y[17])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[6] + 
                  48*A123nuan*c3*qubicTerms*pow(y[15],2)*y[6] + 6*A12nunuee*c1*y[6]*y[9] + 72*A12nuan*c3*y[6]*y[9] - 
                  24*A12nunuee*c1*pow(sw,2)*y[6]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[6]*y[9] - 
                  3*A12nunuee*c1*y[4]*y[10] + 36*A12nuan*c3*y[4]*y[10] + 24*A12nunuee*c1*pow(sw,4)*y[4]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[10] + 48*A123nuan*c3*qubicTerms*y[1]*y[6]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[10] - 3*A12nunuee*c1*y[5]*y[11] + 36*A12nuan*c3*y[5]*y[11] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[11] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[11] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[6]*y[11] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[11] - 
                  3*A12nunuee*c1*y[6]*y[12] - 36*A12nuan*c3*y[6]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[6]*y[12] - 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[12] + 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[6]*y[12] + 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[12] - 
                  3*A12nunuee*c1*y[1]*y[13] + 36*A12nuan*c3*y[1]*y[13] + 24*A12nunuee*c1*pow(sw,4)*y[1]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[13] + 48*A123nuan*c3*qubicTerms*y[4]*y[6]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[13] - 3*A12nunuee*c1*y[2]*y[14] + 36*A12nuan*c3*y[2]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[2]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[6]*y[14] - 36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[14] + 
                  6*A134nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(y[6] + y[15]) + 
                  48*A123nuan*c3*qubicTerms*y[6]*y[7]*y[16] - A12nunuee*c1*pow(3,0.5)*y[6]*y[17] - 
                  12*A12nuan*c3*pow(3,0.5)*y[6]*y[17] + 4*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[6]*y[17] - 
                  8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[6]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[6]*y[17] + 48*A123nuan*c3*qubicTerms*y[6]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[9]*y[17] + 
                  2*A123nunuee*c1*(-3*(-1 + 8*pow(sw,4))*(y[4]*y[10] + y[5]*y[11] + y[1]*y[13] + y[2]*y[14]) - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8])*y[15] - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*y[6]*(6*y[9] - 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[15]*(A12nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + 3*y[0]*y[3] + pow(3,0.5)*y[0]*y[8] + 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*y[3]*y[12] + 
                  4*(pow(y[6],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[7]*y[16]) + 
                  4*y[8]*y[17])))))
                  -
                  znu_lhs*y[6];

      f[7]  =     (4*d21*cos(theta13)*sin(2*theta12)*(sin(theta23)*y[1] + cos(theta23)*y[4]) + 
                  2*(-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(-(cos(theta23)*y[1]) + sin(theta23)*y[4]) - 
                  8*e2*y[6] - 4*d21*sin(2*theta12)*sin(theta13)*
                  (2*sin(2*theta23)*y[6] + cos(2*theta23)*(y[3] - pow(3,0.5)*y[8])) + 
                  (d21 - 2*d31 + 3*d21*cos(2*theta12))*(2*cos(2*theta23)*y[6] + 
                  sin(2*theta23)*(-y[3] + pow(3,0.5)*y[8])) - 
                  (-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*
                  (2*cos(2*theta23)*y[6] + sin(2*theta23)*(-y[3] + pow(3,0.5)*y[8])) - 
                  8*vs*(y[4]*y[10] + y[5]*y[11] - y[6]*y[12] - y[1]*y[13] - y[2]*y[14] + y[3]*y[15] - 
                  pow(3,0.5)*y[8]*y[15] + pow(3,0.5)*y[6]*y[17]))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[7] + 
                  48*A123nuan*c3*qubicTerms*pow(y[16],2)*y[7] + 6*A12nunuee*c1*y[7]*y[9] + 72*A12nuan*c3*y[7]*y[9] - 
                  24*A12nunuee*c1*pow(sw,2)*y[7]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[7]*y[9] - 
                  3*A12nunuee*c1*y[5]*y[10] + 36*A12nuan*c3*y[5]*y[10] + 24*A12nunuee*c1*pow(sw,4)*y[5]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[10] + 48*A123nuan*c3*qubicTerms*y[1]*y[7]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[10] + 3*A12nunuee*c1*y[4]*y[11] - 36*A12nuan*c3*y[4]*y[11] - 
                  24*A12nunuee*c1*pow(sw,4)*y[4]*y[11] + 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[11] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[7]*y[11] + 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[11] - 
                  3*A12nunuee*c1*y[7]*y[12] - 36*A12nuan*c3*y[7]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[7]*y[12] - 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[12] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[7]*y[12] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[12] + 
                  3*A12nunuee*c1*y[2]*y[13] - 36*A12nuan*c3*y[2]*y[13] - 24*A12nunuee*c1*pow(sw,4)*y[2]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[13] + 48*A123nuan*c3*qubicTerms*y[4]*y[7]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[13] - 3*A12nunuee*c1*y[1]*y[14] + 36*A12nuan*c3*y[1]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[1]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[7]*y[14] - 36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[6]*y[7]*y[15] + 
                  6*A134nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(y[7] + y[16]) - 
                  A12nunuee*c1*pow(3,0.5)*y[7]*y[17] - 12*A12nuan*c3*pow(3,0.5)*y[7]*y[17] + 
                  4*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[7]*y[17] - 8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[7]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[7]*y[17] + 48*A123nuan*c3*qubicTerms*y[7]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[9]*y[17] + 
                  2*A123nunuee*c1*(-3*(-1 + 8*pow(sw,4))*(y[5]*y[10] - y[4]*y[11] - y[2]*y[13] + y[1]*y[14]) - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8])*y[16] - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*y[7]*(6*y[9] - 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[16]*(A12nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + 3*y[0]*y[3] + pow(3,0.5)*y[0]*y[8] + 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*y[3]*y[12] + 
                  4*(pow(y[7],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15]) + 
                  4*y[8]*y[17])))))
                  -
                  znu_lhs*y[7];

      f[8]  =     pow(3,0.5)*(cos(theta13)*(cos(theta23)*(d31 - d21*pow(sin(theta12),2))*sin(theta13) - 
                  d21*cos(theta12)*sin(theta12)*sin(theta23))*y[5] + 
                  d31*cos(theta23)*pow(cos(theta13),2)*sin(theta23)*y[7] - 
                  d21*(cos(theta12)*cos(2*theta23)*sin(theta12)*sin(theta13) + 
                  cos(theta23)*(pow(cos(theta12),2) - pow(sin(theta12),2)*pow(sin(theta13),2))*sin(theta23))*y[7] + 
                  vs*(-(y[5]*y[13]) + y[4]*y[14] - y[7]*y[15] + y[6]*y[16]))
                  +
                  -0.010416666666666666*(c*pow(3,-0.5)*pow(_PI_,-3)*
                  (-24*A34nunuee*c1*pow(sw,2) - 72*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[9],2)*y[8] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[17],2)*y[8] + 72*A12nuan*c3*pow(3,0.5)*y[8]*y[9] + 
                  72*A12nuan*c3*y[1]*y[10] - 72*A123nuan*c3*qubicTerms*y[0]*y[1]*y[10] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[8]*y[10] - 72*A123nuan*c3*qubicTerms*y[1]*y[9]*y[10] + 
                  72*A12nuan*c3*y[2]*y[11] - 72*A123nuan*c3*qubicTerms*y[0]*y[2]*y[11] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[8]*y[11] - 72*A123nuan*c3*qubicTerms*y[2]*y[9]*y[11] + 
                  72*A12nuan*c3*y[3]*y[12] - 72*A123nuan*c3*qubicTerms*y[0]*y[3]*y[12] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[8]*y[12] - 72*A123nuan*c3*qubicTerms*y[3]*y[9]*y[12] - 
                  36*A12nuan*c3*y[4]*y[13] + 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[13] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[8]*y[13] + 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[13] - 
                  36*A12nuan*c3*y[5]*y[14] + 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[14] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[8]*y[14] + 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[14] - 
                  36*A12nuan*c3*y[6]*y[15] + 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[15] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[8]*y[15] + 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[15] - 
                  36*A12nuan*c3*y[7]*y[16] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[16] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[8]*y[16] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[16] - 
                  72*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[0],2)*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[8],2)*y[17] + 72*A12nuan*c3*pow(3,0.5)*y[0]*y[17] - 
                  72*A12nuan*c3*y[8]*y[17] + 72*A123nuan*c3*qubicTerms*y[0]*y[8]*y[17] + 
                  72*A123nuan*c3*qubicTerms*y[8]*y[9]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[10]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[11]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[12]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[13]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[14]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[15]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[16]*y[17] + 
                  2*A134nunuee*c1*(pow(3,0.5)*(3 - 8*pow(sw,2) + 24*pow(sw,4))*y[8] + 
                  12*pow(sw,2)*(y[0] + y[3] + y[9] + y[12]) + pow(3,0.5)*(3 - 8*pow(sw,2) + 24*pow(sw,4))*y[17])\
                  - 2*A123nunuee*c1*(-6*y[1]*y[10] - 6*y[2]*y[11] + 6*y[3]*y[12] + 3*y[4]*y[13] + 3*y[5]*y[14] - 
                  3*y[6]*y[15] - 3*y[7]*y[16] + 6*pow(3,0.5)*y[0]*y[17] + 
                  2*y[8]*(pow(3,0.5)*(3*y[9] + 24*pow(sw,4)*y[9] + 4*pow(sw,2)*(-2*y[9] + y[12])) + 
                  (-3 + 16*pow(sw,2) - 24*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(6*y[0]*y[9] + 6*y[3]*y[9] + 6*y[0]*y[12] + 3*y[6]*y[15] + 3*y[7]*y[16] - 
                  6*pow(sw,2)*(-2*y[1]*y[10] - 2*y[2]*y[11] - 2*y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16]) + 4*pow(3,0.5)*(-1 + 3*pow(sw,2))*y[0]*y[17] + 
                  2*pow(3,0.5)*y[3]*y[17])) + A12nunuee*c1*
                  (-6*y[1]*y[10] - 6*y[2]*y[11] + 6*y[3]*y[12] + 3*y[4]*y[13] + 3*y[5]*y[14] - 3*y[6]*y[15] - 
                  3*y[7]*y[16] + 6*pow(3,0.5)*y[0]*y[17] + 
                  2*y[8]*(pow(3,0.5)*(3*y[9] + 24*pow(sw,4)*y[9] + 4*pow(sw,2)*(-2*y[9] + y[12])) + 
                  (-3 + 16*pow(sw,2) - 24*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(6*y[0]*y[9] + 6*y[3]*y[9] + 6*y[0]*y[12] + 3*y[6]*y[15] + 3*y[7]*y[16] - 
                  6*pow(sw,2)*(-2*y[1]*y[10] - 2*y[2]*y[11] - 2*y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16]) + 4*pow(3,0.5)*(-1 + 3*pow(sw,2))*y[0]*y[17] + 
                  2*pow(3,0.5)*y[3]*y[17]))))
                  -
                  znu_lhs*y[8];


      /* --------- rbar ---------- */

      f[9]  =     (c*c1*pow(_PI_,-3)*(3*A34nunuee*(3 - 4*pow(sw,2) + 24*pow(sw,4)) - 
                  A134nunuee*(9*y[0] + 9*y[9] + 72*pow(sw,4)*(y[0] + y[9]) + 
                  4*pow(sw,2)*(-3*y[0] + 6*y[3] + 2*pow(3,0.5)*y[8] - 3*y[9] + 6*y[12] + 2*pow(3,0.5)*y[17])) + 
                  (2*A123nunuee - A12nunuee)*(9*y[0]*y[9] - 
                  6*(y[1]*y[10] + y[2]*y[11] - y[3]*y[12] + y[4]*y[13] + y[5]*y[14] - y[6]*y[15] - y[7]*y[16] - 
                  y[8]*y[17]) + 4*pow(sw,2)*(y[0]*(-3*y[9] + 6*y[12] + 2*pow(3,0.5)*y[17]) + 
                  2*(3*y[3]*y[9] + pow(3,0.5)*y[8]*y[9] + pow(3,0.5)*y[8]*y[12] - 3*y[6]*y[15] - 
                  3*y[7]*y[16] + pow(3,0.5)*y[3]*y[17] - 2*y[8]*y[17])) + 
                  24*pow(sw,4)*(3*y[0]*y[9] + 2*(y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16] + y[8]*y[17])))))/144.
                  -
                  znu_lhs*y[9];

      f[10] =     -(e1*y[11]) + e2*y[11] + (4*d21*sin(2*theta12)*sin(theta13)*
                  (sin(2*theta23)*y[11] + cos(2*theta23)*y[14]) - 
                  2*(d21 - 2*d31 + 3*d21*cos(2*theta12))*cos(theta23)*(cos(theta23)*y[11] - sin(theta23)*y[14]) + 
                  (-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*
                  ((-3 + cos(2*theta23))*y[11] - sin(2*theta23)*y[14]) + 
                  2*(d21 - 2*d31 - d21*cos(2*theta12))*cos(theta23)*sin(2*theta13)*y[16] + 
                  4*d21*cos(theta13)*sin(2*theta12)*sin(theta23)*y[16])/8. + 
                  vs*(-2*y[3]*y[11] + 2*y[2]*y[12] - y[7]*y[13] + y[6]*y[14] - y[5]*y[15] + y[4]*y[16])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[1] + 
                  48*A123nuan*c3*qubicTerms*pow(y[10],2)*y[1] - 6*A12nunuee*c1*y[1]*y[9] + 72*A12nuan*c3*y[1]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[1]*y[9] + 
                  6*A134nunuee*c1*(-y[1] + y[10] + 8*pow(sw,4)*(y[1] + y[10])) + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[2]*y[11] + 48*A123nuan*c3*qubicTerms*y[1]*y[3]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[4]*y[13] + 3*A12nunuee*c1*y[6]*y[13] + 36*A12nuan*c3*y[6]*y[13] - 
                  12*A12nunuee*c1*pow(sw,2)*y[6]*y[13] + 24*A12nunuee*c1*pow(sw,4)*y[6]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[13] - 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[13] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[5]*y[14] + 3*A12nunuee*c1*y[7]*y[14] + 36*A12nuan*c3*y[7]*y[14] - 
                  12*A12nunuee*c1*pow(sw,2)*y[7]*y[14] + 24*A12nunuee*c1*pow(sw,4)*y[7]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[14] - 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[14] - 
                  3*A12nunuee*c1*y[4]*y[15] + 36*A12nuan*c3*y[4]*y[15] + 24*A12nunuee*c1*pow(sw,4)*y[4]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[15] + 48*A123nuan*c3*qubicTerms*y[1]*y[6]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[15] - 3*A12nunuee*c1*y[5]*y[16] + 36*A12nuan*c3*y[5]*y[16] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[16] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[16] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[7]*y[16] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[16] - 
                  2*A12nunuee*c1*pow(3,0.5)*y[1]*y[17] + 24*A12nuan*c3*pow(3,0.5)*y[1]*y[17] + 
                  16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[1]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[1]*y[17] + 48*A123nuan*c3*qubicTerms*y[1]*y[8]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[9]*y[17] - 
                  2*A123nunuee*c1*(2*pow(3,0.5)*(1 + 8*pow(sw,4))*y[8]*y[10] + 6*y[0]*(y[10] + 8*pow(sw,4)*y[10]) + 
                  3*(y[6]*y[13] + y[7]*y[14] - 4*pow(sw,2)*(-2*y[3]*y[10] + y[6]*y[13] + y[7]*y[14]) - 
                  y[4]*y[15] - y[5]*y[16] + 8*pow(sw,4)*(y[6]*y[13] + y[7]*y[14] + y[4]*y[15] + y[5]*y[16])) + 
                  2*(-1 + 8*pow(sw,4))*y[1]*(3*y[9] + pow(3,0.5)*y[17])) + 
                  2*y[10]*(-6*A124nuenue*c2 + 6*A12nuenue*c2 + 
                  A12nunuee*c1*(3*(1 + 8*pow(sw,4))*y[0] + 12*pow(sw,2)*y[3] + 
                  pow(3,0.5)*(1 + 8*pow(sw,4))*y[8]) + 
                  12*c3*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[1],2) + y[2]*y[11] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17])))))
                  -
                  znu_lhs*y[10];

      f[11] =     (8*e1*y[10] - 8*e2*y[10] - 4*d21*sin(2*theta12)*sin(theta13)*
                  (sin(2*theta23)*y[10] + cos(2*theta23)*y[13]) + 
                  2*(d21 - 2*d31 + 3*d21*cos(2*theta12))*cos(theta23)*(cos(theta23)*y[10] - sin(theta23)*y[13]) - 
                  (-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*
                  ((-3 + cos(2*theta23))*y[10] - sin(2*theta23)*y[13]) - 
                  2*(-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(-2*sin(theta23)*y[12] + cos(theta23)*y[15]) + 
                  4*d21*cos(theta13)*sin(2*theta12)*(2*cos(theta23)*y[12] + sin(theta23)*y[15]) + 
                  8*vs*(2*y[3]*y[10] - 2*y[1]*y[12] - y[6]*y[13] - y[7]*y[14] + y[4]*y[15] + y[5]*y[16]))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[2] + 
                  48*A123nuan*c3*qubicTerms*pow(y[11],2)*y[2] - 6*A12nunuee*c1*y[2]*y[9] + 72*A12nuan*c3*y[2]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[2]*y[9] + 48*A123nuan*c3*qubicTerms*y[1]*y[2]*y[10] + 
                  6*A134nunuee*c1*(-y[2] + y[11] + 8*pow(sw,4)*(y[2] + y[11])) + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[3]*y[12] + 48*A123nuan*c3*qubicTerms*y[2]*y[4]*y[13] - 
                  3*A12nunuee*c1*y[7]*y[13] - 36*A12nuan*c3*y[7]*y[13] + 12*A12nunuee*c1*pow(sw,2)*y[7]*y[13] - 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[13] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[13] + 48*A123nuan*c3*qubicTerms*y[2]*y[5]*y[14] + 
                  3*A12nunuee*c1*y[6]*y[14] + 36*A12nuan*c3*y[6]*y[14] - 12*A12nunuee*c1*pow(sw,2)*y[6]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[14] - 3*A12nunuee*c1*y[5]*y[15] + 36*A12nuan*c3*y[5]*y[15] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[15] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[15] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[6]*y[15] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[15] + 
                  3*A12nunuee*c1*y[4]*y[16] - 36*A12nuan*c3*y[4]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[4]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[16] + 48*A123nuan*c3*qubicTerms*y[2]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[16] - 2*A12nunuee*c1*pow(3,0.5)*y[2]*y[17] + 
                  24*A12nuan*c3*pow(3,0.5)*y[2]*y[17] + 16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[2]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[2]*y[17] + 48*A123nuan*c3*qubicTerms*y[2]*y[8]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[9]*y[17] - 
                  2*A123nunuee*c1*(2*pow(3,0.5)*(1 + 8*pow(sw,4))*y[8]*y[11] + 6*y[0]*(y[11] + 8*pow(sw,4)*y[11]) + 
                  3*(-(y[7]*y[13]) + y[6]*y[14] + 4*pow(sw,2)*(2*y[3]*y[11] + y[7]*y[13] - y[6]*y[14]) - 
                  y[5]*y[15] + y[4]*y[16] + 8*pow(sw,4)*(-(y[7]*y[13]) + y[6]*y[14] + y[5]*y[15] - y[4]*y[16]))
                  + 2*(-1 + 8*pow(sw,4))*y[2]*(3*y[9] + pow(3,0.5)*y[17])) + 
                  2*y[11]*(-6*A124nuenue*c2 + 6*A12nuenue*c2 + 
                  A12nunuee*c1*(3*(1 + 8*pow(sw,4))*y[0] + 12*pow(sw,2)*y[3] + 
                  pow(3,0.5)*(1 + 8*pow(sw,4))*y[8]) + 
                  12*c3*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[2],2) + y[1]*y[10] + y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17])))))
                  -
                  znu_lhs*y[11];

      f[12] =     cos(theta13)*sin(theta13)*(-2*d31*sin(theta23)*y[11] + 
                  cos(theta23)*(-d31 + d21*pow(sin(theta12),2))*y[14]) - 
                  d21*cos(theta23)*pow(cos(theta12),2)*sin(theta23)*y[16] + 
                  d31*cos(theta23)*pow(cos(theta13),2)*sin(theta23)*y[16] + 
                  d21*pow(sin(theta12),2)*sin(theta23)*(sin(2*theta13)*y[11] + cos(theta23)*pow(sin(theta13),2)*y[16]) + 
                  d21*cos(theta12)*sin(theta12)*(cos(theta13)*(-2*cos(theta23)*y[11] + sin(theta23)*y[14]) - 
                  cos(2*theta23)*sin(theta13)*y[16]) + 
                  vs*(-2*y[2]*y[10] + 2*y[1]*y[11] - y[5]*y[13] + y[4]*y[14] + y[7]*y[15] - y[6]*y[16])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-24*A34nunuee*c1*pow(sw,2) - 
                  72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[3] + 48*A123nuan*c3*qubicTerms*pow(y[12],2)*y[3] + 
                  24*A12nunuee*c1*pow(sw,2)*y[0]*y[9] + 6*A12nunuee*c1*y[3]*y[9] + 72*A12nuan*c3*y[3]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[3]*y[9] + 8*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[8]*y[9] + 
                  48*A123nuan*c3*qubicTerms*y[1]*y[3]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[3]*y[11] + 
                  2*A12nunuee*c1*(3*(1 + 8*pow(sw,4))*y[0] + 12*pow(sw,2)*y[3] + pow(3,0.5)*(1 + 8*pow(sw,4))*y[8])*
                  y[12] - 3*A12nunuee*c1*y[4]*y[13] + 36*A12nuan*c3*y[4]*y[13] + 
                  24*A12nunuee*c1*pow(sw,4)*y[4]*y[13] - 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[13] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[4]*y[13] - 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[13] - 
                  3*A12nunuee*c1*y[5]*y[14] + 36*A12nuan*c3*y[5]*y[14] + 24*A12nunuee*c1*pow(sw,4)*y[5]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[14] + 48*A123nuan*c3*qubicTerms*y[3]*y[5]*y[14] - 
                  36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[14] - 3*A12nunuee*c1*y[6]*y[15] - 36*A12nuan*c3*y[6]*y[15] + 
                  12*A12nunuee*c1*pow(sw,2)*y[6]*y[15] - 24*A12nunuee*c1*pow(sw,4)*y[6]*y[15] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[15] + 48*A123nuan*c3*qubicTerms*y[3]*y[6]*y[15] + 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[15] - 3*A12nunuee*c1*y[7]*y[16] - 36*A12nuan*c3*y[7]*y[16] + 
                  12*A12nunuee*c1*pow(sw,2)*y[7]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[16] + 48*A123nuan*c3*qubicTerms*y[3]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[16] + 8*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[0]*y[17] + 
                  2*A12nunuee*c1*pow(3,0.5)*y[3]*y[17] + 24*A12nuan*c3*pow(3,0.5)*y[3]*y[17] + 
                  16*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[3]*y[17] - 
                  24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[3]*y[17] + 8*A12nunuee*c1*pow(sw,2)*y[8]*y[17] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[8]*y[17] - 24*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[9]*y[17] + 
                  2*A134nunuee*c1*(3*(y[3] + y[12]) + 24*pow(sw,4)*(y[3] + y[12]) + 
                  4*pow(sw,2)*(3*y[0] + pow(3,0.5)*y[8] + 3*y[9] + pow(3,0.5)*y[17])) + 
                  24*c3*y[12]*(A12nuan*(3*y[0] + pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-3*pow(y[0],2) - pow(3,0.5)*y[0]*y[8] - pow(3,0.5)*y[8]*y[9] + 
                  2*(pow(y[3],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + 
                  y[7]*y[16]) + 2*y[8]*y[17])) - 
                  2*A123nunuee*c1*(2*(1 + 8*pow(sw,4))*(3*y[0] + pow(3,0.5)*y[8])*y[12] - 
                  3*(y[4]*y[13] + y[5]*y[14] + y[6]*y[15] + y[7]*y[16]) + 
                  2*y[3]*(3*(1 + 8*pow(sw,4))*y[9] + 12*pow(sw,2)*y[12] + pow(3,0.5)*(1 + 8*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(3*y[6]*y[15] + 3*y[7]*y[16] + 
                  6*pow(sw,2)*(y[4]*y[13] + y[5]*y[14] - y[6]*y[15] - y[7]*y[16]) + 
                  2*y[8]*(pow(3,0.5)*y[9] + y[17]) + 2*y[0]*(3*y[9] + pow(3,0.5)*y[17])))))
                  -
                  znu_lhs*y[12];

      f[13] =     -0.125*(d21*cos(2*theta12)*((-3 + cos(2*theta13))*cos(2*theta23) + 6*pow(cos(theta13),2))*y[14]) - 
                  ((d21 - 2*d31 + 8*e1 - (d21 - 2*d31)*(3*cos(2*theta13) + 2*cos(2*theta23)*pow(cos(theta13),2)) + 
                  4*d21*sin(2*theta12)*sin(theta13)*sin(2*theta23) + 8*vs*(y[3] + pow(3,0.5)*y[8]))*y[14])/8. + 
                  d21*cos(theta12)*sin(theta12)*(cos(2*theta23)*sin(theta13)*y[11] + cos(theta13)*cos(theta23)*y[16]) + 
                  sin(theta23)*(cos(theta23)*(d21*pow(cos(theta12),2) - d31*pow(cos(theta13),2))*y[11] - 
                  d21*cos(theta23)*pow(sin(theta12),2)*pow(sin(theta13),2)*y[11] + 
                  cos(theta13)*(d31 - d21*pow(sin(theta12),2))*sin(theta13)*y[16]) + 
                  vs*(y[7]*y[10] + y[6]*y[11] + y[5]*y[12] - y[2]*y[15] - y[1]*y[16] + pow(3,0.5)*y[5]*y[17])
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[4] + 
                  48*A123nuan*c3*qubicTerms*pow(y[13],2)*y[4] - 6*A12nunuee*c1*y[4]*y[9] + 72*A12nuan*c3*y[4]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[4]*y[9] + 48*A123nuan*c3*qubicTerms*y[1]*y[4]*y[10] + 
                  3*A12nunuee*c1*y[6]*y[10] + 36*A12nuan*c3*y[6]*y[10] - 12*A12nunuee*c1*pow(sw,2)*y[6]*y[10] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[10] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[4]*y[11] - 
                  3*A12nunuee*c1*y[7]*y[11] - 36*A12nuan*c3*y[7]*y[11] + 12*A12nunuee*c1*pow(sw,2)*y[7]*y[11] - 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[11] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[11] + 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[11] - 3*A12nunuee*c1*y[4]*y[12] + 36*A12nuan*c3*y[4]*y[12] + 
                  24*A12nunuee*c1*pow(sw,4)*y[4]*y[12] - 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[4]*y[12] - 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[12] + 
                  6*A134nunuee*c1*(-y[4] + y[13] + 8*pow(sw,4)*(y[4] + y[13])) + 
                  48*A123nuan*c3*qubicTerms*y[4]*y[5]*y[14] - 3*A12nunuee*c1*y[1]*y[15] + 36*A12nuan*c3*y[1]*y[15] + 
                  24*A12nunuee*c1*pow(sw,4)*y[1]*y[15] - 36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[15] + 
                  48*A123nuan*c3*qubicTerms*y[4]*y[6]*y[15] - 36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[15] + 
                  3*A12nunuee*c1*y[2]*y[16] - 36*A12nuan*c3*y[2]*y[16] - 24*A12nunuee*c1*pow(sw,4)*y[2]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[16] + 48*A123nuan*c3*qubicTerms*y[4]*y[7]*y[16] + 
                  36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[16] + A12nunuee*c1*pow(3,0.5)*y[4]*y[17] - 
                  12*A12nuan*c3*pow(3,0.5)*y[4]*y[17] - 8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[4]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[4]*y[17] + 48*A123nuan*c3*qubicTerms*y[4]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[9]*y[17] - 
                  2*A123nunuee*c1*(3*y[6]*y[10] - 3*y[7]*y[11] + 6*y[0]*y[13] + 3*y[3]*y[13] - 
                  pow(3,0.5)*y[8]*y[13] + 12*pow(sw,2)*
                  (-(y[6]*y[10]) + y[7]*y[11] + y[3]*y[13] + pow(3,0.5)*y[8]*y[13]) - 3*y[1]*y[15] + 
                  3*y[2]*y[16] + 8*pow(sw,4)*(3*y[6]*y[10] - 3*y[7]*y[11] + 6*y[0]*y[13] + 3*y[3]*y[13] - 
                  pow(3,0.5)*y[8]*y[13] + 3*y[1]*y[15] - 3*y[2]*y[16]) + 
                  (-1 + 8*pow(sw,4))*y[4]*(6*y[9] + 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[13]*(12*(-A124nuenue + A12nuenue)*c2 + 
                  A12nunuee*c1*(6*(1 + 8*pow(sw,4))*y[0] + 3*(1 + 4*pow(sw,2) + 8*pow(sw,4))*y[3] + 
                  pow(3,0.5)*(-1 + 12*pow(sw,2) - 8*pow(sw,4))*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + y[0]*(-3*y[3] + pow(3,0.5)*y[8]) - 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*(pow(y[4],2) + y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + 
                  y[5]*y[14] + y[6]*y[15] + y[7]*y[16] + y[8]*y[17]))))))
                  -
                  znu_lhs*y[13];

      f[14] =     ((d21 - 2*d31 + 8*e1 + 8*vs*(y[3] + pow(3,0.5)*y[8]))*y[13] - 
                  (d21 - 2*d31 + 3*d21*cos(2*theta12))*(sin(2*theta23)*y[10] + cos(2*theta23)*y[13]) - 
                  (d21 - 2*d31)*cos(2*theta13)*(sin(2*theta23)*y[10] + (3 + cos(2*theta23))*y[13]) - 
                  8*vs*(y[6]*y[10] - y[7]*y[11] + y[4]*y[12] - y[1]*y[15] + y[2]*y[16] + pow(3,0.5)*y[4]*y[17]) - 
                  2*(d21 - 2*d31)*sin(2*theta13)*(-(sin(theta23)*y[15]) + cos(theta23)*(y[12] + pow(3,0.5)*y[17])) + 
                  d21*cos(2*theta12)*(6*pow(cos(theta13),2)*y[13] + 
                  cos(2*theta13)*(sin(2*theta23)*y[10] + cos(2*theta23)*y[13]) + 
                  2*sin(2*theta13)*(-(sin(theta23)*y[15]) + cos(theta23)*(y[12] + pow(3,0.5)*y[17]))) - 
                  4*d21*sin(2*theta12)*(sin(theta13)*(cos(2*theta23)*y[10] - sin(2*theta23)*y[13]) + 
                  cos(theta13)*(cos(theta23)*y[15] + sin(theta23)*(y[12] + pow(3,0.5)*y[17]))))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[5] + 
                  48*A123nuan*c3*qubicTerms*pow(y[14],2)*y[5] - 6*A12nunuee*c1*y[5]*y[9] + 72*A12nuan*c3*y[5]*y[9] + 
                  48*A12nunuee*c1*pow(sw,4)*y[5]*y[9] + 48*A123nuan*c3*qubicTerms*y[1]*y[5]*y[10] + 
                  3*A12nunuee*c1*y[7]*y[10] + 36*A12nuan*c3*y[7]*y[10] - 12*A12nunuee*c1*pow(sw,2)*y[7]*y[10] + 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[10] - 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[10] + 48*A123nuan*c3*qubicTerms*y[2]*y[5]*y[11] + 
                  3*A12nunuee*c1*y[6]*y[11] + 36*A12nuan*c3*y[6]*y[11] - 12*A12nunuee*c1*pow(sw,2)*y[6]*y[11] + 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[11] - 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[11] - 
                  36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[11] - 3*A12nunuee*c1*y[5]*y[12] + 36*A12nuan*c3*y[5]*y[12] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[12] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[5]*y[12] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[4]*y[5]*y[13] + 
                  6*A134nunuee*c1*(-y[5] + y[14] + 8*pow(sw,4)*(y[5] + y[14])) - 3*A12nunuee*c1*y[2]*y[15] + 
                  36*A12nuan*c3*y[2]*y[15] + 24*A12nunuee*c1*pow(sw,4)*y[2]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[15] + 48*A123nuan*c3*qubicTerms*y[5]*y[6]*y[15] - 
                  36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[15] - 3*A12nunuee*c1*y[1]*y[16] + 36*A12nuan*c3*y[1]*y[16] + 
                  24*A12nunuee*c1*pow(sw,4)*y[1]*y[16] - 36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[16] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[7]*y[16] - 36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[16] + 
                  A12nunuee*c1*pow(3,0.5)*y[5]*y[17] - 12*A12nuan*c3*pow(3,0.5)*y[5]*y[17] - 
                  8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[5]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[5]*y[17] + 48*A123nuan*c3*qubicTerms*y[5]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[9]*y[17] - 
                  2*A123nunuee*c1*((6*(1 + 8*pow(sw,4))*y[0] + 3*(1 + 4*pow(sw,2) + 8*pow(sw,4))*y[3] + 
                  pow(3,0.5)*(-1 + 12*pow(sw,2) - 8*pow(sw,4))*y[8])*y[14] + 
                  3*(y[7]*y[10] + y[6]*y[11] - 4*pow(sw,2)*(y[7]*y[10] + y[6]*y[11]) - y[2]*y[15] - y[1]*y[16] + 
                  8*pow(sw,4)*(y[7]*y[10] + y[6]*y[11] + y[2]*y[15] + y[1]*y[16])) + 
                  (-1 + 8*pow(sw,4))*y[5]*(6*y[9] + 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[14]*(12*(-A124nuenue + A12nuenue)*c2 + 
                  A12nunuee*c1*(6*(1 + 8*pow(sw,4))*y[0] + 3*(1 + 4*pow(sw,2) + 8*pow(sw,4))*y[3] + 
                  pow(3,0.5)*(-1 + 12*pow(sw,2) - 8*pow(sw,4))*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] + 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + y[0]*(-3*y[3] + pow(3,0.5)*y[8]) - 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*(pow(y[5],2) + y[1]*y[10] + y[2]*y[11] + y[3]*y[12] + 
                  y[4]*y[13] + y[6]*y[15] + y[7]*y[16] + y[8]*y[17]))))))
                  -
                  znu_lhs*y[14];

      f[15] =     ((-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(cos(theta23)*y[11] + sin(theta23)*y[14]) - 
                  4*e2*y[16] + cos(2*theta23)*(-(d21*cos(2*theta12)*(-3 + cos(2*theta13))) + 
                  2*(d21 - 2*d31)*pow(cos(theta13),2))*y[16] + 
                  4*d21*cos(theta12)*sin(theta12)*(cos(theta13)*cos(theta23)*y[14] - 
                  sin(theta23)*(cos(theta13)*y[11] + 4*cos(theta23)*sin(theta13)*y[16])) + 
                  4*vs*(y[5]*y[10] - y[4]*y[11] - y[7]*y[12] + y[2]*y[13] - y[1]*y[14] + y[3]*y[16] - 
                  pow(3,0.5)*y[8]*y[16] + pow(3,0.5)*y[7]*y[17]))/4.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[6] + 
                  48*A123nuan*c3*qubicTerms*pow(y[15],2)*y[6] + 6*A12nunuee*c1*y[6]*y[9] + 72*A12nuan*c3*y[6]*y[9] - 
                  24*A12nunuee*c1*pow(sw,2)*y[6]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[6]*y[9] - 
                  3*A12nunuee*c1*y[4]*y[10] + 36*A12nuan*c3*y[4]*y[10] + 24*A12nunuee*c1*pow(sw,4)*y[4]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[10] + 48*A123nuan*c3*qubicTerms*y[1]*y[6]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[10] - 3*A12nunuee*c1*y[5]*y[11] + 36*A12nuan*c3*y[5]*y[11] + 
                  24*A12nunuee*c1*pow(sw,4)*y[5]*y[11] - 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[11] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[6]*y[11] - 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[11] - 
                  3*A12nunuee*c1*y[6]*y[12] - 36*A12nuan*c3*y[6]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[6]*y[12] - 
                  24*A12nunuee*c1*pow(sw,4)*y[6]*y[12] + 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[6]*y[12] + 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[12] - 
                  3*A12nunuee*c1*y[1]*y[13] + 36*A12nuan*c3*y[1]*y[13] + 24*A12nunuee*c1*pow(sw,4)*y[1]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[13] + 48*A123nuan*c3*qubicTerms*y[4]*y[6]*y[13] - 
                  36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[13] - 3*A12nunuee*c1*y[2]*y[14] + 36*A12nuan*c3*y[2]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[2]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[6]*y[14] - 36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[14] + 
                  6*A134nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(y[6] + y[15]) + 
                  48*A123nuan*c3*qubicTerms*y[6]*y[7]*y[16] - A12nunuee*c1*pow(3,0.5)*y[6]*y[17] - 
                  12*A12nuan*c3*pow(3,0.5)*y[6]*y[17] + 4*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[6]*y[17] - 
                  8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[6]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[6]*y[17] + 48*A123nuan*c3*qubicTerms*y[6]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[9]*y[17] + 
                  2*A123nunuee*c1*(-3*(-1 + 8*pow(sw,4))*(y[4]*y[10] + y[5]*y[11] + y[1]*y[13] + y[2]*y[14]) - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8])*y[15] - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*y[6]*(6*y[9] - 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[15]*(A12nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + 3*y[0]*y[3] + pow(3,0.5)*y[0]*y[8] + 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*y[3]*y[12] + 
                  4*(pow(y[6],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[7]*y[16]) + 
                  4*y[8]*y[17])))))
                  -
                  znu_lhs*y[15];

      f[16] =     (-4*d21*cos(theta13)*sin(2*theta12)*(sin(theta23)*y[10] + cos(theta23)*y[13]) + 
                  2*(-d21 + 2*d31 + d21*cos(2*theta12))*sin(2*theta13)*(cos(theta23)*y[10] - sin(theta23)*y[13]) + 
                  8*e2*y[15] - 8*vs*(y[4]*y[10] + y[5]*y[11] - y[6]*y[12] - y[1]*y[13] - y[2]*y[14] + y[3]*y[15] - 
                  pow(3,0.5)*y[8]*y[15] + pow(3,0.5)*y[6]*y[17]) + 
                  4*d21*sin(2*theta12)*sin(theta13)*(2*sin(2*theta23)*y[15] + 
                  cos(2*theta23)*(y[12] - pow(3,0.5)*y[17])) - 
                  (d21 - 2*d31 + 3*d21*cos(2*theta12))*(2*cos(2*theta23)*y[15] + 
                  sin(2*theta23)*(-y[12] + pow(3,0.5)*y[17])) + 
                  (-d21 + 2*d31 + d21*cos(2*theta12))*cos(2*theta13)*
                  (2*cos(2*theta23)*y[15] + sin(2*theta23)*(-y[12] + pow(3,0.5)*y[17])))/8.
                  +
                  -0.010416666666666666*(c*pow(_PI_,-3)*(-72*A123nuan*c3*qubicTerms*pow(y[9],2)*y[7] + 
                  48*A123nuan*c3*qubicTerms*pow(y[16],2)*y[7] + 6*A12nunuee*c1*y[7]*y[9] + 72*A12nuan*c3*y[7]*y[9] - 
                  24*A12nunuee*c1*pow(sw,2)*y[7]*y[9] + 48*A12nunuee*c1*pow(sw,4)*y[7]*y[9] - 
                  3*A12nunuee*c1*y[5]*y[10] + 36*A12nuan*c3*y[5]*y[10] + 24*A12nunuee*c1*pow(sw,4)*y[5]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[10] + 48*A123nuan*c3*qubicTerms*y[1]*y[7]*y[10] - 
                  36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[10] + 3*A12nunuee*c1*y[4]*y[11] - 36*A12nuan*c3*y[4]*y[11] - 
                  24*A12nunuee*c1*pow(sw,4)*y[4]*y[11] + 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[11] + 
                  48*A123nuan*c3*qubicTerms*y[2]*y[7]*y[11] + 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[11] - 
                  3*A12nunuee*c1*y[7]*y[12] - 36*A12nuan*c3*y[7]*y[12] + 12*A12nunuee*c1*pow(sw,2)*y[7]*y[12] - 
                  24*A12nunuee*c1*pow(sw,4)*y[7]*y[12] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[12] + 
                  48*A123nuan*c3*qubicTerms*y[3]*y[7]*y[12] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[12] + 
                  3*A12nunuee*c1*y[2]*y[13] - 36*A12nuan*c3*y[2]*y[13] - 24*A12nunuee*c1*pow(sw,4)*y[2]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[0]*y[2]*y[13] + 48*A123nuan*c3*qubicTerms*y[4]*y[7]*y[13] + 
                  36*A123nuan*c3*qubicTerms*y[2]*y[9]*y[13] - 3*A12nunuee*c1*y[1]*y[14] + 36*A12nuan*c3*y[1]*y[14] + 
                  24*A12nunuee*c1*pow(sw,4)*y[1]*y[14] - 36*A123nuan*c3*qubicTerms*y[0]*y[1]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[5]*y[7]*y[14] - 36*A123nuan*c3*qubicTerms*y[1]*y[9]*y[14] + 
                  48*A123nuan*c3*qubicTerms*y[6]*y[7]*y[15] + 
                  6*A134nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(y[7] + y[16]) - 
                  A12nunuee*c1*pow(3,0.5)*y[7]*y[17] - 12*A12nuan*c3*pow(3,0.5)*y[7]*y[17] + 
                  4*A12nunuee*c1*pow(3,0.5)*pow(sw,2)*y[7]*y[17] - 8*A12nunuee*c1*pow(3,0.5)*pow(sw,4)*y[7]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[0]*y[7]*y[17] + 48*A123nuan*c3*qubicTerms*y[7]*y[8]*y[17] + 
                  12*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[9]*y[17] + 
                  2*A123nunuee*c1*(-3*(-1 + 8*pow(sw,4))*(y[5]*y[10] - y[4]*y[11] - y[2]*y[13] + y[1]*y[14]) - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8])*y[16] - 
                  (1 - 4*pow(sw,2) + 8*pow(sw,4))*y[7]*(6*y[9] - 3*y[12] - pow(3,0.5)*y[17])) + 
                  y[16]*(A12nunuee*c1*(1 - 4*pow(sw,2) + 8*pow(sw,4))*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  12*c3*(A12nuan*(6*y[0] - 3*y[3] - pow(3,0.5)*y[8]) + 
                  A123nuan*qubicTerms*(-6*pow(y[0],2) + 3*y[0]*y[3] + pow(3,0.5)*y[0]*y[8] + 3*y[3]*y[9] + 
                  pow(3,0.5)*y[8]*y[9] + 4*y[3]*y[12] + 
                  4*(pow(y[7],2) + y[1]*y[10] + y[2]*y[11] + y[4]*y[13] + y[5]*y[14] + y[6]*y[15]) + 
                  4*y[8]*y[17])))))
                  -
                  znu_lhs*y[16];

      f[17] =     pow(3,0.5)*(cos(theta13)*(cos(theta23)*(-d31 + d21*pow(sin(theta12),2))*sin(theta13) + 
                  d21*cos(theta12)*sin(theta12)*sin(theta23))*y[14] - 
                  d31*cos(theta23)*pow(cos(theta13),2)*sin(theta23)*y[16] + 
                  d21*(cos(theta12)*cos(2*theta23)*sin(theta12)*sin(theta13) + 
                  cos(theta23)*(pow(cos(theta12),2) - pow(sin(theta12),2)*pow(sin(theta13),2))*sin(theta23))*y[16]\
                  + vs*(-(y[5]*y[13]) + y[4]*y[14] - y[7]*y[15] + y[6]*y[16]))
                  +
                  -0.010416666666666666*(c*pow(3,-0.5)*pow(_PI_,-3)*
                  (-24*A34nunuee*c1*pow(sw,2) - 72*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[9],2)*y[8] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[17],2)*y[8] + 72*A12nuan*c3*pow(3,0.5)*y[8]*y[9] + 
                  72*A12nuan*c3*y[1]*y[10] - 72*A123nuan*c3*qubicTerms*y[0]*y[1]*y[10] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[8]*y[10] - 72*A123nuan*c3*qubicTerms*y[1]*y[9]*y[10] + 
                  72*A12nuan*c3*y[2]*y[11] - 72*A123nuan*c3*qubicTerms*y[0]*y[2]*y[11] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[8]*y[11] - 72*A123nuan*c3*qubicTerms*y[2]*y[9]*y[11] + 
                  72*A12nuan*c3*y[3]*y[12] - 72*A123nuan*c3*qubicTerms*y[0]*y[3]*y[12] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[8]*y[12] - 72*A123nuan*c3*qubicTerms*y[3]*y[9]*y[12] - 
                  36*A12nuan*c3*y[4]*y[13] + 36*A123nuan*c3*qubicTerms*y[0]*y[4]*y[13] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[8]*y[13] + 36*A123nuan*c3*qubicTerms*y[4]*y[9]*y[13] - 
                  36*A12nuan*c3*y[5]*y[14] + 36*A123nuan*c3*qubicTerms*y[0]*y[5]*y[14] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[8]*y[14] + 36*A123nuan*c3*qubicTerms*y[5]*y[9]*y[14] - 
                  36*A12nuan*c3*y[6]*y[15] + 36*A123nuan*c3*qubicTerms*y[0]*y[6]*y[15] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[8]*y[15] + 36*A123nuan*c3*qubicTerms*y[6]*y[9]*y[15] - 
                  36*A12nuan*c3*y[7]*y[16] + 36*A123nuan*c3*qubicTerms*y[0]*y[7]*y[16] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[8]*y[16] + 36*A123nuan*c3*qubicTerms*y[7]*y[9]*y[16] - 
                  72*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[0],2)*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*pow(y[8],2)*y[17] + 72*A12nuan*c3*pow(3,0.5)*y[0]*y[17] - 
                  72*A12nuan*c3*y[8]*y[17] + 72*A123nuan*c3*qubicTerms*y[0]*y[8]*y[17] + 
                  72*A123nuan*c3*qubicTerms*y[8]*y[9]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[1]*y[10]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[2]*y[11]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[3]*y[12]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[4]*y[13]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[5]*y[14]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[6]*y[15]*y[17] + 
                  48*A123nuan*c3*qubicTerms*pow(3,0.5)*y[7]*y[16]*y[17] + 
                  2*A134nunuee*c1*(pow(3,0.5)*(3 - 8*pow(sw,2) + 24*pow(sw,4))*y[8] + 
                  12*pow(sw,2)*(y[0] + y[3] + y[9] + y[12]) + pow(3,0.5)*(3 - 8*pow(sw,2) + 24*pow(sw,4))*y[17])\
                  - 2*A123nunuee*c1*(-6*y[1]*y[10] - 6*y[2]*y[11] + 6*y[3]*y[12] + 3*y[4]*y[13] + 3*y[5]*y[14] - 
                  3*y[6]*y[15] - 3*y[7]*y[16] + 6*pow(3,0.5)*y[0]*y[17] + 
                  2*y[8]*(pow(3,0.5)*(3*y[9] + 24*pow(sw,4)*y[9] + 4*pow(sw,2)*(-2*y[9] + y[12])) + 
                  (-3 + 16*pow(sw,2) - 24*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(6*y[0]*y[9] + 6*y[3]*y[9] + 6*y[0]*y[12] + 3*y[6]*y[15] + 3*y[7]*y[16] - 
                  6*pow(sw,2)*(-2*y[1]*y[10] - 2*y[2]*y[11] - 2*y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16]) + 4*pow(3,0.5)*(-1 + 3*pow(sw,2))*y[0]*y[17] + 
                  2*pow(3,0.5)*y[3]*y[17])) + A12nunuee*c1*
                  (-6*y[1]*y[10] - 6*y[2]*y[11] + 6*y[3]*y[12] + 3*y[4]*y[13] + 3*y[5]*y[14] - 3*y[6]*y[15] - 
                  3*y[7]*y[16] + 6*pow(3,0.5)*y[0]*y[17] + 
                  2*y[8]*(pow(3,0.5)*(3*y[9] + 24*pow(sw,4)*y[9] + 4*pow(sw,2)*(-2*y[9] + y[12])) + 
                  (-3 + 16*pow(sw,2) - 24*pow(sw,4))*y[17]) + 
                  4*pow(sw,2)*(6*y[0]*y[9] + 6*y[3]*y[9] + 6*y[0]*y[12] + 3*y[6]*y[15] + 3*y[7]*y[16] - 
                  6*pow(sw,2)*(-2*y[1]*y[10] - 2*y[2]*y[11] - 2*y[3]*y[12] + y[4]*y[13] + y[5]*y[14] + 
                  y[6]*y[15] + y[7]*y[16]) + 4*pow(3,0.5)*(-1 + 3*pow(sw,2))*y[0]*y[17] + 
                  2*pow(3,0.5)*y[3]*y[17]))))
                  -
                  znu_lhs*y[17];

      



      /* -------------- temperature evolution ---------------- */

      // derivative of the trace; this has to be defined after f0, f9 are set/evaluated
      double Trace_dr_drb;
      Trace_dr_drb  = 3. * (f[0] + f[9]);

      // for zgamma evolution
      double rhoe_pe, drhoe, rhog, drhog, rhoe;
      rhog        = pow(_PI_,2.)/15.*pow(Tg,4.);
      drhog       = pow(_PI_,2.)/15.*4.*pow(Tg,3.);
      rhoe        = density.rhoe(Tg);
      rhoe_pe     = rhoe + density.pe(Tg);
      drhoe       = density.drhoe(Tg);

      // for znu evolution
      double rhonu, Hubble;
      rhonu       = 7./8. * pow(_PI_,2.)/30. * Trace_r_rb * pow(Tnu,4.);
      Hubble      = pow(_Tref_/t,2.)/mpl;


      // transfer rate
      double I_transfer;
      double inorm, i1,i2;
      // i1, i2 dummy variables to isolate nuenue and nunuee part
      i1=1; i2=1;
      inorm= pow(_GF_,2.) * pow(_Tref_/t,9.);


      I_transfer  =  (inorm*pow(_PI_,-5)*(2*A34nunueeI1*i1*(3 - 4*pow(sw,2) + 24*pow(sw,4)) - A12nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) + A34nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) - 
                        A12nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) + A34nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) - A12nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) + 
                        A34nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) + 4*A12nunueeI1*i1*y[1]*y[10] + 4*A12nunueeI1*i1*y[2]*y[11] + 4*A12nunueeI1*i1*y[4]*y[13] + 
                        4*A12nunueeI1*i1*y[5]*y[14] - 4*A12nunueeI1*i1*y[6]*y[15] - 4*A12nunueeI1*i1*y[7]*y[16] - A12nuenueI1*i2*(y[9] - 2*pow(3,-0.5)*y[17]) + 
                        A34nuenueI1*i2*(y[9] - 2*pow(3,-0.5)*y[17]) - 2*A12nunueeI1*i1*(y[0] - 2*pow(3,-0.5)*y[8])*(y[9] - 2*pow(3,-0.5)*y[17]) - 
                        A12nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) + A34nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) - 
                        2*A12nunueeI1*i1*(y[0] - y[3] + pow(3,-0.5)*y[8])*(y[9] - y[12] + pow(3,-0.5)*y[17]) - A12nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) + 
                        A34nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) - 2*A12nunueeI1*i1*(y[0] + y[3] + pow(3,-0.5)*y[8])*(y[9] + y[12] + pow(3,-0.5)*y[17]) + 
                        4*pow(sw,2)*(A12nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) - A34nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) + A12nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) - 
                        A34nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) - A12nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) + A34nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) + 
                        4*A12nunueeI1*i1*y[6]*y[15] + 4*A12nunueeI1*i1*y[7]*y[16] + (A12nuenueI1 - A34nuenueI1)*i2*(y[9] - 2*pow(3,-0.5)*y[17]) + 
                        2*A12nunueeI1*i1*(y[0] - 2*pow(3,-0.5)*y[8])*(y[9] - 2*pow(3,-0.5)*y[17]) + A12nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) - 
                        A34nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) + 2*A12nunueeI1*i1*(y[0] - y[3] + pow(3,-0.5)*y[8])*(y[9] - y[12] + pow(3,-0.5)*y[17]) - 
                        A12nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) + A34nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) - 
                        2*A12nunueeI1*i1*(y[0] + y[3] + pow(3,-0.5)*y[8])*(y[9] + y[12] + pow(3,-0.5)*y[17])) - 
                        8*pow(sw,4)*(A12nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) - A34nuenueI1*i2*(y[0] - 2*pow(3,-0.5)*y[8]) + A12nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) - 
                        A34nuenueI1*i2*(y[0] - y[3] + pow(3,-0.5)*y[8]) + A12nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) - A34nuenueI1*i2*(y[0] + y[3] + pow(3,-0.5)*y[8]) + 
                        4*A12nunueeI1*i1*y[1]*y[10] + 4*A12nunueeI1*i1*y[2]*y[11] + 4*A12nunueeI1*i1*y[4]*y[13] + 4*A12nunueeI1*i1*y[5]*y[14] + 4*A12nunueeI1*i1*y[6]*y[15] + 
                        4*A12nunueeI1*i1*y[7]*y[16] + (A12nuenueI1 - A34nuenueI1)*i2*(y[9] - 2*pow(3,-0.5)*y[17]) + 
                        2*A12nunueeI1*i1*(y[0] - 2*pow(3,-0.5)*y[8])*(y[9] - 2*pow(3,-0.5)*y[17]) + A12nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) - 
                        A34nuenueI1*i2*(y[9] - y[12] + pow(3,-0.5)*y[17]) + 2*A12nunueeI1*i1*(y[0] - y[3] + pow(3,-0.5)*y[8])*(y[9] - y[12] + pow(3,-0.5)*y[17]) + 
                        A12nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) - A34nuenueI1*i2*(y[9] + y[12] + pow(3,-0.5)*y[17]) + 
                        2*A12nunueeI1*i1*(y[0] + y[3] + pow(3,-0.5)*y[8])*(y[9] + y[12] + pow(3,-0.5)*y[17]))))/32.;   
      

      
      /* ----- znu, zgamma evolution */      
      f[18] = - zon * ( 1./4.*y[18] * ( Trace_dr_drb/Trace_r_rb - I_transfer/(t*Hubble*rhonu)) );
      f[19] = zon * ( y[19]/t - (4*rhog + 3*rhoe_pe) / (_Tref_*(drhoe+drhog)) - I_transfer/(Hubble*_Tref_*(drhoe+drhog)) );


      

      return GSL_SUCCESS;
}



/* Numerical calculation of the Jacobian matrix of ode_fun() */
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
      // Small constant for step size calculation
      const double epsilon = 1e-8;  
      // Number of equations in the system
      const size_t n = 20;          

      // Jacobian matrix view
      gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, n, n);

      // Original function values at y
      double dydt_orig[n];
      ode_fun(t, y, dydt_orig, params);

      // Compute each partial derivative with respect to y using adaptive central differences
      for (size_t i = 0; i < n; i++) {
            // Calculate the adaptive step size h_i for y[i]
            double h_i = sqrt(epsilon) * std::max(std::abs(y[i]), 1.0);
            // double h_i  = epsilon;

            // Perturb y[i] forward and backward by h_i
            double y_perturbed_forward[n], y_perturbed_backward[n];
            for (size_t j = 0; j < n; j++) {
                  y_perturbed_forward[j] = y[j];  // Copy original y values
                  y_perturbed_backward[j] = y[j];
            }
            y_perturbed_forward[i] += h_i;
            y_perturbed_backward[i] -= h_i;

            // Calculate the function values for the forward and backward perturbed states
            double dydt_forward[n], dydt_backward[n];
            ode_fun(t, y_perturbed_forward, dydt_forward, params);
            ode_fun(t, y_perturbed_backward, dydt_backward, params);

            // Compute partial derivatives for each function with respect to y[i]
            for (size_t k = 0; k < n; k++) {
                  double derivative = (dydt_forward[k] - dydt_backward[k]) / (2 * h_i);
                  gsl_matrix_set(&dfdy_mat.matrix, k, i, derivative);
            }
      }

      // Compute partial derivatives with respect to time t using central differences
      double h_t = sqrt(epsilon) * std::max(std::abs(t), 1.0);
      // double h_t  = epsilon;
      // Calculate function values at t + h_t and t - h_t
      double dydt_time_forward[n], dydt_time_backward[n];
      ode_fun(t + h_t, y, dydt_time_forward, params);
      ode_fun(t - h_t, y, dydt_time_backward, params);

      // Compute time derivatives
      for (size_t i = 0; i < n; i++) {
            dfdt[i] = (dydt_time_forward[i] - dydt_time_backward[i]) / (2 * h_t);
      }

      return GSL_SUCCESS;
}







/* ODE for Xn */
int ode_Xn(double t, const double y[], double f[], void *params){
      auto* p     = static_cast<Xn_ODE_Params*>(params);
      double Tnu  = (*(p->interp_Tnu))(t);
      double Tg   = (*(p->interp_Tg))(t);
      double r11  = (*(p->interp_r11))(t);
      double rb11 = (*(p->interp_rb11))(t);
      double H    = (*(p->interp_H))(t);


      /* PDG24 average neutron lifetime in seconds divided by hbar to bring it into units of MeV^-1 */
      double tau_n            = _taunsec_/(_hbar_*pow(10.,-6.));
      /* get A by matching 0. temp rate formula to experimental value of tau_n */
      double A    = 1./(tau_n*0.05700456948361);


      /* partial rates of neutron destruction */

      /* FD interpolation functions */
      double Gamma_nnu_pe, Gamma_ne_pnu1, Gamma_ne_pnu2, Gamma_n_penu1, Gamma_n_penu2, Gamma_n_penu3, Gamma_n_penu4;
      Gamma_nnu_pe            = r11*weakrates2D.Gamma_nnu_pe(Tnu, Tg);
      Gamma_ne_pnu1           = weakrates1D.Gamma_ne_pnu1(Tg);
      Gamma_ne_pnu2           = rb11*weakrates2D.Gamma_ne_pnu2(Tnu, Tg);
      Gamma_n_penu1           = 1./(A*tau_n);
      Gamma_n_penu2           = weakrates1D.Gamma_n_penu2(Tg);
      Gamma_n_penu3           = rb11*weakrates1D.Gamma_n_penu3(Tnu);
      Gamma_n_penu4           = rb11*weakrates2D.Gamma_n_penu4(Tnu, Tg);
      double Gamma_np;
      Gamma_np                = A*(Gamma_nnu_pe + Gamma_ne_pnu1 - Gamma_ne_pnu2 + Gamma_n_penu1 - Gamma_n_penu2 - Gamma_n_penu3 + Gamma_n_penu4);



      /* partial rates of neutron production (inverse direction) */

      /* FD interpolation functions */
      double Gamma_pe_nnu1, Gamma_pe_nnu2, Gamma_pnu_ne, Gamma_penu_n;
      Gamma_pe_nnu1           = weakrates1D.Gamma_pe_nnu1(Tg);
      Gamma_pe_nnu2           = r11*weakrates2D.Gamma_pe_nnu2(Tnu, Tg);
      Gamma_pnu_ne            = rb11*weakrates2D.Gamma_pnu_ne(Tnu, Tg);
      Gamma_penu_n            = rb11*weakrates2D.Gamma_penu_n(Tnu, Tg);
      double Gamma_pn;
      Gamma_pn                = A*(Gamma_pe_nnu1 - Gamma_pe_nnu2 + Gamma_pnu_ne + Gamma_penu_n);


      f[0] = 1./(H*t)*(weakratescorrection.Gamma_pn_Corr(Tg)*Gamma_pn*(1.-y[0]) - weakratescorrection.Gamma_np_Corr(Tg)*Gamma_np*y[0]);

      return GSL_SUCCESS;
}




/*
analytically calculate the osc. average values for (r - rb)
*/
array<double, 3> get_dn_alpha_osc_ave_from_r_rbar_fb(double normalize, complx r_fb[3][3], complx rbar_fb[3][3], void *params_ave){
      std::vector<double> params_var= *(std::vector<double> *)params_ave;
      double theta12, theta13, theta23;

      double r11,r22,r12re,r13re,r23re,r33;
      double rb11,rb22,rb12re,rb13re,rb23re,rb33;


      r11=r_fb[0][0].real(); r22=r_fb[1][1].real(); r33=r_fb[2][2].real();
      r12re=r_fb[0][1].real(); r13re=r_fb[0][2].real(); r23re=r_fb[1][2].real(); 

      rb11=rbar_fb[0][0].real(); rb22=rbar_fb[1][1].real(); rb33=rbar_fb[2][2].real();
      rb12re=rbar_fb[0][1].real(); rb13re=rbar_fb[0][2].real(); rb23re=rbar_fb[1][2].real(); 


      theta12     = params_var[0];
      theta13     = params_var[1];
      theta23     = params_var[2];


      double nnbar11, nnbar22, nnbar33;


      nnbar11     = (r11 - rb11)*pow(cos(theta13),4)*pow(sin(theta12),4) + (r11 - rb11)*pow(sin(theta13),4) + 
                        2*cos(theta13)*pow(sin(theta13),3)*((r13re - rb13re)*cos(theta23) + (r12re - rb12re)*sin(theta23)) - 
                        2*pow(cos(theta13),3)*pow(sin(theta12),4)*sin(theta13)*
                        ((r13re - rb13re)*cos(theta23) + (r12re - rb12re)*sin(theta23)) + 
                        pow(cos(theta12),2)*pow(cos(theta13),2)*pow(sin(theta12),2)*
                        (r22 + r33 - rb22 - rb33 + (r22 - r33 - rb22 + rb33)*cos(2*theta23) - 
                        2*(r23re - rb23re)*sin(2*theta23)) - 
                        ((11 - 4*cos(2*theta12) + cos(4*theta12))*pow(cos(theta13),2)*pow(sin(theta13),2)*
                        (-r22 - r33 + rb22 + rb33 + (r22 - r33 - rb22 + rb33)*cos(2*theta23) - 
                        2*(r23re - rb23re)*sin(2*theta23)))/16. + 
                        pow(cos(theta12),4)*pow(cos(theta13),2)*
                        ((r11 - rb11)*pow(cos(theta13),2) + (r33 - rb33)*pow(cos(theta23),2)*pow(sin(theta13),2) + 
                        r22*pow(sin(theta13),2)*pow(sin(theta23),2) - rb22*pow(sin(theta13),2)*pow(sin(theta23),2) - 
                        2*r12re*cos(theta13)*sin(theta13)*sin(theta23) + rb12re*sin(2*theta13)*sin(theta23) + 
                        cos(theta23)*((-r13re + rb13re)*sin(2*theta13) - 2*rb23re*pow(sin(theta13),2)*sin(theta23)) + 
                        r23re*pow(sin(theta13),2)*sin(2*theta23)) + 
                        cos(theta12)*pow(cos(theta13),2)*pow(sin(theta12),3)*
                        (2*cos(theta13)*((r12re - rb12re)*cos(theta23) + (-r13re + rb13re)*sin(theta23)) + 
                        sin(theta13)*(-2*(r23re - rb23re)*cos(2*theta23) + (-r22 + r33 + rb22 - rb33)*sin(2*theta23))) + 
                        pow(cos(theta12),3)*pow(cos(theta13),2)*sin(theta12)*
                        (-2*cos(theta13)*((r12re - rb12re)*cos(theta23) + (-r13re + rb13re)*sin(theta23)) + 
                        sin(theta13)*(2*(r23re - rb23re)*cos(2*theta23) + (r22 - r33 - rb22 + rb33)*sin(2*theta23)));


      nnbar22     = (r22 - rb22)*pow(cos(theta23),4)*pow(sin(theta12),4) + 
                        (r33 - rb33)*pow(cos(theta23),2)*(pow(cos(theta13),4) + pow(sin(theta12),4)*(1 + pow(sin(theta13),4)))*
                        pow(sin(theta23),2) + ((r11 - rb11)*(11 - 4*cos(2*theta12) + cos(4*theta12))*pow(sin(2*theta13),2)*
                        pow(sin(theta23),2) + 8*(-8*(r12re - rb12re)*cos(theta13)*pow(sin(theta12),4)*pow(sin(theta13),3)*
                        pow(sin(theta23),3) + 4*(r22 - rb22)*pow(cos(theta13),4)*pow(sin(theta23),4) + 
                        4*(r22 - rb22)*pow(sin(theta12),4)*pow(sin(theta13),4)*pow(sin(theta23),4) + 
                        (3*r22 - 2*r33 - 3*rb22 + 2*rb33)*pow(sin(2*theta12),2)*pow(sin(theta13),2)*
                        pow(sin(2*theta23),2) + 8*(r12re - rb12re)*pow(cos(theta13),3)*pow(sin(theta23),3)*sin(theta13)))
                        /32. - 2*(r23re - rb23re)*pow(cos(theta23),3)*pow(sin(theta12),4)*sin(theta23) + 
                        2*cos(theta23)*pow(sin(theta23),2)*(-((r13re - rb13re)*cos(theta13)*pow(sin(theta12),4)*
                        pow(sin(theta13),3)) + (r13re - rb13re)*pow(cos(theta13),3)*sin(theta13) + 
                        (r23re - rb23re)*pow(cos(theta13),4)*sin(theta23) + 
                        (r23re - rb23re)*pow(sin(theta12),4)*pow(sin(theta13),4)*sin(theta23)) + 
                        pow(cos(theta12),4)*((r22 - rb22)*pow(cos(theta23),4) + 
                        ((r33 - rb33)*(11 - 4*cos(2*theta13) + cos(4*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/8. - 
                        2*(r23re - rb23re)*pow(cos(theta23),3)*sin(theta23) + 
                        2*cos(theta23)*pow(sin(theta13),3)*pow(sin(theta23),2)*
                        ((-r13re + rb13re)*cos(theta13) + (r23re - rb23re)*sin(theta13)*sin(theta23)) + 
                        pow(sin(theta13),2)*pow(sin(theta23),2)*
                        ((r11 - rb11)*pow(cos(theta13),2) + 
                        sin(theta23)*((-r12re + rb12re)*sin(2*theta13) + (r22 - rb22)*pow(sin(theta13),2)*sin(theta23))))\
                        + cos(theta12)*pow(sin(theta12),3)*(-2*cos(theta13)*
                        ((r12re - rb12re)*pow(cos(theta23),3) - 
                        3*(r12re - rb12re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                        (r13re - rb13re)*pow(sin(theta13),2)*pow(sin(theta23),3) + 
                        (r13re - rb13re)*(-2 + cos(2*theta13))*pow(cos(theta23),2)*sin(theta23)) + 
                        sin(theta13)*(2*(r23re - rb23re)*pow(cos(theta23),4) + 
                        3*(r23re - rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2) + 
                        (-2*r22 + 3*r33 + 2*rb22 - 3*rb33 + (2*r22 - r33 - 2*rb22 + rb33)*cos(2*theta13))*cos(theta23)*
                        pow(sin(theta23),3) + 2*(r23re - rb23re)*pow(sin(theta13),2)*pow(sin(theta23),4) + 
                        (4*r22 - 3*r33 - 4*rb22 + 3*rb33 + (r33 - rb33)*cos(2*theta13))*pow(cos(theta23),3)*sin(theta23))\
                        - (r11 - rb11)*pow(cos(theta13),2)*sin(theta13)*sin(2*theta23)) + 
                        pow(cos(theta12),3)*sin(theta12)*(2*cos(theta13)*
                        ((r12re - rb12re)*pow(cos(theta23),3) - 
                        3*(r12re - rb12re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                        (r13re - rb13re)*pow(sin(theta13),2)*pow(sin(theta23),3) + 
                        (r13re - rb13re)*(-2 + cos(2*theta13))*pow(cos(theta23),2)*sin(theta23)) - 
                        sin(theta13)*(2*(r23re - rb23re)*pow(cos(theta23),4) + 
                        3*(r23re - rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2) + 
                        (-2*r22 + 3*r33 + 2*rb22 - 3*rb33 + (2*r22 - r33 - 2*rb22 + rb33)*cos(2*theta13))*cos(theta23)*
                        pow(sin(theta23),3) + 2*(r23re - rb23re)*pow(sin(theta13),2)*pow(sin(theta23),4) + 
                        (4*r22 - 3*r33 - 4*rb22 + 3*rb33 + (r33 - rb33)*cos(2*theta13))*pow(cos(theta23),3)*sin(theta23))\
                        + (r11 - rb11)*pow(cos(theta13),2)*sin(theta13)*sin(2*theta23)) + 
                        2*pow(cos(theta12),2)*pow(sin(theta12),2)*
                        ((r11 - rb11)*pow(cos(theta13),2)*pow(cos(theta23),2) + 
                        (r33 - rb33)*pow(cos(theta23),4)*pow(sin(theta13),2) - 
                        6*(r23re - rb23re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),3) + 
                        r33*pow(sin(theta13),2)*pow(sin(theta23),4) - rb33*pow(sin(theta13),2)*pow(sin(theta23),4) - 
                        3*r12re*pow(cos(theta23),2)*sin(2*theta13)*sin(theta23) + 
                        pow(cos(theta23),3)*((-r13re + rb13re)*sin(2*theta13) + 
                        6*(r23re - rb23re)*pow(sin(theta13),2)*sin(theta23)) + 
                        r13re*sin(2*theta13)*sin(theta23)*sin(2*theta23) + 
                        cos(theta13)*sin(theta13)*(3*rb12re*cos(theta23) - 2*rb13re*sin(theta23))*sin(2*theta23));


      nnbar33     = 2*pow(cos(theta13),3)*pow(cos(theta23),2)*sin(theta13)*
                        ((r13re - rb13re)*cos(theta23) + (r12re - rb12re)*sin(theta23)) + 
                        2*pow(cos(theta12),3)*sin(theta12)*sin(theta13)*
                        ((r23re - rb23re)*pow(cos(theta23),4)*pow(sin(theta13),2) + 
                        (3*(r23re - rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/2. + 
                        ((-3*r22 + 4*r33 + 3*rb22 - 4*rb33 + (r22 - rb22)*cos(2*theta13))*cos(theta23)*pow(sin(theta23),3))/
                        2. + (r23re - rb23re)*pow(sin(theta23),4) + 
                        pow(cos(theta23),3)*(r22 - rb22 + (r22 - 2*r33 - rb22 + 2*rb33)*pow(sin(theta13),2))*sin(theta23)) + 
                        pow(cos(theta12),4)*((r33 - rb33)*pow(cos(theta23),4)*pow(sin(theta13),4) + 
                        ((r22 - rb22)*(11 - 4*cos(2*theta13) + cos(4*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/8. - 
                        2*(r23re - rb23re)*cos(theta23)*pow(sin(theta23),3) + (r33 - rb33)*pow(sin(theta23),4) + 
                        2*(r23re - rb23re)*pow(cos(theta23),3)*pow(sin(theta13),4)*sin(theta23)) + 
                        (4*(r33 - rb33)*pow(cos(theta23),4)*pow(sin(theta12),4)*pow(sin(theta13),4) - 
                        8*(r23re - rb23re)*cos(theta23)*pow(sin(theta12),4)*pow(sin(theta23),3) + 
                        (-2*r22 + 3*r33 + 2*rb22 - 3*rb33)*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2) + 
                        pow(sin(theta12),4)*(4*(r33 - rb33)*pow(sin(theta23),4) + 
                        (r22 - rb22)*(1 + pow(sin(theta13),4))*pow(sin(2*theta23),2)) + 
                        8*(r23re - rb23re)*pow(cos(theta23),3)*pow(sin(theta12),4)*pow(sin(theta13),4)*sin(theta23))/4. + 
                        2*pow(cos(theta12),2)*pow(sin(theta12),2)*
                        ((r22 - rb22)*pow(cos(theta23),4)*pow(sin(theta13),2) - 
                        6*(r23re - rb23re)*pow(cos(theta23),3)*pow(sin(theta13),2)*sin(theta23) + 
                        pow(sin(theta23),3)*(rb12re*sin(2*theta13) + (r22 - rb22)*pow(sin(theta13),2)*sin(theta23)) + 
                        cos(theta23)*pow(sin(theta23),2)*((-r13re + rb13re)*sin(2*theta13) + 
                        6*(r23re - rb23re)*pow(sin(theta13),2)*sin(theta23))) + 
                        pow(cos(theta13),4)*pow(cos(theta23),2)*
                        ((r33 - rb33)*pow(cos(theta23),2) + (r22 - rb22)*pow(sin(theta23),2) + (r23re - rb23re)*sin(2*theta23))
                        + ((r11 - rb11)*pow(cos(theta13),2)*((7 + cos(4*theta12))*pow(cos(theta23),2)*pow(sin(theta13),2) + 
                        2*pow(sin(2*theta12),2)*pow(sin(theta23),2) - sin(4*theta12)*sin(theta13)*sin(2*theta23)))/4. - 
                        2*cos(theta13)*(pow(cos(theta12),4)*pow(cos(theta23),2)*pow(sin(theta13),3)*
                        ((r13re - rb13re)*cos(theta23) + (r12re - rb12re)*sin(theta23)) + 
                        pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta13),3)*
                        ((r13re - rb13re)*cos(theta23) + (r12re - rb12re)*sin(theta23)) - 
                        pow(cos(theta12),2)*pow(sin(theta12),2)*sin(theta13)*sin(theta23)*
                        (r12re - 2*rb12re + (3*r12re - 2*rb12re)*cos(2*theta23) - 2*(r13re - rb13re)*sin(2*theta23)) + 
                        cos(theta12)*pow(sin(theta12),3)*(-((r12re - rb12re)*pow(cos(theta23),3)*pow(sin(theta13),2)) + 
                        (r12re - 2*rb12re + rb12re*cos(2*theta13))*cos(theta23)*pow(sin(theta23),2) + 
                        (-r13re + rb13re)*pow(sin(theta23),3) + 
                        3*(r13re - rb13re)*pow(cos(theta23),2)*pow(sin(theta13),2)*sin(theta23) + 
                        r12re*pow(sin(theta13),2)*sin(theta23)*sin(2*theta23)) + 
                        pow(cos(theta12),3)*sin(theta12)*((r12re - rb12re)*pow(cos(theta23),3)*pow(sin(theta13),2) + 
                        (-2*r12re + rb12re + r12re*cos(2*theta13))*cos(theta23)*pow(sin(theta23),2) + 
                        (r13re - rb13re)*pow(sin(theta23),3) - 
                        3*(r13re - rb13re)*pow(cos(theta23),2)*pow(sin(theta13),2)*sin(theta23) + 
                        rb12re*pow(sin(theta13),2)*sin(theta23)*sin(2*theta23))) + 
                        (cos(theta12)*pow(sin(theta12),3)*sin(theta13)*
                        (2*(r23re - rb23re)*cos(2*theta13 - 4*theta23) + 2*(r23re - rb23re)*cos(2*theta13 - 2*theta23) + 
                        4*r23re*cos(2*theta23) - 4*rb23re*cos(2*theta23) - 12*r23re*cos(4*theta23) + 
                        12*rb23re*cos(4*theta23) + 2*r23re*cos(2*(theta13 + theta23)) - 
                        2*rb23re*cos(2*(theta13 + theta23)) + 2*r23re*cos(2*theta13 + 4*theta23) - 
                        2*rb23re*cos(2*theta13 + 4*theta23) - r22*sin(2*theta13 - 4*theta23) + 
                        r33*sin(2*theta13 - 4*theta23) + rb22*sin(2*theta13 - 4*theta23) - 
                        rb33*sin(2*theta13 - 4*theta23) + 2*r33*sin(2*theta13 - 2*theta23) - 
                        2*rb33*sin(2*theta13 - 2*theta23) - 4*r33*sin(2*theta23) + 4*rb33*sin(2*theta23) - 
                        6*r22*sin(4*theta23) + 6*r33*sin(4*theta23) + 6*rb22*sin(4*theta23) - 6*rb33*sin(4*theta23) - 
                        2*r33*sin(2*(theta13 + theta23)) + 2*rb33*sin(2*(theta13 + theta23)) + 
                        r22*sin(2*theta13 + 4*theta23) - r33*sin(2*theta13 + 4*theta23) - 
                        rb22*sin(2*theta13 + 4*theta23) + rb33*sin(2*theta13 + 4*theta23)))/8.;


      return {normalize*nnbar11 , normalize*nnbar22 , normalize*nnbar33 };

}




/*
analytically calculate the osc. average values for (r + rb)
*/
array<double, 3> get_np_alpha_osc_ave_from_r_rbar_fb(double normalize, complx r_fb[3][3], complx rbar_fb[3][3], void *params_ave){
      std::vector<double> params_var= *(std::vector<double> *)params_ave;
      double theta12, theta13, theta23;

      double r11,r22,r12re,r13re,r23re,r33;
      double rb11,rb22,rb12re,rb13re,rb23re,rb33;


      r11=r_fb[0][0].real(); r22=r_fb[1][1].real(); r33=r_fb[2][2].real();
      r12re=r_fb[0][1].real(); r13re=r_fb[0][2].real(); r23re=r_fb[1][2].real(); 

      rb11=rbar_fb[0][0].real(); rb22=rbar_fb[1][1].real(); rb33=rbar_fb[2][2].real();
      rb12re=rbar_fb[0][1].real(); rb13re=rbar_fb[0][2].real(); rb23re=rbar_fb[1][2].real(); 


      theta12     = params_var[0];
      theta13     = params_var[1];
      theta23     = params_var[2];


      double np11, np22, np33;


      np11   = (r11 + rb11)*pow(cos(theta13),4)*pow(sin(theta12),4) + (r11 + rb11)*pow(sin(theta13),4) + 
                  2*cos(theta13)*pow(sin(theta13),3)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) - 
                  2*pow(cos(theta13),3)*pow(sin(theta12),4)*sin(theta13)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) + 
                  pow(cos(theta12),2)*pow(cos(theta13),2)*pow(sin(theta12),2)*(r22 + r33 + rb22 + rb33 + (r22 - r33 + rb22 - rb33)*cos(2*theta23) - 2*(r23re + rb23re)*sin(2*theta23)) + 
                  ((11 - 4*cos(2*theta12) + cos(4*theta12))*pow(cos(theta13),2)*pow(sin(theta13),2)*
                  (r22 + r33 + rb22 + rb33 + (-r22 + r33 - rb22 + rb33)*cos(2*theta23) + 2*(r23re + rb23re)*sin(2*theta23)))/16. + 
                  pow(cos(theta12),3)*pow(cos(theta13),2)*sin(theta12)*(-2*(r12re + rb12re)*cos(theta13)*cos(theta23) + 2*(r23re + rb23re)*cos(2*theta23)*sin(theta13) + 
                  2*(r13re + rb13re)*cos(theta13)*sin(theta23) + (r22 - r33 + rb22 - rb33)*sin(theta13)*sin(2*theta23)) + 
                  pow(cos(theta12),4)*pow(cos(theta13),2)*((r11 + rb11)*pow(cos(theta13),2) - sin(2*theta13)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) + 
                  pow(sin(theta13),2)*((r33 + rb33)*pow(cos(theta23),2) + (r22 + rb22)*pow(sin(theta23),2) + (r23re + rb23re)*sin(2*theta23))) + 
                  cos(theta12)*pow(cos(theta13),2)*pow(sin(theta12),3)*(2*cos(theta13)*((r12re + rb12re)*cos(theta23) - (r13re + rb13re)*sin(theta23)) - 
                  sin(theta13)*(2*(r23re + rb23re)*cos(2*theta23) + (r22 - r33 + rb22 - rb33)*sin(2*theta23)));

      np22  = (r22 + rb22)*pow(cos(theta23),4)*pow(sin(theta12),4) + r11*pow(cos(theta13),2)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                  rb11*pow(cos(theta13),2)*pow(sin(theta13),2)*pow(sin(theta23),2) + r11*pow(cos(theta13),2)*pow(sin(theta12),4)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                  rb11*pow(cos(theta13),2)*pow(sin(theta12),4)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                  (r33 + rb33)*pow(cos(theta23),2)*(pow(cos(theta13),4) + pow(sin(theta12),4)*(1 + pow(sin(theta13),4)))*pow(sin(theta23),2) - 
                  2*r12re*cos(theta13)*pow(sin(theta12),4)*pow(sin(theta13),3)*pow(sin(theta23),3) - 2*rb12re*cos(theta13)*pow(sin(theta12),4)*pow(sin(theta13),3)*pow(sin(theta23),3) + 
                  r22*pow(cos(theta13),4)*pow(sin(theta23),4) + rb22*pow(cos(theta13),4)*pow(sin(theta23),4) + r22*pow(sin(theta12),4)*pow(sin(theta13),4)*pow(sin(theta23),4) + 
                  rb22*pow(sin(theta12),4)*pow(sin(theta13),4)*pow(sin(theta23),4) + (3*r22*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/4. - 
                  (r33*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/2. + (3*rb22*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/4. - 
                  (rb33*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/2. + 2*r12re*pow(cos(theta13),3)*pow(sin(theta23),3)*sin(theta13) + 
                  2*rb12re*pow(cos(theta13),3)*pow(sin(theta23),3)*sin(theta13) - 2*(r23re + rb23re)*pow(cos(theta23),3)*pow(sin(theta12),4)*sin(theta23) + 
                  2*cos(theta23)*pow(sin(theta23),2)*(-((r13re + rb13re)*cos(theta13)*pow(sin(theta12),4)*pow(sin(theta13),3)) + (r13re + rb13re)*pow(cos(theta13),3)*sin(theta13) + 
                  (r23re + rb23re)*pow(cos(theta13),4)*sin(theta23) + (r23re + rb23re)*pow(sin(theta12),4)*pow(sin(theta13),4)*sin(theta23)) - 
                  2*cos(theta12)*pow(sin(theta12),3)*((r11 + rb11)*cos(theta23)*pow(cos(theta13),2)*sin(theta13)*sin(theta23) + 
                  cos(theta13)*((r12re + rb12re)*pow(cos(theta23),3) - 3*(r12re + rb12re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                  (r13re + rb13re)*pow(sin(theta13),2)*pow(sin(theta23),3) + (r13re + rb13re)*(-2 + cos(2*theta13))*pow(cos(theta23),2)*sin(theta23)) - 
                  (sin(theta13)*(2*(r23re + rb23re)*pow(cos(theta23),4) + 3*(r23re + rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2) + 
                  (-2*r22 + 3*r33 - 2*rb22 + 3*rb33 + (2*r22 - r33 + 2*rb22 - rb33)*cos(2*theta13))*cos(theta23)*pow(sin(theta23),3) + 
                  2*(r23re + rb23re)*pow(sin(theta13),2)*pow(sin(theta23),4) + (4*r22 - 3*r33 + 4*rb22 - 3*rb33 + (r33 + rb33)*cos(2*theta13))*pow(cos(theta23),3)*sin(theta23)))
                  /2.) + pow(cos(theta12),4)*((r22 + rb22)*pow(cos(theta23),4) + 
                  ((r33 + rb33)*(11 - 4*cos(2*theta13) + cos(4*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/8. - 2*(r23re + rb23re)*pow(cos(theta23),3)*sin(theta23) + 
                  2*cos(theta23)*pow(sin(theta13),3)*pow(sin(theta23),2)*(-((r13re + rb13re)*cos(theta13)) + (r23re + rb23re)*sin(theta13)*sin(theta23)) + 
                  pow(sin(theta13),2)*pow(sin(theta23),2)*((r11 + rb11)*pow(cos(theta13),2) + 
                  sin(theta23)*(-((r12re + rb12re)*sin(2*theta13)) + (r22 + rb22)*pow(sin(theta13),2)*sin(theta23)))) + 
                  pow(cos(theta12),3)*sin(theta12)*(2*cos(theta13)*((r12re + rb12re)*pow(cos(theta23),3) - 3*(r12re + rb12re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),2) + 
                  (r13re + rb13re)*pow(sin(theta13),2)*pow(sin(theta23),3) + (r13re + rb13re)*(-2 + cos(2*theta13))*pow(cos(theta23),2)*sin(theta23)) - 
                  sin(theta13)*(2*(r23re + rb23re)*pow(cos(theta23),4) + 3*(r23re + rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2) + 
                  (-2*r22 + 3*r33 - 2*rb22 + 3*rb33 + (2*r22 - r33 + 2*rb22 - rb33)*cos(2*theta13))*cos(theta23)*pow(sin(theta23),3) + 
                  2*(r23re + rb23re)*pow(sin(theta13),2)*pow(sin(theta23),4) + (4*r22 - 3*r33 + 4*rb22 - 3*rb33 + (r33 + rb33)*cos(2*theta13))*pow(cos(theta23),3)*sin(theta23)) + 
                  (r11 + rb11)*pow(cos(theta13),2)*sin(theta13)*sin(2*theta23)) + 
                  2*pow(cos(theta12),2)*pow(sin(theta12),2)*((r11 + rb11)*pow(cos(theta13),2)*pow(cos(theta23),2) + (r33 + rb33)*pow(cos(theta23),4)*pow(sin(theta13),2) - 
                  6*(r23re + rb23re)*cos(theta23)*pow(sin(theta13),2)*pow(sin(theta23),3) + (r33 + rb33)*pow(sin(theta13),2)*pow(sin(theta23),4) - 
                  3*(r12re + rb12re)*pow(cos(theta23),2)*sin(2*theta13)*sin(theta23) - 
                  pow(cos(theta23),3)*((r13re + rb13re)*sin(2*theta13) - 6*(r23re + rb23re)*pow(sin(theta13),2)*sin(theta23)) + 
                  (r13re + rb13re)*sin(2*theta13)*sin(theta23)*sin(2*theta23));

      np33  = r33*pow(cos(theta23),4)*pow(sin(theta12),4)*pow(sin(theta13),4) + rb33*pow(cos(theta23),4)*pow(sin(theta12),4)*pow(sin(theta13),4) + 
                  r22*pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta23),2) + rb22*pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta23),2) + 
                  r22*pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta13),4)*pow(sin(theta23),2) + 
                  rb22*pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta13),4)*pow(sin(theta23),2) - 2*r23re*cos(theta23)*pow(sin(theta12),4)*pow(sin(theta23),3) - 
                  2*rb23re*cos(theta23)*pow(sin(theta12),4)*pow(sin(theta23),3) + r33*pow(sin(theta12),4)*pow(sin(theta23),4) + rb33*pow(sin(theta12),4)*pow(sin(theta23),4) - 
                  (r22*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/2. + (3*r33*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/4. - 
                  (rb22*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/2. + (3*rb33*pow(sin(2*theta12),2)*pow(sin(theta13),2)*pow(sin(2*theta23),2))/4. + 
                  2*r23re*pow(cos(theta23),3)*pow(sin(theta12),4)*pow(sin(theta13),4)*sin(theta23) + 2*rb23re*pow(cos(theta23),3)*pow(sin(theta12),4)*pow(sin(theta13),4)*sin(theta23) + 
                  2*pow(cos(theta13),3)*pow(cos(theta23),2)*sin(theta13)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) + 
                  2*pow(cos(theta12),3)*sin(theta12)*sin(theta13)*((r23re + rb23re)*pow(cos(theta23),4)*pow(sin(theta13),2) + 
                  (3*(r23re + rb23re)*(-3 + cos(2*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/2. + 
                  ((-3*r22 + 4*r33 - 3*rb22 + 4*rb33 + (r22 + rb22)*cos(2*theta13))*cos(theta23)*pow(sin(theta23),3))/2. + (r23re + rb23re)*pow(sin(theta23),4) + 
                  pow(cos(theta23),3)*(r22 + rb22 + (r22 - 2*r33 + rb22 - 2*rb33)*pow(sin(theta13),2))*sin(theta23)) + 
                  pow(cos(theta12),4)*((r33 + rb33)*pow(cos(theta23),4)*pow(sin(theta13),4) + 
                  ((r22 + rb22)*(11 - 4*cos(2*theta13) + cos(4*theta13))*pow(cos(theta23),2)*pow(sin(theta23),2))/8. - 2*(r23re + rb23re)*cos(theta23)*pow(sin(theta23),3) + 
                  (r33 + rb33)*pow(sin(theta23),4) + 2*(r23re + rb23re)*pow(cos(theta23),3)*pow(sin(theta13),4)*sin(theta23)) + 
                  2*pow(cos(theta12),2)*pow(sin(theta12),2)*((r22 + rb22)*pow(cos(theta23),4)*pow(sin(theta13),2) - 
                  6*(r23re + rb23re)*pow(cos(theta23),3)*pow(sin(theta13),2)*sin(theta23) + 
                  pow(sin(theta23),3)*(-((r12re + rb12re)*sin(2*theta13)) + (r22 + rb22)*pow(sin(theta13),2)*sin(theta23)) + 
                  cos(theta23)*pow(sin(theta23),2)*(-((r13re + rb13re)*sin(2*theta13)) + 6*(r23re + rb23re)*pow(sin(theta13),2)*sin(theta23))) + 
                  pow(cos(theta13),4)*pow(cos(theta23),2)*((r33 + rb33)*pow(cos(theta23),2) + (r22 + rb22)*pow(sin(theta23),2) + (r23re + rb23re)*sin(2*theta23)) + 
                  (r11 + rb11)*pow(cos(theta13),2)*(((7 + cos(4*theta12))*pow(cos(theta23),2)*pow(sin(theta13),2))/4. - 
                  2*cos(theta23)*pow(cos(theta12),3)*sin(theta12)*sin(theta13)*sin(theta23) + 
                  cos(theta12)*pow(sin(theta12),2)*(2*cos(theta12)*pow(sin(theta23),2) + sin(theta12)*sin(theta13)*sin(2*theta23))) - 
                  2*cos(theta13)*(pow(cos(theta12),4)*pow(cos(theta23),2)*pow(sin(theta13),3)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) + 
                  pow(cos(theta23),2)*pow(sin(theta12),4)*pow(sin(theta13),3)*((r13re + rb13re)*cos(theta23) + (r12re + rb12re)*sin(theta23)) + 
                  pow(cos(theta12),3)*sin(theta12)*((r12re + rb12re)*pow(cos(theta23),3)*pow(sin(theta13),2) + 
                  (r12re + rb12re)*(-2 + cos(2*theta13))*cos(theta23)*pow(sin(theta23),2) + (r13re + rb13re)*pow(sin(theta23),3) - 
                  3*(r13re + rb13re)*pow(cos(theta23),2)*pow(sin(theta13),2)*sin(theta23)) - 
                  2*pow(cos(theta12),2)*pow(sin(theta12),2)*sin(theta13)*((r12re + rb12re)*cos(theta23) - (r13re + rb13re)*sin(theta23))*sin(2*theta23) + 
                  cos(theta12)*pow(sin(theta12),3)*(-((r12re + rb12re)*pow(cos(theta23),3)*pow(sin(theta13),2)) + (r12re + rb12re)*cos(theta23)*pow(sin(theta23),2) - 
                  (r13re + rb13re)*pow(sin(theta23),3) + 3*(r13re + rb13re)*pow(cos(theta23),2)*pow(sin(theta13),2)*sin(theta23) + 
                  (r12re + rb12re)*pow(sin(theta13),2)*sin(theta23)*sin(2*theta23))) + 
                  (cos(theta12)*pow(sin(theta12),3)*sin(theta13)*(2*(r23re + rb23re)*cos(2*theta13 - 4*theta23) + 2*(r23re + rb23re)*cos(2*theta13 - 2*theta23) + 
                  4*r23re*cos(2*theta23) + 4*rb23re*cos(2*theta23) - 12*r23re*cos(4*theta23) - 12*rb23re*cos(4*theta23) + 2*r23re*cos(2*(theta13 + theta23)) + 
                  2*rb23re*cos(2*(theta13 + theta23)) + 2*r23re*cos(2*theta13 + 4*theta23) + 2*rb23re*cos(2*theta13 + 4*theta23) - r22*sin(2*theta13 - 4*theta23) + 
                  r33*sin(2*theta13 - 4*theta23) - rb22*sin(2*theta13 - 4*theta23) + rb33*sin(2*theta13 - 4*theta23) + 2*r33*sin(2*theta13 - 2*theta23) + 
                  2*rb33*sin(2*theta13 - 2*theta23) - 4*r33*sin(2*theta23) - 4*rb33*sin(2*theta23) - 6*r22*sin(4*theta23) + 6*r33*sin(4*theta23) - 6*rb22*sin(4*theta23) + 
                  6*rb33*sin(4*theta23) - 2*r33*sin(2*(theta13 + theta23)) - 2*rb33*sin(2*(theta13 + theta23)) + r22*sin(2*theta13 + 4*theta23) - r33*sin(2*theta13 + 4*theta23) + 
                  rb22*sin(2*theta13 + 4*theta23) - rb33*sin(2*theta13 + 4*theta23)))/8.;   


      return {normalize*np11 , normalize*np22 , normalize*np33 };

}





/*
get all flavor of neutrino and anti-neutrino density for given r, rb matrices in flavor space
return an array of 6 scalars
*/
array<double, 6> get_rho_alpha_from_r_rbar_fb(double znu, complx r_fb[3][3], complx rbar_fb[3][3]){
      double norm= 7.*pow(_PI_,2.)/240.*pow(znu,4.);

      double r11,r22,r33;
      double rb11,rb22,rb33;

      r11=r_fb[0][0].real(); r22=r_fb[1][1].real(); r33=r_fb[2][2].real();
      rb11=rbar_fb[0][0].real(); rb22=rbar_fb[1][1].real(); rb33=rbar_fb[2][2].real();


      return {norm*r11, norm*r22, norm*r33, norm*rb11, norm*rb22, norm*rb33};
}




/*
Define the 1D evolution of zgamma after neutrino decoupling
*/
int ode_zg(double t, const double y[], double f[], void *params){
      double rhoe_pe, drhoe, rhog, drhog;
      double Tg;
      Tg    = _Tref_*y[0]/t;

      rhog  = pow(_PI_,2.)/15.*pow(Tg,4.);
      drhog = pow(_PI_,2.)/15.*4.*pow(Tg,3.);

      rhoe_pe     = density.rhoe(Tg) + density.pe(Tg);
      drhoe       = density.drhoe(Tg);

      f[0]  = y[0]/t - (4*rhog + 3*rhoe_pe) / (_Tref_*(drhoe+drhog)) ;

      return GSL_SUCCESS;
}


/* get Yp according to our fitting formula */
double YP(double nnbar11, double rho_e, double rhobar_e, double zg, double Neff){
      double Ypref = 0.247046;
      double dNeff = Neff - 3.044;

      double rho_free   = 2.*7./8.*pow(_PI_,2.)/30.;

      double xie        = 6.*nnbar11*pow( (rho_e+rhobar_e)/rho_free,-3./4.);
      double TnueTg     = pow((rho_e+rhobar_e)/rho_free,1./4.)*1./zg;

      double dT         = TnueTg - pow(11./4.,-1./3.);

      return Ypref * exp(-0.961*xie - 0.185*pow(xie,2.) + 0.00102*pow(xie,3.) 
                        +0.016/0.3*dNeff-0.682*dT+0.72*pow(dT,2.)+0.0108809*xie*dNeff-3.472176*xie*dT
                        -0.90547497*pow(xie,2.)*dT+2.22854596976*xie*pow(dT,2.)+0.0160705967*dNeff*dT);
}

/* get DH according to our fitting formula */
double DH(double nnbar11, double rho_e, double rhobar_e, double zg, double Neff){
      double DHref = 2.43674;
      double dNeff = Neff - 3.044;

      double rho_free   = 2.*7./8.*pow(_PI_,2.)/30.;

      double xie        = 6.*nnbar11*pow( (rho_e+rhobar_e)/rho_free,-3./4.);
      double TnueTg     = pow((rho_e+rhobar_e)/rho_free,1./4.)*1./zg;

      double dT         = TnueTg - pow(11./4.,-1./3.);

      return DHref * exp(-0.528*xie+0.22*pow(xie,2.)-0.052*pow(xie,3.)+0.133*dNeff-0.378*dT+0.562*pow(dT,2.)
                        -0.026*xie*dNeff-1.5*xie*dT-0.0096*dNeff*dT);
}


/*
***************************************
---------------------------------------
solve the ODE system
---------------------------------------
This function evolves the kinetic equation
i) Including neutrino oscillations from Tini -> Tave == neutrino decoupling,
ii) Nuetrino decoupled photon-elextron evolution.
input :
void* params =
{
      dm21, dm31,
      theta12, theta13, theta23,
      xi1,  xi2,  xi3,
      zon, vs, T_ave
}
and returns (for T = T_end if not otherwise specified) :
{     
      xi1ini,  xi2ini,  xi3ini,
      (n-\bar{n})_e/T^3, (n-\bar{n})_\mu/T^3, (n-\bar{n})_\tau/T^3, 
      z_nu^dec, z_gamma^final, N_{eff}, 
      \rho_\alpha, \bar{\rho}_\alpha,
      X_n(T_cm = 1.501 MeV), Yp, D/H 
}
***************************************
*/
array<double, 18> ode_evolve(void* params, string filename="", string flush="no"){
      std::vector<double> params_var= *(std::vector<double> *)params;
      double dm21, dm31, theta12, theta13, theta23;
      double xi1, xi2, xi3;
      double zon, vs, T_ave;

      dm21= params_var[0]; dm31= params_var[1];
      theta12= params_var[2]; theta13= params_var[3]; theta23= params_var[4];
      xi1= params_var[5];  xi2= params_var[6];  xi3= params_var[7];
      zon= params_var[8]; vs= params_var[9]; T_ave= params_var.back();


      double elapsed_seconds_ode_evolve_double, time_limit, length_vec_ode_sol, write_ratio;
      //set time limit for one ODE step in seconds; strongly depends on length_vec_ode_sol
      time_limit = 1;
      //set (inverse) step size for the ODE, I think we need >=2e5 to cover the full parameter space
      length_vec_ode_sol = 2.e5;
      // steps written to file (if safe="yes"): length_vec_ode_sol / write_ratio
      write_ratio = 1e2;

      // Set the custom error handler
      gsl_set_error_handler(&custom_gsl_error_handler);

      if(vs== 0){
            cout << "evolve system in adiabatic approximation " << endl;
      }
      else{
            cout << "evolve full system non-adiabatically " << endl;
            vs = 1;
      }

      if (dm31 > 0){
            cout << "Normal Hierarchy " << endl;
      }
      else{
            cout << "Inverted Hierarchy " << endl;
      }

      if (zon == 1){
            cout << "evolve neutrino and photon temperature " << endl;
      }
      else{
            cout << "assume equal neutrino and photon temperature " << endl;
            zon = 0;
      }
      

      /* 
      calculate alternatice input parameters and print them to terminal
      */
      double xitilde1, xitilde2, xitilde3, xitilde_tot;
      xitilde1    = paramtranslate.get_xi_tilde_alpha_from_xi_alpha(xi1);
      xitilde2    = paramtranslate.get_xi_tilde_alpha_from_xi_alpha(xi2);
      xitilde3    = paramtranslate.get_xi_tilde_alpha_from_xi_alpha(xi3);
      double dn1, dn2, dn3, dn, dNeff;
      dn1    = paramtranslate.get_dn_alpha_from_xi_alpha(xi1);
      dn2    = paramtranslate.get_dn_alpha_from_xi_alpha(xi2);
      dn3    = paramtranslate.get_dn_alpha_from_xi_alpha(xi3);
      dn     = paramtranslate.get_dn_from_xi_alpha(xi1,xi2,xi3);
      xitilde_tot = paramtranslate.get_xitilde_tot_from_dn_tot(dn);
      dNeff = paramtranslate.get_delta_Neff_from_xi_alpha(xi1,xi2,xi3);

      cout << "Input parameters and equivalents: " << endl;
      cout << "-----------------" << endl;
      cout << " \u03BE_ini = " << xi1 << " " << xi2 << " " << xi3 << endl;
      cout << " \u03BE\u0303_ini = " << xitilde1 << " " << xitilde2 << " " << xitilde3 << endl;
      cout << " \u0394n_ini= " << dn1 << " " << dn2 << " " << dn3 << endl;
      cout << " \u0394n    = " << dn << endl;
      cout << " \u03BE\u0303_tot = " << xitilde_tot << endl;
      cout << " \u0394Neff = " << dNeff << endl;
      cout << "-----------------" << endl;
      /* ---- finish print ---- */



      /* ---- set intial conditions ---- */
      /* initial conditions in flavor space */
      complx r_flavor_ini[3][3];
      complx rbar_flavor_ini[3][3];
      double znu_ini, zg_ini;
      double xi_zon_ini[4]= {xi1,xi2,xi3,zon};
      ini_correction(xi_zon_ini, r_flavor_ini, rbar_flavor_ini, znu_ini, zg_ini);

      /* initial conditions translated to the GellMann matrix expansion coefficients */
      complx lvec[9][3][3];
      GellMannMatrices(lvec);
      double r_GellMann_ini[9];
      double rbar_GellMann_ini[9];
      get_GellMann_coefficients_from_mat(r_flavor_ini, lvec, r_GellMann_ini);
      get_GellMann_coefficients_from_mat(rbar_flavor_ini, lvec, rbar_GellMann_ini);
      double y[20] = {  r_GellMann_ini[0], r_GellMann_ini[1], r_GellMann_ini[2], 
                        r_GellMann_ini[3], r_GellMann_ini[4], r_GellMann_ini[5], 
                        r_GellMann_ini[6], r_GellMann_ini[7], r_GellMann_ini[8], 
                        rbar_GellMann_ini[0], rbar_GellMann_ini[1], rbar_GellMann_ini[2], 
                        rbar_GellMann_ini[3], rbar_GellMann_ini[4], rbar_GellMann_ini[5], 
                        rbar_GellMann_ini[6], rbar_GellMann_ini[7], rbar_GellMann_ini[8],
                        znu_ini, zg_ini
                        };


      /* 
      define parameters passed to kinetic eq.
      Basically a reduced set of the same parameters as passed to this function ode_evolve
      */
      vector<double> params_ode     = {dm21, dm31, theta12, theta13, theta23, zon, vs};
      void *ptr_params_ode          = &params_ode;


      /* ---- initialize GSL ODE solver ---- */
      const int dim = 20;
      const double epsabs = pow(10,-6);
      const double epsrel = pow(10,-5);
      const double hstart = pow(10,-5);

      gsl_odeiv2_system sys         = {ode_fun, jac, dim, ptr_params_ode};
      gsl_odeiv2_driver * driver    = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msbdf, hstart, epsabs, epsrel);
                                                                              

      /* --- set time span to solve the kinetic eq. ---- */
      double T_ini= _Tref_; 
      double t0 = _Tref_/T_ini; double t_ave = _Tref_/T_ave;
      double t = t0; double ti;



      // initialize flavor basis matrices
      complx r_fb[3][3], rbar_fb[3][3];


      /*
      initialize aux result vector with certain dimension
      factor to translate from r-\bar{r} to (n-\bar{n})_\alpha/cmT^3.
      (n-\bar{n})/Tcm^3 = 1/znu^3 3 \zeta(3)/(4\pi^2) (r-\bar{r})
      */
      double n_alpha_factor  = 3.*_Zeta3_/(4.*pow(_PI_,2.));
      array<double, 3> dn_alpha_res;
      array<double, 3> np_alpha_res;
      vector<double> myres;
      myres.push_back( _Tref_/t0 );
      if (filename != ""){
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[0][0].real() - rbar_flavor_ini[0][0].real()) );
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[1][1].real() - rbar_flavor_ini[1][1].real()) );
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[2][2].real() - rbar_flavor_ini[2][2].real()) );
            myres.push_back( znu_ini );
            myres.push_back( zg_ini );
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[0][0].real() + rbar_flavor_ini[0][0].real()) );
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[1][1].real() + rbar_flavor_ini[1][1].real()) );
            myres.push_back( pow(znu_ini,3.) * n_alpha_factor * (r_flavor_ini[2][2].real() + rbar_flavor_ini[2][2].real()) );
      }
      
      // to build interpolation needed for Xn evolution
      double H_ini      = pow(_Tref_,2.) * 1./(_MPL_/sqrt(11./2. * pow(zg_ini,4.) + 7./8. * (trace(r_flavor_ini) + trace(rbar_flavor_ini)) * pow(znu_ini,4.)));
      vector<double> myres_x, myres_r11, myres_rb11, myres_H, myres_Tnu, myres_Tg;
      myres_x.push_back(t0);
      myres_Tnu.push_back(_Tref_/t0 * znu_ini);
      myres_Tg.push_back(_Tref_/t0 * zg_ini);
      myres_r11.push_back(y[0] + y[3] + y[8]/sqrt(3));
      myres_rb11.push_back(y[9] + y[12] + y[17]/sqrt(3));
      myres_H.push_back(H_ini);


      // dummy variable to check for timeout
      double timeoutT;
      timeoutT    = 0;

      /* -- variables for the progress bar -- */
      double progress, pos;
      int barWidth= 40;
      progress = 0.0;
      /* -- */
      cout << "start the integration ... " << endl;
      for (int i = 1; i < length_vec_ode_sol; i++){
            /* start clock for timeout determination */
            auto start_ode_evolve   = chrono::system_clock::now();

            /* evolve system one timestep */
            ti = t0 + i * (t_ave - t0) / (length_vec_ode_sol-1);
            int status = gsl_odeiv2_driver_apply (driver, &t, ti, y);
            /* check for success/failure of evolution */
            if (status != GSL_SUCCESS) {
                  if (gsl_error_flag) {
                        gsl_odeiv2_driver_free(driver);
                        gsl_set_error_handler(nullptr);  
                        // Handle the singular matrix error
                        cout << "\n   \033[31mFor \u03BE_ini : " << xi1 << " " << xi2 << " " << xi3 << " GSL Singular Matrix Error: " << gsl_strerror(status) << " at T = " << _Tref_/ti << "\033[39;49m" << endl;
                        return {-99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};
                  } 
                  else {
                        // Handle other GSL errors
                        // Use gsl_strerror to print the error message
                        gsl_odeiv2_driver_free(driver);
                        gsl_set_error_handler(nullptr);  
                        cout << "\n   \033[31mFor \u03BE_ini : " << xi1 << " " << xi2 << " " << xi3 << " GSL ODE Error: " << gsl_strerror(status) << " at T = " << _Tref_/ti << "\033[39;49m" << endl;
                        return {-99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};
                  }
            }

            /* check for timeout */
            auto end_ode_evolve     = chrono::system_clock::now();
            chrono::duration<double> elapsed_seconds_ode_evolve   = end_ode_evolve-start_ode_evolve;
            elapsed_seconds_ode_evolve_double                     = elapsed_seconds_ode_evolve.count();
            if (elapsed_seconds_ode_evolve_double > time_limit){
                  cout << "\n   \033[31mFor \u03BE_ini : " << xi1 << " " << xi2 << " " << xi3 << " reached timeout limit of " << time_limit << " seconds per step at T = " << _Tref_/ti << "\033[39;49m" << endl;
                  timeoutT    = _Tref_/ti;
                  break;
            }
                        
            /* translate result of ODE to flavor basis and get final flavor assymmetry */ 
            double coeff_r[9]        = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]};
            double coeff_rbar[9]     = {y[9], y[10], y[11], y[12], y[13], y[14], y[15], y[16], y[17]};
            get_mat_from_GellMann_coeff(coeff_r, lvec, r_fb);
            get_mat_from_GellMann_coeff(coeff_rbar, lvec, rbar_fb);
            double znu  = y[18];
            double zg   = y[19];
            dn_alpha_res = get_dn_alpha_from_r_rbar_fb(n_alpha_factor*pow(znu,3.), r_fb, rbar_fb);
            np_alpha_res = get_np_alpha_from_r_rbar_fb(n_alpha_factor*pow(znu,3.), r_fb, rbar_fb);
            if (filename != "" && fmod(i, write_ratio) == 0){
                  myres.push_back( _Tref_/ti );
                  myres.push_back( dn_alpha_res[0] );
                  myres.push_back( dn_alpha_res[1] );
                  myres.push_back( dn_alpha_res[2] );
                  myres.push_back( znu );
                  myres.push_back( zg );
                  myres.push_back( np_alpha_res[0] );
                  myres.push_back( np_alpha_res[1] );
                  myres.push_back( np_alpha_res[2] );
            }

            double H    = pow(_Tref_/ti,2.) * 1./(_MPL_/sqrt(11./2. * pow(zg,4.) + 7./8. * 3. * (y[0] + y[9]) * pow(znu,4.)));
            if (fmod(i, write_ratio) == 0 ){
                  myres_x.push_back(ti);
                  myres_Tnu.push_back(_Tref_/ti * znu);
                  myres_Tg.push_back(_Tref_/ti * zg);
                  myres_r11.push_back(y[0] + y[3] + y[8]/sqrt(3));
                  myres_rb11.push_back(y[9] + y[12] + y[17]/sqrt(3));
                  myres_H.push_back(H);
            }


            /* progress bar */
            if (flush == "yes"){
                  cout << " [";
                  progress = ((double)i)/((double)length_vec_ode_sol) + 1./(double)length_vec_ode_sol; 
                  pos = barWidth * progress - 1.;
                  for (int k = 0; k < barWidth; ++k) {
                        if (k < int(pos) ) cout << "=";
                        else if (k == int(pos)) cout << ">";
                        else cout << " ";
                  }
                  cout << "] " << int(progress * 100.0) << " %" << " @ " << "T = " << _Tref_/ti << "        " << "\r";
                  cout.flush();
            }
            /* -- print normalization of znu, zg -- */
            if (i == int(length_vec_ode_sol/500)){
                  cout << std::fixed << std::setprecision(6) << "numerical deviation from initial znu, zg equilibrium  : " << 1 - znu/znu_ini << ", " << 1 - zg/zg_ini << "               \n";
                  cout.precision(6);
            }
      }

      /* last time step stored in vectors for the internal Xn evolution */
      myres_x.push_back(t);
      myres_Tnu.push_back(_Tref_/ti * y[18]);
      myres_Tg.push_back(_Tref_/ti * y[19]);
      myres_r11.push_back(y[0] + y[3] + y[8]/sqrt(3));
      myres_rb11.push_back(y[9] + y[12] + y[17]/sqrt(3));
      myres_H.push_back(pow(_Tref_/ti,2.) * 1./(_MPL_/sqrt(11./2. * pow(y[18],4.) + 7./8. * 3. * (y[0] + y[9]) * pow(y[19],4.))));



      /* free the GSL driver when intergaration finished */
      gsl_odeiv2_driver_free (driver);

      /* if code did not timed out then everything OK; else exit */
      if (timeoutT == 0){
            cout << "\033[32m\n...integration finished\033[0m " << endl;
      }
      else{
            return {-99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99};
      }


      cout << "now evolve Xn ... " << endl;

      Interpolator1D interp_Tnu(myres_x, myres_Tnu);
      Interpolator1D interp_Tg(myres_x, myres_Tg);
      Interpolator1D interp_r11(myres_x, myres_r11);
      Interpolator1D interp_rb11(myres_x, myres_rb11);
      Interpolator1D interp_H(myres_x, myres_H);

      Xn_ODE_Params params_Xn = {&interp_Tnu, &interp_Tg, &interp_r11, &interp_rb11, &interp_H};
      double Xn_ini     = 1./(1. + 1./(pow(_Mn_/_Mp_,3./2.)*exp(-_Q_/_Tref_)*exp(-xi1)));
      double y_Xn[1] = {Xn_ini}; 

      gsl_odeiv2_system sys_Xn = {ode_Xn, nullptr, 1, &params_Xn};
  
      gsl_odeiv2_driver* driver_Xn = gsl_odeiv2_driver_alloc_y_new(&sys_Xn, gsl_odeiv2_step_msadams, 1e-10, 1e-10, 1e-10);

      double t0_Xn = _Tref_/T_ini;
      t = t0_Xn;
      double length_vec_ode_Xn_sol = 5.e3;

      vector<double> myres_Xn;
      myres_Xn.push_back(t);
      myres_Xn.push_back(Xn_ini);

      progress = 0.0;
      for (int i = 1; i < length_vec_ode_Xn_sol; i++){
            /* evolve system one timestep */
            ti = t0_Xn + i * (t_ave - t0_Xn) / (length_vec_ode_Xn_sol-1);
            int status = gsl_odeiv2_driver_apply (driver_Xn, &t, ti, y_Xn);
            if (status != GSL_SUCCESS) {
                  std::cerr << "ODE solver error at t = " << t << "\n";
                  exit(0);
                  break;
            }
            myres_Xn.push_back(ti);
            myres_Xn.push_back(y_Xn[0]);
            /* progress bar */
            if (flush == "yes"){
                  cout << " [";
                  progress = ((double)i)/((double)length_vec_ode_Xn_sol) + 1./(double)length_vec_ode_Xn_sol; 
                  pos = barWidth * progress - 1.;
                  for (int k = 0; k < barWidth; ++k) {
                        if (k < int(pos) ) cout << "=";
                        else if (k == int(pos)) cout << ">";
                        else cout << " ";
                  }
                  cout << "] " << int(progress * 100.0) << " %" << " @ " << "T = " << _Tref_/ti << "        " << "\r";
                  cout.flush();
            }
      }

      double Xn_final   = y_Xn[0];
      gsl_odeiv2_driver_free(driver_Xn);
      cout << "\033[32m\n...integration finished\033[0m " << endl;
      cout.precision(8);

      
      /* get final values of znu and zg at neutrino decoupling */
      double znudec     = y[18];
      double zgdec      = y[19];



      /* average oscillations */
      cout << "analytically calculate osc. average solution ...";
      double coeff_r[9]       = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8]};
      double coeff_rbar[9]    = {y[9], y[10], y[11], y[12], y[13], y[14], y[15], y[16], y[17]};
      get_mat_from_GellMann_coeff(coeff_r, lvec, r_fb);
      get_mat_from_GellMann_coeff(coeff_rbar, lvec, rbar_fb);
      vector<double> params_ave     = { theta12, theta13, theta23 };
      void *ptr_params_ave          = &params_ave;
      dn_alpha_res   = get_dn_alpha_osc_ave_from_r_rbar_fb(n_alpha_factor * pow(znudec,3.), r_fb, rbar_fb, ptr_params_ave);
      np_alpha_res   = get_np_alpha_osc_ave_from_r_rbar_fb(n_alpha_factor * pow(znudec,3.), r_fb, rbar_fb, ptr_params_ave);
      if (filename != ""){
            myres.push_back( _Tref_/t_ave );
            myres.push_back( dn_alpha_res[0] ); myres.push_back( dn_alpha_res[1] ); myres.push_back( dn_alpha_res[2] );
            myres.push_back( znudec ); 
            myres.push_back( zgdec );
            myres.push_back( np_alpha_res[0] ); myres.push_back( np_alpha_res[1] ); myres.push_back( np_alpha_res[2] );
      }
      cout << " finished " << endl;



      /* --------- print final results ---------- */
      array<double, 6> rho_res;
      rho_res     = get_rho_alpha_from_r_rbar_fb(znudec, r_fb, rbar_fb);

      cout << "\nX_n(T=1.5 MeV)     : " << Xn_final << endl;
      cout << "\u03BE_ini              : " << xi1 << " " << xi2 << " " << xi3 << endl;
      cout << "(n-n\u0305)/T^3 at BBN   : " << dn_alpha_res[0] <<  " " << dn_alpha_res[1] << " " << dn_alpha_res[2] << endl;
      cout << "z_\u03BD at \u03BD decoupling: " << znudec << endl; 
      cout << "z_\u03B3 at \u03BD decoupling: " << zgdec << endl; 
      cout << "\u03C1_\u03B1 at BBN         : " << rho_res[0] << " " << rho_res[1] << " " << rho_res[2] << endl;
      cout << "\u03C1\u0305_\u03B1 at BBN         : " << rho_res[3] << " " << rho_res[4] << " " << rho_res[5] << endl;
      /* ------- end print ------- */



      /* evolve zgamma after neutrino decoupling from Tave to 0.006 MeV (something smaller than _Me_) */
      const int dim_zg = 1;
      gsl_odeiv2_system sys_zg = {ode_zg, nullptr, dim_zg, nullptr};
      gsl_odeiv2_driver * d_zg = gsl_odeiv2_driver_alloc_y_new (&sys_zg, gsl_odeiv2_step_msadams, 1e-10, 1e-10, 1e-10);
      
      double tini, tend;
      tini  = _Tref_/T_ave;
      tend  = _Tref_/0.006; // some number smaller than _Me_ = 0.511
      t = tini;

      double y_zg[1] = {zgdec};

      if (filename != ""){
            int num=1e3;
            for (int i = 1; i < num; i++){
                  ti = tini + i * (tend - tini) / (num-1);
                  int status_zg = gsl_odeiv2_driver_apply (d_zg, &t, ti, y_zg);
                  if (status_zg != GSL_SUCCESS) {
                        std::cerr << "Error: " << gsl_strerror(status_zg) << std::endl;
                        break;
                  }
                  
                  myres.push_back( _Tref_/ti );
                  myres.push_back( dn_alpha_res[0] );
                  myres.push_back( dn_alpha_res[1] );
                  myres.push_back( dn_alpha_res[2] );
                  myres.push_back( znudec );
                  myres.push_back( y_zg[0] );
                  myres.push_back( np_alpha_res[0] ); 
                  myres.push_back( np_alpha_res[1] ); 
                  myres.push_back( np_alpha_res[2] );
            }
      }
      else{
            int status_zg = gsl_odeiv2_driver_apply (d_zg, &t, tend, y_zg);
            if (status_zg != GSL_SUCCESS) {
                  std::cerr << "Error: " << gsl_strerror(status_zg) << std::endl;
            }    
      }
      gsl_odeiv2_driver_free(d_zg);
      double zg_final   = y_zg[0];
      cout << "z_\u03B3 final          : " << zg_final << endl;

      /* 
      get Neff after e+e- annihilation has finished 
      and normalize our value for \xi_\alpha = 0 to SM accepted value of 3.044
      */
      double NeffNorm, Neff;
      NeffNorm    = 3.044 / 3.00183060; 
      Neff  = NeffNorm * pow(11./4.,4./3.) * 1./2. * 3. * (y[0] + y[9]) * pow(znudec/zg_final,4.);
      cout << "Final Neff         : " << Neff << endl;


      /* get final YP value */
      double yp;
      yp = YP(dn_alpha_res[0], rho_res[0], rho_res[3], zg_final, Neff);
      cout << "YP                 : " << yp << endl; 

      /* get final DH value */
      double dh;
      dh = DH(dn_alpha_res[0], rho_res[0], rho_res[3], zg_final, Neff);
      cout << "D/H                : " << dh << endl; 

      if (filename != ""){
            /* ------------ safe to file ---------- */
            const int dim_vector= 9;
            const int row_size   = int(myres.size()/dim_vector);
            double myres_matrix [row_size][dim_vector];
            for (int i = 0; i < row_size; ++i) {
                  for (int j=0; j < dim_vector; j++){
                        myres_matrix[i][j] =  myres[dim_vector*i+j];
                  }
            }

            ofstream file_asym;
            string file_name_asym;         
            file_name_asym         = _Path_ + "/output/asym_evolve_" + filename + ".dat";
            file_asym.open(file_name_asym);
            file_asym << "# " << "T" << "\t" << "nnbar11" << "\t" << "nnbar22" << "\t" << "nnbar33" << "\t" << 
                  "znu" << "\t" << "zg" 
                  << "\t" << "np11" << "\t" << "np22" << "\t" << "np33"
                  << endl;
            for (int i = 0; i < row_size ; i++){
                  double t          = myres_matrix[i][0]; 
                  double nnbar11    = myres_matrix[i][1];
                  double nnbar22    = myres_matrix[i][2]; 
                  double nnbar33    = myres_matrix[i][3]; 
                  double znu        = myres_matrix[i][4];
                  double zg         = myres_matrix[i][5];
                  double np11       = myres_matrix[i][6];
                  double np22       = myres_matrix[i][7];
                  double np33       = myres_matrix[i][8];

                  file_asym << std::fixed << std::setprecision(16) << std::defaultfloat << 
                        t << "\t"<< nnbar11 << "\t"<< nnbar22 << "\t"<< nnbar33 << "\t"<< 
                        znu << "\t"<< zg 
                        << "\t" << np11 << "\t" << np22 << "\t" << np33 
                        << endl;
            }
            file_asym.close();
            /* ------------ end safe to file ---------- */

            /* ----------- safe Xn evolution ---------- */
            int dim_vector_Xn= 2;
            int row_size_Xn   = int(myres_Xn.size()/dim_vector_Xn);
            double myres_matrix_Xn [row_size_Xn][dim_vector_Xn];
            for (int i = 0; i < row_size_Xn; ++i) {
                  for (int j=0; j < dim_vector_Xn; j++){
                        myres_matrix_Xn[i][j] =  myres_Xn[dim_vector_Xn*i+j];
                  }
            }

            ofstream file_Xn;
            string file_name_Xn;         
            file_name_Xn         = _Path_ + "/output/Xn_evolve_" + filename + ".dat";
            file_Xn.open(file_name_Xn);
            file_Xn << "# " << "x" << "\t" << "Xn" 
                  << endl;
            for (int i = 0; i < row_size_Xn ; i++){
                  double t    = myres_matrix_Xn[i][0]; 
                  double Xn   = myres_matrix_Xn[i][1];

                  file_Xn << std::fixed << std::setprecision(16) << std::defaultfloat << 
                        t << "\t"<< Xn 
                        << endl;
            }
            file_Xn.close();

            /* ----------- safe final result ---------- */
            ofstream file_coflasy;
            string file_name_coflasy;         
            file_name_coflasy = _Path_ + "/output/coflasy_" + filename + ".dat";
            file_coflasy.open(file_name_coflasy);

            file_coflasy << "# " << "1: xi1_ini" << "\t" << "2: xi2_ini"  << "\t" << "3: xi3_ini" 
                        << "\t" << "4: dn1_final" << "\t" << "5: dn2_final"  << "\t" << "6: dn3_final" 
                        << "\t" << "7: znu_final" << "\t" << "8: zg_final"  << "\t" << "9: Neff_final" 
                        << "\t" << "10: rho_e" << "\t" << "11: rho_mu"  << "\t" << "12: rho_tau" 
                        << "\t" << "13: rhobar_e" << "\t" << "14: rhobar_mu"  << "\t" << "15: rhobar_tau" 
                        << "\t" << "16: Xn_dec" << "\t" << "17: Yp"  << "\t" << "18: D/H" 
                        << endl;

            file_coflasy << std::fixed << std::setprecision(8) << std::defaultfloat 
                        << xi1 << "\t" << xi2  << "\t" << xi3
                        << "\t" << dn_alpha_res[0] << "\t" << dn_alpha_res[1]  << "\t" << dn_alpha_res[2]
                        << "\t" << znudec << "\t" << zg_final  << "\t" << Neff
                        << "\t" << rho_res[0] << "\t" << rho_res[1]  << "\t" <<  rho_res[2]
                        << "\t" << rho_res[3] << "\t" << rho_res[4] << "\t" <<  rho_res[5]
                        << "\t" << Xn_final << "\t" << yp  << "\t" << dh
                        << endl;
            file_coflasy.close();
      }


      
      return {xi1, xi2, xi3, dn_alpha_res[0], dn_alpha_res[1], dn_alpha_res[2], znudec, zg_final, Neff, rho_res[0], rho_res[1], rho_res[2], rho_res[3], rho_res[4], rho_res[5], Xn_final, yp, dh};


}








// main
int main(int argc, char* argv[]){
      (void)argc; (void)argv;
      auto start_main = chrono::system_clock::now();
      cout.precision(6);

      /*
      set the electron density and pressure functions
      */
      set_density();

      /*
      set the Collision interpolation function
      */
      set_collision1D();
      set_collision2D();
      set_collision2D_me();

      /*
      set the Transfer rate interpolation function
      */
      set_transfer1D();
      set_transfer2D();

      /*
      set the weak rates for neutron to proton fraction
      */
      set_weakrates1D();
      set_weakrates2D();
      set_weakratescorrection();


      if (argv[1] == NULL){
            cout << "\033[31mYou need to pass at least one argument.\033[39;49m" << endl;
            cout << "exemplary usage: ./coflasy-c.exe ../example.ini" << endl;
            cout << endl;
      }
      else{
            string arg(argv[1]);
            bool found_ini = arg.find(".ini") != string::npos;
            if (found_ini == true){
                  vector<double> vars     = read_ini(argv[1]);
                  double xi1,xi2,xi3;
                  vector<string> filename = read_ini_filename(argv[1]);
                  /* check that filename[0] (filename) is not empty and if yes set it to "results" */
                  if (filename[0] == ""){
                        filename[0] = "results";
                  }
                  cout << "results will be stored in: " << filename[0] << ".dat" << endl;
                  /* check that filename[1] (flush) is either yes or not and if not set it to yes */
                  if (filename[1] != "yes" && filename[1] != "no"){
                        filename[1] = "yes";
                  }
                  cout << "progressbar will be shown: " << filename[1] << endl;
                  if (vars[5] == 0){
                        cout << "passed xi_alpha" << endl;
                        xi1= vars[6]; xi2 = vars[7]; xi3 = vars[8];
                  }
                  else if( vars[5] == 1){
                        cout << "passed tilde{xi}_alpha" << endl;
                        xi1= paramtranslate.get_xi_alpha_from_xi_tilde_alpha(vars[9]);
                        xi2= paramtranslate.get_xi_alpha_from_xi_tilde_alpha(vars[10]);
                        xi3= paramtranslate.get_xi_alpha_from_xi_tilde_alpha(vars[11]);
                  }
                  else if( vars[5] == 2){
                        cout << "passed \u0394n_alpha" << endl;
                        xi1= paramtranslate.get_xi_alpha_from_dn_alpha(vars[12]);
                        xi2= paramtranslate.get_xi_alpha_from_dn_alpha(vars[13]);
                        xi3= paramtranslate.get_xi_alpha_from_dn_alpha(vars[14]);
                  }
                  else{
                        cout << "\033[31mCould not read the intial condition for the chemical potential.\033[39;49m" << endl;
                        return _SUCCESS_;
                  }


                  vector<double> params   = {vars[0], vars[1], vars[2], vars[3], vars[4], xi1, xi2, xi3, vars[15], vars[16], vars[17]};
                  void *ptr_params        = &params;
                  ode_evolve(ptr_params,filename[0],filename[1]);
            }
            else{
                  cout << "\033[31mPassed non-readable argument.\033[39;49m" << endl;
                  cout << "Allowed argument is:" << endl;
                  cout << "file with an .ini extension" << endl;
                  cout << "see example.ini for an example of an .ini file" << endl;
                  cout << endl;
            }
      }


      /***********************************************************
       Done with everything -- do not put anything below this line
      ***********************************************************/
      auto end_main = chrono::system_clock::now();
      chrono::duration<double> elapsed_seconds_main = end_main-start_main;
      cout <<"\ntime: " << scientific << setprecision(4) << elapsed_seconds_main.count() << " seconds" << endl;

      return _SUCCESS_;
}
