/*
Set interpolation function for 1D electron density terms as a function of Tgamma = Tref*zg/x
*/

class Density{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_Tg;

    vector<double> myval_drhoe;
    vector<double> myval_rhoe;
    vector<double> myval_pe;
    vector<double> myval_rhomu;
    vector<double> myval_pmu;


    /* --- */

    vector<double> vec_Tg;

    vector<double> vec_drhoe;
    vector<double> vec_rhoe;
    vector<double> vec_pe;
    vector<double> vec_rhomu;
    vector<double> vec_pmu;
    



    /* --- */

    const double *double_Tg;

    const double *double_drhoe;
    const double *double_rhoe;
    const double *double_pe;
    const double *double_rhomu;
    const double *double_pmu;


    /* --- */
    int num_Tg;


    /* --- */

    gsl_spline *spline_drhoe;
    gsl_spline *spline_rhoe;
    gsl_spline *spline_pe;
    gsl_spline *spline_rhomu;
    gsl_spline *spline_pmu;



    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;


public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    Density(void){
        ifstream myfile_Tg;
        myfile_Tg.open(_Path_ + "/data/rho_p/Tg.dat");
        if ( myfile_Tg.is_open() ) {
            double tmp = 0.0;
            while (myfile_Tg >> tmp) {
                myval_Tg.push_back(tmp);
            }
        }
        myfile_Tg.close();
        for (int i = 0; i < int(myval_Tg.size()); i++){
            vec_Tg.push_back(myval_Tg[i]);
        }



        /* --- */
        ifstream myfile_drhoe;
        myfile_drhoe.open(_Path_ + "/data/rho_p/drhoe.dat");
        if ( myfile_drhoe.is_open() ) {
            double tmp = 0.0;
            while (myfile_drhoe >> tmp) {
                myval_drhoe.push_back(tmp);
            }
        }
        myfile_drhoe.close();
        for (int i = 0; i < int(myval_drhoe.size()); i++){
            vec_drhoe.push_back(myval_drhoe[i]);
        }

        /* --- */
        ifstream myfile_rhoe;
        myfile_rhoe.open(_Path_ + "/data/rho_p/rhoe.dat");
        if ( myfile_rhoe.is_open() ) {
            double tmp = 0.0;
            while (myfile_rhoe >> tmp) {
                myval_rhoe.push_back(tmp);
            }
        }
        myfile_rhoe.close();
        for (int i = 0; i < int(myval_rhoe.size()); i++){
            vec_rhoe.push_back(myval_rhoe[i]);
        }


        /* --- */
        ifstream myfile_pe;
        myfile_pe.open(_Path_ + "/data/rho_p/pe.dat");
        if ( myfile_pe.is_open() ) {
            double tmp = 0.0;
            while (myfile_pe >> tmp) {
                myval_pe.push_back(tmp);
            }
        }
        myfile_pe.close();
        for (int i = 0; i < int(myval_pe.size()); i++){
            vec_pe.push_back(myval_pe[i]);
        }


        /* --- */
        ifstream myfile_rhomu;
        myfile_rhomu.open(_Path_ + "/data/rho_p/rhomu.dat");
        if ( myfile_rhomu.is_open() ) {
            double tmp = 0.0;
            while (myfile_rhomu >> tmp) {
                myval_rhomu.push_back(tmp);
            }
        }
        myfile_rhomu.close();
        for (int i = 0; i < int(myval_rhomu.size()); i++){
            vec_rhomu.push_back(myval_rhomu[i]);
        }

        /* --- */
        ifstream myfile_pmu;
        myfile_pmu.open(_Path_ + "/data/rho_p/pmu.dat");
        if ( myfile_pmu.is_open() ) {
            double tmp = 0.0;
            while (myfile_pmu >> tmp) {
                myval_pmu.push_back(tmp);
            }
        }
        myfile_pmu.close();
        for (int i = 0; i < int(myval_pmu.size()); i++){
            vec_pmu.push_back(myval_pmu[i]);
        }


        /* --- */
        double_Tg       = &vec_Tg[0];

        double_drhoe    = &vec_drhoe[0];
        double_rhoe     = &vec_rhoe[0];
        double_pe       = &vec_pe[0];
        double_rhomu    = &vec_rhomu[0];
        double_pmu      = &vec_pmu[0];


        /* --- */
        num_Tg          = vec_Tg.size()*sizeof(double_Tg) / sizeof(double); 


        /* --- */
        InterType       = gsl_interp_linear;
        acc             = gsl_interp_accel_alloc();

        /* --- */
        spline_drhoe    = gsl_spline_alloc(InterType, num_Tg);
        spline_rhoe     = gsl_spline_alloc(InterType, num_Tg);
        spline_pe       = gsl_spline_alloc(InterType, num_Tg);
        spline_rhomu    = gsl_spline_alloc(InterType, num_Tg);
        spline_pmu      = gsl_spline_alloc(InterType, num_Tg);


    }


    /**************************************************
    initialize interpolation 
    ***************************************************/


    void set_drhoe(void){
        gsl_spline_init(spline_drhoe, double_Tg, double_drhoe, num_Tg);
    }

    void set_rhoe(void){
        gsl_spline_init(spline_rhoe, double_Tg, double_rhoe, num_Tg);
    }

    void set_pe(void){
        gsl_spline_init(spline_pe, double_Tg, double_pe, num_Tg);
    }

    void set_rhomu(void){
        gsl_spline_init(spline_rhomu, double_Tg, double_rhomu, num_Tg);
    }

    void set_pmu(void){
        gsl_spline_init(spline_pmu, double_Tg, double_pmu, num_Tg);
    }


    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/


    double drhoe(double Tg){
        return gsl_spline_eval(spline_drhoe, Tg , acc);
    }

    double rhoe(double Tg){
        return gsl_spline_eval(spline_rhoe, Tg , acc);
    }

    double pe(double Tg){
        return gsl_spline_eval(spline_pe, Tg , acc);
    }

    double rhomu(double Tg){
        return gsl_spline_eval(spline_rhomu, Tg , acc);
    }

    double pmu(double Tg){
        return gsl_spline_eval(spline_pmu, Tg , acc);
    }



};


/*
set collision terms globally, accesible by various modules
*/
Density density;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_density(void){
    density.set_drhoe();
    density.set_rhoe();
    density.set_pe();
    density.set_rhomu();
    density.set_pmu();


}