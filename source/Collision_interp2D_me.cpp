/*
Set the collision term as function of Tg and electron mass (zg, x)
*/

class Collision2D_me{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_zg;
    vector<double> myval_x;
    // nu nu <-> e e
    vector<double> myval_A34_nunu_ee_I0_me;




    vector<double> vec_zg;
    vector<double> vec_x;
    vector<double> vec_A34_nunu_ee_I0_me;


    const double *double_zg;
    const double *double_x;
    const double *double_A34_nunu_ee_I0_me;


    int num_zg;
    int num_x;
    void *za_alloc;
    double* za;


    gsl_spline2d *spline_A34_nunu_ee_I0_me;



    const gsl_interp2d_type *InterType;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;





public:
    /*
    Constructor of the class; read files of the Transfer Rates and set the shared interpolation settings
    */
    Collision2D_me(void){
        ifstream myfile_zg;
        myfile_zg.open(_Path_ + "/data/Collision/z.dat");
        if ( myfile_zg.is_open() ) {
            double tmp = 0.0;
            while (myfile_zg >> tmp) {
                myval_zg.push_back(tmp);
            }
        }
        myfile_zg.close();
        for (int i = 0; i < int(myval_zg.size()); i++){
            vec_zg.push_back(myval_zg[i]);
        }

        ifstream myfile_x;
        myfile_x.open(_Path_ + "/data/Collision/T_MeV.dat");
        if ( myfile_x.is_open() ) {
            double tmp = 0.0;
            while (myfile_x >> tmp) {
                myval_x.push_back(tmp);
            }
        }
        myfile_x.close();
        for (int i = 0; i < int(myval_x.size()); i++){
            vec_x.push_back(myval_x[i]);
        }


        // A34_nunu_ee
        ifstream myfile_A34_nunu_ee_I0_me;
        myfile_A34_nunu_ee_I0_me.open(_Path_ + "/data/Collision/A34_nunu_ee_I0_me.dat");
        if ( myfile_A34_nunu_ee_I0_me.is_open() ) {
            double tmp = 0.0;
            while (myfile_A34_nunu_ee_I0_me >> tmp) {
                myval_A34_nunu_ee_I0_me.push_back(tmp);
            }
        }
        myfile_A34_nunu_ee_I0_me.close();
        for (int i = 0; i < int(myval_A34_nunu_ee_I0_me.size()); i++){
            vec_A34_nunu_ee_I0_me.push_back(myval_A34_nunu_ee_I0_me[i]);
        }



        double_zg      = &vec_zg[0];
        double_x       = &vec_x[0];
        num_zg         = vec_zg.size()*sizeof(double_zg) / sizeof(double); /* x grid points */
        num_x          = vec_x.size()*sizeof(double_x) / sizeof(double); /* y grid points */
        za_alloc        = malloc(num_zg * num_x * sizeof(double));
        za              = (double *)za_alloc;

        InterType       = gsl_interp2d_bilinear;
        xacc            = gsl_interp_accel_alloc();
        yacc            = gsl_interp_accel_alloc();


        spline_A34_nunu_ee_I0_me      = gsl_spline2d_alloc(InterType, num_zg, num_x);

    }




    /**************************************************
    Set interpolation grid
    ***************************************************/


    /*
    set the interpolation grid for the A34_nunu_ee_I0_me rate
    */
    void set_A34_nunu_ee_I0_me(void){
        double_A34_nunu_ee_I0_me    = &vec_A34_nunu_ee_I0_me[0];
        for (int i = 0; i < int(myval_zg.size()); i++){
            for (int j = 0; j < int(myval_x.size()); j++){
                gsl_spline2d_set(spline_A34_nunu_ee_I0_me, za, i, j, double_A34_nunu_ee_I0_me[int(myval_zg.size())*i+j]);
            }
        }
        /* initialize interpolation */
        gsl_spline2d_init(spline_A34_nunu_ee_I0_me, double_zg, double_x, za, num_zg, num_x);
    }




    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/


    /*
    Function returning the interpolated A34_nunu_ee_I0_me rate
    */
    double A34_nunu_ee_I0_me(double zg, double x){
        return gsl_spline2d_eval(spline_A34_nunu_ee_I0_me, zg , x, xacc, yacc);
    }

    
};


/*
set rates globally, accesible by various modules
*/
Collision2D_me collision2D_me;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_collision2D_me(void){
    collision2D_me.set_A34_nunu_ee_I0_me();
}