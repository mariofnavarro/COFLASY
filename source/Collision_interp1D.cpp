/*
Set interpolation function for 1D Collision terms 
*/

class Collision1D{
private:
    /*
    Set variables used within the class
    */
    vector<double> myval_z;

    vector<double> myval_A12_nunu_ee_I0;
    vector<double> myval_A12_nuan_I0;
    vector<double> myval_A123_nuan_I0;


    /* --- */

    vector<double> vec_z;

    vector<double> vec_A12_nunu_ee_I0;
    vector<double> vec_A12_nuan_I0;
    vector<double> vec_A123_nuan_I0;
    



    /* --- */

    const double *double_z;

    const double *double_A12_nunu_ee_I0;
    const double *double_A12_nuan_I0;
    const double *double_A123_nuan_I0;



    /* --- */
    int num_z;


    /* --- */

    gsl_spline *spline_A12_nunu_ee_I0;
    gsl_spline *spline_A12_nuan_I0;
    gsl_spline *spline_A123_nuan_I0;



    /* --- */

    const gsl_interp_type *InterType;
    gsl_interp_accel *acc;

public:
    /*
    Constructor of the class; read files of the rates and set the shared interpolation settings
    */
    Collision1D(void){
        ifstream myfile_z;
        myfile_z.open(_Path_ + "/data/Collision/z.dat");
        if ( myfile_z.is_open() ) {
            double tmp = 0.0;
            while (myfile_z >> tmp) {
                myval_z.push_back(tmp);
            }
        }
        myfile_z.close();
        for (int i = 0; i < int(myval_z.size()); i++){
            vec_z.push_back(myval_z[i]);
        }


        /* --- */
        ifstream myfile_A12_nunu_ee_I0;
        myfile_A12_nunu_ee_I0.open(_Path_ + "/data/Collision/A12_nunu_ee_I0.dat");
        if ( myfile_A12_nunu_ee_I0.is_open() ) {
            double tmp = 0.0;
            while (myfile_A12_nunu_ee_I0 >> tmp) {
                myval_A12_nunu_ee_I0.push_back(tmp);
            }
        }
        myfile_A12_nunu_ee_I0.close();
        for (int i = 0; i < int(myval_A12_nunu_ee_I0.size()); i++){
            vec_A12_nunu_ee_I0.push_back(myval_A12_nunu_ee_I0[i]);
        }


        /* --- */
        ifstream myfile_A12_nuan_I0;
        myfile_A12_nuan_I0.open(_Path_ + "/data/Collision/A12_nuan_I0.dat");
        if ( myfile_A12_nuan_I0.is_open() ) {
            double tmp = 0.0;
            while (myfile_A12_nuan_I0 >> tmp) {
                myval_A12_nuan_I0.push_back(tmp);
            }
        }
        myfile_A12_nuan_I0.close();
        for (int i = 0; i < int(myval_A12_nuan_I0.size()); i++){
            vec_A12_nuan_I0.push_back(myval_A12_nuan_I0[i]);
        }


        /* --- */
        ifstream myfile_A123_nuan_I0;
        myfile_A123_nuan_I0.open(_Path_ + "/data/Collision/A123_nuan_I0.dat");
        if ( myfile_A123_nuan_I0.is_open() ) {
            double tmp = 0.0;
            while (myfile_A123_nuan_I0 >> tmp) {
                myval_A123_nuan_I0.push_back(tmp);
            }
        }
        myfile_A123_nuan_I0.close();
        for (int i = 0; i < int(myval_A123_nuan_I0.size()); i++){
            vec_A123_nuan_I0.push_back(myval_A123_nuan_I0[i]);
        }




        /* --- */
        double_z                = &vec_z[0];

        double_A12_nunu_ee_I0      = &vec_A12_nunu_ee_I0[0];
        double_A12_nuan_I0      = &vec_A12_nuan_I0[0];
        double_A123_nuan_I0      = &vec_A123_nuan_I0[0];


        /* --- */
        num_z                   = vec_z.size()*sizeof(double_z) / sizeof(double); 


        /* --- */
        InterType               = gsl_interp_linear;
        acc                     = gsl_interp_accel_alloc();

        /* --- */
        spline_A12_nunu_ee_I0      = gsl_spline_alloc(InterType, num_z);
        spline_A12_nuan_I0      = gsl_spline_alloc(InterType, num_z);
        spline_A123_nuan_I0      = gsl_spline_alloc(InterType, num_z);


    }


    /**************************************************
    initialize interpolation 
    ***************************************************/

    void set_A12_nunu_ee_I0(void){
        gsl_spline_init(spline_A12_nunu_ee_I0, double_z, double_A12_nunu_ee_I0, num_z);
    }

    void set_A12_nuan_I0(void){
        gsl_spline_init(spline_A12_nuan_I0, double_z, double_A12_nuan_I0, num_z);
    }

    void set_A123_nuan_I0(void){
        gsl_spline_init(spline_A123_nuan_I0, double_z, double_A123_nuan_I0, num_z);
    }



    /**************************************************
    Declare methods returning the interpolated value
    ***************************************************/

    double A12_nunu_ee_I0(double z){
        return gsl_spline_eval(spline_A12_nunu_ee_I0, z , acc);
    }

    double A12_nuan_I0(double z){
        return gsl_spline_eval(spline_A12_nuan_I0, z , acc);
    }

    double A123_nuan_I0(double z){
        return gsl_spline_eval(spline_A123_nuan_I0, z , acc);
    }



};


/*
set collision terms globally, accesible by various modules
*/
Collision1D collision1D;
/*
call this in the int main() (amiqs.cpp) to set the interpolation grid
*/
void set_collision1D(void){
    // nu e -> nu e
    collision1D.set_A12_nunu_ee_I0();
    // nu nu -> nu nu annihilation
    collision1D.set_A12_nuan_I0();
    collision1D.set_A123_nuan_I0();

}