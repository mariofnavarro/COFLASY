/*
Class to translate between different parametrization of initial conditions.
1) the chemical potential \xi = \mu/T
2) the reduced chemical potential \tilde{\xi} = \xi + \xi^3/\pi^2 used in e.g. 2405.06509
it is related to dn = n-\bar{n}/T^3 = 1/6 \tilde{\xi}
*/

class ParamTranslate{
private:
    double xi_tilde_alpha_to_dn_alpha_factor = 1./6.;
    double dn_alpha_to_xi_tilde_alpha_factor = 1./xi_tilde_alpha_to_dn_alpha_factor;
public:
    /*
    function translating from \xi_\alpha to \tilde{\xi}_\alpha
    */
    double get_xi_tilde_alpha_from_xi_alpha(double xi){
        return xi + pow(xi,3.)/pow(_PI_,2.);
    }

    /*
    function translating from \tilde{\xi}_\alpha to L_\alpha
    */
    double get_dn_alpha_from_xi_tilde_alpha(double xi_tilde){
        return this->xi_tilde_alpha_to_dn_alpha_factor * xi_tilde;
    }

    /*
    function translating from \tilde{\xi}_\alpha to \xi_\alpha
    */
    double get_xi_alpha_from_xi_tilde_alpha(double xitilde){

        double nom        = pow(_PI_/6.,2./3.) * ( -2.*pow(3.,1./3.)*pow(_PI_,2./3.) + pow(2.,1./3.)*pow(9.*xitilde + sqrt(12.*pow(_PI_,2.)+81.*pow(xitilde,2.)),2./3.) );
        double denom      = pow(9.*xitilde + sqrt(12.*pow(_PI_,2.)+81.*pow(xitilde,2.)),1./3.);

        return nom/denom;
    }

    /*
    function translating from dn_\alpha to \tilde{\xi}_\alpha
    */
    double get_xi_tilde_alpha_from_dn_alpha(double dn){
        return this->dn_alpha_to_xi_tilde_alpha_factor * dn;
    }



    /*
    function translating from \xi_\alpha to dn_\alpha
    */
    double get_dn_alpha_from_xi_alpha(double xi){
        return this->get_dn_alpha_from_xi_tilde_alpha(this->get_xi_tilde_alpha_from_xi_alpha(xi));
    }


    /*
    function translating from dn_\alpha to \xi_\alpha
    */
    double get_xi_alpha_from_dn_alpha(double dn){
        return this->get_xi_alpha_from_xi_tilde_alpha(this->get_xi_tilde_alpha_from_dn_alpha(dn));
    }




    /*
    function returning the total \Delta n as a function of \xi_\alpha
    */
    double get_dn_from_xi_alpha(double xi1, double xi2, double xi3){
        double dn1, dn2, dn3;
        dn1    = this->get_dn_alpha_from_xi_alpha(xi1);
        dn2    = this->get_dn_alpha_from_xi_alpha(xi2);
        dn3    = this->get_dn_alpha_from_xi_alpha(xi3);  

        return dn1+dn2+dn3;
    }


    /*
    function returning \sum_\alpha \tilde{\xi}_\alpha from total \Delta n input
    */
    double get_xitilde_tot_from_dn_tot(double dn){
        return this->dn_alpha_to_xi_tilde_alpha_factor*dn;
    }



    /*
    calculate the contribution to \Delta N_{eff} for given chmical potentials \xi_\alpha
    */
    double get_delta_Neff_from_xi_alpha(double xi1, double xi2, double xi3){
        return 30./7./pow(_PI_,2.)*(pow(xi1,2.) + pow(xi2,2.) + pow(xi3,2.)) + 15./7./pow(_PI_,4.)*(pow(xi1,4.)+pow(xi2,4.)+pow(xi3,4.));
    }


};

/*
initialize class globally
*/
ParamTranslate paramtranslate;


