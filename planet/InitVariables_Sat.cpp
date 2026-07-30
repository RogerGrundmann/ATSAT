/*
 * Atmosphere General Circulation Modell(ATSAT) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in aa spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and nh3 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
 * 
 * class to prepare the boundary and initial conditions for diverse variables
*/
#include "cSaturnModel.h"
#include "Utils.h"

using namespace std;
using namespace AtomUtils;

void cSaturnModel::init_temperature(){
    cout << endl << "      ATSAT: init_Temperature" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    int j_half = (jm-1)/2;
    double d_j_half = (double)j_half;
    double t_h2_eff = t_pole - t_equator; // non-dimensional
    double d_j = 0.0;
    double height = 0.0;

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            d_j = (double)(j);
            t.x[0][j][k] = t_h2_eff * (d_j * d_j/(d_j_half * d_j_half) 
                - 2.0 * d_j/d_j_half) + t_pole; // non-dimensional    assumption of t_pole < t_equator
            t.x[0][j][k] = t.x[0][j][k]/t_ref;

            for(int i = 0; i < im; i++){
                height = get_layer_height(i);
                t.x[i][j][k] = - gam * height/t_ref + t.x[0][j][k];  // linear temperature decay up to tropopause

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height[km] = " << height << endl
        << "     t_equator[°C] = " << t_equator - t_ref 
        << "     t_pole[°C] = " << t_pole - t_ref << endl
        << "     gam[K/km] = " << gam
        << "     gam * height[K] = " << gam * height << endl
        << "     t_u[°C] = " << t.x[i][j][k] * t_ref - t_ref
        << "     t_u[K] = " << t.x[i][j][k] * t_ref << endl;
*/

            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_Temperature\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_Temperature ended" << endl;
    return;
}
/*
*
*/
void cSaturnModel::init_PressureDynamic(){
    cout << endl << "      ATSAT: init_PressureDynamic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                p_dyn.x[i][j][k] = r_mix/2.0 
                    * (pow((u.x[i][j][k] * u_0), 2.0) 
                    + pow((v.x[i][j][k] * u_0), 2.0) 
                    + pow((w.x[i][j][k] * u_0), 2.0))/3.0 * 1e-5; // in bar = 1e5 Pa


/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height = " << height << endl
        << "     gam = " << gam 
        << "     g = " << g 
        << "     R_ref[J/(g*K)] = " << R_ref 
        << "     p_ref[bar] = " << p_ref 
        << "     exp_pressure = " << exp_pressure << endl
        << "     t_ref[K] = " << t_ref
        << "     t_ref[°C] = " << t_ref - t_ref << endl
        << "     t_u[K] = " << t_u
        << "     t_u[°C] = " << t_u - t_ref
        << "     p_u[bar] = " << p_stat.x[i][j][k] << endl;
*/
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_PressureDynamic ended" << endl;
    return;
}
/*
*
*/
/*
void cSaturnModel::init_Density(){
    cout << endl << "      ATSAT: init_Density" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                rho.x[i][j][k] = 
                    1.0e5 * (p_dyn.x[i][j][k] + p_stat.x[i][j][k])  //  in kg/m³
                    /(R_mix * t.x[i][j][k] * t_ref);


    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     R_mix = " << R_mix << endl
        << "     t[K] = " << t.x[i][j][k] * t_ref - t_ref
        << "     p_dyn[bar] = " << p_dyn.x[i][j][k]
        << "     p_stat[bar] = " << p_stat.x[i][j][k]
        << "     rho[kg/m³] = " << rho.x[i][j][k] << endl;

            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_Density\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_Density ended" << endl;
    return;
}
*/
/*
*
*/
void cSaturnModel::init_PressureStatic(){
    cout << endl << "      ATSAT: init_PressureStatic" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    double exp_pressure = 0.0;
    double t_u = 0.0;

//    R_ref = r_mix;

    for(int k = 0; k < km; k++){
        for(int j = 0; j < jm; j++){
            for(int i = 0; i < im; i++){
                t_u = t.x[i][j][k] * t_ref;
                exp_pressure = g/(gam * R_ref);
                p_stat.x[i][j][k] = p_ref 
                    * pow((t_u/t_ref), exp_pressure); // in bar
/*
                p_stat.x[i][j][k] = p_ref 
                    * pow((1.0 - gam * (height - height_0)/t_ref), 
                    exp_pressure);
*/

//                if(is_land(SeaMount, i, j, k)) p_stat.x[i][j][k] = 0.0;


/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     height = " << height << endl
        << "     gam = " << gam 
        << "     g = " << g 
        << "     R_ref[J/(g*K)] = " << R_ref 
        << "     p_ref[bar] = " << p_ref 
        << "     exp_pressure = " << exp_pressure << endl
        << "     t_ref[K] = " << t_ref
        << "     t_ref[°C] = " << t_ref - t_ref << endl
        << "     t_u[K] = " << t_u
        << "     t_u[°C] = " << t_u - t_ref
        << "     p_u[bar] = " << p_stat.x[i][j][k] << endl;
*/
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_PressureStatic\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_PressureStatic ended" << endl;
    return;
}
/*
*
*/
void cSaturnModel::Forces(){
    cout << endl << "      ATSAT: Forces" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    double rm, sinthe, costhe, rmsinthe;

    for(int i = 1; i < im-1; i++){
        rm = rad.z[i];
        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            costhe = cos(the.z[j]);
            rmsinthe = rm * sinthe;
                for(int k = 1; k < km-1; k++){

// influence of the Coriolis force
                    double Coriolis_rad = - 2.0 * omega * sinthe * w.x[i][j][k];
                    double Coriolis_the = + 2.0 * omega * costhe * w.x[i][j][k];
                    double Coriolis_phi = + 2.0 * omega * (- costhe * v.x[i][j][k] 
                        + sinthe * u.x[i][j][k]);

// influence of the centrifugal force
                    double dpdr = (p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k])/(2.0 * dr);
                    double dpdthe = (p_dyn.x[i][j+1][k] - p_dyn.x[i][j-1][k])/(2.0 * dthe);
                    double dpdphi = (p_dyn.x[i][j][k+1] - p_dyn.x[i][j][k-1])/(2.0 * dphi);

// acting forces
                    CoriolisForce.x[i][j][k] = Coriolis * r_mix 
                        * sqrt((pow(Coriolis_rad, 2) 
                        + pow(Coriolis_the, 2) 
                        + pow(Coriolis_phi, 2))/3.0);

                    // Magnitude of the centrifugal acceleration, times the density.
                    // With a_r = Omega^2*r*sin^2 and a_theta = Omega^2*r*sin*cos (see the
                    // derivation in RHS_Sat_Turb.cpp), the magnitude is
                    //     |a| = Omega^2*r*sin(theta)*sqrt(sin^2 + cos^2) = Omega^2*r*sin(theta),
                    // which is just Omega^2 times the distance from the rotation axis.
                    //
                    // This line used to read Omega^2*r*(1 + |sin(theta)|), a THIRD formula —
                    // neither a component of the force nor its magnitude, and largest nowhere in
                    // particular: it stayed at full strength Omega^2*r at the pole, where the
                    // centrifugal force is zero. The diagnostic and the equation now agree.
                    CentrifugalForce.x[i][j][k] = centrifugal * r_mix
                        * omega * omega * rm * sinthe;

                    BuoyancyForce.x[i][j][k] = buoyancy 
//                        * r_mix * g * (1.0 - (t.x[i][j][k] - 1.0))  //  rho0 * g - rho0 * (t - t0)/t0 * g    for   del_rho << rho0
                         * r_mix * g * (p_stat.x[i][j][k] + p_dyn.x[i][j][k])  //  rho * g
                        /(r_mix * R_mix * t.x[i][j][k] * t_ref) * 1e5;  // in N/m³

                    PresGradForce.x[i][j][k] = 
                        - sqrt((pow(dpdr, 2) 
                        + pow(dpdthe/rm, 2) 
                        + pow(dpdphi/rmsinthe, 2))/3.0)/L_atm * 1.0e5;

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "     i = " << i 
        << "     CoriolisForce[kN/m²] = " << CoriolisForce.x[i][j][k] << endl
        << "     CentrifugalForce[N/m³] = " << CentrifugalForce.x[i][j][k] << endl
        << "     BuoyancyForce[N/m³] = " << BuoyancyForce.x[i][j][k] << endl
        << "     PresGradForce[kN/m²] = " << PresGradForce.x[i][j][k] << endl << endl;
*/
            }
        }
    }
/*
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for Forces\n", elapsed.count() * 1e-9);
*/
    cout << "      ATSAT: Forces ended" << endl;
    return;
}
/*
*
*/
/*
void cSaturnModel::init_vapour(std::string gas, 
    double &c_tropopause, double &coeff_A, double &coeff_B, 
    double &coeff_A_i, double &coeff_B_i, 
    double &t_0, double &t_00, 
    double &ep, double &r, double &m,
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet, double &X,
    Array &c, Array &cloud, Array &ice, Array &cloudiness){

    cout << endl << "      ATSAT: init_vapour of "  << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // initial and boundary conditions of water vapour on water and land surfaces
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula

    double t_u = 0.0;
    double p_u = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;


    double rr = 0.0;

    double mix_ratio = 0.0; // planetary Sciences          volume mixing ratio
    if(gas == "H2")  rr = r_h2;
    if(gas == "HE")  rr = r_he;


// water vapour distribution decreasing approaching tropopause
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
//            int i_mount = get_surface_layer(j, k);
            for(int i = 0; i < im; i++){

                t_u = t.x[i][j][k] * t_ref;
                p_u = p_stat.x[i][j][k];
                E_Rain = cSaturnModel::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);  // saturation water vapour pressure for the water phase at t > 0°C in bar
//                q_Rain = ep * E_Rain/(p_u - E_Rain);  // species vapour amount at saturation with species formation in kg/kg
                q_Rain = ep * E_Rain/p_u;  // species vapour amount at saturation with species formation in kg/kg

                c.x[i][j][k] = r_mix * q_Rain; 
                if(c.x[i][j][k] <= 0.0)  c.x[i][j][k] = 0.0;


    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "init_vapour of " << gas << endl
        << "     i = " << i 
        << "     height = " << get_layer_height(i) 
        << "     i_trop = " << i_trop << endl 
        << "     t_u[K] = " << t.x[i][j][k]
        << "     p_u[bar] = " << p_stat.x[i][j][k] << endl

        << "     E_Rain[bar] = " << E_Rain
        << "     q_Rain[kg/kg] = " << q_Rain << endl
        << "     ep[/] = " << ep << endl

        << "     c[g/kg] = " << c.x[i][j][k] * 1e3 << endl << endl;


                if(is_land(SeaMount, i, j, k)) c.x[i][j][k] = 0.0;
            } // end i
        }// end k
    }// end j


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_vapour\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_vapour of "  << gas << endl;

    return;
}
*/
/*
*
*/
void cSaturnModel::init_vapour_cloud_ice(std::string gas, 
    double &c_tropopause, double &coeff_A, double &coeff_B, 
    double &coeff_A_i, double &coeff_B_i, 
    double &t_0, double &t_00, 
    double &ep, double &r, double &m,
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet,
    Array &c, Array &cloud, Array &ice, Array &cloudiness){

    cout << endl << "      ATSAT: init_vapour_cloud_ice of "  << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // initial and boundary conditions of water vapour on water and land surfaces
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula

    double t_u = 0.0;
    double p_u = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double magnus = 0.0;
    double r_max_equator = 0.0;
    double r_max_pole = 0.0;
    double r_max_add_equator = 0.0;
    double r_max_add_pole = 0.0;
    double t_add_equator = 0.0;
    double t_add_pole = 0.0;

    if(gas == "CH4"){
        r_max_equator = r_ch4;
        r_max_pole = 0.7 * r_ch4;
        magnus = 2.0;
    }
    if(gas == "H2O"){
        r_max_equator = r_h2o;
        r_max_pole = 0.7 * r_h2o;
        magnus = 2.0;
    }
    if(gas == "NH3"){
        r_max_equator = r_nh3;
        r_max_pole = 0.7 * r_nh3;

        r_max_add_equator = r_nh3_add;
        r_max_add_pole = 0.7 * r_nh3_add;

        t_add_equator = t_0_nh3;
        t_add_pole = 1.01 * t_add_equator;

        magnus = 2.0;
    }

// water vapour distribution decreases approaching tropopause,  
//    at p_stat = 1.0 bar = 1013.0 hPa ... E_Rain = 0.0061 bar = 6,1 hPa (t_u = 0°C) ... E_Rain = 0.00564 bar = 56.4 hPa (t_u = -35°C)


    r_max = std::vector<double>(jm, r_max_pole); // radial location of r_max
    double r_max_eff = r_max_pole - r_max_equator;  // coefficient for the zonal parabolic r_max extention

    r_max_add = std::vector<double>(jm, r_max_add_pole); // radial location of r_max
    double r_max_add_eff = r_max_add_pole - r_max_add_equator;  // coefficient for the zonal parabolic r_max extention

    t_add = std::vector<double>(jm, t_add_pole); // radial location of r_max
    double t_add_eff = t_add_pole - t_add_equator;  // coefficient for the zonal parabolic r_max extention

    double d_j_half = (double)(jm-1)/2.0;

    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        r_max[j] = r_max_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_pole;
        r_max_add[j] = r_max_add_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_add_pole;
        t_add[j] = t_add_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + t_add_pole;
    }

    for(int j = 0; j < jm; j++){

        for(int k = 0; k < km; k++){
            for(int i = 0; i <= im-1; i++){
                t_u = t.x[i][j][k] * t_ref;
                p_u = p_stat.x[i][j][k];

//                E_Rain = 1e3 * cSaturnModel::Clausius_Clapeyron(t_u, coeff_A, coeff_B);  // saturation species vapour pressure for the liquid phase at t > 0°C in bar
                E_Rain = cSaturnModel::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);

//                q_Rain = ep * E_Rain/(p_u - E_Rain);  // species vapour amount at saturation with species formation in kg/kg
                q_Rain = ep * E_Rain/p_u;  // species vapour amount at saturation with species formation in kg/kg

                c.x[i][j][k] = magnus * r_mix * q_Rain;            // coefficient is arbitrary

                if(c.x[i][j][k] >= r_max[j])  c.x[i][j][k] = r_max[j];

                if((gas == "NH3")&&(p_u >= p_00_nh3))  c.x[i][j][k] = r_max[j];
                if((gas == "H2O")&&(p_u >= p_0_h2o))  c.x[i][j][k] = r_max[j];
                if((gas == "CH4")&&(p_u >= p_0_ch4))  c.x[i][j][k] = r_max[j];

            }

/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if((j == 90)&&(k == 180)) cout << endl
                    << gas << "  cloudwater and cloudice    °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k << endl
                    << "   height = " << get_layer_height(i) << endl 
                    << "   cloud_loc = " << cloud_loc[j]
                    << "   cloud_loc_equator = " << cloud_loc_equator 
                    << "   cloud_loc_pole = " << cloud_loc_pole << endl
                    << "   x_cloud = " << x_cloud
                    << "   alfa_s = " << alfa_s << endl
                    << "   Humility_critical[/] = " << cSaturnModel::Humility_critical(x_cloud, 1.0, 0.8)
                    << "   del_q_ls[kg/kg] = " << del_q_ls << endl
                    << "   t_u[K] = " << t_u
                    << "   t_u[°C] = " << t_u - t_ref
                    << "   p_u[bar] = " << p_u << endl
                    << "   t_0[K] = " << t_0
                    << "   t_00[K] = " << t_00 << endl
                    << "   h_T[/] = " << h_T << endl
                    << "   E_Rain[bar] = " << E_Rain
                    << "   q_Rain[kg/kg] = " << q_Rain << endl
                    << "   c[kg/kg] = " << c.x[i][j][k]
                    << "   cloud[kg/kg] = " << cloud.x[i][j][k]
                    << "   ice[kg/kg] = " << ice.x[i][j][k] << endl
                    << "   cloudiness[/] = " << cloudiness.x[i][j][k] << endl << endl;
*/
        }// end k
    }// end j


/*
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(SeaMount, i, j, k)){
                    c.x[i][j][k] = 0.0;
                }
            }
        }
    }
*/
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_vapour_cloud_ice\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_vapour_cloud_ice of "  << gas << endl;

    return;
}
/*
*
*/
void cSaturnModel::init_h2s(std::string gas, 
    double &c_tropopause, double &coeff_A, double &coeff_B, 
    double &coeff_A_i, double &coeff_B_i, 
    double &t_0, double &t_00, 
    double &ep, double &r, double &m,
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet, Array &c){

    cout << endl << "      ATSAT: init_vapour_c_i of "  << gas << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // initial and boundary conditions of water vapour on water and land surfaces
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula

    double t_u = 0.0;
    double p_u = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double corr = 4.0e-3; // correction value to match the density values of  Planetary sciences p. 90 2010
    double r_max_equator = r_h2s;
    double r_max_pole = 0.7 * r_h2s;
    double d_j_half = (double)(jm-1)/2.0;

    r_max = std::vector<double>(jm, r_max_pole); // radial location of r_max
    double r_max_eff = r_max_pole - r_max_equator;  // coefficient for the zonal parabolic r_max extention

    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        r_max[j] = r_max_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_pole;
    }

    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_trop = im-1;
            for(int i = 0; i <= i_trop; i++){
                t_u = t.x[i][j][k] * t_ref;
                p_u = p_stat.x[i][j][k];

//                E_Rain = 1e3 * cSaturnModel::Clausius_Clapeyron(t_u, coeff_A, coeff_B);  // saturation species vapour pressure for the liquid phase at t > 0°C in bar
                E_Rain = cSaturnModel::saturation_vapour_pressure
                    (t_u, C, L0, R, del_alf, del_bet);

//                q_Rain = ep * E_Rain/(p_u - E_Rain);  // species vapour amount at saturation with species formation in kg/kg
                q_Rain = ep * E_Rain/p_u;  // species vapour amount at saturation with species formation in kg/kg

                c.x[i][j][k] = corr * r_mix * q_Rain;

                if((c.x[i][j][k]) >= r_max[j])  c.x[i][j][k] = r_max[j];

/*
                cout.precision(8);
                cout.setf(ios::fixed);
                if((j == 90)&&(k == 180)) cout << endl
                    << gas << "  cloudwater and cloudice    °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << endl
                    << "   i = " << i << "   j = " << j << "   k = " << k << endl
                    << "   height = " << get_layer_height(i) << endl 
                    << "   r_mix[Kg/m³] = " << r_mix
                    << "   ep[/] = " << ep << endl
                    << "   t_u[K] = " << t_u
                    << "   t_u[°C] = " << t_u - t_ref
                    << "   p_u[bar] = " << p_u << endl
                    << "   t_0[K] = " << t_0
                    << "   t_00[K] = " << t_00 << endl
                    << "   p_0[bar] = " << p_0_h2s
                    << "   p_00[bar] = " << p_00_h2s << endl
                    << "   E_Rain[bar] = " << E_Rain
                    << "   q_Rain[g/m³] = " << 1e3 * q_Rain << endl
                    << "   c[g/m³] = " << 1e3 * c.x[i][j][k] << endl << endl;
*/
            } // end i cloud
        }// end k
    }// end j

/*
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(is_land(SeaMount, i, j, k)){
                    c.x[i][j][k] = 0.0;
                }
            }
        }
    }
*/
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_vapour_c_i\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_vapour_c_i of "  << gas << endl;

    return;
}
/*
*
*/
void cSaturnModel::init_nh4sh(std::string gas, double &c_tropopause,
    double &coeff_A, double &coeff_B, 
    double &ep, double &r, double &m, 
    double &C, double &L0, double &R, 
    double &del_alf, double &del_bet, Array &c){

    cout << endl << "      ATSAT: init_nh4sh of " << " begin" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    // initial and boundary conditions of water vapour on water and land surfaces
    // value 0.04 stands for the maximum value of 40 g/kg, g water vapour per kg dry air
    // water vapour contents computed by Clausius-Clapeyron-formula

    double t_u = 0.0;
    double p_u = 0.0;
    double E_Rain = 0.0;
    double q_Rain = 0.0;
    double corr = 6.0; // correction value to match the density values of  Planetary sciences p. 90 2010
    double r_max_equator = r_nh4sh;
    double r_max_pole = 0.7 * r_nh4sh;
    double d_j_half = (double)(jm-1)/2.0;

    r_max = std::vector<double>(jm, r_max_pole); // radial location of r_max
    double r_max_eff = r_max_pole - r_max_equator;  // coefficient for the zonal parabolic r_max extention

    for(int j = 0; j < jm; j++){
        double d_j = (double)j;
        r_max[j] = r_max_eff 
            * AtomUtils::parabola((double)d_j/(double)d_j_half) + r_max_pole;
    }

// water vapour distribution decreasing approaching tropopause
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                t_u = t.x[i][j][k] * t_ref;
                p_u = p_stat.x[i][j][k];  // in bar
//                if((t_u <= t_0_nh4sh)&&(t_u >= t_00_nh4sh)){
//                if((p_u <= p_00_nh4sh)&&(p_u >= p_0_nh4sh)){
//                    E_Rain = 1e3 * cSaturnModel::Clausius_Clapeyron(t_u, coeff_A, coeff_B);  // saturation water vapour pressure for the water phase at t > 0°C in bar
                    E_Rain = cSaturnModel::saturation_vapour_pressure
                        (t_u, C, L0, R, del_alf, del_bet);  // saturation water vapour pressure for the water phase at t > 0°C in bar

//                    q_Rain = ep * E_Rain/(p_u - E_Rain);  // water vapour amount at saturation with water formation in kg/kg = non-dimensional
                    q_Rain = ep * E_Rain/p_u;  // water vapour amount at saturation with water formation in kg/kg = non-dimensional

                    c.x[i][j][k] = corr * r_mix * q_Rain; //  non- dimensional

                    if(c.x[i][j][k] <= 0.0)  c.x[i][j][k] = 0.0;

                    if(c.x[i][j][k] >= r_max[j])  c.x[i][j][k] = r_max[j];

                    if(t_u > 1.05 * t_0_nh4sh)  c.x[i][j][k] = 0.0;
//                    if(c.x[i][j][k] > r_max[j])  c.x[i][j][k] = 0.0;
//                    if(p_u > p_0_nh4sh)  c.x[i][j][k] = 0.0;

                    if(is_land(SeaMount, i, j, k))  c.x[i][j][k] = 0.0;

/*
    cout.precision(10);
    cout.setf(ios::fixed);
    if((j==90)&&(k==180)) cout << endl 
        << "init_vapour of " << gas << endl
        << "     i = " << i
        << "     height = " << get_layer_height(i) << endl 
        << "     t_u[K] = " << t_u
        << "     p_u[bar] = " << p_u << endl

        << "     r_mix[Kg/m³] = " << r_mix
        << "     magnus[/] = " << magnus << endl

        << "     E_Rain[bar] = " << E_Rain
        << "     q_Rain[kg/kg] = " << q_Rain << endl

        << "     c[g/kg] = " << c.x[i][j][k] * 1e3 << endl << endl;
*/

            } // end i
        }// end k
    }// end j


    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for init_nh4sh\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: init_nh4sh of "  << gas << endl;

    return;
}
/*
*
*/
