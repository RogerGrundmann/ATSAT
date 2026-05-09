#include "cSaturnModel.h"
#include "VelocityInitializerSat.h"

using namespace std;
using namespace AtomUtils;

void VelocityInitializerSat::compute() { m.SaturnCellStructure(); }
/*
*
*/
void cSaturnModel::SaturnCellStructure(){
    cout << endl << "      Saturn: SaturnCellStructure" << endl;

    // boundary condition for the velocity components in the circulation cells

    // initial velocity components in the northern and southern
    // Pole, Ferrel and Hadley cells

    // equator (at j=90 compares to 0° latitude)
    


// northen hemisphere

// u-component

    init_u(u, 0);  // pole (center)
    init_u(u, 3);

    init_u(u, 5);
    init_u(u, 10);

    init_u(u, 15);
    init_u(u, 16);
 
    init_u(u, 17);
    init_u(u, 21);

    init_u(u, 28);
    init_u(u, 30);

    init_u(u, 32);
    init_u(u, 36);

    init_u(u, 40);
    init_u(u, 42);

    init_u(u, 44);
    init_u(u, 50);

    init_u(u, 57);
    init_u(u, 58);

    init_u(u, 59);
    init_u(u, 60);

    init_u(u, 70);
    init_u(u, 80);

    init_u(u, 90); // equator (center)

// southen hemisphere

// u-component

    init_u(u, 100);
    init_u(u, 110);

    init_u(u, 120);
    init_u(u, 121);

    init_u(u, 122);
    init_u(u, 123);

    init_u(u, 130);
    init_u(u, 136);

    init_u(u, 138);
    init_u(u, 140);

    init_u(u, 144);
    init_u(u, 148);

    init_u(u, 150);
    init_u(u, 152);

    init_u(u, 159);
    init_u(u, 163);

    init_u(u, 164);
    init_u(u, 165);

    init_u(u, 170);
    init_u(u, 175);

    init_u(u, 177);
    init_u(u, 180);  // pole (center)



    double v_max = 40.0;

// equator (center) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    init_v_or_w(v, 90, 0, 0);  // equator (center)
    init_v_or_w(w, 90, 200.0, 0.0);  // equator (center)


// 1. Hadley-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Hadley cells v
    init_v_or_w(v, 80, - v_max, v_max);  // max value
    init_v_or_w(v, 100, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 80, 280.0, 0.0);  // max value
    init_v_or_w(w, 100, 280.0, 0.0);  // max value



// 1. Ferrel-cell &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Ferrel cells v
    init_v_or_w(v, 70, 0, 0);  // center
    init_v_or_w(v, 110, 0, 0);  // center

    init_v_or_w(v, 60, v_max, - v_max);  // min value
    init_v_or_w(v, 120, - v_max, v_max);  // min value

// Ferrel cells w
    init_v_or_w(w, 70, 200, 0);  // center
    init_v_or_w(w, 110, 200, 0);  // center

    init_v_or_w(w, 60, 60.0, 0.0);  // min value
    init_v_or_w(w, 120, 60.0, 0.0);  // min value



// 2. Hadley-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Hadley cells v
    init_v_or_w(v, 59, 0.0, 0.0);  // center
    init_v_or_w(v, 121, 0.0, 0.0);  // center

    init_v_or_w(v, 58, - v_max, v_max);  // max value
    init_v_or_w(v, 122, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 59, 80.0, 0.0);  // center
    init_v_or_w(w, 121, 80.0, 0.0);  // center

    init_v_or_w(w, 58, 90.0, 0.0);  // max value
    init_v_or_w(w, 122, 90.0, 0.0);  // max value



// 2. Ferrel-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Hadley cells v
    init_v_or_w(v, 57, 0.0, 0.0);  // center
    init_v_or_w(v, 123, 0.0, 0.0);  // center

    init_v_or_w(v, 50, v_max, - v_max);  // min value
    init_v_or_w(v, 130, - v_max, v_max);  // min value

// Ferrel cells w
    init_v_or_w(w, 57, 80.0, 0.0);  // center
    init_v_or_w(w, 123, 80.0, 0.0);  // center

    init_v_or_w(w, 50, -20.0, 0.0);  // min value
    init_v_or_w(w, 130, -20.0, 0.0);  // min value



// 3. Hadley-cell &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Hadley cells v
    init_v_or_w(v, 44, 0.0, 0.0);  // center
    init_v_or_w(v, 136, 0.0, 0.0);  // center

    init_v_or_w(v, 42, - v_max, v_max);  // max value
    init_v_or_w(v, 138, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 44, 60.0, 0.0);  // center
    init_v_or_w(w, 136, 60.0, 0.0); // center

    init_v_or_w(w, 42, 140.0, 0.0);  // max value
    init_v_or_w(w, 138, 140.0, 0.0);  // max value



// 3. Ferrel-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Ferrel cells v
    init_v_or_w(v, 40, 0.0, 0.0); //  center
    init_v_or_w(v, 140, 0.0, 0.0); //  center

    init_v_or_w(v, 36, v_max, - v_max);  // min value
    init_v_or_w(v, 144, - v_max, v_max);  // min value

// Ferrel cells w
    init_v_or_w(w, 40, 80.0, 0.0); //  center
    init_v_or_w(w, 140, 80.0, 0.0); //  center

    init_v_or_w(w, 36, 20.0, 0.0);  // min value
    init_v_or_w(w, 144, 20.0, 0.0);  // min value



// 4. Hadley-cell &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Hadley cells v
    init_v_or_w(v, 32, 0.0, 0.0); //  center
    init_v_or_w(v, 148, 0.0, 0.0); //  center

    init_v_or_w(v, 30, - v_max, v_max);  // max value
    init_v_or_w(v, 150, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 32, 80.0, 0.0); //  center
    init_v_or_w(w, 148, 80.0, 0.0); //  center

    init_v_or_w(w, 30, 140.0, 0.0);  // max value
    init_v_or_w(w, 150, 140.0, 0.0);  // max value



// 4. Ferrel-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Ferrel cells v
    init_v_or_w(v, 28, 0.0, 0.0); //  center
    init_v_or_w(v, 152, 0.0, 0.0); //  center

    init_v_or_w(v, 21, v_max, - v_max);  // min value
    init_v_or_w(v, 159, - v_max, v_max);  // min value

// Ferrel cells w
    init_v_or_w(w, 28, 60.0, 0.0); //  center
    init_v_or_w(w, 152, 60.0, 0.0); //  center

    init_v_or_w(w, 21, -20.0, 0.0);  // min value
    init_v_or_w(w, 159, -20, 0.0);  // min value



// 5. Hadley-cell &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Hadley cells v
    init_v_or_w(v, 17, 0.0, 0.0); //  center
    init_v_or_w(v, 163, 0.0, 0.0); //  center

    init_v_or_w(v, 16, - v_max, v_max);  // max value
    init_v_or_w(v, 164, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 17, 40.0, 0.0); //  center
    init_v_or_w(w, 163, 40.0, 0.0); //  center

    init_v_or_w(w, 16, 90.0, 0.0);  // max value
    init_v_or_w(w, 164, 90.0, 0.0);  // max value



// 5. Ferrel-cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Hadley cells v
    init_v_or_w(v, 15, 0.0, 0.0); //  center
    init_v_or_w(v, 165, 0.0, 0.0); //  center

    init_v_or_w(v, 10, v_max, - v_max);  // min value
    init_v_or_w(v, 170, - v_max, v_max);  // min value

// Ferrel cells w
    init_v_or_w(w, 15, 40.0, 0.0); //  center
    init_v_or_w(w, 165, 40.0, 0.0); //  center

    init_v_or_w(w, 10, -20.0, 0.0);  // min value
    init_v_or_w(w, 170, -20, 0.0);  // min value



// 6. Hadley-cell &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// Hadley cells v
    init_v_or_w(v, 5, 0.0, 0.0); //  center
    init_v_or_w(v, 175, 0.0, 0.0); //  center

    init_v_or_w(v, 3, - v_max, v_max);  // max value
    init_v_or_w(v, 177, v_max, - v_max);  // max value

// Hadley cells w
    init_v_or_w(w, 5, 60.0, 0.0); //  center
    init_v_or_w(w, 175, 60.0, 0.0); //  center

    init_v_or_w(w, 3, 150.0, 0.0);  // max value
    init_v_or_w(w, 177, 150.0, 0.0);  // max value



// poles (center) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    init_v_or_w(v, 0, 0.0, 0.0); // 4. Ferrel center
    init_v_or_w(w, 180, 0.0, 0.0); // 4. Ferrel center





// forming diagonals 

// northen hemisphere

// u-velocity
    form_diagonals(u, 80, 90);
    form_diagonals(u, 70, 80);

    form_diagonals(u, 60, 70);
    form_diagonals(u, 59, 60);

    form_diagonals(u, 58, 59);
    form_diagonals(u, 57, 58);

    form_diagonals(u, 50, 57);
    form_diagonals(u, 44, 50);

    form_diagonals(u, 42, 44);
    form_diagonals(u, 40, 42);

    form_diagonals(u, 36, 40);
    form_diagonals(u, 32, 36);

    form_diagonals(u, 30, 32);
    form_diagonals(u, 28, 30);

    form_diagonals(u, 21, 28);
    form_diagonals(u, 17, 21);

    form_diagonals(u, 16, 17);
    form_diagonals(u, 15, 16);

    form_diagonals(u, 10, 15);
    form_diagonals(u, 5, 10);

    form_diagonals(u, 3, 5);
    form_diagonals(u, 0, 3);





// v-velocity
    form_diagonals(v, 80, 90);
    form_diagonals(v, 70, 80);

    form_diagonals(v, 60, 70);
    form_diagonals(v, 59, 60);

    form_diagonals(v, 58, 59);
    form_diagonals(v, 57, 58);

    form_diagonals(v, 50, 57);
    form_diagonals(v, 44, 50);

    form_diagonals(v, 42, 44);
    form_diagonals(v, 40, 42);

    form_diagonals(v, 36, 40);
    form_diagonals(v, 32, 36);

    form_diagonals(v, 30, 32);
    form_diagonals(v, 28, 30);

    form_diagonals(v, 21, 28);
    form_diagonals(v, 17, 21);

    form_diagonals(v, 16, 17);
    form_diagonals(v, 15, 16);

    form_diagonals(v, 10, 15);
    form_diagonals(v, 5, 10);

    form_diagonals(v, 3, 5);
    form_diagonals(v, 0, 3);


// w-velocity
    form_diagonals(w, 80, 90);
    form_diagonals(w, 70, 80);

    form_diagonals(w, 60, 70);
    form_diagonals(w, 59, 60);

    form_diagonals(w, 58, 59);
    form_diagonals(w, 57, 58);

    form_diagonals(w, 50, 57);
    form_diagonals(w, 44, 50);

    form_diagonals(w, 42, 44);
    form_diagonals(w, 40, 42);

    form_diagonals(w, 36, 40);
    form_diagonals(w, 32, 36);

    form_diagonals(w, 30, 32);
    form_diagonals(w, 28, 30);

    form_diagonals(w, 21, 28);
    form_diagonals(w, 17, 21);

    form_diagonals(w, 16, 17);
    form_diagonals(w, 15, 16);

    form_diagonals(w, 10, 15);
    form_diagonals(w, 5, 10);

    form_diagonals(w, 3, 5);
    form_diagonals(w, 0, 3);


    

// southen hemisphere

// u-velocity
    form_diagonals(u, 90, 100);
    form_diagonals(u, 100, 110);

    form_diagonals(u, 110, 120);
    form_diagonals(u, 120, 121);

    form_diagonals(u, 121, 122);
    form_diagonals(u, 122, 123);

    form_diagonals(u, 123, 130);
    form_diagonals(u, 130, 136);

    form_diagonals(u, 136, 138);
    form_diagonals(u, 138, 140);

    form_diagonals(u, 140, 144);
    form_diagonals(u, 144, 148);

    form_diagonals(u, 148, 150);
    form_diagonals(u, 150, 152);

    form_diagonals(u, 152, 159);
    form_diagonals(u, 159, 163);

    form_diagonals(u, 163, 164);
    form_diagonals(u, 164, 165);

    form_diagonals(u, 165, 170);
    form_diagonals(u, 170, 175);

    form_diagonals(u, 175, 177);
    form_diagonals(u, 177, 180);


// v-velocity
    form_diagonals(v, 90, 100);
    form_diagonals(v, 100, 110);

    form_diagonals(v, 110, 120);
    form_diagonals(v, 120, 121);

    form_diagonals(v, 121, 122);
    form_diagonals(v, 122, 123);

    form_diagonals(v, 123, 130);
    form_diagonals(v, 130, 136);

    form_diagonals(v, 136, 138);
    form_diagonals(v, 138, 140);

    form_diagonals(v, 140, 144);
    form_diagonals(v, 144, 148);

    form_diagonals(v, 148, 150);
    form_diagonals(v, 150, 152);

    form_diagonals(v, 152, 159);
    form_diagonals(v, 159, 163);

    form_diagonals(v, 163, 164);
    form_diagonals(v, 164, 165);

    form_diagonals(v, 165, 170);
    form_diagonals(v, 170, 175);

    form_diagonals(v, 175, 177);
    form_diagonals(v, 177, 180);


// w-velocity
    form_diagonals(w, 90, 100);
    form_diagonals(w, 100, 110);

    form_diagonals(w, 110, 120);
    form_diagonals(w, 120, 121);

    form_diagonals(w, 121, 122);
    form_diagonals(w, 122, 123);

    form_diagonals(w, 123, 130);
    form_diagonals(w, 130, 136);

    form_diagonals(w, 136, 138);
    form_diagonals(w, 138, 140);

    form_diagonals(w, 140, 144);
    form_diagonals(w, 144, 148);

    form_diagonals(w, 148, 150);
    form_diagonals(w, 150, 152);

    form_diagonals(w, 152, 159);
    form_diagonals(w, 159, 163);

    form_diagonals(w, 163, 164);
    form_diagonals(w, 164, 165);

    form_diagonals(w, 165, 170);
    form_diagonals(w, 170, 175);

    form_diagonals(w, 175, 177);
    form_diagonals(w, 177, 180);



// changing velocities from the northern hemisphere to the southern
    
    u_trans = std::vector<double>(jm, 0.0);
    v_trans = std::vector<double>(jm, 0.0);
    w_trans = std::vector<double>(jm, 0.0);

    int j_max = jm-1;

    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 0; j <= j_max; j++){
                w_trans[j] = w.x[i][j][k];
            }
            for(int j = j_max; j >= 0; j--){
                w.x[i][j_max-j][k] = w_trans[j];
            }
        }
    }

    // non dimensionalization by u_0
    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
/*
                if(is_land(SeaMount, i, j, k))     
                    u.x[i][j][k] = v.x[i][j][k] = w.x[i][j][k] = 0.0;
                else{
*/
                    u.x[i][j][k] = u.x[i][j][k]/u_0;
                    v.x[i][j][k] = v.x[i][j][k]/u_0;
                    w.x[i][j][k] = w.x[i][j][k]/u_0;

//                    u.x[i][j][k] = un.x[i][j][k] = 0.0;
//                    v.x[i][j][k] = vn.x[i][j][k] = 0.0;
//                }
            }
        }
    }

    cout << "      Saturn: SaturnCellStructure ended" << endl;
    return;
}
/*
*
*/
void cSaturnModel::form_diagonals(Array &a, int start, int end){

//    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int j = start; j < end; j++){
            for(int i = 0; i < im; i++){
                a.x[i][j][k] = (a.x[i][end][k] - a.x[i][start][k]) *
                    (j - start)/(double)(end - start) + a.x[i][start][k];
            }
        }
    }
    return;
}
/*
*
*/
void  cSaturnModel::init_u(Array &u, int j){
//    cout << endl << "      ATSAT: init_u" << endl;

    double u_max = 80.0;
    double u_center = 0.0;


    int tropopause_layer = im_tropopause[j];
    double tropopause_height = get_layer_height(tropopause_layer);

    for(int k = 0; k < km; k++){
        for(int i = 0; i < tropopause_layer; i++){

            double layer_height = get_layer_height(i), 
                half_tropopause_height = tropopause_height/1.5;
            double ratio;
 
            if(layer_height < half_tropopause_height){    
                ratio = layer_height/(half_tropopause_height/0.5);



                switch(j){

                    case 0:  u.x[i][0][k] = - u_max * ratio; break;
                    case 3:  u.x[i][3][k] = u_center * ratio; break;

                    case 5:  u.x[i][5][k] = u_max * ratio; break;
                    case 10:  u.x[i][10][k] = u_center * ratio; break;

                    case 15:  u.x[i][15][k] = - u_max * ratio; break;
                    case 16:  u.x[i][16][k] = u_center * ratio; break;

                    case 17:  u.x[i][17][k] = u_max * ratio; break;
                    case 21:  u.x[i][21][k] = u_center * ratio; break;

                    case 28:  u.x[i][28][k] = - u_max * ratio; break;
                    case 30:  u.x[i][30][k] = u_center * ratio; break;

                    case 32:  u.x[i][32][k] = u_max * ratio; break;
                    case 36:  u.x[i][36][k] = u_center * ratio; break;

                    case 40:  u.x[i][40][k] = - u_max * ratio; break;
                    case 42:  u.x[i][42][k] = u_center * ratio; break;

                    case 44:  u.x[i][44][k] = u_max * ratio; break;
                    case 50:  u.x[i][50][k] = u_center * ratio; break;

                    case 57:  u.x[i][57][k] = - u_max * ratio; break;
                    case 58:  u.x[i][58][k] = u_center * ratio; break;

                    case 59:  u.x[i][59][k] = u_max * ratio; break;
                    case 60:  u.x[i][60][k] = u_center * ratio; break;

                    case 70:  u.x[i][70][k] = - u_max * ratio; break; 
                    case 80:  u.x[i][80][k] = u_center * ratio; break;


                    case 90:  u.x[i][90][k] = u_max * ratio; break;


                    case 100:  u.x[i][100][k] = u_center * ratio; break;
                    case 110:  u.x[i][110][k] = - u_max * ratio; break;

                    case 120:  u.x[i][120][k] = u_center * ratio; break;
                    case 121:  u.x[i][121][k] = u_max * ratio; break;

                    case 122:  u.x[i][122][k] = u_center * ratio; break;
                    case 123:  u.x[i][123][k] = - u_max * ratio; break;

                    case 130:  u.x[i][130][k] = u_center * ratio; break;
                    case 136:  u.x[i][136][k] = u_max * ratio; break;

                    case 138:  u.x[i][138][k] = u_center * ratio; break;
                    case 140:  u.x[i][140][k] = - u_max * ratio; break;

                    case 144:  u.x[i][144][k] = u_center * ratio; break;
                    case 148:  u.x[i][148][k] = u_max * ratio; break;

                    case 150:  u.x[i][150][k] = u_center * ratio; break;
                    case 152:  u.x[i][152][k] = - u_max * ratio; break;

                    case 159:  u.x[i][159][k] = u_center * ratio; break;
                    case 163:  u.x[i][163][k] = u_max * ratio; break;

                    case 164:  u.x[i][164][k] = u_center * ratio; break;
                    case 165:  u.x[i][165][k] = - u_max * ratio; break;

                    case 170:  u.x[i][170][k] = u_center * ratio; break;
                    case 175:  u.x[i][175][k] = u_max * ratio; break;

                    case 177:  u.x[i][177][k] = u_center * ratio; break;
                    case 180:  u.x[i][180][k] = - u_max * ratio; break;

                }
            }
            else{
                ratio = (tropopause_height-layer_height)/half_tropopause_height;
                switch(j){

                    case 0:  u.x[i][0][k] = - u_max * ratio; break;
                    case 3:  u.x[i][3][k] = u_center * ratio; break;

                    case 5:  u.x[i][5][k] = u_max * ratio; break;
                    case 10:  u.x[i][10][k] = u_center * ratio; break;

                    case 15:  u.x[i][15][k] = - u_max * ratio; break;
                    case 16:  u.x[i][16][k] = u_center * ratio; break;

                    case 17:  u.x[i][17][k] = u_max * ratio; break;
                    case 21:  u.x[i][21][k] = u_center * ratio; break;

                    case 28:  u.x[i][28][k] = - u_max * ratio; break;
                    case 30:  u.x[i][30][k] = u_center * ratio; break;

                    case 32:  u.x[i][32][k] = u_max * ratio; break;
                    case 36:  u.x[i][36][k] = u_center * ratio; break;

                    case 40:  u.x[i][40][k] = - u_max * ratio; break;
                    case 42:  u.x[i][42][k] = u_center * ratio; break;

                    case 44:  u.x[i][44][k] = u_max * ratio; break;
                    case 50:  u.x[i][50][k] = u_center * ratio; break;

                    case 57:  u.x[i][57][k] = - u_max * ratio; break;
                    case 58:  u.x[i][58][k] = u_center * ratio; break;

                    case 59:  u.x[i][59][k] = u_max * ratio; break;
                    case 60:  u.x[i][60][k] = u_center * ratio; break;

                    case 70:  u.x[i][70][k] = - u_max * ratio; break; 
                    case 80:  u.x[i][80][k] = u_center * ratio; break;


                    case 90:  u.x[i][90][k] = u_max * ratio; break;


                    case 100:  u.x[i][100][k] = u_center * ratio; break;
                    case 110:  u.x[i][110][k] = - u_max * ratio; break;

                    case 120:  u.x[i][120][k] = u_center * ratio; break;
                    case 121:  u.x[i][121][k] = u_max * ratio; break;

                    case 122:  u.x[i][122][k] = u_center * ratio; break;
                    case 123:  u.x[i][123][k] = - u_max * ratio; break;

                    case 130:  u.x[i][130][k] = u_center * ratio; break;
                    case 136:  u.x[i][136][k] = u_max * ratio; break;

                    case 138:  u.x[i][138][k] = u_center * ratio; break;
                    case 140:  u.x[i][140][k] = - u_max * ratio; break;

                    case 144:  u.x[i][144][k] = u_center * ratio; break;
                    case 148:  u.x[i][148][k] = u_max * ratio; break;

                    case 150:  u.x[i][150][k] = u_center * ratio; break;
                    case 152:  u.x[i][152][k] = - u_max * ratio; break;

                    case 159:  u.x[i][159][k] = u_center * ratio; break;
                    case 163:  u.x[i][163][k] = u_max * ratio; break;

                    case 164:  u.x[i][164][k] = u_center * ratio; break;
                    case 165:  u.x[i][165][k] = - u_max * ratio; break;

                    case 170:  u.x[i][170][k] = u_center * ratio; break;
                    case 175:  u.x[i][175][k] = u_max * ratio; break;

                    case 177:  u.x[i][177][k] = u_center * ratio; break;
                    case 180:  u.x[i][180][k] = - u_max * ratio; break;

                }
            }

        }
    }

//    cout << "      ATSAT: init_u ended" << endl;
    return;
}
/*
*
*/
void  cSaturnModel::init_v_or_w(Array &v_or_w, int j, double coeff_trop, double coeff_sl){
//    cout << endl << "      AGCM: init_v_or_w" << endl;
//    int tropopause_layer = get_tropopause_layer(j);
    int tropopause_layer = im_tropopause[j];
    double tropopause_height = get_layer_height(tropopause_layer);
    for(int k = 0; k < km; k++){
        if(is_ocean_surface(SeaMount, 20, j, k)){ // velocity v changes sign at i=20
            coeff_sl = v_or_w.x[0][j][k];
        }
        for(int i = 0; i < tropopause_layer; i++){
            v_or_w.x[i][j][k] = (coeff_trop - coeff_sl) 
                * get_layer_height(i)/tropopause_height + coeff_sl;
        }
    }
    init_v_or_w_above_tropopause(v_or_w, j, coeff_trop);

//    cout << "      AGCM: init_v_or_w ended" << endl;
    return;
}
/*
*
*/
void  cSaturnModel::init_v_or_w_above_tropopause(Array &v_or_w, int j, double coeff){
//     cout << endl << "      AGCM: init_v_or_w_above_tropopause" << endl;
//    int tropopause_layer = get_tropopause_layer(j);
    int tropopause_layer = im_tropopause[j];
    if(tropopause_layer >= im-1) return;
    double tropopause_height = get_layer_height(tropopause_layer);

//    #pragma omp parallel for
    for(int k = 0; k < km; k++){
        for(int i = tropopause_layer; i < im; i++){
            v_or_w.x[i][j][k] = coeff 
                * (get_layer_height(im-1) - get_layer_height(i))
                /(get_layer_height(im-1) - tropopause_height);
//                v.x[i][j][k] = - v.x[i][j][k];
        }
    }
//     cout << "      AGCM: init_v_or_w_above_tropopause ended" << endl;
    return;
}
