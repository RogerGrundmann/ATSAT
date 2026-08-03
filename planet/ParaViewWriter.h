/*
 * SHARED OUTPUT — the ParaView file machinery. MUST BE BYTE-IDENTICAL IN EVERY MODEL THAT USES
 * IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * The second non-physics shared header, for the same reason as Reporting.h: measured against
 * ATJUP, ParaView_Sat.cpp overlapped 76.5 % line-for-line, the second-highest of any pair, and
 * it is scaffolding rather than science — the code a third and fourth planet copy verbatim.
 *
 * ===== WHAT IS HERE, AND THE ONE THING THAT DECIDED THE SPLIT =====
 *
 * Here: the five dumpers, the .vtk slice machinery (open + close), and the PlotData writer. All
 * of it is geometry and file format — no field of either planet appears in it.
 *
 * NOT here: the field LISTS. Each writer ends in a run of dump_*() calls naming which arrays go
 * into the file with what coefficient, and those lists are where the two models disagree. The
 * disagreement is not cosmetic and not about which fields exist; it is the SAME units question
 * that keeps printMinMax unshared, and ParaView carries more of it:
 *
 *     species        ATSAT 1e3            ATJUP 1.0
 *     nh4sh          ATSAT 1e3            ATJUP r_mix * 1e6
 *     temperature    ATSAT t*t_ref/10.0   ATJUP t*t_ref - 273.15   (K/10 against degC)
 *     PressureDyn    ATSAT 1.0            ATJUP p_dyn_to_bar()*1e3
 *     thermalflux    ATSAT 1e3            ATJUP 1e-3
 *     Centrifugal    ATSAT 1e3            ATJUP 1.0
 *
 * Sharing the lists would mean choosing one model's reading of the units for both, silently, in
 * a refactor. ATJUP has settled its side (its note in PrintMsg_Jup.cpp records that the species
 * factor was wrong by 1.2844x and was corrected); ATSAT has not had that correction. So the
 * lists stay in each planet's own file until the question is settled on its own terms, and this
 * header takes only what has no opinion about units.
 *
 * That split is not a compromise — it is the same line Reporting.h draws, and drawing it twice
 * in the same place is evidence it is the right one. What a new planet inherits is every line
 * that says HOW to write a ParaView file; what it must still write is the one list saying WHICH
 * of its fields to put in one.
 *
 * ===== WHAT IS NOT HERE YET =====
 *
 * paraview_panorama_vts and paraview_sphere_vts keep their own file headers. Both are .vts
 * writers whose header has to name every scalar in one attribute string, so their machinery is
 * bound up with their field list in a way the .vtk slices are not; separating them is a second
 * job. The dumpers and PlotData below are shared by all of them regardless.
 */

#pragma once

#include "ATPhys.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

/*
 * The five dumpers. These were already byte-identical in both models — in the ParaViewSaturn and
 * ParaViewJupiter namespaces — before this header existed, so they move across unchanged.
 *
 * They are free functions in a namespace rather than members of the template below because that
 * is what lets each model keep saying `using namespace ParaViewIO;` at the top of a writer and
 * call `dump_radial(...)` unqualified, exactly as it called its own copy. None of them needs the
 * planet type: an Array knows its own extents.
 */
namespace ParaViewIO{

    inline void dump_array(const std::string &name, Array &a, double multiplier, std::ofstream &f){
        f <<  "    <DataArray type=\"Float32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for(int k = 0; k < a.km; k++){
            for(int j = 0; j < a.jm; j++){
                for(int i = 0; i < a.im; i++){
                    f << (a.x[i][j][k] * multiplier) << std::endl;
                }
                f << "\n";
            }
            f << "\n";
        }
        f << "\n";
        f << "    </DataArray>\n";
    }

    inline void dump_radial(const std::string &desc, Array &a, double multiplier, int i, std::ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << std::endl;
            }
        }
    }

    inline void dump_radial_2d(const std::string &desc, Array_2D &a, double multiplier, std::ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for(int j = 0; j < a.jm; j++){
            for(int k = 0; k < a.km; k++){
                f << (a.y[j][k] * multiplier) << std::endl;
            }
        }
    }

    inline void dump_zonal(const std::string &desc, Array &a, double multiplier, int k, std::ofstream &f){
        f <<  "SCALARS " << desc << " float " << 1 << std::endl;
        f <<  "LOOKUP_TABLE default" << std::endl;
        for(int i = 0; i < a.im; i++){
            for(int j = 0; j < a.jm; j++){
                f << (a.x[i][j][k] * multiplier) << std::endl;
            }
        }
    }

    inline void dump_longal(const std::string &desc, Array &a, double multiplier, int j, std::ofstream &f){
        f << "SCALARS " << desc << " float " << 1 << std::endl;
        f << "LOOKUP_TABLE default" << std::endl;
        for(int i = 0; i < a.im; i++){
            for(int k = 0; k < a.km; k++){
                f << (a.x[i][j][k] * multiplier) << std::endl;
            }
        }
    }
}

template<class Planet>
class ParaViewWriter{

    Planet &m;

public:

    explicit ParaViewWriter(Planet &m_) : m(m_) {}

    /*
     * Open one .vtk slice file and write everything up to and including POINT_DATA — the header,
     * the structured-grid declaration and the point coordinates. The caller then appends its own
     * field list and calls close_slice().
     *
     * The three slice writers differ only in the numbers passed here. Their geometry loops were
     * three copies of one nest, distinguished by which coordinate advances on the inner index:
     *
     *     radial   title "Radial"        n_fast km, n_slow jm, step 0.1,   axis y
     *     zonal    title "Zonal"         n_fast jm, n_slow im, step 0.05,  axis y
     *     longal   title "Longitudinal"  n_fast km, n_slow im, step 0.025, axis z
     *
     * The outer step is 0.1 in all three. `inner_is_z` picks which of y and z advances; the other
     * stays 0, which is what made the longal copy look unlike the other two.
     *
     * Returned by value — ofstream has been movable since C++11 and each writer keeps it as a
     * local under its old name, so the field list below it needs no editing.
     */
    std::ofstream open_slice(const char *kind, const char *title, int idx, int n,
                             int n_fast, int n_slow, double inner_step, bool inner_is_z) const
    {
        const std::string name = m.output_path + "/" + Planet::planet_name() + "_"
            + kind + "_" + std::to_string(idx) + "_" + std::to_string(n) + ".vtk";

        std::ofstream f;
        f.precision(4);
        f.setf(std::ios::fixed);
        f.open(name);
        if(!f.is_open()){
            std::cerr << "ERROR: could not open paraview_vtk file " << __FILE__
                << " at line " << __LINE__ << "\n";
            abort();
        }

        f <<  "# vtk DataFile Version 3.0" << std::endl;
        f <<  title << "_Data_" << Planet::planet_short() << "_Circulation\n";
        f <<  "ASCII" << std::endl;
        f <<  "DATASET STRUCTURED_GRID" << std::endl;
        f <<  "DIMENSIONS " << n_fast << " "<< n_slow << " " << 1 << std::endl;
        f <<  "POINTS " << n_slow * n_fast << " float" << std::endl;

        double x = 0.0, y = 0.0, z = 0.0;
        const double dx = 0.1;
        for(int s = 0; s < n_slow; s++){
            for(int q = 0; q < n_fast; q++){
                double &moving = inner_is_z ? z : y;
                if(q == 0) moving = 0.0;
                else       moving = moving + inner_step;
                f << x << " " << y << " "<< z << std::endl;
            }
            if(inner_is_z) z = 0.0; else y = 0.0;
            x = x + dx;
        }
        f <<  "POINT_DATA " << n_slow * n_fast << std::endl;
        return f;
    }

    /*
     * Close the file and report it. THE ONE PLACE THIS HEADER CHANGES AN EXISTING OUTPUT, and
     * only a log line — no data file is affected.
     *
     * The two models did not print the same thing, and one of them was malformed. ATJUP printed
     *     File:  Jupiter_radial_20_1.vtk  has been written to Directory: ...
     * while ATSAT built the name twice and nested one inside the other:
     *     File:  [Saturn_radial_20_1.vtk]_Sat_radial_20_1.vtk  has been written to Directory: ...
     * — the bracketed part is the whole file name, and everything after it is a second, differently
     * spelled attempt at the same name. There is no reading in which that is intended.
     *
     * So this is ATJUP's form, adopted for both. It is a deliberate, visible change to ATSAT's log
     * and it is the only one in this commit: the .vtk and .xyz files themselves are byte-identical
     * before and after, which is what the acceptance test checks.
     */
    void close_slice(std::ofstream &f, const char *kind, int idx, int n) const {
        f.close();
        std::cout << "   File:  " << Planet::planet_name() << "_" << kind << "_"
            << idx << "_" << n << ".vtk"
            << "  has been written to Directory:  " << m.output_path << std::endl;
    }

    /*
     * The surface .xyz table. This one was 100 % identical between the two models — the only
     * difference was the planet's name in the file name — so nothing here is parameterised
     * beyond that.
     *
     * Note what it writes against what its own header row promises: the header names TEN columns
     * and the loop writes NINE, so every column from "temperature(Celsius)" rightward is read
     * under the wrong heading, and "precipitation(mm)" and "precipitable water(mm)" have no data
     * at all. That is carried across unchanged and deliberately — it is a real defect, it is
     * identical in both models, and correcting it changes an output file something downstream
     * may parse. It wants its own commit and its own decision about which columns were intended.
     */
    void plot_data() const {
        const std::string name = m.output_path + "/PlotData_" + Planet::planet_name() + ".xyz";
        std::ofstream f;
        f.precision(4);
        f.setf(std::ios::fixed);
        f.open(name);
        if(!f.is_open()){
            std::cerr << "ERROR: could not open PlotData file " << __FILE__
                      << " at line " << __LINE__ << "\n";
            abort();
        }
        f << "lons(deg)" << ", " << "lats(deg)" << ", "
            << "topography" << ", " << "v-velocity(m/s)" << ", "
            << "w-velocity(m/s)" << ", " << "velocity-mag(m/s)" << ", "
            << "temperature(Celsius)" << ", " << "water_vapour(g/kg)"
            << ", " << "precipitation(mm)" << ", "
            <<  "precipitable water(mm)" << std::endl;
        double vel_mag;
        for(int k = 0; k < m.km; k++){
            for(int j = 0; j < m.jm; j++){
                vel_mag = sqrt(pow(m.v.x[0][j][k] * m.u_0, 2) + pow(m.w.x[0][j][k] * m.u_0, 2));
                f << k << " " << j << " " << m.SeaMount.x[0][j][k] << " "
                    << m.v.x[0][j][k] * m.u_0 << " " << m.w.x[0][j][k] * m.u_0 << " "
                    << vel_mag << " " << m.t.x[0][j][k] * m.t_ref - m.t_ref << " "
                    << m.h2o.x[0][j][k] << " "<< m.nh3.x[0][j][k] <<  std::endl;
            }
        }
        f.close();
        return;
    }
};
