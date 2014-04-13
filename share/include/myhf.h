/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef HF_H
#define HF_H
/*
 *  hf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */


#include <vector>
#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libmints/vector.h>
#include <libdiis/diismanager.h>
#include <libdiis/diisentry.h>
#include <psi4-dec.h>
#include <libqt/qt.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Matrix;
class Vector;
class SimpleVector;
class TwoBodySOInt;
class JK;
namespace scf {

class myHF : public Wavefunction {
protected:

    SharedMatrix SO2AO_;
    
    /// The kinetic energy matrix
    SharedMatrix T_;
    /// The 1e potential energy matrix
    SharedMatrix V_;
    /// The orthogonalization matrix (symmetric or canonical)
    SharedMatrix X_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_temp_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_F_temp_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_C_temp_;

    /// Old C Alpha matrix (if needed for MOM)
    SharedMatrix Ca_old_;
    /// Old C Beta matrix (if needed for MOM)
    SharedMatrix Cb_old_;

    /// Previous iteration's energy and current energy
    double Eold_;
    double E_;

    /// Table of energy components
    std::map<std::string, double> energies_;

    /// The RMS error in the density
    double Drms_;

    /// Max number of iterations for myHF
    int maxiter_;

    /// Fail if we don't converge by maxiter?
    bool fail_on_maxiter_;

    /// Current Iteration
    int iteration_;

    /// Nuclear repulsion energy
    double nuclearrep_;

    /// Whether DIIS was performed this iteration, or not
    bool diis_performed_;

    /// DOCC vector from input (if found)
    bool input_docc_;

    /// SOCC vector from input (if found)
    bool input_socc_;

    //Initial SAD doubly occupied may be more than ndocc
    int sad_nocc_[8];

    /// Mapping arrays
    int *so2symblk_;
    int *so2index_;

    /// SCF algorithm type
    std::string scf_type_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// How big of a perturbation
    double lambda_;
    /// With what...
    enum perturb { nothing, dipole_x, dipole_y, dipole_z, embpot, dx, sphere };
    perturb perturb_;

    /// The value below which integrals are neglected
    double integral_threshold_;

    /// The soon to be ubiquitous JK object
    boost::shared_ptr<JK> jk_;

    /// Are we to do MOM?
    bool MOM_enabled_;
    /// Are we to do excited-state MOM?
    bool MOM_excited_;
    /// MOM started?
    bool MOM_started_;
    /// MOM performed?
    bool MOM_performed_;

    /// Are we to fractionally occupy?
    bool frac_enabled_;
    /// Frac started? (Same thing as frac_performed_)
    bool frac_performed_;

    /// DIIS manager intiialized?
    bool initialized_diis_manager_;
    /// DIIS manager for all SCF wavefunctions
    boost::shared_ptr<DIISManager> diis_manager_;

    /// How many min vectors for DIIS
    int min_diis_vectors_;
    /// How many max vectors for DIIS
    int max_diis_vectors_;
    /// When do we start collecting vectors for DIIS
    int diis_start_;
    /// Are we even using DIIS?
    int diis_enabled_;

    /// The amount (%) of the previous orbitals to mix in during SCF damping
    double damping_percentage_;
    /// The energy convergence at which SCF damping is disabled
    double damping_convergence_;
    /// Whether to use SCF damping
    bool damping_enabled_;
    /// Whether damping was actually performed this iteration
    bool damping_performed_;

    // parameters for hard-sphere potentials
    double radius_; // radius of spherical potential
    double thickness_; // thickness of spherical barrier
    int r_points_; // number of radial integration points
    int theta_points_; // number of colatitude integration points
    int phi_points_; // number of azimuthal integration points

public:
    /// Nuclear contributions
    Vector nuclear_dipole_contribution_;
    Vector nuclear_quadrupole_contribution_;

    /// The number of iterations needed to reach convergence
    int iterations_needed() {return iterations_needed_;}

    /// The JK object (or null if it has been deleted)
    boost::shared_ptr<JK> jk() const { return jk_; }

    /// The RMS error in the density
    double rms_density_error() {return Drms_;}
protected:

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    // Temporarily converting to virtual function for testing embedding
    // potentials.  TDC, 5/23/12.
    virtual void form_H();
    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();

    /// Prints the orbital occupation
    void print_occupation();

    /// Perform casting of C from old basis to new basis if desired.
    SharedMatrix BasisProjection(SharedMatrix Cold, int* napi, boost::shared_ptr<BasisSet> old_basis, boost::shared_ptr<BasisSet> new_basis);

    /// Common initializer
    void common_init();

    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    //virtual void finalize();

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation();

    /// Maximum overlap method for prevention of oscillation/excited state SCF
    void MOM();
    /// Start the MOM algorithm (requires one iteration worth of setup)
    void MOM_start();

    /// Fractional occupation UHF/UKS
    void frac();
    /// Renormalize orbitals to 1.0 before saving to chkpt
    void frac_renormalize();

    /// Check the stability of the wavefunction, and correct (if requested)
    virtual void stability_analysis();
    void print_stability_analysis(std::vector<std::pair<double, int> > &vec);


    /// Determine how many core and virtual orbitals to freeze
    void compute_fcpi();
    void compute_fvpi();

    /// Prints the orbitals energies and symmetries (helper method)
    void print_orbitals(const char* header, std::vector<std::pair<double,
                        std::pair<const char*, int> > > orbs);

    /// Prints the orbitals in arbitrary order (works with MOM)
    void print_orbitals();

    /// Prints the energy breakdown from this SCF
    void print_energies();

    /// Prints some opening information
    void print_header();

    /// Prints some details about nsopi/nmopi, and initial occupations
    void print_preiterations();

    /// Do any needed integral setup
    virtual void integrals();

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number_;

    /// The number of electrons
    int nelectron_;

    /// The charge of the system
    int charge_;

    /// The multiplicity of the system (specified as 2 Ms + 1)
    int multiplicity_;

    /// The number of iterations need to reach convergence
    int iterations_needed_;

    /// Compute energy for the iteration.
    virtual double compute_E() = 0;

    /// Save the current density and energy.
    virtual void save_density_and_energy() = 0;

    /// Check MO phases
    void check_phases();

    /// SAD Guess and propagation
    void compute_SAD_guess();

    /// Reset to regular occupation from the fractional occupation
    void reset_SAD_occupation();

    /// Form the guess (gaurantees C, D, and E)
    virtual void guess();

    /** Computes the density matrix (D_) */
    virtual void form_D() =0;

    /** Applies damping to the density update */
    virtual void damp_update();

    /** Compute the MO coefficients (C_) */
    virtual void form_C() =0;

    /** Transformation, diagonalization, and backtransform of Fock matrix */
    virtual void diagonalize_F(const SharedMatrix& F, SharedMatrix& C, boost::shared_ptr<Vector>& eps);

    /** Computes the Fock matrix */
    virtual void form_F() =0;

    /** Computes the initial MO coefficients (default is to call form_C) */
    virtual void form_initial_C() { form_C(); }

    /** Forms the G matrix */
    virtual void form_G() =0;

    /** Computes the initial energy. */
    virtual double compute_initial_E() { return 0.0; }

    /** Test convergence of the wavefunction */
    virtual bool test_convergency() { return false; }

    /** Compute/print spin contamination information (if unrestricted) **/
    virtual void compute_spin_contamination();

    /** Saves information to the checkpoint file */
    virtual void save_information() {}

    /** Compute the orbital gradient */
    virtual void compute_orbital_gradient(bool save_diis) {}

    /** Performs DIIS extrapolation */
    virtual bool diis() { return false; }

    /** Form Fia (for DIIS) **/
    virtual SharedMatrix form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi);

    /** Form X'(FDS - SDF)X (for DIIS) **/
    virtual SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso);

    /** Save orbitals to use later as a guess **/
    virtual void save_orbitals();

    /** Load orbitals from previous computation, projecting if needed **/
    virtual void load_orbitals();

    /** Save SAPT info (TODO: Move to Python driver **/
    virtual void save_sapt_info() {}

    /** Saves all wavefunction information to the checkpoint file*/
    void dump_to_checkpoint();

public:
    myHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    myHF(Options& options, boost::shared_ptr<PSIO> psio);
    
    SharedMatrix sotoao() { return SO2AO_;} // wavefunction already has aotoso() 
    double EHF() { return E_; }

    virtual ~myHF();
    
    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    virtual void finalize();

    virtual double compute_energy();
};

}} // Namespaces

#endif
/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 *  hf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include <libmints/mints.h>

#include <libfunctional/superfunctional.h>
#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <liboptions/liboptions_python.h>
#include <psifiles.h>
#include <libfock/jk.h>

//#include "hf.h"

#include <psi4-dec.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

myHF::myHF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6)
{
    common_init();
}

myHF::myHF(Options& options, boost::shared_ptr<PSIO> psio)
    : Wavefunction(options, psio),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6)
{
    common_init();
}

myHF::~myHF()
{
}

void myHF::common_init()
{
    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
    SO2AO_ = pet->sotoao();
    //SO2AO_->identity();
    //AO2SO_->identity();

    attempt_number_ = 1;

    // This quantity is needed fairly soon
    nirrep_ = factory_->nirrep();

    integral_threshold_ = options_.get_double("INTS_TOLERANCE");

    scf_type_ = options_.get_str("SCF_TYPE");

    H_.reset(factory_->create_matrix("One-electron Hamiltonion"));
    X_.reset(factory_->create_matrix("X"));

    nmo_ = 0;
    nso_ = 0;
    int* dimpi = factory_->colspi();
    for (int h = 0; h< factory_->nirrep(); h++){
        nsopi_[h] = dimpi[h];
        nmopi_[h] = nsopi_[h]; //For now
        nso_ += nsopi_[h];
        nmo_ += nmopi_[h]; //For now
    }

    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Should we continue if we fail to converge?
    fail_on_maxiter_ = options_.get_bool("FAIL_ON_MAXITER");

    // Read in DOCC and SOCC from memory
    int nirreps = factory_->nirrep();
    int ndocc = 0, nsocc = 0;
    input_docc_ = false;
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        // Map the symmetry of the input DOCC, to account for displacements
        boost::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            int full_nirreps = old_pg->char_table().nirrep();
            if(options_["DOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            int *temp_docc = new int[full_nirreps];
            for(int h = 0; h < full_nirreps; ++h)
                temp_docc[h] = options_["DOCC"][h].to_integer();
            map_irreps(temp_docc);
            doccpi_ = temp_docc;
            delete[] temp_docc;
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["DOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h)
                doccpi_[h] = options_["DOCC"][h].to_integer();
        }
        for (int i=0; i<nirreps; ++i)
            ndocc += 2*doccpi_[i];
    } else {
        for (int i=0; i<nirreps; ++i)
            doccpi_[i] = 0;
    }

    if(options_.get_str("REFERENCE") == "RKS" || options_.get_str("REFERENCE") == "UKS")
        name_ = "DFT";
    else
        name_ = "SCF";

    input_socc_ = false;
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        // Map the symmetry of the input SOCC, to account for displacements
        boost::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            int full_nirreps = old_pg->char_table().nirrep();
            if(options_["SOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            int *temp_socc = new int[full_nirreps];
            for(int h = 0; h < full_nirreps; ++h)
                temp_socc[h] = options_["SOCC"][h].to_integer();
            map_irreps(temp_socc);
            soccpi_ = temp_socc;
            delete[] temp_socc;
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["SOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h)
                soccpi_[h] = options_["SOCC"][h].to_integer();
        }
        for (int i=0; i<nirreps; ++i)
            nsocc += soccpi_[i];
    } else {
        for (int i=0; i<nirreps; ++i)
            soccpi_[i] = 0;
    }

    if (input_socc_ || input_docc_) {
        for (int h = 0; h < nirrep_; h++) {
            nalphapi_[h] = doccpi_[h] + soccpi_[h];
            nbetapi_[h]  = doccpi_[h];
        }
    }


    // Read information from checkpoint
    nuclearrep_ = molecule_->nuclear_repulsion_energy();

    // Determine the number of electrons in the system
    charge_ = molecule_->molecular_charge();
    nelectron_  = 0;
    for (int i=0; i<molecule_->natom(); ++i)
        nelectron_ += (int)molecule_->Z(i);
    nelectron_ -= charge_;

    // If the user told us the multiplicity, read it from the input
    if(molecule_->multiplicity_specified()){
        multiplicity_ = molecule_->multiplicity();
    }else{
        if(nelectron_%2){
            multiplicity_ = 2;
            molecule_->set_multiplicity(2);
            if (WorldComm->me() == 0) {
            // There are an odd number of electrons
                fprintf(outfile,"\tThere are an odd number of electrons - assuming doublet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");
            }
        }else{
            multiplicity_ = 1;
            // There are an even number of electrons
            if (WorldComm->me() == 0) {
                fprintf(outfile,"\tThere are an even number of electrons - assuming singlet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");
            }
        }
    }

    // Make sure that the multiplicity is reasonable
    if(multiplicity_ - 1 > nelectron_){
        char *str = new char[100];
        sprintf(str, "There are not enough electrons for multiplicity = %d, \n"
                     "please check your input and use the MULTP keyword", multiplicity_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }
    if(multiplicity_%2 == nelectron_%2){
        char *str = new char[100];
        sprintf(str, "A multiplicity of %d with %d electrons is impossible.\n"
                     "Please check your input and use the MULTP and/or CHARGE keywords",
                     multiplicity_, nelectron_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }

    nbeta_  = (nelectron_ - multiplicity_ + 1)/2;
    nalpha_ = nbeta_ + multiplicity_ - 1;

    perturb_h_ = false;
    perturb_h_ = options_.get_bool("PERTURB_H");
    perturb_ = nothing;
    lambda_ = 0.0;
    if (perturb_h_) {
        string perturb_with;

        lambda_ = options_.get_double("PERTURB_MAGNITUDE");

        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X")
                perturb_ = dipole_x;
            else if (perturb_with == "DIPOLE_Y")
                perturb_ = dipole_y;
            else if (perturb_with == "DIPOLE_Z")
                perturb_ = dipole_z;
            else if (perturb_with == "EMBPOT") {
                perturb_ = embpot;
                lambda_ = 1.0;
            }
            else if (perturb_with == "DX") {
                perturb_ = dx;
                lambda_ = 1.0;
            }
            else if (perturb_with == "SPHERE") {
                perturb_ = sphere;
                lambda_ = 1.0;
            }
            else {
                if (WorldComm->me() == 0) {
                    fprintf(outfile, "Unknown PERTURB_WITH. Applying no perturbation.\n");
                }
            }
        } else {
            if (WorldComm->me() == 0) {
                fprintf(outfile, "PERTURB_H is true, but PERTURB_WITH not found, applying no perturbation.\n");
            }
        }
    }

    // How much stuff shall we echo to the user?
    if(options_["PRINT"].has_changed())
        print_ = options_.get_int("PRINT");

    if(options_["DAMPING_PERCENTAGE"].has_changed()){
        // The user has asked for damping to be turned on
        damping_enabled_ = true;
        damping_percentage_ = options_.get_double("DAMPING_PERCENTAGE") / 100.0;
        if(damping_percentage_ < 0.0 || damping_percentage_ > 1.0)
            throw PSIEXCEPTION("DAMPING_PERCENTAGE must be between 0 and 100.");
        damping_convergence_ = options_.get_double("DAMPING_CONVERGENCE");
    }else{
        damping_enabled_ = false;
    }

    // Handle common diis info
    diis_enabled_ = true;
    min_diis_vectors_ = 4;

    // Allocate memory for DIISmin_diis_vectors_
    //  First, did the user request a different number of diis vectors?
    min_diis_vectors_ = options_.get_int("DIIS_MIN_VECS");
    max_diis_vectors_ = options_.get_int("DIIS_MAX_VECS");
    diis_start_ = options_.get_int("DIIS_START");
    diis_enabled_ = options_.get_bool("DIIS");

    // Don't perform DIIS if less than 2 vectors requested, or user requested a negative number
    if (min_diis_vectors_ < 2) {
        // disable diis
        diis_enabled_ = false;
    }

    initialized_diis_manager_ = false;

    MOM_enabled_ = (options_.get_int("MOM_START") != 0);
    MOM_excited_ = (options_["MOM_OCC"].size() != 0 && MOM_enabled_);
    MOM_started_ = false;
    MOM_performed_ = false;

    frac_enabled_ = (options_.get_int("FRAC_START") != 0);
    frac_performed_ = false;

    print_header();
}

void myHF::damp_update()
{
    throw PSIEXCEPTION("Sorry, damping has not been implemented for this "
                       "type of SCF wavefunction yet.");
}

void myHF::integrals()
{
    if (print_ && WorldComm->me() == 0)
        fprintf(outfile, "  ==> Integral Setup <==\n\n");

    // Build the JK from options, symmetric type
    try {
        jk_ = JK::build_JK();
    }
    catch(const BasisSetNotFound& e) {
        if (options_.get_str("SCF_TYPE") == "DF" || options_.get_int("DF_SCF_GUESS") == 1) {
            fprintf(outfile, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            fprintf(outfile, "%s\n", e.what());
            fprintf(outfile, "   Turning off DF and switching to PK method.\n");
            fprintf(outfile, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            options_.set_str("SCF", "SCF_TYPE", "PK");
            options_.set_bool("SCF", "DF_SCF_GUESS", false);
            jk_ = JK::build_JK();
        }
        else
            throw; // rethrow the error
    }

    // Tell the JK to print
    jk_->set_print(print_);
    // Give the JK 75% of the memory
    jk_->set_memory((ULI)(options_.get_double("SCF_MEM_SAFETY_FACTOR")*(Process::environment.get_memory() / 8L)));

    // DFT sometimes needs custom stuff
    if ((options_.get_str("REFERENCE") == "UKS" || options_.get_str("REFERENCE") == "RKS")) {

        // Need a temporary functional
        boost::shared_ptr<SuperFunctional> functional = 
            SuperFunctional::current(options_);
        
        // K matrices
        jk_->set_do_K(functional->is_x_hybrid());
        // wK matrices 
        jk_->set_do_wK(functional->is_x_lrc());
        // w Value
        jk_->set_omega(functional->x_omega());
    }

    // Initialize
    jk_->initialize(); 
    // Print the header
    jk_->print_header();
}

void myHF::finalize()
{
    // This will be the only one
    if (!options_.get_bool("SAVE_JK")) {
        jk_.reset();
    }

    // Clean up after DIIS
    if(initialized_diis_manager_)
        diis_manager_->delete_diis_file();
    diis_manager_.reset();
    initialized_diis_manager_ = false;

    // Figure out how many frozen virtual and frozen core per irrep
    compute_fcpi();
    compute_fvpi();
    energy_ = E_;

    dump_to_checkpoint();

    //Sphalf_.reset();
    X_.reset();
    T_.reset();
    V_.reset();
    diag_temp_.reset();
    diag_F_temp_.reset();
    diag_C_temp_.reset();

    // Close the chkpt
    if(psio_->open_check(PSIF_CHKPT))
        psio_->close(PSIF_CHKPT, 1);
}

void myHF::find_occupation()
{
    // Don't mess with the occ, MOM's got it!
    if (MOM_started_) {
        MOM();
    } else {
        std::vector<std::pair<double, int> > pairs_a;
        std::vector<std::pair<double, int> > pairs_b;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs_a.push_back(make_pair(epsilon_a_->get(h, i), h));
        }
        for (int h=0; h<epsilon_b_->nirrep(); ++h) {
            for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
                pairs_b.push_back(make_pair(epsilon_b_->get(h, i), h));
        }
        sort(pairs_a.begin(),pairs_a.end());
        sort(pairs_b.begin(),pairs_b.end());

        if(!input_docc_ && !input_socc_){
            memset(nalphapi_, 0, sizeof(int) * epsilon_a_->nirrep());
            for (int i=0; i<nalpha_; ++i)
                nalphapi_[pairs_a[i].second]++;
        }
        if(!input_docc_ && !input_socc_){
            memset(nbetapi_, 0, sizeof(int) * epsilon_b_->nirrep());
            for (int i=0; i<nbeta_; ++i)
                nbetapi_[pairs_b[i].second]++;
        }

        int old_socc[8];
        int old_docc[8];
        for(int h = 0; h < nirrep_; ++h){
            old_socc[h] = soccpi_[h];
            old_docc[h] = doccpi_[h];
        }

        for (int h = 0; h < nirrep_; ++h) {
            soccpi_[h] = std::abs(nalphapi_[h] - nbetapi_[h]);
            doccpi_[h] = std::min(nalphapi_[h] , nbetapi_[h]);
        }

        bool occ_changed = false;
        for(int h = 0; h < nirrep_; ++h){
            if( old_socc[h] != soccpi_[h] || old_docc[h] != doccpi_[h]){
                occ_changed = true;
                break;
            }
        }

        // If print > 2 (diagnostics), print always
        if((print_ > 2 || (print_ && occ_changed)) && iteration_ > 0){
            if (WorldComm->me() == 0)
                fprintf(outfile, "\tOccupation by irrep:\n");
            print_occupation();
        }
        // Start MOM if needed (called here because we need the nocc
        // to be decided by Aufbau ordering prior to MOM_start)
        MOM_start();
    }
    // Do fractional orbital normalization here.
    frac();
}

void myHF::print_header()
{
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif

    if (WorldComm->me() == 0) {
        fprintf(outfile, "\n");
        fprintf(outfile, "         ---------------------------------------------------------\n");
        fprintf(outfile, "                                   SCF\n");
        fprintf(outfile, "            by Justin Turney, Rob Parrish, and Andy Simmonett\n");
        fprintf(outfile, "                             %4s Reference\n", options_.get_str("REFERENCE").c_str());
        fprintf(outfile, "                      %3d Threads, %6ld MiB Core\n", nthread, memory_ / 1000000L);
        fprintf(outfile, "         ---------------------------------------------------------\n");
        fprintf(outfile, "\n");
        fprintf(outfile, "  ==> Geometry <==\n\n");
    }

    molecule_->print();

    if (WorldComm->me() == 0) {
        fprintf(outfile, "  Running in %s symmetry.\n\n", molecule_->point_group()->symbol().c_str());

        fprintf(outfile, "  Nuclear repulsion = %20.15f\n\n", nuclearrep_);
        fprintf(outfile, "  Charge       = %d\n", charge_);
        fprintf(outfile, "  Multiplicity = %d\n", multiplicity_);
        fprintf(outfile, "  Electrons    = %d\n", nelectron_);
        fprintf(outfile, "  Nalpha       = %d\n", nalpha_);
        fprintf(outfile, "  Nbeta        = %d\n\n", nbeta_);

        fprintf(outfile, "  ==> Algorithm <==\n\n");
        fprintf(outfile, "  SCF Algorithm Type is %s.\n", options_.get_str("SCF_TYPE").c_str());
        fprintf(outfile, "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");
        if (MOM_excited_)
            fprintf(outfile, "  Excited-state MOM enabled.\n");
        else
            fprintf(outfile, "  MOM %s.\n", MOM_enabled_ ? "enabled" : "disabled");
        fprintf(outfile, "  Fractional occupation %s.\n", frac_enabled_ ? "enabled" : "disabled");
        fprintf(outfile, "  Guess Type is %s.\n", options_.get_str("GUESS").c_str());
        fprintf(outfile, "  Energy threshold   = %3.2e\n", energy_threshold_);
        fprintf(outfile, "  Density threshold  = %3.2e\n", density_threshold_);
        fprintf(outfile, "  Integral threshold = %3.2e\n\n", integral_threshold_);
        fflush(outfile);

        fprintf(outfile, "  ==> Primary Basis <==\n\n");
    }
    basisset_->print_by_level(outfile, print_);
    fflush(outfile);
}
void myHF::print_preiterations()
{
    CharacterTable ct = molecule_->point_group()->char_table();

    if (WorldComm->me() == 0) {
        fprintf(outfile, "   -------------------------------------------------------\n");
        fprintf(outfile, "    Irrep   Nso     Nmo     Nalpha   Nbeta   Ndocc  Nsocc\n");
        fprintf(outfile, "   -------------------------------------------------------\n");
        for (int h= 0; h < nirrep_; h++) {
            fprintf(outfile, "     %-3s   %6d  %6d  %6d  %6d  %6d  %6d\n", ct.gamma(h).symbol(), nsopi_[h], nmopi_[h], nalphapi_[h], nbetapi_[h], doccpi_[h], soccpi_[h]);
        }
        fprintf(outfile, "   -------------------------------------------------------\n");
        fprintf(outfile, "    Total  %6d  %6d  %6d  %6d  %6d  %6d\n", nso_, nmo_, nalpha_, nbeta_, nbeta_, nalpha_-nbeta_);
        fprintf(outfile, "   -------------------------------------------------------\n\n");
    }
}

void myHF::form_H()
{
    
    T_ = SharedMatrix(factory_->create_matrix(PSIF_SO_T));
    V_ = SharedMatrix(factory_->create_matrix(PSIF_SO_V));

    // Assumes these have already been created and stored
    T_->load(psio_, PSIF_OEI);
    V_->load(psio_, PSIF_OEI);

    if (debug_ > 2)
        T_->print(outfile);

    if (debug_ > 2)
        V_->print(outfile);

    if (perturb_h_) {
    
      if(perturb_ == embpot || perturb_ == sphere || perturb_ == dx) { // embedding potential read from file
        if(nirrep_ > 1)
          throw PSIEXCEPTION("RHF_embed: embedding, dx, and spherical potentials require 'symmetry c1'.");
        int nso = 0;
        for(int h=0; h < nirrep_; h++) nso += nsopi_[h];
        int nao = basisset_->nao();

        // Set up AO->SO transformation matrix (u)
        MintsHelper helper(options_, 0);
        SharedMatrix aotoso = helper.petite_list(true)->aotoso();
        int *col_offset = new int[nirrep_];
        col_offset[0] = 0;
        for(int h=1; h < nirrep_; h++)
          col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);

        double **u = block_matrix(nao, nso);
        for(int h=0; h < nirrep_; h++)
          for(int j=0; j < aotoso->coldim(h); j++)
            for(int i=0; i < nao; i++)
              u[i][j+col_offset[h]] = aotoso->get(h, i, j);
        delete[] col_offset;
        

        double *phi_ao, *phi_so, **V_eff;
        phi_ao = init_array(nao);
        phi_so = init_array(nso);
        V_eff = block_matrix(nso, nso);

        if(perturb_ == embpot) {

          FILE* input = fopen("EMBPOT", "r");
          int npoints;
          fscanf(input, "%d", &npoints);
          fprintf(outfile, "  npoints = %d\n", npoints);
          double x, y, z, w, v;
          double max = 0;
          for(int k=0; k < npoints; k++) {
            fscanf(input, "%lf %lf %lf %lf %lf", &x, &y, &z, &w, &v);
            if(fabs(v) > max) max = fabs(v);

            basisset_->compute_phi(phi_ao, x, y, z);
            // Transform phi_ao to SO basis
            C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1, 0.0, &(phi_so[0]), 1);
            for(int i=0; i < nso; i++)
              for(int j=0; j < nso; j++)                
                V_eff[i][j] += w * v * phi_so[i] * phi_so[j];
          } // npoints

          fprintf(outfile, "  Max. embpot value = %20.10f\n", max);
          fclose(input);

        } // embpot
        else if(perturb_ == dx) {
          dx_read(V_eff, phi_ao, phi_so, nao, nso, u);             

        } // dx file
        else if(perturb_ == sphere) {
          radius_ = options_.get_double("RADIUS");
          thickness_ = options_.get_double("THICKNESS");
          r_points_ = options_.get_int("R_POINTS");
          theta_points_ = options_.get_int("THETA_POINTS");
          phi_points_ = options_.get_int("PHI_POINTS");
          fprintf(outfile, "  Hard spherical potential radius         = %3.2f bohr\n", radius_);
          fprintf(outfile, "  Spherical potential thickness           = %3.2f bohr\n", thickness_);
          fprintf(outfile, "  Number of radial integration points     = %d\n", r_points_);
          fprintf(outfile, "  Number of colatitude integration points = %d\n", theta_points_);
          fprintf(outfile, "  Number of azimuthal integration points  = %d\n", phi_points_);

          double r_step = thickness_/r_points_; // bohr
          double theta_step = 2*pc_pi/theta_points_; // 1 degree in radians
          double phi_step = 2*pc_pi/phi_points_; // 1 degree in radians
          double weight = r_step * theta_step * phi_step;
          for(double r=radius_; r < radius_+thickness_; r += r_step) {
            for(double theta=0.0; theta < pc_pi; theta += theta_step) {  /* colatitude */
              for(double phi=0.0; phi < 2*pc_pi; phi += phi_step) { /* azimuthal */

                double x = r * sin(theta) * cos(phi);
                double y = r * sin(theta) * sin(phi);
                double z = r * cos(theta);

                double jacobian = weight * r * r * sin(theta);

                basisset_->compute_phi(phi_ao, x, y, z);

                C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1,
                        0.0, &(phi_so[0]), 1);

                for(int i=0; i < nso; i++)
                  for(int j=0; j < nso; j++)
                    V_eff[i][j] += jacobian * 1.0e6 * phi_so[i] * phi_so[j];
              }
            }
          }
        } // sphere

        if(WorldComm->me() == 0) {
          fprintf(outfile, "  Perturbing H by %f V_eff.\n", lambda_);
          if(options_.get_int("PRINT") > 3) mat_print(V_eff, nso, nso, outfile);
        }

        if(perturb_ == dx) {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, V_eff[i][j]); // ignore nuclear potential
        }
        else {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, (V_eff[i][j] + V_->get(i,j)));
        }

        free(phi_ao);
        free(phi_so);
        free_block(V_eff);
      }  // embpot or sphere
      else {
        OperatorSymmetry msymm(1, molecule_, integral_, factory_);
        vector<SharedMatrix> dipoles = msymm.create_matrices("Dipole");
        OneBodySOInt *so_dipole = integral_->so_dipole();
        so_dipole->compute(dipoles);

        if (perturb_ == dipole_x && (WorldComm->me() == 0)) {
            if (msymm.component_symmetry(0) != 0){
                fprintf(outfile, "  WARNING: You requested mu(x) perturbation, but mu(x) is not symmetric.\n");
            }
            else {
                if (WorldComm->me() == 0)
                    fprintf(outfile, "  Perturbing H by %f mu(x).\n", lambda_);
                dipoles[0]->scale(lambda_);
                V_->add(dipoles[0]);
            }
        } else if (perturb_ == dipole_y) {
            if (msymm.component_symmetry(1) != 0){
                if (WorldComm->me() == 0)
                    fprintf(outfile, "  WARNING: You requested mu(y) perturbation, but mu(y) is not symmetric.\n");
            }
            else {
                if (WorldComm->me() == 0)
                    fprintf(outfile, "  Perturbing H by %f mu(y).\n", lambda_);
                dipoles[1]->scale(lambda_);
                V_->add(dipoles[1]);
            }
        } else if (perturb_ == dipole_z) {
            if (msymm.component_symmetry(2) != 0){
                if (WorldComm->me() == 0)
                    fprintf(outfile, "  WARNING: You requested mu(z) perturbation, but mu(z) is not symmetric.\n");
            }
            else {
                if (WorldComm->me() == 0)
                    fprintf(outfile, "  Perturbing H by %f mu(z).\n", lambda_);
                dipoles[2]->scale(lambda_);
                V_->add(dipoles[2]);
            }
        }

      } // end dipole perturbations
    } // end perturb_h_

    
    // If an external field exists, add it to the one-electron Hamiltonian
    //boost::python::object pyExtern = dynamic_cast<PythonDataType*>(options_["EXTERN"].get())->to_python();
    //boost::shared_ptr<ExternalPotential> external = boost::python::extract<boost::shared_ptr<ExternalPotential> >(pyExtern);
    //cout<<"here!!"<<endl;
    /* if (external) {
        if (H_->nirrep() != 1)
            throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");
        if (print_) {
            external->set_print(print_);
            external->print();
        }
        SharedMatrix Vprime = external->computePotentialMatrix(basisset_);
        if (print_ > 3)
            Vprime->print();
        V_->add(Vprime);


        // Extra nuclear repulsion
        double enuc2 = external->computeNuclearEnergy(molecule_);
        if (print_) {
            fprintf(outfile, "  Old nuclear repulsion        = %20.15f\n", nuclearrep_);
            fprintf(outfile, "  Additional nuclear repulsion = %20.15f\n", enuc2);
            fprintf(outfile, "  Total nuclear repulsion      = %20.15f\n\n", nuclearrep_ + enuc2);
        }
        nuclearrep_ += enuc2;

    }  // end external */

    // Save perturbed V_ for future (e.g. correlated) calcs
    V_->save(psio_, PSIF_OEI);

    H_->copy(T_);
    H_->add(V_);

    if (print_ > 3)
        H_->print(outfile);
}

void myHF::form_Shalf()
{
    // ==> SYMMETRIC ORTHOGONALIZATION <== //

    // S_ is computed by wavefunction

    SharedMatrix eigvec= factory_->create_shared_matrix("L");
    SharedMatrix eigtemp= factory_->create_shared_matrix("Temp");
    SharedMatrix eigtemp2= factory_->create_shared_matrix("Temp2");
    SharedMatrix eigvec_store= factory_->create_shared_matrix("TempStore");
    SharedVector eigval(factory_->create_vector());
    SharedVector eigval_store(factory_->create_vector());

    //Used to do this 3 times, now only once
    S_->diagonalize(eigvec, eigval);
    eigvec_store->copy(eigvec);
    eigval_store->copy(eigval.get());

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    int *dimpi = eigval->dimpi();
    double min_S = fabs(eigval->get(0,0));
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi[h]; ++i) {
            if (min_S > eigval->get(h,i))
                min_S = eigval->get(h,i);
            double scale = 1.0 / sqrt(eigval->get(h, i));
            eigval->set(h, i, scale);
        }
    }
    if (print_ && (WorldComm->me() == 0))
        fprintf(outfile,"  Minimum eigenvalue in the overlap matrix is %14.10E.\n",min_S);
    // Create a vector matrix from the converted eigenvalues
    eigtemp2->set_diagonal(eigval);

    eigtemp->gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    X_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);

    // ==> CANONICAL ORTHOGONALIZATION <== //

    // Decide symmetric or canonical
    double S_cutoff = options_.get_double("S_TOLERANCE");
    if (min_S > S_cutoff && options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC") {

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile,"  Using Symmetric Orthogonalization.\n");

    } else {

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile,"  Using Canonical Orthogonalization with cutoff of %14.10E.\n",S_cutoff);

        //Diagonalize S (or just get a fresh copy)
        eigvec->copy(eigvec_store.get());
        eigval->copy(eigval_store.get());
        int delta_mos = 0;
        for (int h=0; h<nirrep_; ++h) {
            //in each irrep, scale significant cols i  by 1.0/sqrt(s_i)
            int start_index = 0;
            for (int i=0; i<dimpi[h]; ++i) {
                if (S_cutoff  < eigval->get(h,i)) {
                    double scale = 1.0 / sqrt(eigval->get(h, i));
                    eigvec->scale_column(h, i, scale);
                } else {
                    start_index++;
                    nmopi_[h]--;
                    nmo_--;
                }
            }
            if (print_>2 && (WorldComm->me() == 0))
                fprintf(outfile,"  Irrep %d, %d of %d possible MOs eliminated.\n",h,start_index,nsopi_[h]);

            delta_mos += start_index;
        }

        X_->init(nirrep_,nsopi_,nmopi_,"X (Canonical Orthogonalization)");
        for (int h=0; h<eigval->nirrep(); ++h) {
            //Copy significant columns of eigvec into X in
            //descending order
            int start_index = 0;
            for (int i=0; i<dimpi[h]; ++i) {
                if (S_cutoff  < eigval->get(h,i)) {
                } else {
                    start_index++;
                }
            }
            for (int i=0; i<dimpi[h]-start_index; ++i) {
                for (int m = 0; m < dimpi[h]; m++)
                    X_->set(h,m,i,eigvec->get(h,m,dimpi[h]-i-1));
            }
        }

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile,"  Overall, %d of %d possible MOs eliminated.\n\n",delta_mos,nso_);

        // Refreshes twice in RHF, no big deal
        epsilon_a_->init(nmopi_);
        Ca_->init(nirrep_,nsopi_,nmopi_,"MO coefficients");
        epsilon_b_->init(nmopi_);
        Cb_->init(nirrep_,nsopi_,nmopi_,"MO coefficients");
    }

    // Temporary variables needed by diagonalize_F
    diag_temp_   = SharedMatrix(new Matrix(nirrep_, nmopi_, nsopi_));
    diag_F_temp_ = SharedMatrix(new Matrix(nirrep_, nmopi_, nmopi_));
    diag_C_temp_ = SharedMatrix(new Matrix(nirrep_, nmopi_, nmopi_));

    if (print_ > 3) {
        S_->print(outfile);
        X_->print(outfile);
    }
    fflush(outfile);
}

void myHF::compute_fcpi()
{
    // FROZEN_DOCC takes precedence, FREEZE_CORE directive has second priority
    if (options_["FROZEN_DOCC"].has_changed()) {
        if (options_["FROZEN_DOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_DOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_integer();
        }
    } else {

        int nfzc = 0;
        if (options_.get_int("NUM_FROZEN_DOCC") != 0) {
            nfzc = options_.get_int("NUM_FROZEN_DOCC");
        } else {
            nfzc = molecule_->nfrozen_core(options_.get_str("FREEZE_CORE"));
        }
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
            frzcpi_[h] = 0;
        }
        sort(pairs.begin(),pairs.end());

        for (int i=0; i<nfzc; ++i)
            frzcpi_[pairs[i].second]++;
    }
}

void myHF::compute_fvpi()
{
    // FROZEN_UOCC takes precedence, FREEZE_UOCC directive has second priority
    if (options_["FROZEN_UOCC"].has_changed()) {
        if (options_["FROZEN_UOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_UOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_integer();
        }
    } else {
        int nfzv = options_.get_int("NUM_FROZEN_UOCC");
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
            frzvpi_[h] = 0;
        }
        sort(pairs.begin(),pairs.end(), greater<std::pair<double, int> >());

        for (int i=0; i<nfzv; ++i)
            frzvpi_[pairs[i].second]++;
    }
}

void myHF::print_orbitals(const char* header, std::vector<std::pair<double, std::pair<const char*, int> > > orbs)
{
    if (WorldComm->me() == 0) {
        fprintf(outfile, "\t%-70s\n\n\t", header);
        int count = 0;
        for (int i = 0; i < orbs.size(); i++) {
            fprintf(outfile, "%4d%-4s%11.6f  ", orbs[i].second.second, orbs[i].second.first, orbs[i].first);
            if (count++ % 3 == 2 && count != orbs.size())
                fprintf(outfile, "\n\t");
        }
        fprintf(outfile, "\n\n");
    }
}

void myHF::print_orbitals()
{
    char **labels = molecule_->irrep_labels();

    if (WorldComm->me() == 0)
        fprintf(outfile, "\tOrbital Energies (a.u.)\n\t-----------------------\n\n");

    std::string reference = options_.get_str("REFERENCE");
    if((reference == "RHF") || (reference == "RKS")){

        std::vector<std::pair<double, std::pair<const char*, int> > > occ;
        std::vector<std::pair<double, std::pair<const char*, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_e.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occ.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                vir.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));

        }
        std::sort(occ.begin(), occ.end());
        std::sort(vir.begin(), vir.end());

        print_orbitals("Doubly Occupied:", occ);
        print_orbitals("Virtual:", vir);

    }else if((reference == "UHF") || (reference == "UKS") ||
        (reference == "CUHF")){

        std::vector<std::pair<double, std::pair<const char*, int> > > occA;
        std::vector<std::pair<double, std::pair<const char*, int> > > virA;
        std::vector<std::pair<double, std::pair<const char*, int> > > occB;
        std::vector<std::pair<double, std::pair<const char*, int> > > virB;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_eA;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_eA.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_eA.begin(), orb_eA.end());

            std::vector<int> orb_orderA(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderA[orb_eA[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occA.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_orderA[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                virA.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_orderA[a] + 1)));

            std::vector<std::pair<double, int> > orb_eB;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_eB.push_back(make_pair(epsilon_b_->get(h,a), a));
            std::sort(orb_eB.begin(), orb_eB.end());

            std::vector<int> orb_orderB(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderB[orb_eB[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                occB.push_back(make_pair(epsilon_b_->get(h,a), make_pair(labels[h],orb_orderB[a] + 1)));
            for (int a = nbetapi_[h]; a < nmopi_[h]; a++)
                virB.push_back(make_pair(epsilon_b_->get(h,a), make_pair(labels[h],orb_orderB[a] + 1)));

        }
        std::sort(occA.begin(), occA.end());
        std::sort(virA.begin(), virA.end());
        std::sort(occB.begin(), occB.end());
        std::sort(virB.begin(), virB.end());

        print_orbitals("Alpha Occupied:", occA);
        print_orbitals("Alpha Virtual:", virA);
        print_orbitals("Beta Occupied:", occB);
        print_orbitals("Beta Virtual:", virB);

    }else if(reference == "ROHF"){

        std::vector<std::pair<double, std::pair<const char*, int> > > docc;
        std::vector<std::pair<double, std::pair<const char*, int> > > socc;
        std::vector<std::pair<double, std::pair<const char*, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_e.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                docc.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nbetapi_[h] ; a < nalphapi_[h]; a++)
                socc.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h] ; a < nmopi_[h]; a++)
                vir.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));

        }
        std::sort(docc.begin(), docc.end());
        std::sort(socc.begin(), socc.end());
        std::sort(vir.begin(), vir.end());

        print_orbitals("Doubly Occupied:", docc);
        print_orbitals("Singly Occupied:", socc);
        print_orbitals("Virtual:", vir);

    }else{
        throw PSIEXCEPTION("Unknown reference in myHF::print_orbitals");
    }

    for(int h = 0; h < nirrep_; ++h)
        free(labels[h]);
    free(labels);

    if (WorldComm->me() == 0)
        fprintf(outfile, "\tFinal Occupation by Irrep:\n");
    print_occupation();
}

void myHF::guess()
{
    // don't save guess energy as "the" energy because we need to avoid
    // a false positive test for convergence on the first iteration (that
    // was happening before in tests/scf-guess-read before I removed
    // the statements putting this into E_).  -CDS 3/25/13
    double guess_E; 

    //What does the user want?
    //Options will be:
    // "READ"-try to read MOs from guess file, projecting if needed
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties
    string guess_type = options_.get_str("GUESS");
    if (guess_type == "READ" && !psio_->exists(PSIF_SCF_MOS)) {
        fprintf(outfile, "  SCF Guess was Projection but file not found.\n");
        fprintf(outfile, "  Switching over to SAD guess.\n\n");
        guess_type = "SAD";
    }

    if (guess_type == "READ") {

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile, "  SCF Guess: Projection.\n\n");

        load_orbitals(); // won't save the energy from here
        form_D();

    } else if (guess_type == "SAD") {

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile, "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF.\n\n");

        //Superposition of Atomic Density (RHF only at present)
        compute_SAD_guess();
        guess_E = compute_initial_E();

    } else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile, "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");

        Fa_->zero(); //Try Fa_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        int *opi = S_->rowspi();
        int nirreps = S_->nirrep();
        for (h=0; h<nirreps; ++h) {
            for (i=0; i<opi[h]; ++i) {
                Fa_->set(h,i,i,H_->get(h,i,i));
                for (j=0; j<i; ++j) {
                    Fa_->set(h,i,j,0.875*S_->get(h,i,j)*(H_->get(h,i,i)+H_->get(h,j,j)));
                    Fa_->set(h,j,i,Fa_->get(h,i,j));
                }
            }
        }
        Fb_->copy(Fa_);
        form_initial_C();
        find_occupation();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "CORE") {

        if (print_ && (WorldComm->me() == 0))
            fprintf(outfile, "  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        Fa_->copy(H_); //Try the core Hamiltonian as the Fock Matrix
        Fb_->copy(H_);

        form_initial_C();
        find_occupation();
        form_D();
        guess_E = compute_initial_E();
    }
    if (print_ > 3) {
        Ca_->print();
        Cb_->print();
        Da_->print();
        Db_->print();
        Fa_->print();
        Fb_->print();
    }

    // This is confusing the user and valgrind.
    //if (print_ && (WorldComm->me() == 0))
    //    fprintf(outfile, "  Guess energy: %20.14f\n\n", guess_E);

    E_ = 0.0; // don't use this guess in our convergence checks
}

void myHF::save_orbitals()
{
    psio_->open(PSIF_SCF_MOS,PSIO_OPEN_NEW);

    if (print_ && (WorldComm->me() == 0))
        fprintf(outfile,"\n  Saving occupied orbitals to File %d.\n", PSIF_SCF_MOS);

    psio_->write_entry(PSIF_SCF_MOS,"SCF ENERGY",(char *) &(E_),sizeof(double));
    psio_->write_entry(PSIF_SCF_MOS,"NIRREP",(char *) &(nirrep_),sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NSOPI",(char *) &(nsopi_[0]),nirrep_*sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NALPHAPI",(char *) &(nalphapi_[0]),nirrep_*sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NBETAPI",(char *) &(nbetapi_[0]),nirrep_*sizeof(int));

    char *basisname = strdup(options_.get_str("BASIS").c_str());
    int basislength = strlen(options_.get_str("BASIS").c_str()) + 1;

    psio_->write_entry(PSIF_SCF_MOS,"BASIS NAME LENGTH",(char *)(&basislength),sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"BASIS NAME",basisname,basislength*sizeof(char));

    // upon loading, need to know what value of puream was used
    int old_puream = (basisset_->has_puream() ? 1 : 0);
    psio_->write_entry(PSIF_SCF_MOS,"PUREAM",(char *)(&old_puream),sizeof(int));

    SharedMatrix Ctemp_a(new Matrix("ALPHA MOS", nirrep_, nsopi_, nalphapi_));
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nalphapi_[h]; i++)
                Ctemp_a->set(h,m,i,Ca_->get(h,m,i));
    Ctemp_a->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);

    SharedMatrix Ctemp_b(new Matrix("BETA MOS", nirrep_, nsopi_, nbetapi_));
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nbetapi_[h]; i++)
                Ctemp_b->set(h,m,i,Cb_->get(h,m,i));
    Ctemp_b->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);

    psio_->close(PSIF_SCF_MOS,1);
    free(basisname);
}

void myHF::load_orbitals()
{
    psio_->open(PSIF_SCF_MOS,PSIO_OPEN_OLD);

    int basislength, old_puream;
    psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME LENGTH",
        (char *)(&basislength),sizeof(int));
    char *basisnamec = new char[basislength];
    psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME",basisnamec,
        basislength*sizeof(char));
    psio_->read_entry(PSIF_SCF_MOS,"PUREAM",(char *)(&old_puream),
        sizeof(int));
    bool old_forced_puream = (old_puream) ? true : false;
    std::string basisname(basisnamec);

    if (basisname == "")
        throw PSIEXCEPTION("SCF::load_orbitals: Custom basis sets not allowed for projection from a previous SCF");

    if (print_) {
        if (basisname != options_.get_str("BASIS")) {
            if (WorldComm->me() == 0) {
                fprintf(outfile,"  Computing basis set projection from %s to %s.\n", \
                    basisname.c_str(),options_.get_str("BASIS").c_str());
            }
        } else {
            if (WorldComm->me() == 0)
                fprintf(outfile,"  Using orbitals from previous SCF, no projection.\n");
        }
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(old_forced_puream));
    molecule_->set_basis_all_atoms(basisname, "DUAL_BASIS_SCF");
    boost::shared_ptr<BasisSet> dual_basis = BasisSet::construct(parser, molecule_, "DUAL_BASIS_SCF");

    psio_->read_entry(PSIF_SCF_MOS,"SCF ENERGY",(char *) &(E_),sizeof(double));

    int old_nirrep, old_nsopi[8];
    psio_->read_entry(PSIF_SCF_MOS,"NIRREP",(char *) &(old_nirrep),sizeof(int));

    if (old_nirrep != nirrep_)
        throw PSIEXCEPTION("SCF::load_orbitals: Projection of orbitals between different symmetries is not currently supported");

    psio_->read_entry(PSIF_SCF_MOS,"NSOPI",(char *) (old_nsopi),nirrep_*sizeof(int));
    psio_->read_entry(PSIF_SCF_MOS,"NALPHAPI",(char *) &(nalphapi_[0]),nirrep_*sizeof(int));
    psio_->read_entry(PSIF_SCF_MOS,"NBETAPI",(char *) &(nbetapi_[0]),nirrep_*sizeof(int));

    for (int h = 0; h < nirrep_; h++) {
        doccpi_[h] = nbetapi_[h];
        soccpi_[h] = nalphapi_[h] - nbetapi_[h];
    }

    SharedMatrix Ctemp_a(new Matrix("ALPHA MOS", nirrep_, old_nsopi, nalphapi_));
    Ctemp_a->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
    SharedMatrix Ca;
    if (basisname != options_.get_str("BASIS")) {
        Ca = BasisProjection(Ctemp_a, nalphapi_, dual_basis, basisset_);
    } else {
        Ca = Ctemp_a;
    }
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nalphapi_[h]; i++)
                Ca_->set(h,m,i,Ca->get(h,m,i));

    SharedMatrix Ctemp_b(new Matrix("BETA MOS", nirrep_, old_nsopi, nbetapi_));
    Ctemp_b->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
    SharedMatrix Cb;
    if (basisname != options_.get_str("BASIS")) {
        Cb = BasisProjection(Ctemp_b, nbetapi_, dual_basis, basisset_);
    } else {
        Cb = Ctemp_b;
    }
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nbetapi_[h]; i++)
                Cb_->set(h,m,i,Cb->get(h,m,i));

    psio_->close(PSIF_SCF_MOS,1);
    delete[] basisnamec;
}


void myHF::check_phases()
{
    for (int h=0; h<nirrep_; ++h) {
        for (int p = 0; p < Ca_->colspi(h); ++p) {
            for (int mu = 0; mu < Ca_->rowspi(h); ++mu) {
                if (fabs(Ca_->get(h, mu, p)) > 1.0E-3) {
                    if (Ca_->get(h, mu, p) < 1.0E-3) {
                        Ca_->scale_column(h, p, -1.0);
                    }
                    break;
                }
            }
        }
    }

    if (Ca_ != Cb_) {
        for (int h=0; h<nirrep_; ++h) {
            for (int p = 0; p < Cb_->colspi(h); ++p) {
                for (int mu = 0; mu < Cb_->rowspi(h); ++mu) {
                    if (fabs(Cb_->get(h, mu, p)) > 1.0E-3) {
                        if (Cb_->get(h, mu, p) < 1.0E-3) {
                            Cb_->scale_column(h, p, -1.0);
                        }
                        break;
                    }
                }
            }
        }
    }
}

void myHF::dump_to_checkpoint()
{
    if(!psio_->open_check(PSIF_CHKPT))
        psio_->open(PSIF_CHKPT, PSIO_OPEN_OLD);
    chkpt_->wt_nirreps(nirrep_);
    char **labels = molecule_->irrep_labels();
    chkpt_->wt_irr_labs(labels);
    for(int h = 0; h < nirrep_; ++h)
        free(labels[h]);
    free(labels);
    chkpt_->wt_nmo(nmo_);
    chkpt_->wt_nso(nso_);
    chkpt_->wt_nao(basisset_->nao());
    chkpt_->wt_ref(0);
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_enuc(nuclearrep_);
    chkpt_->wt_orbspi(nmopi_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_openpi(soccpi_);
    chkpt_->wt_phase_check(0);
    chkpt_->wt_sopi(nsopi_);
    // Figure out total number of frozen docc/uocc orbitals
    int nfzc = 0;
    int nfzv = 0;
    for (int h = 0; h < nirrep_; h++) {
        nfzc += frzcpi_[h];
        nfzv += frzvpi_[h];
    }
    chkpt_->wt_nfzc(nfzc);
    chkpt_->wt_nfzv(nfzv);
    // These were computed by myHF::finalize()
    chkpt_->wt_frzcpi(frzcpi_);
    chkpt_->wt_frzvpi(frzvpi_);

    int m = 0;
    for(int h = 0; h < nirrep_; ++h)
        if(soccpi_[h]) ++m;
    chkpt_->wt_iopen(m*(m+1)/2);

    if(options_.get_str("REFERENCE") == "UHF" ||
        options_.get_str("REFERENCE") == "CUHF"){

        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_alpha_evals(values);
        delete[] values;
        values = epsilon_b_->to_block_vector();
        chkpt_->wt_beta_evals(values);
        delete[] values;
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_alpha_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
        vectors = Cb_->to_block_matrix();
        chkpt_->wt_beta_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
    }else{
        // All other reference type yield restricted orbitals
        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_evals(values);
        delete[] values;
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
        double *ftmp = Fa_->to_lower_triangle();
        chkpt_->wt_fock(ftmp);
        delete[] ftmp;
    }
    psio_->close(PSIF_CHKPT, 1);
}

double myHF::compute_energy()
{
    std::string reference = options_.get_str("REFERENCE");

    bool converged = false;
    MOM_performed_ = false;
    diis_performed_ = false;
    // Neither of these are idempotent
    if (options_.get_str("GUESS") == "SAD" || options_.get_str("GUESS") == "READ")
        iteration_ = -1;
    else
        iteration_ = 0;

    if (print_ && (WorldComm->me() == 0))
        fprintf(outfile, "  ==> Pre-Iterations <==\n\n");

    if (print_)
        print_preiterations();

    // Andy trick 2.0
    std::string old_scf_type = options_.get_str("SCF_TYPE");
    if (options_.get_bool("DF_SCF_GUESS") && !(old_scf_type == "DF" || old_scf_type == "CD")) {
         fprintf(outfile, "  Starting with a DF guess...\n\n");
         if(!options_["DF_BASIS_SCF"].has_changed()) {
             // TODO: Match Dunning basis sets 
             molecule_->set_basis_all_atoms("CC-PVDZ-JKFIT", "DF_BASIS_SCF");
         }
         scf_type_ = "DF";
         options_.set_str("SCF","SCF_TYPE","DF"); // Scope is reset in proc.py. This is not pretty, but it works
    }

    if(attempt_number_ == 1){
        boost::shared_ptr<MintsHelper> mints (new MintsHelper(options_, 0));
        mints->one_electron_integrals();
        

        integrals();
        

        timer_on("Form H");
        
        form_H(); //Core Hamiltonian
        
        timer_off("Form H");
        

        timer_on("Form S/X");
        form_Shalf(); //S and X Matrix
        timer_off("Form S/X");

        timer_on("Guess");
        guess(); // Guess
        timer_off("Guess");

    }else{
        // We're reading the orbitals from the previous set of iterations.
        form_D();
        E_ = compute_initial_E();
    }

    bool df = (options_.get_str("SCF_TYPE") == "DF");

    if (WorldComm->me() == 0) {
        fprintf(outfile, "  ==> Iterations <==\n\n");
        fprintf(outfile, "%s                        Total Energy        Delta E     RMS |[F,P]|\n\n", df ? "   " : "");
    }
    fflush(outfile);

    // SCF iterations
    do {
        iteration_++;

        save_density_and_energy();

        // Call any preiteration callbacks
        call_preiteration_callbacks();

        timer_on("Form G");
        form_G();
        timer_off("Form G");

        // Reset fractional SAD occupation
        if (iteration_ == 0 && options_.get_str("GUESS") == "SAD")
            reset_SAD_occupation();

        timer_on("Form F");
        form_F();
        timer_off("Form F");

        if (print_>3) {
            Fa_->print(outfile);
            Fb_->print(outfile);
        }

	//bool has_efp = options.get("HAS_EFP");

	// XXX
	//if (has_efp) {
	//efp_get_multipole_count
	// allocate arrays
	//efp_get_multipoles
	// compute 1e contributions
	//}

        E_ = compute_E();

	// XXX
	//if (has_efp) {
	//double efp_energy;
	//efp_scf_update(efp, &efp_energy);
	//E_ += efp_energy;
	//}

        timer_on("DIIS");
        bool add_to_diis_subspace = false;
        if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
            add_to_diis_subspace = true;

        compute_orbital_gradient(add_to_diis_subspace);

        if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
            diis_performed_ = diis();
        } else {
            diis_performed_ = false;
        }
        timer_off("DIIS");

        if (print_>4 && diis_performed_ && (WorldComm->me() == 0)) {
            fprintf(outfile,"  After DIIS:\n");
            Fa_->print(outfile);
            Fb_->print(outfile);
        }

        // If we're too well converged, or damping wasn't enabled, do DIIS
        damping_performed_ = (damping_enabled_ && iteration_ > 1 && Drms_ > damping_convergence_);

        std::string status = "";
        if(diis_performed_){
            if(status != "") status += "/";
            status += "DIIS";
        }
        if(MOM_performed_){
            if(status != "") status += "/";
            status += "MOM";
        }
        if(damping_performed_){
            if(status != "") status += "/";
            status += "DAMP";
        }
        if(frac_performed_){
            if(status != "") status += "/";
            status += "FRAC";
        }



        timer_on("Form C");
        form_C();
        timer_off("Form C");
        timer_on("Form D");
        form_D();
        timer_off("Form D");

        Process::environment.globals["SCF ITERATION ENERGY"] = E_;

        // After we've built the new D, damp the update if
        if(damping_performed_) damp_update();

        if (print_ > 3){
            Ca_->print(outfile);
            Cb_->print(outfile);
            Da_->print(outfile);
            Db_->print(outfile);
        }

        converged = test_convergency();

        df = (options_.get_str("SCF_TYPE") == "DF");

        if (WorldComm->me() == 0) {
            fprintf(outfile, "   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n", df ? "DF-" : "",
                              reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, status.c_str());
            fflush(outfile);
        }

        // If a an excited MOM is requested but not started, don't stop yet
        if (MOM_excited_ && !MOM_started_) converged = false;

        // If a fractional occupation is requested but not started, don't stop yet
        if (frac_enabled_ && !frac_performed_) converged = false;

        // If a DF Guess environment, reset the JK object, and keep running
        if (converged && options_.get_bool("DF_SCF_GUESS") && !(old_scf_type == "DF" || old_scf_type == "CD")) {
            fprintf(outfile, "\n  DF guess converged.\n\n"); // Be cool dude. 
            converged = false;
            if(initialized_diis_manager_)
                diis_manager_->reset_subspace();
            scf_type_ = old_scf_type;
            options_.set_str("SCF","SCF_TYPE",old_scf_type);
            old_scf_type = "DF";
            integrals();
        }

        // Call any postiteration callbacks
        call_postiteration_callbacks();

    } while (!converged && iteration_ < maxiter_ );

    if (WorldComm->me() == 0)
        fprintf(outfile, "\n  ==> Post-Iterations <==\n\n");

    check_phases();
    compute_spin_contamination();
    frac_renormalize();

    if (converged || !fail_on_maxiter_) {
        // Need to recompute the Fock matrices, as they are modified during the SCF interation
        // and might need to be dumped to checkpoint later
        form_F();

        // Print the orbitals
        if(print_)
            print_orbitals();

        if (WorldComm->me() == 0 && converged) {
            fprintf(outfile, "  Energy converged.\n\n");
        }
        if (WorldComm->me() == 0 && !converged) {
            fprintf(outfile, "  Energy did not converge, but proceeding anyway.\n\n");
        }
        if (WorldComm->me() == 0) {
            fprintf(outfile, "  @%s%s Final Energy: %20.14f", df ? "DF-" : "", reference.c_str(), E_);
            if (perturb_h_) {
                fprintf(outfile, " with %f perturbation", lambda_);
            }
            fprintf(outfile, "\n\n");
            print_energies();
        }

        // Properties
        if (print_) {
            boost::shared_ptr<OEProp> oe(new OEProp());
            oe->set_title("SCF");
            oe->add("DIPOLE");

            if (print_ >= 2) {
                oe->add("QUADRUPOLE");
                oe->add("MULLIKEN_CHARGES");
            }

            if (print_ >= 3) {
                oe->add("LOWDIN_CHARGES");
                oe->add("MAYER_INDICES");
                oe->add("WIBERG_LOWDIN_INDICES");
            }

            if (WorldComm->me() == 0)
                fprintf(outfile, "  ==> Properties <==\n\n");
            oe->compute();

//  Comments so that autodoc utility will find these PSI variables
//
//  Process::environment.globals["SCF DIPOLE X"] =
//  Process::environment.globals["SCF DIPOLE Y"] =
//  Process::environment.globals["SCF DIPOLE Z"] =
//  Process::environment.globals["SCF QUADRUPOLE XX"] =
//  Process::environment.globals["SCF QUADRUPOLE XY"] =
//  Process::environment.globals["SCF QUADRUPOLE XZ"] =
//  Process::environment.globals["SCF QUADRUPOLE YY"] =
//  Process::environment.globals["SCF QUADRUPOLE YZ"] =
//  Process::environment.globals["SCF QUADRUPOLE ZZ"] =

        }

        save_information();
    } else {
        if (WorldComm->me() == 0) {
            fprintf(outfile, "  Failed to converged.\n");
            fprintf(outfile, "    NOTE: MO Coefficients will not be saved to Checkpoint.\n");
        }
        E_ = 0.0;
        if(psio_->open_check(PSIF_CHKPT))
            psio_->close(PSIF_CHKPT, 1);

        // Throw if we didn't converge?
        die_if_not_converged();
    }

    // Orbitals are always saved, in case an MO guess is requested later
    save_orbitals();
    if (options_.get_str("SAPT") != "FALSE") //not a bool because it has types
        save_sapt_info();

    WorldComm->sync();

    // Perform wavefunction stability analysis
    if(options_.get_str("STABILITY_ANALYSIS") != "NONE")
        stability_analysis();

    // Clean memory off, handle diis closeout, etc
    //finalize();


    //fprintf(outfile,"\nComputation Completed\n");
    fflush(outfile);
    return E_;
}

void myHF::print_energies()
{
    fprintf(outfile, "   => Energetics <=\n\n");
    fprintf(outfile, "    Nuclear Repulsion Energy =        %24.16f\n", energies_["Nuclear"]);
    fprintf(outfile, "    One-Electron Energy =             %24.16f\n", energies_["One-Electron"]);
    fprintf(outfile, "    Two-Electron Energy =             %24.16f\n", energies_["Two-Electron"]);
    fprintf(outfile, "    DFT Exchange-Correlation Energy = %24.16f\n", energies_["XC"]); 
    fprintf(outfile, "    Empirical Dispersion Energy =     %24.16f\n", energies_["-D"]);
    fprintf(outfile, "    Total Energy =                    %24.16f\n", energies_["Nuclear"] + 
        energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"] + energies_["-D"]); 
    fprintf(outfile, "\n");
    
    Process::environment.globals["NUCLEAR REPULSION ENERGY"] = energies_["Nuclear"];
    Process::environment.globals["ONE-ELECTRON ENERGY"] = energies_["One-Electron"];
    Process::environment.globals["TWO-ELECTRON ENERGY"] = energies_["Two-Electron"];
    if (fabs(energies_["XC"]) > 1.0e-14) {
        Process::environment.globals["DFT XC ENERGY"] = energies_["XC"];
        Process::environment.globals["DFT FUNCTIONAL TOTAL ENERGY"] = energies_["Nuclear"] + 
            energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"];
        Process::environment.globals["DFT TOTAL ENERGY"] = energies_["Nuclear"] + 
            energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"] + energies_["-D"];
    } else {
        Process::environment.globals["myHF TOTAL ENERGY"] = energies_["Nuclear"] + 
            energies_["One-Electron"] + energies_["Two-Electron"];
    }
    if (fabs(energies_["-D"]) > 1.0e-14) {
        Process::environment.globals["DISPERSION CORRECTION ENERGY"] = energies_["-D"];
    }
//  Comment so that autodoc utility will find this PSI variable
//     It doesn't really belong here but needs to be linked somewhere
//  Process::environment.globals["DOUBLE-HYBRID CORRECTION ENERGY"]
}

void myHF::print_occupation()
{
    if (WorldComm->me() == 0) {
        char **labels = molecule_->irrep_labels();
        std::string reference = options_.get_str("REFERENCE");
        fprintf(outfile, "\t      ");
        for(int h = 0; h < nirrep_; ++h) fprintf(outfile, " %4s ", labels[h]); fprintf(outfile, "\n");
        fprintf(outfile, "\tDOCC [ ");
        for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", doccpi_[h]);
        fprintf(outfile, " %4d ]\n", doccpi_[nirrep_-1]);
        if(reference != "RHF" && reference != "RKS"){
            fprintf(outfile, "\tSOCC [ ");
            for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", soccpi_[h]);
            fprintf(outfile, " %4d ]\n", soccpi_[nirrep_-1]);
        }
        if (MOM_excited_) {
            // Also print nalpha and nbeta per irrep, which are more physically meaningful
            fprintf(outfile, "\tNA   [ ");
            for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", nalphapi_[h]);
            fprintf(outfile, " %4d ]\n", nalphapi_[nirrep_-1]);
            fprintf(outfile, "\tNB   [ ");
            for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", nbetapi_[h]);
            fprintf(outfile, " %4d ]\n", nbetapi_[nirrep_-1]);
        }

        for(int h = 0; h < nirrep_; ++h) free(labels[h]); free(labels);
        fprintf(outfile,"\n");
    }
}

void myHF::diagonalize_F(const SharedMatrix& Fm, SharedMatrix& Cm, boost::shared_ptr<Vector>& epsm)
{
    //Form F' = X'FX for canonical orthogonalization
    diag_temp_->gemm(true, false, 1.0, X_, Fm, 0.0);
    diag_F_temp_->gemm(false, false, 1.0, diag_temp_, X_, 0.0);

    //Form C' = eig(F')
    diag_F_temp_->diagonalize(diag_C_temp_, epsm);

    //Form C = XC'
    Cm->gemm(false, false, 1.0, X_, diag_C_temp_, 0.0);
}

void myHF::reset_SAD_occupation()
{
    // RHF style for now
    for (int h = 0; h < Da_->nirrep(); h++) {
        nalphapi_[h] = sad_nocc_[h];
        nbetapi_[h]  = sad_nocc_[h];
        doccpi_[h]   = sad_nocc_[h];
        soccpi_[h]   = 0;
    }
}
SharedMatrix myHF::form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi)
{
    int* nsopi = Cso->rowspi();
    int* nmopi = Cso->colspi();
    int* nvirpi = new int[nirrep_];

    for (int h = 0; h < nirrep_; h++)
        nvirpi[h] = nmopi[h] - noccpi[h];

    SharedMatrix Fia(new Matrix("Fia (Some Basis)", nirrep_, noccpi, nvirpi));

    // Hack to get orbital e for this Fock
    SharedMatrix C2(new Matrix("C2", nirrep_, nsopi, nmopi));
    boost::shared_ptr<Vector> E2(new Vector("E2", nirrep_, nmopi));
    diagonalize_F(Fso, C2, E2);

    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi[h];
        int nso = nsopi[h];
        int nvir = nvirpi[h];
        int nocc = noccpi[h];

        if (nmo == 0 || nso == 0 || nvir == 0 || nocc == 0) continue;

        //double** C = Cso->pointer(h);
        double** C = C2->pointer(h);
        double** F = Fso->pointer(h);
        double** Fiap = Fia->pointer(h);

        double** Temp = block_matrix(nocc, nso);

        C_DGEMM('T','N',nocc,nso,nso,1.0,C[0],nmo,F[0],nso,0.0,Temp[0],nso);
        C_DGEMM('N','N',nocc,nvir,nso,1.0,Temp[0],nso,&C[0][nocc],nmo,0.0,Fiap[0],nvir);

        free_block(Temp);

        //double* eps = E2->pointer(h);
        //for (int i = 0; i < nocc; i++)
        //    for (int a = 0; a < nvir; a++)
        //        Fiap[i][a] /= eps[a + nocc] - eps[i];

    }

    //Fia->print();

    delete[] nvirpi;

    return Fia;
}
SharedMatrix myHF::form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso)
{
    SharedMatrix FDSmSDF(new Matrix("FDS-SDF", nirrep_, nsopi_, nsopi_));
    SharedMatrix DS(new Matrix("DS", nirrep_, nsopi_, nsopi_));

    DS->gemm(false,false,1.0,Dso,S_,0.0);
    FDSmSDF->gemm(false,false,1.0,Fso,DS,0.0);

    SharedMatrix SDF(FDSmSDF->transpose());
    FDSmSDF->subtract(SDF);

    DS.reset();
    SDF.reset();

    SharedMatrix XP(new Matrix("X'(FDS - SDF)", nirrep_, nmopi_, nsopi_));
    SharedMatrix XPX(new Matrix("X'(FDS - SDF)X", nirrep_, nmopi_, nmopi_));
    XP->gemm(true,false,1.0,X_,FDSmSDF,0.0);
    XPX->gemm(false,false,1.0,XP,X_,0.0);

    //XPX->print();

    return XPX;
}

void myHF::print_stability_analysis(std::vector<std::pair<double, int> > &vec)
{
    std::sort(vec.begin(), vec.end());
    std::vector<std::pair<double, int> >::const_iterator iter = vec.begin();
    fprintf(outfile, "\t");
    char** irrep_labels = molecule_->irrep_labels();
    int count = 0;
    for(; iter != vec.end(); ++iter){
        ++count;
        fprintf(outfile, "%4s %-10.6f", irrep_labels[iter->second], iter->first);
        if(count == 4){
            fprintf(outfile, "\n\t");
            count = 0;
        }else{
            fprintf(outfile, "    ");
        }
    }
    if(count)
        fprintf(outfile, "\n\n");
    else
        fprintf(outfile, "\n");

    for(int h = 0; h < nirrep_; ++h)
        free(irrep_labels[h]);
    free(irrep_labels);
}
void myHF::stability_analysis()
{
    throw PSIEXCEPTION("Stability analysis hasn't been implemented yet for this wfn type.");
}
}}
/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 *  hf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

//#include "hf.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

void myHF::MOM_start()
{
    // Perhaps no MOM (or at least no MOM_start())? 
    if (iteration_ != options_.get_int("MOM_START") || iteration_ == 0 || MOM_started_)  return;
    
    // If we're here, we're doing MOM of some kind
    MOM_started_ = true;
    MOM_performed_ = true; // Gets printed next iteration   
    //
    // Build Ca_old_ matrices
    Ca_old_ = SharedMatrix(new Matrix("C Alpha Old (SO Basis)", nirrep_, nsopi_, nmopi_));
    if (!same_a_b_orbs()) {
        Cb_old_ = SharedMatrix(new Matrix("C Beta Old (SO Basis)", nirrep_, nsopi_, nmopi_));
    } else {
        Cb_old_ = Ca_old_;
    }
    
    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    // If no excitation requested, it's a stabilizer MOM, nothing fancy, don't print 
    if (!options_["MOM_OCC"].size()) return;      
   
    // If we're here, its an exciting MOM
    fprintf(outfile, "\n");
    print_orbitals(); 
    fprintf(outfile, "\n  ==> MOM Excited-State Iterations <==\n\n");
    
    // Reset iterations and DIIS (will automagically restart)
    iteration_ = 0;
    if (initialized_diis_manager_) {
        diis_manager_->delete_diis_file();
        diis_manager_.reset();
        initialized_diis_manager_ = false;
    }
 
    // Find out which orbitals are where
    std::vector<std::pair<double, std::pair<int, int> > > orbs_a;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_a_->pointer(h);
        for (int a = 0; a < nmo; a++)
            orbs_a.push_back(make_pair(eps[a], make_pair(h, a)));
    }
    std::vector<std::pair<double, std::pair<int, int> > > orbs_b;
    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi_[h];
        if (nmo == 0) continue;
        double* eps = epsilon_b_->pointer(h);
        for (int a = 0; a < nmo; a++)
            orbs_b.push_back(make_pair(eps[a], make_pair(h, a)));
    }
    sort(orbs_a.begin(),orbs_a.end());
    sort(orbs_b.begin(),orbs_b.end());

    if (options_["MOM_OCC"].size() != options_["MOM_VIR"].size())
        throw PSIEXCEPTION("SCF: MOM_OCC and MOM_VIR are not the same size");
    
    std::set<int> tot_check;

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        tot_check.insert(options_["MOM_OCC"][n].to_integer());
        tot_check.insert(options_["MOM_VIR"][n].to_integer());
    }

    if (tot_check.size() != 2*options_["MOM_OCC"].size())
        throw PSIEXCEPTION("SCF::MOM_start: Occupied/Virtual index array elements are not unique");

    CharacterTable ct = molecule_->point_group()->char_table();

    fprintf(outfile, "  Excitations:\n");

    for (int n = 0; n < options_["MOM_OCC"].size(); n++) {
        int i = options_["MOM_OCC"][n].to_integer(); 
        int a = options_["MOM_VIR"][n].to_integer();
       
        bool si = (i > 0);
        bool sa = (a > 0);
    
        i = abs(i) - 1;
        a = abs(a) - 1;

        if (options_.get_str("REFERENCE") == "RHF" || options_.get_str("REFERENCE") == "RKS") {

            if (!si || !sa) 
                throw PSIEXCEPTION("SCF::MOM_start: RHF/RKS requires + -> + in input, as only double excitations are performed");

            int hi = orbs_a[i].second.first;
            int ha = orbs_a[a].second.first;
 
            if (hi == ha) {
                // Same irrep
                int pi = orbs_a[i].second.second;
                int pa = orbs_a[a].second.second;
  
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pa];
                eps[pa] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
            
                delete[] Ct;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            } else {
                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
            
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[hi]--;
                nbetapi_[hi]--;
                doccpi_[hi]--; 

                // Vir -> Occ 
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_a_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
            
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[ha]++;
                nbetapi_[ha]++;
                doccpi_[ha]++; 

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "AB -> AB", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            }

        } else if (options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS") {

            if (si && sa) {

                // Alpha - alpha
                int hi = orbs_a[i].second.first;
                int ha = orbs_a[a].second.first;
 
                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_a[i].second.second;
                    int pa = orbs_a[a].second.second;
  
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_a_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
                
                    delete[] Ct;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_a[i].second.second;
                    int pi2 = nalphapi_[hi] - 1;  
 
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Ca_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_a_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nalphapi_[hi]--;

                    // Vir -> Occ 
                    int pa = orbs_a[a].second.second;
                    int pa2 = nalphapi_[ha];  
 
                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Ca_->pointer(ha);
                    Ct = new double[nso]; 
                    eps = epsilon_a_->pointer(ha);
     
                    // Swap eigvals 
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nalphapi_[ha]++;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                }
            } else if (!si && !sa) {
                // Beta->Beta
                int hi = orbs_b[i].second.first;
                int ha = orbs_b[a].second.first;
 
                if (hi == ha) {
                    // Same irrep
                    int pi = orbs_b[i].second.second;
                    int pa = orbs_b[a].second.second;
  
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_b_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pa];
                    eps[pa] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa], nmo);
                
                    delete[] Ct;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                } else {
                    // Different irrep
                    // Occ -> Vir
                    int pi = orbs_b[i].second.second;
                    int pi2 = nbetapi_[hi] - 1;  
 
                    int nso = nsopi_[hi];
                    int nmo = nmopi_[hi];

                    double** Ca = Cb_->pointer(hi);
                    double* Ct = new double[nso]; 
                    double* eps = epsilon_b_->pointer(hi);
                    double  epst;
     
                    // Swap eigvals 
                    epst = eps[pi];
                    eps[pi] = eps[pi2];
                    eps[pi2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nbetapi_[hi]--;

                    // Vir -> Occ 
                    int pa = orbs_b[a].second.second;
                    int pa2 = nbetapi_[ha];  
 
                    nso = nsopi_[ha];
                    nmo = nmopi_[ha];

                    Ca = Cb_->pointer(ha);
                    Ct = new double[nso]; 
                    eps = epsilon_b_->pointer(ha);
     
                    // Swap eigvals 
                    epst = eps[pa];
                    eps[pa] = eps[pa2];
                    eps[pa2] = epst; 

                    // Swap eigvecs
                    C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                    C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                    C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                    delete[] Ct;
                    
                    // Redo indexing
                    nbetapi_[ha]++;

                    fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
                }
            } else if (!si && sa) {
                // Beta->Alpha
                int hi = orbs_b[i].second.first;
                int ha = orbs_a[a].second.first;
 
                // Different irrep
                // Occ -> Vir
                int pi = orbs_b[i].second.second;
                int pi2 = nbetapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Cb_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_b_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nbetapi_[hi]--;
                nbeta_--;

                // Vir -> Occ 
                int pa = orbs_a[a].second.second;
                int pa2 = nalphapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Ca_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_a_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[ha]++;
                nalpha_++;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "B  -> A ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            } else if (si && !sa) {
                // Alpha->Beta
                int hi = orbs_a[i].second.first;
                int ha = orbs_b[a].second.first;
 
                // Different irrep
                // Occ -> Vir
                int pi = orbs_a[i].second.second;
                int pi2 = nalphapi_[hi] - 1;  
 
                int nso = nsopi_[hi];
                int nmo = nmopi_[hi];

                double** Ca = Ca_->pointer(hi);
                double* Ct = new double[nso]; 
                double* eps = epsilon_a_->pointer(hi);
                double  epst;
     
                // Swap eigvals 
                epst = eps[pi];
                eps[pi] = eps[pi2];
                eps[pi2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pi], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pi2], nmo, &Ca[0][pi], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pi2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nalphapi_[hi]--;
                nalpha_--;

                // Vir -> Occ 
                int pa = orbs_b[a].second.second;
                int pa2 = nbetapi_[ha];  
 
                nso = nsopi_[ha];
                nmo = nmopi_[ha];

                Ca = Cb_->pointer(ha);
                Ct = new double[nso]; 
                eps = epsilon_b_->pointer(ha);
     
                // Swap eigvals 
                epst = eps[pa];
                eps[pa] = eps[pa2];
                eps[pa2] = epst; 

                // Swap eigvecs
                C_DCOPY(nso, &Ca[0][pa], nmo, Ct, 1);
                C_DCOPY(nso, &Ca[0][pa2], nmo, &Ca[0][pa], nmo);
                C_DCOPY(nso, Ct, 1, &Ca[0][pa2], nmo);
                
                delete[] Ct;
                
                // Redo indexing
                nbetapi_[ha]++;
                nbeta_++;

                fprintf(outfile, "   %8s: %4d%-4s -> %4d%-4s \n", "A  -> B ", pi + 1, ct.gamma(hi).symbol(), pa + 1, ct.gamma(ha).symbol()); 
            }
            if (nalpha_ < nbeta_) throw PSIEXCEPTION("PSI::MOM_start: Nbeta ends up being less than Nalpha, this is not supported");
           
            // Fix doccpi/soccpi. 
            for (int h = 0; h < nirrep_; h++) {
                std::vector<std::pair<double,std::pair<int, bool> > > alphas;    
                std::vector<std::pair<double,std::pair<int, bool> > > betas;

                for (int i = 0; i < nalphapi_[h]; i++) {
                    alphas.push_back(make_pair(epsilon_a_->get(h,i), make_pair(i, true)));                                
                }
                for (int i = 0; i < nbetapi_[h]; i++) {
                    betas.push_back(make_pair(epsilon_b_->get(h,i), make_pair(i, true)));                                
                }
                for (int i = nalphapi_[h]; i < nmopi_[h]; i++) {
                    alphas.push_back(make_pair(epsilon_a_->get(h,i), make_pair(i, false)));                                
                }
                for (int i = nbetapi_[h]; i < nmopi_[h]; i++) {
                    betas.push_back(make_pair(epsilon_b_->get(h,i), make_pair(i, false)));                                
                }
                sort(alphas.begin(),alphas.end());       
                sort(betas.begin(),betas.end());       

                doccpi_[h] = 0;
                soccpi_[h] = 0;

                for (int i = 0; i < nmopi_[h]; i++) {
                    bool alpha_occ = alphas[i].second.second;
                    bool beta_occ = betas[i].second.second;
                    if (alpha_occ && beta_occ)
                        doccpi_[h]++;
                    else if (alpha_occ || beta_occ) // Careful, could be beta occ
                        soccpi_[h]++;
                }
            }

        } else if (options_.get_str("REFERENCE") == "ROHF") {
            throw PSIEXCEPTION("SCF::MOM_start: MOM excited states are not implemented for ROHF");
        }
    }
    Ca_old_->copy(Ca_);
    Cb_old_->copy(Cb_);

    fprintf(outfile, "\n                        Total Energy        Delta E      Density RMS\n\n");
}
void myHF::MOM()
{
    // Go MOM go!
    // Alpha
    for (int h = 0; h < nirrep_; h++) {
    
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nalpha = nalphapi_[h];
    
        if (nso == 0 || nmo == 0 || nalpha == 0) continue;
        
        double** Cold = Ca_old_->pointer(h);
        double** Cnew = Ca_->pointer(h);
        double*  eps  = epsilon_a_->pointer(h);
        double** S = S_->pointer(h);

        double* c = new double[nso];
        double* d = new double[nso];
        double* p = new double[nmo]; 
        
        memset(static_cast<void*>(c), '\0', sizeof(double)*nso);
        
        for (int i = 0; i < nalpha; i++)
            C_DAXPY(nso,1.0,&Cold[0][i],nmo,c,1);

        C_DGEMV('N',nso,nso,1.0,S[0],nso,c,1,0.0,d,1);
        C_DGEMV('T',nso,nmo,1.0,Cnew[0],nmo,d,1,0.0,p,1);

        //fprintf(outfile,"  P_a:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++)
            pvec[a] = make_pair(fabs(p[a]), a);
        sort(pvec.begin(),pvec.end(), greater<std::pair<double, int> >());

        //fprintf(outfile,"  P_a sorted:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nalpha);
        for (int a = 0; a < nalpha; a++)
            occvec[a] = make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(),occvec.end());

        //fprintf(outfile,"  P_a_occ sorted:\n");
        //for (int a = 0; a < nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nalpha);
        for (int a = 0; a < nmo - nalpha; a++)
            virvec[a] = make_pair(eps[pvec[a + nalpha].second], pvec[a + nalpha].second);
        sort(virvec.begin(),virvec.end());

        //fprintf(outfile,"  P_a_vir sorted:\n");
        //for (int a = 0; a < nmo - nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso,nmo);

        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]),static_cast<void*>(Cnew[0]),sizeof(double)*nso*nmo);
        memcpy(static_cast<void*>(p),      static_cast<void*>(eps),    sizeof(double)*nmo);

        for (int a = 0; a < nalpha; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso,&Ct[0][occvec[a].second],nmo,&Cnew[0][a],nmo);    
        } 

        for (int a = 0; a < nmo - nalpha; a++) {
            eps[a + nalpha] = virvec[a].first;
            C_DCOPY(nso,&Ct[0][virvec[a].second],nmo,&Cnew[0][a + nalpha],nmo);    
        } 

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }   
   
    // Beta
    if (same_a_b_orbs()) return;
  
    for (int h = 0; h < nirrep_; h++) {
    
        // Indexing
        int nso = nsopi_[h];
        int nmo = nmopi_[h];
        int nbeta = nbetapi_[h];
    
        if (nso == 0 || nmo == 0 || nbeta == 0) continue;
        
        double** Cold = Cb_old_->pointer(h);
        double** Cnew = Cb_->pointer(h);
        double*  eps  = epsilon_b_->pointer(h);
        double** S = S_->pointer(h);

        double* c = new double[nso];
        double* d = new double[nso];
        double* p = new double[nmo]; 
        
        memset(static_cast<void*>(c), '\0', sizeof(double)*nso);
        
        for (int i = 0; i < nbeta; i++)
            C_DAXPY(nso,1.0,&Cold[0][i],nmo,c,1);

        C_DGEMV('N',nso,nso,1.0,S[0],nso,c,1,0.0,d,1);
        C_DGEMV('T',nso,nmo,1.0,Cnew[0],nmo,d,1,0.0,p,1);

        //fprintf(outfile,"  P_a:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: %14.10f\n", a + 1, p[a]);

        // Find the largest contributions
        std::vector<std::pair<double, int> > pvec;
        pvec.resize(nmo);
        for (int a = 0; a < nmo; a++)
            pvec[a] = make_pair(fabs(p[a]), a);
        sort(pvec.begin(),pvec.end(), greater<std::pair<double, int> >());

        //fprintf(outfile,"  P_a sorted:\n");
        //for (int a = 0; a < nmo; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, pvec[a].second, pvec[a].first);

        // Now order the mos in each group
        std::vector<std::pair<double, int> > occvec;
        occvec.resize(nbeta);
        for (int a = 0; a < nbeta; a++)
            occvec[a] = make_pair(eps[pvec[a].second], pvec[a].second);
        sort(occvec.begin(),occvec.end());

        //fprintf(outfile,"  P_a_occ sorted:\n");
        //for (int a = 0; a < nbeta; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, occvec[a].second, occvec[a].first);

        std::vector<std::pair<double, int> > virvec;
        virvec.resize(nmo - nbeta);
        for (int a = 0; a < nmo - nbeta; a++)
            virvec[a] = make_pair(eps[pvec[a + nbeta].second], pvec[a + nbeta].second);
        sort(virvec.begin(),virvec.end());

        //fprintf(outfile,"  P_a_vir sorted:\n");
        //for (int a = 0; a < nmo - nalpha; a++)
        //    fprintf(outfile,"   a = %3d: Index = %3d, %14.10f\n", a + 1, virvec[a].second, virvec[a].first);

        double** Ct = block_matrix(nso,nmo);
        
        // Use Cold and p as a buffer
        memcpy(static_cast<void*>(Ct[0]),static_cast<void*>(Cnew[0]),sizeof(double)*nso*nmo);
        memcpy(static_cast<void*>(p),      static_cast<void*>(eps),    sizeof(double)*nmo);

        for (int a = 0; a < nbeta; a++) {
            eps[a] = occvec[a].first;
            C_DCOPY(nso,&Ct[0][occvec[a].second],nmo,&Cnew[0][a],nmo);    
        } 

        for (int a = 0; a < nmo - nbeta; a++) {
            eps[a + nbeta] = virvec[a].first;
            C_DCOPY(nso,&Ct[0][virvec[a].second],nmo,&Cnew[0][a + nbeta],nmo);    
        } 

        free_block(Ct);

        delete[] c;
        delete[] d;
        delete[] p;
    }   
    Cb_old_->copy(Cb_);
}

}}
/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 * frac.cc. So close to a great Battlestar Galactica joke. 
 * The way this works:
 * 1) form_C returns 1's normalized C matrices.
 * 2) find_occupation determines the occupations of said matrix via either Aufbau or MOM selection
 * 3) find_occupation then calls this, which renormalizes the C matrices with \sqrt(val) for each frac occ
 * 4) find_D is then computed with the renormalized C matrices, and everything is transparent until 1) on the next cycle
 *
 * Some executive decisions:
 *  -DIIS: Upon FRAC start, the old DIIS info is nuked, and DIIS_START is incremented by the current iteration count. 
 *         Thus, DIIS begins again on the next iteration. The exception is if you set FRAC_DIIS to false,
 *         in which case DIIS will cease for good upon FRAC start. 
 *  -MOM:  To use MOM with FRAC, set MOM_START to either FRAC_START or FRAC_START+1 (I think I prefer the latter).
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <libmints/mints.h>
#include <libqt/qt.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

//#include "hf.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

void myHF::frac()
{
    // Perhaps no frac?
    if (iteration_ < options_.get_int("FRAC_START") || options_.get_int("FRAC_START") == 0)  return;
    frac_performed_ = true;    

    // First frac iteration, blow away the diis and print the frac task
    if (iteration_ == options_.get_int("FRAC_START")) {

        // Throw unless UHF/UKS
        if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS"))
            throw PSIEXCEPTION("Fractional Occupation SCF is only implemented for UHF/UKS");

        // Throw if no frac tasks
        if (!options_["FRAC_OCC"].size())
            throw PSIEXCEPTION("Fractional Occupation SCF requested, but empty FRAC_OCC/FRAC_VAL vector");

        // Throw if inconsistent size
        if (options_["FRAC_OCC"].size() != options_["FRAC_VAL"].size()) 
            throw PSIEXCEPTION("Fractional Occupation SCF: FRAC_OCC/FRAC_VAL are of different dimensions");

        // Throw if the user is being an idiot with docc/socc
        if (input_docc_ || input_socc_) 
            throw PSIEXCEPTION("Fractional Occupation SCF: Turn off DOCC/SOCC");

        // Throw if the user is trying to start MOM before FRAC
        if (options_.get_int("MOM_START") <= options_.get_int("FRAC_START") && options_.get_int("MOM_START") != 0)
            throw PSIEXCEPTION("Fractional Occupation SCF: MOM must start after FRAC");

        // Throw if the use is just way too eager
        if (MOM_excited_) 
            throw PSIEXCEPTION("Fractional Occupation SCF: Don't try an excited-state MOM");
    
        // Close off a previous burn-in SCF
        fprintf(outfile, "\n");
        print_orbitals(); 
        
        // frac header
        fprintf(outfile, "\n  ==> Fractionally-Occupied SCF Iterations <==\n\n");
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double(); 
            
            // Throw if user requests frac occ above nalpha/nbeta 
            int max_i = (i > 0 ? nalpha_: nbeta_);
            if(abs(i) > max_i) {
                if (i > 0)
                    nalpha_++;
                else 
                    nbeta_++;
            }

            // Throw if the user is insane
            if (val < 0.0)
                throw PSIEXCEPTION("Fractional Occupation SCF: PSI4 is not configured for positrons. Please annihilate and start again");             

            fprintf(outfile, "    %-5s orbital %4d will contain %11.3E electron.\n", (i > 0 ? "Alpha" : "Beta"), abs(i), val);
        }
        fprintf(outfile, "\n");

        // Make sure diis restarts correctly/frac plays well with MOM
        if (initialized_diis_manager_) {
            diis_manager_->delete_diis_file();
            diis_manager_.reset();
            initialized_diis_manager_ = false;
            diis_start_ += iteration_ + 1;
        }

        // Turn yonder DIIS off if requested
        if (!options_.get_bool("FRAC_DIIS")) {
            diis_enabled_ = false;
        }

        // Load the old orbitals in if requested
        if (options_.get_bool("FRAC_LOAD")) {
            fprintf(outfile, "    Orbitals reloaded from file, your previous iterations are garbage.\n\n");
            load_orbitals();
        }

        // Keep the printing nice
        fprintf(outfile, "                        Total Energy        Delta E      Density RMS\n\n");
        fflush(outfile);

        // Prevent spurious convergence (technically this iteration comes from the N-electron system anyways)
        frac_performed_ = false;
    }
    
    // Every frac iteration: renormalize the Ca/Cb matrices

    // Sort the eigenvalues in the usual manner
    std::vector<boost::tuple<double,int,int> > pairs_a;
    std::vector<boost::tuple<double,int,int> > pairs_b;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs_a.push_back(boost::tuple<double,int,int>(epsilon_a_->get(h, i), h, i));
    }
    for (int h=0; h<epsilon_b_->nirrep(); ++h) {
        for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
            pairs_b.push_back(boost::tuple<double,int,int>(epsilon_b_->get(h, i), h, i));
    }
    sort(pairs_a.begin(),pairs_a.end());
    sort(pairs_b.begin(),pairs_b.end());
    
    // Renormalize the C matrix entries
    for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
        int i = options_["FRAC_OCC"][ind].to_integer();
        double val = options_["FRAC_VAL"][ind].to_double(); 
        bool is_alpha = (i > 0);
        i = abs(i) - 1; // Back to C ordering

        int i2  = ((is_alpha) ? get<2>(pairs_a[i]) : get<2>(pairs_b[i])); 
        int h   = ((is_alpha) ? get<1>(pairs_a[i]) : get<1>(pairs_b[i])); 

        int nso = Ca_->rowspi()[h];
        int nmo = Ca_->colspi()[h];

        double** Cp = ((is_alpha) ? Ca_->pointer(h) : Cb_->pointer(h));

        // And I say all that to say this
        C_DSCAL(nso, sqrt(val), &Cp[0][i], nmo); 
    } 
}
void myHF::frac_renormalize()
{
    if (!options_.get_bool("FRAC_RENORMALIZE") || !frac_enabled_) return;

    // Renormalize the fractional occupations back to 1, if possible before storage 
    fprintf(outfile, "    FRAC: Renormalizing orbitals to 1.0 for storage.\n\n");

    // Sort the eigenvalues in the usual manner
    std::vector<boost::tuple<double,int,int> > pairs_a;
    std::vector<boost::tuple<double,int,int> > pairs_b;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs_a.push_back(boost::tuple<double,int,int>(epsilon_a_->get(h, i), h, i));
    }
    for (int h=0; h<epsilon_b_->nirrep(); ++h) {
        for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
            pairs_b.push_back(boost::tuple<double,int,int>(epsilon_b_->get(h, i), h, i));
    }
    sort(pairs_a.begin(),pairs_a.end());
    sort(pairs_b.begin(),pairs_b.end());
    
    // Renormalize the C matrix entries
    for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
        int i = options_["FRAC_OCC"][ind].to_integer();
        double val = options_["FRAC_VAL"][ind].to_double(); 
        bool is_alpha = (i > 0);
        i = abs(i) - 1; // Back to C ordering

        int i2  = ((is_alpha) ? get<2>(pairs_a[i]) : get<2>(pairs_b[i])); 
        int h   = ((is_alpha) ? get<1>(pairs_a[i]) : get<1>(pairs_b[i])); 

        int nso = Ca_->rowspi()[h];
        int nmo = Ca_->colspi()[h];

        double** Cp = ((is_alpha) ? Ca_->pointer(h) : Cb_->pointer(h));

        // And I say all that to say this: TODO: This destroys FMP2 computations if val == 0
        if (val != 0.0) 
            C_DSCAL(nso, 1.0 / sqrt(val), &Cp[0][i], nmo); 
    } 
}

void myHF::compute_spin_contamination()
{
    if (!(options_.get_str("REFERENCE") == "UHF" || options_.get_str("REFERENCE") == "UKS" || options_.get_str("REFERENCE") == "CUHF"))
        return;

    double nalpha = (double) nalpha_;
    double nbeta  = (double) nbeta_;

    // Adjust for fractional occupation
    if (frac_performed_) {
        for (int ind = 0; ind < options_["FRAC_OCC"].size(); ind++) {
            int i = options_["FRAC_OCC"][ind].to_integer();
            double val = options_["FRAC_VAL"][ind].to_double(); 
            bool is_alpha = (i > 0);
            if (is_alpha) {
                nalpha -= (1.0 - val);
            } else {
                nbeta  -= (1.0 - val);
            }
        }    
    }

    SharedMatrix S = SharedMatrix(factory_->create_matrix("S (Overlap)"));
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basisset_,basisset_, basisset_,basisset_));
    boost::shared_ptr<OneBodySOInt> so_overlap(fact->so_overlap());
    so_overlap->compute(S);

    double dN = 0.0;

    for (int h =0; h < S->nirrep(); h++) {
        int nbf = S->colspi()[h];
        int nmo = Ca_->colspi()[h];
        int na = nalphapi_[h];
        int nb = nbetapi_[h];
        if (na == 0 || nb == 0 || nbf == 0 || nmo == 0)
            continue;

        SharedMatrix Ht (new Matrix("H Temp", nbf, nb));
        SharedMatrix Ft (new Matrix("F Temp", na, nb));

        double** Sp = S->pointer(h);
        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Htp = Ht->pointer(0);
        double** Ftp = Ft->pointer(0);

        C_DGEMM('N','N',nbf,nb,nbf,1.0,Sp[0],nbf,Cbp[0],nmo,0.0,Htp[0],nb);
        C_DGEMM('T','N',na,nb,nbf,1.0,Cap[0],nmo,Htp[0],nb,0.0,Ftp[0],nb);

        dN += C_DDOT(na*(long int)nb, Ftp[0], 1, Ftp[0], 1);
    }

    double nmin = (nbeta < nalpha ? nbeta : nalpha);
    double dS = nmin - dN;
    double nm = (nalpha - nbeta) / 2.0;
    double S2 = fabs(nm) * (fabs(nm) + 1.0);

    fprintf(outfile, "   @Spin Contamination Metric: %17.9E\n", dS);
    fprintf(outfile, "   @S^2 Expected:              %17.9E\n", S2);
    fprintf(outfile, "   @S^2 Observed:              %17.9E\n", S2 + dS);
    fprintf(outfile, "   @S   Expected:              %17.9E\n", nm);
    fprintf(outfile, "   @S   Observed:              %17.9E\n", nm);

    if (frac_performed_) {
        fprintf(outfile, "   @Nalpha:                    %17.9E\n", nalpha);
        fprintf(outfile, "   @Nbeta:                     %17.9E\n", nbeta);
    }
    fprintf(outfile, "\n");
}

}}
/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#ifndef LIBSCF_SAD_H
#define LIBSCF_SAD_H

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class Molecule;
class Matrix;

namespace scf {

class mySADGuess {

protected:

    int print_;
    int debug_;

    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
    SharedMatrix AO2SO_;

    int nalpha_;
    int nbeta_;

    Options& options_;

    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Ca_;
    SharedMatrix Cb_;

    void common_init();

    SharedMatrix form_D_AO();
    void getUHFAtomicDensity(boost::shared_ptr<BasisSet> atomic_basis, int n_electrons, int multiplicity, double** D);
    void atomicUHFHelperFormCandD(int nelec, int norbs,double** Shalf, double**F, double** C, double** D);

    void form_D();
    void form_C();

public:

    mySADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options);
    virtual ~mySADGuess();

    void compute_guess();
    
    SharedMatrix Da() const { return Da_; } 
    SharedMatrix Db() const { return Db_; } 
    SharedMatrix Ca() const { return Ca_; } 
    SharedMatrix Cb() const { return Cb_; } 

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }

}; 

}}

#endif
/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*
 *  sad.cc
 *
 * Routines for the high-meantenance SAD guess
 * and daul-basis projections
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

//#include "hf.h"
//#include "sad.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

mySADGuess::mySADGuess(boost::shared_ptr<BasisSet> basis, int nalpha, int nbeta, Options& options) :
    basis_(basis), nalpha_(nalpha), nbeta_(nbeta), options_(options)
{
    common_init();
}
mySADGuess::~mySADGuess()
{
}
void mySADGuess::common_init()
{
    molecule_ = basis_->molecule();

    boost::shared_ptr<IntegralFactory> ints(new IntegralFactory(basis_));
    boost::shared_ptr<PetiteList> petite(new PetiteList(basis_,ints));
    AO2SO_ =  petite->aotoso();

    print_ = options_.get_int("SAD_PRINT");
    debug_ = options_.get_int("DEBUG");
}
void mySADGuess::compute_guess()
{
    form_D();
    form_C();
}
void mySADGuess::form_D()
{
    // Build Neutral D in AO basis (block diagonal)
    SharedMatrix DAO = form_D_AO();

    // Transform Neutral D from AO to SO basis
    Da_ = SharedMatrix(new Matrix("Da SAD",AO2SO_->colspi(),AO2SO_->colspi()));

    double* temp = new double[AO2SO_->rowspi()[0] * (ULI) AO2SO_->max_ncol()];
    for (int h = 0; h < Da_->nirrep(); h++) {
        int nao = AO2SO_->rowspi()[h];
        int nso = AO2SO_->colspi()[h];
        if (!nao || !nso) continue;

        double** DAOp = DAO->pointer();
        double** DSOp = Da_->pointer(h);
        double** Up = AO2SO_->pointer(h);

        C_DGEMM('N','N',nao,nso,nao,1.0,DAOp[0],nao,Up[0],nso,0.0,temp,nso);
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,temp,nso,0.0,DSOp[0],nso);
    }
    delete[] temp;

    // Scale Da to true electron count
    double npair = 0.0;
    for (int A = 0; A < molecule_->natom(); A++) {
        npair += 0.5 * molecule_->Z(A);
    }
    Da_->scale(((double) nalpha_) / npair);

    // Build/Scale Db if needed
    if (nalpha_ == nbeta_) {
        Db_ = Da_;
    } else {
        Db_ = SharedMatrix(Da_->clone());
        Db_->set_name("Db SAD");
        Db_->scale(((double) nbeta_) / ((double) nalpha_));
    }

    if (debug_) {
        Da_->print();
        Db_->print();
    }
}
void mySADGuess::form_C()
{
    Ca_ = Da_->partial_cholesky_factorize(options_.get_double("SAD_CHOL_TOLERANCE"));
    Ca_->set_name("Ca SAD");
    if (nalpha_ == nbeta_) {

        Cb_ = Ca_;
    } else {
        Cb_ = SharedMatrix(Ca_->clone());
        Cb_->set_name("Cb SAD");
        Cb_->scale(sqrt(((double)nbeta_)/((double)nalpha_)));
    }

    if (debug_) {
        Ca_->print();
        Cb_->print();
    }
}
SharedMatrix mySADGuess::form_D_AO()
{
    std::vector<boost::shared_ptr<BasisSet> > atomic_bases;

    if (print_ > 6) {
        fprintf(outfile,"\n  Constructing atomic basis sets\n  Molecule:\n");
        molecule_->print();
    }

    //Build the atomic basis sets for libmints use in UHF
    for (int A = 0; A<molecule_->natom(); A++) {
        atomic_bases.push_back(basis_->atomic_basis_set(A));
        if (print_>6) {
            fprintf(outfile,"  SAD: Atomic Basis Set %d\n", A);
            atomic_bases[A]->molecule()->print();
            fprintf(outfile,"\n");
            atomic_bases[A]->print(outfile);
            fprintf(outfile,"\n");
        }
    }

    //Spin occupations per atom, to be determined by Hund's Rules
    //or user input
    int* nalpha = init_int_array(molecule_->natom());
    int* nbeta = init_int_array(molecule_->natom());
    int* nelec = init_int_array(molecule_->natom());
    int* nhigh = init_int_array(molecule_->natom());
    int tot_elec = 0;

    //Ground state high spin occupency array, atoms 0 to 36 (see Giffith's Quantum Mechanics, pp. 217)
    const int reference_S[] = {0,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,2,1,0,1,0,1,2,3,6,5,4,3,2,1,0,1,2,3,2,1,0};
    const int MAX_Z = 36;

    if (print_>1)
        fprintf(outfile,"  Determining Atomic Occupations\n");
    for (int A = 0; A<molecule_->natom(); A++) {
        int Z = molecule_->Z(A);
        if (Z>MAX_Z) {
            throw std::domain_error(" Only Atoms up to 36 (Kr) are currently supported with SAD Guess");
        }
        nhigh[A] = reference_S[Z];
        nelec[A] = Z;
        tot_elec+= nelec[A];
        nbeta[A] = (nelec[A]-nhigh[A])/2;
        nalpha[A] = nelec[A]-nbeta[A];
        if (print_>1)
            fprintf(outfile,"  Atom %d, Z = %d, nelec = %d, nhigh = %d, nalpha = %d, nbeta = %d\n",A,Z,nelec[A],nhigh[A],nalpha[A],nbeta[A]);
    }
    fflush(outfile);

    // Determine redundant atoms
    int* unique_indices = init_int_array(molecule_->natom()); // All atoms to representative unique atom
    int* atomic_indices = init_int_array(molecule_->natom()); // unique atom to first representative atom
    int* offset_indices = init_int_array(molecule_->natom()); // unique atom index to rank
    int nunique = 0;
    for (int l = 0; l < molecule_->natom(); l++) {
        unique_indices[l] = l;
        atomic_indices[l] = l;
    }
    for (int l = 0; l < molecule_->natom() - 1; l++) {
        for (int m = l + 1; m < molecule_->natom(); m++) {
            if (unique_indices[m] != m)
                continue; //Already assigned
            if (molecule_->Z(l) != molecule_->Z(m))
                continue;
            if (nalpha[l] != nalpha[m])
                continue;
            if (nbeta[l] != nbeta[m])
                continue;
            if (nhigh[l] != nhigh[m])
                continue;
            if (nelec[l] != nelec[m])
                continue;
            if (atomic_bases[l]->nbf() != atomic_bases[m]->nbf())
                continue;
            if (atomic_bases[l]->nshell() != atomic_bases[m]->nshell())
                continue;
            if (atomic_bases[l]->nprimitive() != atomic_bases[m]->nprimitive())
                continue;
            if (atomic_bases[l]->max_am() != atomic_bases[m]->max_am())
                continue;
            if (atomic_bases[l]->max_nprimitive() != atomic_bases[m]->max_nprimitive())
                continue;
            if (atomic_bases[l]->has_puream() !=  atomic_bases[m]->has_puream())
                continue;

            // Semi-Rigorous match obtained
            unique_indices[m] = l;
        }
    }
    for (int l = 0; l < molecule_->natom(); l++) {
        if (unique_indices[l] == l) {
            atomic_indices[nunique] = l;
            offset_indices[l] = nunique;
            nunique++;
        }
    }

    //Atomic D matrices within the atom specific AO basis
    double*** atomic_D = (double***)malloc(nunique*sizeof(double**));
    for (int A = 0; A<nunique; A++) {
        atomic_D[A] = block_matrix(atomic_bases[atomic_indices[A]]->nbf(),atomic_bases[atomic_indices[A]]->nbf());
    }

    if (print_ > 1)
        fprintf(outfile,"\n  Performing Atomic UHF Computations:\n");
    for (int A = 0; A<nunique; A++) {
        int index = atomic_indices[A];
        if (print_ > 1)
            fprintf(outfile,"\n  UHF Computation for Unique Atom %d which is Atom %d:",A, index);
        getUHFAtomicDensity(atomic_bases[index],nelec[index],nhigh[index],atomic_D[A]);
    }
    if (print_)
        fprintf(outfile,"\n");

    fflush(outfile);

    //Add atomic_D into D (scale by 1/2, we like effective pairs)
    SharedMatrix DAO = SharedMatrix(new Matrix("D_SAD (AO)", basis_->nbf(), basis_->nbf()));
    for (int A = 0, offset = 0; A < molecule_->natom(); A++) {
        int norbs = atomic_bases[A]->nbf();
        int back_index = unique_indices[A];
        for (int m = 0; m<norbs; m++)
            for (int n = 0; n<norbs; n++)
                DAO->set(0,m+offset,n+offset,0.5*atomic_D[offset_indices[back_index]][m][n]);
        offset+=norbs;
    }

    for (int A = 0; A<nunique; A++)
        free_block(atomic_D[A]);
    free(atomic_D);

    free(nelec);
    free(nhigh);
    free(nalpha);
    free(nbeta);

    free(atomic_indices);
    free(unique_indices);
    free(offset_indices);

    if (debug_) {
        DAO->print();
    }

    return DAO;
}
void mySADGuess::getUHFAtomicDensity(boost::shared_ptr<BasisSet> bas, int nelec, int nhigh, double** D)
{
    boost::shared_ptr<Molecule> mol = bas->molecule();

    int nbeta = (nelec-nhigh)/2;
    int nalpha = nelec-nbeta;
    int natom = mol->natom();
    int norbs = bas->nbf();

    if (print_>1) {
        fprintf(outfile,"\n");
        bas->print(outfile);
        fprintf(outfile,"  Occupation: nalpha = %d, nbeta = %d, norbs = %d\n",nalpha,nbeta,norbs);
        fprintf(outfile,"\n  Atom:\n");
        mol->print();
    }

    if (natom != 1) {
        throw std::domain_error("SAD Atomic UHF has been given a molecule, not an atom");
    }

    double** Dold = block_matrix(norbs,norbs);
    double **Shalf = block_matrix(norbs, norbs);
    double** Ca = block_matrix(norbs,norbs);
    double** Cb = block_matrix(norbs,norbs);
    double** Da = block_matrix(norbs,norbs);
    double** Db = block_matrix(norbs,norbs);
    double** Fa = block_matrix(norbs,norbs);
    double** Fb = block_matrix(norbs,norbs);
    double** Fa_old = block_matrix(norbs,norbs);
    double** Fb_old = block_matrix(norbs,norbs);
    double** Ga = block_matrix(norbs,norbs);
    double** Gb = block_matrix(norbs,norbs);

    IntegralFactory integral(bas, bas, bas, bas);
    MatrixFactory mat;
    mat.init_with(1,&norbs,&norbs);
    OneBodyAOInt *S_ints = integral.ao_overlap();
    OneBodyAOInt *T_ints = integral.ao_kinetic();
    OneBodyAOInt *V_ints = integral.ao_potential();
    TwoBodyAOInt *TEI = integral.eri();

    //Compute Shalf;
    //Fill S
    SharedMatrix S_UHF(mat.create_matrix("S_UHF"));
    S_ints->compute(S_UHF);
    double** S = S_UHF->to_block_matrix();

    if (print_>6) {
        fprintf(outfile,"  S:\n");
        print_mat(S,norbs,norbs,outfile);
    }
    // S^{-1/2}

    // First, diagonalize S
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(norbs);
    int lwork = norbs * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',norbs,S[0],norbs,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    double **S_copy = block_matrix(norbs, norbs);
    C_DCOPY(norbs*norbs,S[0],1,S_copy[0],1);

    // Now form S^{-1/2} = U(T)*s^{-1/2}*U,
    // where s^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of S
    for (int i=0; i<norbs; i++) {
        if (eigval[i] < 1.0E-10)
            eigval[i] = 0.0;
        else
            eigval[i] = 1.0 / sqrt(eigval[i]);

        // scale one set of eigenvectors by the diagonal elements s^{-1/2}
        C_DSCAL(norbs, eigval[i], S[i], 1);
    }
    free(eigval);

    // Smhalf = S_copy(T) * S
    C_DGEMM('t','n',norbs,norbs,norbs,1.0,
            S_copy[0],norbs,S[0],norbs,0.0,Shalf[0],norbs);

    free_block(S);
    free_block(S_copy);

    if (print_>6) {
        fprintf(outfile,"  S^-1/2:\n");
        print_mat(Shalf,norbs,norbs,outfile);
    }

    //Compute H
    SharedMatrix T_UHF(mat.create_matrix("T_UHF"));
    T_ints->compute(T_UHF);
    SharedMatrix V_UHF(mat.create_matrix("V_UHF"));
    V_ints->compute(V_UHF);
    SharedMatrix H_UHF(mat.create_matrix("H_UHF"));
    H_UHF->zero();
    H_UHF->add(T_UHF);
    H_UHF->add(V_UHF);
    double** H = H_UHF->to_block_matrix();

    delete S_ints;
    delete T_ints;
    delete V_ints;

    if (print_>6) {
        fprintf(outfile,"  H:\n");
        print_mat(H,norbs,norbs,outfile);
    }

    //Compute initial Ca and Da from core guess
    atomicUHFHelperFormCandD(nalpha,norbs,Shalf,H,Ca,Da);
    //Compute initial Cb and Db from core guess
    atomicUHFHelperFormCandD(nbeta,norbs,Shalf,H,Cb,Db);
    //Compute intial D
    C_DCOPY(norbs*norbs,Da[0],1,D[0],1);
    C_DAXPY(norbs*norbs,1.0,Db[0],1,D[0],1);
    if (print_>6) {
        fprintf(outfile,"  Ca:\n");
        print_mat(Ca,norbs,norbs,outfile);

        fprintf(outfile,"  Cb:\n");
        print_mat(Cb,norbs,norbs,outfile);

        fprintf(outfile,"  Da:\n");
        print_mat(Da,norbs,norbs,outfile);

        fprintf(outfile,"  Db:\n");
        print_mat(Db,norbs,norbs,outfile);

        fprintf(outfile,"  D:\n");
        print_mat(D,norbs,norbs,outfile);
    }

    //Compute inital E for reference
    double E = C_DDOT(norbs*norbs,D[0],1,H[0],1);
    E += C_DDOT(norbs*norbs,Da[0],1,Fa[0],1);
    E += C_DDOT(norbs*norbs,Db[0],1,Fb[0],1);
    E *= 0.5;

    const double* buffer = TEI->buffer();

    double E_tol = options_.get_double("SAD_E_CONVERGENCE");
    double D_tol = options_.get_double("SAD_D_CONVERGENCE");
    int maxiter = options_.get_int("SAD_MAXITER");
    int f_mixing_iteration = options_.get_int("SAD_F_MIX_START");

    double E_old;
    int iteration = 0;

    bool converged = false;
    if (print_>1) {
        fprintf(outfile, "\n  Initial Atomic UHF Energy:    %14.10f\n\n",E);
        fprintf(outfile, "                                         Total Energy            Delta E              Density RMS\n\n");
        fflush(outfile);
    }
    do {

        iteration++;

        //Copy the old values over for error analysis
        E_old = E;
        //I'm only going to use the total for now, could be expanded later
        C_DCOPY(norbs*norbs,D[0],1,Dold[0],1);
        //And old Fock matrices for level shift
        C_DCOPY(norbs*norbs,Fa[0],1,Fa_old[0],1);
        C_DCOPY(norbs*norbs,Fb[0],1,Fb_old[0],1);

        //Form Ga and Gb via integral direct
        memset((void*) Ga[0], '\0',norbs*norbs*sizeof(double));
        memset((void*) Gb[0], '\0',norbs*norbs*sizeof(double));

        //At the moment this is 8-fold slower than it could be, we'll see if it is signficant
        for (int MU = 0; MU < bas->nshell(); MU++) {
        int numMU = bas->shell(MU).nfunction();
        for (int NU = 0; NU < bas->nshell(); NU++) {
        int numNU = bas->shell(NU).nfunction();
        for (int LA = 0; LA < bas->nshell(); LA++) {
        int numLA = bas->shell(LA).nfunction();
        for (int SI = 0; SI < bas->nshell(); SI++) {
        int numSI = bas->shell(SI).nfunction();
        TEI->compute_shell(MU,NU,LA,SI);
        for (int m = 0, index = 0; m < numMU; m++) {
        int omu = bas->shell(MU).function_index() + m;
        for (int n = 0; n < numNU; n++) {
        int onu = bas->shell(NU).function_index() + n;
        for (int l = 0; l < numLA; l++) {
        int ola = bas->shell(LA).function_index() + l;
        for (int s = 0; s < numSI; s++, index++) {
        int osi = bas->shell(SI).function_index() + s;
             //fprintf(outfile,"  Integral (%d, %d| %d, %d) = %14.10f\n",omu,onu,ola,osi,buffer[index]);
             Ga[omu][onu] += D[ola][osi]*buffer[index];
             //Ga[ola][osi] += D[omu][onu]*buffer[index];
             Ga[omu][osi] -= Da[onu][ola]*buffer[index];
             Gb[omu][onu] += D[ola][osi]*buffer[index];
             //Gb[ola][osi] += D[omu][onu]*buffer[index];
             Gb[omu][osi] -= Db[onu][ola]*buffer[index];
        } } } } } } } }

        //Form Fa and Fb
        C_DCOPY(norbs*norbs,H[0],1,Fa[0],1);
        C_DAXPY(norbs*norbs,1.0,Ga[0],1,Fa[0],1);
        C_DCOPY(norbs*norbs,H[0],1,Fb[0],1);
        C_DAXPY(norbs*norbs,1.0,Gb[0],1,Fb[0],1);

        //Compute E
        E = C_DDOT(norbs*norbs,D[0],1,H[0],1);
        E += C_DDOT(norbs*norbs,Da[0],1,Fa[0],1);
        E += C_DDOT(norbs*norbs,Db[0],1,Fb[0],1);
        E *= 0.5;

        //Perform any required convergence stabilization
        //F-mixing (should help with oscillation)
        //20% old, 80% new Fock matrix for now
        if (iteration >= f_mixing_iteration) {
            C_DSCAL(norbs*norbs,0.8,Fa[0],1);
            C_DSCAL(norbs*norbs,0.8,Fb[0],1);
            C_DAXPY(norbs*norbs,0.2,Fa_old[0],1,Fa[0],1);
            C_DAXPY(norbs*norbs,0.2,Fb_old[0],1,Fb[0],1);
        }

        //Diagonalize Fa and Fb to from Ca and Cb and Da and Db
        atomicUHFHelperFormCandD(nalpha,norbs,Shalf,Fa,Ca,Da);
        atomicUHFHelperFormCandD(nbeta,norbs,Shalf,Fb,Cb,Db);

        //Form D
        C_DCOPY(norbs*norbs,Da[0],1,D[0],1);
        C_DAXPY(norbs*norbs,1.0,Db[0],1,D[0],1);

        //Form delta D and Drms
        C_DAXPY(norbs*norbs,-1.0,D[0],1,Dold[0],1);
        double Drms = sqrt(1.0/(1.0*norbs*norbs)*C_DDOT(norbs*norbs,Dold[0],1,Dold[0],1));

        double deltaE = fabs(E-E_old);

        if (print_>6) {
            fprintf(outfile,"  Fa:\n");
            print_mat(Fa,norbs,norbs,outfile);

            fprintf(outfile,"  Fb:\n");
            print_mat(Fb,norbs,norbs,outfile);

            fprintf(outfile,"  Ga:\n");
            print_mat(Ga,norbs,norbs,outfile);

            fprintf(outfile,"  Gb:\n");
            print_mat(Gb,norbs,norbs,outfile);

            fprintf(outfile,"  Ca:\n");
            print_mat(Ca,norbs,norbs,outfile);

            fprintf(outfile,"  Cb:\n");
            print_mat(Cb,norbs,norbs,outfile);

            fprintf(outfile,"  Da:\n");
            print_mat(Da,norbs,norbs,outfile);

            fprintf(outfile,"  Db:\n");
            print_mat(Db,norbs,norbs,outfile);

            fprintf(outfile,"  D:\n");
            print_mat(D,norbs,norbs,outfile);
        }
        if (print_>1)
            fprintf(outfile, "  @Atomic UHF iteration %3d energy: %20.14f    %20.14f %20.14f\n", iteration, E, E-E_old, Drms);
        if (iteration > 1 && deltaE < E_tol && Drms < D_tol)
            converged = true;

        if (iteration > maxiter) {
            fprintf(outfile, "\n WARNING: Atomic UHF is not converging! Try casting from a smaller basis or call Rob at CCMST.\n");
            break;
        }

        //Check convergence
    } while (!converged);
    if (converged && print_ > 1)
        fprintf(outfile, "  @Atomic UHF Final Energy for atom %s: %20.14f\n", mol->symbol(0).c_str(),E);

    delete TEI;
    free_block(Dold);
    free_block(Ca);
    free_block(Cb);
    free_block(Da);
    free_block(Db);
    free_block(Fa);
    free_block(Fb);
    free_block(Fa_old);
    free_block(Fb_old);
    free_block(Ga);
    free_block(Gb);
    free_block(H);
    free_block(Shalf);
}
void mySADGuess::atomicUHFHelperFormCandD(int nelec, int norbs,double** Shalf, double**F, double** C, double** D)
{
    //Forms C in the AO basis for SAD Guesses
    double **Temp = block_matrix(norbs,norbs);
    double **Fp = block_matrix(norbs,norbs);
    double **Cp = block_matrix(norbs,norbs);

    //Form F' = X'FX = XFX for symmetric orthogonalization
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Shalf[0],norbs,F[0],norbs,0.0,Temp[0],norbs);
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Temp[0],norbs,Shalf[0],norbs,0.0,Fp[0],norbs);

    //Form C' = eig(F')
    double *eigvals = init_array(norbs);
    sq_rsp(norbs, norbs, Fp,  eigvals, 1, Cp, 1.0e-14);
    free(eigvals);

    //Form C = XC'
    C_DGEMM('N','N',norbs,norbs,norbs,1.0,Shalf[0],norbs,Cp[0],norbs,0.0,C[0],norbs);

    //Form D = Cocc*Cocc'
    C_DGEMM('N','T',norbs,norbs,nelec,1.0,C[0],norbs,C[0],norbs,0.0,D[0],norbs);

    free_block(Temp);
    free_block(Cp);
    free_block(Fp);
}

void myHF::compute_SAD_guess()
{
    boost::shared_ptr<mySADGuess> guess(new mySADGuess(basisset_,nalpha_,nbeta_,options_));
    guess->compute_guess();

    Da_->copy(guess->Da());
    Db_->copy(guess->Db());

    for (int h = 0; h < Da_->nirrep(); h++) {

        int nso = guess->Ca()->rowspi()[h];
        int nmo = guess->Ca()->colspi()[h];
        if (nmo > X_->colspi()[h])
            nmo = X_->colspi()[h];

        sad_nocc_[h] = nmo;

        if (!nso || !nmo) continue;

        double** Cap = Ca_->pointer(h);
        double** Cbp = Cb_->pointer(h);
        double** Ca2p = guess->Ca()->pointer(h);
        double** Cb2p = guess->Cb()->pointer(h);

        for (int i = 0; i < nso; i++) {
            ::memcpy((void*) Cap[i], (void*) Ca2p[i], nmo*sizeof(double));
            ::memcpy((void*) Cbp[i], (void*) Cb2p[i], nmo*sizeof(double));
        }
    }

    int temp_nocc;
    for (int h = 0 ; h < Da_->nirrep(); h++) {
        temp_nocc = sad_nocc_[h];
        sad_nocc_[h] = doccpi_[h];
        nalphapi_[h] = temp_nocc;
        nbetapi_[h]  = temp_nocc;
        doccpi_[h]   = temp_nocc;
        soccpi_[h]   = 0;
    }

    fflush(outfile);

    E_ = 0.0; // This is the -1th iteration
}
SharedMatrix myHF::BasisProjection(SharedMatrix C_A, int* noccpi, boost::shared_ptr<BasisSet> old_basis, boost::shared_ptr<BasisSet> new_basis)
{

    //Based on Werner's method from Mol. Phys. 102, 21-22, 2311
    boost::shared_ptr<IntegralFactory> newfactory(new IntegralFactory(new_basis,new_basis,new_basis,new_basis));
    boost::shared_ptr<IntegralFactory> hybfactory(new IntegralFactory(old_basis,new_basis,old_basis,new_basis));
    boost::shared_ptr<OneBodySOInt> intBB(newfactory->so_overlap());
    boost::shared_ptr<OneBodySOInt> intAB(hybfactory->so_overlap());

    boost::shared_ptr<PetiteList> pet(new PetiteList(new_basis, newfactory));
    SharedMatrix AO2USO(pet->aotoso());

    SharedMatrix SAB(new Matrix("S_AB", C_A->nirrep(), C_A->rowspi(), AO2USO->colspi()));
    SharedMatrix SBB(new Matrix("S_BB", C_A->nirrep(), AO2USO->colspi(), AO2USO->colspi()));

    intAB->compute(SAB);
    intBB->compute(SBB);

    //SAB->print();
    //SBB->print();

    newfactory.reset();
    hybfactory.reset();
    intAB.reset();
    intBB.reset();
    pet.reset();

    // Constrained to the same symmetry at the moment, we can relax this soon
    SharedMatrix C_B(new Matrix("C_B", C_A->nirrep(),AO2USO->colspi(), noccpi));

    // Block over irreps (soon united irreps)
    for (int h = 0; h < C_A->nirrep(); h++) {

        int nocc = noccpi[h];
        int na = C_A->rowspi()[h];
        int nb = AO2USO->colspi()[h];

        if (nocc == 0 || na == 0 || nb == 0) continue;

        double** Ca = C_A->pointer(h);
        double** Cb = C_B->pointer(h);
        double** Sab = SAB->pointer(h);
        double** Sbb = SBB->pointer(h);

        int CholError = C_DPOTRF('L',nb,Sbb[0],nb);
        if (CholError !=0 )
            throw std::domain_error("S_BB Matrix Cholesky failed!");

        //Inversion (in place)
        int IError = C_DPOTRI('L',nb,Sbb[0],nb);
        if (IError !=0 )
            throw std::domain_error("S_BB Inversion Failed!");

        //LAPACK is smart and all, only uses half of the thing
        for (int m = 0; m<nb; m++)
            for (int n = 0; n<m; n++)
                Sbb[m][n] = Sbb[n][m];

        //Form T
        double** Temp1 = block_matrix(nb,nocc);
        C_DGEMM('T','N',nb,nocc,na,1.0,Sab[0],nb,Ca[0],nocc,0.0,Temp1[0],nocc);

        //fprintf(outfile," Temp1:\n");
        //print_mat(Temp1,nb,nocc,outfile);

        double** Temp2 = block_matrix(nb,nocc);
        C_DGEMM('N','N',nb,nocc,nb,1.0,Sbb[0],nb,Temp1[0],nocc,0.0,Temp2[0],nocc);

        //fprintf(outfile," Temp2:\n");
        //print_mat(Temp2,nb,nocc,outfile);

        double** Temp3 = block_matrix(na,nocc);
        C_DGEMM('N','N',na,nocc,nb,1.0,Sab[0],nb,Temp2[0],nocc,0.0,Temp3[0],nocc);

        //fprintf(outfile," Temp3:\n");
        //print_mat(Temp3,na,nocc,outfile);

        double** T = block_matrix(nocc,nocc);
        C_DGEMM('T','N',nocc,nocc,na,1.0,Ca[0],nocc,Temp3[0],nocc,0.0,T[0],nocc);

        //fprintf(outfile," T:\n");
        //print_mat(T,nocc,nocc,outfile);

        //Find T^-1/2
        // First, diagonalize T
        // the C_DSYEV call replaces the original matrix T with its eigenvectors
        double* eigval = init_array(nocc);
        int lwork = nocc * 3;
        double* work = init_array(lwork);
        int stat = C_DSYEV('v','u',nocc,T[0],nocc,eigval, work,lwork);
        if (stat != 0) {
            fprintf(outfile, "C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }
        free(work);

        // Now T contains the eigenvectors of the original T
        // Copy T to T_copy
        double **T_mhalf = block_matrix(nocc, nocc);
        double **T_copy = block_matrix(nocc, nocc);
        C_DCOPY(nocc*nocc,T[0],1,T_copy[0],1);

        // Now form T^{-1/2} = U(T)*t^{-1/2}*U,
        // where t^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of T
        for (int i=0; i<nocc; i++) {
            if (eigval[i] < 1.0E-10)
                eigval[i] = 0.0;
            else
                eigval[i] = 1.0 / sqrt(eigval[i]);

            // scale one set of eigenvectors by the diagonal elements t^{-1/2}
            C_DSCAL(nocc, eigval[i], T[i], 1);
        }
        free(eigval);

        // T_mhalf = T_copy(T) * T
        C_DGEMM('t','n',nocc,nocc,nocc,1.0,
                T_copy[0],nocc,T[0],nocc,0.0,T_mhalf[0],nocc);

        //Form CB
        C_DGEMM('N','N',nb,nocc,nocc,1.0,Temp2[0],nocc,T_mhalf[0],nocc,0.0,Cb[0],nocc);

        free_block(Temp1);
        free_block(Temp2);
        free_block(Temp3);
        free_block(T);
        free_block(T_copy);
        free_block(T_mhalf);

    }
    return C_B;
}

}}
