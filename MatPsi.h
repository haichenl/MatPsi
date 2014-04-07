#include <libmints/mints.h>
#include <libfock/jk.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include "mex.h"
#include "class_handle.hpp"
#include <boost/shared_array.hpp>

using namespace std;
using namespace psi;
using namespace boost;

namespace psi {
    FILE* outfile = stdout;
    char* psi_file_prefix = NULL;
    std::string outfile_name = "";
    extern int read_options(const std::string &name, Options & options, bool suppress_printing = false);
}

class MatPsi {
protected:
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
	boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
	boost::shared_ptr<MatrixFactory> matfac_;
    boost::shared_ptr<DirectJK> directjk_;
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi(std::string mol_string, std::string basis_name);
	
    // copy constructor 
	MatPsi(boost::shared_ptr<MatPsi> inputMatPsi);
	
	~MatPsi() {}
    
    // Molecule properties 
    // number of atoms 
    int natom() { return molecule_->natom(); }
    
    // atom coordinates 
    SharedMatrix coord() { return (molecule_->geometry()).clone(); }
    
    // nuclear repulsion energy 
    double Enuc() { return molecule_->nuclear_repulsion_energy(); }
    
    // Z list 
    SharedVector Zlist();
    
    // number of electrons 
    int nelec();
    
    // Basis set properties 
    // number of basis functions 
    int nbasis() { return basis_->nbf(); }
    
    // map basis number to the number of atom it is centred on 
    SharedVector func2center();
    
    // map basis number to its angular momentum 
    SharedVector func2am();
    
	// One-electron integrals 
    // compute the overlap matrix S 
	SharedMatrix overlap();
    
    // compute the kinetic energy matrix KE 
    SharedMatrix kinetic();
    
    // compute the total potential energy matrix EN 
    SharedMatrix potential();
    
    // compute the atom-separated potential energy matrix ENI 
    boost::shared_array<SharedMatrix> potential_sep();
    
    // compute from a given point charge the environment potential energy matrix ENVI
    SharedMatrix potential_zxyz(const double* Zxyz_array);
    
    // Two-electron integrals 
    // compute the 4-indexed two-electron integral H2(i, j, k, l) 
    double tei_ijkl(int i, int j, int k, int l);
    
    // compute all unique two-electron integrals and put them in a vector; be careful as it costs a huge amount of memory 
    SharedVector tei_alluniq();
    
    // compute all unique two-electron integrals and pre-arrange them for the forming of J and K 
    boost::shared_array<SharedVector> tei_alluniqJK();
    
    // SCF related 
    // for Hartree Fock, compute 2-electron G matrix from occupied molecular orbital coefficient matrix, direct algorithm, consider no geometrical symmetry 
    SharedMatrix HFnosymmMO2G(SharedMatrix coeff, long int memory = 1000000000L, double cutoff = 1.0E-12);
    
};