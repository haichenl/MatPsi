#include <libmints/mints.h>
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
    int natom_;
    int nbasis_;
    int nelec_;
    boost::shared_ptr<Molecule> molecule_;
    boost::shared_ptr<BasisSet> basis_;
	boost::shared_ptr<IntegralFactory> intfac_;
    boost::shared_ptr<TwoBodyAOInt> eri_;
	boost::shared_ptr<MatrixFactory> matfac_;
    
public:
    // constructor; takes in 2 strings and parse them 
    MatPsi(std::string mol_string, std::string basis_name);
	
    // copy constructor 
	MatPsi(boost::shared_ptr<MatPsi> inputMatPsi);
	
	~MatPsi() {
	}
    
    // number of atoms 
    int natom() {
        return natom_;
    }
    
    // number of basis functions 
    int nbasis() {
        return nbasis_;
    }
    
    // number of electrons 
    int nelec() {
        return nelec_;
    }
    
    // atom coordinates 
    SharedMatrix coord() {
        SharedMatrix coordMat = (molecule_->geometry()).clone();
        return coordMat;
    }
	
    // compute the overlap matrix S 
	SharedMatrix overlap();
    
    // compute the kinetic energy matrix KE 
    SharedMatrix kinetic();
    
    // compute the total potential energy matrix EN 
    SharedMatrix potential();
    
    // compute the atom-separated potential energy matrix ENI 
    boost::shared_array<SharedMatrix> potential_sep();
    
    // compute the 4-indexed two electron integral H2(i, j, k, l) 
    double tei_ijkl(int i, int j, int k, int l);
    
};