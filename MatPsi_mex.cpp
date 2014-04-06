#include "mex.h"
#include "class_handle.hpp"
#include "MatPsi.h"
#include "MatPsi.cc"

using namespace std;
using namespace psi;
using namespace boost;

SharedMatrix InputMatrix(mxArray*& Mat_m) {
	int nrow = mxGetM(Mat_m);
	int ncol = mxGetN(Mat_m);
    SharedMatrix Mat_c(new Matrix(nrow, ncol));
	double* Mat_m_pt = mxGetPr(Mat_m);
	double* Mat_c_pt = Mat_c->get_pointer();
	for(int i = 0; i < nrow * ncol; i++) {
		*Mat_c_pt++ = *Mat_m_pt++;
	}
    return Mat_c;
}

void OutputMatrix(mxArray*& Mat_m, SharedMatrix Mat_c) {
	int nrow = Mat_c->nrow();
	int ncol = Mat_c->ncol();
	Mat_m = mxCreateDoubleMatrix( nrow, ncol, mxREAL);
	double* Mat_m_pt = mxGetPr(Mat_m);
	double* Mat_c_pt = Mat_c->get_pointer();
	for(int i = 0; i < nrow * ncol; i++) {
		*Mat_m_pt++ = *Mat_c_pt++;
	}
}

void OutputVector(mxArray*& Mat_m, SharedVector Vec_c) {
	int dim = Vec_c->dim();
	Mat_m = mxCreateDoubleMatrix( dim, 1, mxREAL);
	double* Mat_m_pt = mxGetPr(Mat_m);
	double* Vec_c_pt = Vec_c->pointer();
	for(int i = 0; i < dim; i++) {
		*Mat_m_pt++ = *Vec_c_pt++;
	}
}

void OutputScalar(mxArray*& Mat_m, double scalar) {
	Mat_m = mxCreateDoubleMatrix( 1, 1, mxREAL);
	double* Mat_m_pt = mxGetPr(Mat_m);
    *Mat_m_pt = scalar;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {	
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");
    
    // Constructor
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("Constructor: One output expected.");
        if (nrhs!=2 || !mxIsCell(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 2)
            mexErrMsgTxt("Constructor: Cell array {mol_string, basis_name} input expected.");
        try {
            // Return a handle to a new C++ instance
            mxArray *tmp;
            tmp = mxGetCell(prhs[1], 0);
            std::string mol_string = mxArrayToString(tmp);
            tmp = mxGetCell(prhs[1], 1);
            std::string basis_name = mxArrayToString(tmp);
            plhs[0] = convertPtr2Mat<MatPsi>(new MatPsi(mol_string, basis_name));
        } 
        catch(...) {
            mexErrMsgTxt("Constructor: Cell array {mol_string, basis_name} input expected.");
        }
        return;
    }
    
    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<MatPsi>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }
    
    // Get the class instance pointer from the second input
    MatPsi* MatPsi_obj = convertMat2Ptr<MatPsi>(prhs[1]);
    
    // Call the various class methods 
	// Copy Constructor
	if (!strcmp("MatPsiCopy", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("Copy Constructor: One output expected.");
        try {
            // Return a handle to a new C++ instance
			boost::shared_ptr<MatPsi> MatPsi_str(MatPsi_obj);
            plhs[0] = convertPtr2Mat<MatPsi>(new MatPsi(MatPsi_str));
        } 
        catch(...) {
            mexErrMsgTxt("Copy Constructor failed");
        }
        return;
    }
	
    // Molecule properties
    // natom
    if (!strcmp("natom", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("natom: Unexpected arguments.");
        // Call the method
        OutputScalar(plhs[0], (double)MatPsi_obj->natom());
        return;
    }
    
    // nelec
    if (!strcmp("nelec", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("nelec: Unexpected arguments.");
        // Call the method
        OutputScalar(plhs[0], (double)MatPsi_obj->nelec());
        return;
    }
    
    // coord(natom, 3)
    if (!strcmp("coord", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("coord: Unexpected arguments.");
        // Call the method
		OutputMatrix(plhs[0], MatPsi_obj->coord());
        return;
    }
    
    // Zlist(natom, 1) 
    if (!strcmp("Zlist", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Zlist: Unexpected arguments.");
        // Call the method
		OutputVector(plhs[0], MatPsi_obj->Zlist());
        return;
    }
    
    // Enuc
    if (!strcmp("Enuc", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("Enuc: Unexpected arguments.");
        // Call the method
        OutputScalar(plhs[0], MatPsi_obj->Enuc());
        return;
    }
    
    // Basis set properties 
    // nbasis
    if (!strcmp("nbasis", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("nbasis: Unexpected arguments.");
        // Call the method
        OutputScalar(plhs[0], (double)MatPsi_obj->nbasis());
        return;
    }
    
    // func2center 
    if (!strcmp("func2center", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("func2center: Unexpected arguments.");
        // Call the method
        SharedVector func2centerVec = MatPsi_obj->func2center();
        for(int i = 0; i < func2centerVec->dim(); i++)
            func2centerVec->add(i, 1.0); // + 1 convert C++ convention to Matlab convention 
        OutputVector(plhs[0], func2centerVec);
        return;
    }
    
    // func2am 
    if (!strcmp("func2am", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("func2am: Unexpected arguments.");
        // Call the method
        OutputVector(plhs[0], MatPsi_obj->func2am());
        return;
    }
	
    // One-electron integrals 
	// overlap(nbasis, nbasis)
    if (!strcmp("overlap", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("overlap: Unexpected arguments.");
        // Call the method
		OutputMatrix(plhs[0], MatPsi_obj->overlap());
        return;
    }
    
    // kinetic(nbasis, nbasis)
    if (!strcmp("kinetic", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("kinetic: Unexpected arguments.");
        // Call the method
		OutputMatrix(plhs[0], MatPsi_obj->kinetic());
        return;
    }
    
    // potential(nbasis, nbasis)
    if (!strcmp("potential", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("potential: Unexpected arguments.");
        // Call the method
		OutputMatrix(plhs[0], MatPsi_obj->potential());
        return;
    }
    
    // potential_sep(nbasis, nbasis, natom)
    if (!strcmp("potential_sep", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("potential_sep: Unexpected arguments.");
        // Call the method
        boost::shared_array<SharedMatrix> viMatArray = MatPsi_obj->potential_sep();
        int ncol = viMatArray[0]->ncol();
        int nrow = viMatArray[0]->nrow();
        int natom = MatPsi_obj->natom();
        int dims[3] = {ncol, nrow, natom};
        plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
        double* matlab_pt = mxGetPr(plhs[0]);
        for(int iatom = 0; iatom < natom; iatom++) {
            double* tmp_pt = viMatArray[iatom]->get_pointer();
            for(int i = 0; i < ncol * nrow; i++) {
                matlab_pt[iatom*ncol*nrow + i] = tmp_pt[i];
            }
        }
        return;
    }
    
    // potential_zxyz(nbasis, nbasis)
    if (!strcmp("potential_zxyz", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("potential_zxyz: Unexpected arguments.");
        if (nrhs!=3 || !mxIsCell(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 4)
            mexErrMsgTxt("potential_zxyz: Cell array {Z, x, y, z} input expected.");
        // Call the method
        double Zxyz_array[4];
        for(int i = 0; i < 4; i++) {
            Zxyz_array[i] = (double)mxGetScalar(mxGetCell(prhs[2], i));
        }
		OutputMatrix(plhs[0], MatPsi_obj->potential_zxyz(Zxyz_array));
        return;
    }
    
    // Two-electron integrals 
    // tei_ijkl
    if (!strcmp("tei_ijkl", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("tei_ijkl: Unexpected arguments.");
        if (nrhs!=3 || !mxIsCell(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 4)
            mexErrMsgTxt("tei_ijkl: Cell array {i, j, k, l} input expected.");
        // Call the method
        mxArray *tmp;
        int ind[4];
        for(int i = 0; i < 4; i++) {
            tmp = mxGetCell(prhs[2], i);
            ind[i] = (int)mxGetScalar(tmp) - 1; // -1 convert Matlab convention to C++ convention
            if(ind[i] < 0 || ind[i] >= MatPsi_obj->nbasis())
                mexErrMsgTxt("tei_ijkl: Required index not within scale.");
        }
        OutputScalar(plhs[0], MatPsi_obj->tei_ijkl(ind[0], ind[1], ind[2], ind[3]));
        return;
    }
    
    // tei_alluniq 
    if (!strcmp("tei_alluniq", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("tei_alluniq: Unexpected arguments.");
        // Call the method
        OutputVector(plhs[0], MatPsi_obj->tei_alluniq());
        return;
    }
    
    // tei_alluniqForK 
    if (!strcmp("tei_alluniqForK", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("tei_alluniqForK: Unexpected arguments.");
        // Call the method
        OutputVector(plhs[0], MatPsi_obj->tei_alluniqForK());
        return;
    }
    
    // HFnosymmMO2G(nbasis, nbasis) 
    if (!strcmp("HFnosymmMO2G", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("HFnosymmMO2G: Unexpected arguments.");
        if (nrhs!=3 || !mxIsCell(prhs[2]) || mxGetNumberOfElements(prhs[2]) != 1)
            mexErrMsgTxt("HFnosymmMO2G: Cell array containing a (nbasis by noccupy) matrix {MOmat} input expected.");
        if ( mxGetM(mxGetCell(prhs[2], 0)) != MatPsi_obj->nbasis() )
            mexErrMsgTxt("HFnosymmMO2G: MO matrix dimension does not agree.");
        // Call the method
        mxArray* tmp = mxGetCell(prhs[2], 0);
		OutputMatrix(plhs[0], MatPsi_obj->HFnosymmMO2G(InputMatrix(tmp)));
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

