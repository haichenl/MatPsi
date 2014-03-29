#include "mex.h"
#include "class_handle.hpp"
#include "MatPsi.h"
#include "MatPsi.cc"

using namespace std;
using namespace psi;
using namespace boost;

void MatOut(SharedMatrix Mat_c, mxArray*& Mat_m) {
	int nrow = Mat_c->nrow();
	int ncol = Mat_c->ncol();
	Mat_m = mxCreateDoubleMatrix( nrow, ncol, mxREAL);
	double* Mat_m_pt = mxGetPr(Mat_m);
	double* Mat_c_pt = Mat_c->get_pointer();
	for(int i = 0; i < nrow * ncol; i++) {
		*Mat_m_pt++ = *Mat_c_pt++;
	}
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
	
    // natom
    if (!strcmp("natom", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("natom: Unexpected arguments.");
        // Call the method
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        double* natom_pt = mxGetPr(plhs[0]);
        *natom_pt = (double)MatPsi_obj->natom();
        return;
    }
    
    // nbasis
    if (!strcmp("nbasis", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("nbasis: Unexpected arguments.");
        // Call the method
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        double* nbasis_pt = mxGetPr(plhs[0]);
        *nbasis_pt = (double)MatPsi_obj->nbasis();
        return;
    }
    
    // nelec
    if (!strcmp("nelec", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("nelec: Unexpected arguments.");
        // Call the method
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        double* nelec_pt = mxGetPr(plhs[0]);
        *nelec_pt = (double)MatPsi_obj->nelec();
        return;
    }
    
    // coord(natom, 3)
    if (!strcmp("coord", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("coord: Unexpected arguments.");
        // Call the method
        SharedMatrix zMat = MatPsi_obj->coord();
		MatOut(zMat, plhs[0]);
        return;
    }
	
	// overlap(nbasis, nbasis)
    if (!strcmp("overlap", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("overlap: Unexpected arguments.");
        // Call the method
        SharedMatrix sMat = MatPsi_obj->overlap();
		MatOut(sMat, plhs[0]);
        return;
    }
    
    // kinetic(nbasis, nbasis)
    if (!strcmp("kinetic", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("kinetic: Unexpected arguments.");
        // Call the method
        SharedMatrix tMat = MatPsi_obj->kinetic();
		MatOut(tMat, plhs[0]);
        return;
    }
    
    // potential(nbasis, nbasis)
    if (!strcmp("potential", cmd)) {
        // Check parameters
        if (nlhs < 0 || nrhs < 2)
            mexErrMsgTxt("potential: Unexpected arguments.");
        // Call the method
        SharedMatrix vMat = MatPsi_obj->potential();
		MatOut(vMat, plhs[0]);
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
        int nbasis = MatPsi_obj->nbasis();
        for(int i = 0; i < 4; i++) {
            tmp = mxGetCell(prhs[2], i);
            ind[i] = (int)mxGetScalar(tmp) - 1; // -1 convert Matlab convention to C++ convention
            if(ind[i] < 0 || ind[i] >= nbasis)
                mexErrMsgTxt("tei_ijkl: Required index not within scale.");
        }
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        double* tei_ijkl_pt = mxGetPr(plhs[0]);
        *tei_ijkl_pt = MatPsi_obj->tei_ijkl(ind[0], ind[1], ind[2], ind[3]);
        return;
    }
    
    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}

