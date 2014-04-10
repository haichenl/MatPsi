MatPsi: An interface between Matlab and Psi4
======

##Compilation

Please refer to complink.sh 

##Usage 

Usage is (nearly) all of Matlab convention, except that all the function __input arguments__ have to be in __one cell array__: 

    >> [output1, output2, ...] = matpsi.WhatEverFunction( {input1, input2, ...} );

This is assumed by _Example MATLAB class wrapper for a C++ class_ 's developer Oliver Woodford. 
See http://www.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class for more details. 

####Constructor 

Usually, we construct a MatPsi object using 2 strings, one describes the molecule's geometry (as well as charge and spin), 
and one sets the name of the basis set we use. 

```
mol_string = 

O
H 1 R
H 1 R 2 A

R = .9
A = 104.5

basis_name = 

6-31g*

>> matpsi = MatPsi( {mol_string, basis_name} );
```

Notice that both strings are in a cell array. 

####Copy Constructor

Construct from an existing MatPsi object. 

    >> matpsi2 = matpsi.MatPsiCopy();

####Molecule properties 

1. natom: Number of atoms. 

    ```
    >> matpsi.natom(); 
    ```

2. nelec: Number of electrons. 

    ```
    >> matpsi.nelec(); 
    ```

3. coord: Atom Cartesian coordinates. 

    ```
    >> matpsi.coord(); 
    ```

4. Enuc:  Nuclear repulsion energy. 

    ```
    >> matpsi.Enuc(); 
    ```

####Basis set properties 

1. nbasis: Number of basis functions. 

    ```
    >> matpsi.nbasis(); 
    ```

2. func2center: A vector mapping basis function number to the number of atom it is centred on, a.k.a "basisAtom" in MSQC. 

    ```
    >> matpsi.func2center(); 
    ```

3. func2am(): A vector mapping basis function number to the corresponding angular momentum number, a.k.a "basisType" in MSQC. 

    ```
    >> matpsi.func2am(); 
    ```

####One-electron integrals 

1. overlap: Atomic orbital overlap matrix (S). 

    ```
    >> matpsi.overlap(); 
    ```

2. kinetic: Kinetic energy matrix (KE). 

    ```
    >> matpsi.kinetic(); 
    ```

3. potential: One-electron potential energy matrix (EN). Summed over all atoms. 

    ```
    >> matpsi.potential(); 
    ```

4. potential_sep: Atom-separated One-electron potential energy matrices (ENI). Return a 3-D, (nbasis by nbasis by natom) array. 

    ```
    >> matpsi.potential_sep(); 
    ```

5. potential_zxyz: Environment potential energy matrix (ENVI) for a given point charge in the format of {Z, x, y, z}; Z stands for 
charge magnitude; x, y, z are Cartesian coordinates of the point charge. 

    ```
    >> matpsi.potential_zxyz( {Z, x, y, z} ); 
    ```

####Two-electron integrals 

1. tei_ijkl: 4-indexed two-electron interaction integral. Need four indices as input arguments. Return only one integral value. 

    ```
    >> matpsi.tei_ijkl( {i, j, k, l} ); 
    ```

    Notice that all indices are in a cell array. 

2. tei_uniqN: Total number of unique two-electron interaction integrals. No geometrical symmetry is considered. 

    ```
    >> matpsi.tei_uniqN(); 
    ```

3. tei_alluniq: All unique two-electron interaction integrals. No geometrical symmetry is considered. Be careful as it consumes 
a huge amount of memory. 

    ```
    >> matpsi.tei_alluniq(); 
    ```

4. tei_alluniqJK: Still all unique two-electron interaction integrals but ordered and summed in quick forming of 
the exchange energy matrix K. See Developer's Note for detailed discussions. 

    ```
    >> [Jvec, Kvec] = matpsi.tei_alluniqJK(); 
    ```

5. tei_allfull: Full, 4-D two-electron interaction integrals. No geometrical symmetry is considered. Be careful as it consumes 
a really huge amount of memory. 

    ```
    >> matpsi.tei_allfull(); 
    ```

####SCF related 

1. OccMO2J: For restricted Hartree Fock theory, form the 2-electron Coulomb interaction J matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2J( { OccMO } ); 
    ```
    
    Notice that OccMO is in a cell array and must have nbasis number of rows. 

2. OccMO2K: For restricted Hartree Fock theory, form the 2-electron exchange interaction K matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2K( { OccMO } ); 
    ```
    
    Notice that OccMO is in a cell array and must have nbasis number of rows. 

3. OccMO2G: For restricted Hartree Fock theory, form the 2-electron G = 2 * J - K matrix from occupied molecular orbital coefficients matrix. 
Direct algorithm. No geometrical symmetry is considered. 

    ```
    >> matpsi.OccMO2G( { OccMO } ); 
    ```
    
    Notice that OccMO is in a cell array and must have nbasis number of rows. 

4. DirectRHF: Solve restricted Hartree-Fock functions and returns the final Hartree-Fock energy 

    ```
    >> matpsi.DirectRHF(); 
    ```

5. ERHF: Get the final restricted Hartree-Fock energy. Executable after DirectRHF. 

    ```
    >> matpsi.ERHF(); 
    ```

6. orbital: Restricted Hartree-Fock molecular orbital coefficients. Executable after DirectRHF. 

    ```
    >> matpsi.orbital(); 
    ```

7. Eorb: Restricted Hartree-Fock molecular orbital energies (eigenvalues). Executable after DirectRHF. 

    ```
    >> matpsi.Eorb(); 
    ```

8. density: Restricted Hartree-Fock density matrix. Executable after DirectRHF. 

    ```
    >> matpsi.density(); 
    ```

9. H1Matrix: One-electron (core) Hamiltonian matrix. Executable after DirectRHF. 

    ```
    >> matpsi.H1Matrix(); 
    ```

10. JMatrix: 2-electron Coulomb interaction J matrix. Executable after DirectRHF. 

    ```
    >> matpsi.JMatrix(); 
    ```

11. KMatrix: 2-electron exchange interaction K matrix. Executable after DirectRHF. 

    ```
    >> matpsi.KMatrix(); 
    ```

12. FockMatrix: Fock matrix. Executable after DirectRHF. 

    ```
    >> matpsi.FockMatrix(); 
    ```

##TODO 



##Developer's Note 

Mar. 29: To make two-electron integrals work, the original Psi4 source code has been slightly changed. Below is where and why. 

From "$psi4dir/src/lib/libmins/twobody.h", we know that the virtual class `TwoBodyAOInt` has some properties probably related with 
Psi4's python interface: 

    class TwoBodyAOInt
    {
    protected:
    ...
        /// The PyBuffer object used for sharing the target_ buffer without copying data
        PyBuffer<double> target_pybuffer_;
        /// Whether or not to use the PyBuffer
        bool enable_pybuffer_;
    ...
    }; // class end

and in "$psi4dir/src/lib/libmins/twobody.cc", `TwoBodyAOInt` class virtual constructor: 

    TwoBodyAOInt::TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv) :
        integral_(intsfactory),
        original_bs1_(integral_->basis1()),
        original_bs2_(integral_->basis2()),
        original_bs3_(integral_->basis3()),
        original_bs4_(integral_->basis4()),
        deriv_(deriv),
        target_(0),
        target_pybuffer_(&target_, true)

the last line shown, `target_pybuffer_(&target_, true)`, is known causing some strange segmentation fault errors. 
Eliminating this line fixes it, but as a result, we probably need to set `bool enable_pybuffer_` as `false` forever. 

Apr. 02: A new method, void remove_symmetry(), has been added to JK class. It allows us to get rid of the geometrical symmetry automatically 
imposed (but not used in real computation for some weird reasons) by JK constructor. 

    class JK {
        ...
    public:
        ...
        void remove_symmetry() {
            AO2USO_ = SharedMatrix(new Matrix("AO->SO matrix", primary_->nbf(), primary_->nbf()));
            AO2USO_->identity();
        }
    };

Apr. 05: About forming the exchange energy matrix K 

In order to (efficiently) form K, we loop over all unique two-electron integrals and store them as: 

```
I = i*(i-1)/2 + j;
J = k*(k-1)/2 + l;
Jvec( I*(I-1)/2 + J ) = H2(i, j, k, l);
Kvec( I*(I-1)/2 + J ) = H2(i, l, k, j) + H2(i, k, j, l);
```

where i >= j, k >= l, and I >= J. Then we can form K from teiVecForK and the density matrix in the same manner as forming J, that is, 
restore this vector into an (upper case) N by N (N = n*(n+1)/2, n is the number of basis functions) Hermitian matrix, contract 
the density matrix into an N-length column vector and scale it by C(J) = 2 or C(K) = -1/2 and scale again all of the diagonal elements by 1/2, 
then time the N by N matrix by this density vector. In the real algorithm we can avoid the "restore", or Hermitianize, step and 
save some memory as well as sacrifice some speed, though. 


