MatPsi: An interface between Matlab and Psi4
======

#Usage 

Usage is (nearly) all of Matlab convention, except that all the function __input arguments__ have to be in __one cell array__: 

    >> [output1, output2, ...] = matpsi.WhatEverFunction( {input1, input2, ...} );

This is assumed by _Example MATLAB class wrapper for a C++ class_ 's developer Oliver Woodford. See http://www.mathworks.com/matlabcentral/fileexchange/38964-example-matlab-class-wrapper-for-a-c++-class for more details. 

###Constructor 

Usually, we construct a MatPsi object using 2 strings, one describing the molecule's geometry, and one set the name of the basis set we are going to use. 

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

###Copy Constructor

Construct from an existing MatPsi object. 

    >> matpsi2 = matpsi.MatPsiCopy();

###Molecule and basis set properties 

1. natom: number of atoms. 

    ```
    >> matpsi.natom(); 
    ```

2. nbasis: number of basis functions. 

    ```
    >> matpsi.nbasis(); 
    ```

3. nelec: number of electrons. 

    ```
    >> matpsi.nelec(); 
    ```

###One electron integrals 

1. overlap: atomic orbital overlap matrix (S). 

    ```
    >> matpsi.overlap(); 
    ```

2. kinetic: kinetic energy matrix (KE). 

    ```
    >> matpsi.kinetic(); 
    ```

3. potential: 1-electron potential energy matrix (EN). Summed over all atoms. 

    ```
    >> matpsi.potential(); 
    ```

4. potential_sep: atom-separated 1-electron potential energy matrices (ENI). Returns a 3-D, (nbasis by nbasis by natom) array. 

    ```
    >> matpsi.potential_sep(); 
    ```

###Two electron integrals 

1. tei_ijkl: 4-indexed two electron interaction integral. Needs four indices as input arguments. Returns only one integral value. 

    ```
    >> matpsi.tei_ijkl( {i, j, k, l} ); 
    ```

    Notice that all indices are in a cell array. 

#TODO 

1. Add some more two electron integral methods. 

2. Hartree Fock SCF. 

3. Add a function that takes a density matrix as input and outputs the G = J - 1/2 * K matrix. 

4. Environment. Need to know exactly how we implement. 

#Developer's Note 

To make two electron integrals work, the original Psi4 source code has been slightly changed. Below is where and why. 

From "$psi4dir/src/lib/libmins/twobody.h", we know that the virtual class `TwoBodyAOInt` has the properties below which are probably related with their python interface; 

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

the last line, `target_pybuffer_(&target_, true)`, is known causing some strange segmentation fault errors. Eliminating this line fixes it, but as a result, we probably need to set `bool enable_pybuffer_` as `false` for ever. 




