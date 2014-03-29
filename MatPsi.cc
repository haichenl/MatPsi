// Please directly include this file in MatPsi_mex.cpp to make compilation easier... 

// Constructor
MatPsi::MatPsi(std::string mol_string, std::string basis_name) {
    // necessary initializations
    Process::environment.initialize();
    WorldComm = initialize_communicator(0, NULL);
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);
    Wavefunction::initialize_singletons();
    
    // create molecule object and set its basis set name 
    molecule_ = psi::Molecule::create_molecule_from_string(mol_string);
    molecule_->set_basis_all_atoms(basis_name);
    
    // create basis object 
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    basis_ = BasisSet::construct(parser, molecule_, "BASIS");
    
    // create integral factory object 
    boost::shared_ptr<IntegralFactory> integral_factory_temp(new IntegralFactory(basis_, basis_, basis_, basis_));
    intfac_ = integral_factory_temp;
    
    // create two electron integral generator
    boost::shared_ptr<TwoBodyAOInt> eri_temp(intfac_->eri());
    eri_ = eri_temp;
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    boost::shared_ptr<MatrixFactory> matrix_factory_temp(new MatrixFactory);
    matfac_ = matrix_factory_temp;
    matfac_->init_with(1, nbf, nbf);
    
    // number of atoms 
    natom_ = molecule_->natom();
    
    // number of basis functions 
    nbasis_ = *nbf;
    
    // number of electrons 
    int charge = molecule_->molecular_charge();
    int nelectron  = 0;
    for(int i = 0; i < natom_; i++)
        nelectron += (int)molecule_->Z(i);
    nelectron -= charge;
    nelec_ = nelectron;
}

// copy constructor 
MatPsi::MatPsi(boost::shared_ptr<MatPsi> inputMatPsi) {
    natom_ = inputMatPsi->natom_;
    nbasis_ = inputMatPsi->nbasis_;
    nelec_ = inputMatPsi->nelec_;
    molecule_ = inputMatPsi->molecule_;
    basis_ = inputMatPsi->basis_;
    intfac_ = inputMatPsi->intfac_;
    eri_ = inputMatPsi->eri_;
    matfac_ = inputMatPsi->matfac_;
}

SharedMatrix MatPsi::overlap() {
    SharedMatrix sMat(matfac_->create_matrix("Overlap"));
    boost::shared_ptr<OneBodyAOInt> sOBI(intfac_->ao_overlap());
    sOBI->compute(sMat);
    return sMat;
}

SharedMatrix MatPsi::kinetic() {
    SharedMatrix tMat(matfac_->create_matrix("Kinetic"));
    boost::shared_ptr<OneBodyAOInt> tOBI(intfac_->ao_kinetic());
    tOBI->compute(tMat);
    return tMat;
}

SharedMatrix MatPsi::potential() {
    SharedMatrix vMat(matfac_->create_matrix("Potential"));
    boost::shared_ptr<OneBodyAOInt> vOBI(intfac_->ao_potential());
    vOBI->compute(vMat);
    return vMat;
}

boost::shared_array<SharedMatrix> MatPsi::potential_sep() {
    boost::shared_array<SharedMatrix> viMatArray(new SharedMatrix [natom_]);
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    SharedMatrix Zxyz = viPtI->charge_field();
    SharedMatrix Zxyz_rowi(new Matrix(1, 4));
    for( int i = 0; i < natom_; i++) {
        SharedVector Zxyz_rowi_vec = Zxyz->get_row(0, i);
        Zxyz_rowi->set_row(0, 0, Zxyz_rowi_vec);
        viPtI->set_charge_field(Zxyz_rowi);
        viMatArray[i] = matfac_->create_shared_matrix("potential_sep");
        viOBI->compute(viMatArray[i]);
    }
    return viMatArray;
}

double MatPsi::tei_ijkl(int i, int j, int k, int l) {
    int ish = basis_->function_to_shell(i);
    int jsh = basis_->function_to_shell(j);
    int ksh = basis_->function_to_shell(k);
    int lsh = basis_->function_to_shell(l);
    int ii = i - basis_->shell_to_basis_function(ish);
    int jj = j - basis_->shell_to_basis_function(jsh);
    int kk = k - basis_->shell_to_basis_function(ksh);
    int ll = l - basis_->shell_to_basis_function(lsh);
    int ni = basis_->shell(ish).nfunction();
    int nj = basis_->shell(jsh).nfunction();
    int nk = basis_->shell(ksh).nfunction();
    int nl = basis_->shell(lsh).nfunction();
    eri_->compute_shell(ish, jsh, ksh, lsh);
    const double *buffer = eri_->buffer();
    return buffer[ll+nl*(kk+nk*(jj+nj*ii))];
}