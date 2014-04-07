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
    
    // create basis object; use pure angular momentum functions (5d or 7f etc.) forever 
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    parser->force_puream_or_cartesian_ = true;
    parser->forced_is_puream_ = true;
    basis_ = BasisSet::construct(parser, molecule_, "BASIS");
    
    // create integral factory object 
    intfac_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basis_, basis_, basis_, basis_));
    
    // create two electron integral generator
    eri_ = boost::shared_ptr<TwoBodyAOInt>(intfac_->eri());
    
    // create matrix factory object 
    int nbf[] = { basis_->nbf() };
    matfac_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory);
    matfac_->init_with(1, nbf, nbf);
    
    // create directJK object
    directjk_ = boost::shared_ptr<DirectJK>(new DirectJK(basis_));
    
}

// copy constructor 
MatPsi::MatPsi(boost::shared_ptr<MatPsi> inputMatPsi) {
    molecule_ = inputMatPsi->molecule_;
    basis_ = inputMatPsi->basis_;
    intfac_ = inputMatPsi->intfac_;
    eri_ = inputMatPsi->eri_;
    matfac_ = inputMatPsi->matfac_;
    directjk_ = inputMatPsi->directjk_;
}

SharedVector MatPsi::Zlist() {
    SharedVector zlistvec(new Vector(molecule_->natom()));
    for(int i = 0; i < molecule_->natom(); i++) {
        zlistvec->set(i, (double)molecule_->Z(i));
    }
    return zlistvec;
}

int MatPsi::nelec() {
    int charge = molecule_->molecular_charge();
    int nelectron  = 0;
    for(int i = 0; i < molecule_->natom(); i++)
        nelectron += (int)molecule_->Z(i);
    nelectron -= charge;
    return nelectron;
}

SharedVector MatPsi::func2center() {
    SharedVector func2centerVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2centerVec->set(i, (double)basis_->function_to_center(i));
    }
    return func2centerVec;
}

SharedVector MatPsi::func2am() {
    SharedVector func2amVec(new Vector(basis_->nbf()));
    for(int i = 0; i < basis_->nbf(); i++) {
        func2amVec->set(i, (double)basis_->shell(basis_->function_to_shell(i)).am());
    }
    return func2amVec;
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
    int natom_ = molecule_->natom();
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

SharedMatrix MatPsi::potential_zxyz(const double* Zxyz_array) {
    boost::shared_ptr<OneBodyAOInt> viOBI(intfac_->ao_potential());
    boost::shared_ptr<PotentialInt> viPtI = boost::static_pointer_cast<PotentialInt>(viOBI);
    SharedMatrix Zxyz_row(new Matrix(1, 4));
    for(int i = 0; i < 4; i++)
        Zxyz_row->set(0, i, Zxyz_array[i]);
    viPtI->set_charge_field(Zxyz_row);
    SharedMatrix vZxyzMat(matfac_->create_matrix("Potential_Zxyz"));
    viOBI->compute(vZxyzMat);
    return vZxyzMat;
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

inline int ij2I(int i, int j) {
    if(i < j) {
        int tmp = i;
        i = j;
        j = tmp;
    }
    return i * ( i + 1 ) / 2 + j;
}

SharedVector MatPsi::tei_alluniq() {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbasis_ = basis_->nbf();
    int nuniq = ( nbasis_ * ( nbasis_ + 1 ) * ( nbasis_ * nbasis_ + nbasis_ + 2 ) ) / 8;
    SharedVector tei_alluniqvec(new Vector(nuniq));
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            tei_alluniqvec->set( ij2I( ij2I(intIter.i(), intIter.j()), ij2I(intIter.k(), intIter.l()) ), buffer[intIter.index()] );
        }
    }
    return tei_alluniqvec;
}

boost::shared_array<SharedVector> MatPsi::tei_alluniqJK() {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbasis_ = basis_->nbf();
    int nuniq = ( nbasis_ * ( nbasis_ + 1 ) * ( nbasis_ * nbasis_ + nbasis_ + 2 ) ) / 8;
    boost::shared_array<SharedVector> JKvecs( new SharedVector [2] );
    JKvecs[0] = SharedVector(new Vector(nuniq));
    JKvecs[1] = SharedVector(new Vector(nuniq));
    SharedVector auxvec(new Vector(nuniq));
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            JKvecs[0]->set(ij2I( ij2I(i, j), ij2I(k, l) ), buffer[intIter.index()] );
        }
    }
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            int i = intIter.i();
            int j = intIter.j();
            int k = intIter.k();
            int l = intIter.l();
            JKvecs[1]->set( ij2I( ij2I(i, j), ij2I(k, l) ), JKvecs[0]->get( ij2I( ij2I(i, l), ij2I(k, j) ) ) );
            auxvec->set( ij2I( ij2I(i, j), ij2I(k, l) ), JKvecs[0]->get( ij2I( ij2I(i, k), ij2I(j, l) ) ) );
        }
    }
    JKvecs[1]->add(auxvec);
    return JKvecs;
}

SharedMatrix MatPsi::HFnosymmMO2G(SharedMatrix MO, long int memory, double cutoff) {
    directjk_->set_memory(memory);
    directjk_->set_cutoff(cutoff);
    directjk_->set_do_J(true);
    directjk_->set_do_K(true);
    directjk_->set_allow_desymmetrization(true);
    directjk_->initialize();
    directjk_->remove_symmetry();
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    C_left.push_back(MO);
    directjk_->compute();
    SharedMatrix Gnew = directjk_->J()[0];
    SharedMatrix Knew = directjk_->K()[0];
    Knew->scale(0.5);
    Gnew->add(Knew);
    directjk_->finalize();
    return Gnew;
}