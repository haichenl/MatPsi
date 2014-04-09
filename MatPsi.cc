// Please directly include this file in MatPsi_mex.cpp to make compilation easier... 

// Constructor
MatPsi::MatPsi(std::string mol_string, std::string basis_name, int ncores, unsigned long int memory) {
    // necessary initializations
    Process::environment.initialize();
    WorldComm = initialize_communicator(0, NULL);
    Process::environment.options.set_read_globals(true);
    read_options("", Process::environment.options, true);
    Process::environment.options.set_read_globals(false);
    
    // set cores and memory 
    Process::environment.set_n_threads(ncores);
    Process::environment.set_memory(memory);
    
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
    
    // initialize Hartree-Fock energy 
    ERHF_ = 0.0;
    
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

int MatPsi::tei_uniqN() {
    return ( basis_->nbf() * ( basis_->nbf() + 1 ) * ( basis_->nbf() * basis_->nbf() + basis_->nbf() + 2 ) ) / 8;
}

void MatPsi::tei_alluniq(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nuniq = tei_uniqN();
    const double *buffer = eri_->buffer();
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // Compute quartet
        eri_->compute_shell(shellIter);
        // From the quartet get all the integrals
        AOIntegralsIterator intIter = shellIter.integrals_iterator();
        for (intIter.first(); intIter.is_done() == false; intIter.next()) {
            matpt[ ij2I( ij2I(intIter.i(), intIter.j()), ij2I(intIter.k(), intIter.l()) ) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi::tei_allfull(double* matpt) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
    int nbf_ = basis_->nbf();
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
            matpt[ l+nbf_*(k+nbf_*(j+nbf_*i)) ] = buffer[intIter.index()];
            matpt[ l+nbf_*(k+nbf_*(i+nbf_*j)) ] = buffer[intIter.index()];
            matpt[ k+nbf_*(l+nbf_*(j+nbf_*i)) ] = buffer[intIter.index()];
            matpt[ k+nbf_*(l+nbf_*(i+nbf_*j)) ] = buffer[intIter.index()];
            matpt[ j+nbf_*(i+nbf_*(l+nbf_*k)) ] = buffer[intIter.index()];
            matpt[ j+nbf_*(i+nbf_*(k+nbf_*l)) ] = buffer[intIter.index()];
            matpt[ i+nbf_*(j+nbf_*(l+nbf_*k)) ] = buffer[intIter.index()];
            matpt[ i+nbf_*(j+nbf_*(k+nbf_*l)) ] = buffer[intIter.index()];
        }
    }
}

void MatPsi::tei_alluniqJK(double* matptJ, double* matptK) {
    AOShellCombinationsIterator shellIter = intfac_->shells_iterator();
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
            matptJ[ij2I( ij2I(i, j), ij2I(k, l) )] = buffer[intIter.index()];
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
            matptK[ij2I( ij2I(i, j), ij2I(k, l) )] = matptJ[ij2I( ij2I(i, l), ij2I(k, j) )] + matptJ[ij2I( ij2I(i, k), ij2I(j, l) )];
        }
    }
}

void MatPsi::init_directjk(double cutoff) {
    directjk_->set_cutoff(cutoff);
    directjk_->initialize();
    directjk_->remove_symmetry();
}

SharedMatrix MatPsi::HFnosymmMO2G(SharedMatrix MO) {
    init_directjk();
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    C_left.clear();
    C_left.push_back(MO);
    directjk_->compute();
    SharedMatrix Gnew = directjk_->J()[0];
    SharedMatrix Knew = directjk_->K()[0];
    Knew->scale(-0.5);
    Gnew->add(Knew);
    directjk_->finalize();
    return Gnew;
}

void MatPsi::form_density()
{
    for(int p = 0; p < nbf_; ++p){
        for(int q = 0; q < nbf_; ++q){
            double val = 0.0;
            for(int i = 0; i < ndocc_; ++i){
                val += C_occ_->get(p, i) * C_occ_->get(q, i);
            }
            D_->set(p, q, val);
        }
    }   
}

double MatPsi::compute_electronic_energy()
{
    Matrix HplusF;
    HplusF.copy(H_);
    HplusF.add(F_);
    return D_->vector_dot(HplusF);    
}

double MatPsi::DirectRHF() {

    // Initializations 
    maxiter_ = 100;
    e_convergence_ = 1e-8;
    d_convergence_ = 1e-8;
    
    nbf_ = basis_->nbf();
    int nelec_ = nelec();
    if(nelec_ % 2)
        throw PSIEXCEPTION("This is only an RHF code, but you gave it an odd number of electrons.  Try again!");
    ndocc_ = nelec_ / 2;
    e_nuc_ = molecule_->nuclear_repulsion_energy();
    
    S_ = SharedMatrix(new Matrix("Overlap matrix", nbf_, nbf_));
    H_ = SharedMatrix(new Matrix("Core Hamiltonian matrix", nbf_, nbf_));
    boost::shared_ptr<OneBodyAOInt> sOBI(intfac_->ao_overlap());
    boost::shared_ptr<OneBodyAOInt> tOBI(intfac_->ao_kinetic());
    boost::shared_ptr<OneBodyAOInt> vOBI(intfac_->ao_potential());
    SharedMatrix V = SharedMatrix(new Matrix("Potential integrals matrix", nbf_, nbf_));
    sOBI->compute(S_);
    tOBI->compute(H_);
    vOBI->compute(V);
    H_->add(V);
    
    // Allocate some matrices
    X_  = SharedMatrix(new Matrix("S^-1/2", nbf_, nbf_));
    F_  = SharedMatrix(new Matrix("Fock matrix", nbf_, nbf_));
    Ft_ = SharedMatrix(new Matrix("Transformed Fock matrix", nbf_, nbf_));
    C_  = SharedMatrix(new Matrix("MO Coefficients_", nbf_, nbf_));
    Eorb_  = SharedVector(new Vector("MO Energies_", nbf_));
    C_occ_  = SharedMatrix(new Matrix("Occupied MO Coefficients_", nbf_, ndocc_));
    D_  = SharedMatrix(new Matrix("The Density Matrix", nbf_, nbf_));
    SharedMatrix Temp1(new Matrix("Temporary Array 1", nbf_, nbf_));
    SharedMatrix Temp2(new Matrix("Temporary Array 2", nbf_, nbf_));
    SharedMatrix FDS(new Matrix("FDS", nbf_, nbf_));
    SharedMatrix SDF(new Matrix("SDF", nbf_, nbf_));
    SharedMatrix Evecs(new Matrix("Eigenvectors", nbf_, nbf_));
    SharedVector Evals(new Vector("Eigenvalues", nbf_));

    // Form the X_ matrix (S^-1/2)
    S_->diagonalize(Evecs, Evals, ascending);
    for(int p = 0; p < nbf_; ++p){
        double val = 1.0 / sqrt(Evals->get(p));
        Evals->set(p, val);
    }
    Temp1->set_diagonal(Evals);
    Temp2->gemm(false, true, 1.0, Temp1, Evecs, 0.0);
    X_->gemm(false, false, 1.0, Evecs, Temp2, 0.0);
    
    F_->copy(H_);
    Ft_->transform(F_, X_);
    Ft_->diagonalize(Evecs, Evals, ascending);

    C_occ_->gemm(false, false, 1.0, X_, Evecs, 0.0);
    form_density();

    int iter = 1;
    bool converged = false;
    double e_old;
    double e_new = compute_electronic_energy();
    
    init_directjk();
    std::vector<SharedMatrix>& C_left = directjk_->C_left();
    
    while(!converged && iter < maxiter_){
        e_old = e_new;

        // Add the core Hamiltonian term to the Fock operator
        F_->copy(H_);
        
        // Add the two electron terms to the Fock operator
        C_left.clear();
        C_left.push_back(C_occ_);
        directjk_->compute();
        Temp1 = directjk_->J()[0];
        Temp1->scale(2);
        F_->add(Temp1);
        F_->subtract(directjk_->K()[0]);

        // Transform the Fock operator and diagonalize it
        Ft_->transform(F_, X_);
        Ft_->diagonalize(Evecs, Evals, ascending);

        // Form the orbitals from the eigenvectors of the transformed Fock matrix 
        C_occ_->gemm(false, false, 1.0, X_, Evecs, 0.0);

        // Rebuild the density using the new orbitals
        form_density();

        // Compute the energy
        e_new = compute_electronic_energy();
        double dE = e_new - e_old;
 
        // Compute the orbital gradient, FDS-SDF
        Temp1->gemm(false, false, 1.0, D_, S_, 0.0);
        FDS->gemm(false, false, 1.0, F_, Temp1, 0.0);
        Temp1->gemm(false, false, 1.0, D_, F_, 0.0);
        SDF->gemm(false, false, 1.0, S_, Temp1, 0.0);
        Temp1->copy(FDS);
        Temp1->subtract(SDF);
        double dRMS = Temp1->rms();
     
        converged = (fabs(dE) < e_convergence_) && (dRMS < d_convergence_);

        iter++;
    }
    directjk_->finalize();
    C_->gemm(false, false, 1.0, X_, Evecs, 0.0);
    Eorb_->copy(Evals.get());
    
    if(!converged)
        throw PSIEXCEPTION("The SCF iterations did not converge.");
    
    ERHF_ = e_nuc_ + e_new;
    return ERHF_;
}

double MatPsi::ERHF() { 
    if(F_ == NULL) { // a trick to determine whether Hartree-Fock has been performed 
        throw PSIEXCEPTION("ERHF: Hartree-Fock calculation has not been done.");
    }
    return ERHF_; 
}

SharedMatrix MatPsi::orbital() { 
    if(C_ == NULL) {
        throw PSIEXCEPTION("orbital: Hartree-Fock calculation has not been done.");
    }
    return C_; 
}

SharedVector MatPsi::Eorb() { 
    if(Eorb_ == NULL) {
        throw PSIEXCEPTION("Eorb: Hartree-Fock calculation has not been done.");
    }
    return Eorb_; 
}

SharedMatrix MatPsi::density() { 
    if(D_ == NULL) {
        throw PSIEXCEPTION("density: Hartree-Fock calculation has not been done.");
    }
    return D_; 
}

SharedMatrix MatPsi::H1Matrix() { 
    if(H_ == NULL) {
        throw PSIEXCEPTION("H1Matrix: Hartree-Fock calculation has not been done.");
    }
    return H_; 
}

SharedMatrix MatPsi::FockMatrix() { 
    if(F_ == NULL) {
        throw PSIEXCEPTION("FockMatrix: Hartree-Fock calculation has not been done.");
    }
    return F_; 
}

