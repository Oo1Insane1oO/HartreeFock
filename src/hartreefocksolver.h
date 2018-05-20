#ifndef HARTREEFOCKSOLVER_H
#define HARTREEFOCKSOLVER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <experimental/filesystem>

#include <Eigen/Dense>
#include <yaml-cpp/yaml.h>
#include <mpi.h>

#include <string>

template<class Integrals>
class HartreeFockSolver {
    private:
        Integrals* m_I;

        unsigned int totalSize, iterations;

        int myRank, numProcs;

        double energy;

        static constexpr double mixingFactor = 0.6;

        bool interaction;

        std::string dirpath;

        Eigen::ArrayXd twoBodyNonAntiSymmetrizedElements;
        Eigen::MatrixXd oneBodyElements, overlapElements;

        Eigen::MatrixXd FockMatrix, densityMatrix, coefficients;
        
        inline void setDensityMatrix() {
            /* set density matrix in HartreeFock */
            for (unsigned int c = 0; c < coefficients.rows(); ++c) {
                for (unsigned int d = 0; d < coefficients.cols(); ++d) {
                    densityMatrix(c,d) = 0;
                    for (unsigned int i = 0; i < m_numParticles/2; ++i) {
                        densityMatrix(c,d) += coefficients(c,i) *
                            coefficients(d,i);
                    } // end fori
                } // end ford
            } // end forc
        } // end function setDensityMatrix

        inline void setFockMatrix() {
            /* set Hartree-Fock matrix */
            FockMatrix.setZero();
            for (unsigned int p = 0; p < m_numStates; ++p) {
                for (unsigned int q = p; q < m_numStates; ++q) {
                    FockMatrix(p,q) = oneBodyElements(p,q);
                    for (unsigned int r = 0; r < m_numStates; ++r) {
                        for (unsigned int s = 0; s < m_numStates; ++s) {
                            FockMatrix(p,q) += densityMatrix(r,s) * (2 *
                                    twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates,
                                            p,r,q,s)) -
                                    twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates,
                                            p,r,s,q)));
                        } // end fors
                    } // end forr

                    // matrix is symmetric by definition
                    FockMatrix(q,p) = FockMatrix(p,q);
                } // end forq
            } // end forp
        } // end function sethartreeFockMatrix

        unsigned int dIndex(const unsigned int& N, const unsigned int& i, const
                unsigned int& j, const unsigned int& k, const unsigned int& l)
            const {
            /* calculate offset for 4d-matrix (square case) for indices
             * (i,j,k,l) */
            return i + N * (j + N * (k + N*l));
        } // end function dIndex

        void setNonAntiSymmetrizedElements(const Eigen::ArrayXd& pqrsElements,
                const Eigen::ArrayXd& pqsrElements) {
            /* set symmetric values in full two-body matrix used in
             * Hartree-Fock algorithm */
            int pqrs = 0;
            for (unsigned int p = 0; p < m_numStates; ++p)
            for (unsigned int q = p; q < m_numStates; ++q)
            for (unsigned int r = 0; r < m_numStates; ++r)
            for (unsigned int s = r; s < m_numStates; ++s)
            {
                double value = pqrsElements(pqrs);
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, p,q,r,s))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, r,q,p,s))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, r,s,p,q))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, p,s,r,q))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, q,p,s,r))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, s,p,q,r))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, s,r,q,p))
                    = value;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, q,r,s,p))
                    = value;

                double asvalue = pqsrElements(pqrs);
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, p,q,s,r))
                    = asvalue;
                twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates, q,p,r,s))
                    = asvalue;

                pqrs++;
            } // end for p,q,r,s
        } // end function setNonAntiSymmetrizedElements

        void readTwoBodyMatrix() {
            /* read in two-body non-antisymmetrized elements and set *
             * anti-symmetrized elements */
            std::ifstream twoBodyFile;
            twoBodyFile.open(twoBodyFileName);
            double buf = 0;
            if (twoBodyFile.is_open()) {
                for (unsigned int p = 0; p < totalSize; ++p)
                for (unsigned int q = 0; q < totalSize; ++q)
                for (unsigned int r = 0; r < totalSize; ++r)
                for (unsigned int s = 0; s < totalSize; ++s)
                {
                    if ((p < m_numStates) && 
                        (q < m_numStates) && 
                        (r < m_numStates) && 
                        (s < m_numStates)) {
                        twoBodyFile >> twoBodyNonAntiSymmetrizedElements(
                                dIndex(m_numStates, p,q,r,s));
                    } else {
                        twoBodyFile >> buf;
                    } // end ifelse
                } // end for p,q,r,s
                twoBodyFile.close();
            } // end if
        } //  end function readTwoBodyMatrix

        bool checkAndPrependFileName() {
            /* Prepend path to filename, check if file with larger matrix than
             * given cutoff exists, use that one if it exists. */
            namespace fs = std::experimental::filesystem;

            // prepend path
            twoBodyFileName = dirpath + "/" + twoBodyFileName;

            for (auto& x : fs::directory_iterator(dirpath)) {
                std::string Lstring = "";
                std::string filex = x.path().filename().replace_extension("");
                std::string filew = "";
                for (unsigned int i = filex.find_first_of("w")+1; i <
                        filex.find_first_of("_")+1; ++i) {
                    filew += filex.at(i);
                } // end fori
                if (std::stod(filew) != omega) {
                    continue;
                } // end if
                std::string fileD = "";
                for (unsigned int i = filex.find_first_of("D")+1; i <
                        filex.find_first_of("_")+1; ++i) {
                    fileD += filex.at(i);
                } // end fori
                if (std::stod(fileD) != m_dim) {
                    continue;
                } // end if
                for (unsigned int i = filex.find_first_of("L")+1; i <
                        filex.size(); ++i) {
                    Lstring += filex.at(i);
                } // end fori
                int fileL = 0;
                try {
                    fileL = std::stoi(Lstring);
                } catch (std::string) {
                    continue;
                } // end try-catch
                if (fileL > m_basisSize) {
                    /* return once suitable file is found */
                    totalSize = fileL/2;
                    twoBodyFileName = x.path();
                    return true;
                } // end if
            } // end forx

            return fs::exists(twoBodyFileName);
        } // end function checkFileName

        void weightSizesAndDispl(Eigen::ArrayXi& sizes, Eigen::ArrayXi& displ,
                const Eigen::ArrayXXi& pqMap) {
            /* Weight the sizes based on the sum of the product of
             * (n_p,n_q,n_r,n_s) for each process. make sure total sum of pqrs
             * within process chunk is within originalMean */
            Eigen::ArrayXi sums = Eigen::ArrayXi::Zero(numProcs);
            for (unsigned int i = 0; i < sums.size(); ++i) {
                const Eigen::Ref<const Eigen::ArrayXXi> pqBlock =
                    pqMap.block(displ(i), 0, sizes(i), 2);
                for (unsigned int k = 0; k < pqBlock.rows(); ++k) {
                    sums(i) +=
                        (m_I->GaussianBasis::getnStates(pqBlock(k,0)).array() +
                         1).prod() *
                        (m_I->GaussianBasis::getnStates(pqBlock(k,1)).array() +
                         1).prod();
                } // end fork
            } // end fori

            // find mean value for reference
            sizes.setZero();
            int originalMean = ceil(sums.mean());
            int proc = 0;
            int jsum = 0;
            int k = 0;
            int maxnFactor = originalMean + 3 *
                pow(m_I->GaussianBasis::getn().maxCoeff()+1, 2);
            for (unsigned int pq = 0; pq < pqMap.rows(); ++pq) {
                // iterate over rows in pqMap and add the sum product for each
                // row as above. Break when jsum is larger than original mean,
                // at which size for process proc is set, the sum is reset and
                // next process is taken in.
                auto resetAndGoToNext = [&]() {
                    /* reset values for next process */
                    k=0;
                    jsum = 0;
                    proc++;
                }; // end lambda resetAndGoToNext

                jsum += (m_I->GaussianBasis::getnStates(pqMap(pq,0)).array() +
                        1).prod() *
                    (m_I->GaussianBasis::getnStates(pqMap(pq,1)).array() +
                     1).prod();
                k++;

                if (jsum >= maxnFactor) {
                    /* dont count if overshot */
                    sizes(proc) = k-1;
                    pq -= 1;
                    resetAndGoToNext();
                } else if (jsum > originalMean) {
                    /* go to next process */
                    sizes(proc) = k;
                    resetAndGoToNext();
                } else if ((pq == pqMap.rows()-1) && (jsum <= originalMean)) {
                    sizes(proc) = k;
                } // end ifeif
            } // end for pq

            int n2 = m_numStates * (m_numStates + 1);
            n2 /= 2;
            sizes(numProcs-1) = n2 - sizes.head(numProcs-1).sum();

            // recalculate relative displacements
            for (unsigned int p = 0; p < numProcs; ++p) {
                displ(p) = sizes.head(p).sum();
            } // end forp
        } // end function weightSizesAndDispl

    protected:
        unsigned int m_dim, m_numStates, m_numParticles, m_basisSize;

        double omega;

        std::string twoBodyFileName;

        inline void assemble(unsigned int progressDivider=0) {
            /* assemble integral elements (with symmetries) */

            // class Cartesian (in derived basis) creates full-shell with both
            // spins(each one-body function is created for both spin up and
            // down), divide by two to take the restricted form of the
            // Hartree-Fock equation.

            m_numStates = m_basisSize / 2;
            totalSize = m_numStates;

            // matrix containing elements <i|h|j> (one-body elements) and
            // elements <i|j> (overlap elements)
            if (myRank == 0) {
                /* only calculate oneBody and overlap elements in root */
                oneBodyElements = Eigen::MatrixXd::Zero(m_numStates,
                        m_numStates);
                overlapElements = Eigen::MatrixXd::Zero(m_numStates,
                        m_numStates);

                // set one-body (uncoupled) elements and overlap elements
                for (unsigned int p = 0; p < m_numStates; ++p) {
                    for (unsigned int q = p; q < m_numStates; ++q) {
                        overlapElements(p,q) = m_I->overlapElement(p,q);
                        oneBodyElements(p,q) = m_I->oneBodyElement(p,q);
                        if (p != q) {
                            /* only off-diagonal symmetric elements need to be
                             * set */
                            overlapElements(q,p) = overlapElements(p,q);
                            oneBodyElements(q,p) = oneBodyElements(p,q);
                        } // end if
                    } // end forq
                } // end forp
            } // end if

            // set two-body coupled (Coulomb) integral elements
            if (interaction) {
                // let root check if file of two-body elements exists
                twoBodyNonAntiSymmetrizedElements =
                    Eigen::ArrayXd::Zero(m_numStates * m_numStates *
                            m_numStates * m_numStates);
                bool fileExists = false;
                if (myRank == 0) {
                    if (twoBodyFileName.compare("")) {
                        fileExists = checkAndPrependFileName();
                        if (fileExists) {
                            readTwoBodyMatrix();
                        } // end if
                    } // end if
                } // end if

                // let slaves know of the file-check clarity (#religous?)
                MPI_Bcast(&fileExists, 1, MPI_INT, 0, MPI_COMM_WORLD);

                if (fileExists == 0) {
                    /* only calculate if file does not exist */
                    // create matrix containing pairs (p,q)
                    int subSize = m_numStates * (m_numStates+1);
                    subSize /= 2;
                    Eigen::ArrayXd pqrsElements, pqsrElements;
                    Eigen::ArrayXi displ(numProcs);
                    Eigen::ArrayXi sizes(numProcs);

                    auto calculateDispl = [&,this]() {
                        for (int p = 0; p < numProcs; ++p) {
                            displ(p) = sizes.head(p).sum();
                        } // end forp
                    };
                    using EigenRowArrayXXi = Eigen::Array<int, Eigen::Dynamic,
                          Eigen::Dynamic, Eigen::RowMajor>;
                    EigenRowArrayXXi pqMap;

                    if (myRank == 0) {
                        pqMap = EigenRowArrayXXi::Zero(subSize, 2);
                        int j = 0;
                        for (unsigned int p = 0; p < m_numStates; ++p)
                        for (unsigned int q = p; q < m_numStates; ++q)
                        {
                            pqMap(j,0) = p;
                            pqMap(j,1) = q;
                            j++;
                        } // end for p,q,r,s
                        
                        // distribute sizes
                        for (int p = 0; p < numProcs; ++p) {
                            sizes(p) = Methods::divider(p, subSize, numProcs);
                            displ(p) = sizes.head(p).sum();
                        } // end forp

                        weightSizesAndDispl(sizes, displ, pqMap);
                    } // end if

                    int mySize;
                    MPI_Scatter(sizes.data(), 1, MPI_INT, &mySize, 1, MPI_INT,
                            0, MPI_COMM_WORLD);

                    EigenRowArrayXXi mypqMap = EigenRowArrayXXi(mySize, 2);

                    if (myRank == 0) {
                        sizes *= 2;
                        calculateDispl();
                        std::cout << sizes << std::endl;
                    } // end if
                    MPI_Scatterv(pqMap.data(), sizes.data(), displ.data(),
                            MPI_INT, mypqMap.data(), mySize*2, MPI_INT, 0,
                            MPI_COMM_WORLD);

                    // array containing two-body elements <ij|1/r_12|kl> for
                    // subset (r,s) of set of range of (p,q) in each process
                    int pSize = mySize * subSize;
                    Eigen::ArrayXd myTmpTwoBody(pSize), myTmpTwoBodyAS(pSize);

                    // save first part of progress bar
                    std::string progressPosition, progressBuffer;
                    if (progressDivider) {
                        progressDivider = (int)exp(fmod(4.5, pSize));
                        progressPosition = Methods::stringPos(myRank, 3) +
                            "Progress: [";
                    } // end if

                    int rs = 0;
                    for (int pq = 0; pq < mySize; ++pq) {
                        const int& p = mypqMap(pq,0);
                        const int& q = mypqMap(pq,1);
                        for (unsigned int r = 0; r < m_numStates; ++r) {
                            for (unsigned int s = r; s < m_numStates; ++s) {
                                myTmpTwoBody(rs) = m_I->coulombElement(p, q, r,
                                        s);
                                myTmpTwoBodyAS(rs) = m_I->coulombElement(p, q,
                                        s, r);
                                
                                // print progress
                                if (progressDivider) {
                                    /* show progress if given */
                                    if (!(static_cast<int>(fmod(rs,
                                                        Methods::divider( rs,
                                                            pSize,
                                                            progressDivider)))))
                                    {
                                        /* print only a few times */
                                        progressBuffer = progressPosition;
                                        Methods::printProgressBar(progressBuffer,
                                                (float)((rs==pSize-1) ?  rs :
                                                    (rs+1)) / pSize, 23,
                                                "Two-Body");
                                    } // end if
                                } // end if

                                rs++;
                            } // end fors
                        } // end forr
                    } // end forpq

                    // gather subresults from slaves into complete matrix in
                    // root
                    if (myRank == 0) {
                        pqrsElements = Eigen::ArrayXd::Zero(subSize*subSize);
                        pqsrElements = Eigen::ArrayXd::Zero(subSize*subSize);
                        sizes /= 2;
                        sizes *= subSize;
                        calculateDispl();
                    } // end if
                    MPI_Gatherv(myTmpTwoBody.data(), pSize, MPI_DOUBLE,
                            pqrsElements.data(), sizes.data(), displ.data(),
                            MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    MPI_Gatherv(myTmpTwoBodyAS.data(), pSize, MPI_DOUBLE,
                            pqsrElements.data(), sizes.data(), displ.data(),
                            MPI_DOUBLE, 0, MPI_COMM_WORLD);
                
                    if (myRank == 0) {
                        setNonAntiSymmetrizedElements(pqrsElements,
                                pqsrElements);
                        writeTwoBodyElementsToFile();
                    } // end if
                } // end if

                MPI_Bcast(twoBodyNonAntiSymmetrizedElements.data(),
                        twoBodyNonAntiSymmetrizedElements.size(), MPI_DOUBLE,
                        0, MPI_COMM_WORLD);
            } // end if
        } // end function assemble

    public:
        HartreeFockSolver(Integrals* IIn, const unsigned int dimension,
                unsigned int cut, const unsigned int numParticles) {
            /* set dimensions, cutoff and number of particles and initialize
             * basis and integrals */
            dirpath = "src/integrals/inputs";

            twoBodyFileName = "";

            m_I = IIn;

            m_basisSize = cut;

            m_dim = dimension;
            m_numParticles = numParticles;
            interaction = true;

            // grab info from default communicator 
            MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
            MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
        } // end constructor

        virtual ~HartreeFockSolver () {
        } // end deconstructor

        double iterate(const unsigned int& maxIterations, const double& eps,
                const unsigned int progressDivider) {
            /* run Hartree-Fock algorithm for finding coefficients and energy
             * until threshold convergence or until maxIterations is reached */

            // pre-calculate one- and two-body matrix-elements and set initial
            // density matrix with coefficient matrix set to identity
            assemble(progressDivider);

            if (myRank == 0) {
                coefficients = Eigen::MatrixXd::Identity(m_numStates,
                        m_numStates);
                densityMatrix = Eigen::MatrixXd::Zero(m_numStates,
                        m_numStates);
                setDensityMatrix();
                FockMatrix = Eigen::MatrixXd::Zero(m_numStates, m_numStates);
                Eigen::VectorXd previousEnergies =
                    Eigen::VectorXd::Zero(m_numStates);
                Eigen::MatrixXd oldCoefficients = coefficients;
                Eigen::VectorXd energies = Eigen::VectorXd::Zero(m_numStates);

                // initialize eigenvalue/vector solver for hermitian matrix
                // (Fock matrix is build to be hermitian)
                Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
                    eigenSolver;

                // run Hartree-Fock algorithm
                std::string progressPosition, progressBuffer;
                progressPosition = Methods::stringPos(myRank, 3) + "Progress:"
                    " [";
                for (unsigned int count = 0; count < maxIterations; ++count) {
                    /* run for maxIterations or until convergence is reached */

                    // set HF-matrix with current coefficients
                    setFockMatrix();

                    // find eigenvalues and eigenvector (HartreeFock-energies
                    // and coefficients respectively)
                    eigenSolver.compute(FockMatrix, overlapElements);

                    // find eigenvalues and eigenvectors
                    energies = eigenSolver.eigenvalues();
                    coefficients = eigenSolver.eigenvectors();

                    // set density matrix with new coefficients
                    setDensityMatrix();
                    
                    // perform mixing (for convergence, dont ask why...)
                    densityMatrix = mixingFactor*densityMatrix + (1-mixingFactor)
                        * oldCoefficients;

                    // check for convergence with RMS of difference between
                    // previous and current energies 
                    double diff = sqrt((energies - previousEnergies).norm() /
                            m_numStates);
                    energy = groundStateEnergy(energies);

                    if (diff < eps) {
                        iterations = count;
                        break;
                    } // end if

                    // keep old values
                    oldCoefficients = densityMatrix;
                    previousEnergies = energies;

                    // print progress
                    if (progressDivider) {
                        /* show progress if given */
                        if (!(static_cast<int>(fmod(count,
                                            Methods::divider(count,
                                                maxIterations,
                                                progressDivider))))) {
                            /* print only a few times */
                            progressBuffer = progressPosition;
                            Methods::printProgressBar(progressBuffer,
                                    (float)((count==maxIterations-1) ? count :
                                        (count+1)) / maxIterations, 23, "HF");
                        } // end if
                    } // end if
                } // end forcount

                return energy;
            } // end if

            return 0;
        } // end function iterate

        double groundStateEnergy(const Eigen::VectorXd& eigVals) {
            // find estimate for ground state energy for m_numParticles
            double E = 2*eigVals.segment(0, m_numParticles/2).sum();
            for (unsigned int a = 0; a < m_numStates; ++a)
            for (unsigned int b = 0; b < m_numStates; ++b)
            for (unsigned int c = 0; c < m_numStates; ++c)
            for (unsigned int d = 0; d < m_numStates; ++d)
            {
                E -= densityMatrix(a,b) * densityMatrix(c,d) *
                    (2*twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates,
                                                                a,c,b,d)) -
                     twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates,
                             a,c,d,b)));
            } // end for a,b,c,d

            return E;
        }  // end function groundStateEnergy

         double getTwoBodyElement(const unsigned int& i, const unsigned int& j,
                 const unsigned int& k, const unsigned int& l) {
            /* return two-body non-antisymmetrized elements */
            return twoBodyNonAntiSymmetrizedElements(dIndex(m_numStates,
                        i,j,k,l));
        } // end function getTwoBodyElement 

        const unsigned int& getIterations() const {
            return iterations;
        } // end function getIterations

        void setInteraction(bool a) {
            /* set interaction on (if a=true) or false (if a=false) */
            interaction = a;
        } // end function setInteraction

        void writeTwoBodyElementsToFile() {
            /* write elements to file */
            std::ofstream twoBodyFile(twoBodyFileName);
            if (twoBodyFile.is_open()) {
                for (unsigned int p = 0; p < m_numStates; ++p)
                for (unsigned int q = 0; q < m_numStates; ++q)
                for (unsigned int r = 0; r < m_numStates; ++r)
                for (unsigned int s = 0; s < m_numStates; ++s)
                {
                    twoBodyFile << std::fixed << std::setprecision(16) <<
                        twoBodyNonAntiSymmetrizedElements( dIndex(m_numStates,
                                    p,q,r,s));
                    twoBodyFile << " ";
                } // end for p,q,r,s
            } // end if
            twoBodyFile.close();
        } // end function writeTwoBodyElementsToFile

        void writeCoefficientsToFile(const std::string& filename, const
                std::string& omega) {
            /* write coefficients (in order basis functions were given in
             * integral object) to YAML file */
            std::ofstream outFile(filename + ".yaml");
            YAML::Node info;
            info["E0"] = energy;
            info["I"] = iterations;
            info["omega"] = omega;
            info["dim"] = m_dim;
            info["numbasis"] = m_numStates;
            info["numparticles"] = m_numParticles;
            std::vector<double> tmpCol(m_numStates);
            for (unsigned int p = 0; p < m_numParticles/2; ++p) {
                for (unsigned int i = 0; i < m_numStates; ++i) {
                    tmpCol[i] = coefficients(i,p);
                } // end fori
                info["coeffs"].push_back(tmpCol);
            } // end forp
            if (outFile.is_open()) {
                outFile << info;
            } else {
                std::cout << "File : " + filename + "could not be opened." <<
                    std::endl;
            } // end ifelse
        } // end function writeCoefficientsToFile
};

#endif /* HARTREEFOCKSOLVER_H */
