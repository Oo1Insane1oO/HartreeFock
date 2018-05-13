#ifndef CARTESIAN_H
#define CARTESIAN_H

#include <string>

#include <eigen3/Eigen/Dense>

class Cartesian {
    using EigenIntPtrMat = Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic>;
    using EigenIntPtrVec = Eigen::Matrix<int*, Eigen::Dynamic, 1>;
    private:
        int s;
        unsigned int m_dim, m_size, m_numStates;
        Eigen::VectorXi n, ms, E, M, m;

        Eigen::Matrix<int*, Eigen::Dynamic, Eigen::Dynamic> states;

        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int, const unsigned int);
        void addState(unsigned int&, const unsigned int, const unsigned int,
                const unsigned int);
        void findPrincipal(const unsigned int& , Eigen::MatrixXi&);

        unsigned int calculateDegeneracy(const unsigned int&);

        void setStates(const unsigned int&);
    
    public:
        Cartesian ();
        virtual ~Cartesian ();

        void setup(unsigned int, const unsigned int);

        const EigenIntPtrMat &getStates() const;
        const Eigen::Ref<const EigenIntPtrVec> getStates(const unsigned int&)
            const;
        const int& getn(const unsigned int&, const unsigned int&) const;
        const Eigen::VectorXi &getn() const;
        const int &getn(const unsigned int&) const;
        const Eigen::VectorXi &getE() const;
        const int &getE(const unsigned int&) const;
        const Eigen::VectorXi &getMagic() const;
        const int &getMagic(const unsigned int&) const;
        const unsigned int& getSize() const;
        const unsigned int& getNumberOfStates() const;

        void restructureStates();

        void printStates();

        void writeToFile(std::string);
};

#endif /* CARTESIAN_H */
