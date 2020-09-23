#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include "EigenSystem.h"
#include "Msg.h"



EigenSystem::EigenSystem(void) {

    activeVals = 0;

    ccIjk[0] = new std::complex<double>[64];
    ccIjk[1] = new std::complex<double>[64];
    pi[0].resize(4);
    pi[1].resize(4);
}

EigenSystem::~EigenSystem(void) {

    delete [] ccIjk[0];
    delete [] ccIjk[1];
}

void EigenSystem::calculateEigenSystem(NucleotideSquareMatrix_t& Q) {

    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    Eigen::EigenSolver< NucleotideSquareMatrix_t > eigSolver;
    eigSolver.compute( Q, true );
    eigenValues[activeVals] = eigSolver.eigenvalues();
    Eigen::Matrix<std::complex<double>, 4, 4> eigenVectors = eigSolver.eigenvectors();
    Eigen::Matrix<std::complex<double>, 4, 4> inverseEigenVectors = eigenVectors.inverse();

    // calculate cc_ijk
    std::complex<double>* pc = ccIjk[activeVals];
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                 *(pc++) = eigenVectors(i,k) * inverseEigenVectors(k,j);
}

std::vector<double> EigenSystem::calulateStationaryFrequencies(NucleotideSquareMatrix_t& Q) {
    
    Eigen::VectorXd f = Q.transpose().fullPivLu().kernel();
    f = f / f.sum();
    std::vector<double> stationaryFrequencies(4);
    for (int i=0; i<4; i++)
        stationaryFrequencies[i] = f(i);
    
    for (int i=0; i<4; i++)
        {
        if (f[i] < 0.0)
            {
            std::cout << f << std::endl;
            Msg::error("Negative stationary frequency");
            }
        }
    
    return stationaryFrequencies;
}

void EigenSystem::flipActiveValues(void) {

    if (activeVals == 0)
        activeVals = 1;
    else
        activeVals = 0;
}

std::complex<double>* EigenSystem::getCijk(void) {

    return ccIjk[activeVals];
}

Eigen::Matrix<std::complex<double>, 4, 1>& EigenSystem::getEigenValues(void) {

    return eigenValues[activeVals];
}

double EigenSystem::sumVector(std::vector<double>& v) {

    double sum = 0.0;
    for (int i=0; i<4; i++)
        sum += v[i];
    return sum;
}

void EigenSystem::testStationaryFrequencies(NucleotideSquareMatrix_t& Q) {

}
