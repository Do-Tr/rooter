#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include "EigenSystem.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "ParameterAsrv.h"
#include "TransitionProbabilities.h"




TransitionProbabilities::TransitionProbabilities(int k, Model* m) {

    numStates = 4;
    numGammaCats = k;
    myModel = m;
    
    // dynamically allocate a vector of matrices
    p = new double**[numGammaCats];
    for (int c=0; c<numGammaCats; c++)
        {
        p[c] = new double*[numStates];
        p[c][0] = new double[numStates * numStates];
        for (int i=1; i<numStates; i++)
            p[c][i] = p[c][i-1] + numStates;
        for (int i=0; i<numStates; i++)
            for (int j=0; j<numStates; j++)
                p[c][i][j] = 0.0;
        }
}

TransitionProbabilities::~TransitionProbabilities(void) {

    delete [] p[0][0];
    delete [] p[0];
    delete [] p;
}

void TransitionProbabilities::print(void) {

    std::vector<double>& rates = myModel->getAsrv()->getRates();
    for (int c=0; c<numGammaCats; c++)
        {
        std::cout << "Gamma Category " << c << ": " << rates[c] << std::endl;
        for (int i=0; i<numStates; i++)
            {
            for (int j=0; j<numStates; j++)
                {
                std::cout << std::fixed << std::setprecision(6);
                std::cout << p[c][i][j] << " ";
                }
            std::cout << std::endl;
            }
        std::cout << std::endl;
        }
}

void TransitionProbabilities::tiProbs(double v) {

    EigenSystem& eigs = EigenSystem::eigenSystem();
    Eigen::Matrix<std::complex<double>, 4, 1>& ceigenvalue = eigs.getEigenValues();
    std::complex<double>* ccIjk = eigs.getCijk();
    std::vector<double>& rates = myModel->getAsrv()->getRates();
    
    for (int c=0; c<numGammaCats; c++)
        {
        double r = rates[c];

        std::complex<double> ceigValExp[4];
        for (int s=0; s<4; s++)
            ceigValExp[s] = exp(ceigenvalue[s] * v * r);

        std::complex<double>* ptr = ccIjk;
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                std::complex<double> sum = std::complex<double>(0.0, 0.0);
                for(int s=0; s<4; s++)
                    sum += (*ptr++) * ceigValExp[s];
                p[c][i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
                }
            }
        }
}

