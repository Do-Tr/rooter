#include <cmath>
#include <iomanip>
#include <iostream>
#include "ParameterAsrv.h"
#include "RandomVariable.h"




ParameterAsrv::ParameterAsrv(RandomVariable* rp, Model* mp, std::string nm, int n, double alp) : Parameter(rp, mp, nm) {

    numCategories = n;
    alphaPrior = alp;
    alpha = 1.0 / alphaPrior; // set to the expected value
    rates.resize(numCategories);
    rv->discretizeGamma(rates, alpha, alpha, numCategories, false);
}

ParameterAsrv::ParameterAsrv(ParameterAsrv& b) : Parameter(b.rv, b.modelPtr, b.name) {

    clone(b);
}

Parameter& ParameterAsrv::operator=(Parameter& b) {

    if (this != &b)
        {
        ParameterAsrv* dc = dynamic_cast<ParameterAsrv*>(&b);
        clone(*dc);
        }
    return *this;
}

ParameterAsrv& ParameterAsrv::operator=(ParameterAsrv& b) {

    if (this != &b)
        {
        clone(b);
        }
    return *this;
}

double& ParameterAsrv::operator[](int idx) {

    if (idx >= rates.size())
        {
        std::cout << "Gamma rates index out of bounds";
        exit(0);
        }
    return rates[idx];
}

void ParameterAsrv::clone(ParameterAsrv& b) {

    numCategories = b.numCategories;
    alphaPrior = b.alphaPrior;
    alpha = b.alpha;
    rates = b.rates;
}

std::string ParameterAsrv::getParmHeader(int n) {

    std::string str = "Alphas";
    if  ( n != -1 )
        str += "[" + std::to_string(n) + "]";
    str += '\t';
    return str;
}

std::string ParameterAsrv::getParmString(int n) {

    std::string str = std::to_string(alpha) + '\t';
    return str;
}

double ParameterAsrv::lnPriorProb(void) {

    return rv->lnExponentialPdf( alphaPrior, alpha );
}

void ParameterAsrv::print(void) {

    std::cout << "Shape parameter = " << std::fixed << std::setprecision(3) << alpha;
    std::cout << " ( ";
    double sum = 0.0;
    for (int i=0; i<rates.size(); i++)
        {
        sum += rates[i];
        std::cout << std::fixed << std::setprecision(5) << rates[i] << " ";
        }
    std::cout << ") Average = " << sum / rates.size() << std::endl;
}

void ParameterAsrv::setAsrv(double x) {

    alpha = x;
    rv->discretizeGamma(rates, alpha, alpha, numCategories, false);
}
