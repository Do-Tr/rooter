#include <cmath>
#include <iomanip>
#include <iostream>
#include "Msg.h"
#include "ParameterExchangabilityRates.h"
#include "RandomVariable.h"



ParameterExchangabilityRates::ParameterExchangabilityRates(RandomVariable* rp, Model* mp, std::string nm, bool isTimeReversible) : Parameter(rp, mp, nm) {

    int n = 6;
    if (isTimeReversible == false)
        n = 12;
    f.resize( n );
    a.resize( n );
    for (int i=0; i<a.size(); i++)
        a[i] = 1.0;
    rv->dirichletRv(a, f);
    for (int i=0; i<f.size(); i++)
        f[i] = a[i] / (double)n; // expected value, for now
}

ParameterExchangabilityRates::ParameterExchangabilityRates(ParameterExchangabilityRates& b) : Parameter(b.rv, b.modelPtr, b.name) {

    f.resize( b.f.size() );
    a.resize( b.a.size() );
    clone(b);
}

Parameter& ParameterExchangabilityRates::operator=(Parameter& b) {

    if (this != &b)
        {
        ParameterExchangabilityRates* dc = dynamic_cast<ParameterExchangabilityRates*>(&b);
        clone(*dc);
        }
    return *this;
}

ParameterExchangabilityRates& ParameterExchangabilityRates::operator=(ParameterExchangabilityRates& b) {

    if (this != &b)
        {
        clone(b);
        }
    return *this;
}

double& ParameterExchangabilityRates::operator[](int idx) {

    if (idx >= f.size())
        {
        std::cout << "Exchangability index out of bounds";
        exit(0);
        }
    return f[idx];
}

void ParameterExchangabilityRates::clone(ParameterExchangabilityRates& b) {

    for (int i=0; i<b.a.size(); i++)
        {
        a[i] = 1.0;
        f[i] = b.f[i];
        }
}

void ParameterExchangabilityRates::print(void) {

    std::cout << "Exchangability Rates = ";
    for (int i=0; i<f.size(); i++)
        std::cout << std::fixed << std::setprecision(4) << f[i] << " ";
    std::cout << '\n';
}

double ParameterExchangabilityRates::lnPriorProb(void) {

    return rv->lnDirichletPdf(a, f);
}

bool ParameterExchangabilityRates::sumToWithin(double x, double tol, std::vector<double>& vec) {

    double sum = 0.0;
    for (int i=0; i<vec.size(); i++)
        sum += vec[i];
    
    if (fabs(sum-x) < tol)
        return true;
    return false;
}

std::string ParameterExchangabilityRates::getParmString(int n) {

    std::string tempStr = "";
    for (int i=0; i<f.size(); i++)
        {
        char temp[20];
        sprintf(temp, "%1.4lf\t", f[i]);
        tempStr += temp;
        }
    return tempStr;
}

std::string ParameterExchangabilityRates::getParmHeader(int n) {

    std::string tempString = "";
    char nuc[4] = { 'A', 'C', 'G', 'T' };
    if (f.size() == 6)
        {
        for (int i=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                char temp[20];
                sprintf(temp, "R[%c%c]\t", nuc[i], nuc[j]);
                tempString += temp;
                }
            }
        }
    else
        {
        for (int i=0; i<4; i++)
            {
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    {
                    char temp[20];
                    sprintf(temp, "R[%c%c]\t", nuc[i], nuc[j]);
                    tempString += temp;
                    }
                }
            }
        }
    return tempString;
}
