#ifndef ParameterExchangabilityRates_H
#define ParameterExchangabilityRates_H

#include <vector>
#include "Parameter.h"



class ParameterExchangabilityRates : public Parameter {

    public:
                                        ParameterExchangabilityRates(RandomVariable* rp, Model* mp, std::string nm, bool isTimeReversible);
                                        ParameterExchangabilityRates(ParameterExchangabilityRates& b);
        Parameter&                      operator=(Parameter& b);
        ParameterExchangabilityRates&   operator=(ParameterExchangabilityRates& b);
        double&                         operator[](int idx);
        std::string                     getParmString(int n);
        std::string                     getParmHeader(int n);
        std::vector<double>&            getExchangabilityRates(void) { return f; }
        std::vector<double>&            getExchangabilityAlpha(void) { return a; }
        double                          lnPriorProb(void);
        void                            print(void);
        void                            setExchangabilityRates(std::vector<double> x) { f = x; }
        int                             size(void) { return (int)f.size(); }

    protected:
        void                            clone(ParameterExchangabilityRates& b);
        bool                            sumToWithin(double x, double tol, std::vector<double>& vec);
        std::vector<double>             a;
        std::vector<double>             f;
};

#endif
