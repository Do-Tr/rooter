#ifndef ParameterBaseFrequencies_H
#define ParameterBaseFrequencies_H

#include <vector>
#include "Parameter.h"



class ParameterBaseFrequencies : public Parameter {

    public:
                                        ParameterBaseFrequencies(RandomVariable* rp, Model* mp, std::string nm);
                                        ParameterBaseFrequencies(ParameterBaseFrequencies& b);
        Parameter&                      operator=(Parameter& b);
        ParameterBaseFrequencies&       operator=(ParameterBaseFrequencies& b);
        double&                         operator[](int idx);
        std::string                     getParmString(int n);
        std::string                     getParmHeader(int n);
        std::vector<double>&            getBaseFrequencies(void) { return f; }
        std::vector<double>&            getBaseFrequencyAlpha(void) { return a; }
        double                          lnPriorProb(void);
        void                            print(void);
        void                            setBaseFrequencies(std::vector<double> x) { f = x; };

    protected:
        void                            clone(ParameterBaseFrequencies& b);
        std::vector<double>             a;
        std::vector<double>             f;
};

#endif
