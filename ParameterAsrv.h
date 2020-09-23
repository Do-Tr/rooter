#ifndef ParameterAsrv_H
#define ParameterAsrv_H

#include <string>
#include <vector>
#include "Parameter.h"



class ParameterAsrv : public Parameter {

    public:
                                ParameterAsrv(RandomVariable* rp, Model* mp, std::string nm, int n, double alp);
                                ParameterAsrv(ParameterAsrv& b);
        Parameter&              operator=(Parameter& b);
        ParameterAsrv&          operator=(ParameterAsrv& b);
        double&                 operator[](int idx);
        std::string             getParmHeader(int n);
        std::string             getParmString(int n);
        double&                 getAsrv(void) { return alpha; }
        std::vector<double>&    getRates(void) { return rates; }
        double                  lnPriorProb(void);
        void                    print(void);
        void                    setAsrv(double x);

    protected:
        void                    clone(ParameterAsrv& b);
        double                  alpha;
        double                  alphaPrior;
        int                     numCategories;
        std::vector<double>     rates;
};

#endif
