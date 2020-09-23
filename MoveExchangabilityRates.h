#ifndef MoveExchangabilityRates_H
#define MoveExchangabilityRates_H

#include <vector>
#include "Move.h"
class ParameterExchangabilityRates;



class MoveExchangabilityRates : public Move {

    public:
                                        MoveExchangabilityRates(RandomVariable* r, Model* m, std::string nm, ParameterExchangabilityRates* f[2], bool isTimeReversible);
        void                            accept(void);
        void                            reject(void);
        void                            restore(void);
        void                            tune(void);
        double                          update(void);
        std::vector<double>             values(void);

    protected:
        ParameterExchangabilityRates*   myParameter[2];
        int                             numRates;
        double                          tuning;
        std::vector<double>             oldValue;
};

#endif
