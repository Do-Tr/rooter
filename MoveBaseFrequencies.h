#ifndef MoveBaseFrequencies_H
#define MoveBaseFrequencies_H

#include "Move.h"
class ParameterBaseFrequencies;
class RandomVariable;



class MoveBaseFrequencies : public Move {

    public:
                                    MoveBaseFrequencies(RandomVariable* r, Model* m, std::string nm, ParameterBaseFrequencies* f[2]);
        void                        accept(void);
        void                        reject(void);
        void                        restore(void);
        void                        tune(void);
        double                      update(void);
        std::vector<double>         values(void);

    protected:
        ParameterBaseFrequencies*   myParameter[2];
        double                      tuning;
        std::vector<double>         oldValue;
};

#endif
