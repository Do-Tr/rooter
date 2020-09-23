#ifndef MoveTreeLength_H
#define MoveTreeLength_H

#include "Move.h"
class ParameterTree;
class RandomVariable;



class MoveTreeLength : public Move {

    public:
                                MoveTreeLength(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]);
        void                    accept(void);
        void                    reject(void);
        void                    restore(void);
        void                    tune(void);
        double                  update(void);
        std::vector<double>     values(void);

    protected:
        ParameterTree*          myParameter[2];
        double                  tuning;
        double                  oldValue;
};

#endif
