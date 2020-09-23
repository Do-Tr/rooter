#ifndef MoveNodePosition_H
#define MoveNodePosition_H

#include "Move.h"
class Branch;
class ParameterTree;
class RandomVariable;



class MoveNodePosition : public Move {

    public:
                                MoveNodePosition(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]);
        void                    accept(void);
        void                    reject(void);
        void                    restore(void);
        void                    tune(void);
        double                  update(void);
        std::vector<double>     values(void);
    
    protected:
        ParameterTree*          myParameter[2];
        double                  tuning;
        Branch*                 oldBranch[10];
        double                  oldValue[10];
};

#endif
