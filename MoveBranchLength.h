#ifndef MoveBranchLength_H
#define MoveBranchLength_H

#include "Move.h"
class Branch;
class ParameterTree;
class RandomVariable;



class MoveBranchLength : public Move {

    public:
                                MoveBranchLength(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]);
        void                    accept(void);
        void                    reject(void);
        void                    restore(void);
        void                    tune(void);
        double                  update(void);
        std::vector<double>     values(void);
    
    protected:
        ParameterTree*          myParameter[2];
        double                  tuning;
        Branch*                 oldBranch1;
        double                  oldValue;
};

#endif
