#ifndef MoveRoot_H
#define MoveRoot_H

#include "Move.h"
class Node;
class ParameterTree;



class MoveRoot : public Move {

    public:
                            MoveRoot(RandomVariable* r, Model* m, std::string nm, ParameterTree* t[2]);
        void                accept(void);
        void                reject(void);
        void                restore(void);
        void                tune(void);
        double              update(void);
        std::vector<double> values(void);

    protected:
        ParameterTree*      myParameter[2];
        double              tuning;
};

#endif
