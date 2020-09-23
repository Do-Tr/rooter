#ifndef MoveAsrv_H
#define MoveAsrv_H

#include "Move.h"
class ParameterAsrv;
class RandomVariable;



class MoveAsrv : public Move {

    public:
                                MoveAsrv(RandomVariable* r, Model* m, std::string nm, ParameterAsrv* gs[2]);
        void                    accept(void);
        void                    reject(void);
        void                    restore(void);
        void                    tune(void);
        double                  update(void);
        std::vector<double>     values(void);
    
    protected:
        ParameterAsrv*          myParameter[2];
        double                  tuning;
        double                  oldValue;
};

#endif
