#include "Move.h"



Move::Move(RandomVariable* r, Model* m, std::string nm) {

    rv = r;
    model = m;
    moveName = nm;
    numTries = 0;
    numAccepted = 0;
}

void Move::clear(void) {

    numTries = 0;
    numAccepted = 0;
}
