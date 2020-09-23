#include <iostream>
#include "Alignment.h"
#include "Mcmc.h"
#include "Model.h"
#include "RandomVariable.h"
#include "Settings.h"

void printHeader(void);



int main(int argc, char* argv[]) {

    // print greeting
    printHeader();

    // instantiate the random number generator
    RandomVariable rv;
    
    // read the user settings
    Settings settings(argc, argv);

    // read the alignment
    Alignment alignment(settings.getInputFileName());
    
    // set up the model
    Model model(&rv, &settings, &alignment);

    // run chain
    Mcmc mcmc(&rv, &model, &settings);
    mcmc.run();

    return 0;
}

void printHeader(void) {

    std::cout << std::endl;
    std::cout << "   M. C. Rooter, version 1.0" << std::endl;
    std::cout << std::endl;
    std::cout << "   * John Huelsenbeck, Rasmus Nielsen, and Hongru Wang" << std::endl;
    std::cout << "   * University of California, Berkeley" << std::endl;
    std::cout << "   * " << std::endl;
}
