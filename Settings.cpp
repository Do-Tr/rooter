#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include "Settings.h"



Settings::Settings(int argc, char* argv[]) {

    std::vector<std::string> commandString;

#   if 0 
    commandString.push_back( "-i" );
    commandString.push_back( "/Users/johnh/Desktop/rasmus/humans.phymlAln" );
    commandString.push_back( "-o" );
    commandString.push_back( "/Users/johnh/Desktop/rasmus/test" );
    commandString.push_back( "-t" );
    commandString.push_back( "/Users/johnh/Desktop/rasmus/raxmlRootedTree.nw" );
    commandString.push_back( "-c" );
    commandString.push_back( "/Users/johnh/Desktop/rasmus/tipDates.txt" );
    commandString.push_back( "-tr" );
    commandString.push_back( "yes" );
    commandString.push_back( "-l" );
    commandString.push_back( "1000000" );
    commandString.push_back( "-p" );
    commandString.push_back( "100" );
    commandString.push_back( "-s" );
    commandString.push_back( "100" );
    commandString.push_back( "-lambda" );
    commandString.push_back( "1000.0" );
    commandString.push_back( "-g" );
    commandString.push_back( "4" );
    commandString.push_back( "-e" );
    commandString.push_back( "2.0" );

#   else

    for (int i=1; i<argc; i++)
        {
        std::string cmd = argv[i];
        commandString.push_back(cmd);
        }
#   endif

    enum Mode { DATA_FILE, TREE_FILE, OUTPUT_FILE, TIPTIMES_FILE, CHAIN_LENGTH, PRINT_FREQ, SAMPLE_FREQ, BRLEN_LAMBDA, NUM_GAMMA_CATS, ASRV_LAMBDA, IS_REVERSIBLE, NONE };

    /* set default values for parameters */
    inputFileName          = "";
    treeFileName           = "";
    outPutFileName         = "";
    tipTimesFileName       = "";
    chainLength            = 1000000;
    printFrequency         = 1;
    sampleFrequency        = 100;
    brlenLambda            = 10.0;
    numGammaCats           = 4;
    asrvLambda             = 2.0;
    isReversible           = true;
    
    if (commandString.size() > 1)
        {
        if (commandString.size() % 2 != 0)
            {
            printCommandString(commandString);
            printUsage();
            }
            
        /* read the command-line arguments */
        int status = NONE;
        for (int i=0; i<commandString.size(); i++)
            {
            std::string cmd = commandString[i];
            //std::cout << cmd << std::endl;
            if (status == NONE)
                {
                /* read the parameter specifier */
                if ( cmd == "-i" )
                    status = DATA_FILE;
                else if ( cmd == "-t" )
                    status = TREE_FILE;
                else if ( cmd == "-o" )
                    status = OUTPUT_FILE;
                else if ( cmd == "-c" )
                    status = TIPTIMES_FILE;
                else if ( cmd == "-l" )
                    status = CHAIN_LENGTH;
                else if ( cmd == "-p" )
                    status = PRINT_FREQ;
                else if ( cmd == "-s" )
                    status = SAMPLE_FREQ;
                else if ( cmd == "-tr" )
                    status = IS_REVERSIBLE;
                else if ( cmd == "-lambda" )
                    status = BRLEN_LAMBDA;
                else if ( cmd == "-g" )
                    status = NUM_GAMMA_CATS;
                else if ( cmd == "-e" )
                    status = ASRV_LAMBDA;
                else
                    {
                    std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
                    exit(1);
                    }
                }
            else
                {
                /* read the parameter */
                if ( status == DATA_FILE )
                    inputFileName = commandString[i];
                else if ( status == OUTPUT_FILE )
                    outPutFileName = commandString[i];
                else if ( status == TREE_FILE )
                    treeFileName = commandString[i];
                else if ( status == TIPTIMES_FILE )
                    tipTimesFileName = commandString[i];
                else if ( status == CHAIN_LENGTH )
                    chainLength = atoi(commandString[i].c_str());
                else if ( status == PRINT_FREQ )
                    printFrequency = atoi(commandString[i].c_str());
                else if ( status == SAMPLE_FREQ )
                    sampleFrequency = atoi(commandString[i].c_str());
                else if ( status == IS_REVERSIBLE )
                    {
                    if (commandString[i][0] == 'Y' || commandString[i][0] == 'y')
                        isReversible = true;
                    else if (commandString[i][0] == 'N' || commandString[i][0] == 'n')
                        isReversible = false;
                    else
                        printUsage();
                    }
                else if ( status == BRLEN_LAMBDA )
                    brlenLambda = atof(commandString[i].c_str());
                else if ( status == NUM_GAMMA_CATS )
                    numGammaCats = atoi(commandString[i].c_str());
                else if ( status == ASRV_LAMBDA )
                    asrvLambda = atof(commandString[i].c_str());
                else
                    {
                    std::cerr << "Unknown status reading command line information" << std::endl;
                    exit(1);
                    }
                status = NONE;
                }
            }
        }
    else
        {
        printCommandString(commandString);
        printUsage();
        }
}

void Settings::printCommandString(std::vector<std::string> cs) {

    for (int i=0; i<cs.size(); i++)
        std::cout << "command[" << i << "] = \"" << cs[i] << "\"" << std::endl;
}

void Settings::printUsage(void) {

    std::cout << "Usage:" << std::endl;
    std::cout << "     -i : Input file name" << std::endl;
    std::cout << "     -t : Tree file name (for constraining the analysis to a fixed tree)" << std::endl;
    std::cout << "     -C : File name with tip dates" << std::endl;
    std::cout << "     -o : Output file name" << std::endl;
    std::cout << "     -l : Number of MCMC cycles" << std::endl;
    std::cout << "     -p : Print frequency" << std::endl;
    std::cout << "     -s : Sample frequency" << std::endl;
    std::cout << "    -tr : Is the model time reversible (yes or no)" << std::endl;
    std::cout << "-lambda : Exponential parameter for branch lengths" << std::endl;
    std::cout << "     -g : Number of gamma rate categories" << std::endl;
    std::cout << "     -e : Exponential parameter for shape parameter describing ASRV" << std::endl;
    std::cout << std::endl;
    std::cout << "Example:" << std::endl;
    std::cout << "   ./rooter -i <input file> -o <output file>" << std::endl;
    exit(0);
}

