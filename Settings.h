#ifndef Settings_H
#define Settings_H

#include <string>
#include <vector>



class Settings {

    public:
                        Settings(int argc, char* argv[]);
        double          getAsrvLambda(void) { return asrvLambda; }
        double          getBrlenLambda(void) { return brlenLambda; }
        int             getChainLength(void) { return chainLength; }
        std::string     getInputFileName(void) { return inputFileName; }
        bool            getIsReversible(void) { return isReversible; }
        int             getNumGammaCats(void) { return numGammaCats; }
        std::string     getOutPutFileName(void) { return outPutFileName; }
        int             getPrintFrequency(void) { return printFrequency; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        std::string     getTipTimesFileName(void) { return tipTimesFileName; }
        std::string     getTreeFileName(void) { return treeFileName; }
        void            setAsrvLambda(double x) { asrvLambda = x; }
        void            setInputFileName(std::string s) { inputFileName = s; }
        void            setIsReversible(bool tf) { isReversible = tf; }
        void            setNumGammaCats(int x) { numGammaCats = x; }
        void            setOutPutFileName(std::string s) { outPutFileName = s; }
        void            setTipTimesFileName(std::string s) { tipTimesFileName = s; }
        void            setTreeFileName(std::string s ) { treeFileName = s; }

    private:
        void            printCommandString(std::vector<std::string> cs);
        void            printUsage(void);
        double          brlenLambda;
        std::string     inputFileName;
        std::string     outPutFileName;
        std::string     treeFileName;
        std::string     tipTimesFileName;
        double          asrvLambda;
        int             chainLength;
        int             printFrequency;
        int             sampleFrequency;
        int             numGammaCats;
        bool            isReversible;
};

#endif
