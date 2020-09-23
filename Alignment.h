#ifndef Alignment_H
#define Alignment_H

#include <string>
#include <vector>


class Alignment {

    public:
                                    Alignment(void) = delete;
                                    Alignment(std::string fileName);
                                   ~Alignment(void);
        int                         getNumTaxa(void) { return numTaxa; }
        int                         getNumSites(void) { return numSites; }
        int                         getNumPatterns(void) { return numPatterns; }
        void                        getPossibleNucs (int nucCode, int* nuc);
        int                         getNucleotide(size_t i, size_t j);
        char                        getNucleotideChar(size_t i, size_t j);
        int                         getNucleotideForPattern(size_t i, size_t j);
        char                        getNucleotideCharForPattern(size_t i, size_t j);
        int                         getNumSitesOfPattern(int idx) { return numSitesOfPattern[idx]; }
        int                         getTaxonIndex(std::string ns);
        std::vector<std::string>    getTaxonNames(void);
        std::string                 getTaxonName(int i);
        void                        listTaxa(void);
        void                        print(void);

    protected:
        bool                        areSitesOfSamePattern(int idx1, int idx2);
        void                        compress(void);
        char                        nucChar(int nucCode);
        int                         nucID(char nuc);
        int                         numberUniquePatterns(void);
        std::vector<std::string>    taxonNames;
        int**                       matrix;
        int**                       compressedMatrix;
        int*                        numSitesOfPattern;
        int                         numTaxa;
        int                         numSites;
        int                         numPatterns;
};

#endif
