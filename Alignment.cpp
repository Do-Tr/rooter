#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include "Alignment.h"
#include "Msg.h"



Alignment::Alignment(std::string fileName) {

	// open the file
	std::ifstream seqStream(fileName.c_str());
	if (seqStream.is_open() == true)
        std::cout << "   * Reading data file \"" << fileName << "\"" << std::endl;
    else
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	int taxonNum = 0;
	matrix = NULL;
	numTaxa = numSites = 0;
    while ( getline (seqStream, linestring) )
		{
		std::istringstream linestream(linestring);
        //std::cout << line << " -- \"" << linestring << "\"" << std::endl;
		int ch;
		std::string word = "";
		int wordNum = 0;
		int siteNum = 0;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
            //std::cout << "word:" << wordNum << "\"" << word << "\"" << std::endl;
			if (line == 0)
				{
				// read the number of taxa/chars from the first line
				if (wordNum == 1)
					numTaxa = atoi(word.c_str());
				else
					numSites = atoi(word.c_str());
				if (numTaxa > 0 && numSites > 0 && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numSites];
					for (size_t i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numSites;
					for (size_t i=0; i<numTaxa; i++)
						for (size_t j=0; j<numSites; j++)
							matrix[i][j] = 0;
                    std::cout << "   * Alignment has " << numTaxa << " taxa and " << numSites << " sites" << std::endl;
					}
				}
			else
				{
				if (wordNum == 1)
					{
                    taxonNames.push_back(word);
                    taxonNum++;
					}
				else
					{
                    for (int i=0; i<word.length(); i++)
                        {
                        char site = word.at(i);
                        matrix[taxonNum-1][siteNum++] = nucID(site);
                        }
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		line++;
		}	
	
	// close the file
	seqStream.close();
 
    // compress the data matrix
    compress();
}

Alignment::~Alignment(void) {

	delete [] matrix[0];
	delete [] matrix;
    delete [] compressedMatrix[0];
    delete [] compressedMatrix;
    delete [] numSitesOfPattern;
}

bool Alignment::areSitesOfSamePattern(int idx1, int idx2) {
    
    for (int i=0; i<numTaxa; i++)
        {
        if (matrix[i][idx1] != matrix[i][idx2])
            return false;
        }
    return true;
}

void Alignment::compress(void) {

    // count the number of unique patterns in the data matrix
    numPatterns = numberUniquePatterns();
    std::cout << "   * Compressed alignment has " << numPatterns << " unique site patterns" << std::endl;

    // allocate the compressed matrix
    compressedMatrix = new int*[numTaxa];
    compressedMatrix[0] = new int[numTaxa * numPatterns];
    for (size_t i=1; i<numTaxa; i++)
        compressedMatrix[i] = compressedMatrix[i-1] + numPatterns;
    for (size_t i=0; i<numTaxa; i++)
        for (size_t j=0; j<numPatterns; j++)
            compressedMatrix[i][j] = 0;
    numSitesOfPattern = new int[numPatterns];
    for (int i=0; i<numPatterns; i++)
        numSitesOfPattern[i] = 0;
    
    // fill in the matrix
    std::vector<bool> touched(numSites, false);
    for (int j=0, k=0; j<numSites; j++)
        {
        if (touched[j] == false)
            {
            touched[j] = true;
            numSitesOfPattern[k]++;
            for (int i=0; i<numTaxa; i++)
                compressedMatrix[i][k] = matrix[i][j];
            for (int m=j+1; m<numSites; m++)
                {
                if ( areSitesOfSamePattern(j, m) == true )
                    {
                    touched[m] = true;
                    numSitesOfPattern[k]++;
                    }
                }
            k++;
            }
        }
}

int Alignment::getNucleotide(size_t i, size_t j) {

    return matrix[i][j];
}

char Alignment::getNucleotideChar(size_t i, size_t j) {

    int nucCode = getNucleotide(i,j);
    char nc = nucChar(nucCode);
    return nc;
}

int Alignment::getNucleotideForPattern(size_t i, size_t j) {

    return compressedMatrix[i][j];
}

char Alignment::getNucleotideCharForPattern(size_t i, size_t j) {

    int nucCode = getNucleotideForPattern(i,j);
    char nc = nucChar(nucCode);
    return nc;
}

std::vector<std::string> Alignment::getTaxonNames(void) {

    return taxonNames;
}

/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Alignment::getPossibleNucs (int nucCode, int* nuc) {

	if (nucCode == 1)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 2)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 3)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 0;
		}
	else if (nucCode == 4)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 5)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 6)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 7)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 0;
		}
	else if (nucCode == 8)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 9)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 10)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 11)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 0;
		nuc[3] = 1;
		}
	else if (nucCode == 12)
		{
		nuc[0] = 0;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 13)
		{
		nuc[0] = 1;
		nuc[1] = 0;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 14)
		{
		nuc[0] = 0;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 15)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
	else if (nucCode == 16)
		{
		nuc[0] = 1;
		nuc[1] = 1;
		nuc[2] = 1;
		nuc[3] = 1;
		}
}

void Alignment::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << '\n';
}

std::string Alignment::getTaxonName(int i) {

	return taxonNames[i];
}

int Alignment::getTaxonIndex(std::string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;
}

char Alignment::nucChar(int nucCode) {

    switch (nucCode)
        {
        case 1: return 'A';
        case 2: return 'C';
        case 3: return 'M';
        case 4: return 'G';
        case 5: return 'R';
        case 6: return 'S';
        case 7: return 'V';
        case 8: return 'T';
        case 9: return 'W';
        case 10: return 'Y';
        case 11: return 'H';
        case 12: return 'K';
        case 13: return 'D';
        case 14: return 'B';
        case 15: return 'N';
        default: return '?';
        }
}

/*-------------------------------------------------------------------
|
|   NucID: 
|
|   Take a character, nuc, and return an integer:
|
|       nuc        returns
|        A            1 
|        C            2     
|        G            4      
|        T U          8     
|        R            5      
|        Y           10       
|        M            3      
|        K           12   
|        S            6     
|        W            9      
|        H           11      
|        B           14     
|        V            7      
|        D           13  
|        N - ?       15       
|
-------------------------------------------------------------------*/
int Alignment::nucID(char nuc) {

	char		n;
	
	if (nuc == 'U' || nuc == 'u')
		n = 'T';
	else
		n = nuc;

	if (n == 'A' || n == 'a')
		{
		return 1;
		}
	else if (n == 'C' || n == 'c')
		{
		return 2;
		}
	else if (n == 'G' || n == 'g')
		{
		return 4;
		}
	else if (n == 'T' || n == 't')
		{
		return 8;
		}
	else if (n == 'R' || n == 'r')
		{
		return 5;
		}
	else if (n == 'Y' || n == 'y')
		{
		return 10;
		}
	else if (n == 'M' || n == 'm')
		{
		return 3;
		}
	else if (n == 'K' || n == 'k')
		{
		return 12;
		}
	else if (n == 'S' || n == 's')
		{
		return 6;
		}
	else if (n == 'W' || n == 'w')
		{
		return 9;
		}
	else if (n == 'H' || n == 'h')
		{
		return 11;
		}
	else if (n == 'B' || n == 'b')
		{
		return 14;
		}
	else if (n == 'V' || n == 'v')
		{
		return 7;
		}
	else if (n == 'D' || n == 'd')
		{
		return 13;
		}
	else if (n == 'N' || n == 'n')
		{
		return 15;
		}
	else if (n == '-')
		{
		return 15;
		}
	else if (n == '?')
		{
		return 15;
		}
	else
		return -1;
}

int Alignment::numberUniquePatterns(void) {

    if (matrix == NULL)
        return 0;

    std::vector<bool> touched(numSites, false);
    int numPat = 0;
    for (int i=0; i<numSites; i++)
        {
        if (touched[i] == false)
            {
            touched[i] = true;
            numPat++;
            for (int j=i+1; j<numSites; j++)
                {
                if ( areSitesOfSamePattern(i, j) == true )
                    touched[j] = true;
                }
            }
        }

    return numPat;
}

void Alignment::print(void) {

	int** x = matrix;
		
	std::cout << "        ";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << std::setw(3) << i;
	std::cout << '\n';
	std::cout << "------------------------";
	for (size_t i=0; i<numTaxa; i++)
		std::cout << "---";
	std::cout << '\n';
	for (size_t j=0; j<numSites; j++)
		{
		std::cout << std::setw(4) << j+1 << " -- ";
		for (size_t i=0; i<numTaxa; i++)
			{
			std::cout << std::setw(3) << x[i][j];
			}
		std::cout << '\n';
		}
    std::cout << std::endl;
    
    x = compressedMatrix;
    
    std::cout << "        ";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << std::setw(3) << i;
    std::cout << '\n';
    std::cout << "------------------------";
    for (size_t i=0; i<numTaxa; i++)
        std::cout << "---";
    std::cout << '\n';
    for (size_t j=0; j<numPatterns; j++)
        {
        std::cout << std::setw(4) << j+1 << " -- ";
        for (size_t i=0; i<numTaxa; i++)
            {
            std::cout << std::setw(3) << x[i][j];
            }
        std::cout << " (" << numSitesOfPattern[j] << ")";
        std::cout << '\n';
        }
}


