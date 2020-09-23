#include <algorithm>
#include <iomanip>
#include <iostream>

#include "BipartitionSamples.h"
#include "RbBitSet.h"

bool sortByVal(const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {

    return (a.second > b.second);
}


BipartitionSamples::BipartitionSamples(void) {

    nunSamples = 0;
}

void BipartitionSamples::addBipartition(std::string bp) {

    nunSamples++;
    
    // has this bipartition already been sampled?
    bool alreadySampled = false;
    for (int i=0; i<bitsets.size(); i++)
        {
        if (bitsets[i] == bp)
            {
            alreadySampled = true;
            numVisits[i]++;
            break;
            }
        }
    
    // add it if it hasn't been sampled
    if (alreadySampled == false)
        {
        bitsets.push_back(bp);
        numVisits.push_back(1);
        }
}

void BipartitionSamples::print(void) {
    
    std::vector<std::pair<std::string, int> > vec;

    // copy key-value pairs from the map to the vector
    for (int i=0; i<bitsets.size(); i++)
        {
        vec.push_back( std::make_pair(bitsets[i], numVisits[i]) );
        }
    
    // // sort the vector by increasing order of its pair's second value
    sort(vec.begin(), vec.end(), sortByVal);

    // print the vector
    std::cout << std::endl;
    std::cout << "   * Sampled Root Bipartitions" << std::endl;
    double prob = 0.0;
    for (int i=0; i<vec.size(); i++)
        {
        prob += (double)vec[i].second / nunSamples;
        
        std::cout << "   * " << std::setw(4) << i+1 << " -- ";
        std::cout << vec[i].first << " ";
        std::cout << std::setw(6) << vec[i].second << " ";
        std::cout << std::fixed << std::setprecision(4) << (double)vec[i].second / nunSamples << " ";
        std::cout << std::fixed << std::setprecision(4) << prob << " ";
        std::cout << std::endl;
        }
    std::cout << std::endl;
}
