#ifndef EigenSystem_H
#define EigenSystem_H

#include <complex>
#include <vector>
#include <Eigen/Core>

typedef Eigen::Matrix<double, 4, 1> NucleotideVector_t;
typedef Eigen::Matrix<double, 4, 4> NucleotideSquareMatrix_t;



class EigenSystem {

    public:
        static EigenSystem&                         eigenSystem(void)
                                                        {
                                                        static EigenSystem singleEigenSystem;
                                                        return singleEigenSystem;
                                                        }
        void                                        calculateEigenSystem(NucleotideSquareMatrix_t& Q);
        std::vector<double>                         calulateStationaryFrequencies(NucleotideSquareMatrix_t& Q);
        std::vector<double>&                        getStationaryFrequencies(void) { return pi[activeVals]; }
        void                                        flipActiveValues(void);
        std::complex<double>*                       getCijk(void);
        Eigen::Matrix<std::complex<double>, 4, 1>&  getEigenValues(void);
        void                                        setStationaryFrequencies(std::vector<double> f) { pi[activeVals] = f; }
        void                                        testStationaryFrequencies(NucleotideSquareMatrix_t& Q);

    private:
                                                    EigenSystem(void);
                                                    EigenSystem(const EigenSystem&);
                                                    EigenSystem& operator=(const EigenSystem&);
                                                   ~EigenSystem(void);
        double                                      sumVector(std::vector<double>& v);
        int                                         activeVals;
        Eigen::Matrix<std::complex<double>, 4, 1>   eigenValues[2];
        std::complex<double>*                       ccIjk[2];
        std::vector<double>                         pi[2];
};

#endif
