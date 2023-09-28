// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Random.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "Utility/IpplTimings.h"
#include "Ippl.h"
#include "Random/Distribution.h"
#include "Random/InverseTransformSampling.h"
#include "Random/NormalDistribution.h"

struct custom_cdf{
       KOKKOS_INLINE_FUNCTION double operator()(double x, unsigned int d, const double *params) const {
           if(d==0){
               return ippl::random::uniform_cdf_func<double>(x);
           }
           else{
               return x + (params[2] / params[3]) * Kokkos::sin(params[3] * x);
           }
       }
};
struct custom_pdf{
       KOKKOS_INLINE_FUNCTION double operator()(double x, unsigned int d, double const *params) const {
           if(d==0){
               return ippl::random::uniform_pdf_func<double>();;
           }
           else{
               return  1.0 + params[2] * Kokkos::cos(params[3] * x);
           }
       }
};
struct custom_estimate{
        KOKKOS_INLINE_FUNCTION double operator()(double u, unsigned int d, double const *params) const {
            if(d==0){
                return ippl::random::uniform_estimate_func<double>(u+params[0]*0);
            }
            else{
                return u;
            }
        }
};

KOKKOS_FUNCTION unsigned int doublefactorial(unsigned int n)
{
    if (n == 0 || n==1)
      return 1;
    return n*doublefactorial(n-2);
}

KOKKOS_FUNCTION double NormDistCentMom(double stdev, unsigned int p){
    if(p%2==0){
        return pow(stdev, p)*doublefactorial(p-1);
    }
    else{
        return 0.;
    }
}

KOKKOS_FUNCTION void NormDistCentMoms(double stdev, const int P, double *moms){
    for(int p=0; p<P; p++){
        moms[p] = NormDistCentMom(stdev, p+1);
    }
}

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
        const int Dim = 2;

        using Mesh_t = ippl::UniformCartesian<double, Dim>;

        ippl::Vector<int, 2> nr   = {std::atoi(argv[1]), std::atoi(argv[2])};
        const unsigned int ntotal = std::atol(argv[3]);

        ippl::NDIndex<2> domain;
        for (unsigned i = 0; i < Dim; i++) {
            domain[i] = ippl::Index(nr[i]);
        }

        ippl::e_dim_tag decomp[Dim];
        for (unsigned d = 0; d < Dim; ++d) {
            decomp[d] = ippl::PARALLEL;
        }

        ippl::Vector<double, Dim> rmin   = -4.;
        ippl::Vector<double, Dim> rmax   = 4.;
        ippl::Vector<double, Dim> length = rmax - rmin;

        ippl::Vector<double, Dim> hr     = length / nr;
        ippl::Vector<double, Dim> origin = rmin;

        const bool isAllPeriodic = true;

        Mesh_t mesh(domain, hr, origin);

        ippl::FieldLayout<Dim> fl(domain, decomp, isAllPeriodic);

        ippl::detail::RegionLayout<double, Dim, Mesh_t> rlayout(fl, mesh);

        using view_type  = typename ippl::detail::ViewType<ippl::Vector<double, Dim>, 1>::view_type;
        int seed = 42;
        using size_type = ippl::detail::size_type;
        unsigned int nlocal;
        Kokkos::Random_XorShift64_Pool<> rand_pool64((size_type)(seed + 100 * ippl::Comm->rank()));

        // example of sampling normal in both dimensions
        const double mu1 = 0.8;
        const double sd1 = 0.3;
        const double mu2 = -mu1;
        const double sd2 = sqrt(2.-( sd1*sd1 + 2.*mu1*mu1 ) );
        const double par[4] = {mu1, sd1, mu2, sd2};
        using Dist_t = ippl::random::NormalDistribution<double, Dim>;
        using sampling_t = ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace, Dist_t>;

        
        Dist_t dist(par);
        sampling_t sampling(dist, rmax, rmin, rlayout, ntotal);
        nlocal = sampling.getLocalNum();
        view_type position("position", nlocal);
        sampling.generate(position, rand_pool64);
        
        const int P = 4;
        double moms1[P], moms2[P];
        NormDistCentMoms(sd1, P, moms1);
        NormDistCentMoms(sd2, P, moms2);
        for(int i=0; i<P; i++)
            std::cout << moms1[i] << std::endl;

        for(int i=0; i<P; i++)
            std::cout << moms2[i] << std::endl;
        /*
        // example of sampling normal/uniform in one and harmonic in another with custom functors
        const int DimP = 4;
        double pi    = Kokkos::numbers::pi_v<double>;
        using DistH_t = ippl::random::Distribution<double, Dim, DimP, custom_pdf, custom_cdf, custom_estimate>;
        using samplingH_t = ippl::random::InverseTransformSampling<double, Dim, Kokkos::DefaultExecutionSpace, DistH_t>;
        const double parH[DimP] = {rmin[0], rmax[0], 0.5, 2.*pi/(rmax[1]-rmin[1])*4.0};

        DistH_t distH(parH);
        samplingH_t samplingH(distH, rmax, rmin, rlayout, ntotal);
        nlocal = samplingH.getLocalNum();
        view_type positionH("positionH", nlocal);
        samplingH.generate(positionH, rand_pool64);
        */
    }
    ippl::finalize();
    return 0;
}
