// Landau Damping Test
//   Usage:
//     srun ./LandauDamping
//                  <nx> [<ny>...] <Np> <Nt> <stype>
//                  <lbthres> --overallocate <ovfactor> --info 10
//     nx       = No. cell-centered points in the x-direction
//     ny...    = No. cell-centered points in the y-, z-, ...-direction
//     Np       = Total no. of macro-particles in the simulation
//     Nt       = Number of time steps
//     stype    = Field solver type (FFT and CG supported)
//     lbthres  = Load balancing threshold i.e., lbthres*100 is the maximum load imbalance
//                percentage which can be tolerated and beyond which
//                particle load balancing occurs. A value of 0.01 is good for many typical
//                simulations.
//     ovfactor = Over-allocation factor for the buffers used in the communication. Typical
//                values are 1.0, 2.0. Value 1.0 means no over-allocation.
//     Example:
//     srun ./LandauDamping 128 128 128 10000 10 FFT 0.01 LeapFrog --overallocate 2.0 --info 10
//
// Copyright (c) 2023, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of IPPL.
//
// IPPL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with IPPL. If not, see <https://www.gnu.org/licenses/>.
//

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Random.hpp>
#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "Ippl.h"
#include "Utility/IpplTimings.h"
#include "Manager/PicManager.h"
#include "datatypes.h"
#include "Random/InverseTransformSampling.h"
#include "LandauDampingManager.h"

int main(int argc, char* argv[]) {
    ippl::initialize(argc, argv);
    {
        Inform msg("LandauDamping");
        Inform msg2all("LandauDamping", INFORM_ALL_NODES);
        
        // Create an instance of a manger for the considered application
        LandauDampingManager manager;

        // Read input parameters, assign them to the corresponding memebers of manager
        int arg = 1;
        for (unsigned d = 0; d < Dim; d++) {
            manager.nr[d] = std::atoi(argv[arg++]);
        }
        manager.totalP = std::atoll(argv[arg++]);
        manager.nt  = std::atoi(argv[arg++]);
        manager.solver = argv[arg++];
        manager.lbt = std::atof(argv[arg++]);
        manager.step_method = argv[arg++];
        
        // Perform pre-run operations, including creating mesh, particles,...
       manager.pre_run();
       
       // Loop over time
        manager.time_m = 0.0;
        msg << "Starting iterations ..." << endl;
        for (manager.it = 0; manager.it < manager.nt; manager.it++) {
            // Perfom one step of time integration
            manager.advance();
            
            // Save the output, increment time                       
            manager.post_step();
            
            // Print time
            msg << "Finished time step: " << manager.it + 1 << " time: " << manager.time_m << endl;
        }
        msg << "LandauDamping: End." << endl;
    }
    ippl::finalize();

    return 0;
}
