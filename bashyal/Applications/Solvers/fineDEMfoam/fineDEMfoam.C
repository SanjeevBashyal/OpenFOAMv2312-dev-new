/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Group
    grpMultiphaseSolvers

Description
    Solver for two incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include <fstream>
#include <iostream>
#include <iomanip>

// DEM Includes
#include "particle.H"
#include "particlesGenerator.H"
#include "timeRegistry.H"
#include "gjkContact.H"
#include "contactForce.H"
#include "wall.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "createPimpleControl.H"
    #include "createFields.H"
    #include "createAlphaFluxes.H"
    
    bool correctPhi = false;
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // DEM Initialization
    Info<< "Initializing DEM..." << endl;
    
    // --- Pre-Run Settling Phase ---
    Info << "\n========================================\n"
         << "Pre-Run Settling Phase (DEM Only)\n"
         << "========================================\n" << endl;

    // 1. Setup DEM Objects
    Bashyal::gjkContact contactDetector;
    scalar k_dem = 1e5;
    scalar mu_dem = 0.5;
    scalar gammaN_dem = 50.0;
    scalar crustWidth = 1e-4;
    Bashyal::contactForce forceComputer = Bashyal::contactForce(k_dem, k_dem, k_dem, gammaN_dem, 0.0, mu_dem);
    
    // Walls
    Bashyal::wall floor(point(0,0,0), vector(0,0,1), "floor");
    // Side walls (optional, but good for containment)
    // For this channel flow, maybe just floor is enough if we rely on cyclic/slip for sides?
    // But user asked for "loosely packed 50 particles", they need a container to stack.
    // Let's add side walls for the settling phase, then maybe remove them or ignore them for CFD?
    // Actually, for the CFD part, the mesh has boundaries. For DEM, we need explicit walls if particles hit them.
    // Let's add a floor and maybe periodic boundaries or walls matching the domain.
    // Domain is 0 to 2.0 in X, -0.1 to 0.1 in Y (width), 0 to 0.3 in Z (height).
    // Wait, previous code used Z as vertical? No, OpenFOAM usually Y is vertical.
    // Let's check gravity. constant/g had (0 -9.81 0) usually.
    // In `incipientMotion`, we set gravity.
    // Let's assume Y is vertical (gravity direction).
    // Domain: X (0 to 2), Y (0 to 0.3), Z (0 to 0.2).
    // Floor at Y=0.
    
    Bashyal::wall wallFloor(point(0,0,0), vector(0,1,0), "floor");
    
    // 2. Particles
    Foam::PtrList<Bashyal::particle> particles;
    particles.setSize(51); // 1 Wall + 50 Dynamic
    
    // Wall Particle (Index 0) - Kept for visualization/consistency with previous plan
    particles.set(0, new Bashyal::particle());
    particles[0].setPosition(point(1.0, -0.05, 0.1)); // Below floor
    particles[0].setDC(0.1);
    particles[0].setMass(1000.0);
    particles[0].setBodyOperation(0); // Static
    
    // 3. Settling Loop
    scalar demDt = 1e-4;
    scalar demTime = 0.0;
    scalar settleTime = 2.0; // Max time
    scalar lastEmitTime = -0.01;
    scalar emitInterval = 0.02;
    int pCount = 1;
    int maxParticles = 50;
    
    scalar D = 0.02;
    scalar rho_s = 2500.0;
    scalar mass = rho_s * pow3(D);
    
    Foam::Random rnd(1234);
    
    Info << "Pouring 50 particles..." << endl;
    
    while (demTime < settleTime)
    {
        // Emission
        if (pCount <= maxParticles && (demTime - lastEmitTime) >= emitInterval)
        {
            scalar x = 1.0 + (rnd.sample01<scalar>() - 0.5) * 0.2; // Center X=1.0, width 0.2
            scalar z = 0.1 + (rnd.sample01<scalar>() - 0.5) * 0.15; // Center Z=0.1, width 0.15
            scalar y = 0.2; // Drop from height
            
            point pos(x, y, z);
            
            // Check overlap
            bool overlap = false;
            for (int i=1; i<pCount; ++i)
            {
                if (mag(particles[i].position() - pos) < D) overlap = true;
            }
            
            if (!overlap)
            {
                particles.set(pCount, new Bashyal::particle());
                particles[pCount].setPosition(pos);
                particles[pCount].setDC(D);
                particles[pCount].setMass(mass);
                particles[pCount].setBodyOperation(5); // Free motion
                
                // Random orientation? (Optional, particle class handles it?)
                
                pCount++;
                lastEmitTime = demTime;
                // Info << "Emitted " << pCount-1 << endl;
            }
        }
        
        // Forces
        for (int i=1; i<pCount; ++i) particles[i].clearForceAndTorque();
        
        // Gravity
        vector gVec(0, -9.81, 0);
        for (int i=1; i<pCount; ++i) particles[i].applyForce(mass * gVec, particles[i].position());
        
        // Contacts
        for (int i=1; i<pCount; ++i)
        {
            // Floor
            Bashyal::ContactInfo c = contactDetector.detectContactParticleWall(particles[i], wallFloor, crustWidth);
            if (c.isContact) forceComputer.computeForceParticleWall(particles[i], wallFloor, c, demDt);
            
            // Particle-Particle
            for (int j=i+1; j<pCount; ++j)
            {
                 Bashyal::ContactInfo cp = contactDetector.detectContact(particles[i], particles[j], crustWidth, crustWidth);
                 if (cp.isContact) forceComputer.computeForce(particles[i], particles[j], cp, demDt);
            }
        }
        
        // Update
        scalar totalKE = 0.0;
        for (int i=1; i<pCount; ++i)
        {
            particles[i].update(demDt);
            totalKE += 0.5 * mass * magSqr(particles[i].velocity());
        }
        
        demTime += demDt;
        
        // Exit condition
        if (pCount > maxParticles && totalKE < 1e-4)
        {
            Info << "Settled at t=" << demTime << " KE=" << totalKE << endl;
            break;
        }
        
        if (int(demTime/demDt) % 1000 == 0)
        {
             Info << "DEM Time: " << demTime << " Particles: " << pCount-1 << " KE: " << totalKE << endl;
        }
    }
    
    // Resize to actual count if not full (unlikely)
    particles.setSize(pCount);
    
    Info << "Settling complete. Starting CFD..." << endl;

    std::ofstream dataFile("sedimentation_data.csv");
    dataFile << "Time,X,Y,Z,Vx,Vy,Vz" << std::endl;

    // Initialize Time Registry
    // Note: runTime.endTime().value() might be large, but we can update it
    Bashyal::timeRegistry timeReg(runTime.deltaT().value(), runTime.endTime().value(), "VTK");
    timeReg.setWriteInterval(runTime.controlDict().get<scalar>("writeInterval")); // Sync with OpenFOAM write interval
    
    if (particles.size() > 0)
    {
        timeReg.addParticle(particles[0], "particle1");
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "alphaCourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // Calculate IBM fields
        #include "calculateIBM.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter())
            {
                // mesh.update(); // Removed for static mesh

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                }
            }

            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
            
            // Coupling
            if (pimple.finalIter())
            {
                #include "coupling.H"
            }
        }

        runTime.write();

        // Update time registry
        timeReg.setCurrentTime(runTime.value());
        timeReg.advanceTime(); // This handles writing VTPs if interval is met
        timeReg.writeTimeSeries(); // Update PVD file

        if (particles.size() > 0)
        {
             dataFile << runTime.timeName() << "," 
                      << particles[0].position().x() << "," 
                      << particles[0].position().y() << "," 
                      << particles[0].position().z() << ","
                      << particles[0].velocity().x() << "," 
                      << particles[0].velocity().y() << "," 
                      << particles[0].velocity().z() 
                      << std::endl;
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
