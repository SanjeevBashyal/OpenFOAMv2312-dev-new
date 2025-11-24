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
    Bashyal::particlesGenerator pGen;
    // Create 1 particle for testing
    Foam::PtrList<Bashyal::particle> particles = pGen.createCubical(1, 0.05, 0.05);
    
    if (particles.size() > 0)
    {
        // Test Case 5: Floating Cube
        // Domain: 0.3 x 1.0 x 0.3
        // Particle: 0.1m Cube
        // Density: 500 kg/m3
        // Mass = 500 * 0.1^3 = 0.5 kg
        // Position: (0.15, 0.62, 0.15) (Bottom at 0.57, Water at 0.5)
        particles[0].setPosition(point(0.15, 0.62, 0.15));
        particles[0].setDC(0.1); 
        particles[0].setMass(0.5);
    }

    std::ofstream dataFile("sedimentation_data.csv");
    dataFile << "Time,Y,Vy" << std::endl;

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
                      << particles[0].position().y() << "," 
                      << particles[0].velocity().y() << std::endl;
        }

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
