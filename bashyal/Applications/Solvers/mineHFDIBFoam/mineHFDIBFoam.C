#include "introHeader.H"

// OpenFOAM Default include
#include "fvCFD_rest.H"
// bashyal include
#include "sediment_rest.H"

using namespace Bashyal;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "openfoamObjectsCreate.H"

    turbulence->validate();

    #include "readDynMeshDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nInitializing HFDIBDEM\n" << endl;
    openHFDIBDEM  HFDIBDEM(mesh);
    HFDIBDEM.initialize(lambda,U,refineF,maxRefinementLevel,runTime.timeName());
    #include "initialMeshRefinement.H"
    
    if(HFDIBDEM.getRecordFirstTime())
    {
        HFDIBDEM.setRecordFirstTime(false);
        HFDIBDEM.writeBodiesInfo();
    }

    Info<< "\nStarting time loop\n" << endl;

    scalar CFDTime_(0.0);
    scalar DEMTime_(0.0);
    scalar suplTime_(0.0);

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        clockTime createBodiesTime; // OS time efficiency testing
        HFDIBDEM.createBodies(lambda,refineF);
        suplTime_ += createBodiesTime.timeIncrement(); // OS time efficiency testing
        
        clockTime preUpdateBodiesTime; // OS time efficiency testing
        HFDIBDEM.preUpdateBodies(lambda,f);
        suplTime_ += preUpdateBodiesTime.timeIncrement(); // OS time efficiency testing

        clockTime pimpleRunClockTime; // OS time efficiency testing
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.controlledUpdate();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    lambda *= 0.0;
                    HFDIBDEM.recreateBodies(lambda,refineF);
                }

                f *= lambda;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        CFDTime_ += pimpleRunClockTime.timeIncrement();
        Info << "updating HFDIBDEM" << endl;
        clockTime postUpdateBodiesTime;
        HFDIBDEM.postUpdateBodies(lambda,f);
        suplTime_ += postUpdateBodiesTime.timeIncrement();


        clockTime addRemoveTime;
        HFDIBDEM.addRemoveBodies(lambda,U,refineF);
        suplTime_ += addRemoveTime.timeIncrement();

        clockTime updateDEMTime;
        HFDIBDEM.updateDEM(lambda,refineF);
        DEMTime_ += updateDEMTime.timeIncrement();
        Info << "updated HFDIBDEM" << endl;


        runTime.write();

        clockTime writeBodiesInfoTime;
        if(runTime.outputTime())
        {
            HFDIBDEM.writeBodiesInfo();
        }
        suplTime_ += writeBodiesInfoTime.timeIncrement();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    Info<< " CFDTime_                 = " << CFDTime_             << " s \n" <<
           " Solver suplementary time = " << suplTime_            << " s \n" << 
           " DEMTime_                 = " << DEMTime_             << " s \n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
