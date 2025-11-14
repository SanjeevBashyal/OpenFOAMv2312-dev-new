#include "fvCFD.H"
#include "DEMSimulation.H"
#include "createMeshes.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // --- Simulation Parameters ---
    const Foam::scalar endTime = 5.0;       // seconds
    const Foam::scalar dt = 0.005;          // seconds
    const Foam::scalar writeInterval = 0.05; // seconds
    const Foam::vector gravity(0, 0, -9.81);

    // --- Setup the Simulation ---
    Bashyal::DEMSimulation simulation(gravity);

    // Add particles as defined in the helper file
    addInitialParticles(simulation);

    // --- Run the Simulation ---
    Foam::Info << "\nStarting DEM Simulation" << Foam::endl;
    simulation.run(endTime, dt, writeInterval);

    Foam::Info << "End of application" << Foam::endl;
    return 0;
}