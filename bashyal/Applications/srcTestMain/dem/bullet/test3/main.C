#include "fvCFD.H"

#include "particlesGenerator.H"
#include "bulletIntegrator.H"
#include "timeRegistry.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    // #include "setRootCase.H"
    // #include "createTime.H"

    // Simulation parameters
    const scalar endTime = 2.0;        // seconds
    const scalar dt = 0.01;            // seconds
    const scalar writeInterval = 0.05; // seconds
    const vector gravity(0, 0, -9.81);

    // Output directory (absolute) inside this test folder
    const std::string outDir = "/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/dem/bullet/test3/VTK_output";

    // Generators and registries
    Bashyal::particlesGenerator gen;
    Bashyal::timeRegistry tr(dt, endTime, outDir);
    tr.setWriteInterval(writeInterval);

    // Create particles of various shapes/sizes
    PtrList<Bashyal::particle> spheres = gen.createSpherical(6, 0.08, 12);
    PtrList<Bashyal::particle> cubes   = gen.createCubical(6, 0.08, 0.08);
    PtrList<Bashyal::particle> cuboids = gen.createCuboid(6, vector(0.06, 0.08, 0.10), vector(0.06, 0.08, 0.10));
    PtrList<Bashyal::particle> rounds  = gen.createWellRounded(6, 0.07, 0.07, 0.7, 12);

    // Place particles in a simple grid (non-overlapping by construction)
    const scalar spacing = 0.25; // grid spacing
    const point origin(-0.5, -0.5, 0.6);
    auto placeGrid = [&](PtrList<Bashyal::particle>& list, const point& base, const word& prefix)
    {
        label idx = 0;
        forAll(list, i)
        {
            label gx = idx % 3;
            label gy = (idx / 3) % 3;
            label gz = idx / 9;
            point p = base + vector(gx * spacing, gy * spacing, gz * spacing);
            list[i].setPosition(p);
            tr.addParticle(list[i], prefix + "_" + name(i));
            ++idx;
        }
    };

    placeGrid(spheres, origin + vector(0.0, 0.0, 0.0), "sphere");
    placeGrid(cubes,   origin + vector(0.0, 0.8, 0.0), "cube");
    placeGrid(cuboids, origin + vector(0.8, 0.0, 0.0), "cuboid");
    placeGrid(rounds,  origin + vector(0.8, 0.8, 0.0), "rounded");

    // Bullet integrator setup
    Bashyal::bulletIntegrator integrator(gravity);
    integrator.registerFrom(tr, 0.5, 0.5);

    Info << "Starting DEM Bullet simulation with "
         << (spheres.size() + cubes.size() + cuboids.size() + rounds.size())
         << " particles" << endl;

    // Time loop
    scalar t = 0.0;
    while (t < endTime - SMALL)
    {
        integrator.step(dt, 5);
        tr.advanceTime();
        t += dt;
    }

    // Finalize: write PVD collection
    tr.writeTimeSeries();

    Info << "Finished. Outputs in " << outDir << endl;
    return 0;
}


