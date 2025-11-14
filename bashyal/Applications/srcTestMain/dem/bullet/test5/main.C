#include "fvCFD.H"

#include "particlesGenerator.H"
#include "bulletIntegrator.H"
#include "timeRegistry.H"
#include "packer.H"

using namespace Foam;

int main(int argc, char* argv[])
{
    // Standalone app (no setRootCase / createTime)

    const scalar endTime = 6.0;
    const scalar dt = 0.003;
    const scalar writeInterval = 0.03;
    const vector gravity(0, 0, -9.81);

    const std::string outDir = "/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/dem/bullet/test5/VTK_output";

    // Define packing domain as an axis-aligned box
    point boxMin(-0.5, -0.5, 0.0);
    point boxMax( 0.5,  0.5, 0.8);

    // Use packer to place particles without overlap inside the box
    Bashyal::packer pack(boxMin, boxMax);

    const label nParticles = 60;       // increase count for dense packing
    const scalar radius = 0.035;       // slightly smaller radius to reduce initial stress
    const int resolution = 16;
    PtrList<Bashyal::particle> parts = pack.packSpheres(nParticles, radius, resolution, 2000);

    // Register for outputs
    Bashyal::timeRegistry tr(dt, endTime, outDir);
    tr.setWriteInterval(writeInterval);
    forAll(parts, i)
    {
        tr.addParticle(parts[i], "particle_" + name(i));
    }

    // Build Bullet world
    Bashyal::bulletIntegrator integrator(gravity);

    // Add static box walls (not affected by gravity)
    vector half((boxMax.x()-boxMin.x())/2.0, (boxMax.y()-boxMin.y())/2.0, (boxMax.z()-boxMin.z())/2.0);
    point center((boxMin.x()+boxMax.x())/2.0, (boxMin.y()+boxMax.y())/2.0, (boxMin.z()+boxMax.z())/2.0);
    scalar thickness = 0.05; // thicker walls to reduce tunneling

    // Floor and ceiling with higher friction and lower restitution
    integrator.addStaticBox(point(center.x(), center.y(), boxMin.z()+thickness), vector(half.x(), half.y(), thickness), quaternion::I, 0.05, 0.95);
    integrator.addStaticBox(point(center.x(), center.y(), boxMax.z()-thickness), vector(half.x(), half.y(), thickness), quaternion::I, 0.05, 0.95);
    // +X / -X walls
    integrator.addStaticBox(point(boxMax.x()-thickness, center.y(), center.z()), vector(thickness, half.y(), half.z()), quaternion::I, 0.05, 0.95);
    integrator.addStaticBox(point(boxMin.x()+thickness, center.y(), center.z()), vector(thickness, half.y(), half.z()), quaternion::I, 0.05, 0.95);
    // +Y / -Y walls
    integrator.addStaticBox(point(center.x(), boxMax.y()-thickness, center.z()), vector(half.x(), thickness, half.z()), quaternion::I, 0.05, 0.95);
    integrator.addStaticBox(point(center.x(), boxMin.y()+thickness, center.z()), vector(half.x(), thickness, half.z()), quaternion::I, 0.05, 0.95);

    // Register particles with low restitution to avoid bouncing; higher friction for faster settling
    integrator.registerFrom(tr, 0.1, 0.9);

    Info << "Starting DEM packing simulation (test5) with " << nParticles << " particles" << endl;
    scalar t = 0.0;
    while (t < endTime - SMALL)
    {
        integrator.step(dt, 5);
        tr.advanceTime();
        t += dt;
    }

    tr.writeTimeSeries();
    Info << "Finished test5. Outputs in " << outDir << endl;
    return 0;
}


