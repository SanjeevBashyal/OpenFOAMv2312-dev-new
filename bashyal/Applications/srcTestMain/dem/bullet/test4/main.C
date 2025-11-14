#include "fvCFD.H"

#include "particlesGenerator.H"
#include "bulletIntegrator.H"
#include "timeRegistry.H"

using namespace Foam;

static Foam::quaternion randomQuaternion(Foam::Random& rng)
{
    scalar u1 = rng.sample01<scalar>();
    scalar u2 = rng.sample01<scalar>();
    scalar u3 = rng.sample01<scalar>();
    scalar sq1 = Foam::sqrt(1 - u1);
    scalar sq2 = Foam::sqrt(u1);
    scalar theta1 = 2 * constant::mathematical::pi * u2;
    scalar theta2 = 2 * constant::mathematical::pi * u3;
    scalar w = sq2 * Foam::cos(theta2);
    vector v(sq1 * Foam::cos(theta1), sq1 * Foam::sin(theta1), sq2 * Foam::sin(theta2));
    Foam::quaternion q(w, v);
    q.normalise();
    return q;
}

int main(int argc, char* argv[])
{
    // #include "setRootCase.H"
    // #include "createTime.H"

    const scalar endTime = 5.0;
    const scalar dt = 0.0025;
    const scalar writeInterval = 0.02;
    const vector gravity(0, 0, -9.81);

    const std::string outDir = "/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/dem/bullet/test4/VTK_output";

    Bashyal::particlesGenerator gen;
    Bashyal::timeRegistry tr(dt, endTime, outDir);
    tr.setWriteInterval(writeInterval);

    // Generate 20 near-by particles as spheres for robust collision response
    const label nParticles = 20;
    const scalar radius = 0.05;
    const int resolution = 16;
    PtrList<Bashyal::particle> parts = gen.createSpherical(nParticles, radius, resolution);

    // Initialize positions in a compact cluster and give initial velocities
    Foam::Random rng(42);
    point clusterCenter(0.0, 0.0, 0.6);
    scalar spread = 0.15;
    forAll(parts, i)
    {
        vector offset((rng.sample01<scalar>() - 0.5) * spread,
                      (rng.sample01<scalar>() - 0.5) * spread,
                      (rng.sample01<scalar>() - 0.5) * spread);
        parts[i].setPosition(clusterCenter + offset);

        // Initial linear velocity in random direction
        vector v0((rng.sample01<scalar>() - 0.5) * 2.0,
                  (rng.sample01<scalar>() - 0.5) * 2.0,
                  (rng.sample01<scalar>() - 0.5) * 2.0);
        parts[i].setVelocity(v0);

        // Initial angular velocity
        vector w0((rng.sample01<scalar>() - 0.5) * 20.0,
                  (rng.sample01<scalar>() - 0.5) * 20.0,
                  (rng.sample01<scalar>() - 0.5) * 20.0);
        parts[i].setAngularVelocity(w0);

        // Random orientation
        parts[i].setOrientation(randomQuaternion(rng));

        tr.addParticle(parts[i], "particle_" + name(i));
    }

    // Build Bullet integrator and add static container walls (fixed)
    Bashyal::bulletIntegrator integrator(gravity);

    // Register dynamic particles (with initial velocities)
    integrator.registerFrom(tr, 0.4, 0.6);

    // Add static walls: box container centered at origin
    point c(0, 0, 0.3);
    vector half(0.8, 0.8, 0.3);
    scalar thickness = 0.02;

    // Floor
    integrator.addStaticBox(point(0,0,0.0), vector(half.x(), half.y(), thickness), quaternion::I, 0.4, 0.8);
    // Ceiling
    integrator.addStaticBox(point(0,0,0.6), vector(half.x(), half.y(), thickness), quaternion::I, 0.4, 0.8);
    // +X wall
    integrator.addStaticBox(point( half.x(), 0, 0.3), vector(thickness, half.y(), half.z()), quaternion::I, 0.4, 0.8);
    // -X wall
    integrator.addStaticBox(point(-half.x(), 0, 0.3), vector(thickness, half.y(), half.z()), quaternion::I, 0.4, 0.8);
    // +Y wall
    integrator.addStaticBox(point(0,  half.y(), 0.3), vector(half.x(), thickness, half.z()), quaternion::I, 0.4, 0.8);
    // -Y wall
    integrator.addStaticBox(point(0, -half.y(), 0.3), vector(half.x(), thickness, half.z()), quaternion::I, 0.4, 0.8);

    Info << "Starting DEM Bullet simulation (test4) with " << nParticles << " particles" << endl;

    scalar t = 0.0;
    while (t < endTime - SMALL)
    {
        integrator.step(dt, 5);
        tr.advanceTime();
        t += dt;
    }

    tr.writeTimeSeries();
    Info << "Finished test4. Outputs in " << outDir << endl;
    return 0;
}


