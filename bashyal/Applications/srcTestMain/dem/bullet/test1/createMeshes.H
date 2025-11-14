#include "polyMesh.H"
#include "particle.H"

// Helper function to create a simple box-shaped boundary
Bashyal::boundary createBoxBoundary(
    Foam::scalar lx, Foam::scalar ly, Foam::scalar lz
)
{
    Foam::pointField points(8);
    points[0] = Foam::point(-lx/2, -ly/2, -lz/2);
    points[1] = Foam::point( lx/2, -ly/2, -lz/2);
    points[2] = Foam::point( lx/2,  ly/2, -lz/2);
    points[3] = Foam::point(-lx/2,  ly/2, -lz/2);
    points[4] = Foam::point(-lx/2, -ly/2,  lz/2);
    points[5] = Foam::point( lx/2, -ly/2,  lz/2);
    points[6] = Foam::point( lx/2,  ly/2,  lz/2);
    points[7] = Foam::point(-lx/2,  ly/2,  lz/2);

    Foam::faceList faces(6);
    faces[0] = Foam::face({0, 3, 2, 1}); // bottom
    faces[1] = Foam::face({4, 5, 6, 7}); // top
    faces[2] = Foam::face({0, 4, 7, 3}); // left
    faces[3] = Foam::face({1, 2, 6, 5}); // right
    faces[4] = Foam::face({0, 1, 5, 4}); // front
    faces[5] = Foam::face({3, 7, 6, 2}); // back

    return Bashyal::boundary(points, faces);
}

// Function to add the particles to the simulation
void addInitialParticles(Bashyal::DEMSimulation& sim)
{
    Foam::Info << "Creating and adding initial particles..." << Foam::endl;

    // --- Particle 1: A static floor ---
    {
        Foam::scalar mass = 0.0; // 0 mass makes it static in Bullet
        Foam::point position(0, 0, 0);
        Bashyal::boundary box = createBoxBoundary(20, 20, 0.5);
        // MOI doesn't matter for a static body, but we provide identity
        sim.addParticle(Bashyal::particle(box, mass, Foam::tensor::I, position));
    }

    // --- Particle 2: A falling box ---
    {
        Foam::scalar mass = 2.0;
        Foam::point position(0.5, 0.5, 10.0);
        Bashyal::boundary box = createBoxBoundary(1, 1, 1);
        // Simple MOI for a box
        Foam::tensor moi = Foam::tensor::zero;
        moi.xx() = (1.0/12.0) * mass * (1*1 + 1*1);
        moi.yy() = (1.0/12.0) * mass * (1*1 + 1*1);
        moi.zz() = (1.0/12.0) * mass * (1*1 + 1*1);
        sim.addParticle(Bashyal::particle(box, mass, moi, position));
    }

    // --- Particle 3: Another falling box ---
    {
        Foam::scalar mass = 5.0; // Heavier
        Foam::point position(-0.5, -0.5, 15.0);
        Bashyal::boundary box = createBoxBoundary(1.5, 1.5, 1.5);
        Foam::tensor moi = Foam::tensor::zero;
        moi.xx() = (1.0/12.0) * mass * (1.5*1.5 + 1.5*1.5);
        moi.yy() = (1.0/12.0) * mass * (1.5*1.5 + 1.5*1.5);
        moi.zz() = (1.0/12.0) * mass * (1.5*1.5 + 1.5*1.5);
        sim.addParticle(Bashyal::particle(box, mass, moi, position));
    }
}