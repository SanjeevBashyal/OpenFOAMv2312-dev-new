#include "particle.H"
#include "spaceRegistry.H"
#include "timeRegistry.H"
#include "OFstream.H"
#include "quaternion.H" 
#include "point.H"
#include "vector.H"
#include "faceList.H"
#include "pointField.H"
#include "constants.H"
#include <cmath>

using namespace Bashyal;
using namespace Foam;

int main(int argc, char *argv[])
{
    // Simulation parameters
    scalar endTime = 2.0;
    scalar dt = 0.1;
    
    // Create particles with different initial positions and velocities
    vector linearVelocity1(0.5, 0.0, 0.0);  // Moving in +x direction
    vector linearVelocity2(-0.3, 0.0, 0.0); // Moving in -x direction
    vector angularVelocity(0.0, 0.0, 1.57); // Rotating around Z axis
    
    point initialPosition1(0, 0, 0);
    point initialPosition2(2, 0, 0);  // Start 2 units away in x direction

    // Create a time registry
    timeRegistry registry(dt, endTime, "VTK_spaceRegistry");
    registry.setWriteInterval(0.2);
    registry.setWriteAtStart(true);
    registry.setWriteAtEnd(true);

    // Create particles (cubes)
    pointField vertices(8);
    vertices[0] = point(-0.5, -0.5, -0.5);
    vertices[1] = point(0.5, -0.5, -0.5);
    vertices[2] = point(0.5, 0.5, -0.5);
    vertices[3] = point(-0.5, 0.5, -0.5);
    vertices[4] = point(-0.5, -0.5, 0.5);
    vertices[5] = point(0.5, -0.5, 0.5);
    vertices[6] = point(0.5, 0.5, 0.5);
    vertices[7] = point(-0.5, 0.5, 0.5);

    faceList faces(6);
    faces[0] = face({0, 3, 2, 1});
    faces[1] = face({4, 5, 6, 7});
    faces[2] = face({0, 1, 5, 4});
    faces[3] = face({2, 3, 7, 6});
    faces[4] = face({1, 2, 6, 5});
    faces[5] = face({3, 0, 4, 7});

    // Create boundary objects and particles
    boundary cubeBoundary1(vertices, faces);
    boundary cubeBoundary2(vertices, faces);
    
    particle p1(cubeBoundary1, 1.0, tensor::I, initialPosition1);
    particle p2(cubeBoundary2, 1.0, tensor::I, initialPosition2);

    // Set velocities
    p1.setVelocity(linearVelocity1);
    p1.setAngularVelocity(angularVelocity);
    p2.setVelocity(linearVelocity2);
    p2.setAngularVelocity(angularVelocity);

    // Add particles to time registry
    registry.addParticle(p1, "particle1");
    registry.addParticle(p2, "particle2");

    // Create space registry for collision detection
    boundBox domain(point(-5, -5, -5), point(5, 5, 5));
    vector cellSize(1.0, 1.0, 1.0);
    spaceRegistry spaceReg(domain, cellSize, 1.2);

    // Add particles to space registry
    spaceReg.addParticle(&p1, "particle1");
    spaceReg.addParticle(&p2, "particle2");

    Info << "Starting DEM simulation with collision detection" << endl;
    Info << "Particle 1: position=" << p1.position() << ", velocity=" << p1.velocity() << endl;
    Info << "Particle 2: position=" << p2.position() << ", velocity=" << p2.velocity() << endl;
    
    for (scalar t = 0; t <= endTime + SMALL; t += dt)
    {
        // Update space registry with current particle states
        spaceReg.updateParticles(dt);
        
        // Find potential collision pairs
        auto collisionPairs = spaceReg.findCollisionPairs();
        
        if (!collisionPairs.empty())
        {
            Info << "Time: " << t << "s - Potential collisions detected:" << endl;
            for (const auto& pair : collisionPairs)
            {
                Info << "  " << pair.first << " <-> " << pair.second << endl;
            }
        }
        else
        {
            Info << "Time: " << t << "s - No collisions detected" << endl;
        }
        
        // Advance time (this will update particles and write VTP files)
        registry.advanceTime();
        
        Info << "  Particle 1 position: " << p1.position() << endl;
        Info << "  Particle 2 position: " << p2.position() << endl;
    }
    
    // Write final time series
    registry.writeTimeSeries();
    Info << "Simulation finished. Output written to VTK_spaceRegistry/timeSeries.pvd" << endl;
    return 0;
} 