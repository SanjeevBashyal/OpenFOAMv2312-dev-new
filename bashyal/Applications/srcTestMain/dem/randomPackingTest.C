/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Test Case 5: Random Loose Packing (Void Fraction Test)
    Validates multi-body contacts, geometric exclusion, and energy dissipation.
    
\*---------------------------------------------------------------------------*/

#include "particle.H"
#include "wall.H"
#include "gjkContact.H"
#include "contactForce.H"
#include "timeRegistry.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "Random.H"
#include <fstream>
#include <vector>

using namespace Foam;
using namespace Bashyal;

// Helper to create a cube particle
particle createCubeParticle(point position, scalar L, scalar m, const quaternion& q)
{
    // Define vertices of a cube centered at position
    pointField vertices(8);
    scalar s = L / 2.0;
    
    vertices[0] = point(-s, -s, -s);
    vertices[1] = point( s, -s, -s);
    vertices[2] = point( s,  s, -s);
    vertices[3] = point(-s,  s, -s);
    vertices[4] = point(-s, -s,  s);
    vertices[5] = point( s, -s,  s);
    vertices[6] = point( s,  s,  s);
    vertices[7] = point(-s,  s,  s);

    // Define faces (standard cube connectivity)
    faceList faces(6);
    faces[0] = face{0, 4, 7, 3}; // Left
    faces[1] = face{1, 2, 6, 5}; // Right
    faces[2] = face{0, 1, 5, 4}; // Bottom
    faces[3] = face{3, 7, 6, 2}; // Top
    faces[4] = face{0, 3, 2, 1}; // Back
    faces[5] = face{4, 5, 6, 7}; // Front

    // Create boundary object
    boundary b(vertices, faces);
    
    // Calculate Inertia for Cube: I = m*L^2/6
    scalar Ival = m * L * L / 6.0;
    tensor I = tensor::I * Ival;
    
    // Create particle
    particle p(b, m, I, position);
    p.setOrientation(q);
    
    return p;
}

// Helper function to write wall VTP
void writeWallVtp(const fileName& filename, const wall& w, scalar sizeX, scalar sizeY)
{
    OFstream vtpFile(filename);
    if (!vtpFile.good())
    {
        FatalErrorIn("writeWallVtp") << "Cannot open file " << filename << exit(FatalError);
    }

    // Assume wall is flat on Z=0 for floor, or vertical for sides.
    // But here we construct specific geometry.
    // Let's just write a quad based on position and normal?
    // For the box, we'll manually construct the VTP for the whole box in main or separate helper.
    // This function is generic for a single plane.
    
    vector n = w.normal();
    n /= mag(n);
    
    // Arbitrary tangents
    vector t1(1, 0, 0);
    if (mag(t1 ^ n) < SMALL) t1 = vector(0, 1, 0);
    vector u = (n ^ t1); u /= mag(u);
    vector v = (n ^ u); v /= mag(v);
    
    point center = w.position();
    scalar sx = sizeX / 2.0;
    scalar sy = sizeY / 2.0;
    
    point p0 = center - sx*u - sy*v;
    point p1 = center + sx*u - sy*v;
    point p2 = center + sx*u + sy*v;
    point p3 = center - sx*u + sy*v;

    vtpFile << "<?xml version=\"1.0\"?>" << Foam::endl;
    vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << Foam::endl;
    vtpFile << "  <PolyData>" << Foam::endl;
    vtpFile << "    <Piece NumberOfPoints=\"4\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"1\">" << Foam::endl;
    
    vtpFile << "      <Points>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          " << p0.x() << " " << p0.y() << " " << p0.z() << Foam::endl;
    vtpFile << "          " << p1.x() << " " << p1.y() << " " << p1.z() << Foam::endl;
    vtpFile << "          " << p2.x() << " " << p2.y() << " " << p2.z() << Foam::endl;
    vtpFile << "          " << p3.x() << " " << p3.y() << " " << p3.z() << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "      </Points>" << Foam::endl;

    vtpFile << "      <Polys>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          0 1 2 3" << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << Foam::endl;
    vtpFile << "          4" << Foam::endl;
    vtpFile << "        </DataArray>" << Foam::endl;
    vtpFile << "      </Polys>" << Foam::endl;

    vtpFile << "    </Piece>" << Foam::endl;
    vtpFile << "  </PolyData>" << Foam::endl;
    vtpFile << "</VTKFile>" << Foam::endl;
}

int main(int argc, char *argv[])
{
    Info << "\n========================================\n"
         << "Test Case 5: Random Loose Packing\n"
         << "========================================\n" << Foam::endl;

    // 1. Parameters
    scalar dt = 1e-4;
    scalar endTime = 5.0;       // Allow enough time to settle
    vector gVec(0, 0, -9.81);
    
    // Particles
    label nParticles = 200;     // User requested 200
    scalar L = 0.05;            // 5cm cubes
    scalar m = 1.0;             // 1kg
    
    // Container (10L x 10L base)
    scalar boxW = 10.0 * L;     // 0.5m
    scalar boxH = 20.0 * L;     // Tall enough
    
    // Contact
    scalar k = 1e5;             // Stiffness
    scalar mu = 0.5;            // Friction (Crucial for packing)
    scalar gammaN = 50.0;       // Damping
    scalar crustWidth = 1e-4;
    
    // 2. Setup Objects
    mkDir("VTK_Packing");
    
    // Walls: Floor + 4 Sides
    // Floor at z=0
    wall floor(point(0,0,0), vector(0,0,1), "floor");
    writeWallVtp("VTK_Packing/floor.vtp", floor, boxW, boxW);
    
    // Side Walls (Centered at 0,0)
    // Left (x = -boxW/2) -> Normal (1,0,0)
    wall wallL(point(-boxW/2, 0, boxH/2), vector(1,0,0), "wallL");
    writeWallVtp("VTK_Packing/wallL.vtp", wallL, boxH, boxW); // Height is Z size, Width is Y size

    // Right (x = boxW/2) -> Normal (-1,0,0)
    wall wallR(point(boxW/2, 0, boxH/2), vector(-1,0,0), "wallR");
    writeWallVtp("VTK_Packing/wallR.vtp", wallR, boxH, boxW);

    // Back (y = -boxW/2) -> Normal (0,1,0)
    wall wallB(point(0, -boxW/2, boxH/2), vector(0,1,0), "wallB");
    writeWallVtp("VTK_Packing/wallB.vtp", wallB, boxW, boxH); // Width is X size

    // Front (y = boxW/2) -> Normal (0,-1,0)
    wall wallF(point(0, boxW/2, boxH/2), vector(0,-1,0), "wallF");
    writeWallVtp("VTK_Packing/wallF.vtp", wallF, boxW, boxH);
    
    List<wall> walls;
    walls.append(floor);
    walls.append(wallL);
    walls.append(wallR);
    walls.append(wallB);
    walls.append(wallF);
    
    // Generate Particles
    std::vector<particle> particles;
    particles.reserve(nParticles); // Reserve memory to avoid pointer invalidation
    Random rnd(1234);
    
    // Pouring parameters
    scalar emitInterval = 0.01; // Emit every 0.01s
    scalar lastEmitTime = -emitInterval; // Ensure first one emits at t=0
    scalar emitHeight = boxH - 2.0*L; // Emit near top
    
    // 3. Simulation Components
    gjkContact contactDetector;
    contactForce forceComputer(k, k, k, gammaN, 0.0, mu);
    
    timeRegistry timeReg(dt, endTime, "VTK_Packing");
    timeReg.setWriteInterval(0.05);
    
    // Register walls
    timeReg.addStaticObject("floor", "VTK_Packing/floor.vtp");
    timeReg.addStaticObject("wallL", "VTK_Packing/wallL.vtp");
    timeReg.addStaticObject("wallR", "VTK_Packing/wallR.vtp");
    timeReg.addStaticObject("wallB", "VTK_Packing/wallB.vtp");
    timeReg.addStaticObject("wallF", "VTK_Packing/wallF.vtp");
    
    std::ofstream dataFile("packing_data.csv");
    dataFile << "Time,KE,CoordNum,MaxZ,SolidFraction" << std::endl;
    
    // 4. Time Loop
    scalar t = 0.0;
    label step = 0;
    
    Info << "Starting simulation with pouring..." << Foam::endl;
    
    while (t < endTime)
    {
        // --- Emission Logic ---
        if (particles.size() < nParticles && (t - lastEmitTime) >= emitInterval)
        {
            // Try to emit a particle
            // Random position at emitHeight
            scalar rX = rnd.sample01<scalar>() * (boxW - 2.0*L) - (boxW/2.0 - L);
            scalar rY = rnd.sample01<scalar>() * (boxW - 2.0*L) - (boxW/2.0 - L);
            point pos(rX, rY, emitHeight);
            
            // Check overlap with existing particles
            bool overlap = false;
            for (const auto& p : particles)
            {
                if (mag(p.position() - pos) < L * 1.1) // 1.1 safety factor
                {
                    overlap = true;
                    break;
                }
            }
            
            if (!overlap)
            {
                // Random Rotation
                vector axis(rnd.sample01<scalar>(), rnd.sample01<scalar>(), rnd.sample01<scalar>());
                axis /= (mag(axis) + SMALL);
                scalar angle = rnd.sample01<scalar>() * constant::mathematical::twoPi;
                quaternion q(axis, angle);
                
                particle newP = createCubeParticle(pos, L, m, q);
                particles.push_back(newP);
                
                // Register with timeReg
                timeReg.addParticle(particles.back(), "cube_" + std::to_string(particles.size()-1));
                
                lastEmitTime = t;
                // Info << "Emitted particle " << particles.size() << " at t=" << t << Foam::endl;
            }
        }

        // A. Clear Forces
        for (auto& p : particles) p.clearForceAndTorque();
        
        // B. Apply Gravity
        for (auto& p : particles) p.applyForce(m * gVec, p.position());
        
        // C. Contact Detection & Force
        label totalContacts = 0;
        List<string> activeKeys; 
        
        // Particle-Wall
        for (size_t i=0; i<particles.size(); ++i)
        {
            for (const auto& w : walls)
            {
                ContactInfo c = contactDetector.detectContactParticleWall(particles[i], w, crustWidth);
                if (c.isContact)
                {
                    forceComputer.computeForceParticleWall(particles[i], w, c, dt);
                    activeKeys.append("W_" + w.name() + "_" + name(particles[i].id()));
                }
            }
        }
        
        // Particle-Particle
        for (size_t i=0; i<particles.size(); ++i)
        {
            for (size_t j=i+1; j<particles.size(); ++j)
            {
                if (mag(particles[i].position() - particles[j].position()) > 2.0*L) continue;
                
                ContactInfo c = contactDetector.detectContact(particles[i], particles[j], crustWidth, crustWidth);
                if (c.isContact)
                {
                    forceComputer.computeForce(particles[i], particles[j], c, dt);
                    label idA = particles[i].id();
                    label idB = particles[j].id();
                    std::string key = "P_" + name(min(idA, idB)) + "_" + name(max(idA, idB));
                    activeKeys.append(key);
                    totalContacts++;
                }
            }
        }
        
        forceComputer.cleanHistory(activeKeys);
        
        // D. Update
        scalar totalKE = 0.0;
        scalar maxZ = -GREAT;
        
        for (auto& p : particles)
        {
            p.update(dt);
            totalKE += 0.5 * m * magSqr(p.velocity()); 
            maxZ = max(maxZ, p.position().z() + L/2.0);
        }
        
        // E. Analysis
        if (step % 100 == 0)
        {
            scalar coordNum = (particles.size() > 0) ? 2.0 * scalar(totalContacts) / scalar(particles.size()) : 0.0;
            
            scalar V_solids = particles.size() * L * L * L;
            scalar V_total = boxW * boxW * maxZ;
            scalar phi = (V_total > SMALL) ? V_solids / V_total : 0.0;
            
            dataFile << t << "," << totalKE << "," << coordNum << "," << maxZ << "," << phi << std::endl;
            
            if (step % 1000 == 0)
            {
                Info << "Time=" << t << " N=" << particles.size() << " KE=" << totalKE << " Phi=" << phi << " Zmax=" << maxZ << Foam::endl;
            }
            
            // Early exit if settled?
            // Only if all particles emitted and settled
            if (particles.size() == nParticles && t > (lastEmitTime + 1.0) && totalKE < 0.5)
            {
                Info << "Settled at t=" << t << ". Stopping." << Foam::endl;
                break;
            }
        }
        
        // F. Advance Time
        t += dt;
        step++;
        
        timeReg.setCurrentTime(t);
        timeReg.advanceTime();
    }
    
    timeReg.writeTimeSeries();
    Info << "Simulation completed." << Foam::endl;
    
    return 0;
}
