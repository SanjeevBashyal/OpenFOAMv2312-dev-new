#include "particle.H"
#include "wall.H"
#include "gjkContact.H"
#include "contactForce.H"
#include "timeRegistry.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include <fstream>

using namespace Foam;
using namespace Bashyal;

// Helper to create a cube particle
particle createCubeParticle(point position, scalar L, scalar m)
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
    return p;
}

// Helper function to write wall VTP
void writeWallVtp(const fileName& filename, const wall& w, scalar size = 2.0)
{
    OFstream vtpFile(filename);
    if (!vtpFile.good())
    {
        FatalErrorIn("writeWallVtp") << "Cannot open file " << filename << exit(FatalError);
    }

    // We need to generate points for the plane based on its position and normal.
    // Let's create a local coordinate system on the plane.
    vector n = w.normal();
    n /= mag(n);
    
    // Find two tangent vectors
    vector t1 = vector(1, 0, 0);
    if (mag(t1 ^ n) < SMALL) t1 = vector(0, 1, 0);
    
    vector u = (n ^ t1); // First tangent
    u /= mag(u);
    
    vector v = (n ^ u); // Second tangent
    v /= mag(v);
    
    point center = w.position();
    scalar s = size / 2.0;
    
    // 4 corners
    point p0 = center - s*u - s*v;
    point p1 = center + s*u - s*v;
    point p2 = center + s*u + s*v;
    point p3 = center - s*u + s*v;

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
         << "Test Case 2: Sliding Block on Inclined Plane\n"
         << "========================================\n" << Foam::endl;

    // 1. Setup Parameters
    scalar dt = 1e-4;           // Time step (s)
    scalar endTime = 1.0;       // Duration
    vector gVec(0, 0, -9.81);   // Gravity vector
    scalar g = mag(gVec);
    
    // Geometry
    scalar thetaDeg = 30.0;
    scalar theta = thetaDeg * constant::mathematical::pi / 180.0;
    
    // Wall Normal: Tilted by theta.
    // If flat is (0,0,1), tilted 30 deg around Y axis (raising X) would be:
    // n = (-sin(theta), 0, cos(theta))
    vector wallNormal(-Foam::sin(theta), 0, Foam::cos(theta));
    wallNormal /= mag(wallNormal);
    
    // Particle Properties
    scalar L = 0.1;             // Side length (m)
    scalar m = 1.0;             // Mass (kg)
    
    // Contact Parameters
    scalar k1 = 1e4;            // Normal stiffness (softer for stability)
    scalar k2 = 1e4;            // Unloading stiffness
    scalar kt = 1e4;            // Tangential stiffness
    scalar mu = 0.3;            // Friction coefficient
    scalar gammaN = 200.0;      // Critical damping for k=1e4
    scalar gammaT = 0.0;
    scalar crustWidth = 1e-4;

    // 2. Create Objects
    mkDir("VTK_Sliding");
    
    // Wall (Inclined Plane passing through origin)
    wall inclinedPlane(point(0,0,0), wallNormal, "inclinedPlane");
    writeWallVtp("VTK_Sliding/wall.vtp", inclinedPlane, 5.0);

    // Particle
    // Place slightly inside to engage spring immediately (equilibrium is ~8.5e-4 m)
    point initialPos = point(0,0,0) + (0.5 * L - 1.0e-3) * wallNormal; 
    particle p = createCubeParticle(initialPos, L, m);
    
    // Rotate particle
    vector rotAxis(0, 1, 0);
    quaternion q(rotAxis, theta); 
    p.setOrientation(q);
    
    // Contact System
    gjkContact contactDetector;
    contactForce forceComputer(k1, k2, kt, gammaN, gammaT, mu);
    
    // Time Registry
    timeRegistry timeReg(dt, endTime, "VTK_Sliding");
    timeReg.setWriteInterval(0.05);
    timeReg.addParticle(p, "block");
    timeReg.addStaticObject("inclinedPlane", "VTK_Sliding/wall.vtp");

    // Data Logging
    std::ofstream dataFile("sliding_data.csv");
    dataFile << "Time,VelMag,TangentialForce,NormalForce,PosX,PosY,PosZ" << std::endl;

    Info << "Starting simulation..." << Foam::endl;

    // 3. Time Loop
    scalar t = 0.0;
    label step = 0;
    
    while (t < endTime)
    {
        // A. Clear Forces
        p.clearForceAndTorque();
        
        // B. Apply Gravity
        vector Fg = m * gVec;
        p.applyForce(Fg, p.position());
        
        // C. Detect Contact
        ContactInfo contact = contactDetector.detectContactParticleWall(p, inclinedPlane, crustWidth);
        
        scalar FnMag = 0.0;
        scalar FtMag = 0.0;
        
        List<string> activeKeys;

        if (contact.isContact)
        {
            // D. Compute Contact Force
            forceComputer.computeForceParticleWall(p, inclinedPlane, contact, dt);
            
            activeKeys.append("W_" + inclinedPlane.name() + "_" + name(p.id()));
            
            // Analyze Forces for Validation
            // Total contact force applied to particle
            vector F_contact = p.force() - Fg;
            
            // Decompose into Normal and Tangential components relative to Wall Normal
            vector n = wallNormal; // Wall normal
            vector FnVec = (F_contact & n) * n;
            vector FtVec = F_contact - FnVec;
            
            FnMag = mag(FnVec);
            FtMag = mag(FtVec);
        }
        
        forceComputer.cleanHistory(activeKeys);
        
        // E. Update Particle
        p.update(dt);
        
        // F. Advance Time
        t += dt;
        step++;
        
        timeReg.setCurrentTime(t);
        timeReg.advanceTime();
        
        // Log Data
        if (step % 100 == 0)
        {
            scalar vMag = mag(p.velocity());
            point pos = p.position();
            dataFile << t << "," << vMag << "," << FtMag << "," << FnMag << "," << pos.x() << "," << pos.y() << "," << pos.z() << std::endl;
        }
    }

    timeReg.writeTimeSeries();
    Info << "Simulation completed. Data written to sliding_data.csv and VTK_Sliding/" << Foam::endl;

    return 0;
}
