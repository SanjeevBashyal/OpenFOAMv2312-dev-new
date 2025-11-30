// --- START OF FILE particle.C ---

#include "particle.H"
#include "IOstreams.H"
#include "error.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //
    
    Foam::label particle::particleCount_ = 0;

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    particle::particle()
    :   boundary(), // Call base class default constructor
        mass_(1.0),
        momentOfInertia_(Foam::tensor::I),
        invMomentOfInertia_(Foam::tensor::I),
        position_(Foam::point::zero),
        orientation_(Foam::quaternion::I),
        velocity_(Foam::vector::zero),
        angularVelocity_(Foam::vector::zero),
        velocityOld_(Foam::vector::zero),
        angularVelocityOld_(Foam::vector::zero),
        axis_(Foam::vector(0, 0, 1)),
        axisOld_(Foam::vector(0, 0, 1)),
        omega_(0.0),
        omegaOld_(0.0),
        force_(Foam::vector::zero),
        torque_(Foam::vector::zero),
        FCoupling_(),
        FContact_(),
        acceleration_(Foam::vector::zero),
        angularAcceleration_(Foam::vector::zero),
        totRotMatrix_(Foam::tensor::I),
        bodyOperation_(5), // Default: fully coupled
        updateTorque_(true),
        CoNum_(0.0),
        meanCoNum_(0.0),
        dC_(1.0),
        isActive_(true),
        id_(particleCount_++)
    {
    }

    particle::particle(
        const boundary& b,
        const Foam::scalar mass,
        const Foam::tensor& moi,
        const Foam::point& initialPosition
    )
    :   boundary(b), // Use boundary's copy constructor for geometry
        mass_(mass),
        momentOfInertia_(moi),
        invMomentOfInertia_(inv(moi)),
        position_(initialPosition),
        orientation_(Foam::quaternion::I), // Start with no rotation
        velocity_(Foam::vector::zero),
        angularVelocity_(Foam::vector::zero),
        velocityOld_(Foam::vector::zero),
        angularVelocityOld_(Foam::vector::zero),
        axis_(Foam::vector(0, 0, 1)),
        axisOld_(Foam::vector(0, 0, 1)),
        omega_(0.0),
        omegaOld_(0.0),
        force_(Foam::vector::zero),
        torque_(Foam::vector::zero),
        FCoupling_(),
        FContact_(),
        acceleration_(Foam::vector::zero),
        angularAcceleration_(Foam::vector::zero),
        totRotMatrix_(Foam::tensor::I),
        bodyOperation_(5), // Default: fully coupled
        updateTorque_(true),
        CoNum_(0.0),
        meanCoNum_(0.0),
        dC_(1.0), // Should be calculated from boundary
        isActive_(true),
        id_(particleCount_++)
    {
        // Calculate dC from bounding box
        Foam::boundBox bb = b.createBoundBox();
        dC_ = bb.mag();

        if (mass_ <= 0)
        {
            FatalErrorInFunction
                << "Particle mass must be positive."
                << abort(Foam::FatalError);
        }
        invMomentOfInertia_ = inv(momentOfInertia_);
    }

    particle::particle(
        const Foam::fileName& stlFile,
        const Foam::scalar mass,
        const Foam::tensor& moi,
        const Foam::point& initialPosition
    )
    :   boundary(stlFile), // Use boundary's STL constructor
        mass_(mass),
        momentOfInertia_(moi),
        invMomentOfInertia_(inv(moi)),
        position_(initialPosition),
        orientation_(Foam::quaternion::I), // Start with no rotation
        velocity_(Foam::vector::zero),
        angularVelocity_(Foam::vector::zero),
        velocityOld_(Foam::vector::zero),
        angularVelocityOld_(Foam::vector::zero),
        axis_(Foam::vector(0, 0, 1)),
        axisOld_(Foam::vector(0, 0, 1)),
        omega_(0.0),
        omegaOld_(0.0),
        force_(Foam::vector::zero),
        torque_(Foam::vector::zero),
        FCoupling_(),
        FContact_(),
        acceleration_(Foam::vector::zero),
        angularAcceleration_(Foam::vector::zero),
        totRotMatrix_(Foam::tensor::I),
        bodyOperation_(5), // Default: fully coupled
        updateTorque_(true),
        CoNum_(0.0),
        meanCoNum_(0.0),
        dC_(1.0), // Should be calculated from boundary
        isActive_(true),
        id_(particleCount_++)
    {
        // Calculate dC from bounding box
        Foam::boundBox bb = this->createBoundBox();
        dC_ = bb.mag();

        if (mass_ <= 0)
        {
            FatalErrorInFunction
                << "Particle mass must be positive."
                << abort(Foam::FatalError);
        }
        invMomentOfInertia_ = inv(momentOfInertia_);
    }


    // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

    particle::~particle()
    {
    }


    // * * * * * * * * * * * * * * * * Methods  * * * * * * * * * * * * * * * //

    void particle::update(Foam::scalar dt)
    {
        // Check if body is static
        if (bodyOperation_ == 0)
        {
            return;
        }

        // --- Translational Update (Semi-implicit Euler) ---
        // Total Force = External Force + Contact Force + Coupling Force
        Foam::vector totalForce = force_ + FContact_.F + FCoupling_.F;

        // F = m*a  =>  a = F/m
        Foam::vector acceleration = totalForce / mass_;
        velocity_ += acceleration * dt;
        position_ += velocity_ * dt;

        // --- Rotational Update ---
        // Total Torque = External Torque + Contact Torque + Coupling Torque
        Foam::vector totalTorque = torque_ + FContact_.T + FCoupling_.T;

        // T = I*alpha => alpha = I^-1 * T
        // Transform inverse MOI tensor from local to world frame
        Foam::tensor R = orientation_.R(); // Get rotation matrix from quaternion
        Foam::tensor invI_world = R & invMomentOfInertia_ & R.T();

        // Calculate angular acceleration
        Foam::vector angularAcceleration = invI_world & totalTorque;

        // Update angular velocity
        angularVelocity_ += angularAcceleration * dt;

        Foam::scalar omegaMag = Foam::mag(angularVelocity_);

        // Update orientation using the new angular velocity
        if (omegaMag > Foam::SMALL)
        {
            Foam::vector rotationAxis = angularVelocity_ / omegaMag; // Normalize to get the axis
            Foam::scalar rotationAngle = omegaMag * dt;         // angle = speed * time
            Foam::quaternion deltaQ(rotationAxis, rotationAngle); // Create quaternion from rotation vector
            orientation_ = deltaQ * orientation_;
            orientation_.normalise(); // Normalize to prevent numerical drift
        }
    }

    void particle::applyForce(const Foam::vector& force, const Foam::point& applicationPoint)
    {
        // Add to total force
        force_ += force;

        // Add to total torque: T = r x F
        // where r is the vector from the center of mass to the application point
        Foam::vector r = applicationPoint - position_;
        torque_ += (r ^ force);
    }

    void particle::clearForceAndTorque()
    {
        force_ = Foam::vector::zero;
        torque_ = Foam::vector::zero;
        FContact_.clear();
        FCoupling_.clear();
    }

    Foam::pointField particle::getGlobalPoints() const
    {
        // Get the local points from the base class
        const Foam::pointField& localPoints = this->vertices();
        Foam::pointField globalPoints(localPoints.size());

        // Create the transformation matrix from the orientation quaternion
        Foam::tensor R = orientation_.R();

        // Transform each local point to its global position
        forAll(localPoints, i)
        {
            // 1. Rotate the point from local to world orientation
            // 2. Translate the point by the particle's position
            globalPoints[i] = (R & localPoints[i]) + position_;
        }

        return globalPoints;
    }

    Foam::point particle::getSupport(const Foam::vector& direction) const
    {
        // Transform direction to local frame
        // d_local = R^T * d_world
        Foam::vector localDir = orientation_.R().T() & direction;
        
        // Find the vertex in the local frame that maximizes the dot product with localDir
        // The vertices are stored in the boundary (geomObject) base class
        const Foam::pointField& localPoints = vertices();
        
        Foam::scalar maxDot = -Foam::GREAT;
        Foam::label maxIdx = 0;
        
        forAll(localPoints, i)
        {
            Foam::scalar dot = localPoints[i] & localDir;
            if (dot > maxDot)
            {
                maxDot = dot;
                maxIdx = i;
            }
        }
        
        // Transform the support point back to world frame
        // p_world = position + R * p_local
        return position_ + (orientation_.R() & localPoints[maxIdx]);
    }

    void particle::writeVtp(const Foam::fileName& filename) const
    {
        Foam::OFstream vtpFile(filename);
        if (!vtpFile.good())
        {
            FatalErrorIn("writeVtp") << "Cannot open file " << filename << Foam::exit(Foam::FatalError);
        }

        // Get the global points (transformed to world coordinates)
        const Foam::pointField& globalVertices = getGlobalPoints();
        const Foam::faceList& faces = this->faces();

        // VTK XML Header
        vtpFile << "<?xml version=\"1.0\"?>" << Foam::endl;
        vtpFile << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">" << Foam::endl;
        vtpFile << "  <PolyData>" << Foam::endl;

        // Piece defines the geometry. We have one object, so one piece.
        vtpFile << "    <Piece NumberOfPoints=\"" << globalVertices.size()
                << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\""
                << faces.size() << "\">" << Foam::endl;

        // 1. Write Points (Vertices) - using global coordinates
        vtpFile << "      <Points>" << Foam::endl;
        vtpFile << "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << Foam::endl;
        for (const Foam::point& pt : globalVertices)
        {
            vtpFile << "          " << pt.x() << " " << pt.y() << " " << pt.z() << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Points>" << Foam::endl;

        // 2. Write Polygons (Faces)
        vtpFile << "      <Polys>" << Foam::endl;
        // a) Connectivity: a flat list of all vertex indices for all faces
        vtpFile << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << Foam::endl;
        for (const Foam::face& f : faces)
        {
            vtpFile << "          ";
            for (const Foam::label vIdx : f)
            {
                vtpFile << vIdx << " ";
            }
            vtpFile << Foam::endl;
        }
        vtpFile << "        </DataArray>" << Foam::endl;

        // b) Offsets: the cumulative count of vertices per face
        vtpFile << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << Foam::endl;
        vtpFile << "          ";
        Foam::label offset = 0;
        for (const Foam::face& f : faces)
        {
            offset += f.size();
            vtpFile << offset << " ";
        }
        vtpFile << Foam::endl;
        vtpFile << "        </DataArray>" << Foam::endl;
        vtpFile << "      </Polys>" << Foam::endl;

        // VTK XML Footer
        vtpFile << "    </Piece>" << Foam::endl;
        vtpFile << "  </PolyData>" << Foam::endl;
        vtpFile << "</VTKFile>" << Foam::endl;
    }

    //---------------------------------------------------------------------------//
    // Methods adapted from HFDIBDEM immersedBody class
    //---------------------------------------------------------------------------//

    void particle::updateMovement(Foam::scalar deltaT)
    {
        updateMovementComp(deltaT, velocity_, axis_, omega_);
    }

    void particle::updateMovement(
        const Foam::vector& Vel,
        const Foam::vector& Axis,
        Foam::scalar omega
    )
    {
        // Use a default time step - in HFDIBDEM this would come from mesh_.time().deltaT().value()
        Foam::scalar deltaT = 1e-5; // Default time step
        updateMovementComp(deltaT, Vel, Axis, omega);
    }

    void particle::updateMovementComp(
        Foam::scalar deltaT,
        const Foam::vector& Vel,
        const Foam::vector& Axis,
        Foam::scalar omega
    )
    {
        // Lambda function for translation update (adapted from HFDIBDEM)
        auto updateTranslation = [&]()
        {
            // Get gravity vector (if available)
            Foam::vector gVec = Foam::vector::zero;
            // Note: In pure DEM, gravity might be handled differently
            // This is adapted from HFDIBDEM where it checks meshObjects::gravity
            
            // Gravity force
            Foam::vector FG = mass_ * gVec;

            // Total force: coupling + contact + gravity
            Foam::vector F = FCoupling_.F;
            F += FContact_.F;
            F += FG;

            if (mass_ > 0)
            {
                // Compute current acceleration (assume constant over timeStep)
                acceleration_ = F / mass_;
                // Update body linear velocity
                velocity_ = Vel + deltaT * acceleration_;
            }
        };

        // Lambda function for rotation update (adapted from HFDIBDEM)
        auto updateRotation = [&]()
        {
            if (Foam::mag(momentOfInertia_) > 0)
            {
                Foam::vector T = FCoupling_.T;
                T += FContact_.T;

                // Transform moment of inertia to world frame
                Foam::tensor R = orientation_.R();
                Foam::tensor I_world = R & momentOfInertia_ & R.T();
                Foam::tensor invI_world = inv(I_world);

                // Update body angular acceleration
                angularAcceleration_ = invI_world & T;
                
                // Update body angular velocity
                Foam::vector Omega = Axis * omega + deltaT * angularAcceleration_;
                
                // Split Omega into axis_ and omega_
                omega_ = Foam::mag(Omega);

                if (omega_ < Foam::SMALL)
                {
                    axis_ = Foam::vector::one;
                }
                else
                {
                    axis_ = Omega / (omega_ + Foam::SMALL);
                }
                axis_ /= Foam::mag(axis_);
            }
        };

        // Lambda function for fixed axis rotation (adapted from HFDIBDEM)
        auto updateRotationFixedAxis = [&]()
        {
            Foam::vector T = FCoupling_.T;
            T += FContact_.T;

            // Transform moment of inertia to world frame
            Foam::tensor R = orientation_.R();
            Foam::tensor I_world = R & momentOfInertia_ & R.T();
            Foam::tensor invI_world = inv(I_world);

            // Update body angular velocity
            Foam::vector Omega = Axis * omega + deltaT * (invI_world & T);

            // Split Omega into axis_ and omega_
            omega_ = Foam::mag(Omega);

            Foam::vector newAxis = Omega / (omega_ + Foam::SMALL);
            if ((newAxis & axis_) < 0) 
            {
                axis_ *= (-1.0);
            }
        };

        // Handle different body operation modes (adapted from HFDIBDEM)
        if (bodyOperation_ == 0 || bodyOperation_ == 3)
        {
            return;
        }
        else if (bodyOperation_ == 1)
        {
            updateRotation();
            return;
        }
        else if (bodyOperation_ == 2)
        {
            updateTranslation();
            return;
        }
        else if (bodyOperation_ == 4)
        {
            updateRotationFixedAxis();
            return;
        }

        updateTranslation();
        if (updateTorque_)
        {
            updateRotation();
        }
    }

    void particle::moveParticle(Foam::scalar deltaT)
    {
        if (bodyOperation_ == 0) 
        {
            return;
        }

        if (Foam::mag(deltaT + 1.0) < Foam::SMALL)
        {
            // If deltaT is -1.0, use a default time step
            // In HFDIBDEM, this would use mesh_.time().deltaT().value()
            deltaT = 1e-5; // Default time step
        }

        // Incremental rotation angle
        Foam::scalar angle = omega_ * deltaT;

        // Translation increment
        Foam::vector transIncr = velocity_ * deltaT;

        // Rotation matrix (Rodrigues' rotation formula, adapted from HFDIBDEM)
        Foam::tensor rotMatrix = Foam::cos(angle) * Foam::tensor::I;
        rotMatrix += Foam::sin(angle) * Foam::tensor(
            0.0,      -axis_.z(),  axis_.y(),
            axis_.z(), 0.0,       -axis_.x(),
            -axis_.y(), axis_.x(),  0.0
        );
        rotMatrix += (1.0 - Foam::cos(angle)) * (axis_ * axis_);

        // Update total rotation matrix
        totRotMatrix_ = rotMatrix & totRotMatrix_;

        // Update orientation quaternion from rotation matrix
        // Extract quaternion from rotation matrix
        Foam::scalar trace = rotMatrix.xx() + rotMatrix.yy() + rotMatrix.zz();
        if (trace > 0)
        {
            Foam::scalar s = 0.5 / Foam::sqrt(trace + 1.0);
            Foam::scalar w = 0.25 / s;
            Foam::scalar x = (rotMatrix.zy() - rotMatrix.yz()) * s;
            Foam::scalar y = (rotMatrix.xz() - rotMatrix.zx()) * s;
            Foam::scalar z = (rotMatrix.yx() - rotMatrix.xy()) * s;
            Foam::quaternion deltaQ(w, Foam::vector(x, y, z));
            orientation_ = deltaQ * orientation_;
            orientation_.normalise();
        }

        // Move particle position
        position_ += transIncr;

        // Update angular velocity vector from axis and omega
        angularVelocity_ = axis_ * omega_;
    }

    void particle::updateCoupling()
    {
        // This method is adapted from HFDIBDEM's updateCoupling
        // In HFDIBDEM, it computes fluid-solid coupling forces from a volVectorField
        // For pure DEM, this can be a placeholder or used for external forces
        
        // Placeholder: In pure DEM, coupling forces might come from external sources
        // FCoupling_ should be set externally if needed
        // For now, we just ensure it's initialized
    }

    void particle::resetContactForces()
    {
        FContact_.clear();
    }

    void particle::updateContactForces(const forces& F)
    {
        FContact_ += F;
    }

    void particle::computeBodyCoNumber(Foam::scalar deltaT, Foam::scalar charCellSize)
    {
        // Adapted from HFDIBDEM's computeBodyCoNumber
        Foam::scalar VelMag = Foam::mag(velocity_);
        
        meanCoNum_ = 0.0;
        CoNum_ = 0.0;

        // Rotation body courant number
        Foam::scalar rotCoNumB = omega_ * dC_ * 0.5 * deltaT;

        // For a single particle, we use characteristic cell size
        if (charCellSize > Foam::SMALL)
        {
            Foam::scalar CoNumCell = VelMag * deltaT / charCellSize;
            CoNumCell += rotCoNumB / charCellSize;
            
            CoNum_ = CoNumCell;
            meanCoNum_ = CoNumCell;
        }
    }

    void particle::printStats() const
    {
        // Adapted from HFDIBDEM's printStats
        Foam::vector L = momentOfInertia_ & (axis_ * omega_);
        Foam::vector p = mass_ * velocity_;

        Foam::Info << "-- particle linear momentum: " << p
            << " magnitude: " << Foam::mag(p) << Foam::endl;
        Foam::Info << "-- particle angular momentum: " << L
            << " magnitude: " << Foam::mag(L) << Foam::endl;
        Foam::Info << "-- particle position: " << position_ << Foam::endl;
        Foam::Info << "-- particle linear velocity: " << velocity_
            << " magnitude: " << Foam::mag(velocity_) << Foam::endl;
        Foam::Info << "-- particle angular velocity: " << omega_
            << " magnitude: " << Foam::mag(omega_) << Foam::endl;
        Foam::Info << "-- particle rotation axis: " << axis_
            << " magnitude: " << Foam::mag(axis_) << Foam::endl;
    }

    void particle::printMomentum() const
    {
        // Adapted from HFDIBDEM's printMomentum
        Foam::vector L = momentOfInertia_ & (axis_ * omega_);
        Foam::vector p = mass_ * velocity_;

        Foam::Info << "-- particle linear momentum: " << p
            << " magnitude: " << Foam::mag(p) << Foam::endl;
        Foam::Info << "-- particle angular momentum: " << L
            << " magnitude: " << Foam::mag(L) << Foam::endl;
    }

    void particle::updateOldMovementVars()
    {
        velocityOld_ = velocity_;
        axisOld_ = axis_;
        omegaOld_ = omega_;
        angularVelocityOld_ = angularVelocity_;
    }

} // End namespace Bashyal

// --- END OF FILE particle.C ---