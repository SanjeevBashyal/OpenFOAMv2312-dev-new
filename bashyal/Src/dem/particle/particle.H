// --- START OF FILE particle.H ---

#ifndef particle_H
#define particle_H

#include "boundary.H"
#include "quaternion.H"
#include "vector.H"
#include "tensor.H"
#include "scalar.H"

namespace Bashyal
{
    /**
     * @struct forces
     * @brief Structure to store force and torque vectors together
     */
    struct forces
    {
        Foam::vector F = Foam::vector::zero;
        Foam::vector T = Foam::vector::zero;

        forces()
        {}

        forces(const forces& rhs)
        {
            F = rhs.F;
            T = rhs.T;
        }

        forces(const Foam::vector& F_in, const Foam::vector& T_in):
            F(F_in),
            T(T_in)
        {}

        void operator=(const forces& rhs)
        {
            F = rhs.F;
            T = rhs.T;
        }

        void operator+=(const forces& rhs)
        {
            F += rhs.F;
            T += rhs.T;
        }

        void operator-=(const forces& rhs)
        {
            F -= rhs.F;
            T -= rhs.T;
        }

        void clear()
        {
            F = Foam::vector::zero;
            T = Foam::vector::zero;
        }
    };

    /**
     * @class particle
     * @brief Represents a physical particle for Discrete Element Method (DEM) simulations.
     *
     * This class extends the geometric 'boundary' with physical properties
     * such as mass, velocity, and forces. It provides methods to update its
     * state over time based on applied forces and torques, forming the core
     * of a DEM particle.
     * 
     * Logic adapted from HFDIBDEM immersedBody class.
     */
    class particle : public boundary
    {
    private:
        //- Physical properties
        Foam::scalar mass_;                 // Mass of the particle
        Foam::tensor momentOfInertia_;      // Moment of inertia tensor in the body's local frame
        Foam::tensor invMomentOfInertia_;   // Inverse of the moment of inertia tensor

        //- Kinematic properties (State)
        Foam::point position_;              // Position of the center of mass in the world frame
        Foam::quaternion orientation_;      // Orientation of the particle (rotation from local to world frame)
        Foam::vector velocity_;             // Linear velocity of the center of mass
        Foam::vector angularVelocity_;      // Angular velocity in the world frame

        //- Old state variables (for time integration)
        Foam::vector velocityOld_;          // Previous time step velocity
        Foam::vector angularVelocityOld_;   // Previous time step angular velocity
        Foam::vector axis_;                 // Rotation axis (normalized)
        Foam::vector axisOld_;              // Previous rotation axis
        Foam::scalar omega_;                // Angular velocity magnitude
        Foam::scalar omegaOld_;             // Previous angular velocity magnitude

        //- Dynamic properties (Forces and Torques)
        Foam::vector force_;                // Total force accumulated on the particle
        Foam::vector torque_;               // Total torque accumulated on the particle
        
        //- Force and torque from coupling (fluid-solid interaction)
        forces FCoupling_;
        
        //- Force and torque from contact (particle-particle, particle-wall)
        forces FContact_;
        
        //- Acceleration
        Foam::vector acceleration_;         // Linear acceleration
        Foam::vector angularAcceleration_;  // Angular acceleration
        
        //- Total rotation matrix (cumulative rotation)
        Foam::tensor totRotMatrix_;
        
        //- Body operation mode
        //  0: static body
        //  1: prescribed translation only (freely rotating)
        //  2: prescribed rotation only (freely translating)
        //  3: prescribed translation and rotation
        //  4: prescribed translation with fixed axis rotation
        //  5: fully coupled (free movement)
        int bodyOperation_;
        
        //- Update torque flag
        bool updateTorque_;
        
        //- Courant number
        Foam::scalar CoNum_;
        Foam::scalar meanCoNum_;
        
        //- Characteristic diameter
        Foam::scalar dC_;
        
        //- Is particle active
        bool isActive_;
        
        //- Particle ID
        Foam::label id_;
        static Foam::label particleCount_;

    public:
        // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

        /**
         * @brief Default constructor.
         * Initializes a particle with default physical properties (e.g., mass=1, zero velocity).
         */
        particle();

        /**
         * @brief Constructor from a boundary object and physical properties.
         * @param b The boundary object defining the particle's geometry.
         * @param mass The mass of the particle.
         * @param moi The moment of inertia tensor of the particle.
         * @param initialPosition The initial position of the particle's center of mass.
         */
        particle(
            const boundary& b,
            const Foam::scalar mass,
            const Foam::tensor& moi,
            const Foam::point& initialPosition
        );

        // * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

        /**
         * @brief Destructor.
         */
        ~particle();


        // * * * * * * * * * * * * * * * * Methods * * * * * * * * * * * * * * * //

        /**
         * @brief Updates the particle's state over a given time step 'dt'.
         * Uses a semi-implicit Euler integration scheme to update position and velocity.
         * @param dt The time step for the integration.
         */
        void update(Foam::scalar dt);

        /**
         * @brief Updates movement based on acting forces (adapted from HFDIBDEM).
         * @param deltaT The time step for the integration.
         */
        void updateMovement(Foam::scalar deltaT);

        /**
         * @brief Updates movement with specified velocity and angular velocity.
         * @param Vel Linear velocity.
         * @param Axis Rotation axis.
         * @param omega Angular velocity magnitude.
         */
        void updateMovement(
            const Foam::vector& Vel,
            const Foam::vector& Axis,
            Foam::scalar omega
        );

        /**
         * @brief Component-wise movement update (adapted from HFDIBDEM).
         * @param deltaT Time step.
         * @param Vel Linear velocity.
         * @param Axis Rotation axis.
         * @param omega Angular velocity magnitude.
         */
        void updateMovementComp(
            Foam::scalar deltaT,
            const Foam::vector& Vel,
            const Foam::vector& Axis,
            Foam::scalar omega
        );

        /**
         * @brief Moves the particle according to body operation (adapted from HFDIBDEM).
         * @param deltaT Time step.
         */
        void moveParticle(Foam::scalar deltaT = -1.0);

        /**
         * @brief Updates coupling forces and torques (adapted from HFDIBDEM).
         * This method should be implemented to compute fluid-solid coupling forces.
         */
        void updateCoupling();

        /**
         * @brief Applies a force at a specific point on the particle.
         * Automatically calculates and adds the resulting torque.
         * @param force The force vector to apply.
         * @param applicationPoint The point in world coordinates where the force is applied.
         */
        void applyForce(const Foam::vector& force, const Foam::point& applicationPoint);

        /**
         * @brief Clears the accumulated forces and torques.
         * Should be called at the beginning of each time step before calculating new interactions.
         */
        void clearForceAndTorque();

        /**
         * @brief Resets contact forces.
         */
        void resetContactForces();

        /**
         * @brief Updates contact forces.
         * @param F Contact forces structure.
         */
        void updateContactForces(const forces& F);

        /**
         * @brief Computes body Courant number.
         * @param deltaT Time step.
         * @param charCellSize Characteristic cell size.
         */
        void computeBodyCoNumber(Foam::scalar deltaT, Foam::scalar charCellSize = 1e3);

        /**
         * @brief Returns the points of the geometry transformed into the world frame.
         * @return A pointField containing the vertices in their current world positions.
         */
        Foam::pointField getGlobalPoints() const;

        /**
         * @brief Support function for GJK algorithm.
         * Returns the vertex that is furthest in the given direction.
         * @param direction The search direction.
         * @return The support point in world coordinates.
         */
        Foam::point getSupport(const Foam::vector& direction) const;

        /**
         * @brief Writes the particle geometry to a VTP file using global coordinates.
         * @param filename The output VTP file name.
         */
        void writeVtp(const Foam::fileName& filename) const;

        /**
         * @brief Prints particle statistics (adapted from HFDIBDEM).
         */
        void printStats() const;

        /**
         * @brief Prints particle momentum (adapted from HFDIBDEM).
         */
        void printMomentum() const;

        /**
         * @brief Updates old movement variables.
         */
        void updateOldMovementVars();


        // * * * * * * * * * * * * * * Accessors (Getters) * * * * * * * * * * * * * //

        Foam::scalar mass() const { return mass_; }
        const Foam::point& position() const { return position_; }
        const Foam::quaternion& orientation() const { return orientation_; }
        const Foam::vector& velocity() const { return velocity_; }
        const Foam::vector& angularVelocity() const { return angularVelocity_; }
        const Foam::vector& force() const { return force_; }
        const Foam::vector& torque() const { return torque_; }
        const Foam::tensor& momentOfInertia() const { return momentOfInertia_; }
        
        const Foam::vector& getVel() const { return velocity_; }
        const Foam::scalar& getOmega() const { return omega_; }
        const Foam::vector& getAxis() const { return axis_; }
        const Foam::scalar& getCoNum() const { return CoNum_; }
        const Foam::scalar& getMeanCoNum() const { return meanCoNum_; }
        Foam::scalar getDC() const { return dC_; }
        const int& getBodyOperation() const { return bodyOperation_; }
        const bool& getIsActive() const { return isActive_; }
        const forces& getFCoupling() const { return FCoupling_; }
        const forces& getFContact() const { return FContact_; }
        const Foam::tensor& getTotRotMatrix() const { return totRotMatrix_; }
        
        Foam::label id() const { return id_; }


        // * * * * * * * * * * * * * * Modifiers (Setters) * * * * * * * * * * * * * * //

        void setPosition(const Foam::point& p) { position_ = p; }
        void setVelocity(const Foam::vector& v) { velocity_ = v; }
        void setAngularVelocity(const Foam::vector& w) { angularVelocity_ = w; }
        void setOrientation(const Foam::quaternion& q) { orientation_ = q; }
        void setBodyOperation(int op) { bodyOperation_ = op; }
        void setUpdateTorque(bool update) { updateTorque_ = update; }
        void setDC(Foam::scalar dc) { dC_ = dc; }
        void setIsActive(bool active) { isActive_ = active; }
        void setMass(Foam::scalar m) { mass_ = m; }
    };
}

#endif

// --- END OF FILE particle.H ---