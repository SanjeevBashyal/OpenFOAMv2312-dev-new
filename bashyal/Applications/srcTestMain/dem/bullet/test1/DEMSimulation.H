#ifndef DEMSimulation_H
#define DEMSimulation_H

#include "particle.H"
#include <vector>
#include <memory>

// Bullet Physics Headers
#include <btBulletDynamicsCommon.h>

namespace Bashyal
{
    /**
     * @class DEMSimulation
     * @brief Manages a multi-particle simulation using the Bullet Physics library.
     *
     * This class orchestrates the entire simulation. It initializes the
     * Bullet physics world, adds particles (both as Bashyal::particle and
     * btRigidBody), runs the simulation loop, and synchronizes the state
     * back to the Bashyal::particle objects for analysis and output.
     */
    class DEMSimulation
    {
    private:
        // --- Bullet Physics World Components ---
        std::unique_ptr<btBroadphaseInterface> broadphase_;
        std::unique_ptr<btDefaultCollisionConfiguration> collisionConfiguration_;
        std::unique_ptr<btCollisionDispatcher> dispatcher_;
        std::unique_ptr<btSequentialImpulseConstraintSolver> solver_;
        std::unique_ptr<btDiscreteDynamicsWorld> dynamicsWorld_;

        // --- Particle Management ---
        // Store our custom particles
        Foam::List<particle> particles_;
        // Store Bullet's physics objects. Using std::vector for Bullet's non-Foam types.
        std::vector<std::unique_ptr<btRigidBody>> rigidBodies_;
        std::vector<std::unique_ptr<btCollisionShape>> collisionShapes_;
        std::vector<std::unique_ptr<btDefaultMotionState>> motionStates_;

    public:
        // * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * //

        /**
         * @brief Constructor. Initializes the Bullet physics world.
         * @param gravity The gravity vector for the simulation.
         */
        DEMSimulation(const Foam::vector& gravity);

        // * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * //

        /**
         * @brief Destructor. Cleans up all Bullet objects.
         */
        ~DEMSimulation();

        // * * * * * * * * * * * * * * * * Methods * * * * * * * * * * * * * * * //

        /**
         * @brief Adds a particle to the simulation.
         * Creates a corresponding btRigidBody and adds both to the simulation.
         * @param p The particle to add. It will be moved into the simulation.
         */
        void addParticle(particle&& p);

        /**
         * @brief Runs the main simulation loop until the specified end time.
         * @param endTime The total simulation time to run for.
         * @param dt The time step for each physics update.
         * @param writeInterval The time interval at which to write output files.
         */
        void run(Foam::scalar endTime, Foam::scalar dt, Foam::scalar writeInterval);

    private:
        /**
         * @brief Executes a single step of the physics simulation.
         * @param dt The time step for this update.
         */
        void step(Foam::scalar dt);

        /**
         * @brief Synchronizes the state from Bullet back to our particles.
         * After Bullet calculates the new positions and orientations, this
         * method updates the corresponding Bashyal::particle objects.
         */
        void synchronizeFromBullet();

        /**
         * @brief Writes the current state of all particles to VTP files.
         * @param time The current simulation time, used for file naming.
         */
        void writeOutput(Foam::scalar time) const;
    };
} // End namespace Bashyal

#endif