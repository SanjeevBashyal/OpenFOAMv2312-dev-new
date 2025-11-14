#ifndef BULLET_INTEGRATOR_H
#define BULLET_INTEGRATOR_H

#include "particle.H"
#include "timeRegistry.H"
#include "foamBulletConverter.H"

#include <memory>
#include <vector>

// Bullet Physics
#include <btBulletDynamicsCommon.h>

namespace Bashyal {

class bulletIntegrator
{
private:
    // Bullet world components
    std::unique_ptr<btBroadphaseInterface> broadphase_;
    std::unique_ptr<btDefaultCollisionConfiguration> collisionConfig_;
    std::unique_ptr<btCollisionDispatcher> dispatcher_;
    std::unique_ptr<btSequentialImpulseConstraintSolver> solver_;
    std::unique_ptr<btDiscreteDynamicsWorld> world_;

    // Mapped particles and Bullet objects
    Foam::List<particle*> particles_;                 // pointers, not owning
    std::vector<std::unique_ptr<btRigidBody>> bodies_;
    std::vector<std::unique_ptr<btCollisionShape>> shapes_;
    std::vector<std::unique_ptr<btDefaultMotionState>> motionStates_;

public:
    // Construct with gravity
    explicit bulletIntegrator(const Foam::vector& gravity);

    // Non-copyable
    bulletIntegrator(const bulletIntegrator&) = delete;
    bulletIntegrator& operator=(const bulletIntegrator&) = delete;

    // Movable
    bulletIntegrator(bulletIntegrator&&) = default;
    bulletIntegrator& operator=(bulletIntegrator&&) = default;

    ~bulletIntegrator();

    // Register all particles from a timeRegistry into the Bullet world
    void registerFrom(timeRegistry& tr, Foam::scalar restitution = 0.5, Foam::scalar friction = 0.5);

    // Advance physics by dt and synchronize back to particles
    void step(Foam::scalar dt, int maxSubSteps = 10);

    // Remove all bodies and reset the world (keeps gravity/settings)
    void clear();

    // Add a static box collider (mass=0), unaffected by gravity
    void addStaticBox(const Foam::point& center,
                        const Foam::vector& halfExtents,
                        const Foam::quaternion& orientation,
                        Foam::scalar restitution = 0.2,
                        Foam::scalar friction = 0.8);
};

} // namespace Bashyal

#endif // BULLET_INTEGRATOR_H
