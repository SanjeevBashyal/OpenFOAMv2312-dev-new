#include "bulletIntegrator.H"

namespace Bashyal {

bulletIntegrator::bulletIntegrator(const Foam::vector& gravity)
{
    collisionConfig_ = std::make_unique<btDefaultCollisionConfiguration>();
    dispatcher_ = std::make_unique<btCollisionDispatcher>(collisionConfig_.get());
    broadphase_ = std::make_unique<btDbvtBroadphase>();
    solver_ = std::make_unique<btSequentialImpulseConstraintSolver>();

    world_ = std::make_unique<btDiscreteDynamicsWorld>(
        dispatcher_.get(), broadphase_.get(), solver_.get(), collisionConfig_.get());

    world_->setGravity(FoamBulletConverter::toBullet(gravity));
}

bulletIntegrator::~bulletIntegrator()
{
    clear();
}

void bulletIntegrator::clear()
{
    if (!world_) return;

    for (int i = world_->getNumCollisionObjects() - 1; i >= 0; --i)
    {
        btCollisionObject* obj = world_->getCollisionObjectArray()[i];
        btRigidBody* body = btRigidBody::upcast(obj);
        if (body)
        {
            world_->removeRigidBody(body);
        }
        else
        {
            world_->removeCollisionObject(obj);
        }
    }

    // Clear containers (unique_ptr handle deletion)
    bodies_.clear();
    shapes_.clear();
    motionStates_.clear();
    particles_.clear();
}

void bulletIntegrator::registerFrom(timeRegistry& tr, Foam::scalar restitution, Foam::scalar friction)
{
    // Clear any previous registration
    clear();

    for (const auto& kv : tr.particles())
    {
        particle* pp = kv.second;
        if (!pp) continue;

        // Build convex hull from particle's local vertices
        std::unique_ptr<btConvexHullShape> shape(new btConvexHullShape());
        FoamBulletConverter::appendConvexHullPoints(*shape, pp->vertices());
        shape->recalcLocalAabb();
        shape->optimizeConvexHull();

        // Initial transform
        btTransform startTransform = FoamBulletConverter::toBulletTransform(pp->position(), pp->orientation());
        std::unique_ptr<btDefaultMotionState> motionState(new btDefaultMotionState(startTransform));

        // Inertia
        btVector3 localInertia(0, 0, 0);
        if (pp->mass() > 0)
        {
            shape->calculateLocalInertia(pp->mass(), localInertia);
        }

        btRigidBody::btRigidBodyConstructionInfo rbInfo(pp->mass(), motionState.get(), shape.get(), localInertia);
        rbInfo.m_restitution = restitution;
        rbInfo.m_friction = friction;

        std::unique_ptr<btRigidBody> body(new btRigidBody(rbInfo));

        // Set initial velocities from particle state
        body->setLinearVelocity(FoamBulletConverter::toBullet(pp->velocity()));
        body->setAngularVelocity(FoamBulletConverter::toBullet(pp->angularVelocity()));

        // Add some damping to reduce scattering and help settling
        body->setDamping(0.04, 0.06);

        // Add to world
        world_->addRigidBody(body.get());

        // Store
        particles_.append(pp);
        shapes_.push_back(std::move(shape));
        motionStates_.push_back(std::move(motionState));
        bodies_.push_back(std::move(body));
    }
}

void bulletIntegrator::addStaticBox(const Foam::point& center,
                                    const Foam::vector& halfExtents,
                                    const Foam::quaternion& orientation,
                                    Foam::scalar restitution,
                                    Foam::scalar friction)
{
    // Create shape
    std::unique_ptr<btCollisionShape> shape(new btBoxShape(btVector3(halfExtents.x(), halfExtents.y(), halfExtents.z())));

    // Transform
    btTransform t = FoamBulletConverter::toBulletTransform(center, orientation);
    std::unique_ptr<btDefaultMotionState> ms(new btDefaultMotionState(t));

    // Static body (mass=0)
    btVector3 inertia(0,0,0);
    btRigidBody::btRigidBodyConstructionInfo info(0.0, ms.get(), shape.get(), inertia);
    info.m_restitution = restitution;
    info.m_friction = friction;

    std::unique_ptr<btRigidBody> body(new btRigidBody(info));
    world_->addRigidBody(body.get());

    // Keep ownership
    shapes_.push_back(std::move(shape));
    motionStates_.push_back(std::move(ms));
    bodies_.push_back(std::move(body));
}

void bulletIntegrator::step(Foam::scalar dt, int maxSubSteps)
{
    if (!world_) return;
    world_->stepSimulation(dt, maxSubSteps);

    // Synchronize back to particles
    forAll(particles_, i)
    {
        btRigidBody* body = bodies_[i].get();
        btTransform trans;
        body->getMotionState()->getWorldTransform(trans);

        Foam::point pos;
        Foam::quaternion ori;
        FoamBulletConverter::toFoamTransform(trans, pos, ori);

        particles_[i]->setPosition(pos);
        particles_[i]->setOrientation(ori);
        particles_[i]->setVelocity(FoamBulletConverter::toFoam(body->getLinearVelocity()));
        particles_[i]->setAngularVelocity(FoamBulletConverter::toFoam(body->getAngularVelocity()));
    }
}

} // namespace Bashyal
