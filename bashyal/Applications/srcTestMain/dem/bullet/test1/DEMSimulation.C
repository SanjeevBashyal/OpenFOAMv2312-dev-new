#include "DEMSimulation.H"
#include "Time.H"

namespace Bashyal
{
    // Helper function to convert from Foam::vector to btVector3
    inline btVector3 toBullet(const Foam::vector& v)
    {
        return btVector3(v.x(), v.y(), v.z());
    }

    // Helper function to convert from btVector3 to Foam::vector
    inline Foam::vector toFoam(const btVector3& v)
    {
        return Foam::vector(v.x(), v.y(), v.z());
    }

    // Helper function to convert from Foam::quaternion to btQuaternion
    inline btQuaternion toBullet(const Foam::quaternion& q)
    {
        // Note: OpenFOAM quaternion is (w, x, y, z), Bullet is (x, y, z, w)
        return btQuaternion(q.v().x(), q.v().y(), q.v().z(), q.w());
    }

    // Helper function to convert from btQuaternion to Foam::quaternion
    inline Foam::quaternion toFoam(const btQuaternion& q)
    {
        // Note: OpenFOAM quaternion is (w, x, y, z), Bullet is (x, y, z, w)
        return Foam::quaternion(q.w(), q.x(), q.y(), q.z());
    }


    DEMSimulation::DEMSimulation(const Foam::vector& gravity)
    {
        // --- Standard Bullet Physics Initialization ---
        collisionConfiguration_ = std::make_unique<btDefaultCollisionConfiguration>();
        dispatcher_ = std::make_unique<btCollisionDispatcher>(collisionConfiguration_.get());
        broadphase_ = std::make_unique<btDbvtBroadphase>();
        solver_ = std::make_unique<btSequentialImpulseConstraintSolver>();

        dynamicsWorld_ = std::make_unique<btDiscreteDynamicsWorld>(
            dispatcher_.get(),
            broadphase_.get(),
            solver_.get(),
            collisionConfiguration_.get()
        );

        // Set the gravity for the world
        dynamicsWorld_->setGravity(toBullet(gravity));
        Foam::Info << "DEMSimulation: Bullet physics world initialized." << Foam::endl;
    }

    DEMSimulation::~DEMSimulation()
    {
        // Cleanup must be in reverse order of creation
        for (int i = dynamicsWorld_->getNumCollisionObjects() - 1; i >= 0; --i)
        {
            btCollisionObject* obj = dynamicsWorld_->getCollisionObjectArray()[i];
            btRigidBody* body = btRigidBody::upcast(obj);
            if (body && body->getMotionState())
            {
                delete body->getMotionState();
            }
            dynamicsWorld_->removeCollisionObject(obj);
        }
        // Smart pointers will handle the rest of the cleanup automatically.
        Foam::Info << "DEMSimulation: Cleaned up Bullet physics world." << Foam::endl;
    }

    void DEMSimulation::addParticle(particle&& p)
    {
        // 1. Create the Collision Shape from the particle's geometry
        // We use a convex hull, which is very general.
        auto shape = std::make_unique<btConvexHullShape>();
        const Foam::pointField& localVerts = p.vertices();
        for (const auto& pt : localVerts)
        {
            shape->addPoint(toBullet(pt), false); // 'false' to skip recalculating AABB each time
        }
        shape->recalcLocalAabb();
        shape->optimizeConvexHull();

        // 2. Create the Motion State (initial position and orientation)
        btTransform startTransform;
        startTransform.setIdentity();
        startTransform.setOrigin(toBullet(p.position()));
        startTransform.setRotation(toBullet(p.orientation()));
        auto motionState = std::make_unique<btDefaultMotionState>(startTransform);

        // 3. Define the Rigid Body properties
        btVector3 localInertia(0, 0, 0);
        // Bullet can calculate the inertia tensor for its primitive shapes,
        // but for a convex hull, it's better to provide it from our particle.
        // For simplicity, we assume the provided MOI is for a dynamic object.
        if (p.mass() > 0)
        {
            // Note: Bullet expects the inverse of the moment of inertia tensor
            // if you calculate it yourself. Here we let it calculate for a box
            // shape for simplicity, but a real implementation might use the
            // tensor from the 'particle' class.
            shape->calculateLocalInertia(p.mass(), localInertia);
        }

        btRigidBody::btRigidBodyConstructionInfo rbInfo(p.mass(), motionState.get(), shape.get(), localInertia);
        rbInfo.m_restitution = 0.5f; // Bounciness
        rbInfo.m_friction = 0.5f;

        // 4. Create the Rigid Body and add it to the world
        auto body = std::make_unique<btRigidBody>(rbInfo);
        dynamicsWorld_->addRigidBody(body.get());

        // 5. Store all the created objects for later management
        particles_.append(std::move(p));
        collisionShapes_.push_back(std::move(shape));
        motionStates_.push_back(std::move(motionState));
        rigidBodies_.push_back(std::move(body));
    }

    void DEMSimulation::step(Foam::scalar dt)
    {
        // Tell Bullet to advance the simulation by dt seconds
        dynamicsWorld_->stepSimulation(dt, 10); // 10 max sub-steps for accuracy
    }

    void DEMSimulation::synchronizeFromBullet()
    {
        forAll(particles_, i)
        {
            btRigidBody* body = rigidBodies_[i].get();
            btTransform trans;
            // Get the latest transform from the motion state
            body->getMotionState()->getWorldTransform(trans);

            // Update our custom particle's state
            particles_[i].setPosition(toFoam(trans.getOrigin()));
            particles_[i].setOrientation(toFoam(trans.getRotation()));

            // Also update velocities (optional, but good for data analysis)
            particles_[i].setVelocity(toFoam(body->getLinearVelocity()));
            particles_[i].setAngularVelocity(toFoam(body->getAngularVelocity()));
        }
    }

    void DEMSimulation::writeOutput(Foam::scalar time) const
    {
        Foam::Info << "Writing output for time = " << time << Foam::endl;

        // Create a directory for this time step
        Foam::fileName timeDir = Foam::word(Foam::name(time));
        Foam::mkDir(timeDir);

        forAll(particles_, i)
        {
            Foam::fileName fName = timeDir / ("particle_" + Foam::word(Foam::name(i)) + ".vtp");
            particles_[i].writeVtp(fName);
        }
    }

    void DEMSimulation::run(Foam::scalar endTime, Foam::scalar dt, Foam::scalar writeInterval)
    {
        Foam::scalar currentTime = 0;
        Foam::scalar nextWriteTime = 0;

        // Write initial state
        writeOutput(currentTime);
        nextWriteTime += writeInterval;

        while (currentTime < endTime)
        {
            // Perform one simulation step
            step(dt);

            // Get the results from Bullet and update our particles
            synchronizeFromBullet();

            currentTime += dt;

            // Check if it's time to write output
            if (currentTime >= nextWriteTime)
            {
                writeOutput(currentTime);
                nextWriteTime += writeInterval;
            }
            Foam::Info << "Time = " << currentTime << "\r" << std::flush;
        }

        // Final write
        if (abs(currentTime - nextWriteTime) > dt)
        {
             writeOutput(currentTime);
        }
        Foam::Info << "\nSimulation finished at Time = " << currentTime << Foam::endl;
    }

} // End namespace Bashyal