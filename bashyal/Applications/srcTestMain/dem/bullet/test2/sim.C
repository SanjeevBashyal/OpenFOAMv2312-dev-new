// src/sim.cpp
#include <iostream>
#include <vector>
#include <memory>
#include <sstream>
#include <iomanip>

#include <btBulletDynamicsCommon.h>

// Include your particle class (uploaded file)
#include "particle.H"   // adjust relative path as needed

using namespace Bashyal;

// Helper: convert btVector3 -> Foam::point
static Foam::point btVec3ToFoamPoint(const btVector3 &v) {
    return Foam::point(Foam::scalar(v.x()), Foam::scalar(v.y()), Foam::scalar(v.z()));
}

// Helper: convert btVector3 -> Foam::vector
static Foam::vector btVec3ToFoamVector(const btVector3 &v) {
    return Foam::vector(Foam::scalar(v.x()), Foam::scalar(v.y()), Foam::scalar(v.z()));
}

// Helper: convert btQuaternion -> Foam::quaternion (axis-angle constructor)
static Foam::quaternion btQuatToFoamQuat(const btQuaternion &q_in) {
    // Ensure normalized
    btQuaternion q = q_in.normalized();
    btScalar w = q.getW();
    btScalar x = q.getX();
    btScalar y = q.getY();
    btScalar z = q.getZ();

    // angle = 2 * acos(w)
    btScalar angle = 2.0 * std::acos(w);
    btScalar s = std::sqrt(1.0 - w*w);

    Foam::vector axis(Foam::scalar(1.0), Foam::scalar(0.0), Foam::scalar(0.0));
    if (s > 1e-6) {
        axis = Foam::vector(Foam::scalar(x / s), Foam::scalar(y / s), Foam::scalar(z / s));
    }
    // Foam::quaternion(axis, angle) constructor used in particle.C above
    return Foam::quaternion(axis, Foam::scalar(angle));
}

// Wrapper tying a Bashyal::particle to a Bullet rigid body
struct RigidParticle
{
    std::unique_ptr<particle> part;           // the user's particle object
    btCollisionShape* shape = nullptr;        // owned by Bullet world cleanup (we'll delete)
    btDefaultMotionState* motionState = nullptr;
    btRigidBody* body = nullptr;
};

// Simple function to compute a sphere radius from particle's vertices (fallback)
static btScalar computeBoundingSphereRadius(const particle& p)
{
    // We assume boundary::vertices() is accessible. If not, use getGlobalPoints() or provide radius.
    Foam::pointField localPts = p.getGlobalPoints(); // global; but magnitude used to estimate radius
    btScalar r = 0.0;
    forAll(localPts, i)
    {
        const Foam::point &pt = localPts[i];
        btScalar dist = std::sqrt(pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z());
        if (dist > r) r = dist;
    }
    // small safety factor:
    return r > 1e-6 ? r : 0.1;
}

int main(int argc, char *argv[])
{
    std::cout << "Starting OpenFOAM-Bullet multi-particle sim\n";

    // 1) Bullet world setup
    btBroadphaseInterface* broadphase = new btDbvtBroadphase();
    btDefaultCollisionConfiguration* collisionConfiguration = new btDefaultCollisionConfiguration();
    btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);
    btSequentialImpulseConstraintSolver* solver = new btSequentialImpulseConstraintSolver();
    btDiscreteDynamicsWorld* dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration);

    // gravity (z-up? OpenFOAM default typically z up, adjust if needed)
    dynamicsWorld->setGravity(btVector3(0.f, -9.81f, 0.f)); // y-down used here; change axes to match your setup

    // 2) Create static ground plane (y = 0)
    btCollisionShape* groundShape = new btStaticPlaneShape(btVector3(0, 1, 0), 0); // plane normal up along +y
    btDefaultMotionState* groundMotionState = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1), btVector3(0,-0.0,0)));
    btRigidBody::btRigidBodyConstructionInfo groundRbInfo(0.0f, groundMotionState, groundShape, btVector3(0,0,0));
    btRigidBody* groundBody = new btRigidBody(groundRbInfo);
    dynamicsWorld->addRigidBody(groundBody);

    // 3) Create N particles and add them to Bullet world
    const size_t N = 12;
    std::vector<std::unique_ptr<RigidParticle>> items;
    items.reserve(N);

    for (size_t i = 0; i < N; ++i)
    {
        // Create particle geometry: for demo, we construct a boundary-less particle with mass and initial position.
        // In practice, you probably construct from an existing `boundary` describing the mesh.
        Foam::point initPos(Foam::scalar(0.0 + 0.25 * i), Foam::scalar(5.0 + 0.5 * i), Foam::scalar(0.0));
        Foam::tensor moi = Foam::tensor::I; // placeholder
        Foam::scalar mass = 1.0;

        auto rp = std::make_unique<RigidParticle>();
        rp->part = std::make_unique<particle>(/* boundary b - you may pass a boundary object here */ boundary(), mass, moi, initPos);

        // Determine collision shape (sphere using bounding radius)
        btScalar radius = computeBoundingSphereRadius(*rp->part);
        if (radius <= 0) radius = 0.2;
        rp->shape = new btSphereShape(radius);

        // initial transform
        btTransform startTransform;
        startTransform.setIdentity();
        startTransform.setOrigin(btVector3(initPos.x(), initPos.y(), initPos.z()));
        startTransform.setRotation(btQuaternion(0,0,0,1));

        // Create motion state and rigid body
        rp->motionState = new btDefaultMotionState(startTransform);
        btVector3 localInertia(0,0,0);
        rp->shape->calculateLocalInertia(mass, localInertia);
        btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, rp->motionState, rp->shape, localInertia);
        rp->body = new btRigidBody(rbInfo);

        // Give a small initial velocity for collisions
        rp->body->setLinearVelocity(btVector3(0.0f, -1.0f - 0.1f * float(i), 0.0f));

        // Add to world
        dynamicsWorld->addRigidBody(rp->body);

        // record in vector
        items.push_back(std::move(rp));
    }

    // 4) Simulation loop
    const double dt = 1.0/240.0;      // physics time-step
    const double simTime = 10.0;      // seconds to simulate
    const int steps = int(simTime / dt);
    const int writeEvery = 120;       // write VTP every 120 steps

    for (int step = 0; step < steps; ++step)
    {
        // advance Bullet
        dynamicsWorld->stepSimulation(dt, 10);

        // sync Bullet -> your particle objects
        for (size_t i = 0; i < items.size(); ++i)
        {
            RigidParticle &rp = *items[i];
            btTransform trans;
            rp.body->getMotionState()->getWorldTransform(trans);
            btVector3 pos = trans.getOrigin();
            btQuaternion btq = trans.getRotation();

            // Convert and set into particle
            rp.part->setPosition(btVec3ToFoamPoint(pos));
            rp.part->setOrientation(btQuatToFoamQuat(btq));

            // Update velocity fields in particle (optional)
            btVector3 linVel = rp.body->getLinearVelocity();
            rp.part->setVelocity(btVec3ToFoamVector(linVel));
            // Angular velocity -> convert to Foam::vector
            btVector3 angVel = rp.body->getAngularVelocity();
            rp.part->setAngularVelocity(btVec3ToFoamVector(angVel));
        }

        // write snapshot occasionally
        if (step % writeEvery == 0)
        {
            for (size_t i = 0; i < items.size(); ++i)
            {
                std::ostringstream name;
                name << "particle_" << std::setfill('0') << std::setw(3) << i << "_step_" << step << ".vtp";
                items[i]->part->writeVtp(Foam::fileName(name.str()));
            }
            std::cout << "Wrote VTP at step " << step << "\n";
        }
    }

    std::cout << "Simulation finished. Cleaning up.\n";

    // 5) Cleanup Bullet objects (delete in reverse creation order)
    for (auto &rp : items)
    {
        dynamicsWorld->removeRigidBody(rp->body);
        delete rp->body;
        delete rp->motionState;
        delete rp->shape;
        // unique_ptr rp->part will be auto deleted
    }

    dynamicsWorld->removeRigidBody(groundBody);
    delete groundBody;
    delete groundMotionState;
    delete groundShape;

    delete dynamicsWorld;
    delete solver;
    delete dispatcher;
    delete collisionConfiguration;
    delete broadphase;

    return 0;
}
