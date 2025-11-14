#ifndef packer_H
#define packer_H

#include "quickInclude.H"
#include "PtrList.H"
#include "boundBox.H"
#include "Random.H"

#include "particle.H"
#include "particlesGenerator.H"
#include "foamBulletConverter.H"

#include <memory>
#include <vector>

// Bullet
#include <btBulletDynamicsCommon.h>

namespace Bashyal
{
    class packer
    {
    private:
        mutable Foam::Random rng_;
        Foam::boundBox domain_;

        // Bullet collision world (no dynamics needed)
        std::unique_ptr<btBroadphaseInterface> broadphase_;
        std::unique_ptr<btDefaultCollisionConfiguration> collisionConfig_;
        std::unique_ptr<btCollisionDispatcher> dispatcher_;
        std::unique_ptr<btCollisionWorld> world_;

        // Keep ownership of shapes and objects added to the world
        std::vector<std::unique_ptr<btCollisionShape>> shapes_;
        std::vector<std::unique_ptr<btCollisionObject>> objects_;

        // Helper: random in [a,b]
        inline Foam::scalar uniform(Foam::scalar a, Foam::scalar b) const { return a + (b - a) * rng_.sample01<Foam::scalar>(); }

        // Helper: build box AABB-constrained random position given half extents
        Foam::point randomPositionWithin(const Foam::vector& halfExtents) const;

        // Helper: Bullet collision check for a candidate object (already added to world_)
        bool hasCollision(btCollisionObject* obj) const;

        // Create Bullet shapes
        std::unique_ptr<btCollisionShape> makeSphereShape(Foam::scalar radius) const;
        std::unique_ptr<btCollisionShape> makeBoxShape(const Foam::vector& halfExtents) const;
        std::unique_ptr<btCollisionShape> makeHullShape(const Foam::pointField& points) const;

        // Build a Bullet collision object with given shape and transform
        std::unique_ptr<btCollisionObject> makeCollisionObject(btCollisionShape* shape, const Foam::point& pos, const Foam::quaternion& q) const;

    public:
        packer(const Foam::point& domainMin, const Foam::point& domainMax);
        ~packer();

        // Non-overlapping packing methods (definite sizes)
        Foam::PtrList<Bashyal::particle> packSpheres(
            const Foam::label n,
            const Foam::scalar radius,
            const int resolution,
            const int maxTrialsPerParticle = 1000);

        Foam::PtrList<Bashyal::particle> packCubes(
            const Foam::label n,
            const Foam::scalar side,
            const int maxTrialsPerParticle = 1000);

        Foam::PtrList<Bashyal::particle> packCuboids(
            const Foam::label n,
            const Foam::vector& dims,
            const int maxTrialsPerParticle = 1000);

        Foam::PtrList<Bashyal::particle> packWellRounded(
            const Foam::label n,
            const Foam::scalar size,
            const Foam::scalar roundnessFactor,
            const int resolution,
            const int maxTrialsPerParticle = 1000);
    };
}

#endif


