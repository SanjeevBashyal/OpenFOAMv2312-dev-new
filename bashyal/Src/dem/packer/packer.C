#include "packer.H"
#include "mathematicalConstants.H"

namespace Bashyal
{
    packer::packer(const Foam::point& domainMin, const Foam::point& domainMax)
    : rng_(1234), domain_(domainMin, domainMax)
    {
        // Build Bullet collision-only world
        broadphase_ = std::make_unique<btDbvtBroadphase>();
        collisionConfig_ = std::make_unique<btDefaultCollisionConfiguration>();
        dispatcher_ = std::make_unique<btCollisionDispatcher>(collisionConfig_.get());
        world_ = std::make_unique<btCollisionWorld>(dispatcher_.get(), broadphase_.get(), collisionConfig_.get());
    }

    packer::~packer()
    {
        // Collision objects must be removed before their shapes are destroyed
        for (int i = world_->getNumCollisionObjects() - 1; i >= 0; --i)
        {
            btCollisionObject* obj = world_->getCollisionObjectArray()[i];
            world_->removeCollisionObject(obj);
        }
        objects_.clear();
        shapes_.clear();
    }

    Foam::point packer::randomPositionWithin(const Foam::vector& halfExtents) const
    {
        const Foam::point& mn = domain_.min();
        const Foam::point& mx = domain_.max();

        Foam::scalar x = uniform(mn.x() + halfExtents.x(), mx.x() - halfExtents.x());
        Foam::scalar y = uniform(mn.y() + halfExtents.y(), mx.y() - halfExtents.y());
        Foam::scalar z = uniform(mn.z() + halfExtents.z(), mx.z() - halfExtents.z());
        return Foam::point(x, y, z);
    }

bool packer::hasCollision(btCollisionObject* obj) const
{
    struct ContactCB : public btCollisionWorld::ContactResultCallback
    {
        bool hasHit = false;
        btScalar addSingleResult(
            btManifoldPoint&,
            const btCollisionObjectWrapper*, int, int,
            const btCollisionObjectWrapper*, int, int) override
        {
            hasHit = true;
            return 0;
        }
    } cb;

    for (int i = 0; i < world_->getNumCollisionObjects(); ++i)
    {
        btCollisionObject* other = world_->getCollisionObjectArray()[i];
        if (other == obj) continue;
        cb.hasHit = false;
        world_->contactPairTest(obj, other, cb);
        if (cb.hasHit) return true;
    }
    return false;
}

    std::unique_ptr<btCollisionShape> packer::makeSphereShape(Foam::scalar radius) const
    {
        return std::unique_ptr<btCollisionShape>(new btSphereShape(radius));
    }

    std::unique_ptr<btCollisionShape> packer::makeBoxShape(const Foam::vector& halfExtents) const
    {
        return std::unique_ptr<btCollisionShape>(new btBoxShape(btVector3(halfExtents.x(), halfExtents.y(), halfExtents.z())));
    }

    std::unique_ptr<btCollisionShape> packer::makeHullShape(const Foam::pointField& points) const
    {
        std::unique_ptr<btConvexHullShape> hull(new btConvexHullShape());
        FoamBulletConverter::appendConvexHullPoints(*hull, points);
        hull->recalcLocalAabb();
        return hull;
    }

    std::unique_ptr<btCollisionObject> packer::makeCollisionObject(btCollisionShape* shape, const Foam::point& pos, const Foam::quaternion& q) const
    {
        btTransform t = FoamBulletConverter::toBulletTransform(pos, q);
        std::unique_ptr<btCollisionObject> obj(new btCollisionObject());
        obj->setCollisionShape(shape);
        obj->setWorldTransform(t);
        return obj;
    }

    Foam::PtrList<Bashyal::particle> packer::packSpheres(
        const Foam::label n,
        const Foam::scalar radius,
        const int resolution,
        const int maxTrialsPerParticle)
    {
        Bashyal::particlesGenerator gen;
        Foam::PtrList<Bashyal::particle> particles(n);

        // Create and store a single Bullet shape for all instances to share
        std::unique_ptr<btCollisionShape> baseShape = makeSphereShape(radius);
        shapes_.push_back(std::move(baseShape));
        btCollisionShape* shapePtr = shapes_.back().get();
        for (Foam::label i = 0; i < n; ++i)
        {
            bool placed = false;
            for (int trial = 0; trial < maxTrialsPerParticle && !placed; ++trial)
            {
                Foam::vector halfExt(radius, radius, radius);
                Foam::point pos = randomPositionWithin(halfExt);
                Foam::quaternion rot = Foam::quaternion::I; // spheres don't need rotation

                std::unique_ptr<btCollisionObject> obj = makeCollisionObject(shapePtr, pos, rot);
                world_->addCollisionObject(obj.get());

                if (!hasCollision(obj.get()))
                {
                    // Accept position; keep the object
                    objects_.push_back(std::move(obj));

                    // Generate particle at origin then place it
                    Foam::PtrList<Bashyal::particle> one = gen.createSpherical(1, radius, resolution);
                    one[0].setPosition(pos);
                    particles.set(i, one.release(0));
                    placed = true;
                    continue;
                }

                world_->removeCollisionObject(obj.get());
            }

            if (!placed)
            {
                FatalErrorInFunction << "Failed to place sphere " << i << "/" << n << Foam::nl << Foam::abort(Foam::FatalError);
            }
        }
        return particles;
    }

    Foam::PtrList<Bashyal::particle> packer::packCubes(
        const Foam::label n,
        const Foam::scalar side,
        const int maxTrialsPerParticle)
    {
        Bashyal::particlesGenerator gen;
        Foam::PtrList<Bashyal::particle> particles(n);

        // Half extents and Bullet shape
        Foam::vector halfExt(side/2, side/2, side/2);
        std::unique_ptr<btCollisionShape> baseShape = makeBoxShape(halfExt);
        shapes_.push_back(std::move(baseShape));
        btCollisionShape* shapePtr = shapes_.back().get();

        for (Foam::label i = 0; i < n; ++i)
        {
            bool placed = false;
            for (int trial = 0; trial < maxTrialsPerParticle && !placed; ++trial)
            {
                Foam::point pos = randomPositionWithin(halfExt);
                // Random orientation to avoid alignment artifacts
                Foam::quaternion rot(Foam::vector(uniform(0,1), uniform(0,1), uniform(0,1)), uniform(0, Foam::constant::mathematical::twoPi));

                std::unique_ptr<btCollisionObject> obj = makeCollisionObject(shapePtr, pos, rot);
                world_->addCollisionObject(obj.get());

                if (!hasCollision(obj.get()))
                {
                    objects_.push_back(std::move(obj));

                    Foam::PtrList<Bashyal::particle> one = gen.createCubical(1, side, side);
                    one[0].setPosition(pos);
                    one[0].setOrientation(rot);
                    particles.set(i, one.release(0));
                    placed = true;
                    continue;
                }
                world_->removeCollisionObject(obj.get());
            }

            if (!placed)
            {
                FatalErrorInFunction << "Failed to place cube " << i << "/" << n << Foam::nl << Foam::abort(Foam::FatalError);
            }
        }

        return particles;
    }

    Foam::PtrList<Bashyal::particle> packer::packCuboids(
        const Foam::label n,
        const Foam::vector& dims,
        const int maxTrialsPerParticle)
    {
        Bashyal::particlesGenerator gen;
        Foam::PtrList<Bashyal::particle> particles(n);

        Foam::vector halfExt(dims.x()/2, dims.y()/2, dims.z()/2);
        std::unique_ptr<btCollisionShape> baseShape = makeBoxShape(halfExt);
        shapes_.push_back(std::move(baseShape));
        btCollisionShape* shapePtr = shapes_.back().get();

        for (Foam::label i = 0; i < n; ++i)
        {
            bool placed = false;
            for (int trial = 0; trial < maxTrialsPerParticle && !placed; ++trial)
            {
                Foam::point pos = randomPositionWithin(halfExt);
                Foam::quaternion rot(Foam::vector(uniform(0,1), uniform(0,1), uniform(0,1)), uniform(0, Foam::constant::mathematical::twoPi));

                std::unique_ptr<btCollisionObject> obj = makeCollisionObject(shapePtr, pos, rot);
                world_->addCollisionObject(obj.get());

                if (!hasCollision(obj.get()))
                {
                    objects_.push_back(std::move(obj));

                    Foam::PtrList<Bashyal::particle> one = gen.createCuboid(1, dims, dims);
                    one[0].setPosition(pos);
                    one[0].setOrientation(rot);
                    particles.set(i, one.release(0));
                    placed = true;
                    continue;
                }
                world_->removeCollisionObject(obj.get());
            }

            if (!placed)
            {
                FatalErrorInFunction << "Failed to place cuboid " << i << "/" << n << Foam::nl << Foam::abort(Foam::FatalError);
            }
        }

        return particles;
    }

    Foam::PtrList<Bashyal::particle> packer::packWellRounded(
        const Foam::label n,
        const Foam::scalar size,
        const Foam::scalar roundnessFactor,
        const int resolution,
        const int maxTrialsPerParticle)
    {
        Bashyal::particlesGenerator gen;
        Foam::PtrList<Bashyal::particle> particles(n);

        // Build a representative hull shape for collision: start with unit sphere points, scale by roundnessFactor ranges
        Foam::pointField basePts = gen.createSpherical(1, 1.0, resolution)[0].vertices(); // Not available directly; will instead generate points here
        // Since generator returns particles not points, create sphere points directly
        basePts.clear();
        {
            const int pointsPerRing = resolution;
            const int totalPoints = (resolution - 1) * pointsPerRing + 2;
            basePts.setSize(totalPoints);
            basePts[0] = Foam::point(0,0,1);
            for (int i = 1; i < resolution; ++i)
            {
                Foam::scalar theta = Foam::constant::mathematical::pi * i / resolution;
                Foam::scalar sT = Foam::sin(theta), cT = Foam::cos(theta);
                for (int j = 0; j < resolution; ++j)
                {
                    Foam::scalar phi = 2 * Foam::constant::mathematical::pi * j / resolution;
                    int index = (i - 1) * pointsPerRing + j + 1;
                    basePts[index] = Foam::point(Foam::cos(phi)*sT, Foam::sin(phi)*sT, cT);
                }
            }
            basePts[totalPoints-1] = Foam::point(0,0,-1);
        }

        std::unique_ptr<btCollisionShape> baseShape = makeHullShape(basePts);
        shapes_.push_back(std::move(baseShape));
        btCollisionShape* shapePtr = shapes_.back().get();

        for (Foam::label i = 0; i < n; ++i)
        {
            bool placed = false;
            for (int trial = 0; trial < maxTrialsPerParticle && !placed; ++trial)
            {
                // Conservative half extents from size
                Foam::vector halfExt(size/2, size/2, size/2);
                Foam::point pos = randomPositionWithin(halfExt);
                Foam::quaternion rot(Foam::vector(uniform(0,1), uniform(0,1), uniform(0,1)), uniform(0, Foam::constant::mathematical::twoPi));

                std::unique_ptr<btCollisionObject> obj = makeCollisionObject(shapePtr, pos, rot);
                world_->addCollisionObject(obj.get());

                if (!hasCollision(obj.get()))
                {
                    objects_.push_back(std::move(obj));

                    Foam::PtrList<Bashyal::particle> one = gen.createWellRounded(1, size, size, roundnessFactor, resolution);
                    one[0].setPosition(pos);
                    one[0].setOrientation(rot);
                    particles.set(i, one.release(0));
                    placed = true;
                    continue;
                }
                world_->removeCollisionObject(obj.get());
            }

            if (!placed)
            {
                FatalErrorInFunction << "Failed to place wellRounded particle " << i << "/" << n << Foam::nl << Foam::abort(Foam::FatalError);
            }
        }

        return particles;
    }
}


