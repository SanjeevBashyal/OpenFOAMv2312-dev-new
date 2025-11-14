#ifndef FOAM_BULLET_CONVERTER_H
#define FOAM_BULLET_CONVERTER_H

#include "quickInclude.H"
#include "quaternion.H"

// Bullet
#include <btBulletDynamicsCommon.h>

class FoamBulletConverter
{
public:
    // Vector conversions
    static inline btVector3 toBullet(const Foam::vector& v)
    {
        return btVector3(v.x(), v.y(), v.z());
    }
    static inline Foam::vector toFoam(const btVector3& v)
    {
        return Foam::vector(v.x(), v.y(), v.z());
    }

    // Point conversions (maps to btVector3)
    static inline btVector3 toBulletPoint(const Foam::point& p)
    {
        return btVector3(p.x(), p.y(), p.z());
    }
    static inline Foam::point toFoamPoint(const btVector3& v)
    {
        return Foam::point(v.x(), v.y(), v.z());
    }

    // Quaternion conversions
    static inline btQuaternion toBullet(const Foam::quaternion& q)
    {
        return btQuaternion(q.v().x(), q.v().y(), q.v().z(), q.w());
    }
    static inline Foam::quaternion toFoam(const btQuaternion& q)
    {
        return Foam::quaternion(q.w(), Foam::vector(q.x(), q.y(), q.z()));
    }

    // Transform conversions
    static inline btTransform toBulletTransform(const Foam::point& pos, const Foam::quaternion& q)
    {
        btTransform t;
        t.setIdentity();
        t.setOrigin(toBulletPoint(pos));
        t.setRotation(toBullet(q));
        return t;
    }

    static inline void toFoamTransform(const btTransform& t, Foam::point& pos, Foam::quaternion& q)
    {
        pos = toFoamPoint(t.getOrigin());
        q = toFoam(t.getRotation());
    }

    // Append points to a convex hull shape
    static void appendConvexHullPoints(btConvexHullShape& hull, const Foam::pointField& points);
};

#endif // FOAM_BULLET_CONVERTER_H
