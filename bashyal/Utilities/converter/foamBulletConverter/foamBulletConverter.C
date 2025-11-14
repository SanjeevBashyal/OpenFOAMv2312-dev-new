#include "foamBulletConverter.H"

void FoamBulletConverter::appendConvexHullPoints(btConvexHullShape& hull, const Foam::pointField& points)
{
    forAll(points, i)
    {
        const Foam::point& p = points[i];
        hull.addPoint(btVector3(p.x(), p.y(), p.z()), false);
    }
}
