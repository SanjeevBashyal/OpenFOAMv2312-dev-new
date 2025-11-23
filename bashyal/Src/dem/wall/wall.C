#include "wall.H"
#include "constants.H"
#include "scalar.H"
#include <cmath>

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

    wall::wall()
    :
        type_(WallType::PLANE),
        position_(Foam::vector::zero),
        normal_(Foam::vector(0, 0, 1)),
        radius_(0.0),
        height_(0.0),
        name_("defaultWall"),
        kn_(1e5),
        gammaN_(1e3),
        gammaT_(1e5),
        muC_(0.5),
        crustWidth_(1e-4)
    {}

    wall::wall(
        const Foam::point& position,
        const Foam::vector& normal,
        const Foam::word& name
    )
    :
        type_(WallType::PLANE),
        position_(position),
        normal_(normal),
        radius_(0.0),
        height_(0.0),
        name_(name),
        kn_(1e5),
        gammaN_(1e3),
        gammaT_(1e5),
        muC_(0.5),
        crustWidth_(1e-4)
    {
        // Normalize normal vector
        Foam::scalar magVal = Foam::mag(normal_);
        if (magVal > Foam::SMALL)
        {
            normal_ /= magVal;
        }
        else
        {
            normal_ = Foam::vector(0, 0, 1);
        }
    }

    wall::wall(
        const Foam::point& center,
        const Foam::vector& axis,
        const Foam::scalar radius,
        const Foam::scalar height,
        const Foam::word& name
    )
    :
        type_(WallType::CYLINDER),
        position_(center),
        normal_(axis),
        radius_(radius),
        height_(height),
        name_(name),
        kn_(1e5),
        gammaN_(1e3),
        gammaT_(1e5),
        muC_(0.5),
        crustWidth_(1e-4)
    {
        // Normalize axis vector
        Foam::scalar magVal = Foam::mag(normal_);
        if (magVal > Foam::SMALL)
        {
            normal_ /= magVal;
        }
        else
        {
            normal_ = Foam::vector(0, 0, 1);
        }
    }

    wall::wall(
        const Foam::point& center,
        const Foam::scalar radius,
        const Foam::word& name
    )
    :
        type_(WallType::SPHERE),
        position_(center),
        normal_(Foam::vector::zero),
        radius_(radius),
        height_(0.0),
        name_(name),
        kn_(1e5),
        gammaN_(1e3),
        gammaT_(1e5),
        muC_(0.5),
        crustWidth_(1e-4)
    {}

    wall::~wall()
    {}

    // * * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

    void wall::setContactParameters(
        const Foam::scalar kn,
        const Foam::scalar gammaN,
        const Foam::scalar gammaT,
        const Foam::scalar muC,
        const Foam::scalar crustWidth
    )
    {
        kn_ = kn;
        gammaN_ = gammaN;
        gammaT_ = gammaT;
        muC_ = muC;
        crustWidth_ = crustWidth;
    }

    bool wall::isInside(const Foam::point& p) const
    {
        switch (type_)
        {
            case WallType::PLANE:
            {
                // For plane: inside means on negative side of normal
                Foam::vector toPoint = p - position_;
                return (toPoint & normal_) < 0.0;
            }
            case WallType::CYLINDER:
            {
                // For cylinder: inside means distance from axis < radius
                Foam::vector toPoint = p - position_;
                Foam::vector alongAxis = (toPoint & normal_) * normal_;
                Foam::vector radial = toPoint - alongAxis;
                return Foam::mag(radial) < radius_;
            }
            case WallType::SPHERE:
            {
                // For sphere: inside means distance from center < radius
                Foam::vector toPoint = p - position_;
                return Foam::mag(toPoint) < radius_;
            }
            default:
                return false;
        }
    }

    Foam::scalar wall::distanceToSurface(const Foam::point& p) const
    {
        switch (type_)
        {
            case WallType::PLANE:
            {
                Foam::vector toPoint = p - position_;
                return toPoint & normal_;
            }
            case WallType::CYLINDER:
            {
                Foam::vector toPoint = p - position_;
                Foam::vector alongAxis = (toPoint & normal_) * normal_;
                Foam::vector radial = toPoint - alongAxis;
                Foam::scalar distToAxis = Foam::mag(radial);
                return distToAxis - radius_;
            }
            case WallType::SPHERE:
            {
                Foam::vector toPoint = p - position_;
                return Foam::mag(toPoint) - radius_;
            }
            default:
                return Foam::GREAT;
        }
    }

    Foam::point wall::closestPointOnSurface(const Foam::point& p) const
    {
        switch (type_)
        {
            case WallType::PLANE:
            {
                Foam::vector toPoint = p - position_;
                Foam::scalar dist = toPoint & normal_;
                return p - dist * normal_;
            }
            case WallType::CYLINDER:
            {
                Foam::vector toPoint = p - position_;
                Foam::vector alongAxis = (toPoint & normal_) * normal_;
                Foam::vector radial = toPoint - alongAxis;
                Foam::scalar magRadial = Foam::mag(radial);
                if (magRadial > Foam::SMALL)
                {
                    Foam::vector radialDir = radial / magRadial;
                    return position_ + alongAxis + radius_ * radialDir;
                }
                else
                {
                    // Point is on axis, return arbitrary point on surface
                    return position_ + alongAxis + radius_ * Foam::vector(1, 0, 0);
                }
            }
            case WallType::SPHERE:
            {
                Foam::vector toPoint = p - position_;
                Foam::scalar magToPoint = Foam::mag(toPoint);
                if (magToPoint > Foam::SMALL)
                {
                    return position_ + radius_ * (toPoint / magToPoint);
                }
                else
                {
                    // Point is at center, return arbitrary point on surface
                    return position_ + radius_ * Foam::vector(1, 0, 0);
                }
            }
            default:
                return p;
        }
    }

    Foam::vector wall::normalAtPoint(const Foam::point& p) const
    {
        switch (type_)
        {
            case WallType::PLANE:
            {
                return normal_;
            }
            case WallType::CYLINDER:
            {
                Foam::vector toPoint = p - position_;
                Foam::vector alongAxis = (toPoint & normal_) * normal_;
                Foam::vector radial = toPoint - alongAxis;
                Foam::scalar magRadial = Foam::mag(radial);
                if (magRadial > Foam::SMALL)
                {
                    return radial / magRadial;
                }
                else
                {
                    // Point is on axis, return arbitrary normal
                    return Foam::vector(1, 0, 0);
                }
            }
            case WallType::SPHERE:
            {
                Foam::vector toPoint = p - position_;
                Foam::scalar magToPoint = Foam::mag(toPoint);
                if (magToPoint > Foam::SMALL)
                {
                    return toPoint / magToPoint;
                }
                else
                {
                    // Point is at center, return arbitrary normal
                    return Foam::vector(1, 0, 0);
                }
            }
            default:
                return Foam::vector(0, 0, 1);
        }
    }
}

