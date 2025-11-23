#ifndef wall_H
#define wall_H

#include "point.H"
#include "vector.H"
#include "tensor.H"
#include "List.H"
#include "boundary.H"

namespace Bashyal
{
    /**
     * @enum WallType
     * @brief Types of walls supported in the DEM simulation
     */
    enum class WallType
    {
        PLANE,      // Infinite plane wall
        CYLINDER,   // Cylindrical wall
        SPHERE      // Spherical wall (for curved boundaries)
    };

    /**
     * @class wall
     * @brief Represents a wall boundary for DEM simulations
     * 
     * Walls can be planes, cylinders, or spheres. They are used for
     * particle-wall contact detection and force computation.
     * Based on the Grains3D paper implementation.
     */
    class wall
    {
    private:
        //- Wall type
        WallType type_;

        //- Wall position/center
        Foam::point position_;

        //- Wall normal (for plane) or axis (for cylinder)
        Foam::vector normal_;

        //- Wall radius (for cylinder/sphere)
        Foam::scalar radius_;

        //- Wall height (for cylinder)
        Foam::scalar height_;

        //- Wall name
        Foam::word name_;

        //- Contact parameters (same as particles)
        Foam::scalar kn_;          // Normal contact stiffness
        Foam::scalar gammaN_;      // Normal damping coefficient
        Foam::scalar gammaT_;      // Tangential damping coefficient
        Foam::scalar muC_;         // Coulomb friction coefficient
        Foam::scalar crustWidth_;  // Crust width for GJK (homothety)

    public:
        // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

        /**
         * @brief Default constructor
         */
        wall();

        /**
         * @brief Constructor for a plane wall
         * @param position A point on the plane
         * @param normal Normal vector to the plane (pointing outward from wall)
         * @param name Wall name
         */
        wall(
            const Foam::point& position,
            const Foam::vector& normal,
            const Foam::word& name = "planeWall"
        );

        /**
         * @brief Constructor for a cylindrical wall
         * @param center Center point of cylinder
         * @param axis Axis direction of cylinder
         * @param radius Cylinder radius
         * @param height Cylinder height
         * @param name Wall name
         */
        wall(
            const Foam::point& center,
            const Foam::vector& axis,
            const Foam::scalar radius,
            const Foam::scalar height,
            const Foam::word& name = "cylinderWall"
        );

        /**
         * @brief Constructor for a spherical wall
         * @param center Center point of sphere
         * @param radius Sphere radius
         * @param name Wall name
         */
        wall(
            const Foam::point& center,
            const Foam::scalar radius,
            const Foam::word& name = "sphereWall"
        );

        /**
         * @brief Destructor
         */
        ~wall();

        // * * * * * * * * * * * * * * * * Accessors * * * * * * * * * * * * * * * * //

        /**
         * @brief Get wall type
         */
        WallType type() const { return type_; }

        /**
         * @brief Get wall position/center
         */
        const Foam::point& position() const { return position_; }

        /**
         * @brief Get wall normal/axis
         */
        const Foam::vector& normal() const { return normal_; }

        /**
         * @brief Get wall radius
         */
        Foam::scalar radius() const { return radius_; }

        /**
         * @brief Get wall height
         */
        Foam::scalar height() const { return height_; }

        /**
         * @brief Get wall name
         */
        const Foam::word& name() const { return name_; }

        /**
         * @brief Set contact parameters
         */
        void setContactParameters(
            const Foam::scalar kn,
            const Foam::scalar gammaN,
            const Foam::scalar gammaT,
            const Foam::scalar muC,
            const Foam::scalar crustWidth
        );

        /**
         * @brief Get contact parameters
         */
        Foam::scalar kn() const { return kn_; }
        Foam::scalar gammaN() const { return gammaN_; }
        Foam::scalar gammaT() const { return gammaT_; }
        Foam::scalar muC() const { return muC_; }
        Foam::scalar crustWidth() const { return crustWidth_; }

        /**
         * @brief Check if a point is inside the wall (for plane: on negative side of normal)
         */
        bool isInside(const Foam::point& p) const;

        /**
         * @brief Get distance from point to wall surface
         * @param p Point to check
         * @return Signed distance (negative if inside, positive if outside)
         */
        Foam::scalar distanceToSurface(const Foam::point& p) const;

        /**
         * @brief Get closest point on wall surface to given point
         * @param p Point to check
         * @return Closest point on wall
         */
        Foam::point closestPointOnSurface(const Foam::point& p) const;

        /**
         * @brief Get normal vector at a point on the wall surface
         * @param p Point on or near the wall
         * @return Normal vector (pointing outward from wall)
         */
        Foam::vector normalAtPoint(const Foam::point& p) const;
    };
}

#endif

