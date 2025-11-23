#ifndef gjkContact_H
#define gjkContact_H

#include "point.H"
#include "vector.H"
#include "List.H"
#include "DynamicList.H"
#include "particle.H"
#include "wall.H"

namespace Bashyal
{
    /**
     * @struct ContactInfo
     * @brief Stores contact information between two objects
     */
    struct ContactInfo
    {
        bool isContact;              // True if contact detected
        Foam::List<Foam::point> contactPoints; // List of contact points
        Foam::vector normal;         // Normal vector at contact (unit vector)
        Foam::scalar overlap;        // Overlap distance (negative if no contact) - NOT USED much in new model
        Foam::scalar depth;          // Penetration depth (positive)
        
        // Legacy/Debug fields
        Foam::point contactPoint;    // Average contact point
        Foam::point pointA;          // Closest point on object A
        Foam::point pointB;          // Closest point on object B
        Foam::scalar distance;       // Distance between objects (positive if separated)

        ContactInfo()
        :
            isContact(false),
            contactPoints(0),
            normal(Foam::vector(0, 0, 1)),
            overlap(Foam::GREAT),
            depth(0.0),
            contactPoint(Foam::vector::zero),
            pointA(Foam::vector::zero),
            pointB(Foam::vector::zero),
            distance(Foam::GREAT)
        {}
    };

    /**
     * @class gjkContact
     * @brief GJK-based contact detection for convex polyhedra
     * 
     * Implements the Gilbert-Johnson-Keerthi (GJK) algorithm for computing
     * the distance between convex objects and detecting contacts.
     * Based on the Grains3D paper implementation.
     */
    /**
     * @class gjkContact
     * @brief GJK-based contact detection for convex polyhedra
     * 
     * Implements the Gilbert-Johnson-Keerthi (GJK) algorithm for collision detection,
     * Expanding Polytope Algorithm (EPA) for penetration depth,
     * and Sutherland-Hodgman clipping for contact manifold generation.
     */
    class gjkContact
    {
    private:
        //- Maximum number of iterations for GJK
        static const int maxIterations_ = 64;
        
        //- Maximum number of iterations for EPA
        static const int maxEpaIterations_ = 64;

        //- Tolerance for convergence
        static const Foam::scalar tolerance_;

        //- Internal structure for Simplex
        struct Simplex
        {
            Foam::DynamicList<Foam::point> points; // Changed List to DynamicList
            
            void add(const Foam::point& p)
            {
                points.append(p);
            }
            
            void clear()
            {
                points.clear();
            }
            
            Foam::label size() const
            {
                return points.size();
            }
            
            const Foam::point& operator[](Foam::label i) const
            {
                return points[i];
            }
            
            Foam::point& operator[](Foam::label i)
            {
                return points[i];
            }
        };

        /**
         * @brief Support function for Minkowski difference
         * Returns the point in the Minkowski difference A - B that is
         * furthest in the direction d
         */
        Foam::point support(
            const particle& pA,
            const particle& pB,
            const Foam::vector& d
        ) const;

        /**
         * @brief Support function for particle-wall Minkowski difference
         */
        Foam::point supportParticleWall(
            const particle& p,
            const wall& w,
            const Foam::vector& d
        ) const;

        /**
         * @brief Get support point from a particle in direction d
         */
        Foam::point getParticleSupport(
            const particle& p,
            const Foam::vector& d
        ) const;

        /**
         * @brief Get support point from a wall in direction d
         */
        Foam::point getWallSupport(
            const wall& w,
            const Foam::vector& d
        ) const;

        /**
         * @brief Handle Simplex evolution (Line, Triangle, Tetrahedron)
         * @param simplex Current simplex
         * @param direction Output: new search direction
         * @return True if origin is enclosed (collision)
         */
        bool doSimplex(
            Simplex& simplex,
            Foam::vector& direction
        ) const;

        /**
         * @brief EPA Algorithm to find penetration depth and normal
         * @param simplex Final simplex from GJK (tetrahedron)
         * @param pA First particle
         * @param pB Second particle
         * @param depth Output: penetration depth
         * @param normal Output: contact normal (from B to A)
         * @return True if successful
         */
        bool epa(
            const Simplex& simplex,
            const particle& pA,
            const particle& pB,
            Foam::scalar& depth,
            Foam::vector& normal
        ) const;
        
        /**
         * @brief EPA for particle-wall
         */
        bool epaParticleWall(
            const Simplex& simplex,
            const particle& p,
            const wall& w,
            Foam::scalar& depth,
            Foam::vector& normal
        ) const;

        /**
         * @brief Generate contact manifold using Sutherland-Hodgman clipping
         * @param pA First particle
         * @param pB Second particle
         * @param normal Contact normal
         * @param contactPoints Output: list of contact points
         */
        void generateManifold(
            const particle& pA,
            const particle& pB,
            const Foam::vector& normal,
            Foam::List<Foam::point>& contactPoints
        ) const;

        /**
         * @brief Generate manifold for particle-wall
         */
        void generateManifoldParticleWall(
            const particle& p,
            const wall& w,
            const Foam::vector& normal,
            Foam::List<Foam::point>& contactPoints
        ) const;

    public:
        /**
         * @brief Default constructor
         */
        gjkContact();

        /**
         * @brief Destructor
         */
        ~gjkContact();

        /**
         * @brief Detect contact between two particles
         * @param pA First particle
         * @param pB Second particle
         * @param crustWidthA Crust width for particle A (homothety)
         * @param crustWidthB Crust width for particle B (homothety)
         * @return ContactInfo structure with contact details
         */
        ContactInfo detectContact(
            const particle& pA,
            const particle& pB,
            const Foam::scalar crustWidthA,
            const Foam::scalar crustWidthB
        ) const;

        /**
         * @brief Detect contact between particle and wall
         * @param p Particle
         * @param w Wall
         * @param crustWidthP Crust width for particle (homothety)
         * @return ContactInfo structure with contact details
         */
        ContactInfo detectContactParticleWall(
            const particle& p,
            const wall& w,
            const Foam::scalar crustWidthP
        ) const;
    };
}

#endif

