#ifndef contactForce_H
#define contactForce_H

#include "vector.H"
#include "particle.H"
#include "wall.H"
#include "gjkContact.H"

namespace Bashyal
{
    /**
     * @class contactForce
     * @brief Computes contact forces between particles and walls
     * 
     * Implements the contact force model from the Grains3D paper:
     * - Hookean elastic force
     * - Viscous damping force
     * - Tangential friction force (Coulomb limit)
     */
    struct ContactHistory
    {
        Foam::vector tangentialSpring;
        Foam::scalar maxOverlap;
        bool isLoading;
        
        ContactHistory()
        :
            tangentialSpring(Foam::vector::zero),
            maxOverlap(0.0),
            isLoading(true)
        {}
    };

    class contactForce
    {
    private:
        //- Normal stiffness (loading)
        Foam::scalar kn_;
        
        //- Normal stiffness (unloading)
        Foam::scalar knUnloading_;
        
        //- Tangential stiffness
        Foam::scalar kt_;
        
        //- Normal damping
        Foam::scalar gammaN_;
        
        //- Tangential damping
        Foam::scalar gammaT_;
        
        //- Coulomb friction coefficient
        Foam::scalar muC_;
        
        //- Contact history: Key is "ID1_ID2"
        mutable Foam::HashTable<ContactHistory, Foam::string> history_;

    public:
        //- Constructor
        contactForce(
            const Foam::scalar kn,
            const Foam::scalar knUnloading,
            const Foam::scalar kt,
            const Foam::scalar gammaN,
            const Foam::scalar gammaT,
            const Foam::scalar muC
        );

        //- Destructor
        ~contactForce() {}

        //- Compute force between two particles
        void computeForce(
            particle& pA,
            particle& pB,
            const ContactInfo& contact,
            const Foam::scalar dt
        ) const;

        //- Compute force between particle and wall
        void computeForceParticleWall(
            particle& p,
            const wall& w,
            const ContactInfo& contact,
            const Foam::scalar dt
        ) const;
        
        //- Clean up history for non-existing contacts
        void cleanHistory(const Foam::List<Foam::string>& activeKeys) const;

    private:
        /**
         * @brief Compute relative velocity at contact point
         * @param pA First particle
         * @param pB Second particle
         * @param contactPoint Contact point in world coordinates
         * @param normal Contact normal
         * @return Relative velocity (normal and tangential components)
         */
        void computeRelativeVelocity(
            const particle& pA,
            const particle& pB,
            const Foam::point& contactPoint,
            const Foam::vector& normal,
            Foam::vector& vRelN,
            Foam::vector& vRelT
        ) const;

        /**
         * @brief Compute relative velocity for particle-wall contact
         */
        void computeRelativeVelocityParticleWall(
            const particle& p,
            const wall& w,
            const Foam::point& contactPoint,
            const Foam::vector& normal,
            Foam::vector& vRelN,
            Foam::vector& vRelT
        ) const;
    };
}

#endif

