#include "contactForce.H"
#include "constants.H"
#include "scalar.H"

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    contactForce::contactForce(
        const Foam::scalar kn,
        const Foam::scalar knUnloading,
        const Foam::scalar kt,
        const Foam::scalar gammaN,
        const Foam::scalar gammaT,
        const Foam::scalar muC
    )
    :
        kn_(kn),
        knUnloading_(knUnloading),
        kt_(kt),
        gammaN_(gammaN),
        gammaT_(gammaT),
        muC_(muC),
        history_()
    {}

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void contactForce::cleanHistory(const Foam::List<Foam::string>& activeKeys) const
    {
        // Create a set of active keys for fast lookup
        Foam::HashSet<Foam::string> activeSet;
        forAll(activeKeys, i)
        {
            activeSet.insert(activeKeys[i]);
        }
        
        // Identify keys to remove
        Foam::List<Foam::string> keysToRemove;
        for (auto iter = history_.cbegin(); iter != history_.cend(); ++iter)
        {
            if (!activeSet.found(iter.key()))
            {
                keysToRemove.append(iter.key());
            }
        }
        
        // Remove keys
        forAll(keysToRemove, i)
        {
            history_.erase(keysToRemove[i]);
        }
    }

    void contactForce::computeRelativeVelocity(
        const particle& pA,
        const particle& pB,
        const Foam::point& contactPoint,
        const Foam::vector& normal,
        Foam::vector& vRelN,
        Foam::vector& vRelT
    ) const
    {
        Foam::vector rA = contactPoint - pA.position();
        Foam::vector rB = contactPoint - pB.position();
        
        Foam::vector vA = pA.velocity() + (pA.angularVelocity() ^ rA);
        Foam::vector vB = pB.velocity() + (pB.angularVelocity() ^ rB);
        
        Foam::vector vRel = vA - vB;
        
        Foam::scalar vn = vRel & normal;
        vRelN = vn * normal;
        vRelT = vRel - vRelN;
    }

    void contactForce::computeRelativeVelocityParticleWall(
        const particle& p,
        const wall& w,
        const Foam::point& contactPoint,
        const Foam::vector& normal,
        Foam::vector& vRelN,
        Foam::vector& vRelT
    ) const
    {
        Foam::vector r = contactPoint - p.position();
        Foam::vector v = p.velocity() + (p.angularVelocity() ^ r);
        
        // Wall velocity is assumed zero for now
        Foam::vector vRel = v;
        
        Foam::scalar vn = vRel & normal;
        vRelN = vn * normal;
        vRelT = vRel - vRelN;
    }

    void contactForce::computeForce(
        particle& pA,
        particle& pB,
        const ContactInfo& contact,
        const Foam::scalar dt
    ) const
    {
        if (!contact.isContact || contact.contactPoints.empty()) return;

        // Calculate centroid of contact points
        Foam::point contactPoint = Foam::vector::zero;
        forAll(contact.contactPoints, i)
        {
            contactPoint += contact.contactPoints[i];
        }
        contactPoint /= contact.contactPoints.size();
        
        Foam::vector normal = contact.normal; // Points from B to A
        Foam::scalar depth = contact.depth;
        
        // History Management
        Foam::label idA = pA.id();
        Foam::label idB = pB.id();
        Foam::string key = idA < idB ? Foam::name(idA) + "_" + Foam::name(idB) : Foam::name(idB) + "_" + Foam::name(idA);
        
        if (!history_.found(key))
        {
            history_.insert(key, ContactHistory());
        }
        ContactHistory& hist = history_[key];
        
        // Hysteretic Normal Force (Walton-Braun)
        Foam::scalar FnMag = 0.0;
        
        if (depth >= hist.maxOverlap)
        {
            // Loading
            hist.maxOverlap = depth;
            hist.isLoading = true;
            FnMag = kn_ * depth;
        }
        else
        {
            // Unloading / Reloading
            hist.isLoading = false;
            Foam::scalar delta0 = hist.maxOverlap * (1.0 - kn_ / knUnloading_);
            FnMag = knUnloading_ * (depth - delta0);
        }
        
        if (FnMag < 0) FnMag = 0;
        
        // Normal points from A to B (EPA convention)
        // Force on A should be opposite to normal (repulsive)
        Foam::vector Fn = -FnMag * normal;
        
        // Damping (Viscous)
        Foam::vector vRelN, vRelT;
        computeRelativeVelocity(pA, pB, contactPoint, normal, vRelN, vRelT);
        
        Fn -= gammaN_ * vRelN; // Opposes relative normal velocity
        
        // Tangential Force (Spring + Coulomb)
        // Update spring
        hist.tangentialSpring -= vRelT * dt; // Displacement integration
        
        Foam::vector Ft = kt_ * hist.tangentialSpring;
        Ft -= gammaT_ * vRelT; // Damping
        
        // Coulomb Limit
        Foam::scalar FnTotal = Fn & normal;
        if (FnTotal < 0) FnTotal = 0;
        Foam::scalar CoulombLimit = muC_ * FnTotal;
        
        Foam::scalar FtMag = Foam::mag(Ft);
        if (FtMag > CoulombLimit)
        {
            if (FtMag > Foam::SMALL)
            {
                Ft *= (CoulombLimit / FtMag);
                // Adjust spring to match limited force (ignoring damping for simplicity)
                // hist.tangentialSpring = Ft / kt_; 
            }
            else
            {
                Ft = Foam::vector::zero;
            }
        }
        
        Foam::vector totalForce = Fn + Ft;
        
        // Apply to A
        pA.updateContactForces(forces(totalForce, (contactPoint - pA.position()) ^ totalForce));
        
        // Apply to B (Newton's 3rd Law)
        pB.updateContactForces(forces(-totalForce, (contactPoint - pB.position()) ^ (-totalForce)));
    }

    void contactForce::computeForceParticleWall(
        particle& p,
        const wall& w,
        const ContactInfo& contact,
        const Foam::scalar dt
    ) const
    {
        if (!contact.isContact || contact.contactPoints.empty()) return;

        Foam::point contactPoint = Foam::vector::zero;
        forAll(contact.contactPoints, i)
        {
            contactPoint += contact.contactPoints[i];
        }
        contactPoint /= contact.contactPoints.size();
        
        Foam::vector normal = contact.normal; // Points from Wall to Particle
        Foam::scalar depth = contact.depth;
        
        // History
        Foam::label id = p.id();
        Foam::string key = "W_" + w.name() + "_" + Foam::name(id);
        
        if (!history_.found(key))
        {
            history_.insert(key, ContactHistory());
        }
        ContactHistory& hist = history_[key];
        
        // Normal Force
        Foam::scalar FnMag = 0.0;
        if (depth >= hist.maxOverlap)
        {
            hist.maxOverlap = depth;
            hist.isLoading = true;
            FnMag = kn_ * depth;
        }
        else
        {
            hist.isLoading = false;
            Foam::scalar delta0 = hist.maxOverlap * (1.0 - kn_ / knUnloading_);
            FnMag = knUnloading_ * (depth - delta0);
        }
        if (FnMag < 0) FnMag = 0;
        
        // Normal points from P to W
        // Force on P should be opposite to normal
        Foam::vector Fn = -FnMag * normal;
        
        // Damping
        Foam::vector vRelN, vRelT;
        computeRelativeVelocityParticleWall(p, w, contactPoint, normal, vRelN, vRelT);
        
        Fn -= gammaN_ * vRelN;
        
        // Tangential
        hist.tangentialSpring -= vRelT * dt;
        Foam::vector Ft = kt_ * hist.tangentialSpring;
        Ft -= gammaT_ * vRelT;
        
        Foam::scalar FnTotal = Fn & normal;
        if (FnTotal < 0) FnTotal = 0;
        Foam::scalar CoulombLimit = muC_ * FnTotal;
        
        if (Foam::mag(Ft) > CoulombLimit)
        {
             if (Foam::mag(Ft) > Foam::SMALL)
            {
                Ft *= (CoulombLimit / Foam::mag(Ft));
            }
            else
            {
                Ft = Foam::vector::zero;
            }
        }
        
        Foam::vector totalForce = Fn + Ft;
        
        p.updateContactForces(forces(totalForce, (contactPoint - p.position()) ^ totalForce));
    }
}
