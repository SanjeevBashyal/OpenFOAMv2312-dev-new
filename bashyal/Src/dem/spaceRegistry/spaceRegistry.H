#ifndef SPACE_REGISTRY_H
#define SPACE_REGISTRY_H

#include "particle.H"
#include "boundBox.H"
#include <vector>
#include <map>
#include <set>
#include <utility>

namespace Bashyal {

/**
 * @class spaceRegistry
 * @brief Tracks particle positions in space and identifies potential collision pairs.
 * 
 * This class maintains a spatial registry of all particles and provides
 * efficient collision detection by using bounding box checks and spatial
 * partitioning techniques.
 */
class spaceRegistry {
private:
    // Particle registry
    std::map<std::string, particle*> particles_;
    
    // Spatial data structures for collision detection
    struct ParticleInfo {
        particle* p;
        Foam::boundBox bbox;
        Foam::point predictedPosition;
        Foam::boundBox predictedBbox;
        Foam::scalar collisionRadius;
    };
    
    std::map<std::string, ParticleInfo> particleInfo_;
    
    // Collision pairs (set of pairs of particle names)
    std::set<std::pair<std::string, std::string>> collisionPairs_;
    
    // Spatial partitioning (simple grid-based approach)
    struct SpatialCell {
        std::vector<std::string> particleNames;
    };
    
    Foam::vector cellSize_;
    Foam::boundBox domain_;
    std::map<Foam::label, SpatialCell> spatialGrid_;
    
    // Configuration parameters
    Foam::scalar safetyFactor_;  // Multiplier for collision detection radius
    Foam::scalar timeStep_;      // Current time step for prediction
    
public:
    /**
     * @brief Constructor
     * @param domain The bounding box of the simulation domain
     * @param cellSize The size of spatial partitioning cells
     * @param safetyFactor Safety factor for collision detection (default: 1.2)
     */
    spaceRegistry(const Foam::boundBox& domain, 
                  const Foam::vector& cellSize,
                  Foam::scalar safetyFactor = 1.2);
    
    /**
     * @brief Destructor
     */
    ~spaceRegistry();
    
    /**
     * @brief Add a particle to the spatial registry
     * @param p Pointer to the particle
     * @param name Unique name for the particle
     */
    void addParticle(particle* p, const std::string& name);
    
    /**
     * @brief Remove a particle from the spatial registry
     * @param name Name of the particle to remove
     */
    void removeParticle(const std::string& name);
    
    /**
     * @brief Update particle information and rebuild spatial structures
     * @param dt Time step for position prediction
     */
    void updateParticles(Foam::scalar dt);
    
    /**
     * @brief Find potential collision pairs for the next time step
     * @return Set of particle name pairs that might collide
     */
    std::set<std::pair<std::string, std::string>> findCollisionPairs();
    
    /**
     * @brief Get all registered particles
     * @return Map of particle names to particle pointers
     */
    const std::map<std::string, particle*>& getParticles() const { return particles_; }
    
    /**
     * @brief Get particle by name
     * @param name Name of the particle
     * @return Pointer to the particle, or nullptr if not found
     */
    particle* getParticle(const std::string& name) const;
    
    /**
     * @brief Get the number of registered particles
     * @return Number of particles
     */
    size_t size() const { return particles_.size(); }
    
    /**
     * @brief Clear all particles from the registry
     */
    void clear();
    
    /**
     * @brief Set the time step for collision prediction
     * @param dt Time step
     */
    void setTimeStep(Foam::scalar dt) { timeStep_ = dt; }
    
    /**
     * @brief Set the safety factor for collision detection
     * @param factor Safety factor
     */
    void setSafetyFactor(Foam::scalar factor) { safetyFactor_ = factor; }
    
private:
    /**
     * @brief Calculate bounding box for a particle
     * @param p Pointer to the particle
     * @return Bounding box of the particle
     */
    Foam::boundBox calculateBoundingBox(const particle* p) const;
    
    /**
     * @brief Predict particle position for the next time step
     * @param p Pointer to the particle
     * @param dt Time step
     * @return Predicted position
     */
    Foam::point predictPosition(const particle* p, Foam::scalar dt) const;
    
    /**
     * @brief Calculate collision radius for a particle
     * @param p Pointer to the particle
     * @return Collision radius
     */
    Foam::scalar calculateCollisionRadius(const particle* p) const;
    
    /**
     * @brief Update spatial grid with current particle positions
     */
    void updateSpatialGrid();
    
    /**
     * @brief Get grid cell index for a position
     * @param pos Position in space
     * @return Grid cell index
     */
    Foam::label getGridIndex(const Foam::point& pos) const;
    
    /**
     * @brief Check if two particles might collide based on their bounding boxes
     * @param p1 First particle info
     * @param p2 Second particle info
     * @return True if collision is possible
     */
    bool checkCollisionPossibility(const ParticleInfo& p1, const ParticleInfo& p2) const;
    
    /**
     * @brief Get neighboring grid cells for a given cell
     * @param cellIndex Grid cell index
     * @return Vector of neighboring cell indices
     */
    std::vector<Foam::label> getNeighboringCells(Foam::label cellIndex) const;
};

} // namespace Bashyal

#endif // SPACE_REGISTRY_H 