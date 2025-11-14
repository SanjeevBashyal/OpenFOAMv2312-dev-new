#include "spaceRegistry.H"
#include "constants.H"
#include <algorithm>
#include <cmath>

namespace Bashyal {

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

spaceRegistry::spaceRegistry(const Foam::boundBox& domain, 
                           const Foam::vector& cellSize,
                           Foam::scalar safetyFactor)
    : domain_(domain),
      cellSize_(cellSize),
      safetyFactor_(safetyFactor),
      timeStep_(0.0)
{
    // Validate input parameters
    if (safetyFactor <= 0.0)
    {
        FatalErrorInFunction
            << "Safety factor must be positive, got: " << safetyFactor
            << Foam::abort(Foam::FatalError);
    }
    
    if (cellSize.x() <= 0.0 || cellSize.y() <= 0.0 || cellSize.z() <= 0.0)
    {
        FatalErrorInFunction
            << "Cell size must be positive in all dimensions, got: " << cellSize
            << Foam::abort(Foam::FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

spaceRegistry::~spaceRegistry()
{
    clear();
}

// * * * * * * * * * * * * * * * * Methods * * * * * * * * * * * * * * * //

void spaceRegistry::addParticle(particle* p, const std::string& name)
{
    if (!p)
    {
        FatalErrorInFunction
            << "Cannot add null particle pointer"
            << Foam::abort(Foam::FatalError);
    }
    
    if (particles_.find(name) != particles_.end())
    {
        FatalErrorInFunction
            << "Particle with name '" << name << "' already exists"
            << Foam::abort(Foam::FatalError);
    }
    
    // Add to particle registry
    particles_[name] = p;
    
    // Initialize particle info
    ParticleInfo info;
    info.p = p;
    info.bbox = calculateBoundingBox(p);
    info.predictedPosition = p->position();
    info.predictedBbox = info.bbox;
    info.collisionRadius = calculateCollisionRadius(p);
    
    particleInfo_[name] = info;
}

void spaceRegistry::removeParticle(const std::string& name)
{
    auto it = particles_.find(name);
    if (it != particles_.end())
    {
        particles_.erase(it);
        particleInfo_.erase(name);
        
        // Remove from spatial grid
        updateSpatialGrid();
    }
}

void spaceRegistry::updateParticles(Foam::scalar dt)
{
    timeStep_ = dt;
    
    // Update particle information
    for (auto& pair : particleInfo_)
    {
        const std::string& name = pair.first;
        ParticleInfo& info = pair.second;
        particle* p = info.p;
        
        // Update current bounding box
        info.bbox = calculateBoundingBox(p);
        
        // Predict position for next time step
        info.predictedPosition = predictPosition(p, dt);
        
        // Calculate predicted bounding box
        Foam::pointField globalPoints = p->getGlobalPoints();
        Foam::vector velocity = p->velocity();
        Foam::vector angularVelocity = p->angularVelocity();
        
        // Simple prediction: move all points by velocity * dt
        Foam::pointField predictedPoints(globalPoints.size());
        forAll(globalPoints, i)
        {
            predictedPoints[i] = globalPoints[i] + velocity * dt;
        }
        
        info.predictedBbox = Foam::boundBox(predictedPoints);
        info.collisionRadius = calculateCollisionRadius(p);
    }
    
    // Update spatial grid
    updateSpatialGrid();
}

std::set<std::pair<std::string, std::string>> spaceRegistry::findCollisionPairs()
{
    collisionPairs_.clear();
    
    // Get all particle names
    std::vector<std::string> particleNames;
    for (const auto& pair : particles_)
    {
        particleNames.push_back(pair.first);
    }
    
    // Check all pairs of particles
    for (size_t i = 0; i < particleNames.size(); ++i)
    {
        for (size_t j = i + 1; j < particleNames.size(); ++j)
        {
            const std::string& name1 = particleNames[i];
            const std::string& name2 = particleNames[j];
            
            const ParticleInfo& info1 = particleInfo_.at(name1);
            const ParticleInfo& info2 = particleInfo_.at(name2);
            
            if (checkCollisionPossibility(info1, info2))
            {
                // Add to collision pairs (ensure consistent ordering)
                if (name1 < name2)
                {
                    collisionPairs_.insert({name1, name2});
                }
                else
                {
                    collisionPairs_.insert({name2, name1});
                }
            }
        }
    }
    
    return collisionPairs_;
}

particle* spaceRegistry::getParticle(const std::string& name) const
{
    auto it = particles_.find(name);
    return (it != particles_.end()) ? it->second : nullptr;
}

void spaceRegistry::clear()
{
    particles_.clear();
    particleInfo_.clear();
    collisionPairs_.clear();
    spatialGrid_.clear();
}

// * * * * * * * * * * * * * * * * Private Methods * * * * * * * * * * * * * * //

Foam::boundBox spaceRegistry::calculateBoundingBox(const particle* p) const
{
    if (!p)
    {
        return Foam::boundBox();
    }
    
    Foam::pointField globalPoints = p->getGlobalPoints();
    if (globalPoints.empty())
    {
        return Foam::boundBox();
    }
    
    return Foam::boundBox(globalPoints);
}

Foam::point spaceRegistry::predictPosition(const particle* p, Foam::scalar dt) const
{
    if (!p)
    {
        return Foam::point::zero;
    }
    
    // Simple prediction: current position + velocity * dt
    return p->position() + p->velocity() * dt;
}

Foam::scalar spaceRegistry::calculateCollisionRadius(const particle* p) const
{
    if (!p)
    {
        return 0.0;
    }
    
    Foam::boundBox bbox = calculateBoundingBox(p);
    Foam::vector size = bbox.span();
    
    // Use the maximum dimension as collision radius
    Foam::scalar maxDim = Foam::max(Foam::max(size.x(), size.y()), size.z());
    
    // Apply safety factor
    return maxDim * safetyFactor_ * 0.5;
}

void spaceRegistry::updateSpatialGrid()
{
    // Clear current grid
    spatialGrid_.clear();
    
    // Add particles to grid cells
    for (const auto& pair : particleInfo_)
    {
        const std::string& name = pair.first;
        const ParticleInfo& info = pair.second;
        
        // Get grid cell for current position
        Foam::label cellIndex = getGridIndex(info.p->position());
        
        // Add to spatial grid
        spatialGrid_[cellIndex].particleNames.push_back(name);
    }
}

Foam::label spaceRegistry::getGridIndex(const Foam::point& pos) const
{
    // Calculate grid coordinates
    Foam::vector relativePos = pos - domain_.min();
    Foam::label ix = static_cast<Foam::label>(relativePos.x() / cellSize_.x());
    Foam::label iy = static_cast<Foam::label>(relativePos.y() / cellSize_.y());
    Foam::label iz = static_cast<Foam::label>(relativePos.z() / cellSize_.z());
    
    // Clamp to domain bounds
    Foam::label nx = static_cast<Foam::label>((domain_.max().x() - domain_.min().x()) / cellSize_.x());
    Foam::label ny = static_cast<Foam::label>((domain_.max().y() - domain_.min().y()) / cellSize_.y());
    
    ix = Foam::max(0, Foam::min(ix, nx - 1));
    iy = Foam::max(0, Foam::min(iy, ny - 1));
    iz = Foam::max(0, Foam::min(iz, static_cast<Foam::label>((domain_.max().z() - domain_.min().z()) / cellSize_.z()) - 1));
    
    // Convert to linear index
    return iz * nx * ny + iy * nx + ix;
}

bool spaceRegistry::checkCollisionPossibility(const ParticleInfo& p1, const ParticleInfo& p2) const
{
    // Check if predicted bounding boxes overlap
    if (!p1.predictedBbox.overlaps(p2.predictedBbox))
    {
        return false;
    }
    
    // Additional check: distance between centers
    Foam::point center1 = p1.predictedBbox.midpoint();
    Foam::point center2 = p2.predictedBbox.midpoint();
    Foam::scalar distance = Foam::mag(center2 - center1);
    Foam::scalar minDistance = p1.collisionRadius + p2.collisionRadius;
    
    return distance <= minDistance;
}

std::vector<Foam::label> spaceRegistry::getNeighboringCells(Foam::label cellIndex) const
{
    std::vector<Foam::label> neighbors;
    
    // Calculate grid dimensions
    Foam::label nx = static_cast<Foam::label>((domain_.max().x() - domain_.min().x()) / cellSize_.x());
    Foam::label ny = static_cast<Foam::label>((domain_.max().y() - domain_.min().y()) / cellSize_.y());
    Foam::label nz = static_cast<Foam::label>((domain_.max().z() - domain_.min().z()) / cellSize_.z());
    
    // Convert linear index to 3D coordinates
    Foam::label iz = cellIndex / (nx * ny);
    Foam::label iy = (cellIndex % (nx * ny)) / nx;
    Foam::label ix = cellIndex % nx;
    
    // Check all 26 neighboring cells (3x3x3 - 1)
    for (Foam::label dz = -1; dz <= 1; ++dz)
    {
        for (Foam::label dy = -1; dy <= 1; ++dy)
        {
            for (Foam::label dx = -1; dx <= 1; ++dx)
            {
                if (dx == 0 && dy == 0 && dz == 0) continue; // Skip self
                
                Foam::label newIx = ix + dx;
                Foam::label newIy = iy + dy;
                Foam::label newIz = iz + dz;
                
                // Check bounds
                if (newIx >= 0 && newIx < nx && 
                    newIy >= 0 && newIy < ny && 
                    newIz >= 0 && newIz < nz)
                {
                    Foam::label neighborIndex = newIz * nx * ny + newIy * nx + newIx;
                    neighbors.push_back(neighborIndex);
                }
            }
        }
    }
    
    return neighbors;
}

} // namespace Bashyal 