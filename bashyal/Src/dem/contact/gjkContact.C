#include "gjkContact.H"
#include "labelPair.H"
#include "DynamicList.H" // Added as per instruction
#include "quaternion.H"
#include "tensor.H"
#include "constants.H"
#include "scalar.H"
#include <algorithm>
#include <limits>

namespace Bashyal
{
    // * * * * * * * * * * * * * * * * Static Members * * * * * * * * * * * * * * * //

    const Foam::scalar gjkContact::tolerance_ = 1e-6;

    // * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

    gjkContact::gjkContact()
    {}

    gjkContact::~gjkContact()
    {}

    // * * * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * * //

    Foam::point gjkContact::getParticleSupport(
        const particle& p,
        const Foam::vector& d
    ) const
    {
        return p.getSupport(d);
    }

    Foam::point gjkContact::getWallSupport(
        const wall& w,
        const Foam::vector& d
    ) const
    {
        Foam::vector dNorm = d / (Foam::mag(d) + Foam::SMALL);
        
        switch (w.type())
        {
            case WallType::PLANE:
            {
                // For a plane, support point is at infinity in direction d
                // But we need a finite point, so use a large distance
                Foam::point pos = w.position();
                Foam::vector normal = w.normal();
                
                // If direction points away from plane, support is at infinity
                if ((dNorm & normal) > 0)
                {
                    return pos + 1e6 * dNorm;  // Large distance
                }
                else
                {
                    // If direction points into plane, support is on the plane
                    // But strictly speaking, a plane is infinite.
                    // For GJK against a plane, we usually project the point onto the plane.
                    // Here we return a point far away for now, but GJK with infinite plane is special.
                    // Better approach: Use the closest point on plane logic directly for plane-particle.
                    return pos - 1e6 * dNorm;
                }
            }
            case WallType::CYLINDER:
            {
                Foam::point center = w.position();
                Foam::vector axis = w.normal();
                Foam::scalar radius = w.radius();
                
                // Project direction onto plane perpendicular to axis
                Foam::vector dPerp = dNorm - (dNorm & axis) * axis;
                Foam::scalar magDPerp = Foam::mag(dPerp);
                
                if (magDPerp > Foam::SMALL)
                {
                    dPerp /= magDPerp;
                    // Support point is on cylinder surface in direction dPerp
                    Foam::point support = center + radius * dPerp;
                    
                    // Height limits
                    Foam::scalar height = w.height();
                    Foam::scalar zMin = center.z() - 0.5 * height;
                    Foam::scalar zMax = center.z() + 0.5 * height;
                    
                    // For the "support" in the axial direction, we should check the sign of d & axis
                    if ((dNorm & axis) > 0) support += 0.5 * height * axis;
                    else support -= 0.5 * height * axis;

                    return support;
                }
                else
                {
                    // Direction is parallel to axis
                    if ((dNorm & axis) > 0) return center + (0.5 * w.height() * axis);
                    else return center - (0.5 * w.height() * axis);
                }
            }
            case WallType::SPHERE:
            {
                Foam::point center = w.position();
                Foam::scalar radius = w.radius();
                return center + radius * dNorm;
            }
            default:
                return w.position();
        }
    }

    Foam::point gjkContact::support(
        const particle& pA,
        const particle& pB,
        const Foam::vector& d
    ) const
    {
        return getParticleSupport(pA, d) - getParticleSupport(pB, -d);
    }

    Foam::point gjkContact::supportParticleWall(
        const particle& p,
        const wall& w,
        const Foam::vector& d
    ) const
    {
        return getParticleSupport(p, d) - getWallSupport(w, -d);
    }

    // Helper to check if two vectors are in the same direction
    bool sameDirection(const Foam::vector& d, const Foam::vector& a)
    {
        return (d & a) > 0;
    }

    bool gjkContact::doSimplex(
        Simplex& simplex,
        Foam::vector& direction
    ) const
    {
        // A is the last added point
        Foam::point a = simplex[simplex.size() - 1];
        Foam::vector ao = -a; // Vector from A to Origin

        if (simplex.size() == 4) // Tetrahedron
        {
            Foam::point b = simplex[2];
            Foam::point c = simplex[1];
            Foam::point d = simplex[0];

            Foam::vector ab = b - a;
            Foam::vector ac = c - a;
            Foam::vector ad = d - a;
            
            Foam::vector abc = ab ^ ac;
            Foam::vector acd = ac ^ ad;
            Foam::vector adb = ad ^ ab;

            // Check faces
            if (sameDirection(abc, ao))
            {
                // In front of ABC
                simplex.clear();
                simplex.add(c); simplex.add(b); simplex.add(a); // Keep ABC
                // But wait, we need to check edges too? 
                // Actually, for GJK, if we are outside a face, we reduce to Triangle case
                return doSimplex(simplex, direction);
            }
            else if (sameDirection(acd, ao))
            {
                simplex.clear();
                simplex.add(d); simplex.add(c); simplex.add(a);
                return doSimplex(simplex, direction);
            }
            else if (sameDirection(adb, ao))
            {
                simplex.clear();
                simplex.add(b); simplex.add(d); simplex.add(a);
                return doSimplex(simplex, direction);
            }
            else
            {
                // Inside tetrahedron!
                return true;
            }
        }
        else if (simplex.size() == 3) // Triangle
        {
            Foam::point b = simplex[1];
            Foam::point c = simplex[0];

            Foam::vector ab = b - a;
            Foam::vector ac = c - a;
            Foam::vector abc = ab ^ ac;

            if (sameDirection(abc ^ ac, ao))
            {
                if (sameDirection(ac, ao))
                {
                    simplex.clear();
                    simplex.add(c); simplex.add(a);
                    direction = (ac ^ ao) ^ ac;
                }
                else
                {
                    // Star case (AB)
                    if (sameDirection(ab, ao))
                    {
                        simplex.clear();
                        simplex.add(b); simplex.add(a);
                        direction = (ab ^ ao) ^ ab;
                    }
                    else
                    {
                        simplex.clear();
                        simplex.add(a);
                        direction = ao;
                    }
                }
            }
            else
            {
                if (sameDirection(ab ^ abc, ao))
                {
                    if (sameDirection(ab, ao))
                    {
                        simplex.clear();
                        simplex.add(b); simplex.add(a);
                        direction = (ab ^ ao) ^ ab;
                    }
                    else
                    {
                        simplex.clear();
                        simplex.add(a);
                        direction = ao;
                    }
                }
                else
                {
                    if (sameDirection(abc, ao))
                    {
                        direction = abc;
                    }
                    else
                    {
                        // Winding order fix
                        Simplex newSimplex;
                        newSimplex.add(b); newSimplex.add(c); newSimplex.add(a);
                        simplex = newSimplex;
                        direction = -abc;
                    }
                }
            }
        }
        else if (simplex.size() == 2) // Line
        {
            Foam::point b = simplex[0];
            Foam::vector ab = b - a;
            
            if (sameDirection(ab, ao))
            {
                Foam::vector cross = ab ^ ao;
                if (Foam::mag(cross) < Foam::SMALL)
                {
                    // Origin is on the line segment AB (degenerate case)
                    // Expand in a direction perpendicular to AB
                    Foam::vector perp = ab ^ Foam::vector(1, 0, 0);
                    if (Foam::mag(perp) < Foam::SMALL)
                    {
                        perp = ab ^ Foam::vector(0, 1, 0);
                    }
                    direction = perp;
                }
                else
                {
                    direction = cross ^ ab;
                }
            }
            else
            {
                simplex.clear();
                simplex.add(a);
                direction = ao;
            }
        }
        
        return false;
    }

    // * * * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * * * //

    ContactInfo gjkContact::detectContact(
        const particle& pA,
        const particle& pB,
        const Foam::scalar crustWidthA,
        const Foam::scalar crustWidthB
    ) const
    {
        ContactInfo contact;
        
        // Initial direction
        Foam::vector d = pB.position() - pA.position();
        if (Foam::mag(d) < Foam::SMALL) d = Foam::vector(1, 0, 0);
        
        Simplex simplex;
        simplex.add(support(pA, pB, d));
        
        d = -simplex[0]; // Towards origin
        
        bool collision = false;
        
        for (int i = 0; i < maxIterations_; ++i)
        {
            Foam::point p = support(pA, pB, d);
            
            if ((p & d) < 0)
            {
                // No collision
                collision = false;
                break;
            }
            
            simplex.add(p);
            
            if (doSimplex(simplex, d))
            {
                collision = true;
                break;
            }
        }
        
        if (collision)
        {
            contact.isContact = true;
            
            // EPA to find depth and normal
            Foam::scalar depth;
            Foam::vector normal;
            if (epa(simplex, pA, pB, depth, normal))
            {
                contact.depth = depth;
                contact.normal = normal;
                
                // Clipping to find manifold
                generateManifold(pA, pB, normal, contact.contactPoints);
                
                // Fallback if manifold generation fails
                if (contact.contactPoints.empty())
                {
                    Foam::point pA_s = pA.getSupport(-normal);
                    Foam::point pB_s = pB.getSupport(normal);
                    contact.contactPoints.append(0.5 * (pA_s + pB_s));
                }

                // Average contact point for debug/legacy
                contact.contactPoint = Foam::vector::zero;
                if (contact.contactPoints.size() > 0)
                {
                    for (const auto& pt : contact.contactPoints)
                    {
                        contact.contactPoint += pt;
                    }
                    contact.contactPoint /= contact.contactPoints.size();
                }
            }
            else
            {
                // EPA failed, fallback?
                contact.isContact = false;
            }
        }
        
        return contact;
    }

    ContactInfo gjkContact::detectContactParticleWall(
        const particle& p,
        const wall& w,
        const Foam::scalar crustWidthP
    ) const
    {
        ContactInfo contact;
        
        if (w.type() == WallType::PLANE)
        {
            Foam::vector n = w.normal();
            Foam::point p0 = w.position();
            
            // Find the point on particle furthest in direction -n (into the wall)
            Foam::point deepestPt = p.getSupport(-n);
            
            // Distance from plane: (pt - p0) . n
            Foam::scalar dist = (deepestPt - p0) & n;
            
            if (dist < crustWidthP)
            {
                contact.isContact = true;
                contact.depth = -dist; // Penetration depth
                contact.normal = n;    // Normal points out of wall
                contact.contactPoint = deepestPt;
                contact.contactPoints.append(deepestPt);
                
                // TODO: For face-face contact, we should find all vertices within tolerance
            }
            return contact;
        }

        // Fallback to GJK for other wall types (Cylinder, Sphere)
        Foam::vector d = w.position() - p.position(); // Rough direction
        if (Foam::mag(d) < Foam::SMALL) d = Foam::vector(0, 0, 1);
        
        Simplex simplex;
        simplex.add(supportParticleWall(p, w, d));
        d = -simplex[0];
        
        bool collision = false;
        
        for (int i = 0; i < maxIterations_; ++i)
        {
            Foam::point pt = supportParticleWall(p, w, d);
            if ((pt & d) < 0)
            {
                collision = false;
                break;
            }
            simplex.add(pt);
            if (doSimplex(simplex, d))
            {
                collision = true;
                break;
            }
        }
        
        if (collision)
        {
            contact.isContact = true;
            Foam::scalar depth;
            Foam::vector normal;
            if (epaParticleWall(simplex, p, w, depth, normal))
            {
                contact.depth = depth;
                contact.normal = normal;
                generateManifoldParticleWall(p, w, normal, contact.contactPoints);
                 if (contact.contactPoints.size() > 0)
                {
                    contact.contactPoint = Foam::vector::zero;
                    for (const auto& pt : contact.contactPoints)
                    {
                        contact.contactPoint += pt;
                    }
                    contact.contactPoint /= contact.contactPoints.size();
                }
            }
        }
        
        return contact;
    }
    
    // EPA Implementation
    struct EPAFace
    {
        Foam::label a, b, c;
        Foam::vector normal;
        Foam::scalar distance;
        bool valid;
        
        EPAFace() : a(-1), b(-1), c(-1), distance(Foam::GREAT), valid(false) {}

        EPAFace(Foam::label _a, Foam::label _b, Foam::label _c, const Foam::List<Foam::point>& points)
        : a(_a), b(_b), c(_c), valid(true)
        {
            Foam::vector ab = points[b] - points[a];
            Foam::vector ac = points[c] - points[a];
            normal = ab ^ ac;
            Foam::scalar magN = Foam::mag(normal);
            if (magN > Foam::SMALL)
            {
                normal /= magN;
                distance = normal & points[a];
                
                // Ensure normal points away from origin
                if (distance < 0)
                {
                    normal = -normal;
                    distance = -distance;
                    // Swap b and c to maintain winding
                    std::swap(b, c);
                }
            }
            else
            {
                valid = false;
                distance = Foam::GREAT;
            }
        }
    };

    bool gjkContact::epa(
        const Simplex& simplex,
        const particle& pA,
        const particle& pB,
        Foam::scalar& depth,
        Foam::vector& normal
    ) const
    {
        Foam::DynamicList<Foam::point> points(simplex.points);
        Foam::DynamicList<EPAFace> faces;
        
        // Initial tetrahedron faces
        faces.append(EPAFace(0, 1, 2, points));
        faces.append(EPAFace(0, 2, 3, points));
        faces.append(EPAFace(0, 3, 1, points));
        faces.append(EPAFace(1, 3, 2, points));
        
        for (int iter = 0; iter < maxEpaIterations_; ++iter)
        {
            // Find closest face
            Foam::scalar minDist = Foam::GREAT;
            Foam::label minFaceIdx = -1;
            
            forAll(faces, i)
            {
                if (faces[i].valid && faces[i].distance < minDist)
                {
                    minDist = faces[i].distance;
                    minFaceIdx = i;
                }
            }
            
            if (minFaceIdx == -1) break;
            
            Foam::vector n = faces[minFaceIdx].normal;
            Foam::point p = support(pA, pB, n);
            
            Foam::scalar d = p & n;
            
            if (d - minDist < tolerance_)
            {
                // Convergence
                depth = d;
                normal = n;
                return true;
            }
            
            // Expand polytope
            // Remove faces visible from p
            Foam::DynamicList<Foam::labelPair> edges;
            
            forAll(faces, i)
            {
                if (!faces[i].valid) continue;
                
                if ((faces[i].normal & (p - points[faces[i].a])) > 0)
                {
                    // Face is visible, remove it
                    faces[i].valid = false;
                    
                    // Add edges to list
                    // If edge exists, remove it (shared edge between removed faces)
                    // Else add it
                    
                    auto addEdge = [&](Foam::label a, Foam::label b)
                    {
                        bool found = false;
                        forAll(edges, k)
                        {
                            if ((edges[k].first() == b && edges[k].second() == a) || 
                                (edges[k].first() == a && edges[k].second() == b))
                            {
                                // Edge shared with another removed face -> internal edge -> remove
                                edges[k] = edges.last();
                                edges.remove(edges.size()-1);
                                found = true;
                                break;
                            }
                        }
                        if (!found)
                        {
                            edges.append(Foam::labelPair(a, b));
                        }
                    };
                    
                    addEdge(faces[i].a, faces[i].b);
                    addEdge(faces[i].b, faces[i].c);
                    addEdge(faces[i].c, faces[i].a);
                }
            }
            
            // Add new faces
            points.append(p);
            Foam::label pIdx = points.size() - 1;
            
            for (const auto& edge : edges)
            {
                faces.append(EPAFace(edge.first(), edge.second(), pIdx, points));
            }
        }
        
        return false;
    }
    
    bool gjkContact::epaParticleWall(
        const Simplex& simplex,
        const particle& p,
        const wall& w,
        Foam::scalar& depth,
        Foam::vector& normal
    ) const
    {
        // Similar to EPA but using supportParticleWall
        Foam::List<Foam::point> points = simplex.points;
        Foam::List<EPAFace> faces;
        
        faces.append(EPAFace(0, 1, 2, points));
        faces.append(EPAFace(0, 2, 3, points));
        faces.append(EPAFace(0, 3, 1, points));
        faces.append(EPAFace(1, 3, 2, points));
        
        for (int iter = 0; iter < maxEpaIterations_; ++iter)
        {
            Foam::scalar minDist = Foam::GREAT;
            Foam::label minFaceIdx = -1;
            
            forAll(faces, i)
            {
                if (faces[i].valid && faces[i].distance < minDist)
                {
                    minDist = faces[i].distance;
                    minFaceIdx = i;
                }
            }
            
            if (minFaceIdx == -1) break;
            
            Foam::vector n = faces[minFaceIdx].normal;
            Foam::point pt = supportParticleWall(p, w, n);
            
            Foam::scalar d = pt & n;
            
            if (d - minDist < tolerance_)
            {
                depth = d;
                normal = n;
                return true;
            }
            
            Foam::DynamicList<Foam::labelPair> edges;
            forAll(faces, i)
            {
                if (!faces[i].valid) continue;
                if ((faces[i].normal & (pt - points[faces[i].a])) > 0)
                {
                    faces[i].valid = false;
                    auto addEdge = [&](Foam::label a, Foam::label b)
                    {
                        bool found = false;
                        forAll(edges, k)
                        {
                            if ((edges[k].first() == b && edges[k].second() == a) || 
                                (edges[k].first() == a && edges[k].second() == b))
                            {
                                edges[k] = edges.last();
                                edges.remove(edges.size()-1);
                                found = true;
                                break;
                            }
                        }
                        if (!found) edges.append(Foam::labelPair(a, b));
                    };
                    addEdge(faces[i].a, faces[i].b);
                    addEdge(faces[i].b, faces[i].c);
                    addEdge(faces[i].c, faces[i].a);
                }
            }
            
            points.append(pt);
            Foam::label pIdx = points.size() - 1;
            for (const auto& edge : edges)
            {
                faces.append(EPAFace(edge.first(), edge.second(), pIdx, points));
            }
        }
        return false;
    }
    
    // Clipping Implementation
    
    struct ClipVertex
    {
        Foam::point v;
        // Add feature ID if needed for caching
    };

    // Helper to find the best face (most parallel to normal)
    void getBestFace(
        const particle& p,
        const Foam::vector& normal,
        Foam::List<Foam::point>& facePoints,
        Foam::vector& faceNormal
    )
    {
        // Transform normal to local frame
        Foam::vector localNormal = p.orientation().R().T() & normal;
        
        const Foam::faceList& faces = p.faces();
        const Foam::pointField& points = p.vertices();
        
        Foam::scalar maxDot = -Foam::GREAT;
        Foam::label maxIdx = -1;
        
        forAll(faces, i)
        {
            // Compute face normal (simplified, assumes planar)
            const Foam::face& f = faces[i];
            if (f.size() < 3) continue;
            
            Foam::vector n = (points[f[1]] - points[f[0]]) ^ (points[f[2]] - points[f[0]]);
            n /= (Foam::mag(n) + Foam::SMALL);
            
            Foam::scalar dot = n & localNormal;
            if (dot > maxDot)
            {
                maxDot = dot;
                maxIdx = i;
                faceNormal = n;
            }
        }
        
        if (maxIdx != -1)
        {
            const Foam::face& f = faces[maxIdx];
            facePoints.setSize(f.size());
            forAll(f, i)
            {
                // Transform points to world frame
                facePoints[i] = p.position() + (p.orientation().R() & points[f[i]]);
            }
            // Transform normal to world frame
            faceNormal = p.orientation().R() & faceNormal;
        }
    }
    
    // Clip polygon against a plane
    Foam::List<Foam::point> clip(
        const Foam::List<Foam::point>& polygon,
        const Foam::point& planePoint,
        const Foam::vector& planeNormal
    )
    {
        Foam::List<Foam::point> clipped;
        if (polygon.size() == 0) return clipped;
        
        Foam::point v1 = polygon.last();
        Foam::scalar d1 = (v1 - planePoint) & planeNormal;
        
        for (const auto& v2 : polygon)
        {
            Foam::scalar d2 = (v2 - planePoint) & planeNormal;
            
            if (d1 >= 0 && d2 >= 0)
            {
                // Both inside/on plane
                clipped.append(v2);
            }
            else if (d1 >= 0 && d2 < 0)
            {
                // Exit plane
                Foam::scalar t = d1 / (d1 - d2);
                clipped.append(v1 + t * (v2 - v1));
            }
            else if (d1 < 0 && d2 >= 0)
            {
                // Enter plane
                Foam::scalar t = d1 / (d1 - d2);
                clipped.append(v1 + t * (v2 - v1));
                clipped.append(v2);
            }
            
            v1 = v2;
            d1 = d2;
        }
        
        return clipped;
    }

    void gjkContact::generateManifold(
        const particle& pA,
        const particle& pB,
        const Foam::vector& normal,
        Foam::List<Foam::point>& contactPoints
    ) const
    {
        contactPoints.clear();
        
        // 1. Identify Reference and Incident faces
        // Normal points from B to A
        // Reference face on A (normal roughly matches contact normal)
        // Incident face on B (normal roughly matches -contact normal)
        
        Foam::List<Foam::point> faceA;
        Foam::vector normalA;
        getBestFace(pA, normal, faceA, normalA);
        
        Foam::List<Foam::point> faceB;
        Foam::vector normalB;
        getBestFace(pB, -normal, faceB, normalB);
        
        if (faceA.size() == 0 || faceB.size() == 0) return;
        
        // Determine which is reference and which is incident
        // Usually the one more parallel to the collision normal is the reference
        // But here we can just pick A as reference for simplicity or check alignment
        
        // Let's assume A is reference
        Foam::List<Foam::point> incident = faceB;
        Foam::List<Foam::point> reference = faceA;
        
        // Clip incident face against side planes of reference face
        // Side planes are formed by edge x normalA
        
        Foam::point v1 = reference.last();
        for (const auto& v2 : reference)
        {
            Foam::vector edge = v2 - v1;
            Foam::vector sideNormal = edge ^ normalA; // Pointing inward?
            // Check winding order. Assuming CCW and normal pointing out.
            // Edge x Normal -> Tangent. Side normal should point IN.
            // If CCW, Edge x Normal points OUT. So use Normal x Edge?
            // Let's try Normal x Edge.
            sideNormal = normalA ^ edge;
            sideNormal /= (Foam::mag(sideNormal) + Foam::SMALL);
            
            incident = clip(incident, v1, sideNormal);
            
            if (incident.size() == 0) break;
            
            v1 = v2;
        }
        
        // Finally, keep points below the reference face
        // Project points onto reference face plane and check distance
        Foam::point refPoint = reference[0];
        
        for (const auto& p : incident)
        {
            Foam::scalar d = (p - refPoint) & normalA;
            if (d <= tolerance_) // Below or on surface
            {
                contactPoints.append(p);
            }
        }
    }
    
    void gjkContact::generateManifoldParticleWall(
        const particle& p,
        const wall& w,
        const Foam::vector& normal,
        Foam::List<Foam::point>& contactPoints
    ) const
    {
        contactPoints.clear();
        
        // Particle is incident, Wall is reference
        Foam::List<Foam::point> faceP;
        Foam::vector normalP;
        getBestFace(p, -normal, faceP, normalP); // Normal points from Wall to Particle
        
        if (faceP.size() == 0) return;
        
        // Clip against wall
        // Wall is infinite (or large), so just check points against wall surface
        // For Plane wall:
        if (w.type() == WallType::PLANE)
        {
            for (const auto& pt : faceP)
            {
                if (w.distanceToSurface(pt) <= tolerance_)
                {
                    contactPoints.append(pt);
                }
            }
        }
        else
        {
            // For other walls, simplified check
             for (const auto& pt : faceP)
            {
                if (w.distanceToSurface(pt) <= tolerance_)
                {
                    contactPoints.append(pt);
                }
            }
        }
        
        // If no points found (e.g. vertex collision), add closest point
        if (contactPoints.size() == 0)
        {
            contactPoints.append(w.closestPointOnSurface(p.position()));
        }
    }

}

