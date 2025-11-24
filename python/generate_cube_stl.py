import numpy as np

def write_cube_stl(filename, size):
    s = size / 2.0
    # Vertices
    vertices = np.array([
        [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
        [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
    ])
    
    # Faces (normal pointing out)
    faces = [
        [0, 3, 2], [0, 2, 1], # Bottom
        [4, 5, 6], [4, 6, 7], # Top
        [0, 1, 5], [0, 5, 4], # Front
        [1, 2, 6], [1, 6, 5], # Right
        [2, 3, 7], [2, 7, 6], # Back
        [3, 0, 4], [3, 4, 7]  # Left
    ]
    
    with open(filename, 'w') as f:
        f.write("solid cube\n")
        for face in faces:
            # Calculate normal
            v1 = vertices[face[1]] - vertices[face[0]]
            v2 = vertices[face[2]] - vertices[face[0]]
            normal = np.cross(v1, v2)
            normal /= np.linalg.norm(normal)
            
            f.write(f"facet normal {normal[0]} {normal[1]} {normal[2]}\n")
            f.write("outer loop\n")
            for vid in face:
                v = vertices[vid]
                f.write(f"vertex {v[0]} {v[1]} {v[2]}\n")
            f.write("endloop\n")
            f.write("endfacet\n")
        f.write("endsolid cube\n")

if __name__ == "__main__":
    write_cube_stl("/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/cfd_dem/freeSedimentation/constant/triSurface/cube.stl", 0.002)
