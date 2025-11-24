import numpy as np

def write_cube_stl(filename, size):
    half_size = size / 2.0
    # Define vertices of a cube centered at origin
    vertices = np.array([
        [-half_size, -half_size, -half_size],
        [ half_size, -half_size, -half_size],
        [ half_size,  half_size, -half_size],
        [-half_size,  half_size, -half_size],
        [-half_size, -half_size,  half_size],
        [ half_size, -half_size,  half_size],
        [ half_size,  half_size,  half_size],
        [-half_size,  half_size,  half_size]
    ])

    # Define faces (each face is 2 triangles)
    # Normal pointing outwards
    faces = [
        [0, 3, 1], [1, 3, 2], # Bottom (z=-h) - Normal (0,0,-1)
        [4, 5, 7], [5, 6, 7], # Top (z=h) - Normal (0,0,1)
        [0, 1, 5], [0, 5, 4], # Front (y=-h) - Normal (0,-1,0)
        [1, 2, 6], [1, 6, 5], # Right (x=h) - Normal (1,0,0)
        [2, 3, 7], [2, 7, 6], # Back (y=h) - Normal (0,1,0)
        [3, 0, 4], [3, 4, 7]  # Left (x=-h) - Normal (-1,0,0)
    ]

    with open(filename, 'w') as f:
        f.write(f"solid cube_{size}m\n")
        for face in faces:
            v1, v2, v3 = vertices[face[0]], vertices[face[1]], vertices[face[2]]
            # Calculate normal
            normal = np.cross(v2 - v1, v3 - v1)
            norm_len = np.linalg.norm(normal)
            if norm_len > 0:
                normal /= norm_len
            
            f.write(f"facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n")
            f.write("outer loop\n")
            f.write(f"vertex {v1[0]:.6f} {v1[1]:.6f} {v1[2]:.6f}\n")
            f.write(f"vertex {v2[0]:.6f} {v2[1]:.6f} {v2[2]:.6f}\n")
            f.write(f"vertex {v3[0]:.6f} {v3[1]:.6f} {v3[2]:.6f}\n")
            f.write("endloop\n")
            f.write("endfacet\n")
        f.write(f"endsolid cube_{size}m\n")

if __name__ == "__main__":
    # 0.1m cube
    write_cube_stl("constant/triSurface/cube_0.1m.stl", 0.1)
    print("Generated constant/triSurface/cube_0.1m.stl")
