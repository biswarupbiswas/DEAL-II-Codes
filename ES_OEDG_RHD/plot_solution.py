import matplotlib.pyplot as plt
import numpy as np
import glob
import os

def parse_vtk_simple(filename):
    """
    Robust simple VTK parser for deal.II ASCII output.
    """
    with open(filename, 'r') as f:
        content = f.read()

    # Find Points
    pts_start = content.find("POINTS")
    if pts_start == -1: return None
    pts_header = content[pts_start:content.find("\n", pts_start)]
    num_pts = int(pts_header.split()[1])
    
    pts_data_start = content.find("\n", pts_start) + 1
    # Points are (x,y,z), so num_pts * 3 values
    pts = np.fromstring(content[pts_data_start:], sep=' ', count=num_pts*3).reshape(-1, 3)

    # Find Scalar D
    d_start = content.find("SCALARS D")
    if d_start == -1: return None
    # Skip header and lookup table lines
    data_start = content.find("LOOKUP_TABLE default", d_start)
    data_start = content.find("\n", data_start) + 1
    
    d = np.fromstring(content[data_start:], sep=' ', count=num_pts)
    
    return pts, d

def plot_solution(filename):
    print(f"Parsing {filename}...")
    data = parse_vtk_simple(filename)
    if data is None:
        print("Parsing failed.")
        return
    
    points, d = data
    x, y = points[:, 0], points[:, 1]
    
    plt.figure(figsize=(12, 10))
    # Use tricontourf for high-quality unstructured density map
    cnt = plt.tricontourf(x, y, d, levels=200, cmap='magma')
    plt.colorbar(cnt, label='Density D')
    plt.title('RHD 2D Riemann Problem (RP2D1)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    
    out_file = 'density_plot.png'
    plt.savefig(out_file, dpi=300)
    print(f"Plot saved to {out_file}")

if __name__ == "__main__":
    files = sorted(glob.glob("solution-*.vtk"))
    if files:
        plot_solution(files[-1])
    else:
        print("No solution files found.")
