import numpy as np
import pyvista as pv
import glob
import os
import sys

def con2prim(D, Mx, My, E, gas_gamma=5.0/3.0):
    gm1 = gas_gamma - 1.0
    M2 = Mx**2 + My**2
    M = np.sqrt(M2)
    den = np.maximum(gm1**2 * (M2 + D**2), 1e-20)
    
    a3 = -2.0 * gas_gamma * gm1 * M * E / den
    a2 = (gas_gamma**2 * E**2 + 2.0 * gm1 * M2 - gm1**2 * D**2) / den
    a1 = -2.0 * gas_gamma * M * E / den
    a0 = M2 / den
    
    eps = 1e-10
    v = np.zeros_like(D)
    mask = M > eps
    if np.any(mask):
        Dm, Em, Mm = D[mask], E[mask], M[mask]
        a3m, a2m, a1m, a0m = a3[mask], a2[mask], a1[mask], a0[mask]
        vl = (1.0 / (2.0 * Mm * gm1)) * \
             (gas_gamma * Em - np.sqrt(np.maximum(0, (gas_gamma * Em)**2 - 4.0 * gm1 * Mm**2)))
        vu = Mm / Em
        z = np.where(vl < eps, 0.5 * (1.0 - Dm / Em) * (vl - vu), 0.0)
        vm = 0.5 * (vl + vu) + z
        for _ in range(10):
            v2 = vm * vm
            v3 = v2 * vm
            v4 = v3 * vm
            vm = vm - (v4 + a3m * v3 + a2m * v2 + a1m * vm + a0m) / \
                      (4.0 * v3 + 3.0 * a3m * v2 + 2.0 * a2m * vm + a1m)
        v[mask] = vm
    
    v = np.clip(v, 0, 1.0 - 1e-15)
    L = 1.0 / np.sqrt(1.0 - v**2)
    return D / L

def process_file(filename):
    mesh = pv.read(filename)
    rho = con2prim(mesh.point_data['D'], mesh.point_data['Mx'], 
                   mesh.point_data['My'], mesh.point_data['E'])
    mesh.point_data['log10_rho'] = np.log10(np.maximum(rho, 1e-10))
    return mesh

def setup_plotter(plotter, mesh, clim=None):
    # The "Easy way": Use built-in scalar_bar_args
    sbar_kwargs = {
        "title": "log10(rho)",
        "vertical": False,
        "position_x": 0.1,
        "position_y": 0.9,
        "width": 0.8,
        "height": 0.05,
        "color": "black",
        "title_font_size": 18,
        "label_font_size": 14,
        "fmt": "%.2f"
    }
    
    plotter.set_background("white")
    plotter.add_mesh(mesh, scalars="log10_rho", cmap="jet", clim=clim, 
                     scalar_bar_args=sbar_kwargs, show_edges=False)
    
    # Simple 2D view with axes
    plotter.view_xy()
    # Let PyVista handle the camera automatically for a tight fit
    plotter.reset_camera()

def plot_single(filename):
    print(f"Plotting {filename}...")
    mesh = process_file(filename)
    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=[1000, 1000])
    setup_plotter(plotter, mesh)
    output_name = filename.replace('.vtk', '_log_rho.png')
    plotter.screenshot(output_name)
    plotter.close()
    print(f"Saved {output_name}")

def animate(filenames, output_gif='solution_animation.gif'):
    print(f"Creating animation: {output_gif}")
    pv.OFF_SCREEN = True
    
    print("Finding global limits...")
    all_min, all_max = [], []
    for f in filenames[::max(1, len(filenames)//10)]:
        m = process_file(f)
        all_min.append(m.point_data['log10_rho'].min())
        all_max.append(m.point_data['log10_rho'].max())
    clim = [min(all_min), max(all_max)]

    plotter = pv.Plotter(off_screen=True, window_size=[1000, 1000])
    plotter.open_gif(output_gif)
    for i, filename in enumerate(filenames):
        print(f"Frame {i+1}/{len(filenames)}")
        mesh = process_file(filename)
        plotter.clear()
        setup_plotter(plotter, mesh, clim=clim)
        plotter.write_frame()
    plotter.close()
    print(f"Animation saved to {output_gif}")

if __name__ == "__main__":
    files = sorted(glob.glob("solution-*.vtk"))
    if not files:
        sys.exit("No .vtk files found.")
    if "--animate" in sys.argv:
        animate(files)
    else:
        plot_single(sys.argv[1] if len(sys.argv) > 1 else files[-1])
