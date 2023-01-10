import matplotlib.pyplot as plt
import matplotlib
import numpy as np

cm = matplotlib.cm.magma
norm = matplotlib.colors.Normalize(vmin=-0.5, vmax=0.5)

def draw_tri_lattice(Lx, Ly, r0=(0,0), color='lightgray'):
    az = .5*np.sqrt(3) # vertical displacement for equilateral triangles
    x0, y0 = r0
    x, y = x0, y0
    for i in range(Lx):
        for j in range(Ly):
            xa = x + 1
            ya = y
            xb = x + .5
            yb = y + az
            plt.plot((xb, x, xa, xb), (yb, y, ya, yb), color=color, zorder=0, linewidth=1)
            x = xa
        if i%2 == 1:
            x = x0
            y = y0 + (i+1)*az
        else:
            x = x0 - .5
            y = y0 + (i+1)*az

            

    
def draw_cluster(plaq, r0=(0,0), color='purple'):
    x0, y0 = r0
    plt.plot(plaq['outline'][0]+x0, plaq['outline'][1]+y0,
             color=color, zorder=1, linewidth=1)

def plot_spins(plaq, mf0, r0=(0,0)):
    x0, y0 = r0
    zs = np.real(mf0['z'])
    plt.quiver(plaq['rs'][0]+x0, plaq['rs'][1]+y0, mf0['x'], mf0['y'], zorder=10, color=cm(norm(zs)))