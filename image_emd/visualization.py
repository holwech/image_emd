import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


def single_scatter_plot(x, y, z, figsize=(15,10), color='gray'):
    fig = plt.figure(figsize=figsize)
    ax = axes3d.Axes3D(fig)
    ax.scatter3D(x, y, z, c=color ,marker=',')

def double_scatter_plot(x1, y1, z1, x2, y2, z2, figsize=(15,10), color1='blue', color2='red'):
    fig = plt.figure(figsize=figsize)
    ax = axes3d.Axes3D(fig)
    ax.scatter3D(x1, y1, z1, c=color1, hmarker=',')
    ax.scatter3D(x2, y2, z2, c=color2, marker=',')

# X, Y, Z is the grid data attained from the spline.
# x, y, z is the data points used to model the spline.
def single_surface_plot(X, Y, Z, x, y, z, figsize=(15,10), color='gray', show_points=True, show_wireframe=False):
    fig = plt.figure(figsize=figsize)
    ax = axes3d.Axes3D(fig)
    ax.plot_surface(X, Y, Z,alpha=0.2, color=color)

    if show_points:
        ax.scatter3D(x, y, z, c=color)
    if show_wireframe:
        ax.plot_wireframe(X, Y, Z, color=color)


def double_surface_plot(X1, Y1, Z1, X2, Y2, Z2, x1, y1, z1, x2, y2, z2,
                        figsize=(15,10), color1='blue', color2='red',
                        show_points=True, show_wireframe=False):
    fig = plt.figure(figsize=figsize)
    ax = axes3d.Axes3D(fig)

    ax.plot_surface(X1, Y1, Z1,alpha=0.2, color=color1)
    ax.plot_surface(X2, Y2, Z2,alpha=0.2, color=color2)

    if show_points:
        ax.scatter3D(x1, y1, z1, c=color1)
        ax.scatter3D(x2, y2, z2, c=color2)
    if show_wireframe:
        ax.plot_wireframe(X1, Y1, Z1, color=color1)
        ax.plot_wireframe(X2, Y2, Z2, color=color2)


