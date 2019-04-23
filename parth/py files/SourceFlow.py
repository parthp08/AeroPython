import math
import numpy as np
import matplotlib.pyplot as plt
from MeshGrid import mesh_grid

def sourceFlow(x_source, y_source, source_strength, mesh_properties, velocities_output=True, plotting=True, width=10):
    '''
    Source Flow function.
    - A source is a point from which we imagine that fluid is flowing out, 
      uniformly. Thus, all the streamlines radiate from a single point as 
      straight lines and the radial velocity decreases with the distance 
      from the source point.
      
      - uses meshgrid Function to generate 2D grid
      
    Parameters:
    ----------------
    x_source : float, x-coordinate of source,
                must in inside the boundaries of meshgrid
    y_source : float, y-coordiante of source,
                must in inside the boundaries of meshgrid
    source_strength : positive float, strength of the source
    velocities_output : bool, if True, retruns the x and y velocites 
                        of source flow on meshgrid
    plotting : bool, plot the sourceflow velocity  field, default is True
    width : float, width of the plot, default is 10.0
    
    Returns:
    --------------
    velocities_output : if True
    velocityField_plot : if True 
    None
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # compute the velocity field of source flow on the mesh grid
    u_source = (source_strength / (2 * math.pi) *
                (X - x_source) / ((X - x_source)**2 + (Y - y_source)**2))
    v_source = (source_strength / (2 * math.pi) *
                (Y - y_source) / ((X - x_source)**2 + (Y - y_source)**2))
    
    if plotting:
        # plot the streamlines
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_source, v_source,
                          density=2, linewidth=1, arrowsize=2, arrowstyle='->')
        plt.scatter(x_source, y_source,
                       color='#CD2305', s=80, marker='o')
    
    if velocities_output:
        return u_source, v_source
    return None

if __name__=='__main__':
    # Define Mesh Grid
    # change mesh_properties as per requirement
    mesh_properties = mesh_grid(-4.0,4.0,-2.0,2.0,200, grid_properties=True)

    sourceFlow(-1.0, 0.0, 5.0, mesh_properties, velocities_output=False)
    # change source_strength and width to observe diffrent results
    plt.show()

