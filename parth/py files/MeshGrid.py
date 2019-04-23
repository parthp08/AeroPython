import math
import numpy as np
import matplotlib.pyplot as plt

def mesh_grid(x_start, x_end, y_start, y_end, N, grid_properties=False, plotting=False, width=10.0):
    '''
    Generate 2D mesh grid.
    - used by other flow functions.
    
    Parameters:
    -------------------
    x_start, x_end : float, boundaries in the x-direction
    y_start, y_end : float, boundaries in the y-direction
    N : int, number of points in each direction
    grid_properties : bool, if True retruns mesh properties, 
                    x_start,x_end,y_start,y_end,X,Y
    plotting : bool, plot the meshgrid, default is False
    width : float, width of the plot, default is 10.0
    
    Returns:
    --------------
    gridpoints : if True
    mesh_plot : if True 
    None
    '''
    
    x = np.linspace(x_start, x_end, N)    # creates a 1D-array with the x-coordinates
    y = np.linspace(y_start, y_end, N)    # creates a 1D-array with the y-coordinates
    
    X, Y = np.meshgrid(x, y)              # generates a mesh grid
    
    if plotting:
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.scatter(X, Y, s=5, color='#CD2305', marker='o')
    
    if grid_properties:
        return (x_start, x_end, y_start, y_end, X, Y, N)
    return None

if __name__=='__main__':
    mesh_grid(-2.0,2.0,-1.0,1.0,50,plotting=True)
    plt.show()
