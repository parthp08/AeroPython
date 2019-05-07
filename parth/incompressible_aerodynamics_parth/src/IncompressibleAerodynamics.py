## Incompressible Aerodynamics Library

import numpy as np
import matplotlib.pyplot as plt

## TODO:=>  Doc-strings for all functions

#--------------------------------------------------------------------------------------------------------------------------------------------
def Avilable_Functions():
	"""
	Prints all functions avilable in Aerodynamics library
	"""
	print(
		"""
- NACA4digit
- mesh grid
- Source
- Sink
- Source_Sink_Pair
- FreestreamFlow
- Source_plus_Freestream
- SourceSink_plus_Freestream
- Doublet
- Doublet_plus_Freestream
- NonLifting_flow_over_Cylinder
- NonLifting_flow_over_NACA0012
- VortexFlow
- Vortex_Sink
- Vortex_Source
- VortexSheet
- Horizontal_Infinite_VortexRow
- Lifting_flow_over_Cylinder
		""")

#--------------------------------------------------------------------------------------------------------------------------------------------
def NACA4digit(airfoil_number, points=1000, get_data=False, plotting=True, closed_trailing_edge=True, save_to_file=False):
    """Compute and plot coordinates of a NACA 4-digit airfoil.
    
    Arguments:
        number: Name of the requested 4-digt NACA airfoil entered as a string,
                i.e. '2412'
    
    Optional arguments:
        points: Number of desired airfoil coordinates
                (half(N) for upper and half(N) for lower surface)
        closed_trailing_edge: The trailing edge has a small but finite
                            thickness when using equations in [1] unaltered.
        save_to_file: Saves the coordinates in a CSV file.
        get_data: returns x and y coordinates
        plotting: plot the airfoil geometry
    
    Returns:
        A matrix with two columns pertaining to the x and y-coordinates,
        respectively. The sequence of coordinates is clockwise, starting at the
        trailing edge.
    
    Raises:
        ValueError: Airfoil number does not have four digits or N is negative
        
    References:
        https://en.wikipedia.org/wiki/NACA_airfoil
    
    """

    ## Check whether input parameters are valid or not
    if not (0 < int(airfoil_number) < 10000 and points > 0):
        raise ValueError("Invalid input.")

    c = 1       # by default chord(c) is set to 1; so x = c, x/c = 1

    ## 4-digit airfoil properties
    M = int(airfoil_number[0]) * (c/100)    # max camber
    P = int(airfoil_number[1]) * (c/10)     # max camber location from leading edge
    T = int(airfoil_number[2:4]) * (c/100)  # max thickness

    ## coordinates points
    N = points//2

    ## spacing of coordinates along the x-axis
    x = np.linspace(0, c, N)

    ## The trailing edge has a small but finite thickness by default. The gap
    # can be closed by utilizing a slightly different equation.
    if closed_trailing_edge:
        a4 = -0.1036   # Closed trailing edge
    else:
        a4 = -0.1015   # Open trailing edge

    ## constants
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843

    ## Computing the y-coordinates of the camber line.
    # Camber & Gradient
    fwd_x = x[x < P]
    aft_x = x[x >= P]
    # camber calcultation
    if 0 < P < c:
        fwd_camber_yc = M / P**2 * (2 * P * fwd_x - np.power(fwd_x, 2))
        fwd_camber_dyc_dx = 2*M / P**2 * (P - fwd_x)
        aft_camber_yc = M / (1 - P)**2 * ((1 - 2 * P) + 2 * P * aft_x - np.power(aft_x, 2))
        aft_camber_dyc_dx = 2*M / (1-P)**2 * (P - aft_x)
        camber_yc = np.append(fwd_camber_yc, aft_camber_yc)   # complete mean camber line
        dyc_dx_camber = np.append(fwd_camber_dyc_dx, aft_camber_dyc_dx)
    else:
        camber_yc = np.zeros(np.size(x))   # cases where max_camber is located on leading or triling edge
        dyc_dx_camber = np.zeros(np.size(x))
    #gradient calculation
    theta = np.arctan(dyc_dx_camber)
    
    ## Thickness Distribution
    yt = 5*T*((a0*np.sqrt(x)) + (a1*x) + (a2*(x**2)) + (a3*(x**3)) + (a4*(x**4)))

    ## Upper surface points
    xu = x - (yt * np.sin(theta))
    yu = camber_yc + (yt* np.cos(theta))

    ## Lower surface points
    xl = x + (yt * np.sin(theta))
    yl = camber_yc - (yt * np.cos(theta))
    
    ## coordinates
    x_cord = np.append(xu, xl)
    y_cord = np.append(yu, yl)
    coordinates = np.column_stack((x_cord, y_cord))
    
    ## Saving data to CSV file
    if save_to_file:
        np.savetxt(f"NACA{airfoil_number}_plotting_data.csv", coordinates, delimiter=',')

    if plotting:
        ## Plot the airfoil
        plt.plot(x_cord, y_cord)
        plt.grid()
        plt.title(f"NACA{airfoil_number} airfoil")
        plt.axis('equal')
        plt.show()
    
    if get_data:
        return x_cord, y_cord
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
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
        plt.title("Mesh Grid")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.scatter(X, Y, s=5, color='#CD2305', marker='o')
        plt.show()
    
    if grid_properties:
        return (x_start, x_end, y_start, y_end, X, Y, N)
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Source(x_source, y_source, source_strength, mesh_properties, get_data=True, plotting=False, width=10):
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
    mesh_properties : Tuple, mesh grid properties
    get_data : bool, if True, retruns the values of u, v and psi
    plotting : bool, plot the sourceflow velocity  field, default is False
    width : float, width of the plot, default is 10.0
    
    Returns:
    --------------
    u_source, v_source : velocity field
    psi_source: stream-function
    velocityField_plot : if True 
    None
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # compute the velocity field of source flow on the mesh grid
    u_source = (source_strength / (2 * np.pi) *
                (X - x_source) / ((X - x_source)**2 + (Y - y_source)**2))
    v_source = (source_strength / (2 * np.pi) *
                (Y - y_source) / ((X - x_source)**2 + (Y - y_source)**2))

    # compute the stream function
    psi_source = source_strength / (2 * np.pi) * np.arctan2((Y - y_source), (X - x_source))

    
    if plotting:
        # plot the streamlines
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Source Flow")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_source, v_source,
                          density=2, linewidth=1, arrowsize=2, arrowstyle='->')
        plt.scatter(x_source, y_source,
                       color='#CD2305', s=80, marker='o')
        plt.show()

    if get_data:
        return u_source, v_source, psi_source
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Sink(x_sink, y_sink, sink_strength, mesh_properties, get_data=True, plotting=False, width=10):
    '''
    Sink Flow function.
    - In the source flow, the strength was chosen to be positive. 
      a source with a negative strength is called a sink. Instead of 
      radiating from a single point, the straight streamlines are now 
      converging to a single point.
    - The velocity field corresponding to a sink looks similar to that 
      of a source, except for the direction of the flow.
      
    - uses meshgrid Function to generate 2D grid
      
    Parameters:
    ----------------
    x_sink : float, x-coordinate of sink,
                must in inside the boundaries of meshgrid
    y_sink : float, y-coordiante of sink,
                must in inside the boundaries of meshgrid
    sink_strength : negative float, strength of the sink
    mesh_properties : Tuple, mesh grid properties
    get_data : bool, if True, retruns the values of u, v and psi
    plotting : bool, plot the sinkflow velocity field, default is False
    width : float, width of the plot, default is 10.0
    
    Returns:
    --------------
    velocities_output : if True
    velocityField_plot : if True 
    None
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # compute the velocity field of sink flow on the mesh grid
    u_sink = (sink_strength / (2 * np.pi) *
                (X - x_sink) / ((X - x_sink)**2 + (Y - y_sink)**2))
    v_sink = (sink_strength / (2 * np.pi) *
                (Y - y_sink) / ((X - x_sink)**2 + (Y - y_sink)**2))
    
    # compute the stream function
    psi_sink = sink_strength / (2 * np.pi) * np.arctan2((Y - y_sink), (X - x_sink))

    if plotting:
        # plot the streamlines
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Sink Flow")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_sink, v_sink,
                          density=2, linewidth=1, arrowsize=2, arrowstyle='->')
        plt.scatter(x_sink, y_sink,
                       color='#CD2305', s=80, marker='o')
        plt.show()

    if get_data:
        return u_sink, v_sink, psi_sink
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Source_Sink_Pair(x_source, y_source, x_sink, y_sink, sourceSinkStrength, mesh_properties, get_data=True, plotting=False, width=10.0):
    '''
    Source-Sink Pair Function
    - adding source and sink by superpositon.
    
    - uses sourceFlow, sinkFlow and meshgrid Functions
    
    Parameters:
    ----------------
    x_source : float, x-coordinate of source,
                must in inside the boundaries of meshgrid
    y_source : float, y-coordiante of source,
                must in inside the boundaries of meshgrid
    x_sink : float, x-coordinate of sink,
                must in inside the boundaries of meshgrid
    y_sink : float, y-coordiante of sink,
                must in inside the boundaries of meshgrid
    sourceSinkStrength : float, strength of the source and sink
                        +ve for source, -ve for sink.
    mesh_properties : Tuple, mesh grid properties
    get_data : bool, if True, retruns the values of u, v and psi
    plotting : bool, plot the source-sink flow velocity field, default is False
    width : float, width of the plot, default is 10.0
    
    Returns:
    --------------
    u_sink, v_sink : velocity field
    psi_sink: stream-function
    velocityField_plot : if True 
    None
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # getting source and sink velocities
    u_source,v_source,_ = Source(x_source,y_source,sourceSinkStrength,mesh_properties,plotting=False)
    u_sink,v_sink,_ = Sink(x_sink,y_sink,-sourceSinkStrength,mesh_properties,plotting=False)
    
    # compute the velocity of the pair source/sink by superposition
    u_pair = u_source + u_sink
    v_pair = v_source + v_sink
    
    if plotting:
        # plot the streamlines of the pair source/sink
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Source-Sink Pair")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_pair, v_pair,
                          density=2.0, linewidth=1, arrowsize=2, arrowstyle='->')
        plt.scatter([x_source, x_sink], [y_source, y_sink], 
                       color='#CD2305', s=80, marker='o')
        plt.show()
    
    if get_data:
        return u_pair, v_pair
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def FreestreamFlow(u_inf, mesh_properties):
    '''
    Calculate freestream(Uniform Flow) velocity and streamline field
    
    - uses mesh_grid Function
    
    Parameters:
    ------------
    u_inf : float, free-stream speeed
    mesh_properties: Tuple, grid properties
    
    Returns:
    ---------
    u_freestream ,v_freestream : free-stream velocity field
    psi_freestream : stream-function 
    '''
    
    N = mesh_properties[-1]  # number of grid points
    Y = mesh_properties[-2]  # messgrid points array in y-direction
    
    # compute the freestream velocity field
    u_freestream = u_inf * np.ones((N,N), dtype=float)
    v_freestream = np.zeros((N,N), dtype=float)
    
    # compute the stream-function
    psi_freestream = u_inf * Y
    
    return u_freestream, v_freestream, psi_freestream

#--------------------------------------------------------------------------------------------------------------------------------------------
def Source_plus_Freestream(x_source, y_source, strength_source, u_inf, mesh_properties, get_data=False, plotting=False, width=10):
    '''
    Source Flow + Freesream(Uniform Flow) Function
    - usign superposition
    
    Parameters:
    ------------
    source_strength : positive float, strength of the source
    x_source : float, x-coordinate of source,
                must in inside the boundaries of meshgrid
    y_source : float, y-coordiante of source,
                must in inside the boundaries of meshgrid
    u_inf : float, free-stream speeed
    mesh_properties: Tuple, grid properties
    get_data : bool, if True retruns u,v,x_source,y_source,psi,u_inf
    plotting : bool, plot the meshgrid, default is True
    width: positive float, width of the plot, default=10.0
    
    Returns:
    ------------
    plotting : if True returns Velocity Flow field plot for source + freestream flow
    get_velocity : u, v if True
    None
    '''
    
    # getting values from source and  freestream function
    u_source, v_source, psi_source = Source(x_source, y_source, strength_source, mesh_properties)
    u_freestream, v_freestream, psi_freestream = FreestreamFlow(u_inf, mesh_properties)
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # superposition of the source on the freestream
    u = u_freestream + u_source
    v = v_freestream + v_source
    psi = psi_freestream + psi_source

    # calculate the stagnation point(point where velocity is zero)
    x_stagnation = x_source - strength_source/(2 * np.pi * u_inf)
    y_stagnation = y_source
    
    if plotting:
        # plot the streamline
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.grid(True)
        plt.title("Source in Freestream")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(x_source, y_source, color='#CD2305', s=80, marker='o')
        # display stagnation point
        plt.scatter(x_stagnation, y_stagnation, color='g', s=80, marker='o')
        # display dividing streamline
        plt.contour(X, Y, psi, levels=[-strength_source/2, strength_source/2], colors='#CD2305')
        plt.show()

    if get_data:
        return u, v, psi, x_source, y_source, u_inf
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def SourceSink_plus_Freestream(x_source, y_source, x_sink, y_sink, strength_source, strength_sink, u_inf, mesh_properties, get_data=True,  plotting=False, pressure_field_plot=False, width=10):
    '''
    Source Flow + Sink Flow + Freesream(Uniform Flow) Function
    - usign superposition
    
    Parameters:
    ------------
    strength_source : positive float, strength of the source
    strength_sink : positive float, strength of the sink
    x_source : float, x-coordinate of source,
                must in inside the boundaries of meshgrid
    y_source : float, y-coordiante of source,
                must in inside the boundaries of meshgrid
    x_sink : float, x-coordinate of sink,
                must in inside the boundaries of meshgrid
    y_sink : float, y-coordiante of sink,
                must in inside the boundaries of meshgrid
    u_inf : float, free-stream speeed
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u,v,x_source,y_source,psi,u_inf
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    ------------
    plotting : if True returns Velocity Flow field plot for source + sink + freestream flow
    get_velocity : u, v if True
    None
    '''
    
    # getting values from source, sink and  freestream function
    u_source, v_source, psi_source = Source(x_source, y_source, strength_source, mesh_properties)
    u_sink, v_sink, psi_sink = Sink(x_sink, y_sink, strength_sink, mesh_properties)
    u_freestream, v_freestream, psi_freestream = FreestreamFlow(u_inf, mesh_properties)
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # superposition of the source and sink on the freestream
    u = u_freestream + u_source + u_sink
    v = v_freestream + v_source + v_sink
    psi = psi_freestream + psi_source + psi_sink

    # calculate the stagnation point(point where velocity is zero)
    # two stagnation for this case
    x_stagnation_1 = x_source - strength_source/(2 * np.pi * u_inf)
    y_stagnation_1 = y_source
    x_stagnation_2 = x_sink - strength_sink/(2 * np.pi * u_inf)
    y_stagnation_2 = y_sink
    
    if plotting:
        # plot the streamline
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.grid(True)
        plt.title("Source-Sink pair in Freestream")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=1.5, arrowstyle='->')
        plt.scatter([x_source,x_sink], [y_source,y_sink], color='#CD2305', s=80, marker='o')
        # display stagnation points
        plt.scatter([x_stagnation_1,x_stagnation_2], [y_stagnation_1, y_stagnation_2], color='g', s=80, marker='o')
        # display dividing streamline
        plt.contour(X, Y, psi,
                   levels=[0,],
                   colors='#CD2305')
        plt.show()

    if pressure_field_plot:
        # Pressure Coefficient Field
        #Pressure Field Visulization for Source + Sink + Freestream Flow
            
        # compute the pressure coefficient field
        cp = 1.0 - (u**2 + v**2)/(u_inf**2)
        
        # plot the pressure coefficient field
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(1.1 * width, height))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        contf = plt.contourf(X, Y, cp,
                                levels=np.linspace(-2.0, 1.0, 100), extend='both')
        cbar = plt.colorbar(contf)
        cbar.set_label('$C_p$', fontsize=16)
        cbar.set_ticks([-2.0, -1.0, 0.0, 1.0])
        plt.scatter([x_source, x_sink], [y_source, y_sink],
                       color='#CD2305', s=80, marker='o')
        plt.contour(X, Y, psi,
                       levels=[0.], colors='#CD2305', linewidths=2, linestyles='solid')
        plt.show()
        
    if get_data:
        return u, v, psi, x_source, y_source, x_sink, y_sink, u_inf
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Doublet(strength_doublet, x_doublet, y_doublet, mesh_properties, get_data=True, plotting=False, width=10):
    '''
    Doublet Flow function.
    
    If you look from far enough away, the distance between source 
    and sink approaches zero, and the pattern you see is called a doublet.
    
    - uses get_velocity, get_stream_function, mesh_grid Functions
    
    Parameters:
    ----------------
    doublet_strength : positive float, strength of the doublet
    x_doublet : float, x-coordinate of doublet,
                must in inside the boundaries of meshgrid
    y_doublet : float, y-coordiante of doublet,
                must in inside the boundaries of meshgrid
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u,v,x_doublet,y_doublet,psi,u_inf
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    --------------
    u_doublet, v_doublet : velocity field
    psi_doublet: stream-function
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # compute the velocity field
    u_doublet = (- strength_doublet / (2 * np.pi) * ((X - x_doublet)**2 - (Y - y_doublet)**2) / ((X - x_doublet)**2 + (Y - y_doublet)**2)**2)
    v_doublet = (- strength_doublet / (2 * np.pi) * 2 * (X - x_doublet) * (Y - y_doublet) / ((X - x_doublet)**2 + (Y - y_doublet)**2)**2)

    # compute the stream-function
    psi_doublet = - strength_doublet / (2 * np.pi) * (Y - y_doublet) / ((X - x_doublet)**2 + (Y - y_doublet)**2)
    
    if plotting:
        # plot the streamline
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.grid(True)
        plt.title("Doublet Flow")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_doublet, v_doublet, density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(x_doublet, y_doublet, color='#CD2305', s=80, marker='o')
        plt.show()
       
    if get_data:
        return u_doublet, v_doublet, psi_doublet
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Doublet_plus_Freestream(strength_doublet, x_doublet, y_doublet, u_inf, mesh_properties, get_data=True,  plotting=False, pressure_field_plot=False, width=10):
    '''
    Freesream(Uniform Flow) + Doublet Function
    - usign superposition
    # uniform(Freestream) Flow + Doublet Flow    == Non-Lifting, No-Drag Flow Over Cylinder
    
    Parameters:
    ------------
    strength_doublet : positive float, strength of the doublet
    x_doublet : float, x-coordinate of doublet,
                must in inside the boundaries of meshgrid
    y_doublet : float, y-coordiante of doublet,
                must in inside the boundaries of meshgrid
    u_inf : float, free-stream speeed
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u,v,x_doublet,y_doublet,psi,u_inf,stagnation points
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    ------------
    plotting : if True returns Velocity Flow field plot for Doublet + Freestream flow
    get_velocity : u, v if True
    None
    '''
    
    # getting values from Doublet and  freestream function
    u_doublet, v_doublet, psi_doublet = Doublet(strength_doublet, x_doublet, y_doublet, mesh_properties, get_data=True, plotting=False)
    u_freestream, v_freestream, psi_freestream = FreestreamFlow(u_inf, mesh_properties)
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # superposition of the source and sink on the freestream
    u = u_freestream + u_doublet
    v = v_freestream + v_doublet
    psi = psi_freestream + psi_doublet

    # calculate the stagnation point(point where velocity is zero)
    # two stagnation for this case
    kappa = strength_doublet   # k = strength of the doublet
    x_stagn1, y_stagn1 = +np.sqrt(kappa / (2 * np.pi * u_inf)), 0.0
    x_stagn2, y_stagn2 = -np.sqrt(kappa / (2 * np.pi * u_inf)), 0.0

    # calculate the cylinder radius to add cyliner in figure
    R =  np.sqrt(kappa / (2 *  np.pi * u_inf))
    
    if plotting:
        # plot the streamline
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.grid(True)
        plt.title("Doublet in Freestream / NonLifting_flow_over_Cylinder")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=1.5, arrowstyle='->')
        plt.scatter(x_doublet, y_doublet, color='#CD2305', s=80, marker='o')
        # plot diving line
        # plt.contour(X, Y, psi,
        #        levels=[0.], colors='#CD2305', linewidths=2, linestyles='solid')
        # display stagnation points
        plt.scatter([x_stagn1,x_stagn2], [y_stagn1, y_stagn2], color='g', s=80, marker='o')
        # display circle
        circle = plt.Circle((0, 0), radius=R, color='#CD2305', alpha=0.5)
        plt.gca().add_patch(circle)
        plt.show()
    
    if pressure_field_plot:
        # Pressure Coefficient Field
        # compute the pressure coefficient field
        cp = 1.0 - (u**2 + v**2)/(u_inf**2)
        
        # get mesh_properties
        x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
        
        # plot the pressure coefficient field
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(1.1 * width, height))
        plt.title("Pressure Coefficient Field of NonLifting_flow_over_Cylinder")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        contf = plt.contourf(X, Y, cp,
                                levels=np.linspace(-2.0, 1.0, 100), extend='both')
        cbar = plt.colorbar(contf)
        cbar.set_label('$C_p$', fontsize=16)
        cbar.set_ticks([-2.0, -1.0, 0.0, 1.0])
        plt.scatter(x_doublet, y_doublet,
                       color='#CD2305', s=80, marker='o')
        plt.contour(X,Y,psi,
                       levels=[0.], colors='#CD2305', linewidths=2, linestyles='solid')
        plt.scatter([x_stagn1, x_stagn2], [y_stagn1, y_stagn2],
                       color='g', s=80, marker='o')
                # display circle
        circle = plt.Circle((0, 0), radius=R, color='#CD2305', alpha=0.5)
        plt.gca().add_patch(circle)
        plt.show()

    if get_data:
        return u, v, psi, x_doublet, y_doublet, x_stagn1, x_stagn2, y_stagn1, y_stagn2, u_inf
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def NonLifting_flow_over_Cylinder(strength, x_doublet, y, u_inf, mesh_properties, get_data=True,  plotting=False, pressure_field_plot=False, width=10):
    """
    Non-Lifting flow over a Cylinder

    -- uses function Doublet_plus_Freestream
    """
    return Doublet_plus_Freestream(strength, x_doublet, y, u_inf, mesh_properties, get_data=get_data,  plotting=plotting, pressure_field_plot=pressure_field_plot, width=width)

#--------------------------------------------------------------------------------------------------------------------------------------------
def Sum_Source_Velocity(x_airfoil, y_airfoil, sigma_airfoil, mesh_properties):
    '''
    Calculate the sum of all Source sheets on the Airfoil surface

    -- used by NonLifting_flow_over_NACA0012 function
    
    Paramters:
    ------------
    x_airfoil : np.array, x-coordinates of the airfoil
    y_airfoil : np.array, y-coordinates of the airfoil
    sigma_airfoil : np.array, SourceStrength at each coordiantes of airfoil
    mesh_properties : tuple, grid properties
    
    Returns:
    ---------
    u_source : np.array, sum of all sources velocity at every location on grid
                in x-direction
    y_source : np.array, sum of all sources velocity at every location on grid
                in y-direction
    '''
    
    # gettng values from mesh_properties
    X, Y = mesh_properties[-2], mesh_properties[-3]
    
    u_source = np.zeros_like(X)
    v_source = np.zeros_like(Y)
    
    for (x, y, strength) in zip(x_airfoil, y_airfoil, sigma_airfoil):
        u_add, v_add, _ = Source(x, y, strength, mesh_properties, get_data=True, plotting=False)
        u_source += u_add
        v_source += v_add
    
    return u_source, v_source

#--------------------------------------------------------------------------------------------------------------------------------------------
def NonLifting_flow_over_NACA0012(x_airfoil, y_airfoil, sigma_airfoil, u_inf, mesh_properties, get_data=True, plotting=False, pressure_field_plot=False, width=10):
    '''
    Non-lifting, No-Drag Incomressible Flow over Airfoil using
    source distribution method
    
    Paramters:
    ------------
    x_airfoil : np.array, x-coordinates of the airfoil
    y_airfoil : np.array, y-coordinates of the airfoil
    sigma_airfoil : np.array, SourceStrength at each coordiantes of airfoil
    u_inf : positive float, freestream speed
    mesh_properties : tuple, grid properties
    width : float, width of the plot, default is 10.0
    get_data : bool, if True retruns u, v
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    ---------
    Velocity field over NACA0012 Airfoil 
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # getting the velocity values of SourceSet and freestream flow
    u_sourceSet, v_sourceSet = Sum_Source_Velocity(x_airfoil, y_airfoil, sigma_airfoil, mesh_properties)
    u_freestream, v_freestream, _ = FreestreamFlow(u_inf, mesh_properties)
    
    u = u_sourceSet + u_freestream
    v = v_sourceSet + v_freestream
    
    if plotting:
        # plot the streamlines
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, 
                       density=2, linewidth=1, arrowsize=2, arrowstyle='->')
        plt.scatter(x_airfoil, y_airfoil, 
                color='#CD2305', s=20, marker='o')
        plt.show()

    if pressure_field_plot:
        # compute the pressure coefficient field
        cp = 1.0 - (u**2 + v**2)/(u_inf**2)
        
        # Max pressure
        max_pressure = np.unravel_index(np.argmax(cp), cp.shape)
        
        # location of the max pressure in (x,y) coordinates
        def double_index_pair(X, Y, xy):
            '''
            Returns A Tuple of coordiantes of the point on MeshGrid
            - used for plotting
            '''
            x,y = xy
            return (X[x][y], Y[x][y])
        
        # get mesh_properties
        x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
        
        # getting the airfoil coordinates from Airfoil Plotter
        x_cord_airfoil, y_cord_airfoil = NACA4digit('0012',points=1000, get_data=True, plotting=False)
        
        # plot the pressure coefficient field
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(1.1 * width, height))
        plt.title("NonLifting_flow_over_NACA0012")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        contf = plt.contourf(X, Y, cp,
                                levels=np.linspace(-0.6, 0.75, 140), extend='both')
        cbar = plt.colorbar(contf)
        cbar.set_label('$C_p$', fontsize=16)
        cbar.set_ticks([-0.5, -0.25, 0.0, 0.25, 0.5, 0.75])
        # plots airfoil from airfoil plotter
        plt.plot(x_cord_airfoil, y_cord_airfoil, color='r')
        # plots sources on the airfoil surface
        #plt.scatter(x_airfoil, y_airfoil,
        #               color='#CD2305', s=20, marker='o')
        # plots max pressure point
        max_pressure_coords = double_index_pair(X, Y, max_pressure)
        plt.scatter(*max_pressure_coords,
                       color='#0022FF', s=80, marker='o')
        plt.show()

    if get_data:
        return u,v
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def VortexFlow(x_vortex, y_vortex, strength_vortex, mesh_properties, get_data=True, plotting=False, width=10):
    '''
    Vortex Flow function.
    
    - uses get_velocity, get_stream_function, mesh_grid Functions
    
    Parameters:
    ----------------
    vortex_strength : positive float, strength of the vortex == gamma
    x_vortex : float, x-coordinate of vortex,
                must in inside the boundaries of meshgrid
    y_vortex : float, y-coordiante of vortex,
                must in inside the boundaries of meshgrid
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u_vortex, v_vortex, psi_vortex
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    --------------
    vortex streamline plot
    u_vortex, v_vortex : velocity field
    psi_vortex: stream-function
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # compute the velocity field
    u_vortex = +strength_vortex / (2 *np.pi) * (Y - y_vortex) / ((X - x_vortex)**2 + (Y - y_vortex)**2)
    v_vortex = -strength_vortex / (2 * np.pi) * (X - x_vortex) / ((X - x_vortex)**2 + (Y - y_vortex)**2)

    # compute the stream-function
    psi_vortex = strength_vortex / (4 * np.pi) * np.log((X - x_vortex)**2 + (Y - y_vortex)**2)
    
    if plotting :
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Vortex Flow")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_vortex, v_vortex,
                          density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(x_vortex, y_vortex, color='#CD2305', s=80, marker='o')
        plt.show()
        
    if get_data:
        return u_vortex, v_vortex, psi_vortex
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Vortex_Sink(x_vortex, y_vortex, x_sink, y_sink, strength_vortex, strength_sink, mesh_properties, get_data=False,  plotting=True, width=10):
    '''
    Vortex + Sink FLow Function
    - usign superposition
    
    Parameters:
    ------------
    strength_vortex : positive float, strength of the vortex = gamma
    strength_sink : positive float, strength of the sink
    x_vortex : float, x-coordinate of vortex,
                must in inside the boundaries of meshgrid
    y_vortex : float, y-coordiante of vortex,
                must in inside the boundaries of meshgrid
    x_sink : float, x-coordinate of sink,
                must in inside the boundaries of meshgrid
    y_sink : float, y-coordiante of sink,
                must in inside the boundaries of meshgrid
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u, v, psi
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    ------------
    plotting : if True returns Velocity Flow field plot for vortex + sink
    get_data : u, v, psi if True
    None
    '''
    
    # getting values form vortex and sink
    u_vortex, v_vortex, psi_vortex = VortexFlow(x_vortex, y_vortex, strength_vortex, mesh_properties, get_data=True, plotting=False)
    u_sink, v_sink, psi_sink = Sink(x_sink, y_sink, strength_sink, mesh_properties)
    
    # getting values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    # superposition of the vortex and sink
    u = u_vortex + u_sink
    v = v_vortex + v_sink
    psi = psi_vortex + psi_sink

    if plotting:
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Vortex Sink/Source")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(x_vortex, y_vortex, color='#CD2305', s=80, marker='o')
        plt.show()

    if get_data:
        return u, v, psi
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Vortex_Source(x_vortex, y_vortex, x_source, y_source, strength_vortex, strength_source, mesh_properties,get_data=False,  plotting=True, width=10):
    """
    Vortex + Source FLow Function
    - usign superposition

    -- uses Vortex_Sink function
    """
    return Vortex_Sink(x_vortex, y_vortex, x_source, y_source, strength_vortex, strength_source, mesh_properties, get_data=get_data,  plotting=plotting, width=width)

#--------------------------------------------------------------------------------------------------------------------------------------------
def VortexSheet(x_array, y_array, strength, mesh_properties, get_data=True, plotting=False, width=10):
    '''
    Calculate the sum of all Source sheets on the Airfoil surface
    
    Paramters:
    ------------
    x_airfoil : np.array, x-coordinates of the airfoil
    y_airfoil : np.array, y-coordinates of the airfoil
    sigma_airfoil : np.array, SourceStrength at each coordiantes of airfoil
    mesh_properties : tuple, grid properties
    
    Returns:
    ---------
    u_source : np.array, sum of all sources velocity at every location on grid
                in x-direction
    y_source : np.array, sum of all sources velocity at every location on grid
                in y-direction
    '''
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    u_vortexsheet = np.zeros_like(X)
    v_vortexsheet = u_vortexsheet.copy()
    
    for (x, y) in zip(x_array, y_array):
        u_add, v_add, _ = VortexFlow(x, y, strength, mesh_properties, get_data=True, plotting=False)
        u_vortexsheet += u_add
        v_vortexsheet += v_add
    
    if plotting:
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Vortex Sheet")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u_vortexsheet, v_vortexsheet,
                          density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(x_array, y_array, color='#CD2305', s=80, marker='o')
        plt.show()

    if get_data:
        return u_vortexsheet, v_vortexsheet
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Horizontal_Infinite_VortexRow(strength, step, mesh_properties, get_data=True, plotting=False, width=10):
    """
    Returns the velocity field generated by a Infinite row of vortex.
    
    Parameters
    ----------
    strength: float
        Strength of the vortexes row
    step: distance betwween two vortex in infinite row
    mesh_properties: Tuple, grid properties
    width: positive float, width of the plot, default=10.0
    get_data : bool, if True retruns u, v
    plotting : bool, plot the meshgrid, default is True
    
    Returns:
    --------------
    infinite vortex row streamline plot
    u, v : if get_data is True
    """
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    coeff = 2.0* np.pi / step
    shared = strength / (2.0 * step * (np.cosh(Y*coeff) - np.cos(X*coeff)))
    u = +np.sinh(Y * coeff) * shared
    v = -np.sin (X * coeff) * shared
    
    if plotting:
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Horizontal Infinite row of Vortices")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v,
                          density=2, linewidth=1, arrowsize=1, arrowstyle='->')
        # plt.scatter(x_array, y_array, color='#CD2305', s=80, marker='o')
        plt.show()
    
    if get_data:
        return u, v
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------
def Lifting_flow_over_Cylinder(x_doublet, y_doublet, x_vortex, y_vortex, strength_doublet, strength_vortex, u_inf, mesh_properties, get_data=True, plotting=False, pressure_coefficient_plot=True, width=10):
    """    
    Lifting flow over a rotating cylinder
    - Vortex Lift
    """


    # get velocities from freestream, doublet and vortex
    u_freestream, v_freestream, psi_freestream = FreestreamFlow(u_inf, mesh_properties)
    u_doublet, v_doublet, psi_doublet = Doublet(strength_doublet, x_doublet, y_doublet, mesh_properties, get_data=True, plotting=False)
    u_vortex, v_vortex, psi_vortex = VortexFlow(strength_vortex, x_vortex, y_vortex, mesh_properties, get_data=True, plotting=False)
    
    # superposition of the doublet and the vortex on the freestream flow
    u = u_freestream + u_doublet + u_vortex
    v = v_freestream + v_doublet + v_vortex
    psi = psi_freestream + psi_doublet + psi_vortex
    
    # gettng values from mesh_properties
    x_start, x_end, y_start, y_end, X, Y, N = mesh_properties
    
    kappa = strength_doublet   # k = strength of the doublet
    gamma = strength_vortex    # gamma = strength of the vortex
    
    # calculate the cylinder radius
    R = np.sqrt(kappa / (2 *  np.pi * u_inf))

    # calculate the stagnation points
    # for gamma/(4* np.pi*R) < 1 and = 1
    x_stagn1, y_stagn1 = (+ np.sqrt(R**2 - (gamma / (4 *  np.pi * u_inf))**2),
                          -gamma / (4 *  np.pi * u_inf))
    x_stagn2, y_stagn2 = (- np.sqrt(R**2 - (gamma / (4 *  np.pi * u_inf))**2),
                          -gamma / (4 *  np.pi * u_inf))

    if plotting:
        height = (y_end - y_start) / (x_end - x_start) * width
        plt.figure(figsize=(width, height))
        plt.title("Lifting Flow over a Cylinder")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v,
                          density=2, linewidth=1, arrowsize=1.5, arrowstyle='->')
        circle = plt.Circle((0.0, 0.0), radius=R, color='#CD2305', alpha=0.5)
        plt.gca().add_patch(circle)
        plt.scatter(x_vortex, y_vortex, color='#CD2305', s=100, marker='o')
        plt.scatter([x_stagn1, x_stagn2], [y_stagn1, y_stagn2],
                       color='g', s=100, marker='o')
        plt.show()
        
    if pressure_coefficient_plot:
        # calculate the surface tangential velocity on the cylinder
        theta = np.linspace(0.0, 2 * np.pi, 100)
        u_theta = -2 * u_inf * np.sin(theta) - gamma / (2 * np.pi * R)

        # compute the surface pressure coefficient
        cp = 1.0 - (u_theta / u_inf)**2

        # if there was no vortex
        u_theta_no_vortex = -2 * u_inf * np.sin(theta)
        cp_no_vortex = 1.0 - (u_theta_no_vortex / u_inf)**2

        # plot the surface pressure coefficient
        size = 6
        plt.figure(figsize=(size, size))
        plt.grid(True)
        plt.title("Pressure Coefficient(Cp) on Surface of Cylinder with and without Vortex")
        plt.xlabel(r'$\theta$', fontsize=18)
        plt.ylabel('$C_p$', fontsize=18)
        plt.xlim(theta.min(), theta.max())
        plt.plot(theta, cp,
                    label='with vortex', color='#CD2305', linewidth=2, linestyle='-')
        plt.plot(theta, cp_no_vortex,
                    label='without vortex', color='g', linewidth=2, linestyle='-')
        plt.legend(loc='best', prop={'size':16})        
        plt.show()

    if get_data:
        return u, v, psi, u_inf, R
    return None

#--------------------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------------------------------------------
