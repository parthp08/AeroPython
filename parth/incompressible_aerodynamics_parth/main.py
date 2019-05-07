from src.IncompressibleAerodynamics import *

## to printout all the avilable functions in Aerodynamics library
Avilable_Functions()

## NACA 4-digit airfoil plotter
NACA4digit('0012',points=1000, get_data=False, plotting=True)

## Define Mesh Grid
## change mesh_properties as per requirement
## Do-not comment out " mesh_properties " as all below functions depend on it
mesh_properties = mesh_grid(-4.0,4.0,-2.0,2.0,200, grid_properties=True, plotting=True, width=20)

## Source Flow
Source(-1.0, 0.0, 5.0, mesh_properties, get_data=False, plotting=True, width=20)
## change source_strength and width to observe diffrent results

## Sink Flow
Sink(-1.0, 0.0, -5.0, mesh_properties, get_data=False, plotting=True, width=20)

## Source + Sink pair
Source_Sink_Pair(-1.0,0.0,1.0,0.0,5.0,mesh_properties,get_data=False,plotting=True, width=20)

## Source in Freestream flow
Source_plus_Freestream(-1.0, 0.0, 10.0, 2.0, mesh_properties, get_data=False, plotting=True, width=20)

## Source-Sink pair in Freestream flow
SourceSink_plus_Freestream(-1.0,0.0,1.0,0.0,5.0,-5.0,1.0,mesh_properties,get_data=False, plotting=True, pressure_field_plot=True, width=10)

## Doublet Flow
Doublet(1.0, 0.0, 0.0, mesh_properties, get_data=False, plotting=True, width=20)

## Doublet in Freestream flow
Doublet_plus_Freestream(1.0, 0.0, 0.0, 1.0, mesh_properties, get_data=False, plotting=True, width=20)

## Non lifting flow over a cylinder   # non rotating cylinder
NonLifting_flow_over_Cylinder(1.0, 0.0, 0.0, 1.0, mesh_properties, get_data=False, plotting=True, pressure_field_plot=True, width=20)

## Non-lifting flow over NACA0012 airfoil
## loading NACA0012 data from files
mesh_properties2 = mesh_grid(-1.0,2.0,-0.5,0.5,200, grid_properties=True)
x_airfoil = np.loadtxt('resources/NACA0012_x.txt')  # X-corrdinates of the airfoil
y_airfoil = np.loadtxt('resources/NACA0012_y.txt')  # Y-corrdinates of the airfoil
sigma_airfoil = np.loadtxt('resources/NACA0012_sigma.txt')  # strength of the sources located at the surfaces
NonLifting_flow_over_NACA0012(x_airfoil, y_airfoil, sigma_airfoil, 1.0, mesh_properties2, get_data=False, plotting=True, pressure_field_plot=True, width=20)

## Vortex flow
VortexFlow(0.0, 0.0, 5.0, mesh_properties, get_data=False, plotting=True, width=20)

## Vortex Sink & Vortex Source
Vortex_Sink(0.0, 0.0, 0.0, 0.0, 5.0, -1.0, mesh_properties, get_data=False, plotting=True, width=20)
Vortex_Source(0.0, 0.0, 0.0, 0.0, 5.0, 1.0, mesh_properties, get_data=False, plotting=True)

## Vortex Sheet
mesh_properties3 = mesh_grid(-2.0,2.0,-0.5,0.5,50, grid_properties=True)
x_array = np.linspace(-2.0*16,2.0*16,num=300)
y_array = np.zeros(300,dtype=float)
# strength_array = np.zeros(11, dtype=float) + 10
VortexSheet(x_array, y_array, 5.0, mesh_properties3, get_data=False, plotting=True, width=20)

## Horizontal infinite row of vortcies
mesh_properties3 = mesh_grid(-2.0,2.0,-0.5,0.5,50, grid_properties=True)
Horizontal_Infinite_VortexRow(5.0, 0.01, mesh_properties3, get_data=False, plotting=True, width=20)

## Lifting flow over a cyliner   # rotating cylinder
Lifting_flow_over_Cylinder(0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 1.0, mesh_properties, get_data=False, plotting=True, pressure_coefficient_plot=True, width=20)

