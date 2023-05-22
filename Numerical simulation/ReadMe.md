# Numerical simulation
Import the geometry of discrete fiber networks to Abaqus, and analyze their mechanical behaviors.

## Data
Simulation results are saved in Abaqus folder:
- "DFN_output" is the simulation data of the reconstructed discrete fiber network (DFN).
- "HGO output" is the data of the homogenized fiber network using Holzapfel-Gasser-Ogden model (HGO).

### Main scripts

`myNetworkSim.m`: discretize the geometry of the fiber network and create Abaqus input file.

`SimPostprocess.m`: postprocess the simulation results to plot stress-strain curve for analysis.

### Self-defined functions

`createBeamMesh.m`: discretize the fiber segments into equi-distant elements and define the boundary condition.

## Proccedure
### Discretization of the geometry and create Abaqus input file

1. Open `myNetworkSim.m`, and then run the algorithm.
2. Select the folder with the preprocessed geometry of the fiber network at the pop-up window, and then select the folder to save the Abaqus input file.
(Note: a simple fiber network based on Voronoi diagram is saved in "Geometry_test" folder, which can be used to test the algorithm)
3. After the discretization is completed, an Abaqus input file will be created. And also the initial coordinate of the boundary node set will be saved.
(Note: the Abaqus input files of the discrete fiber network are saved in "Input_fiber" for review)

### Postprocessing of the simulation results

1. Open `SimPostprocess.m`, and then run the algorithm.
2. Select the folders of the simulation results subsequently according to the command of the pop-up window.
3. After the postprocessing is completed, the stress-strain curve of the discrete as well as the homogenized fiber network models will be plotted. And the maximum normal and shear stress will be compared.
(Note: the postprocessing results are saved in the folder "Result_sim" for review)
