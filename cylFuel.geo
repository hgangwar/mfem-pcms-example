/*   Cylindirical Fuel Pin Geometry and Mesh Generator in 3D in gmsh format   */

/*        Description of the geometry and meshing parameters                  */
/*       ----------------------------------------------------                 
/*       |Length of the fuel pin = 1.0 cm                       |
/*       |Radius of the fuel pin = 0.5 cm                       |
/*       |Number of elements in the axial direction = 100       |
/*       |Number of elements in the radial direction = 50       |

/*       |There will be 3 surfaces in the geometry               |
/*       |  Top, Bottom and the curved surface                   |

/*       |It will be a 3D volume mesh                           |


/*       |The mesh will be saved in the file "fuel_pin.msh"     |
/*       ----------------------------------------------------   */


/*       --------------The code starts here---------------------   */
// Enable OpenCASCADE kernel
SetFactory("OpenCASCADE");

// 1. Define the parameters of the geometry

// 1.1. Define the length of the fuel pin
L = 100.0;

// 1.2. Define the radius of the fuel pin
R = 50;

// 1.3. Define the number of elements in the axial direction
Nz = 100;

// 1.4. Define the number of elements in the radial direction
Nr = 50;

// 2. Define the geometry
// Create a disk with radius R in the x-y plane and name it as "Bottom" and
// use it to get the curved surface of the fuel pin and name it as "Curved"
// and use it to get the top surface of the fuel pin and name it as "Top"
// and use it to get the volume of the fuel pin and name it as "Fuel"

// 2.1. Create the bottom surface of the fuel pin
Disk(1) = {0, 0, 0, R};

// 2.2. create the top surface of the fuel pin
Disk(2) = {0, 0, L, R};

// 2.3. Create the curved surface of the fuel pin
Cylinder(3) = {0, 0, 0, 0, 0, L, R};

// 2.4. Create the volume of the fuel pin
Surface Loop(4) = {1, 2, 3};
Volume(5) = {4};

// 3. Define the meshing parameters
MeshSize {1, 2, 3, 4} = 1.0;

// 4. Create physical entities
Physical Volume("fuel_volume", 1) = {5};
Physical Surface("top", 2) = {2};
Physical Surface("bottom", 3) = {1};
Physical Surface("curved", 4) = {3};


