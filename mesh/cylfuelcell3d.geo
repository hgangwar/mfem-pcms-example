/*
* This file generates the fuel pin geometry for the 3D full length fuel pin
* Geometry parameters
@param r radius of the fuel pin = 0.74 cm
@param l length of the fuel pin = 140 cm
*/
r = 0.74;
l = 140;

// * Mesh parameters
cl = 0.2;

// * The points

// * The center point
Point(1) = {0, 0, 0, cl};

// * The perpheral points
Point(2) = {r, 0, 0, cl};
Point(3) = {-r, 0, 0, cl};
Point(4) = {0, r, 0, cl};
Point(5) = {0, -r, 0, cl};


// * The circle arcs
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 5};
Circle(3) = {5, 1, 3};
Circle(4) = {3, 1, 4};


// * The circle
Curve Loop(1) = {4, 1, 2, 3};

// * Create the plane surface using the circle
Plane Surface(1) = {1};

// * The extrusion to create the 3D cylinder
Extrude {0, 0, l} {
  Curve{4}; Curve{1}; Curve{2}; Curve{3}; Point{4}; Point{2}; Surface{1}; Point{5}; Point{3}; Layers {l*4}; 
}

// * The physical entities
// * The surfaces
Physical Surface("top", 1) = {42};
Physical Surface("bottom", 2) = {1};
Physical Surface("curved", 3) = {12, 8, 16, 20};

// * The volumes
Physical Volume("fuel", 4) = {1};

// * mesh version 2.2
Mesh.MshFileVersion = 2.2;

/*
* To create the mesh, open Gmsh
* File -> Open -> select the .geo file
* Click the + button near the Mesh and select 3D
* You will see that the mesh is created
* File -> Export
* Write the file name as yourfilename.msh
* Select the format as Version 2 ASCII
* Click Save
* You can now use this mesh file in MFEM
* Or run the following command
* gmsh -3 -order 1 -format msh2 -o youmeshname.msh youscript.geo
*/

