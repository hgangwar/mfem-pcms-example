// Gmsh project created on Mon Oct 30 10:22:05 2023
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {-1, 0, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, -1, 0, 1.0};
//+
Circle(1) = {4, 1, 2};
//+
Circle(2) = {2, 1, 5};
//+
Circle(3) = {5, 1, 3};
//+
Circle(4) = {3, 1, 4};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 2} {
  Curve{4}; Curve{1}; Curve{2}; Curve{3}; Point{4}; Point{2}; Surface{1}; Point{5}; Point{3}; Layers {10}; Recombine;
}
//+
Physical Surface("top", 43) = {42};
//+
Physical Surface("bottom", 44) = {1};
//+
Physical Surface("curved", 45) = {12, 8, 16, 20};
//+
Physical Volume("fuel", 46) = {1};
