SetFactory("OpenCASCADE");

// ------------------------------
// Parameters (similar spirit to planar_dam_10x5.geo)
// ------------------------------
Nx   = 150; // nominal cells across X (used to set target size)
Ny   = 150; // nominal cells across Y (used to set target size)

Lx = 5.0; // domain width  (e.g., [-2.5, 2.5] shifted to [0,5])
Ly = 5.0; // domain height (same)

Cx = 2.5; // circle center x (middle of box)
Cy = 2.5; // circle center y
R  = 0.5; // dam radius

// Derived nominal sizes (like dx, dy)
dx = Lx / Nx;
dy = Ly / Ny;

// Target mesh sizes
lc_far  = Min(dx, dy);     // far-field target
lc_near = 0.5 * lc_far;    // finer near the circle

// Refinement band thickness around the circle
ref_band = 0.3; // tweak as desired

// ------------------------------
// Geometry: rectangle with circular inclusion
// ------------------------------
Rectangle(1) = {0, 0, 0, Lx, Ly, 0};
Disk(2)      = {Cx, Cy, 0, R, R};

// Boolean fragment: split outer box by the disk to get two surfaces
// s_out = outer field (box minus disk), s_in = disk
out[] = BooleanFragments{ Surface{1}; Delete; }{ Surface{2}; Delete; };
s_in  = out[1]; // inner disk surface tag (first new surface)
s_out = out[0]; // outer field surface tag (second new surface)

// ------------------------------
// Mesh size control
// ------------------------------
// Default min/max (safety caps)
Mesh.CharacteristicLengthMin = 0.25 * lc_near;
Mesh.CharacteristicLengthMax = 2.0  * lc_far;

// Distance to the circle edges
Field[1] = Distance;
Field[1].SurfacesList = {s_in}; // distance from the disk boundary

// Threshold: near the circle -> lc_near; away -> lc_far
Field[2] = Threshold;
Field[2].IField  = 1;
Field[2].LcMin   = lc_near;
Field[2].LcMax   = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = ref_band;

// Apply background field
Background Field = 2;

// Optional: try to recombine (quads where possible). Around a circle
// quads are not guaranteed; comment out if you prefer pure triangles.
// Mesh.RecombineAll = 1;

// ------------------------------
// Physical groups (stable IDs for RDycore)
// ------------------------------
// Surfaces
Physical Surface("inner_disk", 1) = {s_in};
Physical Surface("outer_field", 2) = {s_out};

// Identify boundary curves by bounding boxes (robust even if edges get split)
eps = 1e-9;

// Left boundary (x = 0)
left_curves[] = Curve In BoundingBox(-eps, -eps, -eps, eps, Ly+eps, eps);
Physical Curve("left", 1) = {left_curves[]};

// Right boundary (x = Lx)
right_curves[] = Curve In BoundingBox(Lx-eps, -eps, -eps, Lx+eps, Ly+eps, eps);
Physical Curve("right", 2) = {right_curves[]};

// Bottom boundary (y = 0)
bottom_curves[] = Curve In BoundingBox(-eps, -eps, -eps, Lx+eps, eps, eps);
Physical Curve("bottom", 3) = {bottom_curves[]};

// Top boundary (y = Ly)
top_curves[] = Curve In BoundingBox(-eps, Ly-eps, -eps, Lx+eps, Ly+eps, eps);
Physical Curve("top", 4) = {top_curves[]};
