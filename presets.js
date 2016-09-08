/* Model 1: Ring, Border growth */
var RingBorder=new Object({
    preset:"Ring Border",
    
    // geometry
    geometry:"ring",    // geometry, either "ring" or "sphere"
    Ri:2,               // inner radius in mm
    Ro:6,               // outter radius in mm
    th:0.5,             // ring thickness in mm
    d:0.5,              // typical length of a volume elements

    // mechanics
    gamma:0.1,        // dumping factor
    rho:0.0001,        // mass density
    mu:50,            // shear modulus
    K:50,            // bulk modulus
    
    // growth
    growth:"ring border instantaneous", // growth function
    G:2.0,                                // growth
    T:0.0,                                 // duration of growth (in sec). T=0 means instantaneous growth

    // display
    colormap:"normal",    // normal, deformation, etc.
    wireframe:true,
    perspective:false    // set to false for orthographic perspective
});

/* Model 2: Ring, Tangential growth */
var RingTangential=new Object({
    preset:"Ring Tangential",
    
    // geometry
    geometry:"ring",
    Ri:2,                        // inner radius in mm
    Ro:6,                        // outter radius in mm
    th:0.5,                          // ring thickness in mm
    d:0.5,                       // typical length of a volume elements

    // mechanics
    gamma:0.1,
    rho:0.0001,                     // mass density
    mu:5,                         // shear modulus
    K:50,                         // bulk modulus
    
    // growth
    growth:"ring tangential instantaneous",
    G:2,                        // growth factor
    T:0.0,                         // duration of growth (in sec)

    // display
    colormap:"deformation",
    wireframe:false,
    perspective:true    // set to false for orthographic perspective
});

/* Model 3: Sphere, Surface growth */
var SphereSurface=new Object({
    preset:"Sphere Surface",
    
    // geometry
    geometry:"sphere",
    url:"data/sphere-1115.json", // surface mesh URL
    th:0.3,                          // surface thickness in mm
    fibres:true,
    
    collision: true,
    Kfc: 10000,      // Collision response string stiffness
    
    // mechanics
    gamma:0.1,                        // dumping
    rho:0.0001,                        // mass density
    mu:65,                            // shear modulus of tetrahedra
    K:100,                            // bulk modulus of tetrahedra
    Kf:1,                            // elastic constant for fibres

    // growth
    growth:"surface homogeneous instantaneous",    // growth function
    G:1.6,                            // growth
    T:0.0,                             // duration of growth (in sec). T=0 means instantaneous growth
    
    // display
    colormap:"normal",
    wireframe:false,
    perspective:false,
    showVertexNumbers: false,
    showTriangleNumbers: false,
    surfaceOnly: true
});

/* Model 4: Sphere, Surface growth, finer mesh */
var SphereSurfaceFine=new Object({
    preset:"Sphere Surface Fine",
    
    // geometry
    geometry:"sphere",
    url:"data/sphere-2500.json", // surface mesh URL
    th:0.2,                          // surface thickness in mm
    fibres:true,
    
    // mechanics
    gamma:0.1,                        // dumping
    rho:0.0001,                        // mass density
    mu:65,                            // shear modulus of tetrahedra
    K:100,                            // bulk modulus of tetrahedra
    Kf:1,                            // elastic constant for fibres

    // growth
    growth:"surface homogeneous instantaneous",    // growth function
    G:1.6,                            // growth
    T:0.0,                             // duration of growth (in sec). T=0 means instantaneous growth
    
    // display
    colormap:"normal",
    wireframe:false,
    perspective:false
});
/* Model 5: Ellipsoid, Surface growth */
var EllipsoidSurface=new Object({
    preset:"Ellipsoid Surface",
    
    // geometry
    geometry:"ellipsoid",
    url:"data/ellipsoid-2582.json",     // surface mesh URL
    th:0.2,                          // surface thickness in mm
    fibres:true,
    
    // mechanics
    gamma:0.1,                         // dumping
    rho:0.0001,                         // mass density
    mu:65,                             // shear modulus of tetrahedra
    K:100,                             // bulk modulus of tetrahedra
    Kf:1,                             // elastic constant for fibres
    
    // growth
    growth:"surface homogeneous instantaneous",    // growth function
    G:1.6,                             // growth
    T:0.0,                             // duration of growth (in sec). T=0 means instantaneous growth
    
    // display
    colormap:"normal",
    wireframe:false,
    perspective:false
});

/* Model 6: Block, Border growth */
var BlockBorder=new Object({
    preset:"Block Border",
    
    // geometry
    geometry:"block",    // geometry, "block", "ring" or "sphere"
    Width:8,               // width in mm
    Height:3,           // height in mm
    Depth:0.5,             // depth in mm
    d:0.5,              // typical length of a volume elements

    // mechanics
    gamma:0.1,        // dumping factor
    rho:0.0001,        // mass density
    mu:50,            // shear modulus
    K:5,            // bulk modulus
    
    // growth
    growth:"block border instantaneous", // growth function
    G:3.0,                                // growth
    T:0.0,                                 // duration of growth (in sec). T=0 means instantaneous growth

    // display
    colormap:"deformation",    // normal, deformation, etc.
    wireframe:false,
    perspective:false    // set to false for orthographic perspective
});

/* Model 7: U Block, No growth, linear spring forces collision */
var UBlockCollision=new Object({
    preset:"U Block Collision",
    
    // geometry
    geometry:"ublock",  // geometry, "block", "ring" or "sphere"
    Width:10,           // width in mm
    Height:1,           // height in mm
    Depth:1,            // depth in mm
    d:1,                // typical length of a volume elements

    // mechanics
    gamma:0.1,      // dumping factor
    rho:0.0001,     // mass density
    mu:50,          // shear modulus
    K:5,            // bulk modulus
    
    collision: true,    // enable collision detection
    Kfc: 10000,      // Collision response string stiffness
    
    // growth
    growth:"homogeneous", // growth function
    G:3.0,                // growth
    T:0.0,                // duration of growth (in sec). T=0 means instantaneous growth

    // display
    colormap:"normal",  // normal, deformation, etc.
    wireframe:true,
    perspective:true,    // set to false for orthographic perspective
    showVertexNumbers: true,
    showTriangleNumbers: false,
    surfaceOnly: true
});



