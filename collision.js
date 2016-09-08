/**
 * @class collision
 * @brief Collision detection algorithm
 */
 
var Collision = {
    h: 0,          // cell size
    nhash: 9973,   // number of cells allocated for collision detection (a prime number close to 10k)
    hash: [],      // collision detection hash table
    resp: [],      // vector of collision responses
    epsilon: 0     // tolerance
}

/**
 * @function HashCell
 * @memberof collision
 * @description Return a hash cell timestamped, storing a vertex index
 * @param {integer} timeStamp Current iteration number
 * @param {integer} i Vertex index
 */
function HashCell(timeStamp,i) {
    return {
        timeStamp: timeStamp, // time stamp for the cell creation
        i: i,                 // vertex index
        next: null,           // next cell at the same hash
    }
}

/**
 * @function initHash
 * @memberof collision
 * @description Init hash array
 * @param {Geometry} ge Geometry of the model
 */
function initHash(ge) {
    var i,i0,i1,i2;
    
    // Init hash array
    Collision.hash=[];
    for(i=0;i<Collision.nhash;i++)
        Collision.hash.push(HashCell(0,0));
        
    // Set grid size to average edge length
    // Each surface triangle supports 3 tetrahedra, one side of the 2nd
    // tetrahedron corresponds to a surface triangle. This is the
    // triangle used for computing the average edge length
    var avrgEdgeLength=0;
    for(i=1;i<ge.nt;i+=3)
    {
        i0=ge.t[4*i+0];
        i1=ge.t[4*i+1];
        i2=ge.t[4*i+2];
        avrgEdgeLength+=norm(subtract([ ge.p[3*i0+0],ge.p[3*i0+1],ge.p[3*i0+2]],
                                      [ ge.p[3*i1+0],ge.p[3*i1+1],ge.p[3*i1+2]]));
        avrgEdgeLength+=norm(subtract([ ge.p[3*i1+0],ge.p[3*i1+1],ge.p[3*i1+2]],
                                      [ ge.p[3*i2+0],ge.p[3*i2+1],ge.p[3*i2+2]]));
        avrgEdgeLength+=norm(subtract([ ge.p[3*i2+0],ge.p[3*i2+1],ge.p[3*i2+2]],
                                      [ ge.p[3*i0+0],ge.p[3*i0+1],ge.p[3*i0+2]]));
    }
    avrgEdgeLength/=ge.nt; // because there are 3 measurements per triangle
    Collision.h=avrgEdgeLength;
    Collision.epsilon=Collision.h/20;
    console.log("average edge length:",avrgEdgeLength);
}

/**
 * @function vertexInTetra
 * @memberof collision
 * @description Test if vertex is inside tetrahedron
 * @param {Geometry} ge Geometry of the model
 * @param {Vector} v Vertex to test
 * @param {Vector Array} T Array with for vectors representing the tetrahedron vertices
 * @param {Vector} penetration If there is penetration, the penetration vector
 * @param {float} epsilon Tolerance
 * @return {boolean} True if a penetration was detected
 */
function vertexInTetra(v, T, penetration, epsilon) {
    var    coords;  // double3D
    var    a,b,c,x; // double3D
    
    a=subtract(T[1],T[0]);
    b=subtract(T[2],T[0]);
    c=subtract(T[3],T[0]);
    x=subtract(v,T[0]);
    
    coords=mulInvMatVec(transpose(vecs2Mat(a, b, c)),x);
    penetration.x=coords[0];
    penetration.y=coords[1];
    penetration.z=coords[2];
    if(coords[0]>=-epsilon && coords[1]>=-epsilon && coords[2]>=-epsilon &&
       coords[0]+coords[1]+coords[2]<=1+epsilon )
        return true;
    return false;
}

/**
 * @function closestPtPointTriangle
 * @memberof collision
 * @description Given a point find the closest point in a triangle. From Ericson (2005) Real-time collision detection, p. 141
 * @param {Vector} p Point to test
 * @param {Vector} a 1st vertex of the triangle
 * @param {Vector} b 2nd vertex of the triangle
 * @param {Vector} c 3rd vertex of the triangle
 * @return {Vector} Closest point
 */
function closestPtPointTriangle(p, a, b, c) {
    // Check if P in vertex region outside A
    var ab = subtract(b, a);
    var ac = subtract(c, a);
    var ap = subtract(p, a);
    var d1 = dot(ab, ap);
    var d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0)
        return {q:a, a:1, b:0, c:0}; // barycentric coordinates (1,0,0)
    
    // Check if P in vertex region outside B
    var bp = subtract(p, b);
    var d3 = dot(ab, bp);
    var d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3)
        return {q:b,a:0,b:1,c:0}; // barycentric coordinates (0,1,0)
    
    // Check if P in edge region of AB, if so return projection of P onto AB
    var vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        var v = d1 / (d1 - d3);
        return {q:add(a, scale(ab,v)),a:1-v,b:v,c:0}; // barycentric coordinates (1-v,v,0)
    }
    
    // Check if P in vertex region outside C
    var cp = subtract(p, c);
    var d5 = dot(ab, cp);
    var d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6)
        return {q:c,a:0,b:0,c:1}; // barycentric coordinates (0,0,1)
    
    // Check if P in edge region of AC, if so return projection of P onto AC
    var vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        var w = d2 / (d2 - d6);
        return {q:add(a, scale(ac,w)),a:1-w,b:0,c:w}; // barycentric coordinates (1-w,0,w)
    }
    
    // Check if P in edge region of BC, if so return projection of P onto BC
    var va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        var w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {q:add(b, scale(subtract(c, b), w)),a:0,b:1-w,c:w}; // barycentric coordinates (0,1-w,w)
    }
    
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    var denom = 1 / (va + vb + vc);
    var v = vb * denom;
    var w = vc * denom;
    return {q:add(add(a, scale(ab,v)), scale(ac, w)),a:1-v-w,b:v,c:w}; //=u*a+v*b+w*c,u=va*denom=1-v-w
}

/**
 * @function closestPtPointTetrahedron
 * @memberof collision
 * @description Given a point find the closest point in a tetrahedron. From Ericson (2005) Real-time collision detection, p. 143
 * @param {Vector} p Point to test
 * @param {Vector} a 1st vertex of the tetrahedron
 * @param {Vector} b 2nd vertex of the tetrahedron
 * @param {Vector} c 3rd vertex of the tetrahedron
 * @param {Vector} d 4th vertex of the tetrahedron
 * @return {Vector} Closest point
 */
function closestPtPointTetrahedron(p, a, b, c, d) {
    // Start out assuming point inside all halfspaces, so closest to itself
    var closestPt = p;
    //var bestSqDist = FLT_MAX;
    var bestSqDist;
    // If point outside face abc then compute closest point on abc
    if (pointOutsideOfPlane(p, a, b, c)) {
        var q = closestPtPointTriangle(p, a, b, c);
        var sqDist = dot(subtract(q.q, p), subtract(q.q, p));
        // Update best closest point if (squared) distance is less than current best
//         if (sqDist < bestSqDist) {
//             bestSqDist = sqDist;
//             closestPt = q;
//         }
        bestSqDist = sqDist;
        closestPt = q;
    }
    // Repeat test for face acd
    if (pointOutsideOfPlane(p, a, c, d)) {
        var q = closestPtPointTriangle(p, a, c, d);
        var sqDist = dot(subtract(q.q, p), subtract(q.q, p));
        if (sqDist < bestSqDist) {
            bestSqDist = sqDist;
            closestPt = q.q;
        }
    }
    // Repeat test for face adb
    if (pointOutsideOfPlane(p, a, d, b)) {
        var q = closestPtPointTriangle(p, a, d, b);
        var sqDist = dot(subtract(q.q, p), subtract(q.q, p));
        if (sqDist < bestSqDist) {
            bestSqDist = sqDist;
            closestPt = q.q;
        }
    }
    // Repeat test for face bdc
    if (pointOutsideOfPlane(p, b, d, c)) {
        var q = closestPtPointTriangle(p, b, d, c);
        var sqDist = dot(subtract(q.q, p), subtract(q.q, p));
        if (sqDist < bestSqDist) {
            bestSqDist = sqDist;
            closestPt = q.q;
        }
    }
    return {p:closestPt,d:bestSqDist};
}
/**
 * @function pointOutsideOfPlane_unknownWinding
 * @memberof collision
 * @description Test if point p and d lie on opposite sides of plane through abc (winding unknown)
 * @param {Vector} p The vertex to test
 * @param {Vector} a 1st vertex of the triangle
 * @param {Vector} b 2nd vertex of the triangle
 * @param {Vector} c 3rd vertex of the triangle
 * @param {Vector} d A vertex on the side of the triangle opposite to p
 * @return {boolean} True if the vertex is outside the triangle's plane
 */
function pointOutsideOfPlane_unknownWinding(p, a, b, c, d)
{
    var signp=dot(subtract(p,a),cross(subtract(b,a),subtract(c,a))); //[APABAC]
    var signd=dot(subtract(d,a),cross(subtract(b,a),subtract(c,a))); //[ADABAC]
    // Points on opposite sides if expression signs are opposite
    return signp * signd < 0;
}

/**
 * @function pointOutsideOfPlane
 * @memberof collision
 * @description Test if vertex p lies outside plane through abc (winding known)
 * @param {Vector} p The vertex to test
 * @param {Vector} a 1st vertex of the triangle
 * @param {Vector} b 2nd vertex of the triangle
 * @param {Vector} c 3rd vertex of the triangle
 * @return {boolean} True if the vertex is outside the triangle's plane
 */
function pointOutsideOfPlane(p, a, b, c) {
    return dot(subtract(p,a),cross(subtract(b,a),subtract(c,a)))>=0; //[APABAC]>=0
}

/**
 * @function addToHash
 * @memberof collision
 * @description Hash table 'add' function for collision detection
 * @param {integer} hash Hash number obtained from the x, y, z coordinates of the vertex
 * @param {integer} vertexIndex Index of the vertex being added
 * @param {integer} iter Iteration number
 */
function addToHash(hash, vertexIndex, iter) {
    var h = Collision.hash[hash];
    
    if(h.timeStamp) {
        // find a free cell, or add a new one
        while(h.timeStamp === iter)
        {
            if(h.next === null)
                h.next = HashCell(0,0);
            h = h.next;
        }
    }
    h.timeStamp=iter;
    h.i=vertexIndex;
}

/**
 * @function hashXYZ
 * @memberof collision
 * @description Generates a hash corresponding to x, y, z coordinates
 * @param {float} x X coordinate
 * @param {float} y Y coordinate
 * @param {float} z Z coordinate
 */
function hashXYZ(x,y,z) {
    return Math.abs(
                ((x*92837111)|0) ^
                ((y*689287499)|0) ^
                ((z*283923481)|0)
               )% Collision.nhash;
}

/**
 * @function collision_surf
 * @memberof collision
 * @description Test for self-collisions in an objects' surface
 * @param {Geometry} ge Geometry object
 * @param {integer} iter Iteration number
 */
function collision_surf_addConstraints(ge, iter) {
    var    i,j,k,x,y,z;
    var    h,hash;
    var    min,max;
    var    belongsToTriangle;
    var    p,q,T=[];
    
    // 1st pass: assign vertices to hash
    for(i=0;i<ge.np;i++) {     // np/2 because only the vertices of the external surface are used
        hash = hashXYZ(Math.floor(ge.p[3*i+0]/Collision.h),Math.floor(ge.p[3*i+1]/Collision.h),Math.floor(ge.p[3*i+2]/Collision.h));
        // console.log(hash);
        addToHash(hash,i,iter);
    }

    $("#txt .collide").removeClass("collide");
    // 2nd pass: detect collision with surface triangle bounding box
    for(i=0;i<ge.nf;i++) {      
        T=[
            [ge.p[3*ge.f[3*i+0]+0],ge.p[3*ge.f[3*i+0]+1],ge.p[3*ge.f[3*i+0]+2]],
            [ge.p[3*ge.f[3*i+1]+0],ge.p[3*ge.f[3*i+1]+1],ge.p[3*ge.f[3*i+1]+2]],
            [ge.p[3*ge.f[3*i+2]+0],ge.p[3*ge.f[3*i+2]+1],ge.p[3*ge.f[3*i+2]+2]]
        ];  
        // bounding box for t
        min=[T[0][0],T[0][1],T[0][2]];
        max=[min[0],min[1],min[2]];
        for(j=1;j<3;j++)
        {
            if(T[j][0]<min[0]) min[0]=T[j][0];
            if(T[j][1]<min[1]) min[1]=T[j][1];
            if(T[j][2]<min[2]) min[2]=T[j][2];
            
            if(T[j][0]>max[0]) max[0]=T[j][0];
            if(T[j][1]>max[1]) max[1]=T[j][1];
            if(T[j][2]>max[2]) max[2]=T[j][2];
        }
        
        // scan bounding box
        for(x=Math.floor(min[0]/Collision.h);x<=Math.floor(max[0]/Collision.h);x++)
        for(y=Math.floor(min[1]/Collision.h);y<=Math.floor(max[1]/Collision.h);y++)
        for(z=Math.floor(min[2]/Collision.h);z<=Math.floor(max[2]/Collision.h);z++) {
            hash = hashXYZ(x,y,z);
            h=Collision.hash[hash];
        
            // check for collision
            while(h.timeStamp) {
                if(h.timeStamp === iter) {
                    // skip if there's already a vertex in collision
                    if(Collision.resp[h.i])
                        break;

                    // check if the vertex belongs to tetrahedron
                    belongsToTriangle=false;
                    for(k=0;k<3;k++) {
                        if(ge.f[3*i+k] === h.i) {
                            belongsToTriangle = true;
                            break;
                        }
                    }
                    if(belongsToTriangle === false)
                    {
                        // 3rd pass: check distance to triangle
                        p=[ge.p[3*h.i+0],ge.p[3*h.i+1],ge.p[3*h.i+2]];
                        q=closestPtPointTriangle(p,T[0],T[1],T[2]);
                        if(norm(subtract(p,q.q))<Collision.epsilon) {
                            Collision.resp[h.i] = {
                                pIndex: h.i,
                                qfIndex: i,
                                a:q.a,
                                b:q.b,
                                c:q.c
                            };
                            $("#v"+h.i).addClass("collide");
                            $("#t"+i).addClass("collide");
                        }
                    }
                }
                if(h.next)
                    h=h.next;
                else
                    break;
            }
        }
    }

    return false;
}

/**
 * @function collision_surf_removeConstraint
 * @memberof collision
 * @description Determine if a collision is resolved and remove unnecessary constrains
 * @param {Geometry} ge Geometry object
 * @param {integer} iter Iteration number
 */
function collision_surf_removeConstraints(ge, iter) {
    var i;
    var p,q,pq;
    var a, b, c;
    var dist;
    for(ci in Collision.resp) {
        p=vertex(ge,Collision.resp[ci].pIndex);
        T=triangleVertices(ge,Collision.resp[ci].qfIndex);
        a=Collision.resp[ci].a;
        b=Collision.resp[ci].b;
        c=Collision.resp[ci].c;
        q=add(add(scale(T[0],a),scale(T[1],b)),scale(T[2],c));
        
        pq=scale(add(p,q),0.5);
        
        if(vertexOutsideTetras(pq,simulation.ge)) {
            display.scene.remove(Collision.resp[ci].l);
            delete Collision.resp[ci];
            console.log("remove line");
        }
    }
}

/**
 * @function collision_surf_enforceConstraints
 * @memberof collision
 * @description Enforce collision constraints as liner springs
 */
function collision_surf_enforceConstraints(ge, me, iter) {
    var Force=me.Force;
	var Kfc=me.Kfc;

	var	l;
	
    var ci, pi;
    var p,q,pq,T;
    var a, b, c;
    var dist;

    for(ci in Collision.resp) {
        pi=Collision.resp[ci].pIndex;
        p=vertex(ge,pi);
        T=triangleVertices(ge,Collision.resp[ci].qfIndex);
        a=Collision.resp[ci].a;
        b=Collision.resp[ci].b;
        c=Collision.resp[ci].c;
        q=add(add(scale(T[0],a),scale(T[1],b)),scale(T[2],c));
        
        pq=subtract(p,q);
        l=norm(pq);
        
        Force[3*pi+0]+=-Kfc*(pq[0])*l;
        Force[3*pi+1]+=-Kfc*(pq[1])*l;
        Force[3*pi+2]+=-Kfc*(pq[2])*l;
    }
}

/**
 * @function collision_tetra
 * @memberof collision
 * @description Test for self-collisions in an objects' volume
 * @param {Geometry} ge Geometry object
 * @param {integer} iter Iteration number
 */
function collision_tetra(ge, iter) {
    var    i,j,k,x,y,z;
    var    h,hash;
    var    min,max,penetration;
    var    isInTetrahedron;
    var    p,T=[];
    var    penetration={};
    
    // 1st pass: assign vertices to hash
    for(i=0;i<ge.np;i++) {     // np/2 because only the vertices of the external surface are used
        hash = hashXYZ(Math.floor(ge.p[3*i+0]/Collision.h),Math.floor(ge.p[3*i+1]/Collision.h),Math.floor(ge.p[3*i+2]/Collision.h));
        // console.log(hash);
        addToHash(hash,i,iter);
    }

    $("#txt .collide").removeClass("collide");
    // 2nd pass: detect collision with tetrahedra bounding box
    for(i=0;i<ge.nt;i++) {      
        T=[
            [ge.p[3*ge.t[4*i+0]+0],ge.p[3*ge.t[4*i+0]+1],ge.p[3*ge.t[4*i+0]+2]],
            [ge.p[3*ge.t[4*i+1]+0],ge.p[3*ge.t[4*i+1]+1],ge.p[3*ge.t[4*i+1]+2]],
            [ge.p[3*ge.t[4*i+2]+0],ge.p[3*ge.t[4*i+2]+1],ge.p[3*ge.t[4*i+2]+2]],
            [ge.p[3*ge.t[4*i+3]+0],ge.p[3*ge.t[4*i+3]+1],ge.p[3*ge.t[4*i+3]+2]]
        ];  
        // bounding box for t
        min=[T[0][0],T[0][1],T[0][2]];
        max=[min[0],min[1],min[2]];
        for(j=1;j<4;j++)
        {
            if(T[j][0]<min[0]) min[0]=T[j][0];
            if(T[j][1]<min[1]) min[1]=T[j][1];
            if(T[j][2]<min[2]) min[2]=T[j][2];
            
            if(T[j][0]>max[0]) max[0]=T[j][0];
            if(T[j][1]>max[1]) max[1]=T[j][1];
            if(T[j][2]>max[2]) max[2]=T[j][2];
        }
        
        // scan bounding box
        for(x=Math.floor(min[0]/Collision.h);x<=Math.floor(max[0]/Collision.h);x++)
        for(y=Math.floor(min[1]/Collision.h);y<=Math.floor(max[1]/Collision.h);y++)
        for(z=Math.floor(min[2]/Collision.h);z<=Math.floor(max[2]/Collision.h);z++) {
            hash = hashXYZ(x,y,z);
            
            // check for collision
            h=Collision.hash[hash];
            while(h.timeStamp) {
                if(h.timeStamp === iter) {
                    // check if the vertex belongs to tetrahedron
                    isInTetrahedron=false;
                    for(k=0;k<4;k++) {
                        if(ge.t[4*i+k] === h.i) {
                            isInTetrahedron = true;
                            break;
                        }
                    }
                    if(isInTetrahedron === false)
                    {
                        // 3rd pass: check distance to tetrahedron faces
                        p=[ge.p[3*h.i+0],ge.p[3*h.i+1],ge.p[3*h.i+2]];
                        //q=closestPtPointTetrahedron(p,T[0],T[1],T[2],T[3]);
                        //if(q.d<Collision.epsilon) {
                        if(vertexInTetra(vertex(ge,h.i),tetraVertices(ge,i),penetration, Collision.epsilon)) {
                            $("#v"+h.i).addClass("collide");
                            
                            //console.log("Collision. v:",h.i,"t:",i,"p:",penetration.x,penetration.y,penetration.z);
                            if(0) {}
                        }
                    }
                }
                if(h.next)
                    h=h.next;
                else
                    break;
            }
        }
    }

    return false;
}
