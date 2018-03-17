/**
 * @page Geometry
 * Geometry functions
 *
 * The myGeometry object stores the material and rest configurations.
 * There is also code to generate ring and brain geometries
 *
 * Depends on algebra.js 
 */

var tetraTopo=[
    "001 010 100 111",
    "000 010 100 001",
    "011 001 111 010",
    "101 111 001 100",
    "110 100 010 111"
];

function myGeometry() {
    // Variables storing geometry and topology
    this.np=0;        // number of material vertices in the model
    this.nt=0;        // number of elastic tetrahedra in the model
    this.p=0;         // material vertex coordinates
    this.r=0;         // rest tetrahedra geometry
    this.re=0;        // rest edge geometry
    this.t=0;         // tetrahedra vertex indices (topology) in the material configuration
    this.Volume=0;    // nodal volumes of the rest configuration
}

function vertex(ge,i) {
    return [ge.p[3*i+0],ge.p[3*i+1],ge.p[3*i+2]];
}
function triangleIndices(ge,i) {
    return [ge.t[3*i+0],ge.t[3*i+1],ge.t[3*i+2]];
}
function triangleVertices(ge,i) {
    return [[ge.p[3*ge.f[3*i+0]+0], ge.p[3*ge.f[3*i+0]+1], ge.p[3*ge.f[3*i+0]+2]],
            [ge.p[3*ge.f[3*i+1]+0], ge.p[3*ge.f[3*i+1]+1], ge.p[3*ge.f[3*i+1]+2]],
            [ge.p[3*ge.f[3*i+2]+0], ge.p[3*ge.f[3*i+2]+1], ge.p[3*ge.f[3*i+2]+2]]];
}

function tetraVertices(ge, i) {
    return [
        [ge.p[3*ge.t[4*i+0]+0],ge.p[3*ge.t[4*i+0]+1],ge.p[3*ge.t[4*i+0]+2]],
        [ge.p[3*ge.t[4*i+1]+0],ge.p[3*ge.t[4*i+1]+1],ge.p[3*ge.t[4*i+1]+2]],
        [ge.p[3*ge.t[4*i+2]+0],ge.p[3*ge.t[4*i+2]+1],ge.p[3*ge.t[4*i+2]+2]],
        [ge.p[3*ge.t[4*i+3]+0],ge.p[3*ge.t[4*i+3]+1],ge.p[3*ge.t[4*i+3]+2]]
    ];  
}

/**
 * @function tetraVol
 * @return the volume of a tetrahedron
 */
function tetraVol(a,b,c, d,e,f, g,h,i, j,k,l)
{
    var m=new Object();
    var vol;
    
    m.a=d-a;
    m.b=e-b;
    m.c=f-c;
    m.d=g-a;
    m.e=h-b;
    m.f=i-c;
    m.g=j-a;
    m.h=k-b;
    m.i=l-c;
    
    vol=determinant(m)/6;
    return vol;
}
/*
allTetraVol
returns the total volume of a tetrahedral mesh
*/
function allTetraVol(c) {
    var i;
    var V=0,vol;
    var min,max;
    for(i=0;i<nt;i++) {
        vol=tetraVol(
            c[3*t[4*i+0]+0],
            c[3*t[4*i+0]+1],
            c[3*t[4*i+0]+2],
            c[3*t[4*i+1]+0],
            c[3*t[4*i+1]+1],
            c[3*t[4*i+1]+2],
            c[3*t[4*i+2]+0],
            c[3*t[4*i+2]+1],
            c[3*t[4*i+2]+2],
            c[3*t[4*i+3]+0],
            c[3*t[4*i+3]+1],
            c[3*t[4*i+3]+2]);
        V+=vol;
        if(i==0) min=max=vol;
        else if(vol<min) min=vol;
        else if(vol>max) max=vol;
    }
    //console.log("vol:"+V+", min:"+min+", max:"+max);
    return V;
}

/**
 * @function vertexInSurface
 * @memberof geometry
 * @description Determine if a vertex is inside a surface
 * @return {boolean} True if the vertex is inside the surface, false if it is not, NaN if it is over the surface
 */
function vertexOutsideSurface(v,ge) {
    var i;
    var T;
    var W=0;
    
    for(i=0;i<ge.nf;i++) {
        T=triangleVertices(ge,i);
        W+=solidAngle(T[0],T[1],T[2],v);
    }
    return !(isNaN(W) || W>Math.PI*4-0.1);
}
function vertexOutsideTetras(p, ge) {
    var    i,j,k,x,y,z;
    var    p,T=[];
    var    penetration={};
    
    for(i=0;i<ge.nt;i++) {      
        if(vertexInTetra(p,tetraVertices(ge,i),penetration,0)==true) {
            return false;
        }
    }
    return true;
}
/*
configureRestGeometry
initialises a rest geometry to the actual geometry. Topology is different, however,
because all tetrahedra are disconnected.
*/
function configureRestGeometry(ge) {
    var nt=ge.nt;
    var r=ge.r;
    var p=ge.p;
    var t=ge.t;
    var i,j,k;
    
    for(i=0;i<nt;i++)
    for(j=0;j<4;j++)
    for(k=0;k<3;k++)
        r[(4*i+j)*3+k]=p[t[4*i+j]*3+k];
}
/*
computeNodalVolume
computes the volume measured at each point of a tetrahedral mesh as the sum of
1/4 of the volume of each of the tetrahedra to which it belongs.
*/
function computeNodalVolume(ge) {
    var np=ge.np;
    var nt=ge.nt;
    var p=ge.p;
    var t=ge.t;
    var Volume=ge.Volume;
    var i,n1,n2,n3,n4,vol;
    
    for(i=0;i<np;i++)
        Volume[i]=0;
    for(i=0;i<nt;i++) {
        n1=t[4*i+0];
        n2=t[4*i+1];
        n3=t[4*i+2];
        n4=t[4*i+3];
        vol=tetraVol(
          p[3*n1+0], p[3*n1+1], p[3*n1+2],
          p[3*n2+0], p[3*n2+1], p[3*n2+2],
          p[3*n3+0], p[3*n3+1], p[3*n3+2],
          p[3*n4+0], p[3*n4+1], p[3*n4+2]);
        Volume[n1]+=vol/4;
        Volume[n2]+=vol/4;
        Volume[n3]+=vol/4;
        Volume[n4]+=vol/4;
    }
}

/*

           Creation of geometries

*/

/*
Get global index of a tetrahedron within a block based on its width, height and depth
indices.
@param {integer} i width index of the tetrahedron
@param {integer} j height index of the tetrahedron
@param {integer} k depth index of the tetrahedron
*/
function blockVind(i,j,k,nw,nh) {
    return k*nh*nw+j*nw+i;
}
/*
Create a new block tetrahedron topology.
@param {integer} n global index of the tetrahedron
@param {integer} i base width index of the tetrahedron
@param {integer} j base height index of the tetrahedron
@param {integer} k base depth index of the tetrahedron
@param {string} s description of tetrahedron topology relative to the base i, j, k indices
*/
function blockTetra(n,i,j,k,s,ge) {
    var t=ge.t;
    var nw=ge.nw;
    var nh=ge.nh;
    var a=s.split(" ");
    var b=a.map(function(m) {
            return blockVind(
                i+parseInt(m.charAt(0)),
                j+parseInt(m.charAt(1)),
                k+parseInt(m.charAt(2)),
                nw,nh);
        });
    t[4*n+0]=b[0];
    t[4*n+1]=b[1];
    t[4*n+2]=b[2];
    t[4*n+3]=b[3];
}
/*
Get global index of a tetrahedron within a ring based on its angular, radial and depth
indices.
@param {integer} i angular index of the tetrahedron
@param {integer} j radial index of the tetrahedron
@param {index} depth index of the tetrahedron
*/
function ringVind(i,j,k,ntheta,nxy) {
    return k*nxy*ntheta+j*ntheta+i;
}
/*
Create a new ring tetrahedron topology.
@param {integer} n global index of the tetrahedron
@param {integer} i base angular index of the tetrahedron
@param {integer} j base radial index of the tetrahedron
@param {integer} k base depth index of the tetrahedron
@param {string} s description of tetrahedron topology relative to the base i, j, k indices
*/
function ringTetra(n,i,j,k,s,ge) {
    var t=ge.t;
    var ntheta=ge.ntheta;
    var nxy=ge.nxy;
    var a=s.split(" ");
    var b=a.map(function(m) {
            return ringVind(
                (i+parseInt(m.charAt(0)))%ntheta,
                j+parseInt(m.charAt(1)),
                k+parseInt(m.charAt(2)),
                ntheta,nxy);
        });
    t[4*n+0]=b[0];
    t[4*n+1]=b[1];
    t[4*n+2]=b[2];
    t[4*n+3]=b[3];
}
/*
makeBlock
*/
function makeBlock(params) {
    var def=$.Deferred();
    var ge = new myGeometry();
    
    var Width=params.Width;
    var Height=params.Height;
    var Depth=params.Depth;
    var d=params.d;
    var nw=parseInt(Width/d+0.5)+1;     // volume elements in the outter circle (+0.5 used for rounding)
    var nh=parseInt(Height/d+0.5)+1;    // number of vol. elem. rings in x-y plane
    var nd=parseInt(Depth/d+0.5)+1;     // number of vol. elem. rings in z
    var np=nw*nh*nd;
    var nt=5*(nw-1)*(nh-1)*(nd-1);
    var p=new Float32Array(np*3);
    var r=new Float32Array(nt*4*3);
    var t=new Uint16Array(nt*4);
    var Volume=new Float32Array(np);

    var i,j,k,l,n,m,R;

    ge.params=params;

    ge.Width=Width;
    ge.Height=Height;
    ge.Depth=Depth;
    ge.d=d;
    
    ge.nw=nw;
    ge.nh=nh;
    ge.nd=nd;
    ge.np=np;
    ge.nt=nt;
    ge.p=p;
    ge.r=r;
    ge.t=t;
    ge.Volume=Volume;
    
    console.log("Number of vertices: "+np+" ("+nw+","+nh+","+nd+")");
    console.log("Number of tetrahedra: "+nt);

    // create material vertices
    m=0;
    for(i=0;i<nw;i++)
    for(j=0;j<nh;j++)
    for(k=0;k<nd;k++)
    {
        n=k*nh*nw+j*nw+i;

        // material configuration
        p[3*n+0]=d*i-Width/2;
        p[3*n+1]=Height/2-d*j;
        p[3*n+2]=d*k-Depth/2;
        m++;
    }

    // create the topology of material tetrahedra
    n=0;
    for(i=0;i<nw-1;i++)
    for(j=0;j<nh-1;j++)
    for(k=0;k<nd-1;k++)
    for(l=0;l<5;l++)
        blockTetra(n++,i,j,k,tetraTopo[l],ge);

    // create rest vertices, copying the positions of material vertices
    // (the topology of rest tetrahedra is implicitly defined by its indices)
    configureRestGeometry(ge);

    // compute nodal volume
    computeNodalVolume(ge);

    return def.resolve(ge).promise();
}

/*
makeUBlock
*/
function makeUBlock(params) {
    var def=$.Deferred();
    var ge = new myGeometry();
    
    var Width=params.Width;
    var Height=params.Height;
    var Depth=params.Depth;
    var d=params.d;
    var nw=parseInt(Width/d+0.5)+1;
    var nh=parseInt(Height/d+0.5)+1;
    var nd=parseInt(Depth/d+0.5)+1;
    var np=nw*nh*nd;                                        // number of vertices
    var nt=5*(nw-1)*(nh-1)*(nd-1);                          // number of tetrahedra
    var nf=2*2*((nw-1)*(nh-1)+(nh-1)*(nd-1)+(nd-1)*(nw-1)); // number of triangular faces
    var p=new Float32Array(np*3);
    var r=new Float32Array(nt*4*3);
    var t=new Uint16Array(nt*4);
    var f=new Uint16Array(nf*3);
    var Volume=new Float32Array(np);

    var i,j,k,l,n,m,R;

    ge.params=params;

    ge.Width=Width;
    ge.Height=Height;
    ge.Depth=Depth;
    ge.d=d;
    
    ge.nw=nw;
    ge.nh=nh;
    ge.nd=nd;
    ge.np=np;
    ge.nt=nt;
    ge.nf=nf;
    ge.p=p;
    ge.r=r;
    ge.t=t;
    ge.f=f;
    ge.Volume=Volume;
    
    console.log("Number of vertices: "+np+" ("+nw+","+nh+","+nd+")");
    console.log("Number of tetrahedra: "+nt);
    console.log("Number of surface triangular faces: "+nf);

    // create material vertices
    m=0;
    for(i=0;i<nw;i++)
    for(j=0;j<nh;j++)
    for(k=0;k<nd;k++)
    {
        n=k*nh*nw+j*nw+i;

        // material configuration
        var rho,theta,R=1.5*Depth,x,y;
        theta=(d*i-Width/2)/R;
        if(Math.abs(theta)<Math.PI/2) {
            rho=R-d*k;
            x=rho*Math.cos(theta-Math.PI/2);
            y=rho*Math.sin(theta-Math.PI/2);
        } else {
            x=(theta<0)?(-R+d*k):(R-d*k);
            y=Math.abs(d*i-Width/2)-R*Math.PI/2;
        }
            
        p[3*n+0]=x;
        p[3*n+1]=Height/2-d*j;
        p[3*n+2]=y;
        m++;
    }
    /*
    console.log("m:",m);
    var arr=[   -1.9656, 1.0426, 4.7928, -2.2754, 1.0312, 2.7593, -2.5671, 1.0056, 0.8072, -2.7500, 0.9858, -1.1377, -1.8487, 0.9829, -2.8714, -0.0112, 0.9827, -3.5792, 1.8343, 0.9845, -2.8835, 2.7498, 0.9884, -1.1447, 2.5666, 0.9941, 0.8232, 2.2652, 1.0128, 2.7756, 1.9600, 1.0276, 4.7863, -1.9600, -1.0276, 4.7863, -2.2652, -1.0128, 2.7756, -2.5665, -0.9941, 0.8232, -2.7498, -0.9884, -1.1447, -1.8343, -0.9845, -2.8835, 0.0112, -0.9826, -3.5792, 1.8487, -0.9829, -2.8714, 2.7500, -0.9858, -1.1377, 2.5671, -1.0056, 0.8072, 2.2754, -1.0312, 2.7593, 1.9656, -1.0426, 4.7928, 0.0789, 1.0381, 4.4661, -0.2578, 1.0213, 2.4593, -0.5999, 0.9979, 0.5392, -0.8158, 0.9817, -0.8860, -0.5570, 0.9767, -1.4228, -0.0109, 0.9747, -1.6418, 0.5368, 0.9765, -1.4340, 0.7980, 0.9812, -0.9063, 0.5859, 0.9964, 0.5114, 0.2458, 1.0210, 2.4418, -0.0859, 1.0320, 4.4710, 0.0859, -1.0320, 4.4710, -0.2458, -1.0210, 2.4418, -0.5859, -0.9964, 0.5114, -0.7980, -0.9812, -0.9063, -0.5368, -0.9765, -1.4340, 0.0109, -0.9747, -1.6418, 0.5570, -0.9767, -1.4228, 0.8158, -0.9817, -0.8860, 0.5999, -0.9979, 0.5392, 0.2578, -1.0213, 2.4593, -0.0789, -1.0381, 4.4661];
    for(i=0;i<132;i++)
        p[i]=arr[i];
    */

    // create the topology of surface triangles
    m=0;
    for(i=0;i<nw-1;i++)
    for(j=0;j<nh-1;j++) {
        f[3*m+0]=j*nw+i;
        f[3*m+1]=j*nw+(i+1);
        f[3*m+2]=(j+1)*nw+(i+1);
        m++;
        f[3*m+0]=j*nw+i;
        f[3*m+1]=(j+1)*nw+(i+1);
        f[3*m+2]=(j+1)*nw+i;
        m++
        f[3*m+0]=(nd-1)*nh*nw+j*nw+i;
        f[3*m+1]=(nd-1)*nh*nw+(j+1)*nw+(i+1);
        f[3*m+2]=(nd-1)*nh*nw+j*nw+(i+1);
        m++;
        f[3*m+0]=(nd-1)*nh*nw+j*nw+i;
        f[3*m+1]=(nd-1)*nh*nw+(j+1)*nw+i;
        f[3*m+2]=(nd-1)*nh*nw+(j+1)*nw+(i+1);
        m++;
    }
    for(j=0;j<nh-1;j++)
    for(k=0;k<nd-1;k++) {
        f[3*m+0]=k*nh*nw+j*nw;
        f[3*m+1]=k*nh*nw+(j+1)*nw;
        f[3*m+2]=(k+1)*nh*nw+(j+1)*nw;
        m++;
        f[3*m+0]=k*nh*nw+j*nw;
        f[3*m+1]=(k+1)*nh*nw+(j+1)*nw;
        f[3*m+2]=(k+1)*nh*nw+j*nw;
        m++
        f[3*m+0]=k*nh*nw+j*nw+(nw-1);
        f[3*m+1]=(k+1)*nh*nw+(j+1)*nw+(nw-1);
        f[3*m+2]=k*nh*nw+(j+1)*nw+(nw-1);
        m++;
        f[3*m+0]=k*nh*nw+j*nw+(nw-1);
        f[3*m+1]=(k+1)*nh*nw+j*nw+(nw-1);
        f[3*m+2]=(k+1)*nh*nw+(j+1)*nw+(nw-1);
        m++
    }
    for(k=0;k<nd-1;k++)
    for(i=0;i<nw-1;i++) {
        f[3*m+0]=k*nh*nw+i;
        f[3*m+1]=(k+1)*nh*nw+i;
        f[3*m+2]=(k+1)*nh*nw+(i+1);
        m++;
        f[3*m+0]=k*nh*nw+i;
        f[3*m+1]=(k+1)*nh*nw+(i+1);
        f[3*m+2]=k*nh*nw+(i+1);
        m++
        f[3*m+0]=k*nh*nw+(nh-1)*nw+i;
        f[3*m+1]=(k+1)*nh*nw+(nh-1)*nw+(i+1);
        f[3*m+2]=(k+1)*nh*nw+(nh-1)*nw+i;
        m++;
        f[3*m+0]=k*nh*nw+(nh-1)*nw+i;
        f[3*m+1]=k*nh*nw+(nh-1)*nw+(i+1);
        f[3*m+2]=(k+1)*nh*nw+(nh-1)*nw+(i+1);
        m++
    }
    
    // create the topology of material tetrahedra
    n=0;
    for(i=0;i<nw-1;i++)
    for(j=0;j<nh-1;j++)
    for(k=0;k<nd-1;k++)
    for(l=0;l<5;l++)
        blockTetra(n++,i,j,k,tetraTopo[l],ge);

    // create rest vertices, copying the positions of material vertices
    // (the topology of rest tetrahedra is implicitly defined by its indices)
    configureRestGeometry(ge);

    // compute nodal volume
    computeNodalVolume(ge);

    return def.resolve(ge).promise();
}
/*
makeRing
*/
function makeRing(params) {
    var def=$.Deferred();
    var ge = new myGeometry();
    
    var Ri=params.Ri;
    var Ro=params.Ro;
    var th=params.th;
    var d=params.d;
    var ntheta=parseInt(2*Math.PI*Ro/d);    // volume elements in the outter circle
    var nxy=parseInt((Ro-Ri)/d)+1;          // number of vol. elem. rings in x-y plane
    var nz=parseInt(th/d)+1;                // number of vol. elem. rings in z
    var np=ntheta*nxy*nz;
    var nt=5*ntheta*(nxy-1)*(nz-1);
    var p=new Float32Array(np*3);
    var r=new Float32Array(nt*4*3);
    var t=new Uint16Array(nt*4);
    var Volume=new Float32Array(np);

    var i,j,k,l,n,m,R;
    var theta;

    ge.params=params;

    ge.Ri=Ri;
    ge.Ro=Ro;
    ge.th=th;
    ge.d=d;
    
    ge.ntheta=ntheta;
    ge.nxy=nxy;
    ge.nz=nz;
    ge.np=np;
    ge.nt=nt;
    ge.p=p;
    ge.r=r;
    ge.t=t;
    ge.Volume=Volume;
    
    console.log("Number of vertices:"+np+"("+ntheta+","+nxy+","+nz+")");
    console.log("Number of tetrahedra:",nt);

    // create material vertices
    m=0;
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy;j++)
    for(k=0;k<nz;k++)
    {
        n=k*nxy*ntheta+j*ntheta+i;
        theta=2*Math.PI*(i/ntheta);
        R=Ro*j+Ri*(nxy-j);
        R/=nxy;

        // material configuration
        p[3*n+0]=R*Math.cos(theta);
        p[3*n+1]=R*Math.sin(theta);
        p[3*n+2]=d*k-d*th/2;
        m++;
    }

    // create the topology of material tetrahedra
    n=0;
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy-1;j++)
    for(k=0;k<nz-1;k++)
    for(l=0;l<5;l++)
        ringTetra(n++,i,j,k,tetraTopo[l],ge);

    // create rest vertices, copying the positions of material vertices
    // (the topology of rest tetrahedra is implicitly defined by its indices)
    configureRestGeometry(ge);

    // compute nodal volume
    computeNodalVolume(ge);

    return def.resolve(ge).promise();
}
    
/*
makeSurface
Generate a 1 elastic tetrahedra thick shell from a surface mesh encoded in json format,
and a plexus of linear elastic elements linking the internal surface of the mesh to the
center of the model.

One prism:
     .c'
    /|\
   / | \
  /  |  \
a'.---+---.b'
 |   .c  |
 |  / \  |
 | /   \ |
 |/     \|
a.-------.b
 
Can be decomposed into 3 tetrahedra:
a'c',a'b',a'a
ca,cb,cc'
bb',bc',ba

where a,b,c are vertices of the original surface as loaded from the mesh file,
and a',b','c are the same vertices displaced along the normal of the surface.

@param {Object} param Object with model parametres
@param {Object} surf JSON Object with surface data
*/
function makeSurface(params) {
    var def=$.Deferred();
    $.getJSON(params.url,function(surf) {
        var ge = new myGeometry();
        var th=params.th;
        var d=params.d;

        var np;               // number of material vertices
        var nt;               // number of material tetrahedra
        var nf;               // number of surface faces (triangles)
        var t;                // array for material surface topology
        var p;                // array for material surface geometry
        var r;                // array for rest surface geometry
        var re;               // array for rest fibre geometry
        var Volume;           // array for tetrahedral volumes

        var i,j,k;
        var vol,n1,n2,n3,n4;
        var P,T,NP,NT;
        var NO,n,nor,a=[],b=[];

        // Configure vertices from loaded mesh
        P=surf.p;      // vertices in the mesh
        T=surf.t;      // triangles in the mesh
        NP=surf.np;    // number of vertices
        NT=surf.nt;    // number of triangles
    
        console.log("Number of vertices: "+surf.np);
        console.log("Number of surface triangles: "+surf.nt);

        // HACK: change model size
        for(i=0;i<NP*3;i++)
            P[i]*=4;
    
        // 1. Compute normalised normal vectors
        NO=new Float32Array(NP*3);
        for(i=0;i<NT;i++) {
            // per triangle normal obtained by cross product
            a[0]=P[3*T[3*i+1]+0]-P[3*T[3*i+0]+0];
            a[1]=P[3*T[3*i+1]+1]-P[3*T[3*i+0]+1];
            a[2]=P[3*T[3*i+1]+2]-P[3*T[3*i+0]+2];
            b[0]=P[3*T[3*i+2]+0]-P[3*T[3*i+0]+0];
            b[1]=P[3*T[3*i+2]+1]-P[3*T[3*i+0]+1];
            b[2]=P[3*T[3*i+2]+2]-P[3*T[3*i+0]+2];
            n=cross(a,b);
            // normal distributed among vertices
            NO[3*T[3*i+0]+0]+=n[0];
            NO[3*T[3*i+0]+1]+=n[1];
            NO[3*T[3*i+0]+2]+=n[2];
            NO[3*T[3*i+1]+0]+=n[0];
            NO[3*T[3*i+1]+1]+=n[1];
            NO[3*T[3*i+1]+2]+=n[2];
            NO[3*T[3*i+2]+0]+=n[0];
            NO[3*T[3*i+2]+1]+=n[1];
            NO[3*T[3*i+2]+2]+=n[2];
        }
        // normalisation of vertex normals
        for(i=0;i<NP;i++) {
            var val=NO[3*i+0]*NO[3*i+0]+NO[3*i+1]*NO[3*i+1]+NO[3*i+2]*NO[3*i+2];
            //console.log(NO[3*i+0],NO[3*i+1],NO[3*i+2]);
            nor=Math.sqrt(NO[3*i+0]*NO[3*i+0]+NO[3*i+1]*NO[3*i+1]+NO[3*i+2]*NO[3*i+2]);
            NO[3*i+0]/=nor;
            NO[3*i+1]/=nor;
            NO[3*i+2]/=nor;
        }

        // 2. Allocate memory
        np=NP*2;
        nt=3*NT;
        nf=NT;
        t=new Uint16Array(nt*4);             // array for cortex topology
        f=new Uint16Array(nf*3);             // array for surface topology
        p=new Float32Array(np*3);            // array for material surface geometry
        r=new Float32Array(nt*4*3);          // array for rest surface geometry
        re=new Float32Array(NP);             // array for rest fibre geometry
        Volume=new Float32Array(np);         // array for tetrahedral volumes

        ge.params=params;
    
        ge.np=np;
        ge.nt=nt;
        ge.nf=nf;
        ge.p=p;
        ge.t=t;
        ge.f=f;
        ge.r=r;
        ge.re=re;
        ge.Volume=Volume;

        // 3. Configure triangular surface topology
        // external vertices are indexed after internal ones
        for(i=0;i<NT*3;i++)
            f[i]=T[i]+NP;
            
        // 4. Configure cortex tetrahedra topology
        // each mesh triangle is extruded using 3 tetrahedra
        n=0;
        for(i=0;i<NT;i++)
        {
            t[4*n+0]=T[3*i+0]+NP;
            t[4*n+1]=T[3*i+2]+NP;
            t[4*n+2]=T[3*i+1]+NP;
            t[4*n+3]=T[3*i+0];
            n++;

            t[4*n+0]=T[3*i+2];
            t[4*n+1]=T[3*i+0];
            t[4*n+2]=T[3*i+1];
            t[4*n+3]=T[3*i+2]+NP;
            n++;

            t[4*n+0]=T[3*i+1];
            t[4*n+1]=T[3*i+1]+NP;
            t[4*n+2]=T[3*i+2]+NP;
            t[4*n+3]=T[3*i+0];
            n++;
        }

        // 5. Configure material cortex geometry.
        // internal vertices
        for(i=0;i<3*NP;i++)
            p[i]=P[i];
        // external vertices
        for(i=0;i<3*NP;i++)
            p[i+3*NP]=P[i]+params.th*NO[i];

        // 6. Configure rest cortex geometry
        configureRestGeometry(ge);

        // 7. Compute nodal volume
        computeNodalVolume(ge);
    
        // 8. Configure rest fibre length
        for(i=0;i<NP;i++)
            re[i]=Math.sqrt(Math.pow(p[3*i+0],2)+Math.pow(p[3*i+1],2)+Math.pow(p[3*i+2],2));
        
        def.resolve(ge);
    });
    
    return def.promise();
}
