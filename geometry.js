/**
	Hyperelastic growth model, Roberto Toro 2015
	Geometry functions
	
	The myGeometry object stores the material and rest configurations.
	There is also code to generate ring and brain geometries
	
	Depends on algebra.js 
*/

var ringTetraTopo=[
	"001 010 100 111",
	"000 010 100 001",
	"011 001 111 010",
	"101 111 001 100",
	"110 100 010 111"
];

function myGeometry() {
	// Variables storing geometry and topology
	this.np=0;		// number of material vertices in the model
	this.nt=0;		// number of elastic tetrahedra in the model
	this.p=0;		// material vertex coordinates
	this.r=0;		// rest tetrahedra geometry
	this.re=0;		// rest edge geometry
	this.t=0;		// tetrahedra vertex indices (topology) in the material configuration
	this.Volume=0;	// nodal volumes of the rest configuration
}

/**
tetraVol
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
/**
allTetraVol
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
configureRestGeometry
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
/**
computeNodalVolume
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

/**
Get global index of a tetrahedron within a ring based on its angular, radial and depth
indices.
@param {integer} i angular index of the tetrahedron
@param {integer} j radial index of the tetrahedron
@param {index} depth index of the tetrahedron
*/
function ringVind(i,j,k,ntheta,nxy) {
    return k*nxy*ntheta+j*ntheta+i;
}
/**
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
/**
makeRing
*/
function makeRing(params) {
	var ge = new myGeometry();
	
	var Ri=params.Ri;
	var Ro=params.Ro;
	var th=params.th;
	var d=params.d;
	var ntheta=parseInt(2*Math.PI*Ro/d);	// volume elements in the outter circle
	var nxy=parseInt((Ro-Ri)/d)+1;			// number of vol. elem. rings in x-y plane
	var nz=parseInt(th/d)+1;       			// number of vol. elem. rings in z
	var np=ntheta*nxy*nz;
	var nt=5*ntheta*(nxy-1)*(nz-1);
	var p=new Float32Array(np*3);
	var r=new Float32Array(nt*4*3);
	var t=new Uint16Array(nt*4);
	var Volume=new Float32Array(np);

	var i,j,k,l,n,m,R;
	var theta;

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
		p[3*n+2]=th*(k/nz);
		m++;
	}

	// create the topology of material tetrahedra
	n=0;
	for(i=0;i<ntheta;i++)
	for(j=0;j<nxy-1;j++)
	for(k=0;k<nz-1;k++)
	for(l=0;l<5;l++)
		ringTetra(n++,i,j,k,ringTetraTopo[l],ge);

	// create rest vertices, copying the positions of material vertices
	// (the topology of rest tetrahedra is implicitly defined by its indices)
	configureRestGeometry(ge);

	// compute nodal volume
	computeNodalVolume(ge);

	return ge;
}
	
/**
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
function makeSphere(params,surf)
{
	var ge = new myGeometry();
	var th=params.th;
	var d=params.d;

	var	np;			// number of material vertices
	var	nt;			// number of material tetrahedra
	var	t;			// array for material surface topology
	var p;			// array for material surface geometry
	var r;			// array for rest surface geometry
	var re;			// array for rest fibre geometry
	var Volume;		// array for tetrahedral volumes

	var i,j,k;
	var vol,n1,n2,n3,n4;
	var P,T,NP,NT;
	var	NO,n,nor,a=[],b=[];

	console.log("Number of vertices: "+np);
	console.log("Number of tetrahedra: "+nt);

	// Configure vertices from loaded mesh
	P=surf.p;	// vertices in the mesh
	T=surf.t;	// triangles in the mesh
	NP=surf.np;	// number of vertices
	NT=surf.nt;	// number of triangles
	
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
	t=new Uint16Array(nt*4);			// array for surface topology
	p=new Float32Array(np*3);			// array for material surface geometry
	r=new Float32Array(nt*4*3);			// array for rest surface geometry
	re=new Float32Array(NP);			// array for rest fibre geometry
	Volume=new Float32Array(np);		// array for tetrahedral volumes
	ge.params=params;
	ge.np=np;
	ge.nt=nt;
	ge.p=p;
	ge.t=t;
	ge.Volume=Volume;

	// 3. Configure tetrahedra topology
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

	// 4. Configure material surface geometry.
	// internal vertices
	for(i=0;i<3*NP;i++)
		p[i]=P[i];
	// external vertices
	for(i=0;i<3*NP;i++)
		p[i+3*NP]=P[i]+params.th*NO[i];

	// 5. Configure rest surface geometry
	configureRestGeometry(ge);

	// 6. Compute nodal volume
	computeNodalVolume(ge);
	
	// 7. Configure rest fibre length
	for(i=0;i<NP;i++)
		re[i]=Math.sqrt(Math.pow(p[3*i+0],2)+Math.pow(p[3*i+1],2)+Math.pow(p[3*i+2],2));
	
	return ge;
}
