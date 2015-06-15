/**
	Hyperelastic growth model
*/

// variables storing geometry and topology
var np;					// number of vertices in the model
var nt;					// number of elastic tetrahedra in the model
var p;					// material vertex coordinates
var r;					// rest tetrahedra geometry
var	re;					// rest edge geometry
var t;					// tetrahedra vertex indices (topology)
var Volume;				// nodal volumes of the rest configuration
var Velocity;			// velocity per vertex
var Force;				// force per vertex

// global variables
var gamma=0.1;			// damping
var time=0;				// history, actually
var dt;					// time step
var pi=4*Math.atan(1);	// perimeter of a circle divided by its diameter

// variables for ring geometry
var ntheta;
var nxy;
var nz;
var flag_running=false;

var tetraTopo=[
	"001 010 100 111",
	"000 010 100 001",
	"011 001 111 010",
	"101 111 001 100",
	"110 100 010 111"
];

// variables for display
var	maxJ=1,minJ=0;	// min and max deformation used to adjust colormap
var	 renderer,		// three.js renderer
	 scene,
	 mesh=null,
	 camera,
	 trackball,
	 material,
	 geometry,
	 lines;
	
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
@return Deferred promise
*/
function makeSurface(param)
{
	return $.get(param.url)
	.done(function(surf){
		var P,T,NP,NT,i,j;
		var	NO,n,nor,a=[],b=[];

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

		// 2. Configure tetrahedra vertices
		np=NP*2;
		console.log("Number of vertices:",np);
		nt=3*NT;
		console.log("Number of tetrahedra:",nt);
		t=new Uint16Array(nt*4);			// array for surface topology
		p=new Float32Array(np*3);			// array for material surface geometry
		r=new Float32Array(nt*4*3);			// array for rest surface geometry
		re=new Float32Array(NP);			// array for rest fibre geometry
		Volume=new Float32Array(np);		// array for tetrahedral volumes
		Velocity=new Float32Array(np*3);	// array for vertex velocities
		Force=new Float32Array(np*3);		// array for vertex forces
	
		// 3. configure tetrahedra topology
		n=0;
		// each mesh triangle is extruded using 3 tetrahedra
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
		//nt=nt/2;

		// 4. configure material surface geometry.
		// internal vertices
		for(i=0;i<3*NP;i++)
			p[i]=P[i];
		// external vertices
		for(i=0;i<3*NP;i++)
			p[i+3*NP]=P[i]+param.th*NO[i];

		// 5. configure rest surface geometry	
		for(i=0;i<nt;i++)
		for(j=0;j<4;j++)
		for(k=0;k<3;k++)
			r[(4*i+j)*3+k]=p[t[4*i+j]*3+k];

		// 6. compute nodal volume
		for(i=0;i<np;i++)
			Volume[i]=0;
		var vol;
		var n1,n2,n3,n4;
		for(i=0;i<nt;i++) {
			n1=t[4*i+0];
			n2=t[4*i+1];
			n3=t[4*i+2];
			n4=t[4*i+3];
			vol=tetraVol(
				p[3*n1+0], p[3*n1+1], p[3*n1+2],
				p[3*n2+0], p[3*n2+1], p[3*n2+2],
				p[3*n3+0], p[3*n3+1], p[3*n3+2],
				p[3*n4+0], p[3*n4+1], p[3*n4+2]
			);
			Volume[n1]+=vol/4;
			Volume[n2]+=vol/4;
			Volume[n3]+=vol/4;
			Volume[n4]+=vol/4;
		}
		
		// 7. configure rest fibre length
		for(i=0;i<NP;i++)
			re[i]=Math.sqrt(Math.pow(p[3*i+0],2)+Math.pow(p[3*i+1],2)+Math.pow(p[3*i+2],2));
		
		// 8. compute time step
		var a=0;	// average mesh spacing
		var n1,n2;
		for(i=0;i<nt;i++)
		for(j=0;j<4;j++) {
			n1=t[4*i+j];
			n2=t[4*i+(j+1)%4];
			a+=Math.sqrt(
				Math.pow(p[n1*3+0]-p[n2*3+0],2)+
				Math.pow(p[n1*3+1]-p[n2*3+1],2)+
				Math.pow(p[n1*3+2]-p[n2*3+2],2)
			);
		}
		a=a/nt/4;
		dt = 0.1*Math.sqrt(param.rho*a*a/param.K);
	})
  .fail(function(xhr,err) {
    console.log( "error",xhr,err);
  });
}
/**
makeRing
*/
function makeRing(param)
{
	var def = $.Deferred();
	var Ri=param.Ri,
		Ro=param.Ro,
		th=param.th,
		d=param.d;
    var i;
    var theta;
    
    ntheta=parseInt(2*pi*Ro/d); // volume elements in the outter circle
    nxy=parseInt((Ro-Ri)/d)+1;  // number of vol. elem. rings in x-y plane
    nz=parseInt(th/d)+1;        // number of vol. elem. rings in z
    
    np=ntheta*nxy*nz;
    console.log("Number of vertices:"+np+
                "("+ntheta+","+nxy+","+nz+")");

    p=new Float32Array(np*3);
    Volume=new Float32Array(np);
    Velocity=new Float32Array(3*np);
    Force=new Float32Array(np*3);
    
    nt=5*ntheta*(nxy-1)*(nz-1);
    console.log("Number of tetrahedra:",nt);
    
    t=new Uint16Array(nt*4);
    
    // vertices
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy;j++)
    for(k=0;k<nz;k++)
    {
        n=k*nxy*ntheta+j*ntheta+i;
        theta=2*pi*(i/ntheta);
        R=Ro*(j/nxy)+Ri*(1-j/nxy);
        
        // material configuration
        p[3*n+0]=R*Math.cos(theta);
        p[3*n+1]=R*Math.sin(theta);
        p[3*n+2]=th*(k/nz);
    }

    // tetrahedra
    var n=0;
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy-1;j++)
    for(k=0;k<nz-1;k++)
    for(l=0;l<5;l++)
        tetra(n++,i,j,k,tetraTopo[l]);
    nt=n;

	// compute time step	
	var a=0;	// average mesh spacing
	var n1,n2;
	for(i=0;i<nt;i++)
	for(j=0;j<4;j++) {
		n1=t[4*i+j];
		n2=t[4*i+(j+1)%4];
		a+=Math.sqrt(
			Math.pow(p[n1*3+0]-p[n2*3+0],2)+
			Math.pow(p[n1*3+1]-p[n2*3+1],2)+
			Math.pow(p[n1*3+2]-p[n2*3+2],2)
		);
	}
	a=a/nt/4;
	dt = 0.1*Math.sqrt(param.rho*a*a/param.K);

	// init rest configuration	
	r=new Float32Array(nt*4*3);
	for(i=0;i<nt;i++)
	for(j=0;j<4;j++)
	for(k=0;k<3;k++)
		r[(4*i+j)*3+k]=p[t[4*i+j]*3+k];

	// compute nodal volume
	for(i=0;i<np;i++)
		Volume[i]=0;
	var vol;
	var n1,n2,n3,n4;
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
	//nt=6*5*25;
	def.resolve();
	return def.promise();
}
/**
growHomogeneous
*/
function growHomogeneous(param) {
    var i,j,k;
    var H=param.H;    // homogeneous growth factor
    
	// for each tetrahedron
	for(i=0;i<nt;i++) {
		// for each tetra. node
		for(j=0;j<4;j++) {
			n=t[4*i+j];	// material node indices
			m=4*i+j;	// rest node indices
			// growth function: homogeneous growth
			for(k=0;k<3;k++) {
				r[m*3+k]=p[n*3+k]*H;
			}
		}
	}
}
/**
growBorderInstantaneous
*/
function growBorderInstantaneous(param) {
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var H=param.H;
    
    // tetrahedral boxes
	j=nxy-2;
    for(i=0;i<ntheta;i++)
    for(k=0;k<nz-1;k++)
    {
		// box element index
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		// box's tetrahedral element index
		for(l=0;l<5;l++) {
			numtet=numbox*5+l;
			// tetrahedron's node index
			for(m=0;m<4;m++) {
				im=t[numtet*4+m];	// material vertex index
				ir=numtet*4+m;		// rest vertex index
				for(n=0;n<3;n++)
					r[ir*3+n]=p[im*3+n]*H;
			}
		}
    }
}
/**
growBorderProgressive
*/
function growBorderProgressive(param) {
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var H=param.H;
	var T=param.T;
	
	if(time>T)
		return;
    
    // tetrahedral boxes
	j=nxy-2;
    for(i=0;i<ntheta;i++)
    for(k=0;k<nz-1;k++)
    {
		// box element index
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		// box's tetrahedral element index
		for(l=0;l<5;l++) {
			numtet=numbox*5+l;
			// tetrahedron's node index
			for(m=0;m<4;m++) {
				im=t[numtet*4+m];	// material vertex index
				ir=numtet*4+m;		// rest vertex index
				for(n=0;n<3;n++)
					r[ir*3+n]=p[im*3+n]*Math.pow(H,1/(T/dt));
			}
		}
    }
}
/**
growTangential. Each vertex in the original ring has coordinates x,y,z. The coordinates
x,y are in the plane of the ring, z is in the plane of its thickness. The function
growTangential produces a tangential expansion of the ring at rest, i.e., it alters the
x,y coordinates, without changing the z coordinate (thickness), nor the radial size of
each finite element.
WRONG FOLLOWS:
This is achieved by
displacing each x,y point radially. The amount of displacement is such that the angle
supported by each finite element will be multiplied by the dilatation parametre D. The
total perimeter of a ring at a distance R from the central axis is 2*pi*R and will become
2*pi*R*D. Then, each point x,y at distance R has to be displaced to a distance R*D, i.e.,
they have to be multiplied by a factor D such that (x*D)^2+(y*D)^2=(R*D)^2.
*/
function growTangential(param) {
    var i,j,k,l,m,n;
    var numbox,numtet;
    var im,ir;
    var theta,R,z;
	var Ri=param.Ri,
		Ro=param.Ro,
		th=param.th;
    var D=param.D;     // dilatation (mm), used for the growth function
    var a,di,dj,dtheta;
    
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy-1;j++)
    for(k=0;k<nz-1;k++)
    {
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		for(l=0;l<5;l++) {
			numtet=numbox*5+l;
			a=tetraTopo[l].split(" ");
			for(m=0;m<4;m++) {
				ir=numtet*4+m;		// rest vertex indexOf
				
				di=parseInt(a[m].charAt(0));
				dj=parseInt(a[m].charAt(1));
				theta=2*Math.PI*(i+di)/ntheta;
				R=Ro*((j+dj)/nxy)+Ri*(1-(j+dj)/nxy);
				dtheta=(D-1)*(Math.PI/ntheta)*(2*di-1);
				r[ir*3+0]=R*Math.cos(theta+dtheta);
				r[ir*3+1]=R*Math.sin(theta+dtheta);
			}
		}
	}
}
/**
growSurface
*/
function growSurface(param) {
	var	i,j,k,ir;
	var	H=param.H;	// homogeneous growth
	
	for(i=0;i<nt;i++) // each tetrahedron
	for(j=0;j<4;j++) {
		ir=4*i+j;		// rest vertex
		for(k=0;k<3;k++) {
			r[3*ir+k]*=H;
		}
	}
}
/**
Create a new tetrahedron topology.
@param {integer} n global index of the tetrahedron
@param {integer} i base angular index of the tetrahedron
@param {integer} j base radial index of the tetrahedron
@param {integer} k base depth index of the tetrahedron
@param {string} s description of tetrahedron topology relative to the base i, j, k indices
*/
function tetra(n,i,j,k,s) {
    var a=s.split(" ");
    var b=a.map(function(n) {
            return vind(
                (i+parseInt(n.charAt(0)))%ntheta,
                j+parseInt(n.charAt(1)),
                k+parseInt(n.charAt(2)));
        });
    t[4*n+0]=b[0];
    t[4*n+1]=b[1];
    t[4*n+2]=b[2];
    t[4*n+3]=b[3];
}
/**
Get global index of a tetrahedron within a ring based on its angular, radial and depth
indices.
@param {integer} i angular index of the tetrahedron
@param {integer} j radial index of the tetrahedron
@param {index} depth index of the tetrahedron
*/
function vind(i,j,k) {
    return k*nxy*ntheta+j*ntheta+i;
}
/**
Matrix inversion.
@param {Matrix} m A 3x3 matrix represented as a vector
*/
function invert(m) {
    var det;
    var w=new Object();
    
    det=m.b*m.f*m.g + m.c*m.d*m.h + m.a*m.e*m.i - m.c*m.e*m.g - m.a*m.f*m.h - m.b*m.d*m.i;
    
    w.a=(m.e*m.i - m.f*m.h)/det;
    w.b=(m.c*m.h - m.b*m.i)/det;
    w.c=(m.b*m.f - m.c*m.e)/det;
    
    w.d=(m.f*m.g - m.d*m.i)/det;
    w.e=(m.a*m.i - m.c*m.g)/det;
    w.f=(m.c*m.d - m.a*m.f)/det;
    
    w.g=(m.d*m.h - m.e*m.g)/det;
    w.h=(m.b*m.g - m.a*m.h)/det;
    w.i=(m.a*m.e - m.b*m.d)/det;
    
    return w;
}
/**
Matrix determinant.
@param {Matrix} a A 3x3 matrix represented as a vector
*/
function determinant(a)
{
    var det=
    a.b*a.f*a.g +
    a.c*a.d*a.h +
    a.a*a.e*a.i -
    a.c*a.e*a.g -
    a.a*a.f*a.h -
    a.b*a.d*a.i;
    
    return det;
}
/**
Vector addition.
@param {Vector} a 3x1 vector
@param {Vector} b 3x1 vector
*/
function add(a,b) {
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]];
}
/**
subtract
*/
function subtract(a,b) {
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]];
}
/**
addMat
*/
function addMat(a,b) {
    return {
        a:a.a+b.a, b:a.b+b.b, c:a.c+b.c,
        d:a.d+b.d, e:a.e+b.e, f:a.f+b.f,
        g:a.g+b.g, h:a.h+b.h, i:a.i+b.i};
}
/**
subMat
*/
function subMat(a,b) {
    return {
        a:a.a-b.a, b:a.b-b.b, c:a.c-b.c,
        d:a.d-b.d, e:a.e-b.e, f:a.f-b.f,
        g:a.g-b.g, h:a.h-b.h, i:a.i-b.i};
}
/**
cross
*/
function cross(a,b) {
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]];
}
/**
transpose
*/
function transpose(m) {
    return {a:m.a,b:m.d,c:m.g,
            d:m.b,e:m.e,f:m.h,
            g:m.c,h:m.f,i:m.i};
}
/**
trace
*/
function trace(m) {
    return m.a+m.e+m.i;
}
/**
mulMat
*/
function mulMat(a,b) {
    return {
    a:a.a*b.a+a.b*b.d+a.c*b.g,
    b:a.a*b.b+a.b*b.e+a.c*b.h,
    c:a.a*b.c+a.b*b.f+a.c*b.i,
    
    d:a.d*b.a+a.e*b.d+a.f*b.g,
    e:a.d*b.b+a.e*b.e+a.f*b.h,
    f:a.d*b.c+a.e*b.f+a.f*b.i,
    
    g:a.g*b.a+a.h*b.d+a.i*b.g,
    h:a.g*b.b+a.h*b.e+a.i*b.h,
    i:a.g*b.c+a.h*b.f+a.i*b.i };
}
/**
mulMatVec
*/
function mulMatVec(m,a) {
    return [
        m.a*a[0]+m.b*a[1]+m.c*a[2],
        m.d*a[0]+m.e*a[1]+m.f*a[2],
        m.g*a[0]+m.h*a[1]+m.i*a[2]];
}
/**
mulVecSca
*/
function mulVecSca(a,b) {
    return [a[0]*b,a[1]*b,a[2]*b];
}
/**
mulMatSca
*/
function mulMatSca(m,a) {
    return {
        a:m.a*a, b:m.b*a, c:m.c*a,
        d:m.d*a, e:m.e*a, f:m.f*a,
        g:m.g*a, h:m.h*a, i:m.i*a};
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
printMat
*/
function printMat(M,name) {
    console.log(name+":",
                M.a+" "+M.b+" "+M.c,
                M.d+" "+M.e+" "+M.f,
                M.g+" "+M.h+" "+M.i);
}
/**
printVec
*/
function printVec(V,name) {
    console.log(name+": "+V[0]+","+V[1]+","+V[2]);
}
/**
drawForce
*/
/*
function drawForce(n,f) {
        processing.beginShape(processing.LINES);
        processing.stroke(255,0,0);
        processing.vertex(p[3*n+0],p[3*n+1],p[3*n+2]);
        processing.vertex(p[3*n+0]+f[0],
                          p[3*n+1]+f[1],
                          p[3*n+2]+f[2]);
        processing.endShape();
}
*/
/**
growth
*/
function tetraElasticity(param) {
    var mu=param.mu;    // shear modulus
    var K=param.K;   // bulk modulus
    var n1,n2,n3,n4;    // material tetra vertex indices
	var m1,m2,m3,m4;		// rest tetra vertex indices
    var a,b,c;
    var i,j;
    var x1,x2,x3,x4;
    var ii;
    var Ar=new Object();       // rest tetra
    var A=new Object();        // material tetra
    var F;           // deformation tensor
    var B,J;         // F*F^T, det(F)
    var J1,J2,J3,J4;
    var vol;
    var N1,N2,N3,N4; // tetra normal vectors
    var S,Ss,Sv;     // Stress
    var Ue,Us,Uv;    // Elastic energy
    var I=new Object({ a:1,b:0,c:0,
            d:0,e:1,f:0,
            g:0,h:0,i:1});
        
    // integrate elastic forces
    for(i=0;i<3*np;i++)
        Force[i]=0;
    Ue=0;
    for(i=0;i<nt;i++) {
        n1=t[4*i+0];
        n2=t[4*i+1];
        n3=t[4*i+2];
        n4=t[4*i+3];

		m1=4*i+0;
		m2=4*i+1;
		m3=4*i+2;
		m4=4*i+3;
        
        // material tetra
        a=p[3*n1+0];
        b=p[3*n1+1];
        c=p[3*n1+2];
        A.a=p[3*n2+0]-a;
        A.d=p[3*n2+1]-b;
        A.g=p[3*n2+2]-c;
        A.b=p[3*n3+0]-a;
        A.e=p[3*n3+1]-b;
        A.h=p[3*n3+2]-c;
        A.c=p[3*n4+0]-a;
        A.f=p[3*n4+1]-b;
        A.i=p[3*n4+2]-c;
        
        x1=[A.a,A.d,A.g];
        x2=[A.b,A.e,A.h];
        x3=[A.c,A.f,A.i];
        
        // tetra face negative normals
        // (because traction s=-S*n)
        N1 = cross(x3,x1);
        N2 = cross(x2,x3);
        N3 = cross(x1,x2);
        N4 = cross(subtract(x2,x3),subtract(x1,x3));

         // rest tetra
        a=r[3*m1+0];
        b=r[3*m1+1];
        c=r[3*m1+2];
        Ar.a=r[3*m2+0]-a;
        Ar.d=r[3*m2+1]-b;
        Ar.g=r[3*m2+2]-c;
        Ar.b=r[3*m3+0]-a;
        Ar.e=r[3*m3+1]-b;
        Ar.h=r[3*m3+2]-c;
        Ar.c=r[3*m4+0]-a;
        Ar.f=r[3*m4+1]-b;
        Ar.i=r[3*m4+2]-c;

        // deformation tensor
        //printMat(A,"A");
        //printMat(Ar,"Ar");
        F=mulMat(A,invert(Ar));
        //printMat(F,"F");
        
        J=determinant(F);
        //console.log("J: "+J);
        
        /*
        	Shear stress: mu (FF'-I tr(FF')/3) / J^(-5/3)
        */
        B=mulMat(F,transpose(F));
        //printMat(B,"B");
        
        Ss=mulMatSca(subMat(B,mulMatSca(I,trace(B)/3)),
                     mu/Math.pow(J,5/3));
        //printMat(Ss,"Ss");

        /*
        	Bulk stress: KI(J-1)
        */
        // (why not K(1-1/J)I ??)
        Sv = mulMatSca(I,K*(J-1));
        //printMat(Sv,"Sv");
        
        // total
        S = addMat(Ss,Sv);
        //printMat(S,"S");

        // distribute forces among tetra verts
        // (recycling x1, x2, x3, adding x4)
        x1=mulVecSca(mulMatVec(S,add(add(N1,N2),N3)),1/6);
        x2=mulVecSca(mulMatVec(S,add(add(N1,N3),N4)),1/6);
        x3=mulVecSca(mulMatVec(S,add(add(N2,N3),N4)),1/6);
        x4=mulVecSca(mulMatVec(S,add(add(N1,N2),N4)),1/6);
        //drawForce(n1,x1);
        //drawForce(n2,x2);
        //drawForce(n3,x3);
        //drawForce(n4,x4);
        Force[3*n1+0]+=x1[0];
        Force[3*n1+1]+=x1[1];
        Force[3*n1+2]+=x1[2];
        Force[3*n2+0]+=x2[0];
        Force[3*n2+1]+=x2[1];
        Force[3*n2+2]+=x2[2];
        Force[3*n3+0]+=x3[0];
        Force[3*n3+1]+=x3[1];
        Force[3*n3+2]+=x3[2];
        Force[3*n4+0]+=x4[0];
        Force[3*n4+1]+=x4[1];
        Force[3*n4+2]+=x4[2];
		
		Us = 0.5*mu*(trace(B)/Math.pow(J,2/3)-3);
		Uv = K*(J-Math.log(J)-1);
		Ue += Us + Uv;
    }
    console.log(time,Ue);
}
/**
linElasticity
*/
function linElasticity(param) {
	var	i,j,k;
	var	a,b,c;
	var	l,l0;
	var avf,nor;
	var	f=[];
	var K=param.Kf;
	
	// grey/white interface tetrahedra are those multiple of 3 +1
	n=0;
	avf=0;
	for(i=0;i<nt;i+=3)
	for(j=0;j<3;j++) {
		a=t[4*(i+1)+j];
		
		
		l=Math.sqrt(Math.pow(p[3*a+0],2)+Math.pow(p[3*a+1],2)+Math.pow(p[3*a+2],2));
		l0=re[a];
		
		nor=Math.sqrt(Math.pow(p[3*a+0],2)+Math.pow(p[3*a+1],2)+Math.pow(p[3*a+2],2));
		f[0]=K*(p[3*a+0]/nor)*(1-l/l0);
		f[1]=K*(p[3*a+1]/nor)*(1-l/l0);
		f[2]=K*(p[3*a+2]/nor)*(1-l/l0);
		Force[3*a+0]+=f[0];
		Force[3*a+1]+=f[1];
		Force[3*a+2]+=f[2];
		
		n++;
	}
}
/**
Integrate velocity into displacement
*/
function move(param) {
    var i,j,k,l,m,n;
	var rho=param.rho;

	// damped motion equation
	for(i=0;i<3*np;i++) {
		j=parseInt(i/3);
		Force[i]-=Velocity[i]*gamma*Volume[j];
		Velocity[i]+=Force[i]/(Volume[j]*rho)*dt;
		p[i]+=Velocity[i]*dt;
	}
    
	// update mesh
	for(i=0;i<np;i++) {
		geometry.vertices[i].x=p[3*i];
		geometry.vertices[i].y=p[3*i+1];
		geometry.vertices[i].z=p[3*i+2];
	}
	geometry.verticesNeedUpdate=true;
	
	if(param.fibres==true || lines==undefined)
		return;
	// update lines
	// Lines
	m=0;
    for(j=0;j<nxy;j+=2)
    for(k=0;k<nz;k+=2) {
    	l=0;
		for(i=0;i<=ntheta;i++)
		{
			n=k*nxy*ntheta+j*ntheta+(i%ntheta);
			lines[m].geometry.vertices[l]=geometry.vertices[n];
			l++;
		}
		m++;
	}
	for(i=0;i<=ntheta;i+=8) {
		l=0;
		k=0;
		for(j=0;j<nxy;j++) {
			n=k*nxy*ntheta+j*ntheta+(i%ntheta);
			lines[m].geometry.vertices[l]=geometry.vertices[n];
			l++;
		}
		j=nxy-1;
		for(k=0;k<nz;k++) {
			n=k*nxy*ntheta+j*ntheta+(i%ntheta);
			lines[m].geometry.vertices[l]=geometry.vertices[n];
			l++;
		}
		k=nz-1;
		for(j=nxy-1;j>=0;j--) {
			n=k*nxy*ntheta+j*ntheta+(i%ntheta);
			lines[m].geometry.vertices[l]=geometry.vertices[n];
			l++;
		}
		j=0;
		for(k=nz-1;k>=0;k--) {
			n=k*nxy*ntheta+j*ntheta+(i%ntheta);
			lines[m].geometry.vertices[l]=geometry.vertices[n];
			l++;
		}
		m++;
	}
	for(i=0;i<lines.length;i++)
		lines[i].geometry.verticesNeedUpdate=true;
}
/*=================
       Display
  =================*/
/**
Initialises the three.js render
*/
function initRender() {
	var n,i,j,k,v1,v2,v3;
	var indx=[0,2,1,0,3,2,0,1,3,1,2,3];
	var w=window.innerWidth;
	var h=window.innerHeight;

	renderer = new THREE.WebGLRenderer({antialias:true});
	renderer.setSize(w,h);
	renderer.setClearColor('white');
	document.body.appendChild(renderer.domElement);
	scene = new THREE.Scene();
	//camera = new THREE.PerspectiveCamera( 75, w/h, 1, 100);
	var Z=6;
	camera = new THREE.OrthographicCamera( -Z*w/h,Z*w/h,Z,-Z,1, 100);
	camera.position.z = 10;
	camera.position.y=-20;
	camera.lookAt(scene.position);
	scene.add(camera);
	trackball = new THREE.TrackballControls(camera,renderer.domElement);
}
/**
Render and update virtual trackball
*/
function render() {
	renderer.render(scene,camera);
	trackball.update();
}
/**
Animation function.
Request a new frame, call the render function, call the simulation functions and update
model display.
*/
function animate() {
	requestAnimationFrame(animate);
	render();

	if(time>0.1 || flag_running==false)
		return;
		
	time+=dt;
	$("#log").html("H:"+param.H.toFixed(4)+"<br />t:"+time.toFixed(4)+"<br />maxJ:"+maxJ.toFixed(4)+"<br />minJ:"+minJ.toFixed(4));

	if(param.growthFunc) {
		tetraElasticity(param);
		if(param.fibres)
			linElasticity(param);
		param.deformFunc();
		move(param);
	}

	//growBorderProgressive(param);
	geometry.computeFaceNormals();
	geometry.computeVertexNormals();
	geometry.normalsNeedUpdate = true;

	mesh.needsUpdate=true;
	geometry.colorsNeedUpdate=true;
}
/**
Color map: from dark blue to dark red for negative to positive values
@param {float} x Value between -1 and 1. Values outside the boundaries are clamped.
*/
function bdarkb2r(x) {
	var	r,g,b;
	
	if(x<-1) 		{	r=0;		g=0;		b=0.5;	}
	else if(x<-0.5) {	r=0;		g=0;		b=1.5+x;}
	else if(x<0) 	{	r=2*(x+0.5);g=2*(x+0.5);b=1;	}
	else if(x<0.5)	{	r=1;		g=1-2*x;	b=1-2*x;}
	else if(x<1)	{	r=1.5-x;	g=0;		b=0;	}
	else 			{	r=0.5;		g=0;		b=0;	}
/*
red_top     = [0.5 0 0];
red_middle  = [1 0 0];
white_middle= [1 1 1];
blue_middle = [0 0 1];
blue_bottom = [0 0 0.5];
*/
	return {red:r,green:b,blue:g};
}
/**
Color map: from black to white for negative to positive values.
@param {float} x Value between -1 and 1. Values outside the boundaries are clamped.
*/
function black2white(x) {
	var	r,g,b,w;
	
	r=(x+1)/2;
	g=(x+1)/2;
	b=(x+1)/2;
	w=1-Math.abs(x)*0.0;
	
	if(x<0) { g*=w;b*=w;}
	if(x>0) { r*=w;g*=w;}
	return {red:r,green:b,blue:g};
}
/**
Display deformation of a ring geometry as text
*/
function ring_deformationText() {
/*
	Map deformation into color for ring geometries
*/
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var J,vol0,vol1;
	var	tminJ,tmaxJ;
    
    // tetrahedral boxes
	for(j=0;j<nxy-1;j++)
    for(i=0;i<ntheta;i++)
    for(k=0;k<nz-1;k++)
    {
		// box element index
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		// box's tetrahedral element index
		vol0=0;
		vol1=0;
		for(l=0;l<5;l++) {
			n=numbox*5+l;

			vol0+=tetraVol(
						r[3*(4*n+0)+0],r[3*(4*n+0)+1],r[3*(4*n+0)+2],
						r[3*(4*n+1)+0],r[3*(4*n+1)+1],r[3*(4*n+1)+2],
						r[3*(4*n+2)+0],r[3*(4*n+2)+1],r[3*(4*n+2)+2],
						r[3*(4*n+3)+0],r[3*(4*n+3)+1],r[3*(4*n+3)+2]);
			vol1+=tetraVol(
						p[3*t[4*n+0]+0],p[3*t[4*n+0]+1],p[3*t[4*n+0]+2],
						p[3*t[4*n+1]+0],p[3*t[4*n+1]+1],p[3*t[4*n+1]+2],
						p[3*t[4*n+2]+0],p[3*t[4*n+2]+1],p[3*t[4*n+2]+2],
						p[3*t[4*n+3]+0],p[3*t[4*n+3]+1],p[3*t[4*n+3]+2]);
		}
		J=Math.log(vol1/vol0);
		if(i+j+k==0) {
			tminJ=J;
			tmaxJ=J;
		}
		if(J>tmaxJ)
			tmaxJ=J;
		if(J<tminJ)
			tminJ=J;
		console.log(i,j,k,J);
    }
}
/**
Map deformation of a ring geometry into colors
*/
function ring_deformationColor() {
/*
	Map deformation into color for ring geometries
*/
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var J,vol0,vol1;
	var	tminJ,tmaxJ,color;
    
    // tetrahedral boxes
	for(j=0;j<nxy-1;j++)
    for(i=0;i<ntheta;i++)
    for(k=0;k<nz-1;k++)
    {
		// box element index
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		// box's tetrahedral element index
		vol0=0;
		vol1=0;
		for(l=0;l<5;l++) {
			n=numbox*5+l;

			vol0+=tetraVol(
						r[3*(4*n+0)+0],r[3*(4*n+0)+1],r[3*(4*n+0)+2],
						r[3*(4*n+1)+0],r[3*(4*n+1)+1],r[3*(4*n+1)+2],
						r[3*(4*n+2)+0],r[3*(4*n+2)+1],r[3*(4*n+2)+2],
						r[3*(4*n+3)+0],r[3*(4*n+3)+1],r[3*(4*n+3)+2]);
			vol1+=tetraVol(
						p[3*t[4*n+0]+0],p[3*t[4*n+0]+1],p[3*t[4*n+0]+2],
						p[3*t[4*n+1]+0],p[3*t[4*n+1]+1],p[3*t[4*n+1]+2],
						p[3*t[4*n+2]+0],p[3*t[4*n+2]+1],p[3*t[4*n+2]+2],
						p[3*t[4*n+3]+0],p[3*t[4*n+3]+1],p[3*t[4*n+3]+2]);
		}
		J=Math.log(vol1/vol0);
		if(i+j+k==0) {
			tminJ=J;
			tmaxJ=J;
		}
		if(J>tmaxJ)
			tmaxJ=J;
		if(J<tminJ)
			tminJ=J;
		for(l=0;l<5;l++) {
			n=numbox*5+l;
			for(m=0;m<4;m++) {
				color=black2white(2*(J-minJ)/(maxJ-minJ)-1);
				geometry.faces[n*4+m].color.setRGB(color.red,color.green,color.blue);
			}
		}
    }
    minJ=tminJ;
    maxJ=tmaxJ;
}
/**
Map deformation of a surface geometry into colors.
*/
function surface_deformationColor() {
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var J,vol0,vol1;
	var	tmp=0,mintmp=1,color;
    
    // tetrahedral prisms
    for(i=0;i<nt;i+=3)
    {
		vol0=0;
		vol1=0;
		for(l=0;l<3;l++) {
			n=i+l;

			vol0+=tetraVol(
						r[3*(4*n+0)+0],r[3*(4*n+0)+1],r[3*(4*n+0)+2],
						r[3*(4*n+1)+0],r[3*(4*n+1)+1],r[3*(4*n+1)+2],
						r[3*(4*n+2)+0],r[3*(4*n+2)+1],r[3*(4*n+2)+2],
						r[3*(4*n+3)+0],r[3*(4*n+3)+1],r[3*(4*n+3)+2]);
			vol1+=tetraVol(
						p[3*t[4*n+0]+0],p[3*t[4*n+0]+1],p[3*t[4*n+0]+2],
						p[3*t[4*n+1]+0],p[3*t[4*n+1]+1],p[3*t[4*n+1]+2],
						p[3*t[4*n+2]+0],p[3*t[4*n+2]+1],p[3*t[4*n+2]+2],
						p[3*t[4*n+3]+0],p[3*t[4*n+3]+1],p[3*t[4*n+3]+2]);
		}
		J=Math.log(vol1/vol0);
		if(Math.abs(J)>tmp)
			tmp=Math.abs(J);
		if(Math.abs(J)<mintmp)
			mintmp=Math.abs(J);
		color=black2white(Math.sign(J)*(Math.abs(J)-minJ)/(maxJ-minJ));
		for(l=0;l<3;l++) {
			n=i+l;
			for(m=0;m<4;m++) {
				geometry.faces[n*4+m].color.setRGB(color.red,color.green,color.blue);
			}
		}
    }
    maxJ=tmp+0.0001;
    minJ=mintmp;
}
function displayParametres(param) {
	console.log("display parametres",param);
	var i;
	$("#params").html("");
	for(key in param) {
		$("#params").append(key+" <input type='text' id='"+key+"' size='10'/><br/>");		
		$("#"+key).val(param[key].toString());
	}
}
function startStop() {
	if(flag_running==true) {
		flag_running=false;
		$("#startStop").html("Start");
	}
	else {
		flag_running=true;
		animate(param);
		param.growthFunc(param);
		$("#startStop").html("Stop");
	}
}
function selectSimulation() {
	console.log($("#select").val());
	switch($("#select").val()) {
		case "SphereSurface":
			param=SphereSurface;
			break;
		case "EllipsoidSurface":
			param=EllipsoidSurface;
			break;
		case "RingTangential":
			param=RingTangential;
			break;
		case "RingBorder":
			param=RingBorder;
			break;
	}
	console.log(param);
	initSimulation(param).then(render);
}
function initSimulation(param) {

	return param.geomFunc(param).then(function() {
		var n,i,j,k,v1,v2,v3;
		var indx=[0,2,1,0,3,2,0,1,3,1,2,3];

		displayParametres(param);

		// Mesh
		geometry=new THREE.Geometry();
		for(n=0;n<np;n++)
			geometry.vertices.push(new THREE.Vector3(p[3*n],p[3*n+1],p[3*n+2]));
		for(n=0;n<nt;n++) {
			for(i=0;i<12;i+=3) {
				v1=t[n*4+indx[i]];
				v2=t[n*4+indx[i+1]];
				v3=t[n*4+indx[i+2]];
				geometry.faces.push(new THREE.Face3(v1,v2,v3));
			}
		}
		geometry.computeFaceNormals();
		geometry.computeVertexNormals();
		//material=new THREE.MeshNormalMaterial({color:'blue',wireframe:false});
		material=new THREE.MeshBasicMaterial({wireframe:false,linewidth:0.001,shading:THREE.FlatShading,vertexColors:THREE.FaceColors});
		//material = new THREE.ShaderMaterial({vertexShader: "varying vec3 vnormal;void main(){vnormal=normal;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);}", fragmentShader: "varying vec3 vnormal;void main(){vec3 n=normalize(vec3(1,1,1)+vnormal);gl_FragColor=vec4(n,1);}",shading:THREE.SmoothShading});
		
		if(mesh!=null)
			scene.remove(mesh);
		mesh=new THREE.Mesh(geometry,material);
		//mesh.rotation.x=-3.1415927/3;
		//mesh.dynamic = true;
		scene.add(mesh);
	
	
		if(param.fibres==false) {
		// Lines
			lines=[];
			for(j=0;j<nxy;j+=2)
			for(k=0;k<nz;k+=2) {
				var linegeom=new THREE.Geometry();
				for(i=0;i<=ntheta;i++)
				{
					n=k*nxy*ntheta+j*ntheta+(i%ntheta);
					linegeom.vertices.push(geometry.vertices[n]);
				}
				var linemat=new THREE.LineBasicMaterial({color:0,linewidth:1});
				var line=new THREE.Line(linegeom,linemat);
				scene.add(line);
				line.translateZ(0.01);
				line.translateY(-0.01);
				lines.push(line);
			}
			for(i=0;i<=ntheta;i+=8) {
				var linegeom=new THREE.Geometry();
				k=0;
				for(j=0;j<nxy;j++) {
					n=k*nxy*ntheta+j*ntheta+(i%ntheta);
					linegeom.vertices.push(geometry.vertices[n]);
				}
				j=nxy-1;
				for(k=0;k<nz;k++) {
					n=k*nxy*ntheta+j*ntheta+(i%ntheta);
					linegeom.vertices.push(geometry.vertices[n]);
				}
				k=nz-1;
				for(j=nxy-1;j>=0;j--) {
					n=k*nxy*ntheta+j*ntheta+(i%ntheta);
					linegeom.vertices.push(geometry.vertices[n]);
				}
				j=0;
				for(k=nz-1;k>=0;k--) {
					n=k*nxy*ntheta+j*ntheta+(i%ntheta);
					linegeom.vertices.push(geometry.vertices[n]);
				}
				var linemat=new THREE.LineBasicMaterial({color:0,linewidth:1});
				var line=new THREE.Line(linegeom,linemat);
				scene.add(line);
				lines.push(line);
				line.translateZ(0.01);
				line.translateY(-0.01);
			}
		}
	});
}

var param;

/* Model 1: Ring, Border growth */
var RingBorder=new Object({
	rho:0.0001,					 // mass density
	K:100,						 // bulk modulus
	mu:100,						 // shear modulus
	Ri:2,   					 // inner radius in mm
	Ro:6,   					 // outter radius in mm
	th:0.5, 				     // ring thickness in mm
	d:0.25,  					 // typical length of a volume elements
	D:1,    					 // dilatation (mm), used for the growTangential function
	H:1.8,						 // growth
	T:0.05,						 // duration of growth (in sec)
	fibres:false,
	deformFunc:ring_deformationColor,
	geomFunc:makeRing,
	growthFunc:growBorderInstantaneous // growth function
});

/* Model 2: Ring, Tangential growth */
var RingTangential=new Object({
	rho:0.0001,					 // mass density
	K:50,						 // bulk modulus
	mu:1,						 // shear modulus
	Ri:2,   					 // inner radius in mm
	Ro:6,   					 // outter radius in mm
	th:0.5,	 				     // ring thickness in mm
	d:0.25,  					 // typical length of a volume elements
	D:1.5,    					 // dilatation (mm), used for the growTangential function
	H:2,						 // growth
	T:0.05,						 // duration of growth (in sec)
	fibres:false,
	deformFunc:ring_deformationColor,
	geomFunc:makeRing,
	growthFunc:growTangential	 // growth function
});

/* Model 3: Sphere, Surface growth */
var SphereSurface=new Object({
	rho:0.0001,					 // mass density
	K:50,						 // bulk modulus of tetrahedra
	mu:50,						 // shear modulus of tetrahedra
	Kf:1,						 // elastic constant for fibres
	url:"data/sphere-1115.json",	 // surface mesh URL
	th:0.25,	 				     // surface thickness in mm
	H:1.6,						 // growth
	fibres:true,
	deformFunc:surface_deformationColor,
	geomFunc:makeSurface,
	growthFunc:growSurface		// growth function
});

/* Model 4: Ellipsoid, Surface growth */
var EllipsoidSurface=new Object({
	rho:0.0001,					 // mass density
	K:50,						 // bulk modulus of tetrahedra
	mu:50,						 // shear modulus of tetrahedra
	Kf:0,						 // elastic constant for fibres
	url:"data/ellipsoid-2582.json",	 // surface mesh URL
	th:0.25,	 				     // surface thickness in mm
	H:1.1,						 // growth
	fibres:true,
	deformFunc:surface_deformationColor,
	geomFunc:makeSurface,
	growthFunc:growSurface		// growth function
});

param=SphereSurface;
initRender();
animate();
initSimulation(param);

$("#select").change(selectSimulation);
$("#startStop").click(startStop);
