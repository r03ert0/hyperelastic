/**
	Hyperelastic growth model, Roberto Toro 2015
	Growth functions
	
	Grow functions affect the rest configuration of the model's geometry
	
	Depends on algebra.js and geometry.js
*/

function myGrowth() {
	// Variables storing growth parametres
	this.g={};	// growth tensor in material coordinates (a 3x3 matrix)
	this.T=0;	// time constant
}

function initGrowth(params) {
	var gr=new myGrowth();
	gr.G=params.G;
	gr.T=params.T;
	return gr;
}
/**
growHomogeneous
*/
function growHomogeneous(ge,gr) {
	var nt=ge.nt;
	var t=ge.t;
	var p=ge.p;
	var r=ge.r;
    var H=gr.G;    // homogeneous growth factor
    
    var i,j,k,m,n;
    
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
growRingBorderInstantaneous
*/
function growRingBorderInstantaneous(ge,gr) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;
	var t=ge.t;
	var p=ge.p;
	var r=ge.r;
	var H=gr.G;

	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
    
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
growRingBorderProgressive
*/
function growRingBorderProgressive(ge,gr,time) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;
	var t=ge.t;
	var p=ge.p;
	var r=ge.r;
	var H=gr.G;
	var T=gr.T;
	
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
    
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
growRingTangentialInstantaneous. Each vertex in the original ring has coordinates x,y,z. The coordinates
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
function growRingTangentialInstantaneous(ge,gr) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;
	var t=ge.t;
	var r=ge.r;
	var Ro=ge.Ro;
	var Ri=ge.Ri;
	var D=gr.G;

    var i,j,k,l,m,n;
    var numbox,numtet;
    var im,ir;
    var theta,R,z;
    var a,di,dj,dtheta;
    
    for(i=0;i<ntheta;i++)
    for(j=0;j<nxy-1;j++)
    for(k=0;k<nz-1;k++)
    {
		numbox=i*(nxy-1)*(nz-1)+j*(nz-1)+k;
		for(l=0;l<5;l++) {
			numtet=numbox*5+l;
			a=ringTetraTopo[l].split(" ");
			for(m=0;m<4;m++) {
				ir=numtet*4+m;		// index of rest vertex
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
function growSurfaceInstantaneous(ge,gr) {
	var nt=ge.nt;
	var r=ge.r;
	var	H=gr.G;	// homogeneous growth
	
	var	i,j,k,ir;
	
	for(i=0;i<nt;i++) // each tetrahedron
	for(j=0;j<4;j++) {
		ir=4*i+j;		// rest vertex
		for(k=0;k<3;k++) {
			r[3*ir+k]*=H;
		}
	}
}