/**
Hyperelastic growth model, Roberto Toro 2015
Display

Depends on geometry.js, three.js and simulation.js
*/

function myDisplay() {
	this.renderer=null;		// three.js renderer
	this.scene=null;
	this.mesh=null;
	this.camera=null;
	this.trackball=null;
	this.material=null;
	this.geometry=null;
	this.color=null;
	this.lines=null;

	this.colormap=null;
	this.perspective=true;
	this.wireframe=true;	
}

/**
Initialises the three.js render
*/
function initDisplay(params) {
	var di = new myDisplay();
	var n,i,j,k,v1,v2,v3;
	var w=window.innerWidth;
	var h=window.innerHeight;
	
	di.colormap=params.colormap;
	di.perspective=params.perspective;
	di.wireframe=params.wireframe;

	di.renderer = new THREE.WebGLRenderer({antialias:true});
	di.renderer.setSize(w,h);
	di.renderer.setClearColor('white');
	document.body.appendChild(di.renderer.domElement);
	di.scene = new THREE.Scene();
	
	return di;
}

/**
initMesh
*/
function initMesh(di,si,params) {
	var np=si.ge.np;
	var nt=si.ge.nt;
	var p=si.ge.p;
	var t=si.ge.t;
	var w=window.innerWidth;
	var h=window.innerHeight;
	
	var n,i,j,k,v1,v2,v3;
	var indx=[0,2,1,0,3,2,0,1,3,1,2,3];
	
	di.colormap=params.colormap;
	di.perspective=params.perspective;
	di.wireframe=params.wireframe;

	if(di.camera)
		di.scene.remove(di.camera);
	
	if(di.perspective) {
		di.camera = new THREE.PerspectiveCamera( 75, w/h, 1, 100);
	}
	else {
		var Z=6;
		di.camera = new THREE.OrthographicCamera( -Z*w/h,Z*w/h,Z,-Z,1, 100);
	}
	di.camera.position.z = 10;
	di.camera.position.y=-20;
	di.camera.lookAt(di.scene.position);
	di.scene.add(di.camera);
	
	di.trackball = new THREE.TrackballControls(di.camera,di.renderer.domElement);

	di.geometry=new THREE.Geometry();
	for(n=0;n<np;n++)
		di.geometry.vertices.push(new THREE.Vector3(p[3*n],p[3*n+1],p[3*n+2]));
	for(n=0;n<nt;n++) {
		for(i=0;i<12;i+=3) {
			v1=t[n*4+indx[i+0]];
			v2=t[n*4+indx[i+1]];
			v3=t[n*4+indx[i+2]];
			di.geometry.faces.push(new THREE.Face3(v1,v2,v3));
		}
	}
	di.geometry.computeFaceNormals();
	di.geometry.computeVertexNormals();
	
	switch(di.colormap) {
		case "normal":
			di.material=new THREE.MeshNormalMaterial({wireframe:di.wireframe});
			break;
		case "deformation":
			di.material=new THREE.MeshBasicMaterial({wireframe:di.wireframe,vertexColors:THREE.FaceColors});
			break;
	}

	if(di.mesh!=null)
		di.scene.remove(di.mesh);
	di.mesh=new THREE.Mesh(di.geometry,di.material);
	di.scene.add(di.mesh);
}
/**
updateMesh
*/
function updateMesh(di,si) {
	var np=si.ge.np;
	var p=si.ge.p;
	var i;

	for(i=0;i<np;i++) {
		di.geometry.vertices[i].x=p[3*i];
		di.geometry.vertices[i].y=p[3*i+1];
		di.geometry.vertices[i].z=p[3*i+2];
	}
	di.geometry.verticesNeedUpdate=true;
}

/**
updateRingLines
*/
function updateRingLines(ge) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;

	var i,j,k,l,m,n;

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

/**
Animation function.
Request a new frame, call the render function, call the simulation functions and update
model display.
*/
function animate() {
	requestAnimationFrame(animate);
	display.renderer.render(display.scene,display.camera);
	display.trackball.update();

	simulationStep(simulation);
	
	if(simulationParams.colormap=="deformation" && simulationParams.geometry=="ring") {
		ring_deformationColor(simulation.ge);
	}
		
	updateMesh(display,simulation);

	display.geometry.computeFaceNormals();
	display.geometry.computeVertexNormals();
	display.geometry.normalsNeedUpdate = true;

	display.mesh.needsUpdate=true;
	display.geometry.colorsNeedUpdate=true;
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
	return {red:r,green:b,blue:g*0.95};
}
/**
Display deformation of a ring geometry as text
*/
function ring_deformationText(ge) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;
	var p=ge.p;
	var r=ge.r;
	var t=ge.t;
	var i,j,k,l,n;
	var numbox;
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
var minJ=0,maxJ=1;
function ring_deformationColor(ge) {
	var ntheta=ge.ntheta;
	var nxy=ge.nxy;
	var nz=ge.nz;
	var p=ge.p;
	var r=ge.r;
	var t=ge.t;
	var i,j,k,l,m,n;
	var numbox;
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
				display.geometry.faces[n*4+m].color.setRGB(color.red,color.green,color.blue);
				//display.geometry.faces[n*4+m].color.setRGB(1,0,0);
			}
		}
	}
	minJ=tminJ;
	maxJ=tmaxJ;
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
Map deformation of a surface geometry into colors.
*/
function surface_deformationColor() {
	var i,j,k,l,m,n;
	var numbox,numtet,im,ir;
	var J,vol0,vol1;
	var	tmp=0,mintmp,maxtmp,val;
	var vval,vval2,nval;

	vval=new Float32Array(geometry.vertices.length);
	vval2=new Float32Array(geometry.vertices.length);
	nval=new Uint8Array(geometry.vertices.length);
	for(i=0;i<geometry.vertices.length;i++)
		vval[i]=0;
	for(i=0;i<geometry.vertices.length;i++)
		nval[i]=0;

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

		if(i==0) {
			mintmp=J;
			maxtmp=J;
		}
		if(J>maxtmp)
			maxtmp=J;
		if(J<mintmp)
			mintmp=J;
		
		val=(J-minJ)/(maxJ-minJ);
		for(l=0;l<3;l++) {
			n=i+l;
			for(m=0;m<4;m++) {
				vval[geometry.faces[n*4+m].a]+=val;
				vval[geometry.faces[n*4+m].b]+=val;
				vval[geometry.faces[n*4+m].c]+=val;
				nval[geometry.faces[n*4+m].a]+=1;
				nval[geometry.faces[n*4+m].b]+=1;
				nval[geometry.faces[n*4+m].c]+=1;
			}
		}
	}
	maxJ=maxtmp;
	minJ=mintmp;

	// from face deformation to vertex colors
	for(i=0;i<geometry.vertices.length;i++) {
		vval[i]=vval[i]/nval[i];
	}
	for(j=0;j<5;j++) {
		for(i=0;i<geometry.vertices.length;i++) {
			vval2[i]=0;
			nval[i]=0;
		}
		for(i=0;i<geometry.faces.length;i++) {
			vval2[geometry.faces[i].a]+=vval[geometry.faces[i].b]+vval[geometry.faces[i].c];
			vval2[geometry.faces[i].b]+=vval[geometry.faces[i].c]+vval[geometry.faces[i].a];
			vval2[geometry.faces[i].c]+=vval[geometry.faces[i].a]+vval[geometry.faces[i].b];
			nval[geometry.faces[i].a]+=2;
			nval[geometry.faces[i].b]+=2;
			nval[geometry.faces[i].c]+=2;
		}
		for(i=0;i<geometry.vertices.length;i++) {
			vval[i]=vval2[i]/nval[i];
		}
	}
	mintmp=vval[0];
	maxtmp=vval[0];
	for(i=0;i<geometry.vertices.length;i++) {
		if(vval[i]>maxtmp)
			maxtmp=vval[i];
		if(vval[i]<mintmp)
			mintmp=vval[i];
	}
	for(i=0;i<geometry.vertices.length;i++)
		vval[i]=(vval[i]-mintmp)/(maxtmp-mintmp);    
	for(i=0;i<geometry.vertices.length;i++)
		color[i].setRGB(vval[i],vval[i],vval[i]);
	for(i=0;i<geometry.faces.length;i++) {
		geometry.faces[i].vertexColors[0]=color[geometry.faces[i].a];
		geometry.faces[i].vertexColors[1]=color[geometry.faces[i].b];
		geometry.faces[i].vertexColors[2]=color[geometry.faces[i].c];
	}
}
