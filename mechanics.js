/**
	Hyperelastic growth model, Roberto Toro 2015
	Mechanics
	
	The myMechanics object stores mechanical properties of the model, as well
	as vectors for the force and velocity of each material vertex.
*/

function myMechanics() {
	this.gamma=0.1;		// damping
	this.rho=0;			// mass density
	this.K=0;			// bulk modulus
	this.mu=0;			// shear modulus
	this.Kf=0;			// linear spring Young's modulus

	this.Velocity=0;	// velocity per material vertex
	this.Force=0;		// force per material vertex
}

/**
initMechanics
*/
function initMechanics(ge,params) {
	var me = new myMechanics();
	var np=ge.np;
	var Velocity=new Float32Array(np*3);	// array for vertex velocities
	var Force=new Float32Array(np*3);		// array for vertex forces

	me.Velocity=Velocity;
	me.Force=Force;
	
	me.gamma=params.gamma;
	me.rho=params.rho;
	me.K=params.K;
	me.mu=params.mu;
	me.Kf=params.Kf;
	
	return me;
}
/**
tetraElasticity
*/
function tetraElasticity(ge,me) {
    var np=ge.np;
    var nt=ge.nt;
    var p=ge.p;
    var t=ge.t;
    var r=ge.r;
    var Force=me.Force;
    var mu=me.mu;   	// shear modulus
    var K=me.K;   		// bulk modulus
    
    var n1,n2,n3,n4;    // material tetra vertex indices
	var m1,m2,m3,m4;	// rest tetra vertex indices
    var a,b,c;
    var i,j;
    var x1,x2,x3,x4;
    var ii;
    var Ar=new Object();// rest tetra
    var A=new Object(); // material tetra
    var F;           	// deformation tensor
    var B,J;         	// F*F^T, det(F)
    var J1,J2,J3,J4;
    var vol;
    var N1,N2,N3,N4; 	// tetra normal vectors
    var S,Ss,Sv;     	// Stress
    var Ue,Us,Uv;    	// Elastic energy
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
        
        // tetra face negative normals (because traction s=-S*n)
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
        F=mulMat(A,invert(Ar));
        J=determinant(F);
        
        // Shear stress: mu (FF'-I tr(FF')/3) / J^(-5/3)
        B=mulMat(F,transpose(F));
        Ss=mulMatSca(subMat(B,mulMatSca(I,trace(B)/3)),mu/Math.pow(J,5/3));

        // Bulk stress: KI(J-1)   (why not K(1-1/J)I ??)
        Sv = mulMatSca(I,K*(J-1));
        
        // Total
        S = addMat(Ss,Sv);

        // distribute forces among tetra verts (recycling x1, x2, x3, adding x4)
        x1=mulVecSca(mulMatVec(S,add(add(N1,N2),N3)),1/6);
        x2=mulVecSca(mulMatVec(S,add(add(N1,N3),N4)),1/6);
        x3=mulVecSca(mulMatVec(S,add(add(N2,N3),N4)),1/6);
        x4=mulVecSca(mulMatVec(S,add(add(N1,N2),N4)),1/6);
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
    // store elastic energy values
    me.Us=Us;
    me.Uv=Uv;
    me.Ue=Ue;
}
/**
linElasticity.
Issue: This function is now specific to the surface models. It should be written
       in generic terms.
*/
function linElasticity(ge,me) {
    var nt=ge.nt;
    var p=ge.p;
    var t=ge.t;
    var re=ge.re;
    var Force=me.Force;
	var Kf=me.Kf;

	var	i,j,k;
	var	a,l,l0;
	var nor;
	var	f=[];
	
	// grey/white interface tetrahedra are those multiple of 3 +1
	for(i=0;i<nt;i+=3)
	for(j=0;j<3;j++) {
		a=t[4*(i+1)+j];
		l=Math.sqrt(Math.pow(p[3*a+0],2)+Math.pow(p[3*a+1],2)+Math.pow(p[3*a+2],2));
		l0=re[a];
		nor=Math.sqrt(Math.pow(p[3*a+0],2)+Math.pow(p[3*a+1],2)+Math.pow(p[3*a+2],2));
		f[0]=Kf*(p[3*a+0]/nor)*(1-l/l0);
		f[1]=Kf*(p[3*a+1]/nor)*(1-l/l0);
		f[2]=Kf*(p[3*a+2]/nor)*(1-l/l0);
		Force[3*a+0]+=f[0];
		Force[3*a+1]+=f[1];
		Force[3*a+2]+=f[2];
	}
}
/**
Integrate velocity into displacement
*/
function move(ge,me,dt) {
	var np=ge.np;
	var p=ge.p;
	var Volume=ge.Volume;
	var rho=me.rho;
	var gamma=me.gamma;
	var Force=me.Force;
	var Velocity=me.Velocity;
    var i,j,k,l,m,n;

	// damped motion equation
	for(i=0;i<3*np;i++) {
		j=parseInt(i/3);
		Force[i]-=Velocity[i]*gamma*Volume[j];
		Velocity[i]+=Force[i]/(Volume[j]*rho)*dt;
		p[i]+=Velocity[i]*dt;
	}
}