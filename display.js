/**
 * @page Display
 * Hyperelastic growth model, Roberto Toro 2015
 * Display
 *
 * Depends on geometry.js, three.js and simulation.js
 */

function myDisplay() {
    this.renderer=null;        // three.js renderer
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
 * @function initDisplay
 * @description Initialises the three.js render
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
 * @function initMesh
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
        var Z=12;
        di.camera = new THREE.OrthographicCamera( -Z*w/h,Z*w/h,Z,-Z,1, 100);
    }
    di.camera.position.z = 10;
    di.camera.position.y=-20;
    di.camera.lookAt(di.scene.position);
    di.scene.add(di.camera);
    
    di.trackball = new THREE.TrackballControls(di.camera,di.renderer.domElement);

    di.geometry=new THREE.Geometry();
    
    if(params.surfaceOnly) {
        var nf=si.ge.nf;
        var f=si.ge.f;
        for(n=0;n<np;n++)
            di.geometry.vertices.push(new THREE.Vector3(p[3*n+0],p[3*n+1],p[3*n+2]));
        for(n=0;n<nf;n++) {
            v1=f[n*3+0];
            v2=f[n*3+1];
            v3=f[n*3+2];
            di.geometry.faces.push(new THREE.Face3(v1,v2,v3));
        }
    } else {
        for(n=0;n<np;n++)
            di.geometry.vertices.push(new THREE.Vector3(p[3*n+0],p[3*n+1],p[3*n+2]));
        for(n=0;n<nt;n++) {
            for(i=0;i<12;i+=3) {
                v1=t[n*4+indx[i+0]];
                v2=t[n*4+indx[i+1]];
                v3=t[n*4+indx[i+2]];
                di.geometry.faces.push(new THREE.Face3(v1,v2,v3));
            }
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
    
    $("#txt").remove();
    
    if(simulationParams.showVertexNumbers || simulationParams.showTriangleNumbers) {
        $("body").append("<div id='txt' style='position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none'>");
        if(simulationParams.showVertexNumbers) {
            for(i=0;i<simulation.ge.np;i++)
                $("#txt").append("<span class='v' id=v"+i+" style='position:absolute'>"+i+"</span>");
        }
        if(simulationParams.showTriangleNumbers) {
            for(i=0;i<simulation.ge.nf;i++)
                $("#txt").append("<span class='t' id=t"+i+" style='position:absolute'>"+i+"</span>");
        }
    }
}

/**
 * @function updateMesh
 */
function updateMesh(di,si) {
    var np=si.ge.np;
    var p=si.ge.p;
    var i;

    for(i=0;i<np;i++) {
        di.geometry.vertices[i].x=p[3*i+0];
        di.geometry.vertices[i].y=p[3*i+1];
        di.geometry.vertices[i].z=p[3*i+2];
    }
    di.geometry.verticesNeedUpdate=true;
}

/**
 * @function updateRingLines
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
 * @function animate
 * @description Animation function. Request a new frame, call the render function, call the simulation functions and update model display
 */
function animate() {
    requestAnimationFrame(animate);
    display.renderer.render(display.scene,display.camera);
    display.trackball.update();

    simulationStep(simulation);
    
    if(simulationParams.colormap=="deformation") {
        if(simulationParams.geometry=="block")
            block_deformationColor(simulation.ge);
        else if(simulationParams.geometry=="ublock")
            block_deformationColor(simulation.ge);
        else if(simulationParams.geometry=="ring")
            ring_deformationColor(simulation.ge);
        else if(simulationParams.geometry=="surface")
            surface_deformationColor(simulation.ge);
    }

    updateMesh(display,simulation);

    display.geometry.computeFaceNormals();
    display.geometry.computeVertexNormals();
    display.geometry.normalsNeedUpdate = true;

    display.mesh.needsUpdate=true;
    display.geometry.colorsNeedUpdate=true;
    
    $("#log").html("t: "+simulation.time.toFixed(4)+" Ue:"+simulation.me.Ue.toFixed(2));

    if(simulation.iter%10!==0)
        return;
        
    // display collision information
    if(simulationParams.showVertexNumbers || simulationParams.showTetraNumbers) {
        var i,p,P=simulation.ge.p,F=simulation.ge.f;
        
        // show vertex indices
        if(simulationParams.showVertexNumbers) {
            for(i=0;i<simulation.ge.np;i++) {
                p=screenXY(P[3*i+0],P[3*i+1],P[3*i+2]);
                $("#v"+i).css({left:p.x,top:p.y});
            }
        }

        // show triangle indices
        if(simulationParams.showTriangleNumbers) {
            for(i=0;i<simulation.ge.nf;i++) {
                p=[
                    (P[3*F[3*i+0]+0]+P[3*F[3*i+1]+0]+P[3*F[3*i+2]+0])/3,
                    (P[3*F[3*i+0]+1]+P[3*F[3*i+1]+1]+P[3*F[3*i+2]+1])/3,
                    (P[3*F[3*i+0]+2]+P[3*F[3*i+1]+2]+P[3*F[3*i+2]+2])/3
                ];
                p=screenXY(p[0],p[1],p[2]);
                $("#t"+i).css({left:p.x,top:p.y});
            }
        }
    }

    // show collision response constraints
    if(simulationParams.collision) {
        var ci;
        var p, pi;
        var a, b, c;
        var fi, fi0, fi1, fi2;
        var q0, q1, q2, q;
        for(ci in Collision.resp) {
            pi=Collision.resp[ci].pIndex;
            p=[simulation.ge.p[3*pi+0],simulation.ge.p[3*pi+1],simulation.ge.p[3*pi+2]];
            a=Collision.resp[ci].a;
            b=Collision.resp[ci].b;
            c=Collision.resp[ci].c;
            fi=Collision.resp[ci].qfIndex;
            fi0=simulation.ge.f[3*fi+0];
            fi1=simulation.ge.f[3*fi+1];
            fi2=simulation.ge.f[3*fi+2];
            q0=[simulation.ge.p[3*fi0+0],simulation.ge.p[3*fi0+1],simulation.ge.p[3*fi0+2]];
            q1=[simulation.ge.p[3*fi1+0],simulation.ge.p[3*fi1+1],simulation.ge.p[3*fi1+2]];
            q2=[simulation.ge.p[3*fi2+0],simulation.ge.p[3*fi2+1],simulation.ge.p[3*fi2+2]];
            q=add(add(scale(q0,a),scale(q1,b)),scale(q2,c));
            if(Collision.resp[ci].l==undefined) {
                var lg=new THREE.Geometry();
                lg.vertices.push(new THREE.Vector3(p[0],p[1],p[2]));
                lg.vertices.push(new THREE.Vector3(q[0],q[1],q[2]));
                var lm=new THREE.LineBasicMaterial({color:0xff0000,linewidth:1});
                Collision.resp[ci].l=new THREE.Line(lg,lm);
                display.scene.add(Collision.resp[ci].l);
                console.log("new line");
            }
            Collision.resp[ci].l.geometry.vertices[0].set(p[0],p[1],p[2]);
            Collision.resp[ci].l.geometry.vertices[1].set(q[0],q[1],q[2]);
            Collision.resp[ci].l.geometry.verticesNeedUpdate=true;
        }
    }
}
/**
 * @function screenXY
 */
function screenXY(x,y,z) {
var vector = new THREE.Vector3(x,y,z);
var canvas = display.renderer.domElement;

//vector.set( 1, 2, 3 );

// map to normalized device coordinate (NDC) space
vector.project( display.camera );

// map to 2D screen space
    return {
        x: Math.round( (   vector.x + 1 ) * canvas.width  / 2 ),
        y: Math.round( ( - vector.y + 1 ) * canvas.height / 2 )
    }
}
/**
 * @function bdarkb2r
 * @description Color map: from dark blue to dark red for negative to positive values
 * @param {float} x Value between -1 and 1. Values outside the boundaries are clamped.
 */
function bdarkb2r(x) {
    var    r,g,b;

    if(x<-1)         {    r=0;        g=0;        b=0.5;    }
    else if(x<-0.5) {    r=0;        g=0;        b=1.5+x;}
    else if(x<0)     {    r=2*(x+0.5);g=2*(x+0.5);b=1;    }
    else if(x<0.5)    {    r=1;        g=1-2*x;    b=1-2*x;}
    else if(x<1)    {    r=1.5-x;    g=0;        b=0;    }
    else             {    r=0.5;        g=0;        b=0;    }
    return {red:r,green:b,blue:g};
}
/**
 * @function black2white
 * @description Color map: from black to white for negative to positive values.
 * @param {float} x Value between -1 and 1. Values outside the boundaries are clamped.
 */
function black2white(x) {
    var    r,g,b,w;

    r=(x+1)/2;
    g=(x+1)/2;
    b=(x+1)/2;
    w=1-Math.abs(x)*0.0;

    if(x<0) { g*=w;b*=w;}
    if(x>0) { r*=w;g*=w;}
    return {red:r,green:b,blue:g*0.95};
}
/**
 * @function ring_deformationText
 * @description Display deformation of a ring geometry as text
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
    var    tminJ,tmaxJ;

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
 * @function ring_deformationColor
 * @description Map deformation of a ring geometry into colors
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
    var    tminJ,tmaxJ,color;

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
 * @function block_deformationColor
 * @description Map deformation of a block geometry into colors
 */
var minJ=0,maxJ=1;
function block_deformationColor(ge) {
    var nw=ge.nw;
    var nh=ge.nh;
    var nd=ge.nd;
    var p=ge.p;
    var r=ge.r;
    var t=ge.t;
    var i,j,k,l,m,n;
    var numbox;
    var J,vol0,vol1;
    var    tminJ,tmaxJ,color;

    // tetrahedral boxes
    for(i=0;i<nw-1;i++)
    for(j=0;j<nh-1;j++)
    for(k=0;k<nd-1;k++)
    {
        // box element index
        numbox=i*(nh-1)*(nd-1)+j*(nd-1)+k;
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
            }
        }
    }
    minJ=tminJ;
    maxJ=tmaxJ;
}
/*
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
 * @function surface_deformationColor
 * @description Map deformation of a surface geometry into colors.
 */
function surface_deformationColor(ge) {
    var i,j,k,l,m,n;
    var numbox,numtet,im,ir;
    var J,vol0,vol1;
    var    tmp=0,mintmp,maxtmp,val;
    var vval,vval2,nval;
    
    var nsp=ge.np/2;

    vval=new Float32Array(nsp);
    vval2=new Float32Array(nsp);
    nval=new Uint8Array(nsp);
    for(i=0;i<nsp;i++)
        vval[i]=0;
    for(i=0;i<nsp;i++)
        nval[i]=0;

    // tetrahedral prisms
    for(i=0;i<ge.nt;i+=3)
    {
        vol0=0;
        vol1=0;
        for(l=0;l<3;l++) {
            n=i+l;

            vol0+=tetraVol(
                        ge.r[3*(4*n+0)+0],ge.r[3*(4*n+0)+1],ge.r[3*(4*n+0)+2],
                        ge.r[3*(4*n+1)+0],ge.r[3*(4*n+1)+1],ge.r[3*(4*n+1)+2],
                        ge.r[3*(4*n+2)+0],ge.r[3*(4*n+2)+1],ge.r[3*(4*n+2)+2],
                        ge.r[3*(4*n+3)+0],ge.r[3*(4*n+3)+1],ge.r[3*(4*n+3)+2]);
            vol1+=tetraVol(
                        ge.p[3*t[4*n+0]+0],ge.p[3*t[4*n+0]+1],ge.p[3*t[4*n+0]+2],
                        ge.p[3*t[4*n+1]+0],ge.p[3*t[4*n+1]+1],ge.p[3*t[4*n+1]+2],
                        ge.p[3*t[4*n+2]+0],ge.p[3*t[4*n+2]+1],ge.p[3*t[4*n+2]+2],
                        ge.p[3*t[4*n+3]+0],ge.p[3*t[4*n+3]+1],ge.p[3*t[4*n+3]+2]);
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
                vval[ge.f[n*4+m].a]+=val;
                vval[ge.f[n*4+m].b]+=val;
                vval[ge.f[n*4+m].c]+=val;
                nval[ge.f[n*4+m].a]+=1;
                nval[ge.f[n*4+m].b]+=1;
                nval[ge.f[n*4+m].c]+=1;
            }
        }
    }
    maxJ=maxtmp;
    minJ=mintmp;

    // from face deformation to vertex colors
    for(i=0;i<nsp;i++) {
        vval[i]=vval[i]/nval[i];
    }
    for(j=0;j<5;j++) {
        for(i=0;i<nsp;i++) {
            vval2[i]=0;
            nval[i]=0;
        }
        for(i=0;i<ge.f.length;i++) {
            vval2[ge.f[i].a]+=vval[ge.f[i].b]+vval[ge.f[i].c];
            vval2[ge.f[i].b]+=vval[ge.f[i].c]+vval[ge.f[i].a];
            vval2[ge.f[i].c]+=vval[ge.f[i].a]+vval[ge.f[i].b];
            nval[ge.f[i].a]+=2;
            nval[ge.f[i].b]+=2;
            nval[ge.f[i].c]+=2;
        }
        for(i=0;i<nsp;i++) {
            vval[i]=vval2[i]/nval[i];
        }
    }
    mintmp=vval[0];
    maxtmp=vval[0];
    for(i=0;i<nsp;i++) {
        if(vval[i]>maxtmp)
            maxtmp=vval[i];
        if(vval[i]<mintmp)
            mintmp=vval[i];
    }
    for(i=0;i<nsp;i++)
        vval[i]=(vval[i]-mintmp)/(maxtmp-mintmp);    
    for(i=0;i<nsp;i++)
        color[i].setRGB(vval[i],vval[i],vval[i]);
    for(i=0;i<ge.f.length;i++) {
        ge.f[i].vertexColors[0]=color[ge.f[i].a];
        ge.f[i].vertexColors[1]=color[ge.f[i].b];
        ge.f[i].vertexColors[2]=color[ge.f[i].c];
    }
}
