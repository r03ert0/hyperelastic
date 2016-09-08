/**
	Hyperelastic growth model, Roberto Toro 2015
	Main
*/

var simulationParams;	// simulation parametres
var simulation;			// simulation object
var display;			// display
var flag_running=true;

function mySimulation(params) {
    this.iter=0;
	this.time=0;// history, actually
	this.dt=0;	// time step

	this.ge=0;	// geometry
	this.me=0;	// mechanics
	this.gr=0;	// growth
	this.di=0;	// display
}

/**
Main: init simulation
*/
function initSimulation(params) {
	var def=$.Deferred();
	
	flag_running=true;
	$("#startStop").html("Stop");
	var si=new mySimulation();
	var make;
	
	// create geometry
	switch(params.geometry) {
		case "block":
			make=makeBlock;
			break;
		case "ublock":
			make=makeUBlock;
			break;
		case "ring":
			make=makeRing;
			break;
		case "sphere":
			make=makeSphere;
			break;
	}
	
	make(params)
	.then(function(ge) {
		si.ge=ge;

		// init mechanics
		si.me=initMechanics(si.ge,params);
		
		// init growth
		si.gr=initGrowth(params);
		switch(params.growth) {
			case "homogeneous":
				si.growthFunction=growHomogeneous;
				break;
			case "block border instantaneous":
				si.growthFunction=growBlockBorderInstantaneous;
				break;
			case "ring border instantaneous":
				si.growthFunction=growRingBorderInstantaneous;
				break;
			case "ring tangential instantaneous":
				si.growthFunction=growRingTangentialInstantaneous;
				break;
			case "surface homogeneous instantaneous":
				si.growthFunction=growSurfaceHomogeneousInstantaneous;
				break;
		}

		// compute time step for simulation
		si.dt=computeTimeStep(si.ge,si.me);

		if(params.T==0)
			si.growthFunction(si.ge,si.gr);
			
        // init hash array for collision detection
        if(si.ge.params["collision"])
            initHash(si.ge);

		/*
		$("#select").change(selectSimulation);
		$("#startStop").click(startStop);
		*/

		def.resolve(si);
	});
	
	return def.promise();
}

/**
computeTimeStep
*/
function computeTimeStep(ge,me) {
	var a=0;	// average mesh spacing
	var n1,n2;
	var nt=ge.nt;
	var t=ge.t;
	var p=ge.p;
	var K=me.K;
	var rho=me.rho;
	var dt;
	
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
	dt=0.2*0.1*Math.sqrt(rho*a*a/K);
	
	return dt;
}

/**
Simulation step
Request a new frame, call the render function, call the simulation functions and update
model display.
*/
function simulationStep(si) {
	
	if(flag_running==false)
		return;
		
	// advance time
	si.time+=si.dt;
	si.iter+=1;

	var Ftet=0,Flin=0;

	// add elastic forces from tetrahedra
	Ftet=tetraElasticity(si.ge,si.me);
	
	// add elastic forces from linear elements
	if(si.ge.params["fibres"])// && si.ge.params.fibres==true)
		Flin=linElasticity(si.ge,si.me);
		
	// Collision detection
	if(si.ge.params["collision"]) {
	    if(si.ge.params["surfaceOnly"]) {
    	    collision_surf_addConstraints(si.ge, si.iter);
            collision_surf_removeConstraints(si.ge, si.iter);
            collision_surf_enforceConstraints(si.ge, si.me, si.iter);
        } else
    	    collision_tetra(si.ge, si.iter);
    }

	//console.log("Ftet:",Ftet,"Flin:",Flin);
	
	move(si.ge,si.me,si.dt);
}

function pauseResume(di) {
	if(flag_running==true) {
		flag_running=false;
		$("#startStop").html("Continue");
	}
	else {
		flag_running=true;
	//	animate(param);
	//	param.growthFunc(param);
		$("#startStop").html("Stop");
	}
}

