

var rand = function(min,max) {
    return [ Math.random()/10,Math.random()/10,Math.random()];
}

function rotateGeometry(parms,obj){
    
	var geometry = obj;
    
	for(var i = 0 ; i < geometry.vertices.length ; i+=3){
        
		var obj = geometry.vertices[i];
        
		var parms2 = [];
		parms2.x = obj.x;
		parms2.y = obj.y;
		parms2.z = obj.z;
        
		var newMat = jMatMult(makeMatrix(parms2),makeMatrix(parms));
		//console.log(newMat);
        
		obj.x = newMat[12];
		obj.y = newMat[13];
		obj.z = newMat[14];
		
	}
    
	if(rotateParent.children.length>0){
        removal();
	}
}

function map_range(value, low1, high1, low2, high2) {
	return low2 + (high2 - low2) * (value - low1) / (high1 - low1);
}

function msin(t,l,h){
	return map_range(Math.sin(t),-1,1,l,h);
}

function mcos(t,l,h){
	return map_range(Math.cos(t),-1,1,l,h);
}

var Average = function(s){
    
	this.bucket = [];
	this.value = s;
    
	this.add = function(val){
        
		this.bucket.push(val);
        
		while(this.bucket.length>this.value){
			this.bucket.shift();
		}
	}
    
	this.avg=function(){
        
		var num = 0;
        
		for(var i = 0 ; i < this.bucket.length ; i++){
			num+=this.bucket[i];
		}
        
		return num/this.value;
	}
}

var Debug = function(){
    
	//instantiate a new debugger:
	//var d = new Debug();
	//then use d.print(variable); or d.print(variable + " " + variable2)
	//or d.print(["variable",variable,"variable2",variable2]);
    
	this.debug=true;
	this.memory = [];
	this.overflow = 1e4;
	this.flow = 0;
    
	this.printer = function(v){
        
		var st = "";
		this.flow++;
        
		if(this.flow>this.overflow){
			this.debug=false;
		}
        
		if(v instanceof Object){
			for (var i = 0 ; i < v.length ; i++){
				st+=v[i];
				st+="|";
			}
		}
		else{
			st+=v;
		}
        
		this.memory.push(st);
        
		if(this.memory.length>2)
			this.memory.shift();
        
		return st;
	}
    
	this.log = function(v){
		var st = this.printer(v);
		this.output(v);
	}
    
	this.print = function(v){
		var st = this.printer(v);
		this.output(st);
	}
    
	this.output = function(st){
		if(this.debug){
			if(st!=this.memory[0]||this.memory.length<2)
				console.log(st);
		}
	}
}

function RTSmat(x,y,z,t,s){
    
	var m = [];
	m.push(x);
	m.push(y);
	m.push(z);
	m.push(t);
	m.push(s);
    
	var r1 = jMatMult(x,y);
	var r2 = jMatMult(z,r1);
    
	var ts = jMatMult(t,s);
	
	return jMatMult(r2,ts);
}

function STRmat(x,y,z,t,s){
    
	var m = [];
	m.push(x);
	m.push(y);
	m.push(z);
	m.push(t);
	m.push(s);
    
	var r1 = jMatMult(y,x);
	var r2 = jMatMult(z,r1);
    
	var st = jMatMult(s,t);
	
	return jMatMult(st,r2);
}

var a = 0; ///errrr...

function jMatMult(A,B){
    
	var C = [];
    
	
	for(var i = 0; i < 16; i++){
		C.push(0);
	}
	
	
	C[0]=	A[0]* B[0] 	+ A[1] *B[4] 	+ A[2]* B[8] 	+ A[3]* B[12];
	C[4]=	A[4]* B[0] 	+ A[5]* B[4] 	+ A[6]* B[8] 	+ A[7]* B[12];
	C[8]=	A[8]* B[0]	+ A[9]* B[4] 	+ A[10]*B[8] 	+ A[11]*B[12];
	C[12]=	A[12]*B[0]	+ A[13]*B[4] 	+ A[14]*B[8] 	+ A[15]*B[12];
    
	C[1]=	A[0]* B[1] 	+ A[1]* B[5] 	+ A[2]* B[9] 	+ A[3]* B[13];
	C[5]=	A[4]* B[1] 	+ A[5]* B[5] 	+ A[6]* B[9] 	+ A[7]* B[13];
	C[9]=	A[8]* B[1]	+ A[9]* B[5] 	+ A[10]*B[9] 	+ A[11]*B[13];
	C[13]=	A[12]*B[1]	+ A[13]*B[5] 	+ A[14]*B[9] 	+ A[15]*B[13];
    
	C[2]=	A[0] *B[2] 	+ A[1] *B[6]  	+ A[2] *B[10] 	+ A[3] *B[14];
	C[6]=  	A[4] *B[2] 	+ A[5] *B[6]  	+ A[6] *B[10] 	+ A[7] *B[14];
	C[10]=  A[8] *B[2]	+ A[9] *B[6] 	+ A[10]*B[10] 	+ A[11]*B[14];
	C[14]=  A[12]*B[2]	+ A[13]*B[6] 	+ A[14]*B[10] 	+ A[15]*B[14];
    
	C[3]=	A[0] *B[3] 	+ A[1] *B[7] 	+ A[2] *B[11]	+ A[3] *B[15];
	C[7]=	A[4] *B[3] 	+ A[5] *B[7] 	+ A[6] *B[11]	+ A[7] *B[15];
	C[11]=	A[8] *B[3]	+ A[9] *B[7] 	+ A[10]*B[11]	+ A[11]*B[15];
	C[15]=	A[12]*B[3]	+ A[13]*B[7] 	+ A[14]*B[11]	+ A[15]*B[15];
	
	return C;
}

function recurse(ob,mat){
    
    var kidMat = mat;
    
    if(ob.parent){
        
        var tempKidMat = jMatMult(mat,ob.parent.pMatrix);
        kidMat = recurse(ob.parent,tempKidMat);
        
    }
    
    return kidMat;
}

function makeMatrix(params){
    
    var x = params.x || 0;
    var y = params.y || 0;
    var z = params.z || 0;
    var rx = params.rx || 0;
    var ry = params.ry || 0;
    var rz = params.rz || 0;
    var xs = params.sx || 1;
    var ys = params.sy || 1;
    var zs = params.sz || 1;
    
    this.ident = [ 1,0,0,0,
                  0,1,0,0,
                  0,0,1,0,
                  0,0,0,1 ];
	
	this.tMatrix = [ 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					x,y,z,1 ];
    
	this.sMatrix = [ xs,0,0,0,
					0,ys,0,0,
					0,0,zs,0,
					0,0,0,1 ];
	
	this.rXmatrix = [ 	1,0,0,0,
                     0,Math.cos(rx),-Math.sin(rx),0,
                     0,Math.sin(rx),Math.cos(rx),0,
                     0,0,0,1 ];
	this.rYmatrix = [ 	Math.cos(ry),0,Math.sin(ry),0,
                     0,1,0,0,
                     -Math.sin(ry),0,Math.cos(ry),0,
                     0,0,0,1 ];
	this.rZmatrix = [
                     Math.cos(rz),-Math.sin(rz),0,0,
                     Math.sin(rz),Math.cos(rz),0,0,
                     0,0,1,0,
                     0,0,0,1 ];
    
	return RTSmat(this.rXmatrix,this.rYmatrix,this.rZmatrix,this.tMatrix,this.sMatrix);
    
}

function updateMatrix(params){
	var x = params.x || 0;
    var y = params.y || 0;
    var z = params.z || 0;
    var rx = params.rx || 0;
    var ry = params.ry || 0;
    var rz = params.rz || 0;
    var xs = params.sx || 0;
    var ys = params.sy || 0;
    var zs = params.sz || 0;
}

Object3D = function(){};

Object3D.prototype = new Object3D();

Object3D.prototype.makeObj = function(){
    
	this.hasChildren = false;
	this.pts = [];
	this.edges = [];
	this.matrix = [];
	this.pMatrix = [];
	var tempMat = identity();
	this.tMatrix = tempMat;
	this.children = [];
	
	this.add = function(obj){
		this.children.push(obj);
		this.hasChildren = true;
		obj.parent = this;
	}
    
	this.move = function(obj,mat){
        
		obj.matrix = jMatMult(obj.tMatrix,mat);
	    var kidMat = mat;
        
	    if(obj.children.length>0){
            
	    	for(var t = 0 ; t < obj.children.length ; t++){
			    var tempKidMat = jMatMult(obj.matrix,obj.children[t].pMatrix);
			    //kidMat =
			    obj.move(obj.children[t],tempKidMat);
            }
            
	    }
        
	    return kidMat;
        
	}
    
	this.recurse = function(mat){
        
	    var kidMat = mat;
	    this.pMatrix = mat;
        
	    if(this.parent){
            
            var tempKidMat = jMatMult(mat,this.parent.pMatrix);
            kidMat = recurse(this.parent,tempKidMat);
            
	    }
        
	    this.matrix=kidMat;
	}
    
	this.move2 = function(mat){
        
		//jMatMult(this.matrix,mat);
		this.pMatrix = mat;
		var kidMat = mat;
		//obj.pMatrix = mat;
		//obj.matrix = mat;//jMatMult(obj.matrix,obj.pMatrix);
        
		try{
		    if(this.hasChildren){
		    	for(var t = 0 ; t < this.children.length ; t++){
                    //var tempKidMat = jMatMult(obj, obj.matrix);//obj.children[t].pMatrix);
                    this.children[t].move2(jmatMult(mat,this.pMatrix));
		      	}
                
		    }
		}
		catch(e){
			console.log(obj);
		}
        
	    this.matrix = jMatMult(this.matrix,this.pMatrix);
	}
    
	return this;
}

var O3D = new Object3D();
O3D.makeObj();

function makeSpiral(obj){
	obj.pts = [];
	obj.edges = [];
	for(var i = -50 ; i < 50 ; i++ ){
		obj.pts.push([Math.cos(i/Math.PI)*Math.cos(i/30),i/100,Math.sin(i/Math.PI)*Math.cos(i/30)]);
	}
    
    // THE EDGES OF A UNIT CUBE (INDEXING INTO THE VERTICES)
    
	for(var i = 0 ; i < obj.pts.length-1 ; i++ ){
		obj.edges.push([i,i+1]);
        
	}
	for(var i = 0 ; i < obj.pts.length-10 ; i++ ){
		obj.edges.push([i,i+10]);
        
	}
}

function makeSpiralSimple(obj){
	obj.pts = [];
	obj.edges = [];
    
	for(var i = -50 ; i < 50 ; i++ ){
		obj.pts.push([Math.cos(i/Math.PI)*Math.cos(i/30),i/100,Math.sin(i/Math.PI)*Math.cos(i/30)]);
	}
    
    // THE EDGES OF A UNIT CUBE (INDEXING INTO THE VERTICES)
    
	for(var i = 0 ; i < obj.pts.length-1 ; i++ ){
		obj.edges.push([i,i+1]);
        
	}
}

function makeSphere2(obj,div){
	obj.pts = [];
	obj.edges = [];
    
	for(var i = 0 ; i < div ; i++ ){
		for(var j = 0 ; j < div ; j++ ){
            
			var c = Math.sin((j/(div-1))*Math.PI);
			var b = Math.cos((j/(div-1))*Math.PI);
			
			
			obj.pts.push([c,b+.1,0]);
            
			
            
		}
		
		for(var q = 0 ; q < obj.pts.length ; q++){
            var matR = [];
            matR = makeMatrix(0,0,0,1,1,1,0,Math.PI*2/(div-1),0);
            obj.pts[q] = transform(obj.pts[q], matR);
            //console.log(i + " " + matR);
		}
	}
    
    // THE EDGES OF A UNIT CUBE (INDEXING INTO THE VERTICES)
    
	for(var i = 0 ; i < obj.pts.length-1 ; i++ ){
		if((i+1)%div!=0)
            obj.edges.push([i,i+1]);
        
        
	}
	for(var i = 0 ; i < obj.pts.length-div ; i++ ){
		//if((i+1)%20!=0)
		obj.edges.push([i,i+div]);
        
        
	}
}

function updateVerts(obj,off,sc){
	d.log(obj.vertices);
    for(t in obj.vertices){
        obj.vertices[t]+=(Math.random()-.5)*off;//(Math.random()-.5)*0.01;
        obj.vertices[t]*=sc;
    }
}

var sphwNoise = function(u,v,p) {
    
	var rand = (Math.random()*.002)+1.;
	var rad=.1;
	if(p!=undefined)
        var rad = p.r1;
    //U makes a full circle
    var theta = 2 * Math.PI * u,
    //v makes a half circle from -.5 to .5
    phi = Math.PI * (v - .5),
    
    cosT = Math.cos(theta) *rad*rand ,
    cosP = Math.cos(phi)*1,
    sinT = Math.sin(theta) *rad*rand ,
    sinP = Math.sin(phi)*rad *rand ;
    
    return [ cosT * cosP, sinT * cosP, sinP ];
}

var sph = function(u,v,p) {
    
	var rad=.1;
	if(p!=undefined)
        var rad = p.r1;
    //U makes a full circle
    var theta = 2 * Math.PI * u,
    //v makes a half circle from -.5 to .5
    phi = Math.PI * (v - .5),
    
    cosT = Math.cos(theta) *rad,
    cosP = Math.cos(phi)*1 ,
    sinT = Math.sin(theta) *rad,
    sinP = Math.sin(phi)*rad ;
    
    return [ cosT * cosP, sinT * cosP, sinP ];
}

var drop = function(u,v,p) {
    
	/*http://paulbourke.net/geometry/teardrop/
     x=.5*(1-cos(theta))*sin(theta)*cos(phi)
     y=.5*(1-cos(theta))*sin(theta)*sin(phi)
     z=cos(theta)
     where0<=phi<=2PI
     and0<=theta<=PI
     */
    
	var rad=.1;
	if(p!=undefined)
        var rad = p.r1;
    //U makes a full circle
    var theta = 2 * Math.PI * u,
    //v makes a half circle from -.5 to .5
    phi = Math.PI * (v - .5),
    
    cosT = Math.cos(theta) *rad,
    cosP = Math.cos(phi)*1 ,
    sinT = Math.sin(theta) *rad,
    sinP = Math.sin(phi)*rad ;
    
    var z = .5*(1-cosT)*sinT*-cosP;
    var y = 1*(1-cosT)*sinT*sinP;
    // if(phi>0&&phi<Math.PI*2&&theta>0&&theta<Math.PI)
    var x = cosT-(rad);
    //	else
    //	var z=sinT;
    
    //return [ cosT * cosP, sinT * cosP, Math.pow(sinP,.45) ];
    return [x,y,z ];
}

var tor = function(u,v,parms) {
    
	this.mt = 2;//mult theta
	this.mp = 2;//mult phi
	this.at= 2;//add theta
	this.ap = 2;//add phi
    if(parms!=undefined){
        this.R=parms.r1 || 1;
        this.r=parms.r2 || .5;
        this.mt=parms.mt || 2;
        this.mp=parms.mp || 2;
        this.at=parms.at || 0;
        this.ap=parms.ap || 0;
    }
    
    R=this.R;
    r=this.r;
    
    /* from wikipedia
     x(u,v) =  (a + b cos v) cos u,
     y(u,v) =  (a + b cos v) sin u,
     z(u,v) = b sin v.
     */
    
    var theta = (mt * (Math.PI + at)) * u,
    phi = (mp * (Math.PI + ap)) * v,
    
    cosT = Math.cos(theta) ,
    cosP = Math.cos(phi) ,
    sinT = Math.sin(theta) ,
    sinP = Math.sin(phi) ,
    
    x = (R+r*cosP)*cosT,
    y = (R+r*cosP)*sinT,
    z = r*sinP;
    
    return [x,y,z];
}

var tor2 = function(u,v,parms) {
    
	this.mt = 2;//mult theta
	this.mp = 2;//mult phi
	this.at= 2;//add theta
	this.ap = 2;//add phi
    if(parms!=undefined){
        this.R=parms.r1 || 1;
        this.r=parms.r2 || .5;
        this.mt=parms.mt || 2;
        this.mp=parms.mp || 2;
        this.at=parms.at || 0;
        this.ap=parms.ap || 0;
    }
    
    R=this.R;
    r=this.r;
    
    /* from wikipedia
     x(u,v) =  (a + b cos v) cos u,
     y(u,v) =  (a + b cos v) sin u,
     z(u,v) = b sin v.
     */
    
    var theta = (mt * (Math.PI + at)) * u,
    phi = (mp * (Math.PI + ap)) * v,
    
    cosT = Math.cos(theta) , 
    cosP = Math.cos(phi) ,
    sinT = Math.sin(theta) , 
    sinP = Math.sin(phi);
    
    var scx = 1,
    scy = 0;
    
    if(v<.6)
        scx=12;
    
    if(v>.9 || v<.1)
        scy=2;
    
    var x = (R+r*cosP)*cosT*scx,
    y = (R+r*cosP)*sinT*scx,
    z = scy-(r*sinP);
    
    return [x,y,z];
}


var cyl = function(u,v,parms) {
    var r = parms.r1 || 1;
    var h = parms.r2 || 1;
    var theta = 2 * Math.PI * u;
    if(v==0||v==1){
        r=0.001;
    }
    
    return [r*Math.sin(theta), v*h,r*Math.cos(theta)];
}