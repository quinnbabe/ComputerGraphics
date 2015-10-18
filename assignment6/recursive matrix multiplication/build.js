


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

var a = 0;


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
    this.children = [];
    
    this.add = function(obj){
        this.children.push(obj);
        this.hasChildren = true;
        obj.parent = this;
    }
    
    return this;
}

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
