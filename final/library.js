

//*******************************************************************************************************
///all that follows is referenced from other sources
///*************


//math utils from Ken
var PI = Math.PI;
function cos(t) { return Math.cos(t); }
function dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function floor(t) { return Math.floor(t); }

function ik(a, b, C, D) {
    
    var cc = dot(C,C),
    x = (1 + (a*a - b*b)/cc) / 2,
    y = dot(C,D)/cc;
    
    for (var i = 0 ; i < 3 ; i++)
		D[i] -= y * C[i];
    
    y = sqrt(max(0,a*a - cc*x*x) / dot(D,D));
    
    for (var i = 0 ; i < 3 ; i++)
		D[i] = x * C[i] + y * D[i];
    
}

function lerp(t, a, b) { return a + t * (b - a); }
function min(a, b) { return Math.min(a, b); }
function max(a, b) { return Math.max(a, b); }
function sCurve(t) { return t * t * (3 - t - t); }
function sin(t) { return Math.sin(t); }
function sqrt(t) { return Math.sqrt(t); }
function dist( a,b ) {
    
    var dx = a[0] - b[0];
    var dy = a[1] - b[1];
    var dz = a[2] - b[2];
    
    return dx * dx + dy * dy + dz * dz;
    
}

var Spring = function (parms){
    
    this.P = parms.P || [0,0,0],
    this.V = parms.V || [0,0,0],
    this.F = parms.F || [0,0,0], //force
    this.mass = parms.mass || 1.0,
    this.damp = parms.damp || 1.0;
    
    this.update = function(elapsed) {
        
        this.V[0] += (this.F[0] - this.P[0]) / this.mass * elapsed;
        this.P[0]  = (this.P[0] + this.V[0]) * (1 - this.damp * elapsed);
        this.V[1] += (this.F[1] - this.P[1]) / this.mass * elapsed;
        this.P[1]  = (this.P[1] + this.V[1]) * (1 - this.damp * elapsed);
        this.V[2] += (this.F[2] - this.P[2]) / this.mass * elapsed;
        this.P[2]  = (this.P[2] + this.V[2]) * (1 - this.damp * elapsed);
        
        return this.P;
    }
    
}


//-------- SOLVE TWO LINK INVERSE KINEMATICS -------------

// Given a two link joint from [0,0,0] to end effector position P,
// let link lengths be a and b, and let norm |P| = c.  Clearly a+b >= c.
//
// Problem: find a "knee" position Q such that |Q| = a and |P-Q| = b.
//
// In the case of a point on the x axis R = [c,0,0], there is a
// closed form solution S = [d,e,0], where |S| = a and |R-S| = b:
//
//    d2+e2 = a2                  -- because |S| = a
//    (c-d)2+e2 = b2              -- because |R-S| = b
//
//    c2-2cd+d2+e2 = b2           -- combine the two equations
//    c2-2cd = b2 - a2
//    c-2d = (b2-a2)/c
//    d - c/2 = (a2-b2)/c / 2
//
//    d = (c + (a2-b2/c) / 2      -- to solve for d and e.
//    e = sqrt(a2-d2)

function findD(a,b,c) {
    return Math.max(0, Math.min(a, (c + (a*a-b*b)/c) / 2));
}


function findE(a,d) { return Math.sqrt(a*a-d*d); }

// This leads to a solution to the more general problem:
//
//   (1) R = Mfwd(P)         -- rotate P onto the x axis
//   (2) Solve for S
//   (3) Q = Minv(S)         -- rotate back again

var  Mfwd = new Array(3);
for (var i=0;i<3;i++){
    Mfwd[i]=new Array(3);
}

var Minv = new Array(3);
for (var i=0;i<3;i++){
    Minv[i]=new Array(3);
}

function solve(A,B,P1,P2,D,Q) {
    var P=[];
    P[0] = P1[0]-P2[0];
    P[1] = P1[1]-P2[1];
    P[2] = P1[2]-P2[2];
    var R = new Array(3);
    defineM(P,D);
    rot(Minv,P,R);
    var d = findD(A,B,norm(R));
    var e = findE(A,d);
    var S = [d,e,0];
    rot(Mfwd,S,Q);
    Q[0]=Q[0]+P2[0];Q[1]=Q[1]+P2[1];Q[2]=Q[2]+P2[2];
    return d > 0 && d < A;
}


// Solve2 (a variation on ik/solve) returns the coordinate of the missing vertex of a triangle in the XY plane.
// a&b are leg lengths, C is the free (mouse) coorinate, R is the fixed coordinate, Q is normally [0,0,0].
// Nonzero Q inputs allow additional modification to this position which is included based on ik/solve, but I think
// this should be taken out of the argument list altogether and additional modifications should be performed
// outside the function. invbend should be set to -1 or 1 to reverse the sign on the bend angle (ie the joint will
// bow inward or outward). The distance between the fixed point and the mouse is returned.
function solve2(a,b,C,R,Q,invbend){
    var cdiff = []; //vector pointing from R to C
    cdiff[0] = C[0]-R[0]; cdiff[1] = C[1]-R[1]; cdiff[2] = C[2]-R[2];
    var clen  = norm(cdiff);
    if (clen > a+b) { //if a+b is not long enough to reach the mouse, allow legs to stretch
        Q[0] = R[0] + cdiff[0]/2; Q[1] = R[1] + cdiff[1]/2; Q[2] = R[2];
    } else { //leg length fixed to a+b
        var c0  = cdiff[0]/clen, s0  = cdiff[1]/clen; //sine/cosine from triangle
        var cA = (-a*a+b*b+clen*clen)/2/b/clen; //cos(A) => A is angle formed b/t CR and CQ from law of cosines
        var sA = invbend*Math.sqrt(1-cA*cA); //trig identity to find sin(A): cos(A)^2 + sin(A)^2 = 1
        Q[0] += R[0]+(c0*cA-s0*sA)*a; //trig identity: cos(A+B) = cos(A)*cos(B)-sin(A)*sin(B)
        Q[1] += R[1]+(s0*cA+c0*sA)*a; //trig identity: sin(A+B) = sin(A)*cos(B)+cos(A)*sin(B)
        Q[2] += R[2]+0;
    }
    return clen;
}


// If "knee" position Q needs to be as close as possible to some point D,
// then choose M such that M(D) is in the y>0 half of the z=0 plane.
//
// Given that constraint, define the forward and inverse of M as follows:

function defineM(P,D) {
    X = Minv[0]; Y = Minv[1]; Z = Minv[2];
    
    // Minv defines a coordinate system whose x axis contains P, so X = unit(P).
    
    for (var i = 0 ; i < 3 ; i++)
        X[i] = P[i];
    normalize(X);
    
    // The y axis of Minv is perpendicular to P, so Y = unit( D - X(DÃ‚Â·X) ).
    
    var dDOTx = dot(D,X);
    for (var i = 0 ; i < 3 ; i++)
        Y[i] = D[i] - dDOTx * X[i];
    normalize(Y);
    
    // The z axis of Minv is perpendicular to both X and Y, so Z = XÃƒâ€”Y.
    
    cross(X,Y,Z);
    
    // Mfwd = (Minv)T, since transposing inverts a rotation matrix.
    
    for (var i = 0 ; i < 3 ; i++) {
        Mfwd[i][0] = Minv[0][i];
        Mfwd[i][1] = Minv[1][i];
        Mfwd[i][2] = Minv[2][i];
    }
}

//------------ GENERAL VECTOR MATH SUPPORT -----------

function norm(v) { return Math.sqrt( dot(v,v) ); }

function normalize(v) {
    var normm = norm(v);
    for (var i = 0 ; i < 3 ; i++)
        v[i] /= normm;
}

function dot(a,b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

function cross(a,b,c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

function crossT(a,b) {
    var c = new Array(3);
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

function rot(M,src,dst) {
    for (var i = 0 ; i < 3 ; i++)
        dst[i] = dot(M[i],src);
}


/**
 * Spline from Tween.js, slightly optimized (and trashed)
 * http://sole.github.com/tween.js/examples/05_spline.html
 *
 * @author mrdoob / http://mrdoob.com/
 * @author alteredq / http://alteredqualia.com/
 from THREEJS*/

Point = function(a){
    
	this.x = a[0];
	this.y = a[1];
	this.z = a[2];
}

Spline = function ( points ) {
    
	this.points = points;
    
	var c = [], v3 = { x: 0, y: 0, z: 0 },
	point, intPoint, weight, w2, w3,
	pa, pb, pc, pd;
    
    
	this.getPoint = function ( k ) {
        
		point = ( this.points.length - 1 ) * k;
		intPoint = Math.floor( point );
		weight = point - intPoint;
        
		c[ 0 ] = intPoint === 0 ? intPoint : intPoint - 1;
		c[ 1 ] = intPoint;
		c[ 2 ] = intPoint  > this.points.length - 2 ? this.points.length - 1 : intPoint + 1;
		c[ 3 ] = intPoint  > this.points.length - 3 ? this.points.length - 1 : intPoint + 2;
        
		pa = this.points[ c[ 0 ] ];
		pb = this.points[ c[ 1 ] ];
		pc = this.points[ c[ 2 ] ];
		pd = this.points[ c[ 3 ] ];
        
		w2 = weight * weight;
		w3 = weight * w2;
        
		v3.x = interpolate( pa.x, pb.x, pc.x, pd.x, weight, w2, w3 );
		v3.y = interpolate( pa.y, pb.y, pc.y, pd.y, weight, w2, w3 );
		v3.z = interpolate( pa.z, pb.z, pc.z, pd.z, weight, w2, w3 );
        
		return v3;
        
	};
    
	// Catmull-Rom
    
	function interpolate( p0, p1, p2, p3, t, t2, t3 ) {
        
		var v0 = ( p2 - p0 ) * 0.5,
        v1 = ( p3 - p1 ) * 0.5;
        
		return ( 2 * ( p1 - p2 ) + v0 + v1 ) * t3 + ( - 3 * ( p1 - p2 ) - 2 * v0 - v1 ) * t2 + v0 * t + p1;
        
	};
    
};



/**************************************************************
 *	Quadratic Bezier 3D curve from THREE.js
 **************************************************************/



Quad = function ( points) {
    
	this.v0 = points[0];
	this.v1 = points[1];
	this.v2 = points[2];
    
	
    
	this.getPoint = function ( t ) {
        
		var tx, ty, tz;
        
		tx = b2( t, this.v0.x, this.v1.x, this.v2.x );
		ty = b2( t, this.v0.y, this.v1.y, this.v2.y );
		tz = b2( t, this.v0.z, this.v1.z, this.v2.z );
        
		return new Point([tx,ty,tz]);
        
	}
}




function b2p0( t, p ) {
    
    var k = 1 - t;
    return k * k * p;
    
}

function b2p1( t, p ) {return 2 * ( 1 - t ) * t * p;}

function b2p2( t, p ) {return t * t * p;}

function b2( t, p0, p1, p2 ) {return b2p0( t, p0 ) + b2p1( t, p1 ) + b2p2( t, p2 );}

function setSpan(id, str) {
    document.getElementById(id).firstChild.nodeValue = str;
}


/**
 * @fileoverview gl-matrix - High performance matrix and vector operations
 * @author Brandon Jones
 * @author Colin MacKenzie IV
 * @version 2.2.0
 */

/* Copyright (c) 2013, Brandon Jones, Colin MacKenzie IV. All rights reserved.
 
 mat4.perspective = function (out, fovy, aspect, near, far) {
 var f = 1.0 / Math.tan(fovy / 2),
 nf = 1 / (near - far);
 out[0] = f / aspect;
 out[1] = 0;
 out[2] = 0;
 out[3] = 0;
 out[4] = 0;
 out[5] = f;
 out[6] = 0;
 out[7] = 0;
 out[8] = 0;
 out[9] = 0;
 out[10] = (far + near) * nf;
 out[11] = -1;
 out[12] = 0;
 out[13] = 0;
 out[14] = (2 * far * near) * nf;
 out[15] = 0;
 return out;
 };
 */
// Ported from Stefan Gustavson's java implementation
// http://staffwww.itn.liu.se/~stegu/simplexnoise/simplexnoise.pdf
// Read Stefan's excellent paper for details on how this code works.
//
// Sean McCullough banksean@gmail.com

/**
 * You can pass in a random number generator object if you like.
 * It is assumed to have a random() method.
 */

var noise=function(t,s,b){
    var n = new SimplexNoise();
    
    this.makeNoise = function(){
        return map_range(n.noise(t,t),0,1,s,b);
    }
}

var SimplexNoise = function(r) {
	if (r == undefined) r = Math;
    this.grad3 = [[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],
                  [1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],
                  [0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1]];
    this.p = [];
    for (var i=0; i<256; i++) {
        this.p[i] = Math.floor(r.random()*256);
    }
    // To remove the need for index wrapping, double the permutation table length
    this.perm = [];
    for(var i=0; i<512; i++) {
		this.perm[i]=this.p[i & 255];
	}
    
    // A lookup table to traverse the simplex around a given point in 4D.
    // Details can be found where this table is used, in the 4D noise method.
    this.simplex = [
                    [0,1,2,3],[0,1,3,2],[0,0,0,0],[0,2,3,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,2,3,0],
                    [0,2,1,3],[0,0,0,0],[0,3,1,2],[0,3,2,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,3,2,0],
                    [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                    [1,2,0,3],[0,0,0,0],[1,3,0,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,3,0,1],[2,3,1,0],
                    [1,0,2,3],[1,0,3,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,0,3,1],[0,0,0,0],[2,1,3,0],
                    [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
                    [2,0,1,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,0,1,2],[3,0,2,1],[0,0,0,0],[3,1,2,0],
                    [2,1,0,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,1,0,2],[0,0,0,0],[3,2,0,1],[3,2,1,0]];
};

SimplexNoise.prototype.dot = function(g, x, y) {
	return g[0]*x + g[1]*y;
};

SimplexNoise.prototype.noise = function(xin, yin) {
    var n0, n1, n2; // Noise contributions from the three corners
    // Skew the input space to determine which simplex cell we're in
    var F2 = 0.5*(Math.sqrt(3.0)-1.0);
    var s = (xin+yin)*F2; // Hairy factor for 2D
    var i = Math.floor(xin+s);
    var j = Math.floor(yin+s);
    var G2 = (3.0-Math.sqrt(3.0))/6.0;
    var t = (i+j)*G2;
    var X0 = i-t; // Unskew the cell origin back to (x,y) space
    var Y0 = j-t;
    var x0 = xin-X0; // The x,y distances from the cell origin
    var y0 = yin-Y0;
    // For the 2D case, the simplex shape is an equilateral triangle.
    // Determine which simplex we are in.
    var i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
    if(x0>y0) {i1=1; j1=0;} // lower triangle, XY order: (0,0)->(1,0)->(1,1)
    else {i1=0; j1=1;}      // upper triangle, YX order: (0,0)->(0,1)->(1,1)
    // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
    // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
    // c = (3-sqrt(3))/6
    var x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
    var y1 = y0 - j1 + G2;
    var x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
    var y2 = y0 - 1.0 + 2.0 * G2;
    // Work out the hashed gradient indices of the three simplex corners
    var ii = i & 255;
    var jj = j & 255;
    var gi0 = this.perm[ii+this.perm[jj]] % 12;
    var gi1 = this.perm[ii+i1+this.perm[jj+j1]] % 12;
    var gi2 = this.perm[ii+1+this.perm[jj+1]] % 12;
    // Calculate the contribution from the three corners
    var t0 = 0.5 - x0*x0-y0*y0;
    if(t0<0) n0 = 0.0;
    else {
        t0 *= t0;
        n0 = t0 * t0 * this.dot(this.grad3[gi0], x0, y0);  // (x,y) of grad3 used for 2D gradient
    }
    var t1 = 0.5 - x1*x1-y1*y1;
    if(t1<0) n1 = 0.0;
    else {
        t1 *= t1;
        n1 = t1 * t1 * this.dot(this.grad3[gi1], x1, y1);
    }
    var t2 = 0.5 - x2*x2-y2*y2;
    if(t2<0) n2 = 0.0;
    else {
        t2 *= t2;
        n2 = t2 * t2 * this.dot(this.grad3[gi2], x2, y2);
    }
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to return values in the interval [-1,1].
    return 70.0 * (n0 + n1 + n2);
};

// 3D simplex noise
SimplexNoise.prototype.noise3d = function(xin, yin, zin) {
    var n0, n1, n2, n3; // Noise contributions from the four corners
    // Skew the input space to determine which simplex cell we're in
    var F3 = 1.0/3.0;
    var s = (xin+yin+zin)*F3; // Very nice and simple skew factor for 3D
    var i = Math.floor(xin+s);
    var j = Math.floor(yin+s);
    var k = Math.floor(zin+s);
    var G3 = 1.0/6.0; // Very nice and simple unskew factor, too
    var t = (i+j+k)*G3;
    var X0 = i-t; // Unskew the cell origin back to (x,y,z) space
    var Y0 = j-t;
    var Z0 = k-t;
    var x0 = xin-X0; // The x,y,z distances from the cell origin
    var y0 = yin-Y0;
    var z0 = zin-Z0;
    // For the 3D case, the simplex shape is a slightly irregular tetrahedron.
    // Determine which simplex we are in.
    var i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
    var i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
    if(x0>=y0) {
        if(y0>=z0)
        { i1=1; j1=0; k1=0; i2=1; j2=1; k2=0; } // X Y Z order
        else if(x0>=z0) { i1=1; j1=0; k1=0; i2=1; j2=0; k2=1; } // X Z Y order
        else { i1=0; j1=0; k1=1; i2=1; j2=0; k2=1; } // Z X Y order
    }
    else { // x0<y0
        if(y0<z0) { i1=0; j1=0; k1=1; i2=0; j2=1; k2=1; } // Z Y X order
        else if(x0<z0) { i1=0; j1=1; k1=0; i2=0; j2=1; k2=1; } // Y Z X order
        else { i1=0; j1=1; k1=0; i2=1; j2=1; k2=0; } // Y X Z order
    }
    // A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
    // a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and
    // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where
    // c = 1/6.
    var x1 = x0 - i1 + G3; // Offsets for second corner in (x,y,z) coords
    var y1 = y0 - j1 + G3;
    var z1 = z0 - k1 + G3;
    var x2 = x0 - i2 + 2.0*G3; // Offsets for third corner in (x,y,z) coords
    var y2 = y0 - j2 + 2.0*G3;
    var z2 = z0 - k2 + 2.0*G3;
    var x3 = x0 - 1.0 + 3.0*G3; // Offsets for last corner in (x,y,z) coords
    var y3 = y0 - 1.0 + 3.0*G3;
    var z3 = z0 - 1.0 + 3.0*G3;
    // Work out the hashed gradient indices of the four simplex corners
    var ii = i & 255;
    var jj = j & 255;
    var kk = k & 255;
    var gi0 = this.perm[ii+this.perm[jj+this.perm[kk]]] % 12;
    var gi1 = this.perm[ii+i1+this.perm[jj+j1+this.perm[kk+k1]]] % 12;
    var gi2 = this.perm[ii+i2+this.perm[jj+j2+this.perm[kk+k2]]] % 12;
    var gi3 = this.perm[ii+1+this.perm[jj+1+this.perm[kk+1]]] % 12;
    // Calculate the contribution from the four corners
    var t0 = 0.6 - x0*x0 - y0*y0 - z0*z0;
    if(t0<0) n0 = 0.0;
    else {
        t0 *= t0;
        n0 = t0 * t0 * this.dot(this.grad3[gi0], x0, y0, z0);
    }
    var t1 = 0.6 - x1*x1 - y1*y1 - z1*z1;
    if(t1<0) n1 = 0.0;
    else {
        t1 *= t1;
        n1 = t1 * t1 * this.dot(this.grad3[gi1], x1, y1, z1);
    }
    var t2 = 0.6 - x2*x2 - y2*y2 - z2*z2;
    if(t2<0) n2 = 0.0;
    else {
        t2 *= t2;
        n2 = t2 * t2 * this.dot(this.grad3[gi2], x2, y2, z2);
    }
    var t3 = 0.6 - x3*x3 - y3*y3 - z3*z3;
    if(t3<0) n3 = 0.0;
    else {
        t3 *= t3;
        n3 = t3 * t3 * this.dot(this.grad3[gi3], x3, y3, z3);
    }
    // Add contributions from each corner to get the final noise value.
    // The result is scaled to stay just inside [-1,1]
    return 32.0*(n0 + n1 + n2 + n3);
};

var makeTexture = function(gl, obj, name, value) {
    var image = new Image();
    image.onload = function() {
        var gl = this.gl;
        gl.bindTexture(gl.TEXTURE_2D, obj.textures[name] = gl.createTexture());
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, this);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
        gl.generateMipmap(gl.TEXTURE_2D);
        gl.bindTexture(gl.TEXTURE_2D, null);
    }
    image.gl = gl;
    image.obj = obj;
    image.src = value;
}

var circle = function(divs) {
    var verts = [];
    for (var u = 0; u < divs; u++) {
        var angle = Math.PI * 2 / divs * u;
        verts.push(Math.cos(angle), 0, Math.sin(angle));
    }
    return verts;
}