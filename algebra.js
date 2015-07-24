/**
	Hyperelastic growth model, Roberto Toro 2015
	Algebra
*/

/**
Matrix
*/
function matrix(a,b,c,d,e,f,g,h,i) {
	this.a=a;
	this.b=b;
	this.c=c;
	this.d=d;
	this.e=e;
	this.f=f;
	this.g=g;
	this.h=h;
	this.i=i;
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