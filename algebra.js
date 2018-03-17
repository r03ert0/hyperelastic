/**
 * @page Algebra
 * Algebra functions that make your life better.
 */

/**
 * @library hyperelastic
 * @version 0.0.1
 * @brief Hyperelastic growth simulator in javascript
 */

/**
 * @function matrix
 * @description Makes a 3x3 matrix object out of 9 float values
 * @param {float} a 1st row, 1st column
 * @param {float} b
 * @param {float} c
 * @param {float} d 2nd row, 1st column
 * @param {float} e
 * @param {float} f
 * @param {float} g 3rd row, 1st column
 * @param {float} h
 * @param {float} i
 * @return {Matrix} A matrix object
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
 * @function invert
 * @description Matrix inversion.
 * @param {Matrix} m A 3x3 matrix represented as a vector
 * @return {Matrix} The inverse of matrix m
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
 * @function mulInvMatVec
 * @description Multiply matrix inverse by vector
 * @param {Matrix} m A 3x3 matrix represented as a vector, first row first
 * @param {Vector} p 3x1 vector
 */
function mulInvMatVec(m, p)
{
    var    det;
    var    a,b,c,d,e,f,g,h,i;
    
    det = m.b*m.f*m.g + m.c*m.d*m.h + m.a*m.e*m.i - m.c*m.e*m.g - m.a*m.f*m.h - m.b*m.d*m.i;
    
    a=(m.e*m.i - m.f*m.h);
    b=(m.c*m.h - m.b*m.i);
    c=(m.b*m.f - m.c*m.e);
    
    d=(m.f*m.g - m.d*m.i);
    e=(m.a*m.i - m.c*m.g);
    f=(m.c*m.d - m.a*m.f);
    
    g=(m.d*m.h - m.e*m.g);
    h=(m.b*m.g - m.a*m.h);
    i=(m.a*m.e - m.b*m.d);
    
    return [ (p[0]*a + p[1]*d + p[2]*g)/det,
             (p[0]*b + p[1]*e + p[2]*h)/det,
             (p[0]*c + p[1]*f + p[2]*i)/det];
}
/**
 * @function determinant
 * @description Matrix determinant.
 * @param {Matrix} a A 3x3 matrix represented as a vector
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
 * @description Vector addition.
 * @function add
 * @param {Vector} a 3x1 vector
 * @param {Vector} b 3x1 vector
 */
function add(a,b) {
    return [a[0]+b[0],a[1]+b[1],a[2]+b[2]];
}
/**
 * @description subtract
 * @function subtract
 */
function subtract(a,b) {
    return [a[0]-b[0],a[1]-b[1],a[2]-b[2]];
}
/**
 * @description norm
 * @function norm
 */
function norm(a) {
    return Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

function dot(a,b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

function scale(a,t) {
    return [a[0]*t,a[1]*t,a[2]*t];
}
/**
 * @description addMat
 * @function addMat
 */
function addMat(a,b) {
    return {
        a:a.a+b.a, b:a.b+b.b, c:a.c+b.c,
        d:a.d+b.d, e:a.e+b.e, f:a.f+b.f,
        g:a.g+b.g, h:a.h+b.h, i:a.i+b.i};
}
/**
 * @description subMat
 * @function subMat
 */
function subMat(a,b) {
    return {
        a:a.a-b.a, b:a.b-b.b, c:a.c-b.c,
        d:a.d-b.d, e:a.e-b.e, f:a.f-b.f,
        g:a.g-b.g, h:a.h-b.h, i:a.i-b.i};
}
/**
 * @description cross
 * @function cross
 */
function cross(a,b) {
    return [a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]];
}
/**
 * @description transpose
 * @function transpose
 */
function transpose(m) {
    return {a:m.a,b:m.d,c:m.g,
            d:m.b,e:m.e,f:m.h,
            g:m.c,h:m.f,i:m.i};
}
/**
 * @description trace
 * @function trace
 */
function trace(m) {
    return m.a+m.e+m.i;
}
/**
 * @description Make a matrix out of 3 vectors, one per column
 * @function vecs2Mat
 * @param {Vector} a 3x1 vector
 * @param {Vector} b 3x1 vector
 * @param {Vector} c 3x1 vector
 */
function vecs2Mat(a, b, c)
{
    return { a:a[0],b:b[0],c:c[0],
             d:a[1],e:b[1],f:c[1],
             g:a[2],h:b[2],i:c[2] };
}

/**
 * @description mulMat
 * @function mulMat
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
 * @description mulMatVec
 * @function mulMatVec
 */
function mulMatVec(m,a) {
    return [
        m.a*a[0]+m.b*a[1]+m.c*a[2],
        m.d*a[0]+m.e*a[1]+m.f*a[2],
        m.g*a[0]+m.h*a[1]+m.i*a[2]];
}
/**
 * @description mulVecSca
 * @function mulVecSca
 */
function mulVecSca(a,b) {
    return [a[0]*b,a[1]*b,a[2]*b];
}
/**
 * @description mulMatSca
 * @function mulMatSca
 */
function mulMatSca(m,a) {
    return {
        a:m.a*a, b:m.b*a, c:m.c*a,
        d:m.d*a, e:m.e*a, f:m.f*a,
        g:m.g*a, h:m.h*a, i:m.i*a};
}
/**
 * @function solidAngle
 * @description Compute the solid angle of a tetrahedron. Reference Jacobson et al (2013) "Robust Inside-Outside Segmentation using Generalized Winding Numbers"
 * @param {Vector} a 1st vertex of tetrahedron
 * @param {Vector} b 2nd vertex of tetrahedron
 * @param {Vector} c 3rd vertex of tetrahedron
 * @param {Vector} d 4th vertex of tetrahedron
 */
function solidAngle(a, b, c, d) {
    var A=[a[0]-d[0],a[1]-d[1],a[2]-d[2]];
    var B=[b[0]-d[0],b[1]-d[1],b[2]-d[2]];
    var C=[c[0]-d[0],c[1]-d[1],c[2]-d[2]];
    var detABC=A[0]*B[1]*C[2] + B[0]*C[1]*A[2] + C[0]*A[1]*B[2] - A[0]*C[1]*B[2] - B[0]*A[1]*C[2] - C[0]*B[1]*A[2];
    var na=Math.sqrt(A[0]*A[0]+A[1]*A[1]+A[2]*A[2]);
    var nb=Math.sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
    var nc=Math.sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
    var divisor = na*nb*nc + (A[0]*B[0] + A[1]*B[1] + A[2]*B[2])*nc + (B[0]*C[0] + B[1]*C[1] + B[2]*C[2])*na + (C[0]*A[0] + C[1]*A[1] + C[2]*A[2])*nb;
    var sabc=2*Math.atan(detABC/divisor);
    
    return sabc;
}
/**
 * @description printMat
 * @function printMat
 */
function printMat(M,name) {
    console.log(name+":",
                M.a+" "+M.b+" "+M.c,
                M.d+" "+M.e+" "+M.f,
                M.g+" "+M.h+" "+M.i);
}
/**
 * @description printVec
 * @function printVec
 */
function printVec(V,name) {
    console.log(name+": "+V[0]+","+V[1]+","+V[2]);
}