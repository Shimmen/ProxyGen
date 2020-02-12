#pragma once

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*   2020-02-12: rewrote to use glm instead. (S. Moos)  */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

#include "mathkit.h"

#define FINDMINMAX(x0, x1, x2, min, max) \
    min = max = x0;                      \
    if (x1 < min)                        \
        min = x1;                        \
    if (x1 > max)                        \
        max = x1;                        \
    if (x2 < min)                        \
        min = x2;                        \
    if (x2 > max)                        \
        max = x2;

/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb)                 \
    p0 = a * v0.y - b * v0.z;                      \
    p2 = a * v2.y - b * v2.z;                      \
    if (p0 < p2) {                                 \
        min = p0;                                  \
        max = p2;                                  \
    } else {                                       \
        min = p2;                                  \
        max = p0;                                  \
    }                                              \
    rad = fa * boxHalfSize.y + fb * boxHalfSize.z; \
    if (min > rad || max < -rad)                   \
        return false;

#define AXISTEST_X2(a, b, fa, fb)                  \
    p0 = a * v0.y - b * v0.z;                      \
    p1 = a * v1.y - b * v1.z;                      \
    if (p0 < p1) {                                 \
        min = p0;                                  \
        max = p1;                                  \
    } else {                                       \
        min = p1;                                  \
        max = p0;                                  \
    }                                              \
    rad = fa * boxHalfSize.y + fb * boxHalfSize.z; \
    if (min > rad || max < -rad)                   \
        return false;

/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb)                 \
    p0 = -a * v0.x + b * v0.z;                     \
    p2 = -a * v2.x + b * v2.z;                     \
    if (p0 < p2) {                                 \
        min = p0;                                  \
        max = p2;                                  \
    } else {                                       \
        min = p2;                                  \
        max = p0;                                  \
    }                                              \
    rad = fa * boxHalfSize.x + fb * boxHalfSize.z; \
    if (min > rad || max < -rad)                   \
        return false;

#define AXISTEST_Y1(a, b, fa, fb)                  \
    p0 = -a * v0.x + b * v0.z;                     \
    p1 = -a * v1.x + b * v1.z;                     \
    if (p0 < p1) {                                 \
        min = p0;                                  \
        max = p1;                                  \
    } else {                                       \
        min = p1;                                  \
        max = p0;                                  \
    }                                              \
    rad = fa * boxHalfSize.x + fb * boxHalfSize.z; \
    if (min > rad || max < -rad)                   \
        return false;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)                 \
    p1 = a * v1.x - b * v1.y;                      \
    p2 = a * v2.x - b * v2.y;                      \
    if (p2 < p1) {                                 \
        min = p2;                                  \
        max = p1;                                  \
    } else {                                       \
        min = p1;                                  \
        max = p2;                                  \
    }                                              \
    rad = fa * boxHalfSize.x + fb * boxHalfSize.y; \
    if (min > rad || max < -rad)                   \
        return false;

#define AXISTEST_Z0(a, b, fa, fb)                  \
    p0 = a * v0.x - b * v0.y;                      \
    p1 = a * v1.x - b * v1.y;                      \
    if (p0 < p1) {                                 \
        min = p0;                                  \
        max = p1;                                  \
    } else {                                       \
        min = p1;                                  \
        max = p0;                                  \
    }                                              \
    rad = fa * boxHalfSize.x + fb * boxHalfSize.y; \
    if (min > rad || max < -rad)                   \
        return false;

bool planeBoxOverlap(vec3 normal, vec3 vert, vec3 maxbox)
{
    vec3 vmin, vmax;

    for (int q = 0; q <= 2; q++) {
        float v = vert[q];
        if (normal[q] > 0.0f) {
            vmin[q] = -maxbox[q] - v;
            vmax[q] = maxbox[q] - v;
        } else {
            vmin[q] = maxbox[q] - v;
            vmax[q] = -maxbox[q] - v;
        }
    }

    if (dot(normal, vmin) > 0.0f)
        return false;

    if (dot(normal, vmax) >= 0.0f)
        return true;

    return false;
}

bool triangleBoxIntersection(vec3 boxCenter, vec3 boxHalfSize, vec3 triVerts[3])
{
    // use separating axis theorem to test overlap between triangle and box
    // need to test for overlap in these directions:
    //  1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
    //     we do not even need to test these)
    //  2) normal of the triangle
    //  3) crossproduct(edge from tri, {x,y,z}-direction)
    //     this gives 3x3=9 more tests

    // Move everything so that the boxcenter is in (0,0,0)
    vec3 v0 = triVerts[0] - boxCenter;
    vec3 v1 = triVerts[1] - boxCenter;
    vec3 v2 = triVerts[2] - boxCenter;

    // Compute triangle edges
    vec3 e0 = v1 - v0;
    vec3 e1 = v2 - v1;
    vec3 e2 = v0 - v2;

    //
    // Bullet 3:
    //

    float rad;
    float min, max;
    float p0, p1, p2;

    vec3 fe0 = abs(e0);
    AXISTEST_X01(e0.z, e0.y, fe0.z, fe0.y);
    AXISTEST_Y02(e0.z, e0.x, fe0.z, fe0.x);
    AXISTEST_Z12(e0.y, e0.x, fe0.y, fe0.x);

    vec3 fe1 = abs(e1);
    AXISTEST_X01(e1.z, e1.y, fe1.z, fe1.y);
    AXISTEST_Y02(e1.z, e1.x, fe1.z, fe1.x);
    AXISTEST_Z0(e1.y, e1.x, fe1.y, fe1.x);

    vec3 fe2 = abs(e2);
    AXISTEST_X2(e2.z, e2.y, fe2.z, fe2.y);
    AXISTEST_Y1(e2.z, e2.x, fe2.z, fe2.x);
    AXISTEST_Z12(e2.y, e2.x, fe2.y, fe2.x);

    //
    // Bullet 1:
    //

    // first test overlap in the {x,y,z}-directions
    // find min, max of the triangle each direction, and test for overlap in
    // that direction -- this is equivalent to testing a minimal AABB around
    // the triangle against the AABB

    // test in X-direction
    FINDMINMAX(v0.x, v1.x, v2.x, min, max);
    if (min > boxHalfSize.x || max < -boxHalfSize.x)
        return false;

    // test in Y-direction
    FINDMINMAX(v0.y, v1.y, v2.y, min, max);
    if (min > boxHalfSize.y || max < -boxHalfSize.y)
        return false;

    // test in Z-direction
    FINDMINMAX(v0.z, v1.z, v2.z, min, max);
    if (min > boxHalfSize.z || max < -boxHalfSize.z)
        return false;

    //
    // Bullet 2:
    //

    // test if the box intersects the plane of the triangle
    // compute plane equation of triangle: normal*x+d=0
    vec3 normal = cross(e0, e1);
    if (!planeBoxOverlap(normal, v0, boxHalfSize))
        return false;

    return true;
}
