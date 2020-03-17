#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "Sphere.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <bobyqa.h>
#include <fmt/format.h>
#include <nanoflann.hpp>
#include <random>
#include <unordered_set>
#include <vector>

#ifdef WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

struct PointCloud {

    explicit PointCloud(const std::vector<vec3> points)
        : m_points(points)
    {
    }

    // Must return the number of data points
    inline size_t kdtree_get_point_count() const
    {
        return m_points.size();
    }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline float kdtree_get_pt(const size_t idx, const size_t dim) const
    {
        vec3 points = m_points[idx];
        return points[dim];
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template<class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }

private:
    const std::vector<vec3> m_points;
};

using namespace nanoflann;
typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, PointCloud>, PointCloud, 3> MyKdTree;

struct SphereSet {
    std::vector<Sphere> spheres;
    std::vector<float> currentSov;

    std::vector<vec3> points;

    std::unique_ptr<PointCloud> pointCloud;
    std::unique_ptr<MyKdTree> kdTree;

    const VoxelGrid* filledGrid;
    const std::vector<SimpleMesh>* meshes;
};

struct Triangle {
public:
    Triangle(vec3 v0, vec3 v1, vec3 v2)
    {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;

        vec3 crossP = cross(v1 - v0, v2 - v0);
        area = length(crossP) / 2.0f;
        normal = normalize(crossP);
    }

    vec3 v[3];
    vec3 normal;
    float area;
};

void collectSurfacePoints(const std::vector<SimpleMesh>& meshes, std::vector<vec3>& points)
{
    std::vector<Triangle> allTriangles;
    double totalArea = 0.0;

    for (auto& mesh : meshes) {
        size_t triangleCount = mesh.triangleCount();
        for (size_t i = 0; i < triangleCount; ++i) {
            vec3 v0, v1, v2;
            mesh.triangle(i, v0, v1, v2);
            allTriangles.emplace_back(v0, v1, v2);

            totalArea += allTriangles.back().area;
        }
    }

    // Sort largest triangles first
    std::sort(allTriangles.begin(), allTriangles.end(), [](const Triangle& lhs, const Triangle& rhs) {
        return lhs.area > rhs.area;
    });

    // Not an exact science, but we do want some level of coverage
    size_t numSurfacePoints = 1024; //allTriangles.size();
    points.reserve(numSurfacePoints);

    std::mt19937_64 generator;
    for (size_t i = 0; i < numSurfacePoints; ++i) {

        // Select random triangle proportionally to the triangle area
        double randomArea = std::uniform_real_distribution<double> { 0.0, totalArea }(generator);
        double accumArea = 0.0;
        for (auto& triangle : allTriangles) {
            accumArea += triangle.area;
            if (accumArea >= randomArea) {
                // Select random point in triangle
                float b, c;
                do {
                    b = std::uniform_real_distribution<float> { 0.0f, 1.0f }(generator);
                    c = std::uniform_real_distribution<float> { 0.0f, 1.0f }(generator);
                } while (b + c < 1.0f);

                vec3 e1 = triangle.v[1] - triangle.v[0];
                vec3 e2 = triangle.v[2] - triangle.v[0];

                vec3 sample = triangle.v[0] + (e1 * b) + (e2 * c);
                points.push_back(sample);
                break;
            }
        }
    }

    assert(points.size() == numSurfacePoints);
}

void collectInnerPoints(const VoxelGrid& filledGrid, std::vector<vec3>& points)
{
}

SphereSet preprocess(const std::vector<SimpleMesh>& meshes, const VoxelGrid& filledGrid, int numSpheres)
{
    SphereSet set;
    set.filledGrid = &filledGrid;
    set.meshes = &meshes;

    std::vector<vec3> surfacePoints;
    collectSurfacePoints(meshes, surfacePoints);

    // TODO: Needed or not? They describe it in the paper, but I don't know if it actually is needed..
    //std::vector<vec3> innerPoints;
    //collectInnerPoints(filledGrid, innerPoints);

    set.points = surfacePoints;

    set.pointCloud = std::make_unique<PointCloud>(set.points);
    set.kdTree = std::make_unique<MyKdTree>(3 /* dim */, *set.pointCloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    set.kdTree->buildIndex();

    assert(numSpheres > 0);
    assert(numSpheres <= filledGrid.numFilledVoxels());

#if 1
    aabb3 bounds = filledGrid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };
    std::default_random_engine generator {};

    while (set.spheres.size() < numSpheres) {

        float x = uniformInGridX(generator);
        float y = uniformInGridY(generator);
        float z = uniformInGridZ(generator);

        ivec3 gridCoord = filledGrid.remapToGridSpace({ x, y, z }, std::round);
        if (filledGrid.get(gridCoord) > 0) {
            Sphere sphere { vec3(x, y, z), 0.0f };
            set.spheres.push_back(sphere);
            set.currentSov.push_back(0.0f);
        }
    }

#else
    for (const vec3& point : surfacePoints) {
        Sphere sphere { point, 0.01f };
        set.spheres.push_back(sphere);
    }
#endif

    return set;
}

float solidAngleOfTriangleOnSphere(const Sphere& sphere, const vec3& v0, const vec3& v1, const vec3& v2)
{
    // See https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron for information

    // Correct? Yes, I think so..
    vec3 va = v0 - sphere.origin;
    vec3 vb = v1 - sphere.origin;
    vec3 vc = v2 - sphere.origin;

    float a = length(va);
    float b = length(vb);
    float c = length(vc);

    float numerator = abs(dot(va, cross(vb, vc)));
    float denomenator = a * b * c + dot(va, vb) * c + dot(va, vc) * b + dot(vb, vc) * a;

    float tanHalfOmega = numerator / denomenator;
    float halfOmega = atan(tanHalfOmega);
    if (denomenator < 0.0f) {
        halfOmega += mathkit::PI;
    }

    float omega = 2.0f * halfOmega;
    return omega;
}

bool solveQuadratic(const float& a, const float& b, const float& c, float& x0, float& x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else {
        float q = (b > 0)
            ? -0.5f * (b + sqrtf(discr))
            : -0.5f * (b - sqrtf(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1)
        std::swap(x0, x1);

    return true;
}

std::optional<vec3> lineSegmentSphereIntersectionALT(const Sphere& sphere, const vec3& x1, const vec3& x2)
{
    // Consider the line segment as a ray with a max length
    vec3 origin = x1;
    vec3 dir = x2 - x1;
    float maxDist = length(dir);
    dir /= maxDist;

    vec3 L = origin - sphere.origin;
    float a = dot(dir, dir);
    float b = 2.0f * dot(dir, L);
    float c = dot(L, L) - mathkit::square(sphere.radius);

    float t0, t1;
    if (!solveQuadratic(a, b, c, t0, t1)) {
        return {};
    }

    if (t0 > t1) {
        std::swap(t0, t1);
    }

    if (t0 < 0.0f) {
        t0 = t1; // if t0 is negative, let's use t1 instead
        if (t0 < 0) {
            // both t0 and t1 are negative
            return {};
        }
    }

    float t = t0;
    return origin + t * dir;
}

std::optional<vec3> lineSegmentSphereIntersection(const Sphere& sphere, const vec3& a, const vec3& b)
{
    // Consider the line segment as a ray with a max length
    vec3 o = a;
    vec3 l = b - a;
    float maxDist = length(l);
    l /= maxDist;

    const vec3& c = sphere.origin;
    const float& r = sphere.radius;

    float inner = mathkit::square(dot(l, o - r)) - (length2(o - c) - r * r);
    if (inner < 0.0f) {
        return {};
    }

    float t1 = sqrt(inner);
    float t0 = -dot(l, o - c);

    float d0 = t0 - t1;
    float d1 = t0 + t1;

    bool d0valid = d0 > 0.0f && d0 < maxDist;
    bool d1valid = d1 > 0.0f && d1 < maxDist;
    assert(!(d0valid && d1valid));

    if (d0valid)
        return o + d0 * l;
    if (d1valid)
        return o + d1 * l;
    else
        return {};
}

auto sphereTriangleIntersection(const Sphere& sphere, const Triangle& triangle)
{
    auto a = lineSegmentSphereIntersectionALT(sphere, triangle.v[0], triangle.v[1]);
    auto b = lineSegmentSphereIntersectionALT(sphere, triangle.v[0], triangle.v[2]);
    auto c = lineSegmentSphereIntersectionALT(sphere, triangle.v[1], triangle.v[2]);

    struct {
        vec3 p0;
        vec3 p1;
    } returnVal;

    if (a.has_value()) {
        returnVal.p0 = a.value();
        if (b.has_value())
            returnVal.p1 = b.value();
        else if (c.has_value())
            returnVal.p1 = c.value();
        else
            assert(false);
    } else {
        returnVal.p0 = b.value();
        returnVal.p1 = c.value();
    }

    return returnVal;
}

float sotvCaseA(const Sphere& sphere, const Triangle& triangle, float h)
{
    const float& r = sphere.radius;
    float omega = solidAngleOfTriangleOnSphere(sphere, triangle.v[0], triangle.v[1], triangle.v[2]);
    return (1.0f / 3.0f) * ((r * r * r * omega) - (triangle.area * h));
}

float sotvCaseB(const Sphere& sphere, const Triangle& triangle, float h)
{
    return mathkit::PI * (h * h) * (sphere.radius - (h / 3.0f));
}

float sotvCaseC(const Sphere& sphere, const Triangle& triangle, const bool inSphere[3], float h)
{
    auto [p0, p1] = sphereTriangleIntersection(sphere, triangle);
    vec3 p2 = inSphere[0]
        ? triangle.v[0]
        : inSphere[1]
            ? triangle.v[1]
            : inSphere[2]
                ? triangle.v[2]
                : vec3(0.0f);

    Triangle subTriangle { p0, p1, p2 };
    float caseA = sotvCaseA(sphere, subTriangle, h);

    // TODO: Calulate swing volume!
    float swing = 0.0f;

    return caseA + swing;
}

auto findInsideVertices(const Triangle& triangle, const bool inSphere[3])
{
    struct {
        vec3 a, b;
    } returnVal;

    if (inSphere[0]) {
        returnVal.a = triangle.v[0];
        if (inSphere[1])
            returnVal.b = triangle.v[1];
        else if (inSphere[2])
            returnVal.b = triangle.v[2];
        else
            assert(false);
    } else {
        returnVal.a = triangle.v[1];
        returnVal.b = triangle.v[2];
    }

    return returnVal;
}

float sotvCaseD(const Sphere& sphere, const Triangle& triangle, const bool inSphere[3], float h)
{
    auto [p0, p1] = sphereTriangleIntersection(sphere, triangle);
    auto [p2, p3] = findInsideVertices(triangle, inSphere);

    Triangle caseAtri { p0, p2, p3 };
    float caseA = sotvCaseA(sphere, caseAtri, h);

    Triangle caseCtri { p0, p3, p1 }; // (technically not outside sphere, but they will have a swing volume)
    const bool inSphereCaseC[3] = { false, true, false };
    float caseC = sotvCaseC(sphere, caseCtri, inSphereCaseC, h);

    return caseA + caseC;
}

float sphereOutsideTriangleVolume(const Sphere& sphere, const Triangle& triangle)
{
    bool inSphere[3];
    int numVerticesOutsideSphere = 3;
    for (int i = 0; i < 3; ++i) {
        float dist = distance(sphere.origin, triangle.v[i]);
        inSphere[i] = dist < sphere.radius;
        numVerticesOutsideSphere -= inSphere[i];
    }

    // TODO: Correct?! H is the distance between the sphere center and the triangle's plane
    float h = abs(dot(triangle.normal, sphere.origin - triangle.v[0]));

    // Is the triangle fully outside the sphere?
    // (actually, I'm not sure we should special case this..)
    //if (numVerticesOutsideSphere == 3 && h > sphere.radius) {
    //    return 0.0f;
    //}

    switch (numVerticesOutsideSphere) {
    case 0:
        return sotvCaseA(sphere, triangle, h);
    case 3:
        return sotvCaseB(sphere, triangle, h);
    case 2:
        return sotvCaseC(sphere, triangle, inSphere, h);
    case 1:
        return sotvCaseD(sphere, triangle, inSphere, h);
    default:
        assert(false);
    }
}

#define SIGN(x) (x < 0 ? -1 : +1)

float sphereOutsideVolumeAnalytical(const Sphere& sphere, const SphereSet& set)
{
    // SOV is the volume outside the triangle set, but inside the sphere
    //float sov = 0.0f;
    float vIn = 0.0f;

    for (auto& mesh : *set.meshes) {
        size_t triangleCount = mesh.triangleCount();
        for (size_t i = 0; i < triangleCount; ++i) {
            vec3 v0, v1, v2;
            mesh.triangle(i, v0, v1, v2);
            Triangle triangle { v0, v1, v2 };

            float h = dot(triangle.normal, sphere.origin - v0);
            if (abs(h) > sphere.radius) {
                // No intersection between the sphere and the triangle
                continue;
            }

            float sotv = sphereOutsideTriangleVolume(sphere, triangle);

            vec3 p = sphere.origin - h * triangle.normal;
            vIn += SIGN(dot(triangle.normal, p - sphere.origin)) * sotv;
        }
    }

    // TODO: Figure out if the sphere center is inside the object!
    //  Could we not just make sure it always is? Maybe not optimal for some cases..?
    bool centerInsideObject = true;

    float sov = centerInsideObject
        ? vIn
        : sphere.volume() + vIn;
    return sov;
}

float totalSphereOutsideVolumeAnalytical(const SphereSet& set)
{
    float totalSOV = 0.0f;
    for (auto& sphere : set.spheres) {
        totalSOV += sphereOutsideVolumeAnalytical(sphere, set);
    }
    return totalSOV;
}

/*
float sphereOutsideVolumeUsingVoxels(const Sphere& sphere, const VoxelGrid& grid)
{
    // TODO: Is there enough precision when using voxels?

    if (sphere.radius < 1e-6f) {
        return 0.0f;
    }

    // TODO: It's (very?) important that we don't have holes in the filled grid now..
    int numOutsideVoxels = 0;

    vec3 voxelSize = grid.voxelSize();
    vec3 first = sphere.origin - vec3(sphere.radius) - (0.05f * voxelSize);
    vec3 last = sphere.origin + vec3(sphere.radius) + (0.05f * voxelSize);

    first = max(first, grid.gridBounds().min);
    last = min(last, grid.gridBounds().max);

    for (float z = first.z; z <= last.z; z += voxelSize.z) { // NOLINT(cert-flp30-c)
        for (float y = first.y; y <= last.y; y += voxelSize.y) { // NOLINT(cert-flp30-c)
            for (float x = first.x; x <= last.x; x += voxelSize.x) { // NOLINT(cert-flp30-c)

                vec3 samplePoint = vec3(x, y, z);

                // Only count voxels that actually are inside the sphere
                if (distance2(sphere.origin, samplePoint) > sphere.radius * sphere.radius) {
                    continue;
                }

                ivec3 gridSampleCoord = grid.remapToGridSpace(samplePoint, std::round);

                // NOTE: Should be handled by the coord clamping above this loop
                //if (any(lessThan(gridSampleCoord, ivec3(0))) || any(greaterThanEqual(gridSampleCoord, grid.gridDimensions()))) {
                //    continue;
                //}

                // If this voxel (which is considered to be inside the sphere) is empty, we assume that the voxel
                // is actually outside the mesh (since the grid should be filled), and so we count this voxel.
                uint32_t value = grid.get(gridSampleCoord);
                if (value == 0) {
                    numOutsideVoxels += 1;
                }
            }
        }
    }

    float sov = float(numOutsideVoxels) * (voxelSize.x * voxelSize.y * voxelSize.z);
    return sov;
}

float totalSphereOutsideVolumeUsingVoxels(const SphereSet& set, const VoxelGrid& grid)
{
    float totalSOV = 0.0f;
    for (auto& sphere : set.spheres) {
        totalSOV += sphereOutsideVolumeUsingVoxels(sphere, grid);
    }
    return totalSOV;
}
*/

float pointAssignmentFromListOfPoints(SphereSet& set, float previousError)
{
    // TODO: From the paper "the [current] (naive) algorithm for point assignment can be accelerated"
    //  See p. 5 in the paper for more information. It's already slow and the SOV function is just a stub!

    std::vector<bool> visited;
    visited.resize(set.points.size());

    struct PointRef {
        PointRef(vec3 point, int32_t index)
            : point(point)
            , index(index)
        {
        }
        vec3 point;
        int32_t index;
    };

    std::vector<PointRef> stack {};

    for (auto& [sphereCenter, _] : set.spheres) {
        PointRef pointRef { sphereCenter, -1 };
        stack.push_back(pointRef);
    }

    auto oldSpheres = set.spheres;
    auto oldSphereSov = set.currentSov;

    int numAssigned = 0;
    int numToAssign = set.points.size() + set.spheres.size();

    while (!stack.empty()) {
        PointRef pointRef = stack.back();
        stack.pop_back();

        int bestSphereIdx = -1;
        float smallestSovIncrease = INFINITY;

        float newPotentialRadius = INFINITY;
        float newPotentialSov = INFINITY;

        for (int si = 0; si < set.spheres.size(); ++si) {
            Sphere& sphere = set.spheres[si];

            // Create a test sphere with radius required to contain the test point
            Sphere testSphere = sphere;
            testSphere.radius = distance(sphere.origin, pointRef.point);

            // Check the increase in SOV for the test sphere
            float testSov = sphereOutsideVolumeAnalytical(testSphere, set);
            float sovIncrease = testSov - set.currentSov[si];

            if (sovIncrease < smallestSovIncrease) {
                smallestSovIncrease = sovIncrease;
                bestSphereIdx = si;

                // Save these for speed (radius) and accuracy (SOV)
                newPotentialRadius = testSphere.radius;
                newPotentialSov = testSov;
            }
        }

        // Expand the best sphere to contain the new point
        assert(bestSphereIdx != -1);
        set.spheres[bestSphereIdx].radius = newPotentialRadius;
        set.currentSov[bestSphereIdx] = newPotentialSov;

        // Insert unvisited neighbors into the stack
#if 1
        size_t maxNeighboursToAdd = 10;

    stupidGotoForEasyExperimenting:

        std::vector<size_t> retIndex(maxNeighboursToAdd);
        std::vector<float> outDistSqr(maxNeighboursToAdd);

        size_t numActualNeighbours = set.kdTree->knnSearch(value_ptr(pointRef.point), maxNeighboursToAdd, &retIndex[0], &outDistSqr[0]);

        // In case of less points in the tree than requested:
        retIndex.resize(numActualNeighbours);
        outDistSqr.resize(numActualNeighbours);

        int numAdded = 0;
        for (size_t ki = 0; ki < numActualNeighbours; ++ki) {

            size_t pointIdx = retIndex[ki];
            if (!visited[pointIdx]) {
                visited[pointIdx] = true;

                vec3 point = set.points[pointIdx];
                stack.emplace_back(point, pointIdx);

                numAdded += 1;

                assert(stack.size() <= numToAssign);
            }
        }

        // If we can't find any more points to add and we have considered all possible points
        if (numAdded == 0) {
            if (maxNeighboursToAdd < set.points.size()) {
                maxNeighboursToAdd += 10;
                goto stupidGotoForEasyExperimenting;
            }
        }
#else
        float nearestUnvisitedDist2 = 99999.99f;
        int64_t nearestUnvisitedIdx = -1;

        for (size_t pi = 0; pi < set.points.size(); ++pi) {

            if (visited[pi]) {
                continue;
            }

            const vec3& neighbour = set.points[pi];
            float dist2 = distance2(pointRef.point, neighbour);

            if (dist2 < nearestUnvisitedDist2) {
                nearestUnvisitedDist2 = dist2;
                nearestUnvisitedIdx = pi;
            }
        }

        //assert(nearestUnvisitedIdx != -1);
        if (nearestUnvisitedIdx != -1) {
            visited[nearestUnvisitedIdx] = true;
            stack.emplace_back(set.points[nearestUnvisitedIdx], nearestUnvisitedIdx);
        }
#endif

        numAssigned += 1;

        // Hmm, now we have set.points.size() which isn't exactly the number of assigned stuff, right?
        if (numAssigned % 30 == 0 || numAssigned == numToAssign) {
            fmt::print("     assigned {}/{} points ({:.1f}%)\r", numAssigned, numToAssign, 100.0f * float(numAssigned) / float(numToAssign));
            fflush(stdout);
        }
    }
    fmt::print("\n");

    float newError = totalSphereOutsideVolumeAnalytical(set);

    if (newError < previousError) {
        // Success, continue
        return newError;
    } else {
        // Failure, rollback the changes & trigger sphere teleportation
        set.spheres = oldSpheres;
        set.currentSov = oldSphereSov;
        return previousError;
    }
}

float pointAssignment(SphereSet& set, float previousError)
{
    // TODO: From the paper "the [current] (naive) algorithm for point assignment can be accelerated"
    //  See p. 5 in the paper for more information. It's already slow and the SOV function is just a stub!

    const VoxelGrid& grid = *set.filledGrid;
    size_t numFilledVoxels = grid.numFilledVoxels();

    struct PointRef {
        PointRef(vec3 p, ivec3 c)
            : point(p)
            , gridIndex(c)
        {
        }
        vec3 point;
        ivec3 gridIndex;
    };

    std::unordered_set<uint64_t> visitedPoints {};
    visitedPoints.reserve(numFilledVoxels);

    std::vector<PointRef> stack {};

    for (auto& [sphereCenter, _] : set.spheres) {
        ivec3 gridIndex = grid.remapToGridSpace(sphereCenter, std::round);
        PointRef pointRef { sphereCenter, gridIndex };
        stack.push_back(pointRef);

        uint64_t linearIndex = grid.linearIndex(gridIndex);
        visitedPoints.insert(linearIndex);
    }

    auto oldSpheres = set.spheres;
    auto oldSphereSov = set.currentSov;

    int numAssigned = 0;

    // TODO: OpenMP would be nice here! Just make sure to mark some critical sections! However, it seems like the rounding of
    //  coordinates (std::roundf) is the main contributor of self-time according to the profiler. Worth looking into omp anyway!
    while (!stack.empty()) {
        PointRef pointRef = stack.back();
        stack.pop_back();

        int bestSphereIdx = -1;
        float smallestSovIncrease = INFINITY;

        float newPotentialRadius = INFINITY;
        float newPotentialSov = INFINITY;

        for (int si = 0; si < set.spheres.size(); ++si) {
            Sphere& sphere = set.spheres[si];

            // Create a test sphere with radius required to contain the test point
            Sphere testSphere = sphere;
            testSphere.radius = distance(sphere.origin, pointRef.point);

            // Check the increase in SOV for the test sphere
            float testSov = sphereOutsideVolumeAnalytical(testSphere, set);
            float sovIncrease = testSov - set.currentSov[si];

            if (sovIncrease < smallestSovIncrease) {
                smallestSovIncrease = sovIncrease;
                bestSphereIdx = si;

                // Save these for speed (radius) and accuracy (SOV)
                newPotentialRadius = testSphere.radius;
                newPotentialSov = testSov;
            }
        }

        // Expand the best sphere to contain the new point
        assert(bestSphereIdx != -1);
        set.spheres[bestSphereIdx].radius = newPotentialRadius;
        set.currentSov[bestSphereIdx] = newPotentialSov;

        // Insert unvisited neighbors into the stack
        const auto& immediateNeighbors = grid.immediateFilledNeighbors(pointRef.gridIndex);
        for (ivec3 gridCoord : immediateNeighbors) {

            uint64_t linearIndexOfNeighbor = grid.linearIndex(gridCoord);
            if (visitedPoints.find(linearIndexOfNeighbor) != visitedPoints.end()) {
                continue;
            }

            vec3 voxelCenter = grid.voxelCenterPoint(gridCoord);
            stack.emplace_back(voxelCenter, gridCoord);

            uint64_t linearIndex = grid.linearIndex(gridCoord);
            visitedPoints.insert(linearIndex);

            assert(stack.size() <= numFilledVoxels);
        }

        numAssigned += 1;

        if (numAssigned % 100 == 0 || numAssigned == numFilledVoxels) {
            fmt::print("     assigned {}/{} points ({:.1f}%)\r", numAssigned, numFilledVoxels, 100.0f * float(numAssigned) / float(numFilledVoxels));
            fflush(stdout);
        }
    }
    fmt::print("\n");

    float newError = totalSphereOutsideVolumeAnalytical(set);

    if (newError < previousError) {
        // Success, continue
        return newError;
    } else {
        // Failure, rollback the changes & trigger sphere teleportation
        set.spheres = oldSpheres;
        set.currentSov = oldSphereSov;
        return previousError;
    }
}

REAL sphereFittingObjectiveFunction(const INTEGER n, const REAL* x, void* data)
{
    assert(n == 4);

    Sphere sphere;
    sphere.radius = x[3];
    sphere.origin = { x[0], x[1], x[2] };

    const SphereSet& sphereSet = *reinterpret_cast<SphereSet*>(data);
    float sov = sphereOutsideVolumeAnalytical(sphere, sphereSet);

    return sov;
}

void sphereFitting(SphereSet& set)
{
    for (auto& sphere : set.spheres) {

        if (sphere.radius <= 0.0f) {
            //fmt::print("Ignoring sphere (which we probably shouldn't..?\n");
            continue;
        }

        // Origin (3) + radius (1)
        constexpr INTEGER n = 3 + 1;

        // NPT is the number of interpolation conditions.  Its value must be in the interval
        // [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.
        constexpr INTEGER npt = n + 2; // TODO: Find a good value!

        REAL x[n] = {
            sphere.origin.x,
            sphere.origin.y,
            sphere.origin.z,
            sphere.radius
        };

        aabb3 gridBounds = set.filledGrid->gridBounds();
        double gridMaxDistance = length(gridBounds.max - gridBounds.min);

        // TODO: Find good values!
        REAL xLower[n] = {
            gridBounds.min.x,
            gridBounds.min.y,
            gridBounds.min.z,
            0.0
        };

        // TODO: Find good values!
        REAL xUpper[n] = {
            gridBounds.max.x,
            gridBounds.max.y,
            gridBounds.max.z,
            gridMaxDistance
        };

        // "Typically, RHOBEG should be about one tenth of the greatest expected change to a variable"
        // TODO: Find good value!
        REAL rhoBeg = gridMaxDistance / 10.0;

        // RHOEND should indicate the accuracy that is required in the final values of the variables
        // NOTE: Some models are tiny, others are huge.. So an absolute/constant accuracy is probably
        //  not very helpful here. A fraction of the sphere radius could maybe work as a scale-aware metric?
        // TODO: Find good value!
        REAL rhoEnd = gridMaxDistance * 1e-3;

        //The array W will be used for working space.  Its length must be at least
        // (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
        // will be set to the function value at the solution. */
        REAL workingMemory[(npt + 5) * (npt + n) + 3 * n * (n + 5) / 2];

        INTEGER logLevel = 1;
        INTEGER maxObjFunCalls = 100;

        int status = bobyqa(n, npt, sphereFittingObjectiveFunction, (void*)(&set),
                            x, xLower, xUpper, rhoBeg, rhoEnd,
                            logLevel, maxObjFunCalls, workingMemory);

        switch (status) {
        case BOBYQA_SUCCESS:
            // algorithm converged!
            sphere.origin = { x[0], x[1], x[2] };
            sphere.radius = x[3];
            return;
        case BOBYQA_BAD_NPT:
            fmt::print("bobyqa error: NPT is not in the required interval\n");
            assert(false); // npt is hardcoded to a valid value!
            break;
        case BOBYQA_TOO_CLOSE:
            fmt::print("bobyqa error: insufficient space between the bounds\n");
            break;
        case BOBYQA_ROUNDING_ERRORS:
            fmt::print("bobyqa error: too much cancellation in a denominator\n");
            break;
        case BOBYQA_TOO_MANY_EVALUATIONS:
            fmt::print("bobyqa error: maximum number of function evaluations exceeded\n");
            break;
        case BOBYQA_STEP_FAILED:
            fmt::print("bobyqa error: a trust region step has failed to reduce Q\n");
            break;
        default:
            assert(false);
        }
    }
}

float sphereOverlap(const Sphere& a, const Sphere& b)
{
    // See https://math.stackexchange.com/questions/297751/overlapping-spheres

    float d = distance(a.origin, b.origin);
    if (d >= a.radius + b.radius) {
        return 0.0f;
    }

    const float& r1 = a.radius;
    const float& r2 = b.radius;

    using namespace mathkit;
    return PI / (12.0f * d) * square(r1 + r2 - d) * (square(d) + 2.0f * d * (r1 + r2) - 3.0f * square(r1 - r2));
}

void sphereTeleportation(SphereSet& set)
{
    //
    // Find the sphere which is most overlapped by the other spheres
    //  & the sphere with highest SOV value
    //

    float totalVolume = 0.0f;
    for (auto& sphere : set.spheres) {
        totalVolume += sphere.volume();
    }

    float maxOverlapRatio = -INFINITY;
    int maxOverlapRatioIdx = -1;

    float maxSov = -INFINITY;
    int maxSovIdx = -1;

    for (size_t i = 0; i < set.spheres.size(); ++i) {

        Sphere& self = set.spheres[i];

        float sharedVolume = 0.0f;

        for (size_t k = 0; k < set.spheres.size(); ++k) {
            if (i == k) {
                continue;
            }

            Sphere& other = set.spheres[k];
            sharedVolume += sphereOverlap(self, other);
        }

        float overlapRatio = sharedVolume / (self.volume() + sharedVolume);
        if (overlapRatio > maxOverlapRatio) {
            maxOverlapRatio = overlapRatio;
            maxOverlapRatioIdx = i;
        }

        float sov = sphereOutsideVolumeAnalytical(self, set);
        if (sov > maxSov) {
            maxSov = sov;
            maxSovIdx = i;
        }
    }

    assert(maxSovIdx != -1);
    assert(maxOverlapRatioIdx != -1);

    //
    // Split the highest SOV sphere into two, with new sphere centers being
    // the two points furthest away from each other in the initial sphere,
    // and remove the most overlapped sphere
    //

    Sphere& sphereToSplit = set.spheres[maxSovIdx];
    Sphere& sphereToRemove = set.spheres[maxOverlapRatioIdx];

    const auto& splitCandidates = set.filledGrid->filledVoxelsInSphere(sphereToSplit);
    assert(!splitCandidates.empty());

    if (splitCandidates.size() == 1) {
        // Only one voxel, so pick two two corners of the voxel
        sphereToSplit.origin = splitCandidates[0] - 0.5f * set.filledGrid->voxelSize();
        sphereToRemove.origin = splitCandidates[0] + 0.5f * set.filledGrid->voxelSize();
    } else {
        float maxDistance2 = 0.0f;
        for (size_t i = 0, len = splitCandidates.size(); i < len; ++i) {
            const vec3& a = splitCandidates[i];
            for (size_t k = i + 1; k < len; ++k) {
                const vec3& b = splitCandidates[k];
                float dist2 = distance2(a, b);
                if (dist2 > maxDistance2) {
                    maxDistance2 = dist2;
                    sphereToSplit.origin = a;
                    sphereToRemove.origin = b;
                }
            }
        }

        assert(maxDistance2 > 1e-6f);
    }
}

void generateOutput(const SphereSet& sphereSet)
{
    size_t sphereCount = sphereSet.spheres.size();

    fmt::print("const vec4 spheres[] = vec4[](\n");
    for (size_t i = 0; i < sphereCount; ++i) {
        const Sphere& sphere = sphereSet.spheres[i];
        fmt::print("\tvec4({}, {}, {}, {}){}\n",
                   sphere.origin.x, sphere.origin.y, sphere.origin.z,
                   sphere.radius, (i + 1 == sphereCount) ? "" : ",");
    }
    fmt::print(");\n");
}

void setApplicationWorkingDirectory(char* executableName, const std::string& workingDir)
{
#ifdef WIN32
    char fullPathBuf[_MAX_PATH] = {};
    assert(_fullpath(fullPathBuf, executableName, sizeof(fullPathBuf)));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(_chdir(newWorkingDir.c_str()) == 0);
#else
    char fullPathBuf[PATH_MAX] = {};
    assert(realpath(executableName, fullPathBuf));
    std::string fullPath { fullPathBuf };

    size_t startOfWorkingDirName = fullPath.find(workingDir);
    std::string newWorkingDir = fullPath.substr(0, startOfWorkingDirName + workingDir.length() + 1);
    assert(chdir(newWorkingDir.c_str()) == 0);
#endif
}

int main(int argc, char** argv)
{
    setApplicationWorkingDirectory(argv[0], "ProxyGen");

    std::string path = "assets/Bunny/bunny_lowres.gltf";
    //std::string path = "assets/Avocado/Avocado.gltf";
    //std::string path = "assets/Cube/Cube.gltf";
    //std::string path = "assets/BoomBox/BoomBoxWithAxes.gltf";

    constexpr int numSpheres = 60;
    const size_t gridDimensions = 64;

    std::vector<SimpleMesh> simpleMeshes {};

    fmt::print("=> loading model '{}'\n", path);

    auto [basePath, model] = GltfUtil::loadModel(path);
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, simpleMeshes);

    fmt::print("=> shell voxelization\n");

    aabb3 boundsOfAllMeshes = SimpleMesh::calculateBounds(simpleMeshes);
    VoxelGrid shellGrid { glm::ivec3(gridDimensions), boundsOfAllMeshes };

    for (size_t i = 0; i < simpleMeshes.size(); ++i) {
        fmt::print("     mesh {}/{}\n", i + 1, simpleMeshes.size());
        shellGrid.insertMesh(simpleMeshes[i], true);
    }

    fmt::print("=> color quantization\n");

    shellGrid.quantizeColors(10);

    fmt::print("=> volume filling\n");

    VoxelGrid filledGrid = shellGrid;
    filledGrid.fillVolumes(simpleMeshes);

    fmt::print("=> preprocessing sphere set\n");
    SphereSet sphereSet = preprocess(simpleMeshes, filledGrid, numSpheres);

    fmt::print("=> optimizing begin\n");

    bool didJustTeleport = false;
    float previousError = INFINITY;

    constexpr int maxIterations = 100;
    int iteration = 0;
    for (; iteration < maxIterations; ++iteration) {

        fmt::print("==> point assignment\n");
        //float newError = pointAssignment(sphereSet, previousError);
        float newError = pointAssignmentFromListOfPoints(sphereSet, previousError);
        fmt::print("     error: {}\n", newError);

        generateOutput(sphereSet);

        if (newError < 1e-12f) {
            // in practice this probably shouldn't happen, but since our voxel grid resolution is quite low it can happen.
            fmt::print("==> error is zero\n");
            break;
        }

        if (didJustTeleport) {
            didJustTeleport = false;
            if (newError >= previousError) {
                // "if [teleportation] fails to decrease error, the algorithm terminates"
                fmt::print("==> teleportation failed to decrease error\n");
                break;
            }
        }

        if (newError < previousError) {
            fmt::print("===> sphere fitting");
            sphereFitting(sphereSet);
        } else {
            fmt::print("===> sphere teleportation\n");
            sphereTeleportation(sphereSet);
            didJustTeleport = true;
        }

        previousError = newError;
    }

    if (iteration == maxIterations) {
        fmt::print("==> max iterations reached\n");
    }

    fmt::print("=> optimizing done\n");

    fmt::print("=> generating sphere output\n");
    generateOutput(sphereSet);
}
