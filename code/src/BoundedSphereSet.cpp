#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "Sphere.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <Eigen/Core>
#include <bobyqa.h>
#include <chrono>
#include <fmt/format.h>
#include <fstream>
#include <iomanip>
#include <json.hpp>
#include <random>
#include <spherical_harmonics.h>
#include <unordered_set>
#include <vector>

#ifdef WIN32
#define NOMINMAX
#include <direct.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

struct SphericalHarmonics {
    vec3 L00 {};
    vec3 L11 {};
    vec3 L10 {};
    vec3 L1_1 {};
    vec3 L21 {};
    vec3 L2_1 {};
    vec3 L2_2 {};
    vec3 L20 {};
    vec3 L22 {};
};

struct SphereSet {
    SphereSet(const SimpleMesh& mesh)
        : mesh(mesh)
    {
    }

    std::vector<Sphere> spheres;
    std::vector<SphericalHarmonics> sphereSH;

    const SimpleMesh& mesh;

    struct {
        std::unique_ptr<VoxelGrid> shell;
        std::unique_ptr<VoxelGrid> filled;
        std::unique_ptr<VoxelGrid> inside;
    } grids;
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

vec3 getRandomSphereLocationInsideMesh(const VoxelGrid& grid, std::default_random_engine& generator)
{
    aabb3 bounds = grid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };

    vec3 location;
    ivec3 gridCoord;

    do {
        location = { uniformInGridX(generator),
                     uniformInGridY(generator),
                     uniformInGridZ(generator) };
        gridCoord = grid.remapToGridSpace(location, std::round);
    } while (grid.get(gridCoord) == 0);

    return location;
}

void assignSpheresRandomlyInVolume(SphereSet& set, unsigned numSpheres)
{
    auto& grid = *set.grids.inside;
    assert(numSpheres <= grid.numFilledVoxels());

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    for (size_t i = 0; i < numSpheres; ++i) {
        vec3 location = getRandomSphereLocationInsideMesh(grid, generator);
        Sphere sphere { location, 0.0 };
        set.spheres.push_back(sphere);
    }

    aabb3 bounds = grid.gridBounds();
    std::uniform_real_distribution<float> uniformInGridX { bounds.min.x, bounds.max.x };
    std::uniform_real_distribution<float> uniformInGridY { bounds.min.y, bounds.max.y };
    std::uniform_real_distribution<float> uniformInGridZ { bounds.min.z, bounds.max.z };

    while (set.spheres.size() < numSpheres) {

        float x = uniformInGridX(generator);
        float y = uniformInGridY(generator);
        float z = uniformInGridZ(generator);

        ivec3 gridCoord = grid.remapToGridSpace({ x, y, z }, std::round);
        if (grid.get(gridCoord) > 0) {
            Sphere sphere { vec3(x, y, z), 0.0f };
            set.spheres.push_back(sphere);
        }
    }
}

bool sphereTriangleIntersection(const Sphere& sphere, const vec3& v0, const vec3& v1, const vec3& v2)
{
    // From https://realtimecollisiondetection.net/blog/?p=103

    vec3 P = sphere.origin;
    float r = sphere.radius;
    float rr = r * r;

    vec3 A = v0 - P;
    vec3 B = v1 - P;
    vec3 C = v2 - P;

    vec3 V = cross(B - A, C - A);
    float d = dot(A, V);
    float e = dot(V, V);
    bool sep1 = d * d > rr * e;

    float aa = dot(A, A);
    float ab = dot(A, B);
    float ac = dot(A, C);
    float bb = dot(B, B);
    float bc = dot(B, C);
    float cc = dot(C, C);

    bool sep2 = (aa > rr) & (ab > aa) & (ac > aa);
    bool sep3 = (bb > rr) & (ab > bb) & (bc > bb);
    bool sep4 = (cc > rr) & (ac > cc) & (bc > cc);

    vec3 AB = B - A;
    vec3 BC = C - B;
    vec3 CA = A - C;

    float d1 = ab - aa;
    float d2 = bc - bb;
    float d3 = ac - cc;

    float e1 = dot(AB, AB);
    float e2 = dot(BC, BC);
    float e3 = dot(CA, CA);

    vec3 Q1 = A * e1 - d1 * AB;
    vec3 Q2 = B * e2 - d2 * BC;
    vec3 Q3 = C * e3 - d3 * CA;

    vec3 QC = C * e1 - Q1;
    vec3 QA = A * e2 - Q2;
    vec3 QB = B * e3 - Q3;

    bool sep5 = (dot(Q1, Q1) > rr * e1 * e1) && (dot(Q1, QC) > 0);
    bool sep6 = (dot(Q2, Q2) > rr * e2 * e2) && (dot(Q2, QA) > 0);
    bool sep7 = (dot(Q3, Q3) > rr * e3 * e3) && (dot(Q3, QB) > 0);

    bool separated = sep1 | sep2 | sep3 | sep4 | sep5 | sep6 | sep7;
    return !separated;
}

bool sphereInsideMeshVolume(const Sphere& sphere, const SimpleMesh& mesh)
{
    // NOTE: Implicit, but this assumes that the mesh is one solid blob with no islands!

    // NOTE: To avoid expensive(/not-implemented) checks to see if the sphere centers are inside the volume
    //  of the object, we will assume they are and through other means verify that it will be invariant. The
    //  other means are: 1) spawning points inside the volume, and 2) always making steps small enough that
    //  we wont leave the volume when optimizing.

    if (sphere.radius == 0.0) {
        return true;
    }

    volatile bool result = true;
    #pragma omp parallel for shared(result)
    for (int i = 0; i < mesh.triangleCount(); ++i) {
        if (result) {
            // TODO: Maybe don't reconstruct these every query?!
            vec3 v0, v1, v2;
            mesh.triangle(i, v0, v1, v2);

            if (sphereTriangleIntersection(sphere, v0, v1, v2)) {
                result = false;
            }
        }
    }

    return result;
}

bool sphereSetInsideMeshVolume(const SphereSet& set, const SimpleMesh& mesh)
{
    for (const Sphere& sphere : set.spheres) {
        if (!sphereInsideMeshVolume(sphere, mesh)) {
            return false;
        }
    }
    return true;
}

void expandSpheresMaximally(SphereSet& set)
{
    aabb3 gridBounds = set.grids.shell->gridBounds();
    double gridDiameter = distance(gridBounds.min, gridBounds.max);

    #pragma omp parallel for
    for (int i = 0; i < set.spheres.size(); ++i) {
        Sphere& sphere = set.spheres[i];

        Sphere minSphere { sphere.origin, 0.0 };
        Sphere maxSphere { sphere.origin, gridDiameter };

        assert(sphereInsideMeshVolume(minSphere, set.mesh));
        assert(!sphereInsideMeshVolume(maxSphere, set.mesh));

        while (abs(maxSphere.radius - minSphere.radius) > 1e-8) {

            double midpoint = (minSphere.radius + maxSphere.radius) / 2.0;
            Sphere middleSphere { sphere.origin, midpoint };

            bool midOk = sphereInsideMeshVolume(middleSphere, set.mesh);
            if (midOk) {
                minSphere = middleSphere;
            } else {
                maxSphere = middleSphere;
            }
        }

        // Make sure to be conservative with the choice to we don't break the invariant
        sphere.radius = minSphere.radius;
        assert(sphereInsideMeshVolume(sphere, set.mesh));
    }

#if 0
    std::sort(set.spheres.begin(), set.spheres.end(), [](const Sphere& lhs, const Sphere& rhs) {
        return lhs.radius > rhs.radius;
    });
    fmt::print("   grid 'radius' {}\n", gridDiameter / 2.0);
    for (auto& sphere : set.spheres) {
        fmt::print("     sphere radius {}\n", sphere.radius);
    }
#endif
}

double sphereSetIntersectionVolumeOverestimation(const SphereSet& set)
{
    // Except some special cases there are no analytical solutions for volume intersection of many spheres. See
    // https://www.researchgate.net/publication/289097700_Volume_of_intersection_of_six_spheres_A_special_case_of_practical_interest
    // for some more information. In this case we don't really care about the exact intersection, but we do want to avoid as much
    // intersection as possible, so just overcounting it works here, I think..

    double intersection = 0.0;

    size_t sphereCount = set.spheres.size();
    for (size_t i = 0; i < sphereCount; ++i) {
        const Sphere& self = set.spheres[i];
        for (size_t k = i + 1; k < sphereCount; ++k) {
            const Sphere& other = set.spheres[k];
            intersection += self.overlapWith(other);
        }
    }

    return intersection;
}

double sphereSetUnionVolume(const SphereSet& set)
{
    double volume = 0.0;

    size_t sphereCount = set.spheres.size();
    for (const Sphere& sphere : set.spheres) {
        volume += sphere.volume();
    }

    // TODO: Maybe scale down the intersection a bit, since we know it overestimates?
    double intersectionOverestimation = sphereSetIntersectionVolumeOverestimation(set);
    volume -= 0.25 * intersectionOverestimation;

    return volume;
}

struct SphereFittingData {
    SphereSet* set;
    unsigned sphereIndex;
};

REAL sphereFittingObjectiveFunction(const INTEGER n, const REAL* x, void* dataPtr)
{
    assert(n == 4);

    Sphere testSphere;
    testSphere.radius = x[3];
    testSphere.origin = { x[0], x[1], x[2] };

    SphereFittingData& data = *reinterpret_cast<SphereFittingData*>(dataPtr);
    SphereSet& set = *data.set;

    // Temporarily swap in the sphere to test
    Sphere currentSphere = set.spheres[data.sphereIndex];
    set.spheres[data.sphereIndex] = testSphere;

    if (!sphereInsideMeshVolume(testSphere, set.mesh)) {
        return 99999.99; // BOBYQA doesn't seem to like infinities, but this will do
    }

    // Since BOBYQA solves a minimization problem we want to minimize
    // the *negative* volume to get a maximal volume solution.
    double metric = -sphereSetUnionVolume(set);

    // Replace the current sphere back into the sphere set
    set.spheres[data.sphereIndex] = currentSphere;

    return metric;
}

void sphereFitting(SphereSet& set)
{
    // To avoid fitting the spheres in the same order every time.
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(set.spheres.begin(), set.spheres.end(), std::default_random_engine(seed));

    for (size_t sphereIdx = 0; sphereIdx < set.spheres.size(); ++sphereIdx) {

        Sphere& sphere = set.spheres[sphereIdx];
        assert(sphere.radius > 0.0);

        // Origin (3) + radius (1)
        constexpr INTEGER n = 3 + 1;

        // NPT is the number of interpolation conditions.  Its value must be in the interval
        // [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.
        // TODO: Find a good value!
#if 0
        constexpr INTEGER npt = n + 2; // (6)
#else
        constexpr INTEGER npt = (n + 1) * (n + 2) / 2; // (15)
#endif

        REAL x[n] = {
            sphere.origin.x,
            sphere.origin.y,
            sphere.origin.z,
            sphere.radius
        };

        // To inscribe a box = {min, max} in a sphere we want to find the min & max
        // points which lie along the diagonal that goes through the box center.
        //
        //  sqrt(x^2 + y^2 + y^3) = r
        //   since diagonal -> x=y=z = d
        //             sqrt(3d^2) = r
        //                   3d^2 = r^2
        //                      d = sqrt(r^2 / 3)
        //
        double d = std::sqrt(sphere.radius * sphere.radius / 3.0);
        vec3 boxDiagonal = vec3(d, d, d) * 0.99f; // (avoid overstepping!)
        aabb3 boxInSphere { sphere.origin - boxDiagonal, sphere.origin + boxDiagonal };

        REAL xLower[n] = {
            boxInSphere.min.x,
            boxInSphere.min.y,
            boxInSphere.min.z,
            1e-8
        };

        aabb3 gridBounds = set.grids.shell->gridBounds();
        double gridMaxDistance = length(gridBounds.max - gridBounds.min);

        REAL xUpper[n] = {
            boxInSphere.max.x,
            boxInSphere.max.y,
            boxInSphere.max.z,
            gridMaxDistance
        };

        // "Typically, RHOBEG should be about one tenth of the greatest expected change to a variable"
        // TODO: Find good value!
        REAL rhoBeg = (2.0 * d) / 10.0;

        // RHOEND should indicate the accuracy that is required in the final values of the variables
        // NOTE: Some models are tiny, others are huge.. So an absolute/constant accuracy is probably
        //  not very helpful here. A fraction of the sphere radius could maybe work as a scale-aware metric?
        // TODO: Find good value!
        REAL rhoEnd = gridMaxDistance * 1e-8;

        //The array W will be used for working space.  Its length must be at least
        // (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
        // will be set to the function value at the solution. */
        REAL workingMemory[(npt + 5) * (npt + n) + 3 * n * (n + 5) / 2];

        INTEGER logLevel = 0;
        INTEGER maxObjFunCalls = 10'000;

        SphereFittingData data { &set, sphereIdx };

        int status = bobyqa(n, npt, sphereFittingObjectiveFunction, (void*)(&data),
                            x, xLower, xUpper, rhoBeg, rhoEnd,
                            logLevel, maxObjFunCalls, workingMemory);

        switch (status) {
        case BOBYQA_SUCCESS: {
            // algorithm converged!
            vec3 newOrigin { x[0], x[1], x[2] };
            double newRadius = x[3];

#if 0
            fmt::print("      - ({}, {}, {}) r={}\n", newOrigin.x, newOrigin.y, newOrigin.z, newRadius);

            vec3 diffOrigin = newOrigin - sphere.origin;
            double diffRadius = newRadius - sphere.radius;
            //fmt::print("      - ({}, {}, {}) r={}\n", diffOrigin.x, diffOrigin.y, diffOrigin.z, diffRadius);
            //assert(length(diffOrigin) > 1e-6 && abs(diffRadius > 1e-6));
#endif

            sphere.origin = newOrigin;
            sphere.radius = newRadius;
            return;
        }
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

vec3 randomPointInSphere(const Sphere& sphere)
{
    std::uniform_real_distribution<double> randomOffset { -sphere.radius, +sphere.radius };

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    vec3 offset;
    do {
        offset = { randomOffset(generator), randomOffset(generator), randomOffset(generator) };
    } while (length(offset) > sphere.radius);

    return sphere.origin + offset;
}

void sphereTeleportation(SphereSet& set)
{
    if (set.spheres.size() == 1) {
        // TODO: Maybe pick a new starting point?
        // Or let the teleportation fail like this and restart from the main loop?
        return;
    }

    std::sort(set.spheres.begin(), set.spheres.end(), [](const Sphere& lhs, const Sphere& rhs) {
        return lhs.radius < rhs.radius;
    });

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator { seed };

    size_t numToReplace = std::ceil(set.spheres.size() / 8.0);
    for (size_t i = 0; i < numToReplace; ++i) {

        Sphere& sphere = set.spheres[i];
        sphere.origin = getRandomSphereLocationInsideMesh(*set.grids.inside, generator);

        // NOTE: A sphere expansion will happen after this is done, so we don't need to take care of that now
        sphere.radius = 0.0;
    }
}

double aabbVolume(const aabb3& aabb)
{
    glm::vec<3, double> dims = aabb.max - aabb.min;
    return dims.x * dims.y * dims.z;
}

std::optional<vec3> triangleRayIntersection(vec3 v0, vec3 v1, vec3 v2, vec3 rayOrigin, vec3 rayDirection, float& outT)
{
    // From GraphicsCodex

    const vec3& e1 = v1 - v0;
    const vec3& e2 = v2 - v0;

    vec3 N = normalize(cross(e1, e2));

    vec3 q = cross(rayDirection, e2);
    float a = dot(e1, q);

    // Close to the limit of precision?
    if (fabsf(a) <= 1e-4f) {
        return {};
    }

    const vec3& s = (rayOrigin - v0) / a;
    const vec3& r = cross(s, e1);

    // Barycentric coordinates
    float b[3];
    b[0] = dot(s, q);
    b[1] = dot(r, rayDirection);
    b[2] = 1.0f - b[0] - b[1];

    // Intersected inside triangle?
    float t = dot(e2, r);
    if ((b[0] >= 0) && (b[1] >= 0) && (b[2] >= 0) && t >= 0) {

        vec3 hitPoint = rayOrigin + t * rayDirection;

        // TODO: This worked on a previous project... What's up with the order though?
        vec3 barycentric { b[2], b[0], b[1] };

        outT = t;
        return barycentric;
    }

    return {};
}

vec3 raytraceMesh(const SimpleMesh& mesh, vec3 o, vec3 d)
{
    float closestT = INFINITE;
    int64_t closestTriangle = -1;
    vec3 closestBarycentric {};

    for (size_t ti = 0; ti < mesh.triangleCount(); ++ti) {

        vec3 v0, v1, v2;
        mesh.triangle(ti, v0, v1, v2);

        float t;
        std::optional<vec3> bary = triangleRayIntersection(v0, v1, v2, o, d, t);

        if (bary.has_value() && t < closestT) {
            closestBarycentric = bary.value();
            closestTriangle = ti;
            closestT = t;
        }
    }

    if (closestTriangle == -1) {
        return vec3(0, 0, 0);
    }

    vec2 uv0, uv1, uv2;
    mesh.triangleTexcoords(closestTriangle, uv0, uv1, uv2);

    vec3 b = closestBarycentric;
    vec2 uv = b[0] * uv0 + b[1] * uv1 + b[2] * uv2;

    assert(mesh.hasTexture());
    vec3 color = mesh.texture().sample(uv);

    return color;
}

SphericalHarmonics generateSphericalHarmonics(vec3 sampleCenter, size_t numSamples, const SimpleMesh& mesh)
{

    auto rgbFunction = [&](double phi, double theta) -> vec3 {
        Eigen::Vector3d v = sh::ToVector(phi, theta);
        vec3 direction = vec3(v.x(), v.y(), v.z());
        vec3 Lsrgb = raytraceMesh(mesh, sampleCenter, direction);
        return pow(Lsrgb, vec3(2.2f));
    };

    std::vector<double> r = *sh::ProjectFunction(
        2, [&](double phi, double theta) -> double {
            vec3 color = rgbFunction(phi, theta);
            return color.r;
        },
        numSamples);
    std::vector<double> g = *sh::ProjectFunction(
        2, [&](double phi, double theta) -> double {
            vec3 color = rgbFunction(phi, theta);
            return color.g;
        },
        numSamples);
    std::vector<double> b = *sh::ProjectFunction(
        2, [&](double phi, double theta) -> double {
            vec3 color = rgbFunction(phi, theta);
            return color.b;
        },
        numSamples);

    assert(r.size() == 9);
    assert(g.size() == 9);
    assert(b.size() == 9);

    SphericalHarmonics sh;

    sh.L00 = vec3(r[0], g[0], b[0]);

    sh.L1_1 = vec3(r[1], g[1], b[1]);
    sh.L10 = vec3(r[2], g[2], b[2]);
    sh.L11 = vec3(r[3], g[3], b[3]);

    sh.L2_2 = vec3(r[4], g[4], b[4]);
    sh.L2_1 = vec3(r[5], g[5], b[5]);
    sh.L20 = vec3(r[6], g[6], b[6]);
    sh.L21 = vec3(r[7], g[7], b[7]);
    sh.L22 = vec3(r[8], g[8], b[8]);

    return sh;
}

std::vector<SphericalHarmonics> generateSphericalHarmonicsForSphereColors(const SphereSet& set, size_t numSamples)
{
    std::vector<SphericalHarmonics> SHs;

    for (const Sphere& sphere : set.spheres) {
        SphericalHarmonics sh = generateSphericalHarmonics(sphere.origin, numSamples, set.mesh);
        SHs.push_back(sh);
    }

    return SHs;
}

void generateJsonOutput(const SphereSet& sphereSet, const std::string& outPath)
{
    using namespace nlohmann;

    json j;
    j["proxy"] = "sphere-set";
    j["spheres"] = {};

    for (size_t i = 0; i < sphereSet.spheres.size(); ++i) {

        const Sphere& sphere = sphereSet.spheres[i];
        const SphericalHarmonics& sh = sphereSet.sphereSH[i];

        json jsonSphere = {
            { "center", { sphere.origin.x, sphere.origin.y, sphere.origin.z } },
            { "radius", sphere.radius },
            { "sh",
              { { "L00", { sh.L00.x, sh.L00.y, sh.L00.z } },
                { "L1_1", { sh.L1_1.x, sh.L1_1.y, sh.L1_1.z } },
                { "L10", { sh.L10.x, sh.L10.y, sh.L10.z } },
                { "L11", { sh.L11.x, sh.L11.y, sh.L11.z } },
                { "L2_2", { sh.L2_2.x, sh.L2_2.y, sh.L2_2.z } },
                { "L2_1", { sh.L2_1.x, sh.L2_1.y, sh.L2_1.z } },
                { "L20", { sh.L20.x, sh.L20.y, sh.L20.z } },
                { "L21", { sh.L21.x, sh.L21.y, sh.L21.z } },
                { "L22", { sh.L22.x, sh.L22.y, sh.L22.z } } } }
        };
        j["spheres"].push_back(jsonSphere);
    }

    std::ofstream outstream(outPath);
    outstream << std::setw(4) << j << std::endl;
}

void generateOutput(const SphereSet& sphereSet, bool doPrint = true)
{
    size_t sphereCount = sphereSet.spheres.size();

    std::string result = "const vec4 spheres[] = vec4[](\n";
    for (size_t i = 0; i < sphereCount; ++i) {
        const Sphere& sphere = sphereSet.spheres[i];
        std::string line = fmt::format("\tvec4({}, {}, {}, {}){}\n",
                                       sphere.origin.x, sphere.origin.y, sphere.origin.z,
                                       sphere.radius, (i + 1 == sphereCount) ? "" : ",");
        result += line;
    }
    result += ");\n";

    if (doPrint) {
        fmt::print(result);
    }

#ifdef WIN32
    OpenClipboard(0);
    EmptyClipboard();
    HGLOBAL hg = GlobalAlloc(GMEM_MOVEABLE, result.size() + 1);
    if (!hg) {
        CloseClipboard();
        return;
    }
    memcpy(GlobalLock(hg), result.c_str(), result.size());
    GlobalUnlock(hg);
    SetClipboardData(CF_TEXT, hg);
    CloseClipboard();
    GlobalFree(hg);
#endif
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

    std::string path = "assets/Barrel/barrel.gltf";
    std::string outPath = "assets/Barrel/barrel_spheres.json";

    constexpr unsigned numSpheres = 2;
    const size_t gridDimensions = 126;

    constexpr size_t numShSamples = 4096;

    fmt::print("=> loading model '{}'\n", path);
    auto [basePath, model] = GltfUtil::loadModel(path);
    std::vector<SimpleMesh> simpleMeshes {};
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, simpleMeshes);
    if (simpleMeshes.size() > 1) {
        fmt::print("==> model has more than one meshes; only the first will be used!\n");
    }

    const SimpleMesh& mesh = simpleMeshes[0];
    SphereSet set { mesh };

    fmt::print("=> creating voxel grids\n");

    aabb3 meshBounds { vec3(+INFINITY), vec3(-INFINITY) };
    mesh.extendAABB(meshBounds.min, meshBounds.max);

    set.grids.shell = std::make_unique<VoxelGrid>(glm::ivec3(gridDimensions), meshBounds);
    set.grids.shell->insertMesh(mesh, true);
    set.grids.shell->quantizeColors(10);

    set.grids.filled = std::make_unique<VoxelGrid>(*set.grids.shell);
    set.grids.filled->fillVolumes();

    set.grids.inside = std::make_unique<VoxelGrid>(*set.grids.filled);
    set.grids.inside->subtractGrid(*set.grids.shell);

    fmt::print("=> assigning initial spheres\n");
    assignSpheresRandomlyInVolume(set, numSpheres);

    fmt::print("=> optimizing begin\n");

    bool didJustTeleport = false;
    double previousVolume = -std::numeric_limits<double>::infinity();

    double bestVolume = std::numeric_limits<double>::min();
    std::vector<Sphere> bestSolution {};

    constexpr int maxConsecutiveReverts = 3;
    int numConsecutiveReverts = 0;

    constexpr int maxIterations = 10'000;
    int iteration = 0;

    for (; iteration < maxIterations; ++iteration) {

        fmt::print("==> expanding spheres\n");
        expandSpheresMaximally(set);
        double newVolume = sphereSetUnionVolume(set);
        fmt::print("   volume: {}\n", newVolume);

        if (newVolume > bestVolume) {
            bestVolume = newVolume;
            bestSolution = set.spheres;
        }

        if (newVolume > previousVolume) {

            fmt::print("===> sphere fitting\n");
            sphereFitting(set);

            didJustTeleport = false;
            numConsecutiveReverts = 0;

        } else {

            if (didJustTeleport) {
                fmt::print("==> teleportation failed to increase volume, ");
                if (numConsecutiveReverts < maxConsecutiveReverts) {
                    fmt::print("reverting to known good solution\n");
                    set.spheres = bestSolution;
                    numConsecutiveReverts += 1;
                    //continue; (no, we have to teleport from this new stage first!)
                } else {
                    fmt::print("and max reverts reached, so aborting\n");
                    break;
                }
            }

            fmt::print("===> sphere teleportation\n");
            sphereTeleportation(set);
            didJustTeleport = true;
        }

        previousVolume = newVolume;

        if (iteration % 10 == 0) {
            double itDone = 100.0 * double(iteration) / double(maxIterations);
            fmt::print("==> {} iterations done ({:.1f}%)\n", iteration, itDone);
            generateOutput(set, false);
        }
    }

    if (iteration == maxIterations) {
        fmt::print("==> max iterations reached\n");
    }

    if (numConsecutiveReverts == maxConsecutiveReverts) {
        fmt::print("==> max consectutive reverts reached\n");
    }

    fmt::print("=> optimizing done\n");

    fmt::print("=> generating sphere output with volume {} (AABB volume {})\n", bestVolume, aabbVolume(meshBounds));
    set.spheres = bestSolution;

    fmt::print("=> generating SH color representations for spheres ({} samples)\n", numShSamples);
    std::vector<SphericalHarmonics> SHs = generateSphericalHarmonicsForSphereColors(set, numShSamples);
    set.sphereSH = SHs;

    generateOutput(set);
    generateJsonOutput(set, outPath);
}
