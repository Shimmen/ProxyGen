#include "GltfUtil.h"
#include "SimpleMesh.h"
#include "VoxelGrid.h"
#include "mathkit.h"
#include <fmt/format.h>
#include <vector>

int main()
{
    // TODO: Take these as command line parameters!
#if 0
    std::string path = "../assets/Avocado/Avocado.gltf";
    std::string outfile = "../assets/Avocado.vox";
#else
    std::string path = "../assets/BoomBox/BoomBoxWithAxes.gltf";
    std::string outfile = "../assets/BoomBox.vox";
#endif
    size_t gridDimensions = 126;

    auto [basePath, model] = GltfUtil::loadModel(path);

    std::vector<SimpleMesh> simpleMeshes {};
    GltfUtil::bakeDownModelToSimpleMeshes(model, basePath, simpleMeshes);

    aabb3 boundsOfAllMeshes = SimpleMesh::calculateBounds(simpleMeshes);
    VoxelGrid grid { glm::ivec3(gridDimensions), boundsOfAllMeshes };

    fmt::print("= voxelization begin =\n");
    for (size_t i = 0; i < simpleMeshes.size(); ++i) {
        fmt::print("   mesh {}/{}\n", i + 1, simpleMeshes.size());
        grid.insertMesh(simpleMeshes[i], true);
    }
    fmt::print("= voxelization done  =\n");

    fmt::print("= color quantization begin =\n");
    grid.quantizeColors(200);
    fmt::print("= color quantization done  =\n");

    fmt::print("= volume filling begin =\n");
    grid.fillVolumes(simpleMeshes);
    fmt::print("= volume filling done  =\n");

    grid.writeToVox(outfile);
}
