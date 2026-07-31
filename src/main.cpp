#include "polyscope/polyscope.h"
#pragma warning( push )
#pragma warning( disable : 4101)

#include <igl/exact_geodesic.h>
#include <igl/invert_diag.h>
#include "polyscope/image_quantity.h"
#include "polyscope/camera_view.h"

#include <igl/readOBJ.h>
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

// Your function
#pragma warning( pop )
#include <iostream>

#include "sharp-feature/sharp_feature_interface.h"

std::vector<glm::vec3>  GetFaceClustering(const std::vector<int>& temp)
{
    std::vector<glm::vec3> COLOR_MASKS(22);
    COLOR_MASKS[0] = {0.801, 0.078, 0.254};  // Darker Red
    COLOR_MASKS[1] = {0.185, 0.605, 0.254};  // Darker Green
    COLOR_MASKS[2] = {0.8, 0.678, 0.0};      // Darker Yellow
    COLOR_MASKS[3] = {0.162, 0.288, 0.747};  // Darker Blue
    COLOR_MASKS[4] = {0.760, 0.409, 0.152};  // Darker Orange
    COLOR_MASKS[5] = {0.468, 0.097, 0.605};  // Darker Purple
    COLOR_MASKS[6] = {0.174, 0.741, 0.741};  // Darker Cyan
    COLOR_MASKS[7] = {0.741, 0.096, 0.701};  // Darker Magenta
    COLOR_MASKS[8] = {0.537, 0.764, 0.0};    // Darker Lime
    COLOR_MASKS[9] = {0.780, 0.545, 0.545};  // Darker Pink
    COLOR_MASKS[10] = {0.0, 0.401, 0.401};   // Darker Teal
    COLOR_MASKS[11] = {0.701, 0.545, 0.8};   // Darker Lavender
    COLOR_MASKS[12] = {0.503, 0.288, 0.101}; // Darker Brown
    COLOR_MASKS[13] = {0.8, 0.780, 0.584};   // Darker Beige
    COLOR_MASKS[14] = {0.401, 0.0, 0.0};     // Darker Maroon
    COLOR_MASKS[15] = {0.566, 0.8, 0.664};   // Darker Mint
    COLOR_MASKS[16] = {0.401, 0.401, 0.0};   // Darker Olive
    COLOR_MASKS[17] = {0.8, 0.647, 0.494};   // Darker Peach
    COLOR_MASKS[18] = {0.0, 0.0, 0.358};     // Darker Navy
    COLOR_MASKS[19] = {0.401, 0.401, 0.401}; // Darker Gray
    COLOR_MASKS[20] = {1., 1., 1.};       // Dark Gray{"farClipRatio":20.0,"fov":28.8500003814697,"nearClipRatio":0.005,"projectionMode":"Orthographic","viewMat":[1.0,0.0,0.0,-8.44820690155029,0.0,0.997785210609436,-0.0665190368890762,-0.526245176792145,0.0,0.0665190368890762,0.997785210609436,-7.09799098968506,0.0,0.0,0.0,1.0],"windowHeight":1061,"windowWidth":1924}
    COLOR_MASKS[21] = {0.0, 0.0, 0.0};       // Black (was white)

    std::vector<glm::vec3> faceColorCluster(temp.size());
    for (int i = 0; i < temp.size(); i++)
    {
        faceColorCluster[i] = COLOR_MASKS[temp[i]%22];
    }
    return faceColorCluster;
}

// imgui
#include "imfilebrowser.h"
ImGui::FileBrowser fileDialog;

std::string filePath;
Eigen::MatrixXd meshV;
Eigen::MatrixXi meshF;
Mesh surfaceMesh;
std::string coolMeshName;
std::vector<glm::vec3> vertexPositions;
std::vector<std::vector<size_t>> triFaces;

void openMesh(const std::string& filename)
{
    surfaceMesh = Mesh();
    filePath = filename;
    if (!CGAL::IO::read_polygon_mesh(filename, surfaceMesh)) {
        std::cerr << "Invalid input file." << std::endl;

        // try if the input is imge
        //surfaceMesh = load_and_mesh(fileDialog.GetSelected().string(), vertexPositions, triFaces);
    } else {

        // normalize
        //normalize_mesh(surfaceMesh);
        export_surface_mesh_to_vectors<Point>(surfaceMesh, vertexPositions, triFaces);
    }
    polyscope::removeAllStructures();
    polyscope::removeStructure(coolMeshName);
    auto cgalMesh = polyscope::registerSurfaceMesh(coolMeshName, vertexPositions, triFaces);
    cgalMesh->setEdgeWidth(0.7);
    cgalMesh->resetTransform();

    sharp_feature::surfaceMesh = surfaceMesh;
    sharp_feature::coolMeshName = coolMeshName;
}

void callback() {

  ImGui::PushItemWidth(100);
    if (ImGui::Button("Open new mesh"))
        fileDialog.Open();

    fileDialog.Display();

    if (fileDialog.HasSelected()) {
        std::cout << "Selected filename" << fileDialog.GetSelected().string() << std::endl;
        openMesh(fileDialog.GetSelected().string());
        fileDialog.ClearSelected();
    }
    sharp_feature::displayInterface();
    ImGui::PopItemWidth();
}

int main(int argc, char **argv) {
    // Options
    polyscope::options::autocenterStructures = true;
    polyscope::view::windowWidth = 1920;
    polyscope::view::windowHeight = 1080;

    // Initialize polyscope
    polyscope::init();
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    filePath = std::string("D://library//SharpEdgeSmall//data//cubesmooth2.obj");

    fileDialog.SetTitle("open mesh");
    fileDialog.SetDirectory(filePath);
    fileDialog.SetTypeFilters({ ".obj", ".off", ".inr" });

    coolMeshName = polyscope::guessNiceNameFromPath(filePath);
    //openMesh(filePath);

    std::string testPath("D://Thesis//Ange//SharpFeatureRelease//test_data//cube.obj");

    openMesh(testPath);
    filePath = testPath;
    test_graph_optimizer::test(surfaceMesh);

    // Add the callback
    polyscope::state::userCallback = callback;

    // Show the gui
    polyscope::show();

  return 0;
}
