#include "sharp_feature.h"
#include "test_optimizer.h"
#include "polyscope/combining_hash_functions.h"

namespace sharp_feature
{


    static MethodMetricType currentMethod = DIHEDRAL;

    Mesh surfaceMesh;
    std::string coolMeshName;
    SharpFeature sharpFeature;

    FeatureMetricParameters<double> featureMetricParameters;
    FeatureGraphParameters<double> featureGraphParameters;

    float selectionThreshold = static_cast<float>(featureGraphParameters.threshold);


    void display(const Polylines<Point>& polylines, const std::string& name = "") {

        polyscope::Group * polyGroup;
        bool groupExists = polyscope::state::groups.find("PolyLines"+ name) != polyscope::state::groups.end();

        if (!groupExists)
        {
            polyGroup = polyscope::createGroup("PolyLines"+ name);
        }
        else
        {
            polyGroup = polyscope::getGroup("PolyLines"+ name);
        }

        std::vector<int> pointID(polylines.points.size(),0);
        size_t polyLineID = 0;
        for (const auto& polyline : polylines.lines)
        {
            if (polyline.size() < 2)
                continue;

            std::vector<Point> points;
            for (size_t i = 0; i  < polyline.size(); ++i)
            {
                points.push_back(polylines.points[polyline[i]]);
            }
            std::vector<std::array<size_t, 2>> edges;
            for (size_t i = 0; i + 1 < polyline.size(); ++i)
            {
                edges.push_back({
                                        i,
                                        i + 1
                                });
            }

            auto curveSimplified = polyscope::registerCurveNetwork(name+"polyline_"+std::to_string(polyLineID),   points, edges);
            curveSimplified->resetTransform();
            curveSimplified->setRadius(sharpFeature.sharpFeatureParameters.sizeLine);
            curveSimplified->addToGroup(*polyGroup);
            polyLineID++;
        }
        auto pcl = polyscope::registerPointCloud("polylines_points_"+name, polylines.points);
        pcl->resetTransform();

        // use polylines.points directly

    }
    void displayInterface()
    {
      ImGui::SetNextItemOpen(true);
      if (ImGui::TreeNode("Sharp Feature")) {

          ImGui::SetNextItemOpen(true);
          if (ImGui::TreeNode("FeatureMetric"))
          {
              if (ImGui::CollapsingHeader("Sharp Feature Measure Parameters", ImGuiTreeNodeFlags_DefaultOpen))
              {
                  ImGui::Checkbox("useAreaWeight", &featureMetricParameters.useAreaWeight);
                  ImGui::Checkbox("useAngleWeight", &featureMetricParameters.useAngleWeight);
                  ImGui::Separator();
                  ImGui::Checkbox("useGeodesicWeight", &featureMetricParameters.useGeodesicWeight);
                  ImGui::DragFloat("Sigma", &featureMetricParameters.sigma,
                                   0.01f, 0.001f, 10.0f, "%.3f");
                  ImGui::DragInt("nSamples", &featureMetricParameters.nSamples,
                                  1, 128, 1);
                  ImGui::DragFloat("Geodesic Radius", &featureMetricParameters.geodesicRadius,
                                   0.01f, 0.001f, 10.0f, "%.3f");
                  ImGui::Separator();

              }
              if (ImGui::BeginCombo("Feature Metric Method", MethodMetricNames[currentMethod]))
              {
                  for (int i = 0; i < MethodMetricType::METHOD_COUNT; ++i)
                  {
                      bool selected = (currentMethod == i);
                      if (ImGui::Selectable(MethodMetricNames[i], selected))
                      {
                          currentMethod = static_cast<MethodMetricType>(i);
                          featureMetricParameters.currentMethod = currentMethod;
                      }

                      if (selected)
                          ImGui::SetItemDefaultFocus();
                  }
                  ImGui::EndCombo();
              }
              if(ImGui::Button("Run compute metric"))
              {
                  sharpFeature = SharpFeature(surfaceMesh);
                  sharpFeature.featureMetricParameters = featureMetricParameters;
                  auto metric = sharpFeature.computeSharpnessMetric();

                  if(sharpFeature.hasChanged)
                  {
                      std::vector<glm::vec3> vertexPositions;
                      std::vector<std::vector<size_t>> triFaces;

                      surfaceMesh = sharpFeature.getSurfaceMesh();
                      export_surface_mesh_to_vectors<Point>(surfaceMesh, vertexPositions, triFaces);
                      auto cgalMesh = polyscope::registerSurfaceMesh(coolMeshName, vertexPositions, triFaces);
                      cgalMesh->setEdgeWidth(0.7);
                      cgalMesh->resetTransform();
                  }
                  auto m = polyscope::getSurfaceMesh(coolMeshName)->addVertexScalarQuantity("Metric "+ std::string(MethodMetricNames[currentMethod]), metric);
                  m->setEnabled(true);
              }
              ImGui::TreePop();
          }

          ImGui::SetNextItemOpen(true);
          if (ImGui::TreeNode("FeatureGraph"))
          {
              if (ImGui::SliderFloat(
                      "Selection Threshold",
                      &selectionThreshold,
                      0.0f,
                      1.0f))
              {
                  featureGraphParameters.setThreshold(selectionThreshold);
              }
              ImGui::InputInt("maxIteration", &sharpFeature.sharpFeatureParameters.maxIteration);
              if (ImGui::CollapsingHeader("Sharp Feature Display Parameters", ImGuiTreeNodeFlags_DefaultOpen))
              {
                  ImGui::Checkbox("Display Points", &sharpFeature.sharpFeatureParameters.displayPoints);
                  ImGui::Checkbox("Display Original Graph", &sharpFeature.sharpFeatureParameters.displayOriginalGraph);
                  ImGui::Checkbox("Display Corners", &sharpFeature.sharpFeatureParameters.displayCorners);

                  ImGui::Separator();

                  ImGui::DragFloat("Corner Size", &sharpFeature.sharpFeatureParameters.sizeCorner,
                                   0.01f, 0.001f, 10.0f, "%.3f");
                  ImGui::DragFloat("Point Size", &sharpFeature.sharpFeatureParameters.sizePoint,
                                   0.01f, 0.001f, 10.0f, "%.3f");
                  ImGui::DragFloat("Line Size", &sharpFeature.sharpFeatureParameters.sizeLine,
                                   0.01f, 0.001f, 10.0f, "%.3f");
              }
              if(ImGui::Button("Build initial feature graph"))
              {
                  sharpFeature.featureGraphParameters = featureGraphParameters;
                  sharpFeature.computeInitialGraph();

                  polyscope::getSurfaceMesh(coolMeshName)->addFaceScalarQuantity("selected", sharpFeature.persistentValues.selectedFaces);
              }
              if(ImGui::Button("Compute polylines"))
              {
                  sharpFeature.computePolyLines();

                  //display(sharpFeature.persistentValues.polylines);
              }

              ImGui::TreePop();
          }
          ImGui::SetNextItemOpen(true);
          if (ImGui::TreeNode("Graph Optimizer"))
          {
              if(ImGui::Button("Optimize polylines"))
              {

                  sharpFeature.computeFacePatches();
                  display(sharpFeature.persistentValues.polylines);
                  auto vec = sharpFeature.getVertexDegree();
                  polyscope::getSurfaceMesh(coolMeshName)->addVertexScalarQuantity("degree", vec);
                  auto pcl = polyscope::registerPointCloud("optimal", sharpFeature.getOptimal());
                  pcl->resetTransform();
                  std::vector<int> edgeOrdering;
                  std::unordered_set<std::pair<size_t, size_t>, polyscope::hash_combine::hash<std::pair<size_t, size_t>>> seenEdges;

                  for(auto f : surfaceMesh.faces())
                  {
                      for(auto h : surfaceMesh.halfedges_around_face(surfaceMesh.halfedge(f)))
                      {
                          auto i0 = surfaceMesh.source(h).idx();
                          auto i1 = surfaceMesh.target(h).idx();

                          size_t iMin = std::min(i0, i1);
                          size_t iMax = std::max(i0, i1);
                          auto p = std::make_pair(iMin, iMax);
                          if (seenEdges.find(p) == seenEdges.end()) {
                              edgeOrdering.push_back(edgeOrdering.size());
                              seenEdges.insert(p);
                          }
                      }
                  }


                  if(polyscope::getSurfaceMesh(coolMeshName)->edgePerm.empty())
                      polyscope::getSurfaceMesh(coolMeshName)->setEdgePermutation(edgeOrdering);
                  polyscope::getSurfaceMesh(coolMeshName)->addEdgeScalarQuantity("line", sharpFeature.persistentValues.edgeScalar);

                  bool displayPlanes = false;
                  if(displayPlanes)
                  {

                      auto centroids = sharpFeature.getPlaneCentroids();
                      auto planeNormals = sharpFeature.getPlaneNormals();

                      std::vector<Point> centers;
                      std::vector<Vector> normals;
                      /*for(int i = 0 ; i < centroids.size(); i++)
                      {
                          // for each point of polyline
                          auto pcl2 = polyscope::registerPointCloud("planeCentroids_"+std::to_string(i), centroids[i]);
                          pcl2->addVectorQuantity("normal", planeNormals[i]);
                          pcl2->resetTransform();
                      }*/
                      for(int i = 0 ; i < centroids.size(); i++)
                      {
                        centers.insert(centers.end(), centroids[i].begin(), centroids[i].end());
                          normals.insert(normals.end(), planeNormals[i].begin(), planeNormals[i].end());
                      }
                      auto pcl2 = polyscope::registerPointCloud("planeCentroids", centers);
                      pcl2->addVectorQuantity("normal", normals);
                      pcl2->resetTransform();

                  }

                  sharpFeature.noisePolyLines();
                  display(sharpFeature.persistentValues.polylines, "_beforeOpti");
                  sharpFeature.optimizePolyLines();


                  /*sharpFeature.computeFacePatches();
                  sharpFeature.optimizePolyLines();*/
                  display(sharpFeature.persistentValues.polylines, "_afterOpti");

              }
              if (ImGui::CollapsingHeader("Debug", ImGuiTreeNodeFlags_DefaultOpen))
              {
                  if (ImGui::Button("select logmap source")) {

                      long long int pickVert = polyscope::getSurfaceMesh(coolMeshName)->selectVertex();
                      if (pickVert >= 0) {
                          auto selectedPatch = sharpFeature.getSelectedPatch(pickVert);
                          std::cout<<"pickVert "<<pickVert<<"\n";
                          if(!selectedPatch.empty())
                          {
                              auto p = polyscope::getSurfaceMesh(coolMeshName)->addFaceScalarQuantity("facePatch", selectedPatch);
                              p->setEnabled(true);
                          }

                      }
                  }
              }

              ImGui::TreePop();
          }
          ImGui::TreePop();
      }
        ImGui::SetNextItemOpen(true);
        if (ImGui::TreeNode("TEST Graph Optimizer"))
        {
            if(ImGui::Button("run test"))
            {
                test_graph_optimizer::test(surfaceMesh);
            }
            ImGui::TreePop();
        }
    }
};