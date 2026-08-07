#pragma once
#include "type.h"
#include "feature_metric.h"
#include "feature_graph.h"
#include "graph_optimizer.h"
namespace sharp_feature
{
    struct SharpFeatureParameters
    {
        bool displayPoints = true;
        bool displayOriginalGraph = false;
        bool displayCorners = true;

        float sizeCorner = 0.00100;
        float sizePoint = 0.00100;
        float sizeLine = 0.00100;

        int maxIteration = 1;
    };

    template<typename UserType>
    struct PersistentValues
    {
        std::vector<UserType> vertexMeasure;
        std::vector<int> selectedFaces;

        std::vector<Point> featurePoint;
        std::vector<std::vector<int>> featureGraph;

        Polylines<Point> polylines;
        std::vector<double> edgeScalar;
    };
    class SharpFeature
    {
    public:
        SharpFeature() = default;
        SharpFeature(const Mesh& inputMesh)
        {
            surfaceMesh = inputMesh;
            clean(surfaceMesh);
        }
        void computeSharpFeatureGraph()
        {
            auto metric = computeSharpnessMetric();
        }
        std::vector<double> computeSharpnessMetric()
        {
            persistentValues.vertexMeasure = featureMetric.computeSharpMeasure(surfaceMesh, featureMetricParameters);

            return persistentValues.vertexMeasure;
        }

        void computeSubdivision()
        {
            Mesh copy = surfaceMesh;
            Mesh newMesh;
            if(persistentValues.vertexMeasure.empty()) return;
            std::unordered_set<Mesh::Face_index > facesToSubdivise;
            std::vector<int> faceStatus(copy.num_faces(), 0);
            std::unordered_map<Mesh::Edge_index, Mesh::Vertex_index> cache;
            auto getMidPoint = [&](Mesh::halfedge_index h)
            {
                Mesh::Vertex_index midpoint;
                if(!cache.contains(copy.edge(h)))
                {
                    auto p =
                            CGAL::midpoint(copy.point(copy.source(h)), copy.point(copy.target(h)));
                    midpoint = newMesh.add_vertex(p);
                    cache[copy.edge(h)] = midpoint;
                    subdivisionMetric[copy.edge(h)] = 0.;
                }
                else
                {
                    midpoint = cache[copy.edge(h)];
                }
                return midpoint;
            };
            for(auto v :copy.vertices())
            {
                newMesh.add_vertex(copy.point(v));
            }
            for(auto f :copy.faces())
            {
                auto it = vertices_around_face(copy.halfedge(f), copy).begin();

                auto v0 = *it++;
                auto v1 = *it++;
                auto v2 = *it++;

                double m0 = persistentValues.vertexMeasure[v0];
                double m1 = persistentValues.vertexMeasure[v1];
                double m2 = persistentValues.vertexMeasure[v2];

                if(     m0 > featureGraphParameters.threshold &&
                        m1 > featureGraphParameters.threshold &&
                        m2 > featureGraphParameters.threshold )
                {
                    facesToSubdivise.insert(f);
                    faceStatus[f.idx()] = 1;
                }
            }
            for(auto f :copy.faces())
            {
                if(facesToSubdivise.contains(f))
                {
                    // perform triangle subdivision mid point 1 triangle to 4 triangle
                    auto h0 = copy.halfedge(f);
                    auto h1 = copy.next(h0);
                    auto h2 = copy.next(h1);

                    auto v0 = copy.source(h0);
                    auto v1 = copy.source(h1);
                    auto v2 = copy.source(h2);

                    newMesh.add_face(v0,  getMidPoint(h0), getMidPoint(h2));
                    newMesh.add_face(v1,  getMidPoint(h1), getMidPoint(h0));
                    newMesh.add_face(v2,  getMidPoint(h2), getMidPoint(h1));
                    newMesh.add_face(getMidPoint(h0), getMidPoint(h1), getMidPoint(h2));
                }
                else
                {
                    // iterate halfedges if face of opposite has a status of 1 need to split the face in half
                    bool needSplit = false;
                    for(auto h : copy.halfedges_around_face(copy.halfedge(f)))
                    {
                        auto oh = copy.opposite(h);
                        if(!copy.is_border(oh))
                        {
                            if(facesToSubdivise.contains(copy.face(oh)))
                            {
                                auto v0 = copy.source(h);
                                auto v1 = copy.target(h);
                                auto v2 = copy.target(copy.next(h));
                                newMesh.add_face(v0,  getMidPoint(h), v2);
                                newMesh.add_face(v1,  v2, getMidPoint(h));
                                needSplit = true;
                                break;
                            }
                        }
                    }
                    // otherwise push it in new mesh
                    if(!needSplit)
                    {
                        auto it = vertices_around_face(copy.halfedge(f), copy).begin();

                        auto v0 = *it++;
                        auto v1 = *it++;
                        auto v2 = *it++;
                        newMesh.add_face(v0,  v1, v2);
                    }

                }
            }
            std::vector<Point> subPoint;
            for(auto &e : subdivisionMetric)
            {
                e.second = featureMetric.computeSharpMeasureVertex(e.first, surfaceMesh, featureMetricParameters);
                subPoint.push_back(newMesh.point(getMidPoint(surfaceMesh.halfedge(e.first))));
            }
            /*auto subPoints  = polyscope::registerPointCloud("subPoints", subPoint);
            subPoints->resetTransform();

            std::vector<glm::vec3> vertexPositions;
            std::vector<std::vector<size_t>> triFaces;
            export_surface_mesh_to_vectors(newMesh, vertexPositions, triFaces);
            auto sub = polyscope::registerSurfaceMesh("sub", vertexPositions, triFaces);
            sub->resetTransform();*/
        }
        void computeInitialGraph()
        {
            featureGraph = FeatureGraph<double>();
            computeSubdivision();
            featureGraph.computeSelectedSurfels(surfaceMesh, persistentValues.vertexMeasure, featureGraphParameters);
            persistentValues.selectedFaces = featureGraph.getSelectedFaces();


            featureGraph.computeFeatureGraph(surfaceMesh, persistentValues.vertexMeasure, subdivisionMetric);
            if(sharpFeatureParameters.displayOriginalGraph)
            {
                std::vector<Point> points;
                std::vector<std::array<size_t, 2>> edges;
                std::vector<double> metrics;
                featureGraph.exportGraph(points, edges, metrics);
                auto curveSimplified = polyscope::registerCurveNetwork("feature graph", points, edges);
                curveSimplified->addNodeScalarQuantity("metric", metrics);
                curveSimplified->resetTransform();
                curveSimplified->setRadius(sharpFeatureParameters.sizeLine);
            }

            featureGraph.computeThinning(surfaceMesh);
            std::vector<Point> points;
            std::vector<std::array<size_t, 2>> edges;
            std::vector<double> metrics;
            featureGraph.exportGraph(points, edges, metrics);

           /* auto curveSimplified = polyscope::registerCurveNetwork("feature graph simplified", points, edges);
            curveSimplified->resetTransform();
            curveSimplified->setRadius(sharpFeatureParameters.sizeLine);

            if (sharpFeatureParameters.displayPoints) {
                auto pcl = polyscope::registerPointCloud("feature_graph_points", points);
                pcl->resetTransform();
            }
            */

            featureGraph.exportGraphSimple(persistentValues.featurePoint, persistentValues.featureGraph);
        }
        void computePolyLines()
        {
            graphOptimizer.computePolyLines(persistentValues.featurePoint, persistentValues.featureGraph);
            persistentValues.polylines = graphOptimizer.getPolyLines();

            if (sharpFeatureParameters.displayCorners) {
                std::vector<Point> corners = graphOptimizer.getCorners();
                auto pcl = polyscope::registerPointCloud("corners", corners);
                pcl->resetTransform();
            }
            if (sharpFeatureParameters.displayPoints) {
                auto pcl = polyscope::registerPointCloud("polylines_points", persistentValues.polylines.points);
                pcl->resetTransform();
            }
        }
        void computeFacePatches()
        {
            //graphOptimizer.optimizePolyLines();
            graphOptimizer.computeFacePatches(surfaceMesh);
            persistentValues.edgeScalar = graphOptimizer.edgeFrontiers;

        }
        void noisePolyLines()
        {
            graphOptimizer.noisePolyLines(surfaceMesh);
            persistentValues.polylines = graphOptimizer.getPolyLines();
        }
        void optimizePolyLines()
        {

            graphOptimizer.maxIteration = sharpFeatureParameters.maxIteration;
            graphOptimizer.optimizePoints(surfaceMesh);
            persistentValues.polylines = graphOptimizer.getPolyLines();
        }
        std::vector<int> getSelectedPatch(const int& selectedVert)
        {
            std::vector<int> res(surfaceMesh.number_of_faces(), 0);
            if(selectedVert >= graphOptimizer.labelFaces.size())
            {
                return {};
            }
            auto faces = graphOptimizer.labelFaces[selectedVert];
            int i = 0;
            for(auto& faceID : graphOptimizer.faceNeighbors[selectedVert])
            {
                res[faceID.idx()] = 1 + faces[i];
                i++;
            }
            return res;

        }
        std::vector<Point> getOptimal()
        {
            return graphOptimizer.optimalPoints;
        }
        std::vector<int> getVertexDegree()
        {
            return graphOptimizer.vertexDegree;
        }
        std::vector<std::vector<Point>> getPlaneCentroids()
        {
            return graphOptimizer.planeCentroids;
        }
        std::vector<std::vector<Vector >>  getPlaneNormals()
        {
            return graphOptimizer.planeNormals;
        }


        Mesh getSurfaceMesh()
        {
            return surfaceMesh;
        }
        FeatureMetricParameters<double> featureMetricParameters;
        FeatureGraphParameters<double> featureGraphParameters;
        PersistentValues<double> persistentValues;

        SharpFeatureParameters sharpFeatureParameters;

        bool hasChanged = false;
    private:
        void clean(Mesh& sm)
        {
            std::vector<Mesh::Vertex_index> to_remove;
            for (auto v : sm.vertices()) {
                if (sm.halfedge(v) == Mesh::null_halfedge())
                    to_remove.push_back(v);
            }
            for (auto v : to_remove)
            {
                hasChanged = true;
                sm.remove_vertex(v);
            }

            sm.collect_garbage();
        }
        FeatureMetric<double> featureMetric;
        FeatureGraph<double> featureGraph;
        GraphOptimizer graphOptimizer;
        Mesh surfaceMesh;

        std::unordered_map<Mesh::Edge_index , double> subdivisionMetric;
    };
};