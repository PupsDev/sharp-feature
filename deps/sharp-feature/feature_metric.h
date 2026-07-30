#include <vector>
#include "sampling_fibo.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/curvature.h>
namespace PMP = CGAL::Polygon_mesh_processing;
namespace sharp_feature
{
    template<typename UserType>
    class FeatureMetric
    {
    public:
        FeatureMetric() = default;
        UserType computeSharpMeasureVertex(const Mesh::Edge_index& edge, const Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {
            switch (parameters.currentMethod) {
                case MethodMetricType::DIHEDRAL:
                    return computeDihedralEdge(edge, sm, parameters);
                    break;
                default:
                    return computeDihedralEdge(edge, sm, parameters);
                break;
            }

        }
        std::vector<UserType> computeSharpMeasure(Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {
            switch (parameters.currentMethod)
            {
                case MethodMetricType::DIHEDRAL:
                    return computeMeasureDihedral(sm, parameters);
                    break;
                case MethodMetricType::DISCRETE_GAUSSIAN:
                    return computeGaussianCurvature(sm, parameters);
                    break;
                case MethodMetricType::DISCRETE_MEAN:
                    return computeDiscreteMeanCurvature(sm, parameters);
                    break;
                case MethodMetricType::CURVEDNESS:
                    return computeCurvedness(sm, parameters);
                    break;
                case MethodMetricType::NVT:
                    return computeNVT(sm, parameters);
                    break;
                case MethodMetricType::NVT_Geodesic:
                    return computeNVTGeodesic(sm, parameters);
                    break;
                default:
                    return computeMeasureDihedral(sm, parameters);
                    break;
            }
        }
    private:
        UserType computeDihedralEdge(const Mesh::Edge_index& edge, const Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {
            UserType value;
            auto h1 = sm.halfedge(edge);
            if (!sm.is_border(h1) && !sm.is_border(sm.opposite(h1))) {
                const auto &faceNormal = CGAL::Polygon_mesh_processing::compute_face_normal(sm.face(h1), sm);
                if(sm.face(sm.opposite(h1)) == Mesh::null_face()) return 1.;
                const auto &faceNormalOpposite = CGAL::Polygon_mesh_processing::compute_face_normal(
                        sm.face(sm.opposite(h1)), sm);
                const auto &angle =
                        static_cast<UserType>(1.) - abs(CGAL::scalar_product(faceNormal, faceNormalOpposite));
                value = (UserType(angle));
            } else {
                value = (UserType(0.));
            }
            return value;
        }
        UserType computeNVTEdge(const Mesh::Edge_index& edge, const Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {
            /*Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
            for(auto face: sm.faces_around_face(sm.halfedge(edge)))
            {
                const auto faceTensor = computeTensor(normal);
                tensor +=  faceTensor;
            }
            auto metric = solveTensor(tensor);
            vertexMetrics.push_back(static_cast<double>(metric));*/
        }

        std::vector<UserType> computeNVT(
                Mesh& sm,
                const FeatureMetricParameters<UserType>& parameters)
        {
            vertexMetrics.clear();


            std::vector<std::vector<Vector>> normals(sm.number_of_vertices());
            std::vector<std::vector<Mesh::face_index>> faceNeighbors(sm.number_of_vertices());

            for (const auto &vID: sm.vertices()) {
                for (const auto &h: CGAL::halfedges_around_target(sm.halfedge(vID), sm)) {
                    const auto &vn = sm.source(h);
                    const auto &fn = sm.face(h);
                    if (fn == Mesh::null_face()) continue;
                    faceNeighbors[vID.idx()].push_back(fn);
                }
            }

            size_t vID = 0 ;
            for(const auto& neighbors: faceNeighbors)
            {
                for(auto& faceID : neighbors)
                    normals[vID].push_back( CGAL::Polygon_mesh_processing::compute_face_normal(faceID, sm));
                vID++;
            }
            for(vID = 0 ; vID< sm.number_of_vertices(); vID++)
            {
                Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();
                for(auto normal : normals[vID])
                {
                    const auto faceTensor = computeTensor(normal);
                    tensor +=  faceTensor;
                }
                auto metric = solveTensor(tensor);
                vertexMetrics.push_back(static_cast<double>(metric));
            }


            return vertexMetrics;
        }
        std::vector<UserType> computeNVTGeodesic(
                Mesh& sm,
                const FeatureMetricParameters<UserType>& parameters)
        {
            NVTGeodesicFibonacci<double> featureNVTGeodesic;
            featureNVTGeodesic.set(sm);
            vertexMetrics = featureNVTGeodesic.computeNVTGeodesic(parameters);

            auto initSampling = polyscope::registerPointCloud("sampling initial" , featureNVTGeodesic.samplingInitialPerVertex[0]);
            initSampling->resetTransform();

            auto sampling = polyscope::registerPointCloud("sampling" , featureNVTGeodesic.samplingPerVertex[0]);
            sampling->resetTransform();

            return vertexMetrics;
        }
        std::vector<UserType> computeCurvedness(
                Mesh& sm,
                const FeatureMetricParameters<UserType>& parameters)
        {
            vertexMetrics.clear();
            vertexMetrics.resize(sm.number_of_vertices(), UserType(0));

            for (const auto& v : sm.vertices())
            {
                double H = CGAL::to_double(
                        PMP::discrete_mean_curvature(v, sm));

                double K = CGAL::to_double(
                        PMP::discrete_Gaussian_curvature(v, sm));

                // Numerical robustness
                double C2 = 2.0 * H * H - K;
                C2 = std::max(0.0, C2);

                vertexMetrics[v.idx()] =
                        static_cast<UserType>(std::sqrt(C2));
            }

            return vertexMetrics;
        }
        std::vector<UserType> computeDiscreteMeanCurvature(
                Mesh& sm,
                const FeatureMetricParameters<UserType>& parameters)
        {
            vertexMetrics.assign(sm.number_of_vertices(), UserType(0));

            for (const auto& v : sm.vertices())
            {
                vertexMetrics[v.idx()] =
                        static_cast<UserType>(
                                std::abs(CGAL::to_double(
                                        PMP::discrete_mean_curvature(v, sm))));
            }

            return vertexMetrics;
        }
        std::vector<UserType> computeGaussianCurvature(Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {
            vertexMetrics.clear();
            vertexMetrics.resize(sm.number_of_vertices());

            for (const auto& v : sm.vertices())
            {
                vertexMetrics[v.idx()] = static_cast<UserType>(
                        std::abs(CGAL::to_double(
                                PMP::discrete_Gaussian_curvature(v, sm))));
            }

            return vertexMetrics;
        }
        std::vector<UserType> computeMeasureDihedral(Mesh& sm, const FeatureMetricParameters<UserType>& parameters)
        {

            std::vector <UserType> metrics;
            for (const auto &edge: sm.halfedges()) {
                if (!sm.is_border(edge) && !sm.is_border(sm.opposite(edge))) {
                    const auto &faceNormal = CGAL::Polygon_mesh_processing::compute_face_normal( sm.face(edge), sm);
                    const auto &faceNormalOpposite = CGAL::Polygon_mesh_processing::compute_face_normal(sm.face(sm.opposite(edge)), sm);
                    const auto &angle = static_cast<UserType>(1.) - abs(CGAL::scalar_product(faceNormal, faceNormalOpposite));
                    metrics.push_back(UserType(angle));
                }
                else
                {
                    metrics.push_back(UserType(0.));
                }
            }

            // Compute per vertex metric
            vertexMetrics.clear();
            vertexMetrics.resize(sm.number_of_vertices(), static_cast<UserType>(0.));


            for (const auto &vertex: sm.vertices()) {
                UserType vertexMetric = UserType(0.);
                size_t total = 0;
                auto id = vertex.idx();
                for (const auto &h: CGAL::halfedges_around_source(vertex, sm)) {
                    auto idx = h.idx();
                    if(idx > metrics.size()) continue;
                    vertexMetric =  std::max(vertexMetric, metrics[idx]);
                    total++;
                }
                vertexMetrics[id] = vertexMetric;
            }
            return vertexMetrics;
        }
        std::vector<UserType> vertexMetrics;
    };
};