#pragma once
#include "type.h"
#include <queue>
//#include "optimizer.h"
#include "surface_optimizer.h"
// AABB
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/centroid.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABB_traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<AABB_traits>;
namespace sharp_feature
{
    class GraphOptimizer
    {
    public:
        GraphOptimizer()
        {

        }
        void optimizePolyLines()
        {

        }
        void computePolyLines(const std::vector<Point>& points, const std::vector<std::vector<int>>& featureGraph) {
            size_t featurePointCount = points.size();
            corners.clear();
            polylines.points = points;
            // compute Degrees
            for (size_t i = 0; i < featurePointCount; i++) {
                size_t degree = featureGraph[i].size();
                PolyLineAttribute attribute;
                attribute.degree = degree;
                polylines.attributes.push_back(attribute);
                if (degree != 2) {
                    corners.insert(i);
                }
            }
            struct EdgeKey {
                size_t a;
                size_t b;

                EdgeKey(size_t v0, size_t v1) {
                    a = std::min(v0, v1);
                    b = std::max(v0, v1);
                }

                bool operator==(const EdgeKey &) const = default;
            };
            struct EdgeKeyHash {
                size_t operator()(const EdgeKey &e) const {
                    return std::hash < size_t > {}(e.a)
                           ^ (std::hash < size_t > {}(e.b) << 1);
                }
            };

            std::unordered_set <EdgeKey, EdgeKeyHash> visited;

            std::vector <std::vector<size_t>> polylines_lines;

            // computePolyLines
            for (const auto &cornerID: corners) {
                for (auto &neighbor: featureGraph[cornerID]) {
                    std::vector <size_t> polyline;
                    size_t start = cornerID;
                    polyline.push_back(start);

                    size_t prev = start;
                    size_t current = neighbor;
                    EdgeKey edgeKey(start, static_cast<size_t>(neighbor));
                    if (visited.count(edgeKey)) {
                        continue;
                    }
                    visited.insert(edgeKey);

                    while (true) {
                        polyline.push_back(current);

                        if (polylines.attributes[current].degree != 2)
                            break;

                        size_t next = SIZE_MAX;
                        for (auto n: featureGraph[current]) {
                            if (n != prev) {
                                next = n;
                                break;
                            }
                        }
                        visited.insert({current, next});

                        prev = current;
                        current = next;
                    }
                    polylines_lines.push_back(polyline);
                }
            }
            polylines.lines = polylines_lines;
        }
        // output
        [[nodiscard]] Polylines<Point> getPolyLines() const
        {
            return polylines;
        }
        [[nodiscard]] std::vector<Point> getCorners() const
        {
            std::vector<Point> cornerPoints;
            for(const auto& corner : corners) cornerPoints.push_back(polylines.points[corner]);
            return cornerPoints;
        }
        void computeFrontierLines(const Mesh& surfaceMesh)
        {
            Tree tree(faces(surfaceMesh).first, faces(surfaceMesh).second, surfaceMesh);
            auto getClosestPoint = [&](size_t queryID) -> int
            {
                // Get the query point
                const auto& queryPt = polylines.points[queryID];

                // Find closest point on the mesh and the face it lies on
                auto closest = tree.closest_point_and_primitive(queryPt);
                auto faceDesc = closest.second;

                FT min_length = FT(FLT_MAX);
                int closestID = -1;

                // Iterate over the vertices of the closest face
                for(auto vd : vertices_around_face(surfaceMesh.halfedge(faceDesc), surfaceMesh))
                {
                    int vid = static_cast<int>(vd);  // vertex index
                    K::Vector_3 diff(queryPt, surfaceMesh.point(vd));
                    FT sq_dist = diff.squared_length();


                    if(sq_dist < min_length)
                    {
                        min_length = sq_dist;
                        closestID = vid;
                    }
                }

                return closestID;
            };

            //vertexLabel.resize(surfaceMesh.num_vertices());
            vertexDegree.resize(surfaceMesh.num_vertices(), 0);
            std::vector<std::vector<Mesh::Vertex_index>> candidates;
            std::unordered_map<int, int> cornerMap;

            // compute candidate sequences of closest points
            for(auto& line : polylines.lines) {
                std::vector<Mesh::Vertex_index> sequence;
                for (int i = 0; i < line.size() ; i++) {
                    int firstID = line[i];
                    int closest = getClosestPoint(firstID);
                    sequence.push_back(Mesh::Vertex_index (closest));
                    polylines.attributes[line[i]].vertexID = closest;
                    vertexDegree[closest] = std::max(vertexDegree[closest] , static_cast<int>(polylines.attributes[line[i]].degree));
                }
                auto new_end = std::unique(sequence.begin(), sequence.end());
                sequence.erase(new_end, sequence.end());
                candidates.push_back(sequence);
            }

            //
            std::vector<std::vector<Mesh::Halfedge_index >> segmentHEALL;
            computeDijkstra(surfaceMesh, candidates, segmentHEALL);

            edgeFrontiers.resize(surfaceMesh.number_of_edges(), 0.0);
            halfEdgeFrontier.resize(surfaceMesh.number_of_halfedges(), 0);
            for (auto vec : segmentHEALL) {
                for (auto h : vec) {
                    halfEdgeFrontier[h.idx()] = 1;
                    auto e = surfaceMesh.edge(h);
                    //edgeFrontiers[e.idx()] = 1.0;
                    edgeFrontiers[e.idx()] = edgeToPolyline[e]+1;
                }
            }

        }
        void computeFacePatches(const Mesh& surfaceMesh)
        {
            computeFrontierLines(surfaceMesh);
            computeKMeansPatches(surfaceMesh);
            computeOptimalPoint(surfaceMesh);
        }
        void noisePolyLines(const Mesh& surfaceMesh) {
            Tree tree(faces(surfaceMesh).first, faces(surfaceMesh).second, surfaceMesh);
            auto getClosestPointFace = [&](const Point &queryPt, Mesh::Face_index &faceDesc) -> int {
                // Get the query point

                // Find closest point on the mesh and the face it lies on
                auto closest = tree.closest_point_and_primitive(queryPt);
                faceDesc = closest.second;

                FT min_length = FT(FLT_MAX);
                int closestID = -1;

                // Iterate over the vertices of the closest face
                for (auto vd: vertices_around_face(surfaceMesh.halfedge(faceDesc), surfaceMesh)) {
                    int vid = static_cast<int>(vd);  // vertex index
                    K::Vector_3 diff(queryPt, surfaceMesh.point(vd));
                    FT sq_dist = diff.squared_length();

                    if (sq_dist < min_length) {
                        min_length = sq_dist;
                        closestID = vid;
                    }
                }

                return closestID;
            };
            std::vector<Mesh::Face_index> supportingFaces;
            for (auto &point: polylines.points) {
                Mesh::Face_index supportingFace;
                int closest = getClosestPointFace(point, supportingFace);
                supportingFaces.push_back(supportingFace);
            }
            auto computeCentroid = [&](Mesh::Face_index startFace) {
                auto h = surfaceMesh.halfedge(startFace);

                Point p0 = surfaceMesh.point(source(h, surfaceMesh));
                Point p1 = surfaceMesh.point(target(h, surfaceMesh));
                Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

                Point p = CGAL::centroid(p0, p1, p2);
                return p;
            };
            for(auto &l : polylines.lines)
            {
                for(int i = 0 ; i < l.size()-1;i++)
                {
                    auto& point = polylines.points[l[i]];
                    auto point2 = polylines.points[l[i+1]];
                    auto n_ = point2-point;
                    Vector n = n_ / std::sqrt(n_.squared_length());
                    //point += 0.25*n;
                    point = computeCentroid(Mesh::Face_index (supportingFaces[l[i]]));
                }
            }

        }
        void optimizePoints(const Mesh& surfaceMesh) {
            Tree tree(faces(surfaceMesh).first, faces(surfaceMesh).second, surfaceMesh);
            auto getClosestPointFace = [&](const Point &queryPt, Mesh::Face_index &faceDesc) -> int {
                // Get the query point

                // Find closest point on the mesh and the face it lies on
                auto closest = tree.closest_point_and_primitive(queryPt);
                faceDesc = closest.second;

                FT min_length = FT(FLT_MAX);
                int closestID = -1;

                // Iterate over the vertices of the closest face
                for (auto vd: vertices_around_face(surfaceMesh.halfedge(faceDesc), surfaceMesh)) {
                    int vid = static_cast<int>(vd);  // vertex index
                    K::Vector_3 diff(queryPt, surfaceMesh.point(vd));
                    FT sq_dist = diff.squared_length();

                    if (sq_dist < min_length) {
                        min_length = sq_dist;
                        closestID = vid;
                    }
                }

                return closestID;
            };
            std::vector<Mesh::Face_index> supportingFaces;
            for (auto &point: polylines.points) {
                Mesh::Face_index supportingFace;
                int closest = getClosestPointFace(point, supportingFace);
                supportingFaces.push_back(supportingFace);
            }
            float deltaXY = 0.25f;
            int steps = 1;
            optimizer = SurfaceOptimizer(surfaceMesh);
            optimizer.parameters.deltaXY = deltaXY;
            optimizer.max_iteration = maxIteration;
            size_t id = 0;

            for(auto& point : polylines.points)
            {
                //std::cout<<"Point ID "<< id <<"\n";
                if(id != 123 )
                {
                    id++;
                    continue;
                }

                auto startFace = supportingFaces[id];
                optimizer.optimize(point, optimalPoints[id], startFace);

                id++;
            }

        }
        void optimizePoints_old(const Mesh& surfaceMesh)
        {
            Tree tree(faces(surfaceMesh).first, faces(surfaceMesh).second, surfaceMesh);
            auto getClosestPointFace = [&](const Point& queryPt, Mesh::Face_index& faceDesc) -> int
            {
                // Get the query point

                // Find closest point on the mesh and the face it lies on
                auto closest = tree.closest_point_and_primitive(queryPt);
                faceDesc = closest.second;

                FT min_length = FT(FLT_MAX);
                int closestID = -1;

                // Iterate over the vertices of the closest face
                for(auto vd : vertices_around_face(surfaceMesh.halfedge(faceDesc), surfaceMesh))
                {
                    int vid = static_cast<int>(vd);  // vertex index
                    K::Vector_3 diff(queryPt, surfaceMesh.point(vd));
                    FT sq_dist = diff.squared_length();

                    if(sq_dist < min_length)
                    {
                        min_length = sq_dist;
                        closestID = vid;
                    }
                }

                return closestID;
            };
            std::vector<Mesh::Face_index> supportingFaces;
            for(auto& point : polylines.points)
            {
                Mesh::Face_index supportingFace;
                int closest = getClosestPointFace(point, supportingFace);
                supportingFaces.push_back(supportingFace);
            }
            float deltaXY = 0.1f;
            int steps = 1;
            optimizer = SurfaceOptimizer(surfaceMesh);
            optimizer.parameters.deltaXY = deltaXY;
            optimizer.max_iteration = steps;
            size_t id = 0;
            for(auto& point : polylines.points)
            {


                auto point2 = id == 0 ? polylines.points[1] : polylines.points[id -1];
                auto n_ = point-point2;
                Vector n = n_ / std::sqrt(n_.squared_length());
                Vector t;
                if (std::abs(n.x()) < 0.9)
                    t = Vector(1, 0, 0);
                else
                    t = Vector(0, 1, 0);

                // First tangent
                Vector u = CGAL::cross_product(n, t);
                u = u / std::sqrt(u.squared_length());

                // Second tangent
                Vector v = CGAL::cross_product(n, u);
                point = point + 0.1*u;

                /*std::vector<Point> debugOpti;
                point = surfaceMesh.point(Mesh::Vertex_index (30));
                debugOpti.push_back(point);
                supportingFaces[id] = Mesh::Face_index (55);*/
                auto startFace = supportingFaces[id];
                //optimizer.optimizeUV(point, optimalPoints[id], startFace);
                /*debugOpti.push_back(optimalPoints[id]);
                debugOpti.push_back(point);
                //processed[id] = 1;
                auto pcl = polyscope::registerPointCloud("debugOpti", debugOpti);
                pcl->resetTransform();

                {
                    auto pcl = polyscope::registerPointCloud("gradients", optimizer.gradients);
                    pcl->resetTransform();
                }
                {
                    auto pcl = polyscope::registerPointCloud("debugPath", optimizer.debugPath);
                    pcl->resetTransform();
                }*/
                std::cout<<"Optimization id "<<id<<"\n";





                id++;
            }
        }
        int maxIteration=2;
        std::vector<std::vector<int>> labelFaces;
        std::vector<std::vector<Mesh::Face_index >> faceNeighbors;
        std::vector<int> vertexDegree;
        std::vector<Point> optimalPoints;
        std::vector<std::vector<Point>> planeCentroids;
        std::vector<std::vector<Vector >> planeNormals;
        std::vector<double> edgeFrontiers;
    private:
        [[nodiscard]] static Point computeCentroid(const Mesh& sm, const CGAL::SM_Face_index& face)
        {
            Vector centroid = CGAL::NULL_VECTOR;
            size_t total      = 0;
            for (const auto& vertexHandle : sm.vertices_around_face(sm.halfedge(face)))
            {
                centroid += (sm.point(vertexHandle) - CGAL::ORIGIN);
                total++;
            }
            return CGAL::ORIGIN + (centroid / static_cast<double>(total));
        }
        void computeOptimalPoint(const Mesh& sm)
        {
            planeCentroids.resize( polylines.points.size());
            planeNormals.resize( polylines.points.size());

            auto computeQuadricFromPlane = [](const Plane& plane)
            {
                const auto& [p, n_raw] = plane;

                // normalize normal
                Vector n = n_raw / std::sqrt(n_raw.squared_length());

                double a = n.x();
                double b = n.y();
                double c = n.z();

                double d = -(a * p.x() + b * p.y() + c * p.z());

                Eigen::Vector4d v(a, b, c, d);

                Eigen::Matrix4d Q = v * v.transpose(); // outer product

                return Q;
            };
            auto optimizeQuadric = [](const Eigen::Matrix4d& Q)
            {
                // Extract A (top-left 3x3) and b (top-right 3x1)
                Eigen::Matrix3d A = Q.block<3,3>(0,0);
                Eigen::Vector3d b = Q.block<3,1>(0,3);

                // We solve: A * x = -b
                Eigen::Vector3d x;

                // Check if A is invertible
                if (std::abs(A.determinant()) > 1e-10)
                {
                    x = -A.inverse() * b;
                }
                else
                {
                    // fallback: use pseudo-inverse or just return origin-ish
                    x = -A.completeOrthogonalDecomposition().solve(b);
                }

                return x; // Eigen::Vector3d
            };
            auto computeLineQuadric= [](const Eigen::Vector3d& vi,
                                               const Eigen::Vector3d& n)
            {
                Eigen::Vector3d ni = n.normalized();
                // RNG
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<double> dis(0.0, 1.0);

                // random vector
                Eigen::Vector3d xi(dis(gen), dis(gen), dis(gen));

                // project xi onto plane orthogonal to ni
                xi = xi - (xi.dot(ni) / ni.dot(ni)) * ni;
                xi.normalize();

                // second orthogonal direction
                Eigen::Vector3d yi = xi.cross(ni);


                // build plane equations
                Eigen::Vector4d px(xi.x(), xi.y(), xi.z(), -xi.dot(vi));
                Eigen::Vector4d py(yi.x(), yi.y(), yi.z(), -yi.dot(vi));

                // quadrics
                Eigen::Matrix4d Qx = px * px.transpose();
                Eigen::Matrix4d Qy = py * py.transpose();

                auto q =  Qx + Qy;

                return q.eval();
            };
            // for each point of the line compute optimal point
            for(int i = 0 ; i < polylines.points.size(); i++)
            {
                auto closestPoint = polylines.attributes[i].vertexID;
                std::unordered_map<int, Vector > normalPlanes;


                // Plane centroids
                std::unordered_map<int, Vector > positionPlanes;
                std::unordered_map<int, size_t > counts;
                for(int j = 0 ; j < labelFaces[closestPoint].size(); j++)
                {
                    auto face = faceNeighbors[closestPoint][j];
                    auto faceCentroid =  computeCentroid(sm, face);
                    auto normal = faceNormals[face.idx()];
                    auto planeID = labelFaces[closestPoint][j];

                    if(positionPlanes.contains(planeID))
                    {
                        positionPlanes[planeID] += faceCentroid - CGAL::ORIGIN;
                        normalPlanes[planeID] += normal;
                        counts[planeID] ++;
                    }
                    else
                    {
                        positionPlanes[planeID] = faceCentroid - CGAL::ORIGIN;
                        normalPlanes[planeID] = normal;
                        counts[planeID] = 1;
                    }


                }

                for(auto& position : positionPlanes)
                {
                    position.second /= static_cast<double>(counts[position.first]);
                }
                for(auto& normal : normalPlanes)
                {
                    normal.second /= static_cast<double>(counts[normal.first]);
                }

                Eigen::Matrix4d quadric = Eigen::Matrix4d::Zero();
                for(int planeID = 0 ; planeID < normalPlanes.size(); planeID++)
                {
                    Plane plane = std::make_pair(CGAL::ORIGIN+positionPlanes[planeID], normalPlanes[planeID]);
                    planeCentroids[i].push_back(CGAL::ORIGIN+positionPlanes[planeID]);
                    planeNormals[i].push_back(normalPlanes[planeID]);
                    quadric += computeQuadricFromPlane(plane);
                }

                auto inputPoint = sm.point(Mesh::Vertex_index (closestPoint));
                Vector inputNormal = PMP::compute_vertex_normal(Mesh::Vertex_index (closestPoint), sm);

                Eigen::Vector3d point(inputPoint.x(),inputPoint.y(),inputPoint.z());
                Eigen::Vector3d normal(inputNormal.x(),inputNormal.y(),inputNormal.z());
                Eigen::Matrix4d lq = computeLineQuadric(point, normal);

                quadric += 1e-1 * lq;


                auto eigenPoint = optimizeQuadric(quadric);
                auto cgalPoint =  Point (eigenPoint[0], eigenPoint[1], eigenPoint[2]);
                optimalPoints.push_back(cgalPoint);

            }

        }
        void computeNeighborhood(const Mesh& sm) {

            maxRadius *= maxRadius;
            faceNeighbors.resize(sm.num_vertices());
            vertexNeighbors.resize(sm.num_vertices());
            faceDepth.resize(sm.number_of_faces(), -1);

            if(faceNormals.empty())
            {
                for(const auto& faceID : sm.faces())
                {
                    faceNormals.push_back( CGAL::Polygon_mesh_processing::compute_face_normal(faceID, sm));
                }
            }
            struct ElementQ {
                Mesh::Face_index index;
                int depth = 0;
            };
            auto checkFaceCondition = [&](float maxR, Mesh::Face_index faceIndex, Mesh::Vertex_index vertexIndex) {
                int cond = 0;
                for (auto v: CGAL::vertices_around_face(sm.halfedge(faceIndex), sm)) {
                    const Point &q = sm.point(v);
                    auto sqdist = CGAL::squared_distance(q, sm.point(vertexIndex));
                    if (sqdist > maxR) {
                        cond++;
                    }
                }
                if (cond > 2) return false;
                return true;
            };
            for (const auto &vID: sm.vertices()) {
                std::unordered_set<size_t > allowedBoundary;
                std::queue<ElementQ> q;
                std::unordered_set<Mesh::Face_index> visited;
                for (const auto &h: CGAL::halfedges_around_target(sm.halfedge(vID), sm)) {

                    auto e = sm.edge(h);
                    if(edgeToPolyline.contains(e))
                        allowedBoundary.insert(edgeToPolyline.at(e));
                    const auto &vn = sm.source(h);
                    const auto &fn = sm.face(h);
                    if(fn == Mesh::null_face()) continue;
                    if (!checkFaceCondition(maxRadius, fn, vID)) continue;
                    visited.insert(fn);
                    q.push({fn, 0});
                    faceNeighbors[vID].push_back(fn);
                    faceDepth[fn] = 0;
                }

                while (!q.empty()) {
                    const auto &current = q.front();
                    q.pop();
                    if (current.depth == maxDepth) continue;

                    for (const auto &h: CGAL::halfedges_around_face(sm.halfedge(current.index), sm)) {
                        auto fn = sm.face(sm.opposite(h));
                        if(fn == Mesh::null_face()) continue;
                        if (visited.contains(fn))continue;
                        if (!checkFaceCondition(maxRadius, fn, vID)) continue;
                        auto e = sm.edge(h);
                        if(edgeToPolyline.contains(e))
                            if(!allowedBoundary.contains(edgeToPolyline.at(e))) continue;
                        visited.insert(fn);
                        q.push({fn, current.depth + 1});
                        faceNeighbors[vID].push_back(fn);
                        faceDepth[fn] = current.depth + 1;
                    }

                }
            }
        }
        void computePatches(const Mesh& surfaceMesh)
        {

        }
        void computeKMeansPatches(const Mesh& surfaceMesh)
        {
            int nseed = 3;
            if(faceNeighbors.empty())
            {
                computeNeighborhood(surfaceMesh);
            }
            labelFaces.resize(surfaceMesh.num_vertices());
            centroidsVec.resize(surfaceMesh.num_vertices());
            for(const auto& vID : surfaceMesh.vertices())
            {
                /*if(corners.contains(vID))
                {
                    nseed = 3;
                }
                else
                {
                    nseed = 2;
                }*/


                nseed = vertexDegree[vID];
                if(nseed == 0) continue;
                std::vector<Vector > neighborNormals;
                for(const auto& neighbor: faceNeighbors[vID.idx()])
                {
                    neighborNormals.push_back(faceNormals[neighbor.idx()]);

                }
                std::vector<Vector> centroids;
                std::vector<int> idx(neighborNormals.size());
                std::iota(idx.begin(), idx.end(), 0);

                std::random_device rd;
                std::mt19937 gen(rd());
                std::shuffle(idx.begin(), idx.end(), gen);

                for (int k = 0; k < nseed; k++)
                {
                    centroids.push_back(neighborNormals[idx[k]]);
                }

                std::vector<int> labels(neighborNormals.size(), 0);

                for (int iter = 0; iter < 5; iter++)
                {
                    // reset accumulators
                    std::vector<Vector> sum(nseed, Vector(0,0,0));
                    std::vector<int> count(nseed, 0);

                    // ASSIGN STEP
                    for (int i = 0; i < neighborNormals.size(); i++)
                    {
                        int best = 0;
                        double bestDist = std::numeric_limits<double>::max();

                        for (int k = 0; k < nseed; k++)
                        {
                            //double d = (neighborNormals[i] - centroids[k]).squared_length();
                            double d = 1.0 - neighborNormals[i]* centroids[k];
                            if (d < bestDist)
                            {
                                bestDist = d;
                                best = k;
                            }
                        }

                        labels[i] = best;
                        sum[best] += neighborNormals[i];
                        count[best]++;
                    }

                    // UPDATE STEP
                    for (int k = 0; k < nseed; k++)
                    {
                        if (count[k] > 0)
                        {
                            centroids[k] = sum[k] / count[k];
                        }
                    }
                    for (auto& c : centroids)
                    {
                        c = c / std::sqrt(c.squared_length());
                    }
                }
                labelFaces[vID.idx()] = labels;
                centroidsVec[vID.idx()] = centroids;
            }
        }
        void computeDijkstra(const Mesh& surfaceMesh, const std::vector<std::vector<Mesh::Vertex_index>>& candidates,  std::vector<std::vector<Mesh::Halfedge_index >>& segmentHEALL)
        {
            // Priority queue item
            struct Node
            {
                Mesh::Vertex_index v;
                double dist;
                K::Vector_3 dir;

                bool operator>(const Node& other) const
                {
                    return dist > other.dist;
                }
            };
            std::vector<std::vector<Mesh::Vertex_index>> all_paths;

            size_t idLine = 0;
            for (const auto& seq : candidates)
            {
                std::vector<Mesh::Vertex_index> full_path;
                std::vector<Mesh::Halfedge_index> full_path_HE;

                if (seq.size() < 2)
                {
                    all_paths.push_back(full_path);
                    continue;
                }

                for (size_t i = 0; i < seq.size() - 1; ++i)
                {
                    Mesh::Vertex_index source = seq[i];
                    Mesh::Vertex_index target = seq[i + 1];

                    std::unordered_map<Mesh::Vertex_index, double> dist;
                    std::unordered_map<Mesh::Vertex_index, Mesh::Vertex_index> pred;

                    std::unordered_map<Mesh::Vertex_index, Mesh::Halfedge_index > predHE;

                    // --- Dijkstra with early stop ---
                    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> pq;

                    for (Mesh::Vertex_index v : surfaceMesh.vertices())
                        dist[v] = std::numeric_limits<double>::infinity();

                    dist[source] = 0.0;
                    pq.push({source, 0.0});

                    bool found = false;

                    while (!pq.empty())
                    {
                        Node current = pq.top();
                        pq.pop();

                        Mesh::Vertex_index u = current.v;

                        if (current.dist > dist[u])
                            continue;

                        if (u == target)
                        {
                            found = true;
                            break;
                        }

                        Point pu = surfaceMesh.point(u);

                        for (Mesh::Halfedge_index h : CGAL::halfedges_around_target(surfaceMesh.halfedge(u), surfaceMesh))
                        {
                            Mesh::Vertex_index v = surfaceMesh.source(h);
                            Point pv = surfaceMesh.point(v);

                            double weight = std::sqrt(CGAL::squared_distance(pu, pv));
                            double new_dist = dist[u] + weight;

                            if (new_dist < dist[v])
                            {
                                dist[v] = new_dist;
                                pred[v] = u;
                                predHE[v] = h;
                                pq.push({v, new_dist});
                            }
                        }
                    }

                    if (!found)
                    {
                        std::cerr << "No path found in one sequence\n";
                        continue;
                    }

                    // --- Reconstruct segment ---
                    std::vector<Mesh::Vertex_index> segment;
                    std::vector<Mesh::Halfedge_index> segmentHE;

                    Mesh::Vertex_index current = target;
                    Mesh::Halfedge_index currentHE;
                    while (current != source)
                    {
                        segment.push_back(current);
                        currentHE = predHE[current];
                        current = pred[current];

                        segmentHE.push_back(currentHE);
                        edgeToPolyline[surfaceMesh.edge(currentHE)] = idLine;
                    }

                    segment.push_back(source);
                    std::reverse(segment.begin(), segment.end());

                    std::reverse(segmentHE.begin(), segmentHE.end());

                    // Avoid duplication
                    if (!full_path.empty())
                        segment.erase(segment.begin());

                    full_path.insert(full_path.end(), segment.begin(), segment.end());

                    full_path_HE.insert(full_path_HE.end(), segmentHE.begin(), segmentHE.end());
                }
                segmentHEALL.push_back(full_path_HE);
                all_paths.push_back(full_path);
                idLine++;
            }
            //return all_paths;
        }
        std::unordered_set <size_t> corners;
        Polylines<Point> polylines;
        std::unordered_map<int, std::unordered_map<int , std::unordered_set<int >>> cornerPatches;
        std::vector<int> halfEdgeFrontier;

        std::unordered_map<Mesh::Edge_index , size_t> edgeToPolyline;



        std::vector<std::vector<Mesh::Vertex_index>> vertexNeighbors;
        std::vector<int> faceDepth;
        std::vector<Vector > faceNormals;


        float maxRadius = 3.;
        int maxDepth = 3;
        std::vector<std::vector<Vector>> centroidsVec;
        std::vector<std::vector<Vector>> centroidsPos;

        SurfaceOptimizer optimizer;

    };
};