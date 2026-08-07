#pragma once
#include "graph_optimizer.h"
#include "surface_optimizer.h"
#include "coord.h"
namespace test_graph_optimizer
{
    void test1(Optimizer optimizer, const Mesh& surfaceMesh)
    {
        // cube corner
        // Vertex 0 corner

        // Vertex 2 link to 0
        // Vertex 6 link to 0
        // Vertex 1 link to 0
        // Vertex 3 link to 0

        // face 4 top
        // face 5 left
        // face 2 left
        // face 3 right

        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point& p, int id)
        {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]()
        {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        optimizer.parameters.deltaXY = 0.075;
        optimizer.max_iteration = 10;
        size_t id = 0;

        auto startFace =  Mesh::Face_index(5);
        auto v0 = Mesh::Vertex_index (0);
        auto v1 = Mesh::Vertex_index (2);
        auto p =
                CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1));
        auto optimalPoint = surfaceMesh.point(v0);

        add(p,0);
        add(optimalPoint,1);
        optimizer.optimize(p, optimalPoint, startFace);
        directions = optimizer.directions;
        add(p,2);
        display();
    }
    void test2(Optimizer optimizer, const Mesh& surfaceMesh)
    {
        // cube corner
        // Vertex 0 corner

        // Vertex 2 link to 0
        // Vertex 6 link to 0
        // Vertex 1 link to 0
        // Vertex 3 link to 0

        // face 4 top
        // face 5 left
        // face 2 left
        // face 3 right

        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point& p, int id)
        {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]()
        {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        optimizer.parameters.deltaXY = 0.2;
        optimizer.max_iteration = 40;
        size_t id = 0;

        auto startFace =  Mesh::Face_index(5);
        auto v0 = Mesh::Vertex_index (0);
        auto v1 = Mesh::Vertex_index (2);
        /*auto p =
                CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1));*/
        auto h = surfaceMesh.halfedge(startFace);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

        Point p = CGAL::centroid(p0, p1, p2);
        auto optimalPoint = surfaceMesh.point(v0);

        add(p,0);
        add(optimalPoint,1);
        optimizer.optimize(p, optimalPoint, startFace);
        directions = optimizer.directions;
        add(p,2);
        display();
    }
    void test3(Optimizer optimizer, const Mesh& surfaceMesh)
    {
        // cube corner
        // Vertex 0 corner

        // Vertex 2 link to 0
        // Vertex 6 link to 0
        // Vertex 1 link to 0
        // Vertex 3 link to 0

        // face 4 top
        // face 5 left
        // face 2 left
        // face 3 right

        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point& p, int id)
        {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]()
        {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        optimizer.parameters.deltaXY = 0.1;
        optimizer.max_iteration = 40;
        size_t id = 0;

        auto startFace =  Mesh::Face_index(0);
        auto v0 = Mesh::Vertex_index (0);
        auto v1 = Mesh::Vertex_index (2);
        /*auto p =
                CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1));*/
        auto h = surfaceMesh.halfedge(startFace);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

        Point p = CGAL::centroid(p0, p1, p2);
        auto optimalPoint = surfaceMesh.point(v0);

        add(p,0);
        add(optimalPoint,1);
        optimizer.optimize(p, optimalPoint, startFace);
        directions = optimizer.directions;
        add(p,2);
        display();
    }


    void testIntersectionCircle(const Mesh& surfaceMesh) {
        SurfaceOptimizer surfaceOptimizer(surfaceMesh);
        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point &p, int id) {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]() {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };

        auto startFace = Mesh::Face_index(5);
        auto centroid = computeCentroid(startFace);
        auto v0 = Mesh::Vertex_index(7);
        auto v1 = Mesh::Vertex_index(3);
        auto v2 = Mesh::Vertex_index(5);

        auto p0 = surfaceMesh.point(v0);
        auto p1 = surfaceMesh.point(v1);
        auto p2 = surfaceMesh.point(v2);

        auto angleBetweenVectors = [&](Vector v1, Vector v2) -> double
        {
            double cosTheta =
                    CGAL::to_double(v1 * v2) /
                    std::sqrt(CGAL::to_double(v1.squared_length() * v2.squared_length()));

            cosTheta = std::clamp(cosTheta, -1.0, 1.0);

            return std::acos(cosTheta); // radians
        };

        std::vector<double> angleBounds;
        angleBounds.push_back(0.0);

        double totalAngle = 0.0;

        std::vector<std::pair<Vector, Vector>> vectors;
        auto addAngle = [&](const Point& a, const Point& b)
        {
            double angle = angleBetweenVectors(a - centroid, b - centroid);
            vectors.emplace_back(a - centroid, b - centroid);
            totalAngle += angle;
            angleBounds.push_back(totalAngle);
        };
        auto normalize = [](const Vector& x) {
            return x / std::sqrt(CGAL::to_double(x.squared_length()));
        };
        auto rotateTowards = [&](const Vector& u, const Vector& axis, double alpha)
        {
            return std::cos(alpha) * u
                   + std::sin(alpha) * CGAL::cross_product(axis, u)
                   + (1.0 - std::cos(alpha)) * (axis * u) * axis;
        };

        addAngle(p0, p1);
        addAngle(p1, p2);
        addAngle(p2, p0);

        size_t nSamples = 256;
        float r = 1.;

        size_t lineID = 0;
        auto traceLine = [&] (Point p1, Point p2)
        {
            std::vector<Point> pts = {p1 ,p2};
            std::vector<std::array<size_t, 2>> line = {{0 ,1}};
            auto c = polyscope::registerCurveNetwork("line_"+std::to_string(lineID), pts, line);
            c->resetTransform();
            lineID++;
        };

        for (size_t i = 0; i < nSamples; ++i) {
            // theta = angle of i * 2pi / nSample
            double theta = totalAngle * double(i) / double(nSamples);

            auto it = std::upper_bound(
                    angleBounds.begin(),
                    angleBounds.end(),
                    theta
            );

            if (it == angleBounds.begin() || it == angleBounds.end())
                continue;

            int k = static_cast<int>(std::distance(angleBounds.begin(), it)) - 1;
            auto [u, v] = vectors[k];

            Vector uhat = normalize(u);

            Vector n = CGAL::cross_product(u, v);
            n = normalize(n);

            // local angle inside this sector
            double alpha = theta - angleBounds[k];

            // rotate u toward v around triangle normal
            Vector w = rotateTowards(uhat, n, alpha);

            w = normalize(w);

            auto point =  centroid + r * w;
            // debug visualization
            directions.push_back( point);

            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh, computeCentroid(startFace), startFace);
            auto targetPoint = computeUV(surfaceMesh, point, startFace);
            surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);
            auto endPoint = computeXYZ(surfaceMesh, end, startFace);
            add(endPoint, 1);

        }
        display();

    }
    void testOptimizationSurfaceCircle(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer surfaceOptimizer(surfaceMesh);
        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point &p, int id) {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]() {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };

        auto startFace = Mesh::Face_index(5);
        auto centroid = computeCentroid(startFace);
        auto v0 = Mesh::Vertex_index(7);
        auto v1 = Mesh::Vertex_index(3);
        auto v2 = Mesh::Vertex_index(5);

        auto p0 = surfaceMesh.point(v0);
        auto p1 = surfaceMesh.point(v1);
        auto p2 = surfaceMesh.point(v2);

        auto angleBetweenVectors = [&](Vector v1, Vector v2) -> double
        {
            double cosTheta =
                    CGAL::to_double(v1 * v2) /
                    std::sqrt(CGAL::to_double(v1.squared_length() * v2.squared_length()));

            cosTheta = std::clamp(cosTheta, -1.0, 1.0);

            return std::acos(cosTheta); // radians
        };

        std::vector<double> angleBounds;
        angleBounds.push_back(0.0);

        double totalAngle = 0.0;

        std::vector<std::pair<Vector, Vector>> vectors;
        auto addAngle = [&](const Point& a, const Point& b)
        {
            double angle = angleBetweenVectors(a - centroid, b - centroid);
            vectors.emplace_back(a - centroid, b - centroid);
            totalAngle += angle;
            angleBounds.push_back(totalAngle);
        };
        auto normalize = [](const Vector& x) {
            return x / std::sqrt(CGAL::to_double(x.squared_length()));
        };
        auto rotateTowards = [&](const Vector& u, const Vector& axis, double alpha)
        {
            return std::cos(alpha) * u
                   + std::sin(alpha) * CGAL::cross_product(axis, u)
                   + (1.0 - std::cos(alpha)) * (axis * u) * axis;
        };

        addAngle(p0, p1);
        addAngle(p1, p2);
        addAngle(p2, p0);

        size_t nSamples = 256;
        float r = 0.5;

        size_t lineID = 0;
        auto traceLine = [&] (Point p1, Point p2)
        {
            std::vector<Point> pts = {p1 ,p2};
            std::vector<std::array<size_t, 2>> line = {{0 ,1}};
            auto c = polyscope::registerCurveNetwork("line_"+std::to_string(lineID), pts, line);
            c->resetTransform();
            lineID++;
        };

        for (size_t i = 0; i < nSamples; ++i) {
            // theta = angle of i * 2pi / nSample
            double theta = totalAngle * double(i) / double(nSamples);

            auto it = std::upper_bound(
                    angleBounds.begin(),
                    angleBounds.end(),
                    theta
            );

            if (it == angleBounds.begin() || it == angleBounds.end())
                continue;

            int k = static_cast<int>(std::distance(angleBounds.begin(), it)) - 1;
            auto [u, v] = vectors[k];

            Vector uhat = normalize(u);

            Vector n = CGAL::cross_product(u, v);
            n = normalize(n);

            // local angle inside this sector
            double alpha = theta - angleBounds[k];

            // rotate u toward v around triangle normal
            Vector w = rotateTowards(uhat, n, alpha);

            w = normalize(w);

            auto point =  centroid + r * w;
            // debug visualization
            directions.push_back( point);

            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh, computeCentroid(startFace), startFace);
            auto targetPoint = computeUV(surfaceMesh, point, startFace);
            /*surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);*/

            float totalDistance = 0.5;
            Point pointMoved;
            surfaceOptimizer.moveOnSurface(startFace,
                                           startPoint,
                                           targetPoint-startPoint,
                                    pointMoved,
                                    faceEnd);

            //auto endPoint = computeXYZ(surfaceMesh, end, startFace);
            add(pointMoved, 1);

        }
        display();
    }
    void testMoveSurfaceCircle(const Mesh& surfaceMesh) {
        SurfaceOptimizer surfaceOptimizer(surfaceMesh);
        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point &p, int id) {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]() {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };

        auto startFace = Mesh::Face_index(4);
        auto centroid = computeCentroid(startFace);
        auto v0 = Mesh::Vertex_index(2);
        auto v1 = Mesh::Vertex_index(0);
        auto v2 = Mesh::Vertex_index(3);

        auto p0 = surfaceMesh.point(v0);
        auto p1 = surfaceMesh.point(v1);
        auto p2 = surfaceMesh.point(v2);

        auto angleBetweenVectors = [&](Vector v1, Vector v2) -> double
        {
            double cosTheta =
                    CGAL::to_double(v1 * v2) /
                    std::sqrt(CGAL::to_double(v1.squared_length() * v2.squared_length()));

            cosTheta = std::clamp(cosTheta, -1.0, 1.0);

            return std::acos(cosTheta); // radians
        };

        std::vector<double> angleBounds;
        angleBounds.push_back(0.0);

        double totalAngle = 0.0;

        std::vector<std::pair<Vector, Vector>> vectors;
        auto addAngle = [&](const Point& a, const Point& b)
        {
            double angle = angleBetweenVectors(a - centroid, b - centroid);
            vectors.emplace_back(a - centroid, b - centroid);
            totalAngle += angle;
            angleBounds.push_back(totalAngle);
        };
        auto normalize = [](const Vector& x) {
            return x / std::sqrt(CGAL::to_double(x.squared_length()));
        };
        auto rotateTowards = [&](const Vector& u, const Vector& axis, double alpha)
        {
            return std::cos(alpha) * u
                   + std::sin(alpha) * CGAL::cross_product(axis, u)
                   + (1.0 - std::cos(alpha)) * (axis * u) * axis;
        };

        addAngle(p0, p1);
        addAngle(p1, p2);
        addAngle(p2, p0);

        size_t nSamples = 256;
        float r = 1.;

        size_t lineID = 0;
        auto traceLine = [&] (Point p1, Point p2)
        {
            std::vector<Point> pts = {p1 ,p2};
            std::vector<std::array<size_t, 2>> line = {{0 ,1}};
            auto c = polyscope::registerCurveNetwork("line_"+std::to_string(lineID), pts, line);
            c->resetTransform();
            lineID++;
        };

        for (size_t i = 0; i < nSamples; ++i) {
            // theta = angle of i * 2pi / nSample
            double theta = totalAngle * double(i) / double(nSamples);

            auto it = std::upper_bound(
                    angleBounds.begin(),
                    angleBounds.end(),
                    theta
            );

            if (it == angleBounds.begin() || it == angleBounds.end())
                continue;

            int k = static_cast<int>(std::distance(angleBounds.begin(), it)) - 1;
            auto [u, v] = vectors[k];

            Vector uhat = normalize(u);

            Vector n = CGAL::cross_product(u, v);
            n = normalize(n);

            // local angle inside this sector
            double alpha = theta - angleBounds[k];

            // rotate u toward v around triangle normal
            Vector w = rotateTowards(uhat, n, alpha);

            w = normalize(w);

            auto point =  centroid + r * w;
            // debug visualization
            directions.push_back( point);

            if(i > 79 && i < 84)
            {
                surfaceOptimizer.debug = true;
            }
            else
            {
                surfaceOptimizer.debug = false;
            }
            bool intersectiononly = false;
            if(intersectiononly)
            {

                Mesh::Face_index faceEnd;
                K::Point_2 end;
                auto startPoint = computeUV(surfaceMesh, computeCentroid(startFace), startFace);
                auto targetPoint = computeUV(surfaceMesh, point, startFace);
                surfaceOptimizer.intersection2D(startFace,
                                                startPoint,
                                                targetPoint,
                                                faceEnd, end);
                auto endPoint = computeXYZ(surfaceMesh, end, startFace);
                add(endPoint, 1);
            }
            else
            {
                Mesh::Face_index faceEnd;
                float totalDistance = 0.5;
                auto startPoint = computeCentroid(startFace);
                auto targetPoint = point;
                surfaceOptimizer.moveOnSurface(startFace,
                                               startPoint,
                                               point,
                                                faceEnd,
                                                totalDistance);
                //auto endPoint = computeXYZ(surfaceMesh, end, startFace);
                add(point, 1);

                if(surfaceOptimizer.debug)
                {
                    traceLine(startPoint, targetPoint);
                    add(point, 1);
                    auto pcl = polyscope::registerPointCloud("steps", surfaceOptimizer.debugPath);
                    pcl->addVectorQuantity("normals", surfaceOptimizer.debugNormal);
                    pcl->resetTransform();
                }
            }



        }
        display();

    }
    void testIntersection(const Mesh& surfaceMesh) {
        SurfaceOptimizer surfaceOptimizer(surfaceMesh);
        // cube corner
        // Vertex 0 corner

        // Vertex 2 link to 0
        // Vertex 6 link to 0
        // Vertex 1 link to 0
        // Vertex 3 link to 0

        // face 4 top
        // face 5 left
        // face 2 left
        // face 3 right

        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point &p, int id) {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]() {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };
        auto isEqual = [](const Point& p1, const Point& p2 )
        {
            double eps = 1e-12;

            return (CGAL::to_double(CGAL::squared_distance(p1, p2)) < eps);
        };
        size_t lineID = 0;
        auto traceLine = [&] (Point p1, Point p2)
        {
            std::vector<Point> pts = {p1 ,p2};
            std::vector<std::array<size_t, 2>> line = {{0 ,1}};
            auto c = polyscope::registerCurveNetwork("line_"+std::to_string(lineID), pts, line);
            c->resetTransform();
            lineID++;
        };

        size_t id = 0;

        auto startFace = Mesh::Face_index(4);

        auto f1 = Mesh::Face_index(1);
        auto f2 = Mesh::Face_index(5);
        auto f3 = Mesh::Face_index(3);

        auto v0 = Mesh::Vertex_index(0);
        auto v1 = Mesh::Vertex_index(2);
        auto v2 = Mesh::Vertex_index(3);
        auto v3 = Mesh::Vertex_index(1);


        {
            // centroid to point
            add( computeCentroid(startFace), 0);
            traceLine(computeCentroid(startFace), surfaceMesh.point(v0));
            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh, computeCentroid(startFace), startFace);
            auto targetPoint = computeUV(surfaceMesh, surfaceMesh.point(v0), startFace);
            surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);
            auto endPoint = computeXYZ(surfaceMesh, end, startFace);
            add(endPoint, 1);
            assert(isEqual(endPoint,  surfaceMesh.point(v0)));
        }

        {
            // centroid to midPoint
            add( computeCentroid(startFace), 0);
            traceLine(computeCentroid(startFace), CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v2)));
            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh, computeCentroid(startFace), startFace);
            auto targetPoint = computeUV(surfaceMesh,  CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v2)), startFace);
            surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);
            auto endPoint = computeXYZ(surfaceMesh, end, startFace);
            add(endPoint, 1);
            assert(isEqual(endPoint,   CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v2))));
        }
        {
            //add( CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1)), 0);
            //add(surfaceMesh.point(v0), 0);
            // midEdge to point
            traceLine( CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1)), surfaceMesh.point(v1));
            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh,   CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1)), startFace);
            auto targetPoint = computeUV(surfaceMesh, surfaceMesh.point(v1), startFace);
            surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);
            Point endPoint;
            if(faceEnd != Mesh::null_face())
            {
                endPoint = computeXYZ(surfaceMesh, end, startFace);
            }
            add(endPoint, 1);
           assert(isEqual(endPoint,  surfaceMesh.point(v1)));
        }


        {

            // centroid to midPoint
            startFace = Mesh::Face_index(3);
            auto centroid =   computeCentroid(startFace);
            traceLine(surfaceMesh.point(v0), centroid);

            Mesh::Face_index faceEnd;
            K::Point_2 end;
            auto startPoint = computeUV(surfaceMesh,surfaceMesh.point(v0), startFace);
            auto targetPoint = computeUV(surfaceMesh,  centroid, startFace);

            surfaceOptimizer.intersection2D(startFace,
                                            startPoint,
                                            targetPoint,
                                            faceEnd, end);
            auto endPoint = computeXYZ(surfaceMesh, end, startFace);
            add(endPoint, 1);
            assert(isEqual(endPoint,   CGAL::midpoint(surfaceMesh.point(v3), surfaceMesh.point(v2))));

        }

        display();



    }


    void test1(SurfaceOptimizer optimizer, const Mesh& surfaceMesh)
    {
        // cube corner
        // Vertex 0 corner

        // Vertex 2 link to 0
        // Vertex 6 link to 0
        // Vertex 1 link to 0
        // Vertex 3 link to 0

        // face 4 top
        // face 5 left
        // face 2 left
        // face 3 right

        std::vector<Point> directions;
        std::vector<Point> out;
        std::vector<int> ids;
        auto add = [&](const Point& p, int id)
        {
            out.push_back(p);
            ids.push_back(id);
        };
        auto display = [&]()
        {
            auto pcl = polyscope::registerPointCloud("points", out);
            pcl->addScalarQuantity("ids", ids);
            pcl->resetTransform();

            auto pclDirections = polyscope::registerPointCloud("directions", directions);
            pclDirections->resetTransform();
        };

        size_t id = 0;

        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };



        /*auto startFace =  Mesh::Face_index(5);
        auto v0 = Mesh::Vertex_index (7);
        auto v1 = Mesh::Vertex_index (3);
        auto p =
                CGAL::midpoint(surfaceMesh.point(v0), surfaceMesh.point(v1));
        auto optimalPoint = surfaceMesh.point(Mesh::Vertex_index (5));

        add(p,0);
        add(optimalPoint,1);
        optimizer.optimize(p, optimalPoint, startFace);
        directions = optimizer.directions;
        add(p,2);
        display();*/





        auto optimalPoint = surfaceMesh.point(Mesh::Vertex_index (0));
        size_t i = 0;
        for(auto f : surfaceMesh.faces())
        {
            auto centroid = computeCentroid(f);
            SurfaceOptimizer optimizer2(surfaceMesh);
            optimizer2.parameters.deltaXY = 0.1;
            optimizer2.max_iteration = 100;
            optimizer2.optimize(centroid, optimalPoint, f);
            //directions.insert(directions.end(), optimizer2.directions.begin(), optimizer2.directions.end());
            auto pclDirections = polyscope::registerPointCloud("directions_"+std::to_string(i),  optimizer2.directions);
            pclDirections->resetTransform();
            pclDirections->setEnabled(false);

            /*
            auto pclDirections2 = polyscope::registerPointCloud("gradients"+std::to_string(i),  optimizer2.debugPath);
            pclDirections2->addVectorQuantity("normals", optimizer2.debugNormal);
            pclDirections2->resetTransform();


            auto pclGradient = polyscope::registerPointCloud("gradients",  optimizer2.gradients);
            pclGradient->resetTransform();

            auto pclPointDirections = polyscope::registerPointCloud("pointDirections",  optimizer2.pointCandidates);
            pclPointDirections->resetTransform();*/
            i++;
        }
        //display();
    }
    void testOptimize(const Mesh& surfaceMesh)
    {

        auto f = Mesh::Face_index (0);

        auto range = vertices_around_face(surfaceMesh.halfedge(f), surfaceMesh);
        auto it = range.begin();

        auto v0 = *it;      // first vertex
        ++it;               // advance the iterator
        auto v1 = *it;
        Point cv = surfaceMesh.point(v0);
        SurfaceOptimizer optimizer2(surfaceMesh);

        auto opt = Point(0.,1.,0.);

        optimizer2.directions.push_back(cv);
        optimizer2.directions.push_back(opt);
        optimizer2.optimizeUV(cv, opt, f);

        {
            auto pcl = polyscope::registerPointCloud("directions", optimizer2.directions);
            pcl->resetTransform();
        }

        {
            auto pcl = polyscope::registerPointCloud("gradients", optimizer2.gradients);
            pcl->resetTransform();
        }
        for(int i =0 ; i < optimizer2.debug2Ds.size(); i++)
        {
            optimizer2.debug2Ds[i].print(i);
        }
       /* optimizer2.parameters.deltaXY = 0.1;
        optimizer2.max_iteration = 100;
        optimizer2.optimize(cv, surfaceMesh.point(v1), f);

        auto pclGradient = polyscope::registerPointCloud("gradients",  optimizer2.gradients);
        pclGradient->resetTransform();

        optimizer2.debug2D.print();

        {
            auto pcl = polyscope::registerPointCloud("debug_gradients", optimizer2.debug_gradient);
            pcl->resetTransform();
        }*/
    }

    void testDelta(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };

        auto p1 = surfaceMesh.point(Mesh::Vertex_index(1));
        auto p2 = surfaceMesh.point(Mesh::Vertex_index(2));
        auto p3 = surfaceMesh.point(Mesh::Vertex_index(3));



        std::vector<Point> opts = {p1,p2,p3};
        auto faceStart = Mesh::Face_index (6);
        auto centralPoint = computeCentroid(faceStart);

        centralPoint               = CGAL::midpoint(p1, p2);
        std::vector<Point > points;
        std::vector<Point > grads;
        points.push_back(centralPoint);
        float deltaXY = 0.1;




        for(auto& p : opts)
        {
            auto getEnergy = [&](K::Point_2 delta)
            {
                Mesh::Face_index outFace;
                Point target;
                auto uv = computeUV(surfaceMesh, centralPoint, faceStart);
                optimizer.moveDeltaSurface(faceStart, uv, delta - CGAL::ORIGIN, target, outFace);
                auto d= CGAL::squared_distance(p, target);
                std::cout<<"d "<<d<<"\n";
                grads.push_back(centralPoint+ (deltaXYZToP(surfaceMesh, delta- CGAL::ORIGIN, faceStart)-CGAL::ORIGIN));
                return d;
            };

            K::Vector_2 gradient = K::Vector_2(
                    getEnergy(K::Point_2(-deltaXY, 0)) - getEnergy(K::Point_2(deltaXY, 0)),
                    getEnergy(K::Point_2(0, -deltaXY)) - getEnergy(K::Point_2(0, deltaXY)));

            std::cout<<"gradient "<<gradient<<"\n";
            Mesh::Face_index outFace;
            Point target;
            auto uv = computeUV(surfaceMesh, centralPoint, faceStart);
            optimizer.moveDeltaSurface(faceStart, uv, gradient, target, outFace);
            points.push_back(target);
        }

        {
            auto pcl = polyscope::registerPointCloud("debug_points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("debug_gradients", grads);
            pcl->resetTransform();
        }


    }
    void testDeltaOptim(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };
        auto optimal = computeCentroid(Mesh::Face_index (2));

        Mesh::Face_index f = Mesh::Face_index (6);
        auto point =  computeCentroid(Mesh::Face_index (6));
        std::vector<Point> points = {point};
        optimizer.optimizeUV(point, optimal, f);
        points.push_back(point);
        points = optimizer.directions;
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("debug_gradients", optimizer.gradients);
            pcl->resetTransform();
        }
    }
    void testIntersection1(const Mesh& surfaceMesh)
    {


        SurfaceOptimizer optimizer(surfaceMesh);
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };
        //auto pt = computeCentroid(Mesh::Face_index (6));
        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;

        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(0. ,0.25);
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            std::cout<<"end "<<end<<"\n";
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }

        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(0. ,0.75);
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            std::cout<<"end "<<end<<"\n";
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }

        {
            K::Point_2 start (0.5 ,0.);
            K::Point_2 target(0. ,0.);
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            std::cout<<"end "<<end<<"\n";
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }

        {
            K::Point_2 start (0.5 ,0.0);
            K::Point_2 target(1. ,0.);
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            std::cout<<"end "<<end<<"\n";
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
    }
    void testIntersection2(const Mesh& surfaceMesh)
    {


        SurfaceOptimizer optimizer(surfaceMesh);
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };
        //auto pt = computeCentroid(Mesh::Face_index (6));
        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;

        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(-0.25 ,0.5);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
            start,
            target-start,
            pointMoved,
            faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }
        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(0.25 ,0.5);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }

        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(0. ,0.75);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }
        {
            K::Point_2 start (0. ,0.5);
            K::Point_2 target(0. ,0.25);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }

        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
    }
    void testIntersection3(const Mesh& surfaceMesh)
    {


        SurfaceOptimizer optimizer(surfaceMesh);
        auto computeCentroid = [&](Mesh::Face_index startFace) {
            auto h = surfaceMesh.halfedge(startFace);

            Point p0 = surfaceMesh.point(source(h, surfaceMesh));
            Point p1 = surfaceMesh.point(target(h, surfaceMesh));
            Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

            Point p = CGAL::centroid(p0, p1, p2);
            return p;
        };
        //auto pt = computeCentroid(Mesh::Face_index (6));
        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;

        {
            K::Point_2 start (0.0 ,0.0);
            K::Point_2 target(-0.25 ,0.0);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            K::Point_2 end;
            optimizer.moveDeltaSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);

            auto newStart = computeUV(surfaceMesh, pointMoved, faceEnd);
            /*optimizer.moveDeltaSurface(faceEnd,
                                       newStart,
                                       target-newStart,
                                       pointMoved,
                                       faceEnd);*/

            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                            faceEnd, end);
            auto endPoint = computeXYZ(surfaceMesh, end, faceStart);
            std::cout<<"end "<<end<<"\n";
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);

        }
        /*
        {
            K::Point_2 start (0. ,0.0);
            K::Point_2 target(0.25 ,0.0);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }

        {
            K::Point_2 start (0. ,0.0);
            K::Point_2 target(0. ,0.25);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }*/
        /*{
            K::Point_2 start (0. ,0.0);
            K::Point_2 target(0. ,-0.25);
            Point pointMoved;
            points.push_back(computeXYZ(surfaceMesh, start, faceStart));
            points.push_back(computeXYZ(surfaceMesh, target, faceStart));

            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            std::cout<<"faceEnd "<<faceEnd<<"\n";
            points.push_back(pointMoved);
        }*/

        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
    }

    void testIntersection2D(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Point > targets;
        K::Point_2 start (0. ,0.5);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {

            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            assert("end == target" && target == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            // compute intersection
            assert("end == target" && K::Point_2(0.5,0.5) == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            assert("end == target" && target == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && start == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }

    }
    void testIntersection2D_1(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Point > targets;
        K::Point_2 start (0. ,0.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {

            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            assert("end == target" && target == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            // compute intersection
            assert("end == target" && target == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            assert("end == target" && K::Point_2(0.01 ,0.) == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            //assert("end == target" && start == end);
            assert("end == target" && K::Point_2(0., 0.01 ) == end);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }

    }
    void testIntersection2D_2(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Point > targets;
        K::Point_2 start (0. ,1.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {

            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            std::cout << "value" << ("end == target" && target == end)<<"\n";
            assert("end == target" && end == K::Point_2 (0.01, 0.99));
            assert("face" && faceEnd.idx() == 11);
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            // compute intersection
            assert("end == target" && end == K::Point_2 (0.0, 0.99));
            assert("face" && faceEnd.idx() == 1);

            std::cout << "value" << ("end == target" && target == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.0, 0.75));
            assert("face" && faceEnd.idx() == 1);

            std::cout << "value" << ("end == target" && K::Point_2(0.01 ,0.) == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && end == K::Point_2 (0.0, 0.99));
            assert("face" && faceEnd.idx() == 1);
            //std::cout << "value" << ("end == target" && start == end)<<"\n;
            std::cout << "value" << ("end == target" && K::Point_2(0., 0.01 ) == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }

    }
    void testIntersection2D_3(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Point > targets;
        K::Point_2 start (1. ,0.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {

            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );


            assert("end == target" && end == K::Point_2 (0.99, 0.));
            assert("face" && faceEnd.idx() == 0);


            std::cout << "value" << ("end == target" && target == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.01));
            assert("face" && faceEnd.idx() == 11);
            // compute intersection
            std::cout << "value" << ("end == target" && target == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.));
            assert("face" && faceEnd.idx() == 0);

            std::cout << "value" << ("end == target" && K::Point_2(0.01 ,0.) == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && end == K::Point_2 (0.75, 0.));
            assert("face" && faceEnd.idx() == 0);
            //std::cout << "value" << ("end == target" && start == end)<<"\n;
            std::cout << "value" << ("end == target" && K::Point_2(0., 0.01 ) == end)<<"\n";
            std::cout<<"OK end "<<end<<" faceEnd "<<faceEnd<<"\n";
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }

    }

    void testOptimizationSurface1(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Vector> directions;
        std::vector<Point > targets;
        std::vector<Point > ends;
        K::Point_2 start (1. ,0.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {
            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                               start,
                               target-start,
                               pointMoved,
                               faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        {
            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        /*
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.01));
            assert("face" && faceEnd.idx() == 11);
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && end == K::Point_2 (0.75, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }*/
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->addVectorQuantity("normal", directions);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("ends", ends);
            pcl->resetTransform();
        }

    }
    void testOptimizationSurface2(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Vector> directions;
        std::vector<Point > targets;
        std::vector<Point > ends;
        K::Point_2 start (0. ,1.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {
            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        {
            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        /*
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.01));
            assert("face" && faceEnd.idx() == 11);
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && end == K::Point_2 (0.75, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }*/
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->addVectorQuantity("normal", directions);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("ends", ends);
            pcl->resetTransform();
        }

    }
    void testOptimizationSurface3(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);

        auto faceStart = Mesh::Face_index (6);
        auto h = surfaceMesh.halfedge(faceStart);

        Point p0 = surfaceMesh.point(source(h, surfaceMesh));
        Point p1 = surfaceMesh.point(target(h, surfaceMesh));
        Point p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));
        std::vector<Point > points;
        std::vector<Vector> directions;
        std::vector<Point > targets;
        std::vector<Point > ends;
        K::Point_2 start (0. ,0.0);
        //points.push_back(computeXYZ(surfaceMesh, start, faceStart));
        {
            K::Point_2 target = start + K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        {
            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }
        {
            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            Point pointMoved;
            Mesh::Face_index  faceEnd;
            optimizer.moveOnSurface(faceStart,
                                    start,
                                    target-start,
                                    pointMoved,
                                    faceEnd);
            //assert("end == target" && end == K::Point_2 (0.99, 0.));
            //assert("face" && faceEnd.idx() == 0);
            std::cout<<"pointMoved "<<pointMoved<<"\n";
            std::cout<<"faceEnd "<<faceEnd.idx()<<"\n";

            points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
            directions.insert(directions.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());
            ends.push_back(pointMoved);
        }

        /*
        {
            K::Point_2 target = start + K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.01));
            assert("face" && faceEnd.idx() == 11);
            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0. ,0.25);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );

            assert("end == target" && end == K::Point_2 (0.99, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }
        {

            K::Point_2 target = start - K::Vector_2 (0.25 ,0.);
            targets.push_back(computeXYZ(surfaceMesh, target, faceStart));
            K::Point_2 end;
            Mesh::Face_index  faceEnd;
            optimizer.intersection2D(faceStart,
                                     start,
                                     target,
                                     faceEnd,  end );
            // intersect at start ok
            assert("end == target" && end == K::Point_2 (0.75, 0.));
            assert("face" && faceEnd.idx() == 0);

            points.push_back(computeXYZ(surfaceMesh, end, faceStart));
        }*/
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->addVectorQuantity("normal", directions);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("ends", ends);
            pcl->resetTransform();
        }

    }

    void testCube1(const Mesh& surfaceMesh)
    {
        std::vector<Point> targets;
        SurfaceOptimizer optimizer(surfaceMesh);
        std::vector<Point > points;
        std::vector<Point > candidates;
        std::vector<Vector > normals;

        std::vector<Point > gradients;
        auto startFace = Mesh::Face_index (5);
        Point point(-1.,1.,0.);
        Point optimal(-1.,-1.,1.);
        targets.push_back(point);
        targets.push_back(optimal);
        optimizer.optimize(point, optimal, startFace);

        gradients.insert(gradients.end(), optimizer.gradients.begin(), optimizer.gradients.end());

        candidates.insert(candidates.end(), optimizer.pointCandidates.begin(), optimizer.pointCandidates.end());
        points.insert(points.end(), optimizer.debugPath.begin(), optimizer.debugPath.end());
        normals.insert(normals.end(), optimizer.debugNormal.begin(), optimizer.debugNormal.end());


        targets.push_back(point);
        {
            auto pcl = polyscope::registerPointCloud("gradients", gradients);
            pcl->resetTransform();
        }

        {
            auto pcl = polyscope::registerPointCloud("targets", targets);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("candidates", candidates);
            pcl->resetTransform();
        }
        {
            auto pcl = polyscope::registerPointCloud("points", points);
            pcl->addVectorQuantity("normal", normals);
            pcl->resetTransform();
        }
    }
    void test(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);
        //test1(optimizer, surfaceMesh);

        //test2(optimizer, surfaceMesh);
        //test3(optimizer, surfaceMesh);
        //testIntersection(surfaceMesh);

        //testIntersectionCircle(surfaceMesh);
        //testOptimize( surfaceMesh);
        //test1(optimizer, surfaceMesh);


        //testDelta(surfaceMesh);
        //testDeltaOptim(surfaceMesh);
        //testIntersection1(surfaceMesh);
        //testIntersection2(surfaceMesh);
        //testIntersection3(surfaceMesh);


        // OK
        //testIntersection2D(surfaceMesh);
        //testIntersection2D_1(surfaceMesh);
        //testIntersection2D_2(surfaceMesh);
        //testIntersection2D_3(surfaceMesh);

        //testIntersectionCircle(surfaceMesh);

        //testOptimizationSurface1(surfaceMesh);
        //testOptimizationSurface2(surfaceMesh);
        //testOptimizationSurface3(surfaceMesh);

        //testOptimizationSurfaceCircle(surfaceMesh);

        //test1(optimizer, surfaceMesh);
        testCube1(surfaceMesh);
    }


};