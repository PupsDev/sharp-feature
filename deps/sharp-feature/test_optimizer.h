#pragma once
//#include "graph_optimizer.h"
#include "surface_optimizer.h"
namespace test_graph_optimizer
{
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
    void test2(SurfaceOptimizer optimizer, const Mesh& surfaceMesh)
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
    void test3(SurfaceOptimizer optimizer, const Mesh& surfaceMesh)
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
    void test(const Mesh& surfaceMesh)
    {
        SurfaceOptimizer optimizer(surfaceMesh);
        //test1(optimizer, surfaceMesh);
        //test2(optimizer, surfaceMesh);
        test3(optimizer, surfaceMesh);

    }

};