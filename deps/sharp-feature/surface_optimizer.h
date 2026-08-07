#pragma once
#include "coord.h"
#include "optimizer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <CGAL/Ray_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Point_2.h>
struct SurfaceOptimizerParameters
{
    float lineQuadricWeight = 0.001f;
    float deltaXY = 0.025f;
};

struct DEBUG_2D
{
    K::Point_2 start;
    K::Point_2 end;

    void print()
    {
        std::vector<Point > geometry;
        std::vector<Point > tri;
        Point a (0, 0. , 0.);
        Point b (1, 0. , 0.);
        Point c (0, 1. , 0.);

        tri.push_back(a);
        tri.push_back(b);
        tri.push_back(c);

        std::vector<std::array<size_t, 3>> triangle;
        triangle.push_back({0,1,2});

        std::vector<std::array<size_t, 2>> lines;
        geometry.push_back(Point(start.x(), start.y(),0.));
        geometry.push_back(Point(end.x(), end.y(),0.));
        lines.push_back({0,1});
        auto curve = polyscope::registerCurveNetwork("line", geometry, lines);
        curve->resetTransform();

        auto surfaceTriangle2D = polyscope::registerSurfaceMesh("triangle2D", tri, triangle);
        surfaceTriangle2D->resetTransform();
    }
    void print(int i)
    {
        std::vector<Point > geometry;
        std::vector<Point > tri;
        Point a (0, 0. , 0.);
        Point b (1, 0. , 0.);
        Point c (0, 1. , 0.);

        tri.push_back(a);
        tri.push_back(b);
        tri.push_back(c);

        std::vector<std::array<size_t, 3>> triangle;
        triangle.push_back({0,1,2});

        std::vector<std::array<size_t, 2>> lines;
        geometry.push_back(Point(start.x(), start.y(),0.));
        geometry.push_back(Point(end.x(), end.y(),0.));
        lines.push_back({0,1});
        auto curve = polyscope::registerCurveNetwork("line" + std::to_string(i), geometry, lines);
        curve->resetTransform();

        auto surfaceTriangle2D = polyscope::registerSurfaceMesh("triangle2D"+ std::to_string(i), tri, triangle);
        surfaceTriangle2D->resetTransform();
    }
};
class SurfaceOptimizer {
public:

    SurfaceOptimizer() {

    }

    SurfaceOptimizer(const CGAL::Surface_mesh <Point> &inputMesh) {
        surfaceMesh = inputMesh;
    }

    void moveDelta(const Mesh::Face_index & faceStart,
                   const K::Vector_2& start,
                   const K::Vector_2& target,
                   Mesh::Face_index & faceEnd, K::Vector_2& end )
    {

    }
    static K::Point_2 clampUV( const K::Point_2& uv)
    {
        double x = uv.x();
        double y = uv.y();
        float eps = 10e-3;
        if( x < eps) x = eps;
        if( y < eps) y = eps;

        if( x > 1-eps) x = 1-eps;
        if( y > 1-eps) y = 1-eps;
        return {x,y};
    }
    static K::Vector_2  normalize( const K::Vector_2 & uv)
    {
        return uv / std::sqrt(CGAL::squared_length(uv));
    }
    bool intersectionSegmentTriangleEdge(
            const K::Point_2& start,
            const K::Point_2& outside,
            const K::Point_2& a,
            const K::Point_2& b,
            double& t,
            K::Point_2& p)
    {
        constexpr double eps = 1e-12;

        auto cross = [](const K::Vector_2& u, const K::Vector_2& v)
        {
            return u.x() * v.y() - u.y() * v.x();
        };

        K::Vector_2 r = outside - start;
        K::Vector_2 s = b - a;

        double den = cross(r, s);

        // ---------------------------------------------------------
        // Parallel / collinear
        // ---------------------------------------------------------
        if (std::abs(den) < eps)
        {
            // Not collinear
            if (std::abs(cross(a - start, r)) > eps)
                return false;

            // Project edge endpoints onto the path
            double rr = CGAL::to_double(r.squared_length());

            double t0 = CGAL::to_double((a - start) * r) / rr;
            double t1 = CGAL::to_double((b - start) * r) / rr;

            if (t0 > t1)
                std::swap(t0, t1);

            // No overlap
            if (t1 < 0.0 || t0 > 1.0)
                return false;

            // First point encountered
            t = std::max(0.0, t0);

            p = start + t * r;

            return true;
        }

        // ---------------------------------------------------------
        // Regular case
        // ---------------------------------------------------------

        K::Vector_2 ap = a - start;

        double tt = cross(ap, s) / den;
        double uu = cross(ap, r) / den;

        if (tt < -eps || tt > 1.0 + eps)
            return false;

        if (uu < -eps || uu > 1.0 + eps)
            return false;

        t = tt;
        p = start + t * r;

        return true;
    }
    void intersection2D2(const Mesh::Face_index & faceStart,
                        const K::Point_2& start,
                        const K::Point_2& target,
                        Mesh::Face_index & faceEnd, K::Point_2& end )
    {
        double bestT = std::numeric_limits<double>::max();

        K::Vector_2 dir = normalize(target - start);
        K::Point_2 outside = start + 100.0 * dir;

        for (auto h : CGAL::halfedges_around_face(surfaceMesh.halfedge(faceStart), surfaceMesh))
        {
            K::Point_2 v0 = computeUV(surfaceMesh,
                                      surfaceMesh.point(surfaceMesh.source(h)),
                                      faceStart);

            K::Point_2 v1 = computeUV(surfaceMesh,
                                      surfaceMesh.point(surfaceMesh.target(h)),
                                      faceStart);

            double t;
            K::Point_2 p;

            if (!intersectionSegmentTriangleEdge(start, outside, v0, v1, t, p))
                continue;

            // Ignore the starting point
            if (t < 1e-8)
                continue;

            if (t < bestT)
            {
                bestT = t;
                end = p;
                faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
            }
        }
    }
    bool checkBoundaries(const K::Point_2& vec)
    {
        float eps = 10e-5;
        return (vec[0]  < - eps  ) ||
               (vec[1]  < - eps  ) ||
               ( (vec[1]+vec[0]) > (1+eps)    );
    }
    void moveDeltaSurface(const Mesh::Face_index & faceStart,
                          const K::Point_2 start,
                          const K::Vector_2 delta,
                          Point& target,
                          Mesh::Face_index & outFace)
    {
        auto computeQ = [](const Vector &d,
                           const Vector &n1,
                           const Vector &n2) {
            auto n1norm = std::sqrt(n1.squared_length());
            auto n2norm = std::sqrt(n2.squared_length());

            auto a = CGAL::cross_product(n1, n2);
            a = a / std::sqrt(a.squared_length()); // normalize

            auto n1u = n1 / n1norm;
            auto n2u = n2 / n2norm;

            auto sinphiv = CGAL::cross_product(n1u, n2u);
            double sinphi = std::sqrt(sinphiv.squared_length());
            double cosphi = n1u * n2u; // dot product

            return d * cosphi
                   + CGAL::cross_product(a, d) * sinphi
                   + a * (a * d) * (1.0 - cosphi);
        };
        auto computeDirection3D = [&](const Point &globalStartPoint,
                                      const Point &globalEndPoint,
                                      const Vector &initialNormal,
                                      const Vector &nextNormal) {
            auto direction3D = globalEndPoint - globalStartPoint;
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            if (CGAL::to_double(initialNormal * nextNormal) > 0.99) {
                direction3D = globalEndPoint - globalStartPoint;
            } else {
                direction3D = computeQ(direction3D, initialNormal,
                                       nextNormal);
            }
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            return direction3D;
        };


        if(CGAL::squared_length(delta) < 1e-5) return;

        Point startPoint =  computeXYZ(surfaceMesh, start, faceStart);
        K::Vector_2 deltaUVX = deltaXYZToUVOld(surfaceMesh, delta, startPoint, faceStart);
        Point endPoint =  computeXYZ(surfaceMesh, start+deltaUVX, faceStart);

        std::cout<<"delta "<<delta<<"\n";

        std::cout<<"deltaUVX "<<deltaUVX<<"\n";
        // check if endPoint is inside surface
        if (barycentricCoordinates(surfaceMesh, faceStart, endPoint))
        {
            target = endPoint;
            outFace = faceStart;
            return;
        }

        Mesh::Face_index currentFace = faceStart;
        K::Point_2 currentStart2D = start;
        K::Point_2 end2D = start+deltaUVX;

        std::cout<<"currentStart2D "<<currentStart2D<<"\n";
        std::cout<<"end2D "<<end2D<<"\n";
        /*
        currentStart2D = clampUV(currentStart2D);

        double x = end2D.x();
        double y = end2D.y();
        float eps = 10e-3;

        if( x < eps &&  y > 1-eps)
        {
            x = 0.;
            y = 1-eps;
        }
        if( y < eps &&  x > 1-eps)
        {
            y = 0.;
            x = 1-eps;
        }
        end2D = K::Point_2 (x,y);*/


        Mesh::Face_index nextFace;
        K::Point_2 intersectionPoint2D;
        intersection2D(currentFace, currentStart2D, end2D, nextFace, intersectionPoint2D);

        Point intersectionPoint =  computeXYZ(surfaceMesh, intersectionPoint2D, currentFace);

        target = intersectionPoint;
        outFace = nextFace;

        auto deltaXYZ = intersectionPoint - startPoint;
        auto deltaDiff = deltaXYZToDeltaUV(surfaceMesh, deltaXYZ,  faceStart);
        auto newDelta = deltaUVX - deltaDiff;

        return; // moveDeltaSurface(nextFace, intersectionPoint2D, newDelta, target, outFace);
    }
    void move(const Mesh::Face_index & faceStart,
              const K::Point_2 start2D,
              const K::Vector_2 delta,
              Mesh::Face_index & outFace,
              K::Point_2 & outEnd2D
              )
    {
        auto computeQ = [](const Vector &d,
                           const Vector &n1,
                           const Vector &n2) {
            auto n1norm = std::sqrt(n1.squared_length());
            auto n2norm = std::sqrt(n2.squared_length());

            auto a = CGAL::cross_product(n1, n2);
            a = a / std::sqrt(a.squared_length()); // normalize

            auto n1u = n1 / n1norm;
            auto n2u = n2 / n2norm;

            auto sinphiv = CGAL::cross_product(n1u, n2u);
            double sinphi = std::sqrt(sinphiv.squared_length());
            double cosphi = n1u * n2u; // dot product

            return d * cosphi
                   + CGAL::cross_product(a, d) * sinphi
                   + a * (a * d) * (1.0 - cosphi);
        };
        auto computeDirection3D = [&](const Point &globalStartPoint,
                                      const Point &globalEndPoint,
                                      const Vector &initialNormal,
                                      const Vector &nextNormal) {
            auto direction3D = globalEndPoint - globalStartPoint;
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            if (CGAL::to_double(initialNormal * nextNormal) > 0.99) {
                direction3D = globalEndPoint - globalStartPoint;
            } else {
                direction3D = computeQ(direction3D, initialNormal,
                                       nextNormal);
            }
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            return direction3D;
        };





        auto target = computeXYZ(surfaceMesh, start2D + delta, faceStart);
        auto start = computeXYZ(surfaceMesh, start2D, faceStart);

        if (barycentricCoordinates(surfaceMesh, faceStart, target)) {
            outEnd2D = start2D+delta;
            outFace = faceStart;
            return;
        }
        auto directionSurface = target - start;
        double totalDistance = std::sqrt(CGAL::squared_length(directionSurface));
        double remainingDistance = totalDistance;
        directionSurface /= totalDistance;



        Mesh::Face_index currentFace = faceStart;
        Point currentStart = start;
        Point currentTarget = target;
        auto end2D = start2D+delta;



        int it = 0;
        while(remainingDistance > 0.)
        {
            Vector initialNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(currentFace,surfaceMesh);
            auto currentStart2D = computeUV(surfaceMesh, currentStart, currentFace);

            /*K::Vector_2 deltaUVX = deltaXYZToUVOld(surfaceMesh, end2D-currentStart2D, currentStart, faceStart);
            end2D = currentStart2D + deltaUVX;
            */
            /*
            std::cout<<"currentTarget "<<currentTarget<<"\n";
            std::cout<<"currentStart2D "<<currentStart2D<<"\n";
            std::cout<<"end2D "<<end2D<<"\n";*/

            Mesh::Face_index nextFace;
            K::Point_2 intersectionPoint2D;
            intersection2D(currentFace, currentStart2D, end2D, nextFace, intersectionPoint2D);

            Point intersectionPoint =  computeXYZ(surfaceMesh, intersectionPoint2D, currentFace);

            Vector nextNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(nextFace,surfaceMesh);
            Vector direction3D;
            double distanceToIntersection = 0.;
            if(nextFace == faceStart)
            {
                direction3D = currentTarget - intersectionPoint;
                direction3D /= sqrt(CGAL::squared_length(direction3D));
                distanceToIntersection = std::sqrt(CGAL::squared_length(intersectionPoint - currentStart));
            }
            else
            {
                direction3D = computeDirection3D(intersectionPoint, currentTarget, initialNormal, nextNormal);
                distanceToIntersection = std::sqrt(CGAL::squared_length(intersectionPoint - currentStart));
            }
            /*
            std::cout<<"Intersection "<<intersectionPoint2D<<"\n";
            std::cout<<"nextFace "<<nextFace<<"\n";
             */


            auto previousStart = currentStart;
            currentStart = intersectionPoint;

            remainingDistance-= distanceToIntersection;
            if (remainingDistance <= 0)
            {
                //auto d = distanceToIntersection+remainingDistance;
                //target = previousStart + d * direction3D;
                target = currentStart;
                outFace = currentFace;
                return;
            }
            currentTarget = currentStart + remainingDistance * direction3D;
            Point step = currentTarget;

            auto bary = barycentricCoordinates(surfaceMesh, nextFace, currentTarget);

            debugPath.push_back(intersectionPoint);
            debugNormal.push_back(direction3D);

            if(bary)
            {
                target = step;
                outEnd2D = computeUV(surfaceMesh, target, nextFace);
                outFace = nextFace;

                debugPath.push_back(target);
                debugNormal.push_back(direction3D);
                return ;
            }
            //currentTarget = currentStart;

            currentFace = nextFace;
            end2D = computeUV(surfaceMesh, currentTarget, currentFace);

            it++;
        }
    }
    void moveOnSurface(const Mesh::Face_index & faceStart,
                       const Point& start,
                       Point& target,
                       Mesh::Face_index & outFace,
                       float totalDistance) {
        auto computeQ = [](const Vector &d,
                           const Vector &n1,
                           const Vector &n2) {
            auto n1norm = std::sqrt(n1.squared_length());
            auto n2norm = std::sqrt(n2.squared_length());

            auto a = CGAL::cross_product(n1, n2);
            a = a / std::sqrt(a.squared_length()); // normalize

            auto n1u = n1 / n1norm;
            auto n2u = n2 / n2norm;

            auto sinphiv = CGAL::cross_product(n1u, n2u);
            double sinphi = std::sqrt(sinphiv.squared_length());
            double cosphi = n1u * n2u; // dot product

            return d * cosphi
                   + CGAL::cross_product(a, d) * sinphi
                   + a * (a * d) * (1.0 - cosphi);
        };
        auto computeDirection3D = [&](const Point &globalStartPoint,
                                      const Point &globalEndPoint,
                                      const Vector &initialNormal,
                                      const Vector &nextNormal) {
            auto direction3D = globalEndPoint - globalStartPoint;
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            if (CGAL::to_double(initialNormal * nextNormal) > 0.99) {
                direction3D = globalEndPoint - globalStartPoint;
            } else {
                direction3D = computeQ(direction3D, initialNormal,
                                       nextNormal);
            }
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            return direction3D;
        };


        auto directionSurface = target - start;
        directionSurface /= std::sqrt(CGAL::squared_length(directionSurface));
        auto pointMoved = start + totalDistance * directionSurface;
        if (barycentricCoordinates(surfaceMesh, faceStart, pointMoved)) {
            target = pointMoved;
            outFace = faceStart;
            return;
        }


        double remainingDistance = totalDistance;
        Mesh::Face_index currentFace = faceStart;
        Point currentStart = start;
        Point currentTarget = target;
        auto end2D = computeUV(surfaceMesh, currentTarget, faceStart);



        int it = 0;
        while(remainingDistance > 0.)
        {
            Vector initialNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(currentFace,surfaceMesh);
            auto currentStart2D = computeUV(surfaceMesh, currentStart, currentFace);

            /*K::Vector_2 deltaUVX = deltaXYZToUVOld(surfaceMesh, end2D-currentStart2D, currentStart, faceStart);
            end2D = currentStart2D + deltaUVX;
            */
            /*std::cout<<"currentTarget "<<currentTarget<<"\n";

            std::cout<<"currentStart2D "<<currentStart2D<<"\n";
            std::cout<<"end2D "<<end2D<<"\n";*/

            Mesh::Face_index nextFace;
            K::Point_2 intersectionPoint2D;
            intersection2D(currentFace, currentStart2D, end2D, nextFace, intersectionPoint2D);

            Point intersectionPoint =  computeXYZ(surfaceMesh, intersectionPoint2D, currentFace);

            Vector nextNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(nextFace,surfaceMesh);
            Vector direction3D;
            double distanceToIntersection = 0.;
            if(nextFace == faceStart)
            {
                direction3D = currentTarget - intersectionPoint;
                direction3D /= sqrt(CGAL::squared_length(direction3D));
                distanceToIntersection = std::sqrt(CGAL::squared_length(intersectionPoint - currentStart));
            }
            else
            {
                direction3D = computeDirection3D(intersectionPoint, currentTarget, initialNormal, nextNormal);
                distanceToIntersection = std::sqrt(CGAL::squared_length(intersectionPoint - currentStart));
            }
            /*std::cout<<"Intersection "<<intersectionPoint2D<<"\n";
            std::cout<<"nextFace "<<nextFace<<"\n";*/


            auto previousStart = currentStart;
            currentStart = intersectionPoint;

            remainingDistance-= distanceToIntersection;
            if (remainingDistance <= 0)
            {
                //auto d = distanceToIntersection+remainingDistance;
                //target = previousStart + d * direction3D;
                target = currentStart;
                outFace = currentFace;
                return;
            }
            currentTarget = currentStart + remainingDistance * direction3D;
            Point step = currentTarget;

            auto bary = barycentricCoordinates(surfaceMesh, nextFace, currentTarget);

            //debugPath.push_back(intersectionPoint);
            //debugNormal.push_back(direction3D);

            if(bary)
            {
                target = step;
                outFace = nextFace;
                return ;
            }
            //currentTarget = currentStart;

            currentFace = nextFace;
            end2D = computeUV(surfaceMesh, currentTarget, currentFace);

            it++;
        }
    }
    void moveOnSurface(const Mesh::Face_index & faceStart,
                       const K::Point_2 start,
                       const K::Vector_2 delta,
                       Point& pointMoved,
                       Mesh::Face_index & outFace) {
        auto computeQ = [](const Vector &d,
                           const Vector &n1,
                           const Vector &n2) {
            auto n1norm = std::sqrt(n1.squared_length());
            auto n2norm = std::sqrt(n2.squared_length());

            auto a = CGAL::cross_product(n1, n2);
            a = a / std::sqrt(a.squared_length()); // normalize

            auto n1u = n1 / n1norm;
            auto n2u = n2 / n2norm;

            auto sinphiv = CGAL::cross_product(n1u, n2u);
            double sinphi = std::sqrt(sinphiv.squared_length());
            double cosphi = n1u * n2u; // dot product

            return d * cosphi
                   + CGAL::cross_product(a, d) * sinphi
                   + a * (a * d) * (1.0 - cosphi);
        };
        auto computeDirection3D = [&](const Point &globalStartPoint,
                                      const Point &globalEndPoint,
                                      const Vector &initialNormal,
                                      const Vector &nextNormal) {
            auto direction3D = globalEndPoint - globalStartPoint;
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            if (CGAL::to_double(initialNormal * nextNormal) > 0.99) {
                direction3D = globalEndPoint - globalStartPoint;
            } else {
                direction3D = computeQ(direction3D, initialNormal,
                                       nextNormal);
            }
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            return direction3D;
        };



        pointMoved =  computeXYZ(surfaceMesh, start+delta, faceStart);
        debug_gradient.push_back(pointMoved);

        if (barycentricCoordinates(surfaceMesh, faceStart, pointMoved)) {
            outFace = faceStart;
            return;
        }
        else
        {

            auto startPoint =  computeXYZ(surfaceMesh, start, faceStart);
            auto totalDistance = std::sqrt(CGAL::squared_length(pointMoved- startPoint));
            moveOnSurface(  faceStart,
                            startPoint,
                            pointMoved,
                            outFace,
                            totalDistance);
            return;
        }


    }
    void intersection2D(const Mesh::Face_index & faceStart,
                        K::Point_2& start,
                        const K::Point_2& target,
                        Mesh::Face_index & faceEnd, K::Point_2& end )
    {
        for (const auto& he : CGAL::halfedges_around_face(surfaceMesh.halfedge(faceStart), surfaceMesh))
        {
            const auto v = surfaceMesh.source(he);
            const auto uv = computeUV(surfaceMesh, surfaceMesh.point(v), faceStart);

            if(start == uv)
            {
                if(start.x() == 0. && start.y() == 0.)
                {
                    auto x = 0.01;
                    auto y = start.y();
                    auto candidate = K::Point_2(x, y);
                    faceEnd = surfaceMesh.face(surfaceMesh.opposite(he));
                    const K::Segment_2 seg(candidate  , target);
                    if (seg.has_on(uv))
                    {
                        candidate = K::Point_2(start.x(), 0.01);
                        faceEnd = surfaceMesh.face(surfaceMesh.opposite(surfaceMesh.prev(he)));
                    }
                    start = candidate;
                }
                else
                {
                    if(start.x() == 0. && start.y() == 1.)
                    {
                        auto candidate = K::Point_2(start.x(), start.y() - 0.01);
                        faceEnd = surfaceMesh.face(surfaceMesh.opposite(he));
                        const K::Segment_2 seg(candidate  , target);
                        if (seg.has_on(uv))
                        {
                            //candidate = K::Point_2(1., start.y()-0.01);
                            auto next = surfaceMesh.prev(he);
                            const auto uv2 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.target(next)), faceStart);
                            const auto uv1 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.source(next)), faceStart);
                            candidate = start - 0.01*(uv2-uv1);
                            auto f = surfaceMesh.face(surfaceMesh.opposite(next));
                            faceEnd = f;
                        }
                        start = candidate;
                    }
                    else if(start.x() == 1. && start.y() == 0.)
                    {
                        auto candidate = K::Point_2(start.x()-0.01, start.y());
                        faceEnd = surfaceMesh.face(surfaceMesh.opposite(surfaceMesh.prev(he)));
                        const K::Segment_2 seg(candidate  , target);
                        if (seg.has_on(uv))
                        {
                            //candidate = K::Point_2(start.x()-0.01,1.);
                            auto next = he;
                            auto t = surfaceMesh.target(next);
                            auto s = surfaceMesh.source(next);
                            auto pt = surfaceMesh.point(t);
                            auto ps = surfaceMesh.point(s);
                            const auto uv2 = computeUV(surfaceMesh, pt, faceStart);
                            const auto uv1 = computeUV(surfaceMesh, ps, faceStart);
                            auto vec = uv2-uv1;
                            candidate = start + 0.01*(vec);
                            auto f = surfaceMesh.face(surfaceMesh.opposite(he));
                            faceEnd = f;
                        }
                        start = candidate;
                    }

                }
                faceEnd = faceStart;
            }
        }

        K::Vector_2 dir = target - start;
        dir = normalize(dir);

        K::Point_2 outside = start + 100.0 * dir;

        std::optional<K::Point_2> bestPoint;
        double bestDist = std::numeric_limits<double>::max();

        K::Ray_2 ray(outside, start);
        float epsilonLength = 10e-5;
        constexpr double eps2 = 1e-20;

        auto intersectHalfEdge = [&](Mesh::Halfedge_index h)
        {
            const auto v1 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.source(h)), faceStart);
            const auto v2 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.target(h)), faceStart);
            const K::Segment_2 seg(v1  , v2);
            auto update = [&](const double& d,const K::Point_2& p0)
            {
                if (d < bestDist)
                {
                    bestDist = d;
                    bestPoint = p0;
                    faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
                }
            };
            auto result = CGAL::intersection(ray, seg);
            // check if start it belong to the segment ?
            if (seg.has_on(start) && seg.has_on(target))
            {
                // compute t parameter closest to target
                // then the correct point
                K::Vector_2 edge = v2 - v1;

                double len2 = CGAL::to_double(edge.squared_length());
                if (len2 < 1e-20)
                    return ;

                // Parameter of the projection of 'target' on the edge:
                // p = v1 + t * edge
                double t = CGAL::to_double((target - v1) * edge) / len2;

                K::Point_2 p;

                if (t <= 0.0)
                    p = v1;
                else if (t >= 1.0)
                    p = v2;
                else
                    p = v1 + t * edge;

                double d = 0.;//CGAL::to_double(CGAL::squared_distance(ray.source(), p));

                if (d < bestDist)
                {
                    bestDist = d;
                    bestPoint = p;
                    faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
                }
                return;
            }
            if (!result)
            {
                return ;
            }
            if (const auto* p0 = std::get_if<K::Point_2>(&*result))
            {
                if (CGAL::squared_distance(*p0, ray.source()) < eps2)return ;
                // one intersection
                double d = CGAL::to_double(CGAL::squared_distance(ray.source(), *p0));
                update(d, *p0);
            }
            else if (const auto* s = std::get_if<K::Segment_2>(&*result))
            {
                // overlapping segment
                auto d0 = CGAL::to_double(CGAL::squared_distance(ray.source(), s->source()));
                auto d1 = CGAL::to_double(CGAL::squared_distance(ray.source(), s->target()));
                K::Point_2 p = (d0 < d1) ? s->source() : s->target();
                double d = std::min(d0, d1);   // first point encountered
                update(d, p);
            }
        };

        for (const auto& h : CGAL::halfedges_around_face(surfaceMesh.halfedge(faceStart), surfaceMesh)) {
            intersectHalfEdge(h);
        }

        if(bestPoint)
        {
            end = bestPoint.value();
            for (const auto& he : CGAL::halfedges_around_face(surfaceMesh.halfedge(faceStart), surfaceMesh))
            {
                const auto v = surfaceMesh.source(he);
                const auto uv = computeUV(surfaceMesh, surfaceMesh.point(v), faceStart);

                if(end == uv)
                {
                    if(end.x() == 0. && end.y() == 0.)
                    {
                        auto x = 0.01;
                        auto y = end.y();
                        auto candidate = K::Point_2(x, y);
                        faceEnd = surfaceMesh.face(surfaceMesh.opposite(he));
                        const K::Segment_2 seg(candidate  , target);
                        if (seg.has_on(uv))
                        {
                            candidate = K::Point_2(end.x(), 0.01);
                            faceEnd = surfaceMesh.face(surfaceMesh.opposite(surfaceMesh.prev(he)));
                        }
                        end = candidate;
                    }
                    else
                    {
                        if(end.x() == 0. && end.y() == 1.)
                        {
                            auto candidate = K::Point_2(end.x(), end.y() - 0.01);
                            faceEnd = surfaceMesh.face(surfaceMesh.opposite(he));
                            const K::Segment_2 seg(candidate  , target);
                            if (seg.has_on(uv))
                            {
                                //candidate = K::Point_2(1., end.y()-0.01);
                                auto next = surfaceMesh.prev(he);
                                const auto uv2 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.target(next)), faceStart);
                                const auto uv1 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.source(next)), faceStart);
                                candidate = end - 0.01*(uv2-uv1);
                                auto f = surfaceMesh.face(surfaceMesh.opposite(next));
                                faceEnd = f;
                            }
                            end = candidate;
                        }
                        else if(end.x() == 1. && end.y() == 0.)
                        {
                            auto candidate = K::Point_2(end.x()-0.01, end.y());
                            faceEnd = surfaceMesh.face(surfaceMesh.opposite(surfaceMesh.prev(he)));
                            const K::Segment_2 seg(candidate  , target);
                            if (seg.has_on(uv))
                            {
                                //candidate = K::Point_2(end.x()-0.01,1.);
                                auto next = he;
                                auto t = surfaceMesh.target(next);
                                auto s = surfaceMesh.source(next);
                                auto pt = surfaceMesh.point(t);
                                auto ps = surfaceMesh.point(s);
                                const auto uv2 = computeUV(surfaceMesh, pt, faceStart);
                                const auto uv1 = computeUV(surfaceMesh, ps, faceStart);
                                auto vec = uv2-uv1;
                                candidate = end + 0.01*(vec);
                                auto f = surfaceMesh.face(surfaceMesh.opposite(he));
                                faceEnd = f;
                            }
                            end = candidate;
                        }

                    }
                    faceEnd = faceStart;
                }
            }


        }
        else
        {
            std::cout<<"Can't find intersection point \n";
        }

    }
    auto evaluate (const Point& targetPoint, const Point& inputPoint)
    {

        return CGAL::squared_distance(targetPoint, inputPoint) ;
    }
    void optimizeUV(Point &inputPoint, const Point &targetPoint, Mesh::Face_index &inputFace)
    {
        int iteration = 0;
        float distance  = std::sqrt(CGAL::squared_length(targetPoint-inputPoint));
        max_iteration = 40;
        while (iteration < max_iteration &&  distance> 10e-3)
        {

            auto uv = computeUV(surfaceMesh, inputPoint, inputFace);
            uv = clampUV(uv);
            int id = 0;
            gradients.clear();
            auto getEnergy = [&](K::Point_2 delta) {

                Mesh::Face_index outFace;
                Point target;
                moveDeltaSurface(inputFace, uv, delta - CGAL::ORIGIN, target, outFace);
                gradients.push_back(target);
                std::cout<<"id "<<id++<<" evaluate(targetPoint, target) "<<evaluate(targetPoint, target)<<"\n";
                return evaluate(targetPoint, target);
            };

            K::Vector_2 gradient = K::Vector_2(
                    getEnergy(K::Point_2(-parameters.deltaXY, 0)) - getEnergy(K::Point_2(parameters.deltaXY, 0)),
                    getEnergy(K::Point_2(0, -parameters.deltaXY)) - getEnergy(K::Point_2(0, parameters.deltaXY)));


            std::cout<<"gradient "<<gradient<<"\n";
            double gradient_step_sq_len = std::sqrt(CGAL::squared_length(gradient));


            //if(gradient_step_sq_len > parameters.deltaXY)
            gradient = parameters.deltaXY / gradient_step_sq_len * gradient;


            Mesh::Face_index outFace;
            Point target;
            moveDeltaSurface(inputFace, uv, gradient, target, outFace);

            std::cout<<"evaluate(targetPoint, target) "<<evaluate(targetPoint, target)<<"\n";
            inputPoint = target;
            if(outFace != Mesh::null_face())
                inputFace = outFace;
            distance  = std::sqrt(CGAL::squared_length(targetPoint-inputPoint));
            directions.push_back(inputPoint);
            iteration++;
        }

    }
    void optimize(Point &inputPoint, const Point &targetPoint, Mesh::Face_index &inputFace) {

        int iteration = 0;
        debug_steps.clear();
        float distance  = std::sqrt(CGAL::squared_length(targetPoint-inputPoint));
        //max_iteration = 1;
        while (iteration < max_iteration &&  distance> 10e-3) {

            auto uv = computeUV(surfaceMesh, inputPoint, inputFace);

            gradients.clear();
            auto getEnergy = [&](K::Point_2 delta) {

                debug = false;
                Mesh::Face_index outFace;
                Point target;
                moveOnSurface(inputFace, uv, delta-CGAL::ORIGIN, target, outFace);

                gradients.push_back(target);
                debug = true;
                return evaluate(targetPoint, target);
            };
            K::Vector_2 gradient = K::Vector_2(
                    getEnergy(K::Point_2(-parameters.deltaXY, 0)) - getEnergy(K::Point_2(parameters.deltaXY, 0)),
                    getEnergy(K::Point_2(0, -parameters.deltaXY)) - getEnergy(K::Point_2(0, parameters.deltaXY)));

            double gradient_step_sq_len = std::sqrt(CGAL::squared_length(gradient));
            //std::cout<<"gradient energy "<<gradient_step_sq_len<<"\n";

            if(gradient_step_sq_len > parameters.deltaXY)
                gradient = parameters.deltaXY / gradient_step_sq_len * gradient;

            auto pointCandidate = computeXYZ(surfaceMesh, uv+gradient, inputFace);





            auto totalDistance = std::sqrt(CGAL::squared_length(pointCandidate- inputPoint));
            /*moveOnSurface(inputFace,
                          inputPoint,
                          pointCandidate,
                          outFace,
                          totalDistance);*/

            //moveOnSurface(inputFace, uv, gradient, pointCandidate, outFace);


            K::Point_2 end2D;
            Mesh::Face_index outFace;
            move(inputFace, uv, gradient, outFace, end2D);
            pointCandidate = computeXYZ(surfaceMesh, end2D, outFace);
            pointCandidates.push_back(pointCandidate);
            directions.push_back(pointCandidate);
            if(distance < 0.1)
            {
                if(barycentricCoordinates(surfaceMesh,inputFace, targetPoint))
                {
                    inputPoint = targetPoint;
                    directions.push_back(inputPoint);
                    break;
                }
            }
            inputPoint = pointCandidate;
            inputFace = outFace;
            //std::cout<<"iteration "<<iteration << " dist "<<std::sqrt(CGAL::squared_length(targetPoint-inputPoint))<<"\n";
            distance  = std::sqrt(CGAL::squared_length(targetPoint-inputPoint));
            //directions.push_back(pointCandidate);
            iteration++;
            //parameters.deltaXY = std::max(0.0001 , parameters.deltaXY*0.99);
        }
        //inputPoint = directions.back();

    }
    int max_iteration = 1;
    CGAL::Surface_mesh<Point> surfaceMesh;
    OptimizerParameters parameters;

    std::vector<Point> directions;
    std::vector<Point> gradients;

    //
    std::vector<Point> debug_steps;
    bool debug = false;

    std::vector<Point> pointCandidates;
    std::vector<Point> debugPath;
    std::vector<Vector> debugNormal;

    std::vector<Point> debug_gradient;
    DEBUG_2D debug2D;
    std::vector<DEBUG_2D> debug2Ds;
};