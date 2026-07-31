#pragma once
#include "coord.h"
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

        if (debug) {
            debugPath.push_back(start);
            debugNormal.push_back(directionSurface);
        }

        while(remainingDistance > 0.)
        {
            Vector initialNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(currentFace,surfaceMesh);
            auto currentStart2D = computeUV(surfaceMesh, currentStart, currentFace);


            Mesh::Face_index nextFace;
            K::Point_2 intersectionPoint2D;
            intersection2D(currentFace, currentStart2D, end2D, nextFace, intersectionPoint2D);

            Vector nextNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(nextFace,surfaceMesh);

            Point intersectionPoint =  computeXYZ(surfaceMesh, intersectionPoint2D, currentFace);



            double distanceToIntersection = std::sqrt(CGAL::squared_length(intersectionPoint - currentStart));
            Vector direction3D = computeDirection3D(intersectionPoint, currentTarget, initialNormal, nextNormal);

            if( debug )
            {
                debugPath.push_back(intersectionPoint);
                debugNormal.push_back(direction3D);
            }


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

            if(bary)
            {
                target = step;
                outFace = nextFace;
                return ;
            }
            currentFace = nextFace;
            end2D = computeUV(surfaceMesh, currentTarget, currentFace);
            end2D = clampUV(end2D);
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

        if (barycentricCoordinates(surfaceMesh, faceStart, pointMoved)) {
            outFace = faceStart;
            return;
        }
        else
        {

            auto startPoint =  computeXYZ(surfaceMesh, start, faceStart);
            /*for(auto h : surfaceMesh.halfedges_around_face(surfaceMesh.halfedge(faceStart)))
            {
                if (barycentricCoordinates(surfaceMesh, faceStart, pointMoved)) {

                }
            }*/
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
                        const K::Point_2& start,
                        const K::Point_2& target,
                        Mesh::Face_index & faceEnd, K::Point_2& end )
    {
        K::Vector_2 dir = target - start;
        dir = normalize(dir);

        K::Point_2 outside = start + 100.0 * dir;

        std::optional<K::Point_2> bestPoint;
        double bestDist = std::numeric_limits<double>::max();

        K::Ray_2 ray(outside, start);
        float epsilonLength = 10e-5;
        constexpr double eps2 = 1e-20;


        for (const auto& h : CGAL::halfedges_around_face(surfaceMesh.halfedge(faceStart), surfaceMesh)) {
            const auto v1 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.source(h)), faceStart);
            const auto v2 = computeUV(surfaceMesh, surfaceMesh.point(surfaceMesh.target(h)), faceStart);
            const auto v1v2 = v1 - v2;
            const K::Segment_2 seg(v1  , v2);


            auto result = CGAL::intersection(ray, seg);
            if (!result)
            {
                // check if start it belong to the segment ?
                if (seg.has_on(start))
                {
                    // compute t parameter closest to target
                    // then the correct point
                    K::Vector_2 edge = v2 - v1;
                    K::Vector_2 move = target - start;

                    double len2 = CGAL::to_double(edge.squared_length());
                    if (len2 < 1e-20)
                        continue;

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

                    double d = CGAL::to_double(CGAL::squared_distance(ray.source(), p));

                    if (d < bestDist)
                    {
                        bestDist = d;
                        bestPoint = p;
                        faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
                    }
                }

                continue;
            }

            if (const auto* p0 = std::get_if<K::Point_2>(&*result))
            {
                if (CGAL::squared_distance(*p0, ray.source()) < eps2)
                {
                    continue;
                }

                // one intersection
                double d = CGAL::to_double(CGAL::squared_distance(ray.source(), *p0));
                if (d < bestDist)
                {
                    bestDist = d;
                    bestPoint = *p0;
                    faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
                }
            }
            else if (const auto* s = std::get_if<K::Segment_2>(&*result))
            {
                // overlapping segment
                auto d0 = CGAL::to_double(CGAL::squared_distance(ray.source(), s->source()));
                auto d1 = CGAL::to_double(CGAL::squared_distance(ray.source(), s->target()));

                K::Point_2 p = (d0 > d1) ? s->source() : s->target();

                double d = std::min(d0, d1);   // first point encountered

                if (d < bestDist)
                {
                    bestDist = d;
                    bestPoint = p;
                    faceEnd = surfaceMesh.face(surfaceMesh.opposite(h));
                }
            }
        }

        if(bestPoint)
        {
            end = bestPoint.value();
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
    void optimize(Point &inputPoint, const Point &targetPoint, Mesh::Face_index &inputFace) {

        int iteration = 0;
        debug_steps.clear();
        while (iteration < max_iteration && std::sqrt(CGAL::squared_length(targetPoint-inputPoint)) > 10e-3) {

            auto uv = computeUV(surfaceMesh, inputPoint, inputFace);

            auto getEnergy = [&](K::Point_2 delta) {

              /*  auto target = computeXYZ(surfaceMesh, uv+(delta-CGAL::ORIGIN), inputFace);
                directions.push_back(target);
                auto totalDistance = std::sqrt(CGAL::squared_length(target- inputPoint));
                Mesh::Face_index outFace;
                moveOnSurface(inputFace,
                              inputPoint,
                                   target,
                                   outFace,
                                   totalDistance);
            */
                Mesh::Face_index outFace;
                Point target;
                moveOnSurface(inputFace, uv, delta-CGAL::ORIGIN, target, outFace);
                gradients.push_back(target);
                return evaluate(targetPoint, target);
            };
            K::Vector_2 gradient = K::Vector_2(
                    getEnergy(K::Point_2(-parameters.deltaXY, 0)) - getEnergy(K::Point_2(parameters.deltaXY, 0)),
                    getEnergy(K::Point_2(0, -parameters.deltaXY)) - getEnergy(K::Point_2(0, parameters.deltaXY)));

            double gradient_step_sq_len = std::sqrt(CGAL::squared_length(gradient));
            std::cout<<"gradient energy "<<gradient_step_sq_len<<"\n";

            if(gradient_step_sq_len > parameters.deltaXY)
                gradient = parameters.deltaXY / gradient_step_sq_len * gradient;

            auto pointCandidate = computeXYZ(surfaceMesh, uv+gradient, inputFace);
            Mesh::Face_index outFace;
            auto totalDistance = std::sqrt(CGAL::squared_length(pointCandidate- inputPoint));
            moveOnSurface(inputFace,
                          inputPoint,
                          pointCandidate,
                          outFace,
                          totalDistance);



            inputPoint = pointCandidate;
            inputFace = outFace;
            std::cout<<"iteration "<<iteration << " dist "<<std::sqrt(CGAL::squared_length(targetPoint-inputPoint))<<"\n";
            directions.push_back(pointCandidate);
            iteration++;
            parameters.deltaXY = std::max(0.0001 , parameters.deltaXY*0.99);
        }
        //inputPoint = debug_steps.back();

    }
    int max_iteration = 1;
    CGAL::Surface_mesh<Point> surfaceMesh;
    OptimizerParameters parameters;

    std::vector<Point> directions;
    std::vector<Point> gradients;

    //
    std::vector<Point> debug_steps;
    bool debug = false;

    std::vector<Point> debugPath;
    std::vector<Vector> debugNormal;
};