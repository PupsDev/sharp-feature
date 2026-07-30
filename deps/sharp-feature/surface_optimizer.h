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
    void intersection2D(const Mesh::Face_index & faceStart,
                        const K::Vector_2& start,
                        const K::Vector_2& target,
                        Mesh::Face_index & faceEnd, K::Vector_2& end )
    {
        auto iFace = face;
        double x = localP0.x();
        double y = localP0.y();
        float eps = 10e-3;
        if( x < eps) x = eps;
        if( y < eps) y = eps;

        if( x > 1-eps) x = 1-eps;
        if( y > 1-eps) y = 1-eps;
        localP0 = K::Point_2 (x,y);


        std::vector<K::Segment_2 > segments;
        auto rayDirect = localP - localP0;
        bool foudD = false;

    }
    auto evaluate (const Point& targetPoint, const Point& inputPoint)
    {

        return CGAL::squared_distance(targetPoint, inputPoint) ;
    }
    void optimize(Point &inputPoint, const Point &targetPoint, const Mesh::Face_index &inputFace) {

        /*int iteration = 0;
        debug_steps.clear();
        while (iteration < max_iteration) {
            std::vector <Point> pts;
            std::vector <SurfacePatch> patches;
            std::vector<double> errors;

            float distance = evaluate(inputPoint, targetPoint);
            if (distance < 10e-5) {
                return;
            }

            auto getEnergy = [&](K::Point_2 delta) {
                SurfacePatch patchCopy = patch;
                patchCopy = moveDelta(patchCopy, delta.x(), delta.y());
                Point globalPoint = triangleToWorldOrtho(patchCopy.face, patchCopy.base, patchCopy.uvs.x(),
                                                         patchCopy.uvs.y());
                return evaluate(targetPoint, globalPoint);
            };
            K::Vector_2 gradient = K::Vector_2(
                    getEnergy(K::Point_2(-parameters.deltaXY, 0)) - getEnergy(K::Point_2(parameters.deltaXY, 0)),
                    getEnergy(K::Point_2(0, -parameters.deltaXY)) - getEnergy(K::Point_2(0, parameters.deltaXY)));

            double gradient_step_sq_len = std::sqrt(CGAL::squared_length(gradient));

            gradient = parameters.deltaXY / CGAL::sqrt(gradient_step_sq_len) * gradient;

            SurfacePatch patchCopy = patch;
            patchCopy = moveDelta(patchCopy, gradient.x(), gradient.y());
            auto pointCandidate = triangleToWorldOrtho(patchCopy.face, patchCopy.base, patchCopy.uvs.x(),
                                                       patchCopy.uvs.y());
            debug_steps.push_back(pointCandidate);
            updatePatch(patch, pointCandidate);

            iteration++;
        }
        inputPoint = debug_steps.back();*/

    }
    int max_iteration = 1;
    CGAL::Surface_mesh<Point> surfaceMesh;
    OptimizerParameters parameters;

    //
    std::vector<Point> debug_steps;
};