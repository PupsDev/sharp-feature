#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <CGAL/Ray_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Point_2.h>

struct OptimizerParameters
{
    float lineQuadricWeight = 0.001f;
    float deltaXY = 0.025f;


};

using Base = std::pair<Vector, Vector>;

struct SurfacePatch
{
    Base base;
    Base orthoBase;
    Vector normal;
    Mesh::Face_index face;

    K::Point_2 uvs; // coordinate of the point in the base
};

class Optimizer
{
public:


    Optimizer()
    {

    }
    Optimizer(const CGAL::Surface_mesh<Point>& inputMesh)
    {
        surfaceMesh = inputMesh;
    }
    Eigen::Matrix4d computeLineQuadric(const Eigen::Vector3d& vi,
                                       const Eigen::Vector3d& ni)
    {
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

        return Qx + Qy;
    }
    Point computePlaneIntersection(const std::vector<Plane>& planes, const Point& inputPoint,const Vector & inputNormal )
    {
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
        Eigen::Matrix4d quadric = Eigen::Matrix4d::Zero();
        std::vector<Point > pts;
        for(const auto& plane : planes)
        {
            quadric += computeQuadricFromPlane(plane);
            pts.push_back(plane.first);
        }
        Eigen::Vector3d point(inputPoint.x(),inputPoint.y(),inputPoint.z());
        Eigen::Vector3d normal(inputNormal.x(),inputNormal.y(),inputNormal.z());

        auto lineQuadric = computeLineQuadric(point, normal);
        quadric += parameters.lineQuadricWeight* lineQuadric;
        auto eigenPoint = optimizeQuadric(quadric);

        auto cgalPoint =  Point (eigenPoint[0], eigenPoint[1], eigenPoint[2]);

        auto center =  CGAL::centroid(pts[0] ,pts[1] ,pts[2]);

        //std::cout<<"cgalPoint "<<cgalPoint<<"\n";
        return cgalPoint;// + Vector(center.x(), center.y(), center.z());;
    }

    Point getUVs( Mesh::Face_index f,const Vector& b1, const Vector& b2 , const Point& worldPosition )
    {
        const auto& h0 = surfaceMesh.halfedge(f);
        const auto& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        const auto& normal = CGAL::cross_product(b1,b2);
        Vector d = worldPosition - p0;

        double u = d * b1; // dot product
        double v = d * b2;
        double w = d * normal;

        return Point(u, v, w);
    }
    std::pair<Vector,Vector> orthoBasis(Mesh::Face_index f) {
        const auto h0 = surfaceMesh.halfedge(f);
        const auto h1 = next(h0, surfaceMesh);
        const auto h2 = next(h1, surfaceMesh);

        const Point& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        const Point& p1 = surfaceMesh.point(target(h1, surfaceMesh));
        const Point& p2 = surfaceMesh.point(target(h2, surfaceMesh));

        Vector e1 = p1 - p0;
        Vector e2 = p2 - p0;

        // normal (unit)
        Vector n = CGAL::cross_product(e1, e2);
        n = n / std::sqrt(n.squared_length());

        // first tangent (unit)
        Vector t1 = e1 / std::sqrt(e1.squared_length());

        // make second tangent orthogonal to both n and t1
        Vector t2 = CGAL::cross_product(n, t1);
        t2 = t2 / std::sqrt(t2.squared_length());

        return std::make_pair(t1, t2);
    }

    Point worldToTriangle(
            Mesh::Face_index f,
            const Vector& b1,
            const Vector& b2,
            const Point& worldPosition)
    {
        auto h0 = surfaceMesh.halfedge(f);
        const auto& p0 = surfaceMesh.point(target(h0, surfaceMesh));

        Vector d = worldPosition - p0;

        double a = b1 * b1;
        double b = b1 * b2;
        double c = b2 * b2;
        double d1 = d * b1;
        double d2 = d * b2;

        double det = a * c - b * b;

        double u = (d1 * c - d2 * b) / det;
        double v = (d2 * a - d1 * b) / det;

        return Point(u, v,0.);
    }
    Point worldToTriangle(
            Mesh::Face_index f,
            std::pair<Vector ,Vector> & base,
            const Point& worldPosition)
    {
        const auto& h0 = surfaceMesh.halfedge(f);
        const auto& h1 = next(h0, surfaceMesh);
        const auto& h2 = next(h1, surfaceMesh);

        const auto& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        const auto& p1 = surfaceMesh.point(target(h1, surfaceMesh));
        const auto& p2 = surfaceMesh.point(target(h2, surfaceMesh));

        Vector b1 = p1-p0;
        Vector b2 = p2-p0;
        base.first = b1;
        base.second = b2;
        return worldToTriangle(f, b1,b2, worldPosition);
    }
    Point triangleToWorldOrtho(Mesh::Face_index f, const Vector& b1, const Vector& b2, double u, double v)
    {
        const auto& h0 = surfaceMesh.halfedge(f);
        const auto& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        return p0 + b1 * u +b2 * v;
    }
    Point triangleToWorldOrtho(Mesh::Face_index f, const std::pair<Vector ,Vector> & base, double u, double v)
    {
        return triangleToWorldOrtho(f, base.first, base.second, u,v);
    }

    Point triangleToWorld(Mesh::Face_index f, double u, double v)
    {
        const auto& h0 = surfaceMesh.halfedge(f);
        const auto& h1 = next(h0, surfaceMesh);
        const auto& h2 = next(h1, surfaceMesh);

        const auto& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        const auto& p1 = surfaceMesh.point(target(h1, surfaceMesh));
        const auto& p2 = surfaceMesh.point(target(h2, surfaceMesh));

        return p0 + (p1 - p0) * u + (p2 - p0) * v;
    }
    Vector computeQ (const Vector& d,
                     const Vector& n1,
                     const Vector& n2)
    {
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
    }

    auto evaluate (const Point& targetPoint, const Point& inputPoint)
    {

        return CGAL::squared_distance(targetPoint, inputPoint) ;
    }
    Point intersect(Point nextUV)
    {
        K::Point_2 p1(currentUV.x(), currentUV.y());
        K::Point_2 p2(nextUV.x(), nextUV.y());
        K::Ray_2 ray(p1, p2-p1);
        double mindist = FLT_MAX;
        Mesh::Face_index nextFace;
        Point currentPoint;

        auto convert2D = [&](const Point& mp_in)
        {
            auto uvs = worldToTriangle(currentF, currentBase.first, currentBase.second, mp_in);
            return K::Point_2(uvs.x(),uvs.y());
        };
        for (const auto& h : CGAL::halfedges_around_face(surfaceMesh.halfedge(currentF), surfaceMesh)) {
            auto v1 = convert2D(surfaceMesh.point(source(h, surfaceMesh)));
            auto v2 = convert2D(surfaceMesh.point(target(h, surfaceMesh)));
            auto v1v2 = v1 - v2;
            K::Segment_2 seg(v1 + 10e-5 * v1v2, v2 - 10e-5 * v1v2);
            auto result = CGAL::intersection(ray, seg);
            if (result) {
                if (const K::Point_2 *p = std::get_if<K::Point_2>(&*result)) {
                    auto vec = *p - p1;
                    if(vec.squared_length() < mindist)
                    {
                        currentPoint =  triangleToWorld(currentF, p->x(), p->y());;
                        mindist= vec.squared_length();
                        nextFace = surfaceMesh.face(opposite(h, surfaceMesh));
                    }
                }
            }
        }
        currentF = nextFace;
        currentBase = orthoBasis(currentF);
        currentUV = getUVs(currentF, currentBase.first, currentBase.second, currentPoint);
        return currentPoint;
    }
    K::Point_2 intersect(const SurfacePatch& surfacePatch, const   K::Point_2& origin, const K::Point_2 end, SurfacePatch& outSurfacePatch)
    {
        /*K::Point_2 p1(currentUV.x(), currentUV.y());
        K::Point_2 p2(nextUV.x(), nextUV.y());
        K::Ray_2 ray(p1, p2-p1);*/


        double mindist = FLT_MAX;
        Mesh::Face_index nextFace = surfacePatch.face;
        auto convert2D = [&](const Point& mp_in)
        {
            auto uvs = worldToTriangle(surfacePatch.face, surfacePatch.base.first, surfacePatch.base.second, mp_in);

            return K::Point_2(uvs.x(),uvs.y());
        };
        K::Point_2 currentPoint;
        /* Point currentPoint;


         K::Point_2 p1 = convert2D(origin);
         Point endPoint = Point(origin + 10.* direction);
         K::Point_2 p2 = convert2D(endPoint);

         auto debugpt = triangleToWorldOrtho(surfacePatch.face, surfacePatch.base, p2.x(),p2.y());
         pointDebug2.push_back(debugpt);*/

        std::vector<K::Point_2> pointDebug2D;
        K::Point_2 p1 = origin;
        K::Point_2 p2 =  K::Point_2(origin + 10.*(end-origin));
        //K::Point_2 p2 = endPoint;

        for (auto h : CGAL::halfedges_around_face(
                surfaceMesh.halfedge(surfacePatch.face), surfaceMesh))
        {
            auto vd = target(h, surfaceMesh);
            const auto& point = surfaceMesh.point(vd);

            // use point.x(), point.y(), point.z()
            pointDebug2D.push_back(convert2D(point));
            auto p2d = convert2D(point);
            Point globalStartPoint =  triangleToWorldOrtho(surfacePatch.face, surfacePatch.base, p2d.x(), p2d.y());
            pointDebug.push_back(globalStartPoint);
        }

        /*Point globalStartPoint =  triangleToWorldOrtho(surfacePatch.face, surfacePatch.base, end.x(), end.y());
        pointDebug.push_back(globalStartPoint);*/


        K::Ray_2 ray(p2, p1-p2);
        K::Point_2 intersection;

        pointDebug2D.push_back(p1);
        //pointDebug2D.push_back(p2);


        for (const auto& h : CGAL::halfedges_around_face(surfaceMesh.halfedge(surfacePatch.face), surfaceMesh)) {
            auto v1 = convert2D(surfaceMesh.point(source(h, surfaceMesh)));
            auto v2 = convert2D(surfaceMesh.point(target(h, surfaceMesh)));

            //pointDebug2D.push_back(v1);
            //pointDebug2D.push_back(v2);

            auto v1v2 = v1 - v2;
            K::Segment_2 seg(v1 + 10e-5 * v1v2, v2 - 10e-5 * v1v2);
            auto result = CGAL::intersection(ray, seg);
            if (result) {
                if (const K::Point_2 *p = std::get_if<K::Point_2>(&*result)) {
                    //pointDebug2D.push_back(*p);
                    auto vec = *p - p2;
                    if(vec.squared_length() < mindist)
                    {
                        intersection = *p;
                        //currentPoint =  triangleToWorld(surfacePatch.face, p->x(), p->y());;
                        mindist= vec.squared_length();
                        nextFace = surfaceMesh.face(opposite(h, surfaceMesh));
                    }
                }
            }
        }
        outSurfacePatch.face = nextFace;
        currentPoint =  intersection;//triangleToWorld(nextFace, intersection.x(), intersection.y());
        pointDebug2.clear();
        for(auto& p : pointDebug2D)
        {
            pointDebug2.push_back(Point(p.x(), p.y(), 0.));
        }
        return currentPoint;
    }
    Point moveOnSurface(double u, double v)
    {

        Point nextPoint;
        float epsilon = 0.01f;
        Point currentPoint = triangleToWorld(currentF, currentUV.x(), currentUV.y());
        Vector direction;// = currentPoint-nextPoint;

        if( ( u > 1.+epsilon || (u < -0.-epsilon)) ||
            ( v > 1.+epsilon || (v < -0.-epsilon) ))
        {
            K::Point_2 p1(currentUV.x(), currentUV.y());
            auto nextUV = Point(u, v, 1.);

            Mesh::Face_index nextFace;
            Mesh::Face_index startFace = currentF;
            nextPoint = intersect(nextUV);
            direction = currentPoint-nextPoint;
            nextFace = currentF;
            if(nextFace == Mesh::null_face())
            {
                std::cout<<"Found border \n";
                startFace = nextFace;
            }
            direction = computeQ(direction, CGAL::Polygon_mesh_processing::compute_face_normal(startFace,surfaceMesh),
                                 CGAL::Polygon_mesh_processing::compute_face_normal(nextFace,surfaceMesh));
        }
        else
        {
            nextPoint = triangleToWorldOrtho(currentF,currentBase.first,currentBase.second,u,v);
        }
        return nextPoint;
    }
    Point optimizeLocal( const Point& target)
    {
        std::vector<double> errors;
        std::vector<Point> candidates;

        double u = currentUV.x();
        double v = currentUV.y();

        auto pts0 = moveOnSurface(u+parameters.deltaXY, v);
        candidates.push_back(pts0);
        errors.push_back(evaluate(target, pts0));

        auto pts1 = moveOnSurface(u-parameters.deltaXY, v);
        candidates.push_back(pts1);
        errors.push_back(evaluate(target, pts1));

        auto pts2 = moveOnSurface(u, v+parameters.deltaXY);
        candidates.push_back(pts2);
        errors.push_back(evaluate(target, pts2));

        auto pts3 = moveOnSurface(u, v-parameters.deltaXY);
        candidates.push_back(pts3);
        errors.push_back(evaluate(target, pts3));

        auto minError = std::min_element(errors.begin(), errors.end());
        auto pointCandidate = candidates[std::distance( errors.begin(),minError)];


        auto pcl = polyscope::registerPointCloud("candidates", candidates);
        pcl->resetTransform();


        return pointCandidate;
    }

    std::pair<Vector, Vector> getOrthoBasis(Mesh::Face_index f, const Point& worldPosition ) {
        const auto h0 = surfaceMesh.halfedge(f);
        const auto h1 = next(h0, surfaceMesh);
        const auto h2 = next(h1, surfaceMesh);

        const Point& p0 = surfaceMesh.point(target(h0, surfaceMesh));
        const Point& p1 = surfaceMesh.point(target(h1, surfaceMesh));
        const Point& p2 = surfaceMesh.point(target(h2, surfaceMesh));

        Vector e1 = p1 - p0;
        Vector e2 = p2 - p0;

        // normal (unit)
        Vector n = CGAL::cross_product(e1, e2);
        n = n / std::sqrt(n.squared_length());

        // first tangent (unit)
        Vector t1 = e1 / std::sqrt(e1.squared_length());

        // make second tangent orthogonal to both n and t1
        Vector t2 = CGAL::cross_product(n, t1);
        t2 = t2 / std::sqrt(t2.squared_length());

        return std::make_pair(t1, t2);
    };

    bool checkBoundaries(const glm::vec3& vec)
    {
        float eps = 10e-5;
        return (vec[0]  < - eps  ) ||
               (vec[1]  < - eps  ) ||
               ( (vec[1]+vec[0]) > (1+eps)    );
    }
    bool checkBoundaries(const K::Point_2& vec)
    {
        float eps = 10e-5;
        return (vec[0]  < - eps  ) ||
               (vec[1]  < - eps  ) ||
               ( (vec[1]+vec[0]) > (1+eps)    );
    }

    void updatePatch(SurfacePatch& patch, const Point& inputPoint)
    {
        Point uvs =  worldToTriangle(
                patch.face,
                patch.base,
                inputPoint);
        patch.orthoBase = orthoBasis(patch.face);;
        patch.normal = CGAL::Polygon_mesh_processing::compute_face_normal(patch.face,surfaceMesh);
        patch.uvs = K::Point_2(uvs.x(), uvs.y());
    }
    Point test(SurfacePatch& patchOutput, K::Point_2 startPoint2D, K::Point_2 endPoint2D, float totalDistance)
    {
        //float totalDistance;
        /*K::Point_2 startPoint2D =  patchOutput.uvs;
        K::Point_2 endPoint2D = K::Point_2(vec2.x, vec2.y);*/

        Point globalStartPoint =  triangleToWorldOrtho(patchOutput.face, patchOutput.base, startPoint2D.x(), startPoint2D.y());
        Point globalEndPoint =  triangleToWorldOrtho(patchOutput.face, patchOutput.base, endPoint2D.x(), endPoint2D.y());
        debugStep.push_back(globalStartPoint);

        debugStep.push_back(globalEndPoint);

        while(checkBoundaries(endPoint2D) && totalDistance > 0.)
        {

            auto initialNormal =  CGAL::Polygon_mesh_processing::compute_face_normal(patchOutput.face,surfaceMesh);
            SurfacePatch outSurfacePatch;
            auto intersectionPoint2D = intersect(patchOutput, startPoint2D, endPoint2D, outSurfacePatch);
            double x = intersectionPoint2D.x();
            double y = intersectionPoint2D.y();

            if (x > 1.0 - 1e-5)
                x = 1.0 - 1e-5;

            if (y > 1.0 - 1e-5)
                y = 1.0 - 1e-5;

            if (x < 1e-5)
                x = 1e-5;

            if (y < 1e-5)
                y = 1e-5;

            intersectionPoint2D = K::Point_2(x, y);
            Point intersectionPoint = triangleToWorldOrtho(patchOutput.face, patchOutput.base, intersectionPoint2D.x(),
                                                           intersectionPoint2D.y());

            if(outSurfacePatch.face == Mesh::null_face())
            {
                updatePatch(patchOutput, intersectionPoint);
                break;
            }

            patchOutput = outSurfacePatch;
            updatePatch(patchOutput, intersectionPoint);

            auto nextNormal = CGAL::Polygon_mesh_processing::compute_face_normal(patchOutput.face, surfaceMesh);

            Vector direction3D = globalEndPoint - globalStartPoint;
            direction3D /= std::sqrt(
                    CGAL::to_double(CGAL::squared_length(direction3D)));

            if(CGAL::squared_length(nextNormal-initialNormal) > 1e-5)
            {
                direction3D = computeQ(direction3D, initialNormal,
                                       nextNormal);
                direction3D /= std::sqrt(
                        CGAL::to_double(CGAL::squared_length(direction3D)));
            }


            float localDistance = std::sqrt(CGAL::to_double(CGAL::squared_length(intersectionPoint - globalStartPoint)));
            float remaining = totalDistance - localDistance;
            totalDistance -= (localDistance);

            globalEndPoint = intersectionPoint + remaining * direction3D;
            globalStartPoint = intersectionPoint;

            float localDistance2 = std::sqrt(CGAL::to_double(CGAL::squared_length(globalEndPoint - intersectionPoint)));

            debugStep.push_back(globalStartPoint);

            debugStep2.push_back(globalEndPoint);

            Point uvStart = worldToTriangle(
                    patchOutput.face,
                    patchOutput.base,
                    intersectionPoint);
            startPoint2D =  K::Point_2(uvStart.x(), uvStart.y());

            Point uvs = worldToTriangle(
                    patchOutput.face,
                    patchOutput.base,
                    globalEndPoint);
            endPoint2D = K::Point_2(uvs.x(), uvs.y());
            patchOutput.uvs = endPoint2D;
        }
        return globalEndPoint;
    }
    SurfacePatch moveDelta(const SurfacePatch& patch,
                           double u, double v) {
        SurfacePatch patchOutput = patch;
        glm::mat3 transform(
                patch.base.first[0], patch.base.first[1], patch.base.first[2],
                patch.base.second[0], patch.base.second[1], patch.base.second[2],
                patch.normal[0], patch.normal[1], patch.normal[2]
        );
        auto inverted = glm::inverse(transform);

        glm::mat3 transformOrtho(
                patch.orthoBase.first[0], patch.orthoBase.first[1], patch.orthoBase.first[2],
                patch.orthoBase.second[0], patch.orthoBase.second[1], patch.orthoBase.second[2],
                patch.normal[0], patch.normal[1], patch.normal[2]
        );
        auto invertedOrtho = glm::inverse(transformOrtho);

        glm::vec3 uvPosition = glm::vec3(patch.uvs.x(), patch.uvs.y(), 0.);
        glm::vec3 orthoPosition = invertedOrtho * transform * uvPosition;
        orthoPosition += glm::vec3(u, v, 0.);

        glm::vec3 displacedPosition = uvPosition + glm::vec3(u, v, 0.);// inverted * transformOrtho * orthoPosition;


        K::Point_2 startPoint2D = K::Point_2( patchOutput.uvs.x(),  patchOutput.uvs.y());
        K::Point_2 endPoint2D = K::Point_2(displacedPosition.x, displacedPosition.y);
        Point firstStartPoint = triangleToWorldOrtho(patchOutput.face, patchOutput.base, patchOutput.uvs.x(),
                                                     patchOutput.uvs.y());
        Point firstEndPoint = triangleToWorldOrtho(patchOutput.face, patchOutput.base, endPoint2D.x(), endPoint2D.y());
        float totalDistance =    std::sqrt(CGAL::to_double(CGAL::squared_length(firstEndPoint - firstStartPoint)));


        //Point currentStart = firstStartPoint;
        if(checkBoundaries(displacedPosition)) {
            Point out = test(patchOutput, startPoint2D, endPoint2D, totalDistance);
        }
        else
        {
            patchOutput.uvs = K::Point_2(displacedPosition.x, displacedPosition.y);
        }
        return patchOutput;
    }

    void optimize(Point& inputPoint, const Point& targetPoint, const  Mesh::Face_index& inputFace)
    {
        SurfacePatch patch;
        patch.face = inputFace;
        updatePatch(patch, inputPoint);

        int iteration = 0;

        std::vector<K::Point_2> deltas{K::Point_2(parameters.deltaXY, 0.),
                                       K::Point_2(0., parameters.deltaXY),
                                       K::Point_2(-parameters.deltaXY, 0.),
                                       K::Point_2(0., -parameters.deltaXY)
        };
        directions.clear();
        while(iteration < max_iteration)
        {
            std::vector<Point> pts;
            std::vector<SurfacePatch> patches;
            std::vector<double> errors;

            float distance = evaluate(inputPoint,targetPoint);
            if(distance < 10e-5)
            {
                return;
            }


            auto getEnergy = [&](K::Point_2 delta)
            {
                SurfacePatch patchCopy = patch;
                patchCopy = moveDelta(patchCopy, delta.x(), delta.y());
                Point globalPoint =  triangleToWorldOrtho(patchCopy.face, patchCopy.base, patchCopy.uvs.x(), patchCopy.uvs.y());
                return evaluate(targetPoint, globalPoint);
            };
            K::Vector_2 gradient = K::Vector_2 (
                    getEnergy(K::Point_2(-parameters.deltaXY, 0))-getEnergy(K::Point_2(parameters.deltaXY, 0)),
                    getEnergy(K::Point_2(0, -parameters.deltaXY))-getEnergy(K::Point_2(0, parameters.deltaXY)));

            double gradient_step_sq_len = std::sqrt(CGAL::squared_length(gradient));

            gradient = parameters.deltaXY / CGAL::sqrt(gradient_step_sq_len) * gradient;

            SurfacePatch patchCopy = patch;
            patchCopy = moveDelta(patchCopy, gradient.x(), gradient.y());
            auto pointCandidate =  triangleToWorldOrtho(patchCopy.face, patchCopy.base, patchCopy.uvs.x(), patchCopy.uvs.y());

            directions.push_back(pointCandidate);
            updatePatch(patch, pointCandidate);

            iteration ++;
        }
        inputPoint=directions.back();

    }
    void optimize_OK(Point& inputPoint, const Point& targetPoint, const  Mesh::Face_index& inputFace)
    {
        SurfacePatch patch;
        patch.face = inputFace;
        updatePatch(patch, inputPoint);

        int iteration = 0;

        std::vector<K::Point_2> deltas{K::Point_2(parameters.deltaXY, 0.),
                                       K::Point_2(0., parameters.deltaXY),
                                       K::Point_2(-parameters.deltaXY, 0.),
                                       K::Point_2(0., -parameters.deltaXY)
        };
        directions.clear();
        while(iteration < max_iteration)
        {
            std::vector<Point> pts;
            std::vector<SurfacePatch> patches;
            std::vector<double> errors;

            float distance = evaluate(inputPoint,targetPoint);
            if(distance < 10e-5)
            {
                return;
            }

            for(int i = 0; i <  deltas.size(); i++)
            {
                //std::cout<<"i "<<i<<"\n";
                SurfacePatch patchCopy = patch;
                patchCopy = moveDelta(patchCopy, deltas[i].x(), deltas[i].y());
                Point globalPoint =  triangleToWorldOrtho(patchCopy.face, patchCopy.base, patchCopy.uvs.x(), patchCopy.uvs.y());
                pts.push_back(globalPoint);
                patches.push_back(patchCopy);
                errors.push_back(evaluate(targetPoint, globalPoint));
            }

            auto minError = std::min_element(errors.begin(), errors.end());
            auto pointCandidate = pts[std::distance( errors.begin(),minError)];
            directions.push_back(pointCandidate);
            patch =  patches[std::distance( errors.begin(),minError)];
            updatePatch(patch, pointCandidate);

            iteration ++;
        }
        inputPoint=directions.back();

    }

    int max_iteration = 1;
    CGAL::Surface_mesh<Point> surfaceMesh;
    OptimizerParameters parameters;

    std::pair<Vector, Vector> currentBase;
    Point currentUV;
    Mesh::Face_index currentF;


    std::vector<Point> directions;

    std::vector<Point> startPoints;
    std::vector<Point> interSectionPoints;
    std::vector<Point> endPoints;



    std::vector<Point> pointDebug;
    std::vector<Point> pointDebug2;


    std::vector<Point> debugStep;
    std::vector<Point> debugStep2;
};