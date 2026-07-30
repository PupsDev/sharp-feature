#pragma once
#include "nvt_helper.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
template<typename UserType>
class NVTGeodesicFibonacci
{
public:
    NVTGeodesicFibonacci() = default;
    struct DEBUG
    {
        std::vector<Point> inputs;
        std::vector<Point> outputs;
        void display()
        {
            auto pcl = polyscope::registerPointCloud("inputs", inputs);
            pcl->setPointRadius(0.0001);
            pcl->resetTransform();

            auto pclOut = polyscope::registerPointCloud("outputs", outputs);
            pclOut->setPointRadius(0.0001);
            pclOut->resetTransform();
        }
    };
    struct SampleData
    {
        Point point;
        Mesh::Face_index faceIndex;
        double r; // for fibo
    };
    void set(const Mesh& sm)
    {
        surfaceMesh = sm;
    }
    std::vector<UserType> computeNVTGeodesic(const FeatureMetricParameters<UserType>& parameters)
    {


        int vertexIDGlobal;
        std::vector<std::vector<Vector >> normalPerVertex(surfaceMesh.num_vertices());

        samplingPerVertex.resize(surfaceMesh.num_vertices());
        samplingInitialPerVertex.resize(surfaceMesh.num_vertices());

        std::vector<std::vector<double>> geodistance(surfaceMesh.num_vertices());

        for (int i = 0; i < surfaceMesh.num_vertices(); i++) {
            auto vertexID = Mesh::Vertex_index(i);
            auto samples = computeInitialFiboSampling(vertexID, parameters);

            std::vector<Point> samplePoints;
            for (auto &sample: samples) samplePoints.push_back(sample.point);
            samplingInitialPerVertex[i] = samplePoints;

            std::vector<Point> intersection3D;
            std::vector<int> ids;
            int id = 0;
            for (auto &sample: samples) {
                auto direction = sample.point - surfaceMesh.point(vertexID);
                Mesh::Face_index face;
                intersection3D.push_back(movePointGeodesic( sample, direction, sample.r, face));
                geodistance[i].push_back(sample.r);
                normalPerVertex[i].push_back(CGAL::Polygon_mesh_processing::compute_face_normal(face, surfaceMesh));
                id++;
            }
            samplingPerVertex[i] = intersection3D;

        }

        vertexMetrics.clear();
        auto nId = 0;
        for (const auto &normals: normalPerVertex) {
            Eigen::Matrix3d tensor;
            tensor.setZero();
            int sampleID = 0;
            double rTotal = 0.;
            for (const auto &n: normals) {
                auto mat = computeTensor(n);
                if (parameters.useGeodesicWeight) {
                    double w = exp(-geodistance[nId][sampleID] / parameters.sigma);
                    mat = w * mat;
                    rTotal += w;
                }
                tensor += mat;
                sampleID++;
            }
            nId++;
            if (parameters.useGeodesicWeight) {
                tensor = tensor /= rTotal;
            }
            vertexMetrics.push_back(solveTensor(tensor));
        }
        return vertexMetrics;
    }
    std::vector<SampleData> computeInitialFiboSampling( Mesh::Vertex_index vertexID, const FeatureMetricParameters<UserType>& parameters)
    {
        auto angleBetweenHalfedges = [&](CGAL::SM_Halfedge_index h1, CGAL::SM_Halfedge_index h2) -> double
        {
            auto v1 = surfaceMesh.point(surfaceMesh.target(h1)) - surfaceMesh.point(surfaceMesh.source(h1));
            auto v2 = surfaceMesh.point(surfaceMesh.target(h2)) - surfaceMesh.point(surfaceMesh.source(h2));

            double cosTheta =
                    CGAL::to_double(v1 * v2) /
                    std::sqrt(CGAL::to_double(v1.squared_length() * v2.squared_length()));

            cosTheta = std::clamp(cosTheta, -1.0, 1.0);

            return std::acos(cosTheta); // radians
        };
        struct TriangleData
        {
            Mesh::Face_index faceIndex;
            Vector v0, v1;
            Point p0;
        };


        std::vector<float> planeAngles;
        std::vector<float> planeAnglesSum;
        std::vector<TriangleData> triangles;


        std::vector<SampleData> samples;

        std::vector<Point > ptsDe;
        std::vector<Vector > vecDe;

        //std::cout<<"VertexID "<<vertexID<<"\n";
        float totalSum = 0.;
        for (auto hi : CGAL::halfedges_around_source(vertexID, surfaceMesh))
        {
            //auto hi2 = surfaceMesh.opposite(surfaceMesh.prev(hi));
            auto hi2 = surfaceMesh.next(surfaceMesh.opposite(hi));
            float angle = angleBetweenHalfedges(hi, hi2);
            planeAngles.push_back(angle);
            TriangleData triangle;
            triangle.faceIndex = surfaceMesh.face(hi2);
            triangle.v0 = surfaceMesh.point(surfaceMesh.target(hi))- surfaceMesh.point(surfaceMesh.source(hi));
            triangle.v1 = surfaceMesh.point(surfaceMesh.target(hi2))- surfaceMesh.point(surfaceMesh.source(hi));
            triangle.p0 = surfaceMesh.point(surfaceMesh.source(hi));
            triangles.push_back(triangle);
            planeAnglesSum.push_back(totalSum);
            totalSum += angle;
        }
        planeAnglesSum.push_back(totalSum);
        auto planeAnglesSumNormalized = planeAnglesSum;
        for(auto& angle : planeAnglesSumNormalized) angle/= totalSum;
        // found triangle
        auto normalize = [](const Vector& x) {
            return x / std::sqrt(CGAL::to_double(x.squared_length()));
        };
        auto rotateTowards = [&](const Vector& u, const Vector& axis, double alpha)
        {
            return std::cos(alpha) * u
                   + std::sin(alpha) * CGAL::cross_product(axis, u)
                   + (1.0 - std::cos(alpha)) * (axis * u) * axis;
        };

        const double phi = (1.0 + std::sqrt(5.0)) * 0.5;
        const double twoPi = 2.0 * M_PI;


        for (size_t i = 0; i < parameters.nSamples; ++i) {
            double r = parameters.geodesicRadius * std::sqrt((i + 0.5) / static_cast<double>(parameters.nSamples));  // sample in unit disk so let's scale by geodesic radius that we want
            double theta = twoPi * i / phi;
            float r0 = 0.001;
            theta = std::fmod(theta, twoPi);
            if (theta < 0.0)
                theta += twoPi;

            for (int k = 0; k < planeAnglesSum.size() - 1; k++) {
                auto rep = theta / (2.0 * M_PI);

                if( rep >= planeAnglesSumNormalized[k] && rep < planeAnglesSumNormalized[(k+1)] )
                {

                    auto u = triangles[k].v0;
                    auto v = triangles[k].v1;

                    auto uhat = normalize(u);
                    auto vhat = normalize(v);
                    Vector n = CGAL::cross_product(uhat, vhat);
                    n = normalize(n);
                    float alpha = rep - planeAnglesSumNormalized[k];
                    alpha *= totalSum;

                    Vector w = rotateTowards(uhat, n, alpha);

                    auto vi = triangles[k].p0 + r0*w;
                    auto sample = Point(vi.x(), vi.y(), vi.z());
                    auto faceIndex = triangles[k].faceIndex;
                    SampleData sampleData;
                    sampleData.point = sample;
                    sampleData.faceIndex = faceIndex;
                    sampleData.r = r;
                    if(faceIndex == Mesh::null_face()) continue;
                    else
                    {
                        samples.push_back(sampleData);
                        break;
                    }
                }
            }

        }
        return samples;
    }
private:
    K::Point_2 computeUV(
            const Point& p,
            const Point& p0,
            const Vector& e0,
            const Vector& e1)
    {
        Vector d = p - p0;

        double a00 = CGAL::to_double(e0 * e0);
        double a01 = CGAL::to_double(e0 * e1);
        double a11 = CGAL::to_double(e1 * e1);

        double b0 = CGAL::to_double(d * e0);
        double b1 = CGAL::to_double(d * e1);

        double det = a00 * a11 - a01 * a01;

        double u = ( a11 * b0 - a01 * b1) / det;
        double v = (-a01 * b0 + a00 * b1) / det;

        return K::Point_2(u, v);
    }
    K::Point_2 computeUV(const Point& p, Mesh::Face_index faceIndex)
    {
        auto h = surfaceMesh.halfedge(faceIndex);

        auto v0 = source(h, surfaceMesh);
        auto v1 = target(h, surfaceMesh);
        auto v2 = target(next(h, surfaceMesh), surfaceMesh);

        Point p0 = surfaceMesh.point(v0);
        Point p1 = surfaceMesh.point(v1);
        Point p2 = surfaceMesh.point(v2);

        Vector e0 = p1 - p0;
        Vector e1 = p2 - p0;

        return computeUV(p, p0, e0, e1);
    }
    Point computeXYZ(const K::Point_2& uv, Mesh::Face_index face)
    {
        auto h = surfaceMesh.halfedge(face);

        const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
        const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));
        const Point& p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

        Vector e0 = p1 - p0;
        Vector e1 = p2 - p0;

        return p0 + uv.x() * e0 + uv.y() * e1;
    }
    std::optional<std::array<double, 3>> barycentricCoordinates(
            const Mesh& mesh,
            Mesh::Face_index face,
            const Point& p)
    {
        auto h = mesh.halfedge(face);

        const Point& p0 = mesh.point(source(h, mesh));
        const Point& p1 = mesh.point(target(h, mesh));
        const Point& p2 = mesh.point(target(next(h, mesh), mesh));

        Vector e0 = p1 - p0;
        Vector e1 = p2 - p0;
        Vector vp = p  - p0;

        double d00 = CGAL::to_double(e0 * e0);
        double d01 = CGAL::to_double(e0 * e1);
        double d11 = CGAL::to_double(e1 * e1);
        double d20 = CGAL::to_double(vp * e0);
        double d21 = CGAL::to_double(vp * e1);

        double denom = d00 * d11 - d01 * d01;

        constexpr double eps = 1e-12;

        if (std::abs(denom) < eps)
            return std::nullopt;

        double b1 = (d11 * d20 - d01 * d21) / denom;
        double b2 = (d00 * d21 - d01 * d20) / denom;
        double b0 = 1.0 - b1 - b2;

        if (b0 < -eps || b1 < -eps || b2 < -eps)
            return std::nullopt;

        return std::array<double, 3>{b0, b1, b2};
    }

    Point projectTriangleFace(Point p, Mesh::Face_index f )
    {

        // Get three vertices of the face
        auto h = surfaceMesh.halfedge(f);

        Point p0 = surfaceMesh.point(surfaceMesh.target(h));
        h = surfaceMesh.next(h);
        Point p1 = surfaceMesh.point(surfaceMesh.target(h));
        h = surfaceMesh.next(h);
        Point p2 = surfaceMesh.point(surfaceMesh.target(h));

        // Construct the face plane
        K::Plane_3 plane(p0, p1, p2);
        return plane.projection(p);
    }
    Point iterateMoveGeodesic(const SampleData& sample, const Vector& inputDirection, double geodesicDistance, Mesh::Face_index & outFace)
    {
        auto computeQ = [](const Vector& d,
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
        };
        auto computeIntersection = [&](K::Point_2 localP, K::Point_2 localP0, Mesh::Face_index &face) {

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
            auto copyray = rayDirect;
            copyray /= std::sqrt(CGAL::to_double(CGAL::squared_length(copyray)));
            for (const auto &he: CGAL::halfedges_around_face(surfaceMesh.halfedge(face), surfaceMesh)) {
                auto v1 = computeUV(surfaceMesh.point(source(he, surfaceMesh)), face);
                auto v2 = computeUV(surfaceMesh.point(target(he, surfaceMesh)), face);


                auto v3 = computeUV(surfaceMesh.point(target(surfaceMesh.next(he), surfaceMesh)), face);
                K::Vector_2 vec = v2-v1;

                if( abs(vec * rayDirect) > 0.99)
                {
                    K::Vector_2 inward(-vec.y(), vec.x()); // left normal
                    // Flip if necessary
                    if (inward * (v3 - v1) > 0)
                        inward = -inward;

                    inward /= std::sqrt(CGAL::to_double(CGAL::squared_length(inward)));

                    K::Vector_2 inward2(v3.x(), v3.y()); // left normal
                    rayDirect = rayDirect + .1 * inward2;
                    rayDirect /= std::sqrt(CGAL::to_double(CGAL::squared_length(rayDirect)));
                    foudD = true;
                    break;
                }

            }
            auto localP0_exterior = localP0 +  10.* rayDirect;
            auto direction = localP0 - localP0_exterior;
            direction /= std::sqrt(CGAL::to_double(CGAL::squared_length(direction)));
            K::Ray_2 ray(localP0_exterior, direction);

            K::Point_2 currentPoint;
            Mesh::Face_index nextFace = face;
            bool foundRes = false;
            // check halfedges
            float mindist = FLT_MAX;
            bool found = false;
            for (const auto &he: CGAL::halfedges_around_face(surfaceMesh.halfedge(face), surfaceMesh)) {
                auto v1 = computeUV(surfaceMesh.point(source(he, surfaceMesh)), face);
                auto v2 = computeUV(surfaceMesh.point(target(he, surfaceMesh)), face);
                auto v1v2 = v1 - v2;
                K::Segment_2 seg(v1 + 10e-5 * v1v2, v2 - 10e-5 * v1v2);
                segments.push_back(seg);
                auto result = CGAL::intersection(ray, seg);
                if (result) {
                    if (const K::Point_2 *p = std::get_if<K::Point_2>(&*result)) {
                        auto vec = *p - localP0_exterior;
                        if (vec.squared_length() < mindist) {

                            currentPoint = *p;
                            mindist = vec.squared_length();
                            nextFace = surfaceMesh.face(opposite(he, surfaceMesh));
                            found = true;
                        }
                    }
                    foundRes = true;
                }
            }

            if (found) {
                face = nextFace;
                return currentPoint;
            } else {

                std::cout << "ERROR \n";

                auto result1 = CGAL::intersection(ray, segments[0]);
                auto result2 = CGAL::intersection(ray, segments[1]);
                auto result3 = CGAL::intersection(ray, segments[2]);

                auto p = ray.source();
                auto dir = ray.direction();
                if(!foundRes)
                {
                    Point in = Point(localP0.x(), localP0.y(), 0 );
                    Point out = Point(localP0_exterior.x(), localP0_exterior.y(), 0 );

                    std::vector<Point> triangles;
                    auto h = surfaceMesh.halfedge(iFace);

                    auto v0 = source(h, surfaceMesh);
                    auto v1 = target(h, surfaceMesh);
                    auto v2 = target(next(h, surfaceMesh), surfaceMesh);

                    Point p0 = surfaceMesh.point(v0);
                    Point p1 = surfaceMesh.point(v1);
                    Point p2 = surfaceMesh.point(v2);

                    auto lp0 = computeUV(p0, iFace);
                    auto lp1 = computeUV(p1, iFace);
                    auto lp2 = computeUV(p2, iFace);

                    triangles.push_back(Point(lp0.x(), lp0.y(), 0.));
                    triangles.push_back(Point(lp1.x(), lp1.y(), 0.));
                    triangles.push_back(Point(lp2.x(), lp2.y(), 0.));

                    std::vector<std::array<size_t, 2>> lines;
                    lines.push_back({0,1});
                    lines.push_back({1,2});
                    lines.push_back({2,0});

                    auto a = polyscope::registerCurveNetwork("triangles", triangles, lines);
                    a->setRadius(0.0001);
                    a->resetTransform();

                    debug.inputs.push_back(in);
                    debug.outputs.push_back(out);
                    std::vector<Point> line = {in, out};
                    std::vector<std::array<size_t, 2>> edges  = {
                            {0, 1}
                    };
                    auto b = polyscope::registerCurveNetwork("line", line,edges);
                    b->setRadius(0.0001);
                    b->resetTransform();

                    {
                        auto p2 = localP0_exterior + copyray;
                        Point in2 = Point(localP0.x(), localP0.y(), 0 );
                        Point out2 = Point(p2.x(), p2.y(), 0 );

                        std::vector<Point> line2 = {in2, out2};
                        std::vector<std::array<size_t, 2>> edges2  = {
                                {0, 1}
                        };
                        auto b = polyscope::registerCurveNetwork("line2", line2,edges2);
                        b->setRadius(0.0001);
                        b->resetTransform();
                    }
                }
                else
                {
                    std::cout<<"Found line segment \n";

                }
                /*if(bugFound) return currentPoint;
                bugFound = true;*/

            }

        };

        bool isInside = false;

        Mesh::Face_index currentFace = sample.faceIndex;
        Mesh::Face_index nextFace = currentFace;
        Vector direction = inputDirection;
        Point input = sample.point;
        Point output = sample.point + direction;
        double remainingDistance = geodesicDistance;

        std::vector<Point> sequences;
        int iteration = 0;
        bool cont = true;
        while(cont && !bugFound)
        {

            /*auto bary2 = barycentricCoordinates(surfaceMesh, currentFace, input);
            K::Point_2 input2D1, input2D2;
            if(!bary2)
            {
                auto h = surfaceMesh.halfedge(currentFace);
                input2D1 = computeUV(input, currentFace);
                currentFace = surfaceMesh.face(surfaceMesh.opposite(h));
                input2D2 = computeUV(input, currentFace);
            }*/

            sequences.push_back(input);
            auto input2D = computeUV(input, currentFace);
            auto output2D = computeUV(output, currentFace);

            auto normalCurrentFace = CGAL::Polygon_mesh_processing::compute_face_normal(currentFace,surfaceMesh);

            auto intersection2D = computeIntersection(output2D, input2D, nextFace);

            double x = intersection2D.x();
            double y = intersection2D.y();
            float eps = 10e-5;
            if( x < eps) x = eps;
            if( y < eps) y = eps;

            if( x > 1-eps) x = 1-eps;
            if( y > 1-eps) y = 1-eps;
            intersection2D = K::Point_2 (x,y);

            Point intersection3D = computeXYZ(intersection2D, currentFace);

            if(nextFace == Mesh::null_face())
            {
                outFace = currentFace;
                return intersection3D;
            }
            double distanceToIntersection = std::sqrt(CGAL::squared_length(intersection3D - input));

            auto normalNextFace = CGAL::Polygon_mesh_processing::compute_face_normal(nextFace,surfaceMesh);
            Vector newDirection;
            if(CGAL::to_double(normalCurrentFace * normalNextFace)  > 0.99)
                newDirection = direction;
            else
            {
                newDirection = computeQ(direction, normalCurrentFace, normalNextFace);
            }
            newDirection /= std::sqrt(CGAL::squared_length(newDirection));

            input = intersection3D;
            output = input + newDirection;
            remainingDistance-= distanceToIntersection;
            Point target = input + remainingDistance* newDirection;
            auto bary = barycentricCoordinates(surfaceMesh, nextFace, target);

            currentFace = nextFace;
            direction = newDirection;
            outFace = currentFace;

            if(bugFound)
            {
                auto seq = polyscope::registerPointCloud("quences", sequences);
                seq->resetTransform();
                return target;
            }
            if(bary)
            {
                return target;
            }
            iteration++;
            if(iteration > 100 || remainingDistance < 0)
            {
                cont = false;
                return intersection3D;
            }
        }

        return input;

    }
    Point movePointGeodesic( SampleData& sample, const Vector& direction, double geodesicDistance, Mesh::Face_index & outFace)
    {

        auto directionNormalized =  direction / std::sqrt(CGAL::squared_length(direction));
        Point out = sample.point + directionNormalized;
        Point outSurface = projectTriangleFace(out, sample.faceIndex);
        auto faceIndex = sample.faceIndex;
        Vector directionSurface = outSurface-sample.point;
        directionSurface /= std::sqrt(CGAL::squared_length(directionSurface));

        auto pointMoved = sample.point + geodesicDistance*directionSurface;
        auto bary = barycentricCoordinates(surfaceMesh, faceIndex, pointMoved);
        outFace = faceIndex;
        if (bary)
        {
            return pointMoved;
        } else
        {
            auto bary2 = barycentricCoordinates(surfaceMesh, faceIndex,  sample.point);
            if(!bary2)
            {
                auto h = surfaceMesh.halfedge(faceIndex);
                sample.faceIndex = surfaceMesh.face(surfaceMesh.opposite(h));
            }

            return iterateMoveGeodesic(sample, directionSurface, geodesicDistance, outFace);
            //return pointMoved;
        }
    }

    bool bugFound = false;
    DEBUG debug;
    std::vector<UserType> vertexMetrics;
    Mesh surfaceMesh;

public:
    std::vector<std::vector<Point>> samplingPerVertex;
    std::vector<std::vector<Point>> samplingInitialPerVertex;
};
