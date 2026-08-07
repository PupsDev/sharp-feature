#pragma once
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
K::Point_2 computeUV(const Mesh& surfaceMesh, const Point& p, Mesh::Face_index faceIndex)
{
    assert(faceIndex != Mesh::null_face());

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

Eigen::Vector2d coordinatesInBasis(
        const Vector& v,
        const Vector& e1,
        const Vector& e2)
{
    Eigen::Matrix<double,3,2> A;
    A.col(0) << e1.x(), e1.y(), e1.z();
    A.col(1) << e2.x(), e2.y(), e2.z();

    auto normal = cross_product(e1,e2);
    normal /= std::sqrt(CGAL::squared_length(normal));

    Eigen::Vector3d b(v.x(), v.y(), v.z());

    // Since e1,e2 span the plane, solve A*x = b
    return A.colPivHouseholderQr().solve(b);
}
Eigen::Vector2d coordinatesInBasis2D(
        const K::Vector_2& v,
        const Vector& e1,
        const Vector& e2)
{
    Eigen::Matrix<double,3,2> A;
    A.col(0) << e1.x(), e1.y(), e1.z();
    A.col(1) << e2.x(), e2.y(), e2.z();

    auto normal = cross_product(e1,e2);
    normal /= std::sqrt(CGAL::squared_length(normal));

    Eigen::Vector3d b(v.x(), v.y(), 0.);

    // Since e1,e2 span the plane, solve A*x = b
    return A.colPivHouseholderQr().solve(b);
}
K::Vector_2 deltaXYZToUV(const Mesh& surfaceMesh, const K::Vector_2 & deltaXYZ, const Point& p,Mesh::Face_index face)
{
    auto h = surfaceMesh.halfedge(face);

    const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
    const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));
    const Point& p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

    Vector e0 = p1 - p0;
    e0 /= sqrt(CGAL::squared_length(e0));


    Vector t;
    if (std::abs(e0.x()) < 0.9)
        t = Vector(1, 0, 0);
    else
        t = Vector(0, 1, 0);

    // First tangent
    Vector u = CGAL::cross_product(e0, t);
    // Second tangent
    Vector v = CGAL::cross_product(e0, u);
    float dx = deltaXYZ.x();
    float dy = deltaXYZ.y();
    Point d = p + dx * e0 + dy * v;

    return computeUV(surfaceMesh, d, face) - CGAL::ORIGIN;

}
K::Vector_2 deltaXYZToDeltaUV(const Mesh& surfaceMesh, const Vector & deltaXYZ, Mesh::Face_index face)
{
    auto h = surfaceMesh.halfedge(face);

    const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
    const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));
    const Point& p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

    Vector e0 = p1 - p0;
    e0 /= sqrt(CGAL::squared_length(e0));

    Vector t;
    if (std::abs(e0.x()) < 0.9)
        t = Vector(1, 0, 0);
    else
        t = Vector(0, 1, 0);

    // First tangent
    Vector u = CGAL::cross_product(e0, t);
    u /= sqrt(CGAL::squared_length(u));

    auto vec = coordinatesInBasis(deltaXYZ, e0, u);
    return {vec.x(), vec.y()};
}
Point deltaXYZToP(const Mesh& surfaceMesh, const K::Vector_2 & deltaXYZ,Mesh::Face_index face)
{
    auto h = surfaceMesh.halfedge(face);
    const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
    const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));

    Vector e0 = p1 - p0;
    e0 /= sqrt(CGAL::squared_length(e0));

    Vector t;
    if (std::abs(e0.x()) < 0.9)
        t = Vector(1, 0, 0);
    else
        t = Vector(0, 1, 0);

    // First tangent
    Vector u = CGAL::cross_product(e0, t);
    u /= sqrt(CGAL::squared_length(u));


     return CGAL::ORIGIN+ deltaXYZ.x() *e0 +  deltaXYZ.y() * u;
}
K::Vector_2 deltaXYZToUVOld(const Mesh& surfaceMesh, const K::Vector_2 & deltaXYZ, const Point& p,Mesh::Face_index face)
{

    bool useDefaultBase = false;
    if(useDefaultBase)
    {
        auto h = surfaceMesh.halfedge(face);

        const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
        const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));
        const Point& p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

        Vector e0 = p1 - p0;
        Vector e1 = p2 - p0;

        auto l0 = sqrt(e0.squared_length());
        auto l1 = sqrt(e1.squared_length());
        return {deltaXYZ.x() / l0,  deltaXYZ.y() / l1};
    }
    else
    {
        auto h = surfaceMesh.halfedge(face);

        const Point& p0 = surfaceMesh.point(source(h, surfaceMesh));
        const Point& p1 = surfaceMesh.point(target(h, surfaceMesh));
        const Point& p2 = surfaceMesh.point(target(next(h, surfaceMesh), surfaceMesh));

        Vector e0 = p1 - p0;
        e0 /= sqrt(CGAL::squared_length(e0));

        Vector t;
        if (std::abs(e0.x()) < 0.9)
            t = Vector(1, 0, 0);
        else
            t = Vector(0, 1, 0);

        // First tangent
        Vector u = CGAL::cross_product(e0, t);
        u /= sqrt(CGAL::squared_length(u));


        auto pos = p0 + deltaXYZ.x() *e0 +  deltaXYZ.y() * u;
        return computeUV(surfaceMesh,pos, face)-computeUV(surfaceMesh,p0, face);
    }
}
Point computeXYZ(const Mesh& surfaceMesh,const K::Point_2& uv, Mesh::Face_index face)
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