#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
using FT = typename K::FT;
using Point = typename K::Point_3;
using Vector = typename K::Vector_3;
typedef CGAL::Surface_mesh<Point> Mesh;

using Plane = std::pair<Point, Vector>;
static const char* MethodMetricNames[] = {
        "DIHEDRAL",
        "NVT",
        "NVT_Geodesic",
        "DISCRETE_MEAN",
        "DISCRETE_GAUSSIAN",
        "CURVEDNESS"
};
enum MethodMetricType
{
    DIHEDRAL,
    NVT,
    NVT_Geodesic,
    DISCRETE_MEAN,
    DISCRETE_GAUSSIAN,
    CURVEDNESS,
    METHOD_COUNT
};
template<typename UserType>
struct FeatureMetricParameters
{
    UserType threshold = UserType(0.5);
    int maxDepth = 5;
    UserType maxRadius = UserType(2.);
    MethodMetricType currentMethod = DIHEDRAL;

    bool useAngleWeight = false;
    bool useAreaWeight = false;
    bool useNormalize = false;

    bool useGeodesicWeight = false;
    float sigma = 0.1;
    int nSamples = 32;
    float geodesicRadius = 0.1;
};

template<typename UserType>
struct FeatureGraphParameters
{
    UserType threshold= UserType(0.1);
    void setThreshold(const float& selectionThreshold)
    {
        threshold = static_cast<UserType>(selectionThreshold);
    }
};

struct PolyLineAttribute
{
    size_t degree;
    size_t vertexID;
};
template <typename Point>
struct Polylines {
    std::vector <Point> points;
    std::vector <std::vector<std::size_t>> lines;

    std::vector<PolyLineAttribute> attributes;
};

enum GraphElementType
{
    POINTEL,
    LINEL,
    SURFEL
};

template<typename UserType>
struct GraphElement
{
    size_t id;
    Point position;
    UserType weight;
    bool border = false;
    GraphElementType type;

    GraphElement() = default;
    GraphElement (GraphElement<UserType>& graphElement)
    {
        id = graphElement.id;
        position = graphElement.position;
        weight = graphElement.weight;
        border = graphElement.border;
        type = graphElement.type;
    }
    GraphElement (const GraphElement<UserType>& graphElement)
    {
        id = graphElement.id;
        position = graphElement.position;
        weight = graphElement.weight;
        border = graphElement.border;
        type = graphElement.type;
    }

};

template <typename Point>
void export_surface_mesh_to_vectors(
        const CGAL::Surface_mesh<Point>& mesh,
        std::vector<glm::vec3>& vertexPositions_,
        std::vector<std::vector<size_t>>& facesIn)
{
    vertexPositions_.clear();
    facesIn.clear();

    using Vertex_index = typename Mesh::Vertex_index;
    using Face_index = typename Mesh::Face_index;

    std::unordered_map<Vertex_index, size_t> vertexIndexMap;
    vertexIndexMap.reserve(mesh.number_of_vertices());

    // Assign indices and collect vertex positions
    size_t idx = 0;
    for (Vertex_index v : mesh.vertices()) {
        const Point& p = mesh.point(v);
        vertexPositions_.emplace_back(static_cast<float>(p.x()),
                                      static_cast<float>(p.y()),
                                      static_cast<float>(p.z()));
        vertexIndexMap[v] = idx++;
    }

    // Collect face indices
    for (Face_index f : mesh.faces()) {
        std::vector<size_t> face;
        for (Vertex_index v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            face.push_back(vertexIndexMap[v]);
        }
        if (face.size() >= 3) { // ignore degenerate faces
            facesIn.push_back(face);
        }
    }
}