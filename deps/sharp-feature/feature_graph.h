#pragma once
#include "type.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/heap/priority_queue.hpp>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

namespace sharp_feature
{
    template<typename UserType>
    class FeatureGraph
    {
    public:
        typedef GraphElement<UserType> VertexProperty;

        typedef boost::adjacency_list<
                boost::vecS,
                boost::listS,
                boost::undirectedS,
                VertexProperty
        >
                Graph;

        typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
        typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;

        struct GraphPriority
        {
            Vertex vertex;
            UserType priority;
            bool operator<(const GraphPriority& other) const {
                return priority > other.priority; // reverse for min-heap
            }
        };
        FeatureGraph() = default;

        // Select Surfels given an input metric and a threshold
        void computeSelectedSurfels(const Mesh& sm, const std::vector<UserType>& vertexMetrics, const FeatureGraphParameters<UserType>& parameters)
        {
            selectedSurfels.resize(sm.number_of_faces(), 0);
            for (const auto& face : sm.faces())
            {
                for (const auto &vertexHandle : CGAL::vertices_around_face(sm.halfedge(face), sm))
                {
                    if (vertexMetrics[vertexHandle.idx()] > parameters.threshold)
                    {
                        selectedSurfels[face.idx()] = 1;
                    }
                }
            }

        }
        std::unordered_map<CGAL::SM_Halfedge_index, int> computeBorderHalfedges(const Mesh& sm)
        {
            std::unordered_map<CGAL::SM_Halfedge_index, int> halfEdgeBorder;

            // TODO how does Surface Mesh DS works for non manifold ?
            for (const auto& halfedge : sm.halfedges())
            {
                // Get face indices, if they exist
                auto f = sm.face(halfedge);
                auto fo = sm.face(sm.opposite(halfedge));

                if (f == Mesh::null_face()) continue; // skip if current face is null

                if (fo == Mesh::null_face()) {
                    // Border edge: opposite has no face
                    halfEdgeBorder[halfedge] = 1;
                    continue;
                }
                // Compare selections
                auto fIdx = f.idx();
                auto foIdx = fo.idx();

                halfEdgeBorder[halfedge] = (selectedSurfels[fIdx] != selectedSurfels[foIdx]) ? 1 : 0;
            }
            return halfEdgeBorder;
        }
        void computeFeatureGraphNew(const Mesh& sm, const std::vector<UserType>& vertexMetrics, const std::unordered_map<Mesh::Edge_index , double>& subdivisionMetric)
        {
            std::vector<Vertex> pointels(sm.number_of_vertices());
            std::vector<Vertex> linels(sm.number_of_edges());
            std::vector<Vertex> surfels(sm.number_of_faces());

            std::vector<Vector > faceNormals(sm.number_of_faces());
            std::unordered_map<Mesh::Edge_index, Vertex> halfEdgeMap;
            for(const auto& face : sm.faces())
            {
                faceNormals[face.idx()] = CGAL::Polygon_mesh_processing::compute_face_normal(face, sm);
            }
            size_t currentFace = 0;
            size_t vertexId = 0;
            for (const auto& face : sm.faces()) {
                std::cout << "Graph construction "
                          << static_cast<double >(currentFace) * 100. / static_cast<double>(sm.num_faces()) << "\n";
                currentFace++;
                ;
                if (selectedSurfels[face.id()]) {
                    // Compute Surfel point as centroid
                    const auto &centroid = computeCentroid(sm, sm.halfedge(face));
                    UserType centroidMetric = UserType(0.);
                    size_t total      = 0;
                    for (const auto& vertexHandle : sm.vertices_around_face(sm.halfedge(face)))
                    {
                        centroidMetric += vertexMetrics[vertexHandle.idx()];
                        total++;
                    }
                    centroidMetric /= static_cast<UserType>(total);

                    Vertex surfelElement = boost::add_vertex(graph);
                    auto element = GraphElement<UserType>();
                    element.id = vertexId++;
                    element.position =centroid;
                    element.weight = 0.;//centroidMetric;
                    element.border = false;
                    element.type = GraphElementType::SURFEL;
                    graph[surfelElement] = element;
                    surfels[face.id()] = surfelElement;

                    for(auto vertex : sm.vertices_around_face(sm.halfedge(face)))
                    {
                        Vertex pointelElement = boost::add_vertex(graph);
                        auto elementPointel = GraphElement<UserType>();
                        elementPointel.id = vertexId;
                        elementPointel.position =  sm.point(vertex);
                        elementPointel.weight =  vertexMetrics[vertex.idx()];
                        elementPointel.border = false;
                        elementPointel.type = GraphElementType::POINTEL;
                        graph[pointelElement] = elementPointel;
                        pointels[vertex.idx()] = pointelElement;

                    }
                    for(auto halfedge : sm.halfedges_around_face(sm.halfedge(face)))
                    {
                        const auto& sourceHandle = sm.source(halfedge);
                        const auto& targetHandle = sm.target(halfedge);
                        auto edge = sm.edge(halfedge);
                        if(halfEdgeMap.contains(edge))
                        {
                            linels[sm.edge(halfedge).idx()] = halfEdgeMap[edge];
                        }
                        if(!sm.is_border(halfedge))
                        {
                            const auto &point = computeMidPointHalfEdge(sm, halfedge);
                            // Add midpoint Linel
                            Vertex linelElement = boost::add_vertex(graph);
                            auto elementLinel = GraphElement<UserType>();
                            elementLinel.id = vertexId;
                            elementLinel.position = point;
                            if(subdivisionMetric.contains(sm.edge(halfedge)))
                            {
                                elementLinel.weight =  subdivisionMetric.at(sm.edge(halfedge));
                            }
                            else
                            {
                                elementLinel.weight =  static_cast<UserType>(0.5)*(vertexMetrics[sourceHandle.idx()]+vertexMetrics[targetHandle.idx()]);
                                if(sm.face(sm.opposite(halfedge)) != Mesh::null_face())
                                {
                                    const auto &faceNormal = faceNormals[sm.face(halfedge).idx()];
                                    const auto &faceNormalOpposite = faceNormals[sm.face(sm.opposite(halfedge)).idx()];
                                    const auto &angle =
                                            static_cast<UserType>(1.) - abs(CGAL::scalar_product(faceNormal, faceNormalOpposite));
                                    if(angle < 10e-3) elementLinel.weight = 0.;
                                }

                            }
                            elementLinel.border = sm.is_border(halfedge);
                            elementLinel.type = GraphElementType::LINEL;
                            linels[edge.idx()] = linelElement;
                            halfEdgeMap[edge] = linelElement;
                            graph[linelElement] = elementLinel;
                        }
                    }
                }

            }
            for(auto f : sm.faces())
            {
                if(!selectedSurfels[f.idx()])
                    continue;

                auto surfel = surfels[f.idx()];
                for(auto h : halfedges_around_face(sm.halfedge(f), sm))
                {
                    boost::add_edge(surfel, pointels[sm.source(h).idx()], graph);
                    boost::add_edge(surfel, linels[sm.edge(h).idx()], graph);

                    boost::add_edge(linels[sm.edge(h).idx()], pointels[sm.source(h).idx()], graph);
                    boost::add_edge(linels[sm.edge(h).idx()], pointels[sm.target(h).idx()], graph);
                }
            }
        }
        void computeFeatureGraph(const Mesh& sm, const std::vector<UserType>& vertexMetrics, const std::unordered_map<Mesh::Edge_index , double>& subdivisionMetric)
        {
            std::unordered_map<size_t, size_t> vertexMap;
            std::unordered_map<size_t, size_t> halfEdgeMap;
            std::unordered_map<CGAL::SM_Halfedge_index, int> halfEdgeBorder;

            // TODO how does Surface Mesh DS works for non manifold ?
            for (const auto& halfedge : sm.halfedges())
            {
                // Get face indices, if they exist
                auto f = sm.face(halfedge);
                auto fo = sm.face(sm.opposite(halfedge));

                if (f == Mesh::null_face()) continue; // skip if current face is null

                if (fo == Mesh::null_face()) {
                    // Border edge: opposite has no face
                    halfEdgeBorder[halfedge] = 1;
                    continue;
                }
                // Compare selections
                auto fIdx = f.idx();
                auto foIdx = fo.idx();

                halfEdgeBorder[halfedge] = (selectedSurfels[fIdx] != selectedSurfels[foIdx]) ? 1 : 0;
            }

            size_t currentFace = 0;
            size_t vertexId = 0;
            for (const auto& face : sm.faces())
            {
                std::cout<<"Graph construction "<<static_cast<double >(currentFace)*100. / static_cast<double>(sm.num_faces())<<"\n";
                currentFace++;
                if (selectedSurfels[face.id()]) {

                    // Compute Surfel point as centroid
                    const auto &centroid = computeCentroid(sm, sm.halfedge(face));
                    UserType centroidMetric = UserType(0.);
                    size_t total      = 0;
                    for (const auto& vertexHandle : sm.vertices_around_face(sm.halfedge(face)))
                    {
                        centroidMetric += vertexMetrics[vertexHandle.idx()];
                        total++;
                    }
                    centroidMetric /= static_cast<UserType>(total);

                    Vertex surfelElement = boost::add_vertex(graph);
                    auto element = GraphElement<UserType>();
                    element.id = vertexId++;
                    element.position =centroid;
                    //element.weight = 0;//centroidMetric;
                    element.weight = 0;//centroidMetric;
                    element.border = false;
                    element.type = GraphElementType::SURFEL;
                    graph[surfelElement] = element;

                    std::unordered_set<CGAL::SM_Vertex_index> visitedPointels;
                    // Iterate over face halfedges to compute remaining Pointels
                    for (const auto &halfedge: CGAL::halfedges_around_face(sm.halfedge(face), sm))
                    {
                        const auto& sourceHandle = sm.source(halfedge);
                        const auto& targetHandle = sm.target(halfedge);

                        // Add unique vertex for each pointel
                        if(!vertexMap.count(sourceHandle.idx()))
                        {
                            Vertex pointelElement = boost::add_vertex(graph);
                            auto elementPointel = GraphElement<UserType>();
                            elementPointel.id = vertexId;
                            elementPointel.position =  sm.point(sourceHandle);
                            elementPointel.weight =  vertexMetrics[sourceHandle.idx()];
                            elementPointel.border = false;
                            elementPointel.type = GraphElementType::POINTEL;

                            graph[pointelElement] = elementPointel;
                            vertexMap[sourceHandle.idx()] = vertexId++;

                        }
                        if(!vertexMap.count(targetHandle.idx()))
                        {
                            Vertex pointelElement = boost::add_vertex(graph);
                            auto elementPointel = GraphElement<UserType>();
                            elementPointel.id = vertexId;
                            elementPointel.position =  sm.point(targetHandle);
                            elementPointel.weight =  vertexMetrics[targetHandle.idx()];
                            elementPointel.border = false;
                            elementPointel.type = GraphElementType::POINTEL;

                            graph[pointelElement] = elementPointel;
                            vertexMap[targetHandle.idx()] = vertexId++;

                        }

                        // Add unique edges between Surfels and Pointels
                        if(!visitedPointels.count(sourceHandle))
                        {
                            boost::add_edge(surfelElement, boost::vertex(vertexMap[sourceHandle.idx()], graph),  graph);
                            visitedPointels.insert(sourceHandle);
                        }
                        if(!visitedPointels.count(targetHandle))
                        {
                            boost::add_edge(surfelElement, boost::vertex(vertexMap[targetHandle.idx()], graph), graph);
                            visitedPointels.insert(targetHandle);
                        }

                        if(!halfEdgeMap.count(sm.opposite(halfedge)))
                        {

                            const auto &point = computeMidPointHalfEdge(sm, halfedge);
                            // Add midpoint Linel
                            Vertex linelElement = boost::add_vertex(graph);

                            auto elementLinel = GraphElement<UserType>();
                            elementLinel.id = vertexId;
                            elementLinel.position = point;
                            //elementLinel.weight =  static_cast<UserType>(0.5)*(vertexMetrics[sourceHandle.idx()]+vertexMetrics[targetHandle.idx()]);
                            if(subdivisionMetric.contains(sm.edge(halfedge)))
                            {
                                elementLinel.weight =  subdivisionMetric.at(sm.edge(halfedge));
                            }
                            else
                            {
                                elementLinel.weight =  static_cast<UserType>(0.5)*(vertexMetrics[sourceHandle.idx()]+vertexMetrics[targetHandle.idx()]);
                                auto h1 = halfedge;
                                const auto &faceNormal = CGAL::Polygon_mesh_processing::compute_face_normal(sm.face(h1), sm);
                                if(sm.face(sm.opposite(h1)) == Mesh::null_face())
                                {
                                    elementLinel.weight = 0.;
                                }
                                else
                                {
                                    const auto &faceNormalOpposite = CGAL::Polygon_mesh_processing::compute_face_normal(
                                            sm.face(sm.opposite(h1)), sm);
                                    const auto &angle =
                                            static_cast<UserType>(1.) - abs(CGAL::scalar_product(faceNormal, faceNormalOpposite));
                                    if(angle < 10e-3) elementLinel.weight = 0.;
                                }

                            }
                            elementLinel.border =  static_cast<bool>(halfEdgeBorder[halfedge]);
                            elementLinel.type = GraphElementType::LINEL;
                            graph[linelElement] = elementLinel;

                            halfEdgeMap[halfedge] = vertexId;
                            vertexId++;

                            // Add edges for all graph element
                            // Pointel edge
                            boost::add_edge(linelElement, boost::vertex(vertexMap[sm.source(halfedge).idx()], graph), graph);
                            boost::add_edge(linelElement, boost::vertex(vertexMap[sm.target(halfedge).idx()], graph), graph);

                            if(halfEdgeBorder[halfedge])
                            {
                                graph[boost::vertex(vertexMap[sm.source(halfedge).idx()], graph)].border = true;
                                graph[boost::vertex(vertexMap[sm.target(halfedge).idx()], graph)].border = true;
                            }

                            // Surfel edge
                            boost::add_edge(linelElement, surfelElement, graph);



                            /*
                             * // only for quad maybe ? to test
                            boost::add_edge( boost::vertex(vertexMap[sm.source(halfedge).idx()], graph), boost::vertex(vertexMap[sm.target(halfedge).idx()], graph), graph);
                            if(halfEdgeBorder[halfedge])
                            {
                                graph[boost::vertex(vertexMap[sm.source(halfedge).idx()], graph)].border = true;
                                graph[boost::vertex(vertexMap[sm.target(halfedge).idx()], graph)].border = true;
                            }*/
                        }
                        else
                        {
                            auto vertexIndex = halfEdgeMap[sm.opposite(halfedge)];
                            boost::add_edge(boost::vertex(vertexIndex, graph), surfelElement, graph);
                        }

                    }
                }

            }

        }
        void computeThinning(const Mesh& sm)
        {
            boost::heap::priority_queue<GraphPriority> pq;
            for(const auto& vertice : boost::make_iterator_range(boost::vertices(graph)))
            {
                GraphPriority element;
                element.priority = graph[vertice].weight;
                element.vertex= vertice;
                pq.push(element);
            }
            std::set<Vertex> nonRemovable;
            std::set<Vertex> deleted;

            size_t id = 0;
            while(!pq.empty() )
            {
                auto v = pq.top();
                pq.pop();

                if(!isRemovable(v.vertex))
                {
                    nonRemovable.insert(v.vertex);
                    continue;
                }
                if(deleted.count(v.vertex)) continue;
                for(const auto& neighbor : boost::make_iterator_range(boost::adjacent_vertices(v.vertex, graph)))
                {
                    graph[neighbor].border = true;
                    if(nonRemovable.count(neighbor))
                    {
                        GraphPriority element;
                        element.priority = graph[neighbor].weight;
                        element.vertex= neighbor;
                        pq.push(element);
                    }
                }
                deleted.insert(v.vertex);
                boost::clear_vertex(v.vertex, graph);
                boost::remove_vertex(v.vertex, graph);

                id++;
            }
        }



        // output
        void exportGraphSimple(std::vector<Point> &points, std::vector<std::vector<int>>& graphSimple)
        {
            vertexIndexMap.clear();
            points.clear();
            size_t index = 0;
            for(const auto& vertice : boost::make_iterator_range(boost::vertices(graph)))
            {
                vertexIndexMap[vertice] = index++;
                points.push_back(graph[vertice].position);
            }
            graphSimple.resize(vertexIndexMap.size());

            for(const auto& edge : boost::make_iterator_range(boost::edges(graph)))
            {
                const auto &source = boost::source(edge, graph);
                const auto &target = boost::target(edge, graph);

                graphSimple[vertexIndexMap[source]].push_back(vertexIndexMap[target]);
                graphSimple[vertexIndexMap[target]].push_back(vertexIndexMap[source]);
            }
        }
        [[nodiscard]] void exportGraph(std::vector<Point> &points, std::vector<std::array<size_t, 2>>& edges, std::vector<UserType>& nodeMetric)
        {
            vertexIndexMap.clear();
            size_t index = 0;
            for(const auto& vertice : boost::make_iterator_range(boost::vertices(graph)))
            {
                vertexIndexMap[vertice] = index++;
                points.push_back(graph[vertice].position);
                nodeMetric.push_back(graph[vertice].weight);
            }
            for(const auto& edge : boost::make_iterator_range(boost::edges(graph)))
            {
                const auto &source = boost::source(edge, graph);
                const auto &target = boost::target(edge, graph);

                std::array<size_t, 2> tab = {vertexIndexMap[source], vertexIndexMap[target]};
                edges.push_back(tab);

            }
        }
        [[nodiscard]] std::vector<int> getSelectedFaces() const
        {
            return selectedSurfels;
        }
    private:
        [[nodiscard]] bool isSimple(const Vertex& vertex)
        {
            boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>  subgraph;
            std::set<Vertex> visited;

            std::unordered_map<Vertex, unsigned __int64> mapNeighbors;

            if(boost::out_degree(vertex, graph) < 1) return true;

            for(const auto& neighbor : boost::make_iterator_range(boost::adjacent_vertices(vertex, graph)))
            {
                if(!visited.count(neighbor))
                {
                    auto element = boost::add_vertex(subgraph);
                    visited.insert(neighbor);
                    mapNeighbors[neighbor] = element;
                }

                for(const auto& neighbor2 : boost::make_iterator_range(boost::adjacent_vertices(neighbor, graph)))
                {
                    if(visited.count(neighbor2) && neighbor2 != vertex)
                    {
                        boost::add_edge(mapNeighbors[neighbor], mapNeighbors[neighbor2], subgraph);
                    }
                }
            }

            if(visited.empty()) return true;
            std::vector< int > component(num_vertices(subgraph));
            auto num = boost::connected_components(subgraph, &component[0]);

            return num == 1;
        }


        [[nodiscard]] bool isRemovable(const Vertex& vertex)
        {
            if(graph[vertex].weight < 0.1) return true;
            if(!graph[vertex].border) return false;
            bool condition = true;
            if(graph[vertex].type == POINTEL)
            {
                condition = boost::out_degree(vertex, graph) > 1;
            }
            return isSimple(vertex) && condition;
        }
        [[nodiscard]] Point computeMidPointHalfEdge(const Mesh& sm, const CGAL::SM_Halfedge_index& halfEdge)
        {
            const auto& v1                  = sm.point(sm.source(halfEdge));
            const auto& v2                  = sm.point(sm.target(halfEdge));
            const auto& mid                 = CGAL::midpoint(v1, v2);
            return mid;
        }
        [[nodiscard]] static Point computeCentroid(const Mesh& sm, const CGAL::SM_Halfedge_index& face)
        {
            Vector centroid = CGAL::NULL_VECTOR;
            size_t total      = 0;
            for (const auto& vertexHandle : sm.vertices_around_face(face))
            {
                centroid += (sm.point(vertexHandle) - CGAL::ORIGIN);
                total++;
            }
            return CGAL::ORIGIN + (centroid / static_cast<double>(total));
        }

        std::vector<int> selectedSurfels;
        std::unordered_map<Vertex, size_t> vertexIndexMap;
        Graph graph;
    };
};