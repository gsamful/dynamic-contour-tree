#ifndef __TERRAIN_H__
#define __TERRAIN_H__

#include <vector>
#include <algorithm>
#include <set>

#include "vertex.h"

template <class VT>
class StarGridTerrain {
	std::vector<VT> points;

	public:
	typedef VT Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef VertexPair<const Vertex> Edge;
	typedef VertexRefLess<Vertex> VertexHeightComp;
	typedef typename Vertex::HeightType HeightType;

	// height matrix constructor
	// +* adjacency is used
	StarGridTerrain(const std::vector< std::vector<HeightType> > & mat) {
		for (size_t i = 0 ; i < mat.size() ; ++i) {
			for (size_t j = 0 ; j < mat[0].size() ; ++j) {
				points.push_back(Vertex(mat[i][j]));
			}
		}
		for (size_t i = 0 ; i < mat.size() ; ++i) {
			for (size_t j = 0 ; j < mat[0].size() ; ++j) {
				if ((i + j) % 2 == 0) {
					if (j < mat.size()-1) {
						if (i > 0) {
							points[i*mat[0].size()+j].ConnectBack(points[(i-1)*mat[0].size()+j+1]);
							points[(i-1)*mat[0].size()+j+1].ConnectBack(points[i*mat[0].size()+j]);
//							std::cout << points[i*mat[0].size()+j] << " - " << points[(i-1)*mat[0].size()+j+1] << std::endl;
						}
					} 	
				}
				if (j < mat[0].size()-1) {
					points[i*mat[0].size()+j].ConnectBack(points[i*mat[0].size()+j+1]);
					points[i*mat[0].size()+j+1].ConnectFront(points[i*mat[0].size()+j]);
//					std::cout << points[i*mat[0].size()+j] << " - " << points[i*mat[0].size()+j+1] << std::endl;
				}
				if ((i + j) % 2 == 0) {
					if (j < mat.size()-1) {
						if (i < mat.size()-1) {
							points[i*mat[0].size()+j].ConnectBack(points[(i+1)*mat[0].size()+j+1]);
							points[(i+1)*mat[0].size()+j+1].ConnectBack(points[i*mat[0].size()+j]);
//							std::cout << points[i*mat[0].size()+j] << " - " << points[(i+1)*mat[0].size()+j+1] << std::endl;
						} 
					}
				}
				if (i < mat.size()-1) {
					points[i*mat[0].size()+j].ConnectBack(points[(i+1)*mat[0].size()+j]);
					points[(i+1)*mat[0].size()+j].ConnectBack(points[i*mat[0].size()+j]);
//					std::cout << points[i*mat[0].size()+j] << " - " << points[(i+1)*mat[0].size()+j] << std::endl;
				}
			}
		}
//		std::cout << "DONE" << std::endl;
		for (Vertex& vertex : points) vertex.TypeCheck();
	}

	template <class Comp> 
 	std::vector<Edge> EdgeList(Comp comp) const {
		std::vector<Edge> edges;
		for (const Vertex & point : points) {
			for (VertexRef neighbor : point.Neighbors()) {
				if (comp(point, neighbor)) {
					edges.push_back(Edge(point, neighbor.get()));
				}
			}
		}
		std::sort(edges.begin(), edges.end(), comp);
		return edges;
	}
};

template <class VT> 
class EdgeTerrain {
	std::vector<VT> points;
	std::vector<size_t> n;

	public:
	typedef VT Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef VertexPair<const Vertex> Edge;
	typedef VertexRefLess<Vertex> VertexHeightComp;
	typedef typename Vertex::HeightType HeightType;
	
	void FindBoundary(std::vector<VertexRef> &boundary, VertexRef start) {
		VertexRef cur = start, prev = start;
		boundary.push_back(cur);
		cur = cur.get().Neighbors().back();
		while (cur.get() != start.get()) {
			boundary.push_back(cur);
			const std::vector<VertexRef> nei = cur.get().Neighbors();
			for (size_t i = 0 ; i < nei.size() ; ++i) {
				if (nei[i].get() == prev.get()) {
					prev = cur;
					cur = nei[(i-1+nei.size())%nei.size()];
					break;
				}
			}
		}
//		std::cout << "Boundary Size: " << boundary.size() << std::endl;
	}

	EdgeTerrain() {}

	EdgeTerrain(const std::vector<Vertex> & vtx, const std::vector< std::pair<size_t,size_t> > & edges ) {
		Build(vtx, edges);
	}
	
	void Build(const std::vector<Vertex> & vtx, const std::vector< std::pair<size_t,size_t> > & edges ) {
		points = vtx;
		n.clear();
		n.resize(4,0);
		// Insert Inf Vertex
		points.push_back(Vertex::InfVertex());

		// Edge connection
		for (const auto & edge : edges) {
			points[edge.first].Connect(points[edge.second]);
			points[edge.second].Connect(points[edge.first]);
		}
		
		VertexRef v_small = points[0];
		for (auto & point: points) {
			if (point == points.back()) break;
			point.SortNeighbor();
			if (!typename Vertex::PositionComp()(v_small, point)) v_small = point;
		}
		std::vector<VertexRef> boundary;
		FindBoundary(boundary, v_small);
		
		for (auto & point: points) {
			for (size_t i = 0 ; i < boundary.size() ; ++i) {
				if (boundary[i].get() == point) {
					point.ConnectBoundary(points.back(), boundary[(i-1+boundary.size())%boundary.size()], boundary[(i+1)%boundary.size()]);
					points.back().Connect(point);
					break;
				}
			}
			point.TypeCheck();
			n[point.Type()]++;
		}
//		if (n[1]+2 != n[2]+n[3]) std::cerr << "Lemma doesn't hold" << std::endl;
//		for (auto _n : n) std::cerr << _n << " ";
//		std::cerr << std::endl;
	}

	template <class Comp> 
 	std::vector<Edge> EdgeList(Comp comp) const {
		std::vector<Edge> edges;
		for (const Vertex & point : points) {
			for (VertexRef neighbor : point.Neighbors()) {
				if (comp(point, neighbor)) {
					edges.push_back(Edge(point, neighbor.get()));
				}
			}
		}
		std::sort(edges.begin(), edges.end(), comp);
		return edges;
	}

	std::vector<VT> & PointList() {
		return points;
	}

	const std::vector<VT> & PointList() const {
		return points;
	}

};

#endif /* __TERRAIN_H__ */
