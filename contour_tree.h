#ifndef __CONTOUR_TREE_H__
#define __CONTOUR_TREE_H__

#include <set>
#include <unordered_map>
#include <functional>
#include <stack>

#include "union_find.h"
#include "terrain.h"
#include "vertex.h"

template <class T, class Comp = std::less<T> >
struct RefComp {
	bool operator()(const std::reference_wrapper<T> &v1, const std::reference_wrapper<T> &v2) const {
		return Comp()(v1.get(),v2.get());
	}
};

template <class Vertex, class VertexComp = std::less<Vertex> >
class MergeTreeHelper {
	public:
	typedef std::reference_wrapper<const Vertex> VertexRef;
	typedef VertexPair<const Vertex> Edge;
	const Vertex &prev (const Edge &edge) const {
		return edge.u;
	}

	const Vertex &next (const Edge &edge) const {
		return edge.v;
	}

	bool operator()(const Vertex &v1, const Vertex &v2) const {
		return VertexComp()(v1, v2);
	}

	bool operator()(const Edge &e1, const Edge &e2) const {
		if (e1.v.get() == e2.v.get()) return VertexComp()(e1.u, e2.u);
		return VertexComp()(e1.v, e2.v);
	}

	static const std::string Name() {
		return "Merge Tree";
	}

	static const VertexType LeafType = MINIMUM;
};

template <class Vertex, class VertexComp = std::less<Vertex> >
class SplitTreeHelper {
	public:
	typedef std::reference_wrapper<const Vertex> VertexRef;
	typedef VertexPair<const Vertex> Edge;
	const Vertex &prev (const Edge &edge) const {
		return edge.u;
	}

	const Vertex &next (const Edge &edge) const {
		return edge.v;
	}

	bool operator()(const Vertex &v1, const Vertex &v2) const {
		return VertexComp()(v2, v1);
	}

	bool operator()(const Edge &e1, const Edge &e2) const {
		if (e1.v.get() == e2.v.get()) return VertexComp()(e2.u, e1.u);
		return VertexComp()(e2.v, e1.v);
	}

	static const std::string Name() {
		return "Split Tree";
	}
	
	static const VertexType LeafType = MAXIMUM;
};

template <class Terrain, class Helper>
class OneWayTree {
	private:
	typedef typename Terrain::Edge Edge;
	typedef typename Terrain::Vertex Vertex;
	typedef typename Terrain::VertexRef VertexRef;
	typedef typename Terrain::VertexHeightComp VertexHeightComp;

	struct Node {
		VertexRef v;
		typedef typename std::reference_wrapper<Node> NodeRef;
		NodeRef p;
		std::set< NodeRef, RefComp<Node> > children;
		Node (VertexRef v) : v(v), p(std::ref(*this)) {
			children.clear();
			assert(v.get() == p.get().v.get());
		}
		Node (const Node &n) : v(n.v), p(std::ref(*this)), children(n.children) {
			if (!n.IsRoot()) p = n.p;
		}
		bool operator < (const Node & n) const {
			return VertexHeightComp()(v, n.v);
		}
		bool IsRoot() const {
			return p.get().v.get() == v.get();
		}
	};

	typedef typename Node::NodeRef NodeRef;
	typedef typename std::unordered_map<VertexRef, Node, typename Vertex::Hash, VertexRefEqual<Vertex> > VNMap;
	typedef typename VNMap::iterator VNMapIt;
	typedef typename std::set<VertexRef, RefComp<const Vertex, VertexHeightComp> > VertexRefSet;
	typedef typename std::set<VertexRef, RefComp<const Vertex, VertexHeightComp> > VsIt;

	VNMap VertexToNode;
	VertexRefSet leaves;

	Node &MakeNode(VertexRef v) {
		VNMapIt it = VertexToNode.find(v);
		if (it == VertexToNode.end()) {
			VertexToNode.insert(std::make_pair(v,Node(v)));
			it = VertexToNode.find(v);
			assert(it != VertexToNode.end());
			assert(it->second.IsRoot());
		}
		return it->second;
	}
	
	// For construction
	void Link(Node &u, Node &v) {
		assert(u.p.get().v.get() == u.v.get());
		u.p = std::ref(v);
		v.children.insert(std::ref(u));
	}
	
	public:
	OneWayTree(const Terrain &terrain) {
		Helper helper;
		UnionFind<VertexRef, typename Vertex::Hash, VertexRefEqual<Vertex> > uf;
	
		for (auto edge : terrain.EdgeList(Helper())) {
			VertexRef prev = helper.prev(edge);
			VertexRef next = helper.next(edge);
			uf.Add(prev);
			uf.Add(next);
//			std::cout << prev.get() << " (" << prev.get().Type() << ") - " << next.get() << " (" << next.get().Type() << ")" << std::endl;
			assert(uf.Root(prev).get().Type() != REGULAR);
			if (prev.get().Type() != REGULAR) {
				MakeNode(prev);
				if (IsLeaf(prev)) leaves.insert(prev);
			} 

			if (next.get().Type() == REGULAR) {
				uf.Join(next, uf.Root(prev));
//				std::cout << next.get() << " -> " << uf.Root(prev).get() << " : " << uf.Root(next).get() << std::endl;
			} else {
				VertexRef repr = uf.Root(prev);
//				std::cout << repr.get() << " [" << prev.get() << "] " << " -> " << next.get() << std::endl;
				assert(repr.get().Type() != REGULAR);
				assert(next.get().Type() != REGULAR);
				assert(next.get() == uf.Root(next).get());
				uf.Join(repr, next);
				assert(VertexToNode.find(repr) != VertexToNode.end());
				if (repr.get() != next.get()) Link(VertexToNode.at(repr), MakeNode(next));
			}
		}

//		std::cerr << Helper::Name() << " has been constructed " << std::endl;
	}

	OneWayTree(const OneWayTree &tree) : VertexToNode(tree.VertexToNode), leaves(tree.leaves) {
		for (std::pair<const VertexRef, Node> & p : VertexToNode) {
			std::set< NodeRef, RefComp<Node> > new_children;
			Node &node = p.second;
			for (const NodeRef & child : node.children)
				new_children.insert(std::ref(VertexToNode.at(child.get().v)));
			std::swap(p.second.children, new_children);
			node.p = std::ref(VertexToNode.at(node.p.get().v));
		}
	}

	~OneWayTree() {
		VertexToNode.clear();
	}

	// For deconstruction
	void Remove(VertexRef v) {
		Node &u = VertexToNode.at(v);
		if (!u.IsRoot()) {
			assert(u.children.size() < 2); // should have at most 1 child
			Node &p = u.p;
			p.children.erase(std::ref(u)); 
			if (!u.children.empty()) { // the case where v is degree 2
				Node &c = *u.children.begin();
				c.p = p;
				p.children.insert(c);
			} else { // the case where v is a leaf
				leaves.erase(v);
				if (p.children.empty()) leaves.insert(p.v);
			}
		} else { // This shouldn't be happen
			// the case where v is a root
			// make all children be roots
			if (!u.children.empty())
				for (Node &node : u.children) 
					node.p = node;
			else 
				leaves.erase(v);
		}
		VertexToNode.erase(v);
	}

	bool Empty() const {
		return Size() == 0;
	}

	bool IsRoot(VertexRef v) const {
		return VertexToNode.at(v).IsRoot();
	}

	bool IsLeaf(VertexRef v) const {
		return VertexToNode.at(v).children.empty();
	}

//	std::vector<VertexRef> Leaves() const {
//		std::vector<VertexRef> V;
//		for (const auto &leaf : leaves)
//			V.push_back(*leaf);
//		return V;
//	}

	const VertexRefSet & Leaves() const {
		return leaves;
	}

	VertexRef Parent(VertexRef v) const {
		return VertexToNode.at(v).p.get().v;
	}

	size_t Size() const {
		return VertexToNode.size();
	}

	size_t ChildrenSize(VertexRef v) const {
		return VertexToNode.at(v).children.size();
	}

	const std::vector< VertexRef > Children(VertexRef v) {
		std::vector <VertexRef> V;
		const std::set< NodeRef, RefComp<Node> > & c = VertexToNode.at(v).children;
		V.reserve(c.size());
		for (const Node &n : c)
			V.push_back(n.v);
		return V;
	}

	const std::vector< VertexRef > NodeList() const {
		std::vector <VertexRef> V;
		V.reserve(VertexToNode.size());
		for (const std::pair<VertexRef, Node> & v : VertexToNode)
			V.push_back(v.first);
		return V;
	}

	void Print() const {
		std::cout << "Size: " << Size() << std::endl;
		if (!Size()) return;
		std::reference_wrapper<const Node> n = VertexToNode.begin()->second;
		while (!n.get().IsRoot()) n = std::cref(n.get().p);
		std::stack< std::pair<std::reference_wrapper<const Node>,int> > S;
		S.push(std::make_pair(n,0));
		while(!S.empty()) {
			auto t = S.top();
			S.pop();
			for (int i = 0 ; i < t.second ; ++i) std::cout << "   ";
			std::cout << t.first.get().v.get() << " (" << t.first.get().children.size() << ")" << std::endl;
			for (auto child : t.first.get().children) 
				S.push(std::make_pair(child, t.second+1));
		}
	}

	static const std::string Name() {
		return Helper::Name();
	}
};

//Alias
template <class Terrain>
using MergeTree = OneWayTree<Terrain, MergeTreeHelper<typename Terrain::Vertex, typename Terrain::VertexHeightComp> >;

template <class Terrain>
using SplitTree = OneWayTree<Terrain, SplitTreeHelper<typename Terrain::Vertex, typename Terrain::VertexHeightComp> >;

template <class DT, class Terrain>
class ContourTree {
	DT tree;
	public:
	typedef typename Terrain::Vertex Vertex;
	typedef typename Terrain::VertexRef VertexRef;
	typedef typename Terrain::VertexHeightComp VertexHeightComp;
	typedef typename DT::Node Node;
	typedef typename std::unordered_map<VertexRef, Node *, typename Vertex::Hash, VertexRefEqual<Vertex> > VNMap;
	typedef typename VNMap::iterator VNMapIt;
	typedef typename std::set<VertexRef, RefComp<const Vertex, VertexHeightComp> > LeafSet;

	VNMap VertexToNode;
	Node *gmn;

	template <class T1, class T2>
	void ExtractPossibleLeaves(LeafSet & leaves, const T1 & t1, const T2 & t2) {
		for (const auto & leaf : t1.Leaves()) {
			size_t size = t2.ChildrenSize(leaf);
			assert(size != 0 || t1.Size() == 1);
			if (size == 1) leaves.insert(leaf);
		}
	}

	template <class T1, class T2>
	void RemoveLeaf(LeafSet & leaves, T1 & t1, T2 & t2) {
		VertexRef leaf = *leaves.begin();
		// Assume leaf is a leaf of t1
		assert(t1.IsLeaf(leaf));
		bool t1_root = t1.IsRoot(leaf), t2_root = t2.IsRoot(leaf);
		std::vector<VertexRef> cand;
		if (!t1_root) {
			cand.push_back(t1.Parent(leaf));
			tree.Link(VertexToNode.at(leaf), VertexToNode.at(cand.back()));
			if (!t2_root)
				cand.push_back(t2.Parent(leaf));
		} else 
			gmn = VertexToNode.at(leaf);

		t1.Remove(leaf);
		t2.Remove(leaf);
	
		leaves.erase(leaves.begin());
		for (auto p : cand) {
			if ((t1.IsLeaf(p) && t2.ChildrenSize(p) == 1) 
				|| (t2.IsLeaf(p) && t1.ChildrenSize(p) == 1)) {
				leaves.insert(p);
			}
		}
	}

	ContourTree (MergeTree<Terrain> mt, SplitTree<Terrain> st) {
		// Add Nodes as singletons
		for	(VertexRef v : mt.NodeList()) 
			VertexToNode.insert(std::make_pair(v, tree.Add(v)));
		
		// Prepare leaf extraction
		LeafSet leaves;
		ExtractPossibleLeaves(leaves, mt, st);
		ExtractPossibleLeaves(leaves, st, mt);
		
		while(!mt.Empty()) {
			assert(!leaves.empty());
			if (mt.IsLeaf(*leaves.begin()))
				RemoveLeaf(leaves, mt, st); 
			else 
				RemoveLeaf(leaves, st, mt);
		}
		
//		std::cerr << "Contour Tree has been constructed" << std::endl;
	}

	ContourTree(const Terrain & terrain) : ContourTree(MergeTree<Terrain>(terrain), SplitTree<Terrain>(terrain)) {}

	void Print() {
		std::cout << "Contour Tree Size: " << VertexToNode.size() << std::endl;
		tree.Evert(gmn);
		for (const std::pair<const VertexRef, Node *> & p : VertexToNode) {
			std::cout << p.first.get();
			Node *parent = tree.Parent(p.second);
			if (parent) std::cout << " -> " << tree.Parent(p.second)->key.get() << std::endl;
			else std::cout << " is root" << std::endl;
		}
	}

	void Evert(VertexRef v) {
		tree.Evert(VertexToNode.at(v));
	}

	VertexRef Parent(VertexRef v){
		Node *p = tree.Parent(VertexToNode.at(v));
		if (p) return p->key;
		return v;
	}

	const VNMap & VertexNodeMap() {
		return VertexToNode;
	}
};

#endif /* __CONTOUR_TREE_H__ */
