#ifndef __DYNAMIC_CONTOUR_TREE_INTERFACES_H__
#define __DYNAMIC_CONTOUR_TREE_INTERFACES_H__

#include "./dynamic-tree/euler_tree.h"
#include "./dynamic-tree/link_cut_tree.h"
#include "./dynamic-tree/navigator.h"

template <class CT>
struct BasicContourTree {
	typedef CT Tree;
	typedef typename Tree::ItemType ItemType;
	typedef typename Tree::Node Node;
	typedef typename std::pair<Node *, Node *> Edge;
	Tree tree;

	virtual void Disconnect(Node *, Node *) = 0;
	virtual void Connect(Node *, Node *) = 0;
	virtual Node * Add(ItemType) = 0;
	virtual void Remove(Node *) = 0;

	virtual void Shift(Node *v, ItemType to) {
		v->key = to;
	}

	virtual Node *FindLCA(Node *a, Node *b) {
		return tree.FindLCA(a, b);
	}

	virtual Node *Parent(Node * v) {
		return tree.Parent(v);
	}

	virtual void Evert(Node * v) {
		tree.Evert(v);
	}

	virtual size_t Size() {
		return tree.Size();
	}

	// returns true if t is on the path between f and t in contour tree
	// Unnecessary for dynamic contour tree implementation
//	virtual bool IsOn(Node *d, Node *f, Node *t) = 0;

	// Return Edge of Contour Tree carrying a Regular vertex v
	// If v is not regular then behavior is undefined
	virtual Node* NextNodeTo(Node *, Node *) = 0;

	// Return Edge of Contour Tree carrying a Regular vertex v
	// If v is not regular then return edge is undefined
//	virtual Edge FindEdge(Node *, Node *, ItemType) = 0;

	virtual bool SanityCheck() { return true; }
};

template <class EulerTreeType>
struct ETContourTreeBasic : public BasicContourTree< EulerTreeType > {
	typedef EulerTreeType ET;
	typedef typename ET::ItemType VertexType;
	typedef typename ET::Edge ETEdge;
	typedef BasicContourTree<ET> CT;
	typedef typename CT::Node Node;
	typedef typename CT::Edge Edge;
	typedef typename std::unordered_map<Node *, ETEdge> NE_Type;
	std::unordered_map < Node *, NE_Type > edges;
	
	void Disconnect(Node *a, Node *b) {
		assert(edges.at(a).count(b));
		if (CT::tree.EVERTABLE) CT::tree.Evert(b);
		CT::tree.Cut(edges.at(a).at(b));
		edges.at(a).erase(b);
		edges.at(b).erase(a);
		assert(SanityCheck());
	}

	void Connect(Node *a, Node *b) {
		assert(CT::tree.FindRoot(a) == a);
		ETEdge e = CT::tree.Link(a, b);
		assert(edges.at(a).count(b) == 0);
		assert(edges.at(b).count(a) == 0);
		edges.at(a).insert(std::make_pair(b, e));
		edges.at(b).insert(std::make_pair(a, e));
		assert(SanityCheck());
	}

	Node * Add(VertexType v) {
		Node *n = CT::tree.Add(v);
		edges.insert(std::make_pair(n, NE_Type(0)));
		return n;
	}

	void Remove(Node *v) {
		CT::tree.Remove(v);
		edges.erase(v);
	}

	Node * NextNodeTo(Node *f, Node *t) {
		CT::Evert(t);
		return CT::tree.Parent(f);
	}

	template <typename Comp>
	Edge FindEdge(Node *v_max, Node *v_min, VertexType v, const Comp & comp = Comp()) {
		// This is very Naive Algorithm O(log n + height of the tree)
		CT::Evert(v_max);
		Node *cur = v_min, *p = CT::tree.Parent(v_min);
		while (p != v_max) {
			if (comp(v, p->key)) break;
			cur = p;
			p = CT::tree.Parent(cur);
		}
		return Edge(cur, p);
	}

	bool SanityCheck() {
		for (auto v : edges) {
			for (auto e : v.second) {
				assert(e.second->key.node == v.first || e.second->key.node == e.first);
			}
		}
		return true;
	}
};

template <typename VertexType>
using ETContourTree = ETContourTreeBasic < EulerTree< VertexType > >;

template <class VertexType>
struct LCAETContourTree : public ETContourTreeBasic < LCAEulerTree<VertexType> > {
	typedef LCAEulerTree<VertexType> ET;
	typedef BasicContourTree<ET> CT;
	typedef ETContourTreeBasic < LCAEulerTree<VertexType> > BT;
	typedef typename CT::Node Node;
	typedef typename CT::Edge Edge;
	typedef typename ET::Edge ETEdge;
	typedef typename ET::STNode STNode;
	typedef typename ET::STKey STKey;
	typedef typename ET::ST ST;
	typedef typename ET::Stat Stat;

	bool orientation = true;

	void Disconnect(Node *a, Node *b) {
		orientation = (b == CT::tree.Parent(a));
		BT::Disconnect(a, b);
	}

	void Connect(Node *a, Node *b) {
		if (CT::tree.Parent(a) || (!CT::tree.Parent(a) && !CT::tree.Parent(b) && !orientation)) {
			assert(CT::tree.FindRoot(b) == b);
			std::swap(a, b);
		}
		BT::Connect(a, b);
	}

	struct FindIntervalCondition {
		const Node * t;
		FindIntervalCondition (const Node * t) : t(t) {}
		bool operator()(const Stat & test_stat, const Stat & stat) {
			return ((STKey *)test_stat.key)->node == t;
		}
	};

	ETEdge FindInterval(Node *anc, Node *des) {
		PredNavigator<FindIntervalCondition, STNode> nav((FindIntervalCondition(anc)));
		ETEdge cur = ST::NavigatorSearch(nav, des->repr, anc->repr);
		// Amortization
		// amortized cost of cur is paid by navigator search function
		ST::SplayNode(des->repr); 
		assert(ST::InOrder(des->repr, cur->key.next));
		return cur;
	}

	template <bool ASC, class Comp>
	struct FindEdgeCondition {
		const VertexType v;
		const Comp & comp;
		FindEdgeCondition () {}
		FindEdgeCondition (VertexType v, const Comp & comp) : v(v), comp(comp) {}
		bool operator()(const Stat & test_stat, const Stat & stat) const {
			if (ASC) return (test_stat.Dist(stat) > 0 && 
								(comp(v, ((STKey *)test_stat.key)->node->key)));
			else return (test_stat.Dist(stat) > 0 && 
							(comp(((STKey *)test_stat.key)->node->key, v)));
		}
	};

	Node * NextNodeTo(Node *f, Node *t) {
		Node * lca = CT::tree.FindLCA(f, t);
		if (lca != f) return CT::tree.Parent(f);
		return (ST::Succ(FindInterval(f, t))->key.node);
	}

	template <class Comp>
	Edge FindEdge(Node *v_max, Node *v_min, VertexType v, const Comp & comp = Comp()) {
		Node *lca = CT::tree.FindLCA(v_min, v_max);
		if (comp(lca->key, v.get())) { // target is between lca to v_max
			PredNavigator< FindEdgeCondition<false, Comp>, STNode > nav((FindEdgeCondition<false, Comp>(v, comp)));
			ETEdge v = ST::NavigatorSearch(nav, v_max->repr, FindInterval(lca, v_max));
			ST::SplayNode(v_max->repr); // Amortization
			return Edge(v->key.node, NextNodeTo(v->key.node, v_max));
		} else { // target is between v_min and lca
			PredNavigator< FindEdgeCondition<true, Comp>, STNode > nav((FindEdgeCondition<true, Comp>(v, comp)));
			ETEdge u = ST::NavigatorSearch(nav, v_min->repr, FindInterval(lca, v_min));
			ST::SplayNode(v_min->repr); // Amortization
			return Edge(NextNodeTo(u->key.node, v_min), u->key.node);
		}
	}
};

template <class VertexType>
struct STContourTree : public BasicContourTree< LinkCutTree<VertexType, Statistic, true> > {
	typedef BasicContourTree< LinkCutTree<VertexType, Statistic, true> > CT;
	typedef typename CT::Node Node;
	typedef typename CT::Edge Edge;

	void Disconnect(Node *a, Node *b) {
		CT::tree.Evert(b);
		CT::tree.Cut(a);
	}
	
	void Connect(Node *f, Node *t) {
		assert(CT::tree.FindRoot(f) == f);
		CT::tree.Link(f, t);
	}

	Node * Add(VertexType v) {
		return CT::tree.Add(v);
	}
	
	void Remove(Node * v) {
		CT::tree.Remove(v);		
	}
	
	bool IsOn(Node *d, Node *f, Node *t) {
		// With Evert
		CT::tree.Evert(f);
		return CT::tree.FindLCA(d, t) != f;

		// Without Evert
		/*
		Node *lca = tree.FindLCA(f, t);
		Node *lcaf = tree.FindLCA(d, f);
		Node *lcat = tree.FindLCA(d, t);
		return (lcat == lca && lcaf == d) || (lcaf == lca && lcat == d);
		*/
	}
	Node* NextNodeTo(Node *f, Node *t) {
		if (f == t) return f;
		CT::tree.Evert(t);
		return CT::tree.Parent(f);
	}

	template <class Comp>
	Edge FindEdge(Node *v_max, Node *v_min, VertexType v, const Comp & comp = Comp()) {
		Node *cur = v_min, *low = v_min;
		CT::tree.Evert(v_max);
		assert(CT::tree.FindRoot(v_min) == v_max);
		// Binary Search to Find The Vertex and Make Edge
		CT::tree.Access(low);
		while(cur) {
			CT::tree.PushReverse(cur);
			if (cur->key.get() == v.get()) {
				low = cur;
				break;
			}
			if (comp(cur->key, v)) { // Left is higher
				assert(comp(low->key, cur->key) || low->key.get() == cur->key.get());
				low = cur;
				cur = cur->Left();
			} else {
				cur = cur->Right();
			} 
		}
		assert(low->key.get() != v.get());
		assert(comp(v, NextNodeTo(low, v_max)->key));
		return Edge(low, NextNodeTo(low, v_max));
	}
};

#endif /* __DYNAMIC_CONTOUR_TREE_INTERFACES_H__ */
