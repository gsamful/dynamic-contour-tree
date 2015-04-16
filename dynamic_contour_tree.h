#ifndef __DYNAMIC_CONTOUR_TREE_H__
#define __DYNAMIC_CONTOUR_TREE_H__

#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>

#include "vertex.h"
#include "contour_tree.h"

#include "dynamic_contour_tree_interfaces.h"

// TODO: remove get() usage by creating proper operators
/*
 * Assumption:
 * Vertex has Neighbors(), Successor(v), Predecessor(v)
 * Vertex's can be compared by (<, > :height) (==, != :identity) operators
 * Vertex's Hash Function doesn't use height attribute
 * Dynamic Tree data structure for contour tree has Evert() only for test
 */

template <class DynamicTreeForContourTree, class DynamicTreeForAuxTree, class Terrain> 
class DynamicContourTreeBasic {
public:
	DynamicTreeForContourTree contour_tree;
	DynamicTreeForAuxTree asc_tree, des_tree;

	typedef typename DynamicTreeForContourTree::Node Node;
	typedef typename DynamicTreeForContourTree::Edge Edge;
	typedef typename DynamicTreeForAuxTree::Node AuxNode;

	typedef typename Terrain::Vertex Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef typename std::reference_wrapper<Vertex> MutableVertexRef;
	typedef VertexRefEqual<Vertex> EqualVertexRef;
	typedef VertexRefLess<Vertex> LessVertexRef;

	typedef typename std::unordered_set<VertexRef, typename Vertex::Hash, EqualVertexRef> LinkType;
	typedef typename std::unordered_set<Node *> ConnectionType;
//	typedef typename std::vector<VertexRef> LinkType;
//	typedef typename std::vector<Node *> ConnectionType;
	typedef DynamicContourTreeBasic<DynamicTreeForContourTree, DynamicTreeForAuxTree, Terrain> DCT;

	typedef typename Vertex::HeightType HeightType;

	DynamicContourTreeBasic() { // need a constructor accepting contour tree and terrain
	}
	
	DynamicContourTreeBasic(Terrain & terrain) { Build(terrain); ResetEventCounter(); }
	~DynamicContourTreeBasic() {}

	void Build(const Terrain & terrain) {
		ContourTree<DynamicTreeForAuxTree, Terrain> ct(terrain);
		for (auto vn : ct.VertexNodeMap()) {
			Add(vn.first);
		}
		for (auto vn : VertexToNode) {
			VertexRef p = ct.Parent(vn.first);
			if (!EqualVertexRef()(p, vn.first)) Connect(vn.second, VertexToNode.at(p));
		}
		BuildAuxData(terrain);
//		std::cerr << "Dynamic Contour Tree has been constructed" << std::endl;
	}

	void ResetEventCounter() {
		LocalEventCnt = InterchangeEventCnt = BirthEventCnt = DeathEventCnt = ShiftEventCnt = AuxEventCnt = 0;
	}

	/***************************
	 * Functions for test case *
	 ***************************/

	void Evert(VertexRef v) {
		contour_tree.Evert(VertexToNode.at(v));
	}

	VertexRef Parent(VertexRef v) {
		Node *p = contour_tree.Parent(VertexToNode.at(v));
		if (p) return p->key;
		return v;
	}


	/***********************
	 *  Testing Functions  *
	 ***********************/

	bool EdgeTest(Edge e, VertexRef v) {
		bool ret = true;
//		std::cout << "Testing : " << v.get() << " : ";
		if (IsRegular(v)) {
			assert(e.first && e.second);
//			std::cout << " (Regular) ";
			Node *v_min = VertexToNode.at(DesRoot(v));
			Node *v_max = VertexToNode.at(AscRoot(v));
			Node *lca = contour_tree.FindLCA(v_min, v_max);
			Node *lcaf = contour_tree.FindLCA(v_min, e.first);
			Node *lcat = contour_tree.FindLCA(v_max, e.first);
			ret = ((lca == lcaf && e.first == lcat) || (lca == lcat && e.first == lcaf));
//			std::cout << e.first->key.get() << " " << ret << " -> ";
			lcaf = contour_tree.FindLCA(v_min, e.second);
			lcat = contour_tree.FindLCA(v_max, e.second);
			ret &= ((lca == lcaf && e.second == lcat) || (lca == lcat && e.second == lcaf));
//			std::cout << e.second->key.get() << " " << ret;
			assert(LessVertexRefComp(e.first->key, e.second->key));
		} else {
//			std::cout << " (Saddle) ";
//			if (e.first) std::cout << e.first->key.get() << " -> ";
//			if (e.second) std::cout << e.second->key.get();
			ret &= ((e.first && EqualVertexRef()(e.first->key, v)) || (e.second && EqualVertexRef()(e.second->key, v)));
		}
//		std::cout <<std::endl;
		return ret;
	}

	bool SanityCheck() {
		size_t e = 0;
		for (auto c : Connection) {
			e += c.second.size();
		}
		assert(contour_tree.Size() * 2 == (e + 2));
		for (auto c : Connection) {
			contour_tree.Evert(c.first);
			for (auto n : c.second) {
				assert(contour_tree.Parent(n) == c.first);
			}
		}
		return true;
	}

	void PrintPath(VertexRef v) {
		Node *v_min = VertexToNode.at(DesRoot(v)), *v_max = VertexToNode.at(AscRoot(v));
		while (v_min != v_max) {
//			std::cout << v_min->key.get() << " -> ";
			assert(LessVertexRefComp(v_min->key, NextNodeTo(v_min, v_max)->key));
			v_min = NextNodeTo(v_min, v_max);
		}
//		std::cout << v_max->key.get() << std::endl;
	}

	const ConnectionType Neighbors(VertexRef v) {
		return Neighbors(VertexToNode.at(v));
	}

	//TODO: Tools for traverse contour tree


public:
	/*******************************
	 * Operations for contour tree *
	 *******************************/
	void Disconnect(Node *a, Node *b) {
		size_t size1 = Connection.at(a).size();
		size_t size2 = Connection.at(b).size();
		contour_tree.Disconnect(a, b);
		Remove(Connection.at(a), b);
		Remove(Connection.at(b), a);
		assert(Connection.at(a).size() == size1 - 1);
		assert(Connection.at(b).size() == size2 - 1);
	}
	
	void Connect(Node *f, Node *t) {
		for (auto c : Connection.at(f)) assert(c != t);
		for (auto c : Connection.at(t)) assert(c != f);
		contour_tree.Connect(f, t);
		Insert(Connection.at(f), t);
		Insert(Connection.at(t), f);
	}
	
	Node * Add(VertexRef v) {
		Node *ret = contour_tree.Add(v);
		VertexToNode.insert(std::make_pair(v, ret));		
		Connection.insert(std::make_pair(ret, ConnectionType()));
		return ret;
	}
	
	void Remove(Node * v) {
		assert(Neighbors(v).size() == 0);
		Connection.erase(v);
		VertexToNode.erase(v->key);
		contour_tree.Remove(v);		
	}

	bool IsMonkey(Node *v) {
		return Neighbors(v).size() > 3;
	}

	bool IsExtremal(Node *v) {
		return Neighbors(v).size() == 1;
	}

	bool IsSaddle(Node *v) {
		return !IsExtremal(v);
	}

	bool IsRegular(VertexRef v) {
		return VertexToNode.count(v) == 0;
	}

	// Returns all neighbors of a node v
	const ConnectionType Neighbors(Node *v) {
		return Connection.at(v);
	}

	AuxNode * AscNode(VertexRef v) {
		return VertexToAux.at(v).asc;
	}

	AuxNode * DesNode(VertexRef v) {
		return VertexToAux.at(v).des;
	}

	VertexRef AscRoot(VertexRef v) {
		return asc_tree.FindRoot(AscNode(v))->key;
	}

	VertexRef DesRoot(VertexRef v) {
		return des_tree.FindRoot(DesNode(v))->key;
	}

	LinkType & UpperLink(Vertex v) {
		return VertexToAux.at(v).uk;
	}

	LinkType & LowerLink(Vertex v) {
		return VertexToAux.at(v).lk;
	}

	/***************************************
	* Contour Tree Related Operations      *
	* Dynamic tree implementation specific *
	***************************************/
	Node* NextNodeTo(Node *f, Node *t) {
		return contour_tree.NextNodeTo(f, t);
	}

	Edge FindEdge(VertexRef v) {
		return contour_tree.FindEdge(VertexToNode.at(AscRoot(v)), VertexToNode.at(DesRoot(v)), v, LessVertexRefComp);
	}

protected:
	/*****************
	 * Vertex Helper *
	 *****************/
	/* This Helper is for continuous change of a vertex. */
	template <bool UP>
	struct VertexHelper {
		VertexRef v;
		LinkType & link, & opposite_link;
		// Flip history for changing height of a vertex
		std::unordered_set<VertexRef, typename Vertex::Hash, EqualVertexRef> flip_history;
		LessVertexRef comp;

		VertexHelper(VertexRef v, LinkType & link, LinkType & opposite_link, LessVertexRef comp) 
			: v(v), link(link), opposite_link(opposite_link), flip_history(0), comp(comp) {}

		VertexType Type() {
			if (UP) return DCT::Type(link.size(), opposite_link.size());
			return DCT::Type(opposite_link.size(), link.size());
		}

		size_t ComponentSize() const {
			return link.size();
		}

		bool Dir(VertexRef u) {
			bool flip = (flip_history.count(u));
			bool dir = flip?!UP:UP;
			return comp(v,u) == dir;
		}

		void Flip(VertexRef u) {
			assert(Dir(u)); // u will become false direction by this flip
			VertexRef succ_u = v.get().Successor(u), pred_u = v.get().Predecessor(u);
			if (opposite_link.empty()) {
				link.clear();
				Insert(opposite_link, u);
				Insert(link, succ_u);
			} else {
				if (Dir(pred_u)) Insert(opposite_link, u);
				else Remove(link, u, EqualVertexRef());
				if (!Dir(succ_u)) Remove(opposite_link, succ_u, EqualVertexRef());
				else Insert(link, succ_u);
			}
			if (link.empty() && opposite_link.empty()) Insert(opposite_link, u);
			assert((link.size() == opposite_link.size() && link.size() > 0) || (link.size() + opposite_link.size() == 1));
		}

		VertexRef FindAlternative(VertexRef cur) {
			if (EqualVertexRef()(cur, v)) return v;
			for (auto u : link) {
				if (!EqualVertexRef()(u, cur)) return u;
				VertexRef succ = v.get().Successor(u);
				if (Dir(succ)) return succ;
			}
			return v;
		}

	};
	
	template <bool UP>
	inline VertexRef AuxParent(VertexRef v) {
		if (UP) {
			AuxNode *asc_node = AscNode(v);
			VertexRef ret = (asc_tree.IsRoot(asc_node))?v:asc_tree.Parent(asc_node)->key;
			return ret;
		} else {
			AuxNode *des_node = DesNode(v);
			VertexRef ret = (des_tree.IsRoot(des_node))?v:des_tree.Parent(des_node)->key;
			return ret;
		}
	}

	template <bool UP>
	inline VertexHelper<UP> CreateVertexHelper(VertexRef v) {
		if (UP) {
			return VertexHelper<UP>(v, UpperLink(v), LowerLink(v), LessVertexRefComp);
		} else {
			return VertexHelper<UP>(v, LowerLink(v), UpperLink(v), LessVertexRefComp);
		}
	}


	/*****************
	 * Event Handler *
	 *****************/
	void BirthEvent(VertexRef cri, VertexRef sad/*, Edge &ce*/) {
		// This should be done by saddle
		// since the adjacency of critical has been changed.
//		std::cout << " >> Birth Cri:" << cri.get() << " Sad:" << sad.get() << std::endl;
		Node *critical = Add(cri);
		Node *saddle;
		BirthEventCnt++;
		if (!IsRegular(sad)) { // sad is not regular so it becomes a monkey saddle
			saddle = VertexToNode.at(sad);
//			std::cout << "     I'm make a Monkey saddle!" << std::endl;
			Connect(critical, saddle);
//			assert(ce.first == saddle || ce.second == saddle);
		} else {
			Edge ce = FindEdge(sad);
			saddle = Add(sad);
			Node *eu = ce.first;
			Node *ev = ce.second;
			Disconnect(eu, ev);
			Connect(saddle, ev);
			Connect(eu, saddle);
			Connect(critical, saddle);
		}
	}
	
	void DeathEvent(VertexRef cri, VertexRef sad/*, Edge &ce*/) {
//		std::cout << " >> Death Cri:" << cri.get() << " Sad:" << sad.get() << std::endl;
		Node *critical = VertexToNode.at(cri);
		Node *saddle = VertexToNode.at(sad);
		bool monkey = IsMonkey(saddle);
		DeathEventCnt++;
		assert(contour_tree.Parent(critical) == saddle || contour_tree.Parent(saddle) == critical);
		if (!monkey) {
			std::vector<Node *> e;
			for (auto neighbor : Neighbors(saddle)) {
				if (neighbor != critical) e.push_back(neighbor);
			}
			assert(e.size() == 2);
			assert(IsExtremal(critical));
			assert(Neighbors(saddle).size() == 3);
			Disconnect(saddle, e[0]);
			Disconnect(e[1], saddle);					
			Disconnect(critical, saddle);
			Connect(e[1], e[0]);
			Remove(saddle);
			Remove(critical);
		} else {
//			std::cout << "     Monkey Death.." << std::endl;
			assert(IsExtremal(critical));
			Disconnect(critical, saddle);
			Remove(critical);
//			assert(ce.first == saddle || ce.second == saddle);
		}
	}

	// TODO: Make this simpler
	void ShiftEvent(VertexRef cri, VertexRef reg/*, Edge &ce */) {
//		std::cout << " >> Shift Cri:" << cri.get() << " Reg:" << reg.get() << std::endl;
		Node *from = VertexToNode.at(cri);
		size_t deg_f = Neighbors(from).size();
		bool is_reg = IsRegular(reg);
		ShiftEventCnt++;
		// cri = extremal point or cri = saddle and reg = regular
		if (IsExtremal(from) || (!IsMonkey(from) && is_reg)) { 
			contour_tree.Shift(from, reg);
			VertexToNode.insert(std::make_pair(reg, from));
			VertexToNode.erase(cri);
//			assert(from == ce.first || from == ce.second);
		} else if (is_reg) { // cri is a monkey saddle and reg is regular
			Node *eu, *ev;
			if (LessVertexRefComp(cri, reg)) { // Increasing case
				eu = VertexToNode.at(cri);
				ev = NextNodeTo(eu, VertexToNode.at(AscRoot(reg)));
			} else {
				ev = VertexToNode.at(cri);
				eu = NextNodeTo(ev, VertexToNode.at(DesRoot(reg)));
			}
			Node *to = Add(reg);
			Disconnect(eu, ev);
			Connect(eu, to);
			Connect(to, ev);
//			assert(ce.first == from || ce.second == from);
//			assert(LessVertexRefComp(ce.first->key, ce.second->key));
			// Rest of the connectivity will be handled in interchange event
		} else if (!IsMonkey(from) && !is_reg) { // cri becomes regular and target becomes a monkey
			Node *to = VertexToNode.at(reg);
			assert(IsSaddle(to));
			Disconnect(from, to);
			for (auto n : Neighbors(from)) {
				Disconnect(n, from);
				Connect(n, to);
			}
			Remove(from);
		} else { // cri is a saddle and regular is also saddle
			// Should be handled by interchange event
			assert(IsMonkey(from) && IsSaddle(VertexToNode.at(reg)));
//			std::cout << "     From Monkey Saddle to Monkey Saddle. " << Neighbors(from).size() << " " << Neighbors(VertexToNode.at(reg)).size() << std::endl;
		}
	}

	/********************************
	 * Abstraction of Event Handler *
	 ********************************/

	template <bool UP = true>
	void LocalEventHandler(VertexHelper<UP> & ve, VertexHelper<!UP> & ue/*, Edge  &e */) {
		VertexRef v = ve.v, u = ue.v;
		VertexRef v_extreme = AuxParent<UP>(v), v_alter = ve.FindAlternative(v_extreme);
		VertexRef u_extreme = AuxParent<!UP>(u), u_alter = ue.FindAlternative(u_extreme);
		VertexType v_type = ve.Type(), u_type = ue.Type();
		bool birth_death = false, aux_change = false;
		LocalEventCnt++;

		if (!UP) std::swap(v_extreme, u_extreme), std::swap(v_alter, u_alter), std::swap(v_type, u_type), std::swap(v, u);

//		std::cout << "*************** Local Event " << v.get() << " " << v_type << " with " << u.get() << " " << u_type << std::endl;
		assert(LessVertexRefComp(v, u));
		// Update Ascent Tree
		if (EqualVertexRef()(v_extreme, u)) {
			aux_change = true;
			AuxNode * nv = AscNode(v); 
			asc_tree.Cut(nv);
			VertexRef w = v_alter;
			if (!EqualVertexRef()(w, v)) {
				asc_tree.Link(nv, AscNode(w));
			} else if (u_type != MAXIMUM) { // v becomes a maximum and u is not a maximum
				// Birth Event
//				std::cout << "Max birth event" << std::endl;
				BirthEvent(v, u/*, e*/);
				birth_death = true;
			}
		}

		if (u_type == MAXIMUM) {
			aux_change = true;
			asc_tree.Link(AscNode(u), AscNode(v));
			VertexRef w = v_alter;
			if (EqualVertexRef()(w, v)) { // (v, u) is the only up edge
				// Shift Event
				assert(v_type == REGULAR);
//				std::cout << "Max shift event" << std::endl;
				ShiftEvent(u, v/*, e*/);
			} else { // Death Event
				assert(u_type == MAXIMUM && v_type == SADDLE);
//				std::cout << "Max death" << std::endl;
				DeathEvent(u, v/*, e*/);
				birth_death = true;
			}
		}

		// Update Descent Tree
		if (EqualVertexRef()(u_extreme, v)) {
			aux_change = true;
			AuxNode * nu = DesNode(u); 
			des_tree.Cut(nu);
			VertexRef w = u_alter;
			if (!EqualVertexRef()(w, u)) {
				des_tree.Link(nu, DesNode(w));
			} else if (v_type != MINIMUM) { // v becomes a minimum and u is not a maximum
				// Birth Event
//				std::cout << "Min birth event" << std::endl;
				BirthEvent(u, v/*, e*/);
				birth_death = true;
			}
		}

		if (v_type == MINIMUM) {
			aux_change = true;
			des_tree.Link(DesNode(v), DesNode(u));
			VertexRef w = u_alter;
			if (EqualVertexRef()(w, u)) { // (v, u) is the only up edge
				// Shift Event
				assert(u_type == REGULAR);
//				std::cout << "Min shift event" << std::endl;
				ShiftEvent(v, u/*, e*/);
			} else { // Death Event
				assert(v_type == MINIMUM && u_type == SADDLE);
//				std::cout << "Min death" << std::endl;
				DeathEvent(v, u/*, e*/);
				birth_death = true;
			}
		}

		if (!UP) std::swap(v, u);
		size_t c_s = ve.ComponentSize(), n_c_s;
		ve.Flip(u), ue.Flip(v);
		n_c_s = ve.ComponentSize();

		// Saddle shift
		if (!birth_death && c_s != n_c_s && std::max(c_s, n_c_s) >= 2) {
//			std::cout << "    Saddle shift " << c_s << " " << n_c_s << " " << ve.Type() << " " << ue.ComponentSize() << std::endl;
			if (c_s < n_c_s) ShiftEvent(u, v/*, e*/);
			else ShiftEvent(v, u/*, e*/);
		} 

		if (aux_change) AuxEventCnt++;

		assert(VertexToNode.count(AscRoot(ve.v)));
		assert(VertexToNode.count(DesRoot(ve.v)));
		assert((ve.Type() == REGULAR) == IsRegular(ve.v) && (ve.Type() != SADDLE || IsSaddle(VertexToNode.at(ve.v))));
	}

	// Among the neighbors of b find new neighbors of a after interchange event
	std::vector<Node *> ExtractChildren(Node *a, Node *b) {
		bool UP = LessVertexRefComp(a->key, b->key);
		LinkType & link = UP?UpperLink(a->key):LowerLink(a->key);
		std::vector<Node *> ret(0);
		for (auto v : link) {
			VertexRef to = UP?AscRoot(v):DesRoot(v);
			assert(!IsRegular(to));
			Node *child = NextNodeTo(b, VertexToNode.at(to));
			if (child != a && std::find(ret.begin(), ret.end(), child) == ret.end())
				ret.push_back(child);
		}

		if (UP) {
			Node *naa = NextNodeTo(b, VertexToNode.at(AscRoot(a->key)));
			assert(naa == a || std::find(ret.begin(), ret.end(), naa) != ret.end());
		} else {
			Node *nad = NextNodeTo(b, VertexToNode.at(DesRoot(a->key)));
			assert(nad == a || std::find(ret.begin(), ret.end(), nad) != ret.end());
		} 

		return ret;
	}
	
	void InterChangeEventHandler(VertexRef a, VertexRef b, Edge &e) {
//		std::cout << "*************** InterChangeEvent!! " << a.get() << " " << b.get() << std::endl;
		Node *na = VertexToNode.at(a), *nb = VertexToNode.at(b);
		std::vector<Node *> va = ExtractChildren(na, nb);
		std::vector<Node *> vb = ExtractChildren(nb, na);
		InterchangeEventCnt++;

		// Disconnect and Reconnect na nb is unnecessary if dynamic tree is evertable
		Disconnect(na, nb);
		for (auto n : vb) {
			Disconnect(n, na);
			Connect(n, nb);
		}
		for (auto n : va) {
			Disconnect(n, nb);
			Connect(n, na);
		}
		Connect(nb, na);

        // Sanity Check
		assert(Neighbors(na).size() > 2 && Neighbors(nb).size() > 2);
	}

	void BuildAuxData(const Terrain & terrain) {
		for (const auto & v : terrain.PointList()) {
			VertexToAux.insert(std::make_pair(std::cref(v), VertexAux()));
			VertexAux &aux = VertexToAux.at(std::cref(v));
			aux.asc = asc_tree.Add(v), aux.des = des_tree.Add(v);
			auto prev = v.Neighbors().back();
			for (const auto neighbor : v.Neighbors()) {
				if (LessVertexRefComp(v, neighbor) && LessVertexRefComp(prev, v)) Insert(aux.uk, neighbor);
				if (LessVertexRefComp(neighbor, v) && LessVertexRefComp(v, prev)) Insert(aux.lk, neighbor);
				prev = neighbor;
			}
			if (aux.uk.empty() && aux.lk.empty()) {
				if (LessVertexRefComp(v, prev)) Insert(aux.uk, prev);
				else Insert(aux.lk, prev);
			}
		}

		for (const auto & v : terrain.PointList()) {
			VertexAux &aux = VertexToAux.at(std::cref(v));
			if (!aux.uk.empty()) asc_tree.Link(aux.asc, AscNode(*aux.uk.begin()));
			if (!aux.lk.empty()) des_tree.Link(aux.des, DesNode(*aux.lk.begin()));
		}
	}

	VertexType Type(VertexRef v) {
		size_t u_s = UpperLink(v).size(), l_s = LowerLink(v).size();
		return Type(u_s, l_s);
	}

	static VertexType Type(size_t u_s, size_t l_s) {
		if (u_s >= 2) return SADDLE;
		else if (l_s == 0) return MINIMUM;
		else if (u_s == 0) return MAXIMUM;
		return REGULAR;
	}

	template <class T>
	static void Insert(std::vector<T> & vec, T target) {
		vec.push_back(target);
	}

	template <class T, class EQ = std::equal_to<T> >
	static void Remove(std::vector<T> & vec, T target, EQ comp = EQ()) {
		for (auto & e : vec) {
			if (comp(e, target)) {
				std::swap(e,vec.back());
				vec.pop_back();
				break;
			}
		}
	}

	template <class Container, class T>
	static void Insert(Container & con, T target) {
		con.insert(target);
	}

	template <class Container, class T, class EQ = std::equal_to<T> >
	static void Remove(Container & con, T target, EQ eq = EQ()) {
		con.erase(target);
	}

	struct VertexAux {
		AuxNode *asc, *des;
		LinkType uk, lk;
	};

public:
	// Mapping between Critical Vertex and Node
	std::unordered_map<VertexRef, Node *, typename Vertex::Hash, EqualVertexRef > VertexToNode; 
	// Unified Mapping between vertex and the node of aux trees and link pointers
	std::unordered_map<VertexRef, VertexAux, typename Vertex::Hash, EqualVertexRef> VertexToAux;
	// Connectivity for queries
	std::unordered_map< Node *, ConnectionType > Connection;
	LessVertexRef LessVertexRefComp;
	size_t LocalEventCnt, InterchangeEventCnt, BirthEventCnt, DeathEventCnt, ShiftEventCnt, AuxEventCnt;
};

template <class DynamicTreeForContourTree, class DynamicTreeForAuxTree, class Terrain> 
class DynamicContourTree : public DynamicContourTreeBasic<DynamicTreeForContourTree, DynamicTreeForAuxTree, Terrain>{
public:
	typedef DynamicContourTreeBasic<DynamicTreeForContourTree, DynamicTreeForAuxTree, Terrain> DCTB;

	typedef typename DCTB::Node Node;
	typedef typename DCTB::Edge Edge;
	typedef typename DCTB::AuxNode AuxNode;

	typedef typename DCTB::Vertex Vertex;
	typedef typename DCTB::VertexRef VertexRef;
	typedef typename std::reference_wrapper<Vertex> MutableVertexRef;
	typedef VertexRefEqual<Vertex> EqualVertexRef;
	typedef VertexRefLess<Vertex> LessVertexRef;

	typedef typename DCTB::LinkType LinkType;
	typedef typename DCTB::ConnectionType ConnectionType;
	typedef typename DCTB::DCT DCT;
	typedef typename DCTB::HeightType HeightType;
	template <bool UP>
	using VertexHelper = typename DCTB::template VertexHelper<UP>;

	DynamicContourTree() { // need a constructor accepting contour tree and terrain
	}
	
	// -1-base terrain should forms a terrain containing only one min/max
	DynamicContourTree(Terrain & terrain, bool fast_build = true) {
		if (!fast_build) Build_0base(terrain);
		else DCTB::Build(terrain);
	}

	void Build_0base(Terrain & terrain) {
		std::vector<HeightType> R;
		std::vector<VertexRef> cri;

		for (auto & v : terrain.PointList()) {
			if (v == terrain.PointList().back()) break;
			R.push_back(v.Height());
//			v.ChangeHeight(std::numeric_limits<HeightType>::max());
			v.ChangeHeight(HeightType());
		}

		DCTB::BuildAuxData(terrain);
		for (auto & v : terrain.PointList()) {
			v.TypeCheck();
			assert(v.Type() != SADDLE);
			if (v.Type() != REGULAR) {
				cri.push_back(v);
			}
		}

		assert(cri.size() == 2);
		Node *crin1 = DCTB::Add(cri[0]), *crin2 = DCTB::Add(cri[1]);
		if (LessVertexRef()(cri[0], cri[1])) DCTB::Connect(crin1, crin2);
		else DCTB::Connect(crin2, crin1);
//		std::cout << "Initial Contour Tree: ";
//		std::cout << cri[0].get() << " -> " << cri[1].get() << std::endl;

		for (size_t i = 0; i < R.size(); ++i) {
			ChangeHeight(terrain.PointList()[i], R[i]);
		}
//		std::cerr << "Dynamic Contour Tree has been constructed" << std::endl;
	}

	void ChangeHeight(MutableVertexRef v, HeightType r) {
//		std::cout << "******************************************** Change Height: " << v.get() << " >> " << r << std::endl;
		if (v.get().Height() < r) { // Raising
			Update<true>(std::cref(v), r);
		} else if (v.get().Height() > r) { // Lower
			Update<false>(std::cref(v), r);
		}
		v.get().ChangeHeight(r);
	}

private:
	template<bool UP>
	struct VertexComp {
		bool operator()(const VertexRef &v, const VertexRef &w) const {
			if (UP) return LessVertexRef()(w, v);
			else return LessVertexRef()(v, w);
		}
	};

	/********************************
	 * Abstraction of Event Handler *
	 ********************************/

	template <bool UP = true>
	Edge NextEdgeOnContourTree(Node *cur, VertexRef to) { // For a critical point u
		Node *ret = nullptr;
		ret = DCTB::NextNodeTo(cur, (UP)?DCTB::VertexToNode.at(DCTB::AscRoot(to)):DCTB::VertexToNode.at(DCTB::DesRoot(to)));
		if (ret == cur) ret = nullptr;
		return (UP)?Edge(cur, ret):Edge(ret, cur);
	}
	template <bool UP = true>
	Edge NextEdgeOnContourTree(VertexRef v) { // For the critical point v
		assert(DCTB::Type(v) != REGULAR);
		LinkType & link = UP?DCTB::UpperLink(v):DCTB::LowerLink(v);
		Node *nv = DCTB::VertexToNode.at(v);
		Node *t = nullptr;
		for (auto n : link) {
			Node *t2 = DCTB::NextNodeTo(nv, DCTB::VertexToNode.at((UP)?DCTB::AscRoot(n):DCTB::DesRoot(n)));
			if (!t || ((LessVertexRef()(t2->key, t->key)) == UP)) t = t2; 
		}   
		if (UP) {
			if (DCTB::Type(v) != MAXIMUM) return Edge(nv, t); 
			else return Edge(t, nv);
		} else {
			if (DCTB::Type(v) != MINIMUM) return Edge(t, nv); 
			else return Edge(nv, t);
		}
	}

	template <bool UP = true>
	bool CanCross(VertexRef v, VertexRef u, HeightType r) {
		// TODO : More simpler code
		if (UP && (LessVertexRef()(u, v) || u.get().Height() > r)) return false;
		if (!UP && (LessVertexRef()(v, u) || u.get().Height() < r)) return false;
		if (u.get().Height() == r) {
			Vertex v_image = v.get();
			HeightType t = v_image.Height();
			v_image.ChangeHeight(r);
			if (UP && LessVertexRef()(v_image, u)) return false;
			if (!UP && LessVertexRef()(u, v_image)) return false;
		}
		return true;
	}

	template <bool UP = true>
	void Update(VertexRef v, HeightType r) {
		// Insert all neighbors of height > v into event queue
		std::priority_queue<VertexRef, std::vector<VertexRef>, VertexComp<UP> > vq;
		for (auto neighbor : v.get().Neighbors()) 
			if (CanCross<UP>(v,neighbor,r))	vq.push(neighbor);

		// Handle the node in contour tree
		VertexHelper<UP> ve = DCTB::template CreateVertexHelper<UP>(v);
		ve.ComponentSize();
//		std::cout << UP << " Type: " << ve.Type() << " VQ Size: " << vq.size() << std::endl;
		Edge e = (ve.Type() == REGULAR)?DCTB::FindEdge(ve.v):NextEdgeOnContourTree<UP>(ve.v);
		while (true) {
			Node* next = (UP)?e.second:e.first;
//			std::cout << "v = " << v.get() << " " << " next = " << ((!next)?v.get():next->key.get()) << " vq top " << (!vq.empty()?vq.top().get():v.get()) << " " << e.first << " "<< e.second  << std::endl;
			assert(DCTB::EdgeTest(e, v));
			assert(DCTB::contour_tree.SanityCheck());

			// End loop condition
			if ( UP && (!next || (next->key.get().Height() > r)) && vq.empty() ) break;
			if (!UP && (!next || (next->key.get().Height() < r)) && vq.empty() ) break;

			// Local Event
			if (!vq.empty() && (!next
						    || ( UP && !(LessVertexRef()(next->key, vq.top())))
						    || (!UP && !(LessVertexRef()(vq.top(), next->key))))) {
				VertexRef u = vq.top();
				vq.pop();

				bool before_cri = (UP)?EqualVertexRef()(e.first->key, v):EqualVertexRef()(e.second->key, v);
				VertexHelper<!UP> ue = DCTB::template CreateVertexHelper<!UP>(u);
				DCTB::template LocalEventHandler<UP>(ve, ue/*, e*/);
				if (!DCTB::IsRegular(v)) {
					Node *nv = DCTB::VertexToNode.at(v);
					// The case where Death event happened
					if (!DCTB::IsRegular(u) && DCTB::IsSaddle(nv) && DCTB::IsSaddle(DCTB::VertexToNode.at(u)) ) {
						if (UP) e = Edge(nv, DCTB::VertexToNode.at(u));
						else e = Edge(DCTB::VertexToNode.at(u), nv);
					} else {
						e = NextEdgeOnContourTree<UP>(ve.v);
					}
				} else {
					// In this block FindEdge on u is safe since v is unstable
					// Death / Shift / Nothing
					if (before_cri && DCTB::IsRegular(u)) { // Death
						e = DCTB::FindEdge(u);
					} else if (!DCTB::IsRegular(u)) { // Shift to u
						Node *nu = DCTB::VertexToNode.at(u);
						if (!DCTB::IsExtremal(nu)) e = NextEdgeOnContourTree<!UP>(nu, v);
						else { // if u is an extremal there is only one edge connected to u in the contour tree
							e = Edge(nu, *(DCTB::Neighbors(nu).begin()));
							if (LessVertexRef()(e.second->key, e.first->key)) std::swap(e.first, e.second);
						}
					} else {  // Nothing Happen
//						std::cout << "Just nothing happen.. " << e.first->key.get() << " != " << v.get()<< std::endl;
					}
				}

				// Add Flip Histroy
				ve.flip_history.insert(u);
			} else {
				VertexRef u = next->key;
				assert(!DCTB::IsRegular(u));
//				std::cout << "*************** Event with an endpoint of an edge: " << v.get() << " " << ve.Type()  << " with " << u.get() << " " << Type(u) << std::endl;

				// End loop conditions for generate case
				if (!CanCross<UP>(v, u, r)) break;

				// Interchange Event
				if (ve.Type() == SADDLE && DCTB::IsSaddle(next)) {
					assert(!EqualVertexRef()(u, v));
					DCTB::InterChangeEventHandler(ve.v, u, e);
					e = NextEdgeOnContourTree<UP>(ve.v);
				} else if (DCTB::IsExtremal(next)) {
					if (UP) e = Edge(next, nullptr);
					else e = Edge(nullptr, next);
				} else {
//					std::cout << u.get() << " nothing happen " << ve.Type() << " " << UP << std::endl;
					e = NextEdgeOnContourTree<UP>(next, v);
				}

				// Add Flip History
				ve.flip_history.insert(u);
			}
		}
	}
};

template <class DynamicTreeForContourTree, class DynamicTreeForAuxTree, class Terrain> 
class TimeVaryingDynamicContourTree : public DynamicContourTreeBasic<DynamicTreeForContourTree, DynamicTreeForAuxTree, Terrain>{
public:
	typedef DynamicContourTreeBasic<DynamicTreeForContourTree, DynamicTreeForAuxTree, Terrain> DCTB;

	typedef typename DCTB::Node Node;
	typedef typename DCTB::Edge Edge;
	typedef typename DCTB::AuxNode AuxNode;

	typedef typename DCTB::Vertex Vertex;
	typedef typename DCTB::VertexRef VertexRef;
	typedef typename std::reference_wrapper<Vertex> MutableVertexRef;
	typedef VertexRefEqual<Vertex> EqualVertexRef;
	typedef VertexRefLess<Vertex> LessVertexRef;

	typedef typename DCTB::LinkType LinkType;
	typedef typename DCTB::ConnectionType ConnectionType;
	typedef typename DCTB::DCT DCT;
	typedef typename DCTB::HeightType HeightType;
	typedef typename HeightType::InstType HeightInstType;
	typedef typename HeightType::AugType HeightAugType;
private:

	template <bool UP>
	using VertexHelper = typename DCTB::template VertexHelper<UP>;
	struct event {
		VertexRef v, u;
		HeightAugType t;
		size_t vi, ui;
		event(VertexRef v, VertexRef u, HeightAugType t, size_t vi, size_t ui) : v(v), u(u), t(t), vi(vi), ui(ui) {}
		bool operator < (const event & e) const {
			if (e.t == t) return v.get() < e.v.get();
			return e.t < t;
		}
	};
public:

	HeightAugType CrossingTime(VertexRef a, VertexRef b) {
		HeightType diff = a.get().Height() - b.get().Height();
		if (a.get().Height() < b.get().Height()) return std::numeric_limits<HeightAugType>::lowest();
		if (diff.a == 0) return std::numeric_limits<HeightAugType>::lowest();
		return diff.solve();
	}

	TimeVaryingDynamicContourTree() { // need a constructor accepting contour tree and terrain
	}

	// -1-base terrain should forms a terrain containing only one min/max
	TimeVaryingDynamicContourTree(Terrain & terrain) {
		DCTB::Build(terrain);
		DCTB::ResetEventCounter();
		sortedList.reserve(terrain.PointList().size());
		for (const auto & v : terrain.PointList()) {
			if (v == Vertex::InfVertex()) continue;
			sortedList.push_back(v);
//			points.push_back(v);
		}
		DCTB::LessVertexRefComp.t = 0;
//		std::sort(points.begin(), points.end());
		std::sort(sortedList.begin(), sortedList.end(), DCTB::LessVertexRefComp);
		for (size_t i = 1; i < sortedList.size() ; ++i) {
			HeightAugType ct = CrossingTime(sortedList[i-1],sortedList[i]);
			if (ct >= 0) Q.push(event(sortedList[i-1],sortedList[i],ct,i-1,i));
		}

/*
		size_t n = sortedList.size();
		size_t sqrn = sqrt(n);
		assert(sqrn*sqrn == n);
		for (size_t i = 0 ; i < n; i++) {
			size_t r = i/sqrn, c = i%sqrn;
			size_t t;
			if (i % 2 == 0) { // diagonal
				t = i - (sqrn+1);
				if (t < n && t/sqrn+1 == r) edges.push_back(std::make_pair(i, t));
				t = i - (sqrn-1);
				if (t < n && t/sqrn != r) edges.push_back(std::make_pair(i, t));
				t = i + (sqrn-1);
				if (t < n && t/sqrn != r) edges.push_back(std::make_pair(i, t));
				t = i + (sqrn+1);
				if (t < n && t/sqrn-1 == r) edges.push_back(std::make_pair(i, t));
			}
			t = i - sqrn;
			if (t < n) edges.push_back(std::make_pair(i, t));
			t = i - 1;
			if (t < n && t/sqrn == r) edges.push_back(std::make_pair(i, t));
			t = i + 1;
			if (t < n && t/sqrn == r) edges.push_back(std::make_pair(i, t));
			t = i + sqrn;
			if (t < n) edges.push_back(std::make_pair(i, t));
		}   
*/
	}

	void run() {
		size_t n_aug = 0, n_local = 0, n_interchange = 0;
		// Enough constant to handle the event at 0
		HeightAugType time_before_event = -0.000001;
		while(!Q.empty()){
			event e = Q.top();
			Q.pop();
			if (!(EqualVertexRef()(e.v, sortedList[e.vi]) && EqualVertexRef()(e.u, sortedList[e.ui]))) continue;
//			std::cout << "Time t = " << e.t << " " << e.v.get() << "(" << e.v.get().Height()(time_before_event) << ") " << DCTB::Type(e.v) << " " << e.vi << " " << e.u.get() << "(" << e.u.get().Height()(time_before_event) << ")" << " " << DCTB::Type(e.u) << " " << e.ui << std::endl;
			assert(time_before_event < e.t);
			DCTB::LessVertexRefComp.t = time_before_event;
			if (e.v.get().IsNeighbor(e.u)) {
				VertexHelper<true> ve = DCTB::template CreateVertexHelper<true>(e.v);
				VertexHelper<false> ue = DCTB::template CreateVertexHelper<false>(e.u);
				DCTB::template LocalEventHandler<true>(ve, ue);
				n_local++;
			} else if (DCTB::Type(e.v) == REGULAR || DCTB::Type(e.u) == REGULAR) {
				if (DCTB::Type(e.v) == REGULAR && DCTB::Type(e.u) == REGULAR) {
					if (DCTB::FindEdge(e.v) == DCTB::FindEdge(e.u)) n_aug++;
				} else if (DCTB::Type(e.v) == REGULAR) {
					if (DCTB::FindEdge(e.v).second->key.get() == e.u.get()) n_aug++;
				} else {
					if (DCTB::FindEdge(e.u).first->key.get() == e.v.get()) n_aug++;
				}
			}

			if (DCTB::Type(e.v) == SADDLE && DCTB::Type(e.u) == SADDLE) {
				Node *nv = DCTB::VertexToNode.at(e.v), *nu = DCTB::VertexToNode.at(e.u);
				if (DCTB::NextNodeTo(nv, nu) == nu) {
					Edge edge(nv, nu);
					DCTB::InterChangeEventHandler(e.v, e.u, edge);
					n_interchange++;
				}
			} 

			assert(DCTB::contour_tree.SanityCheck());
			std::swap(sortedList[e.vi], sortedList[e.ui]);
			if (e.vi > 0) {
				HeightAugType ct = CrossingTime(sortedList[e.vi-1],sortedList[e.vi]);
				if (ct >= 0) Q.push(event(sortedList[e.vi-1],sortedList[e.vi],ct,e.vi-1,e.vi));
			}
			if (e.ui+1 < sortedList.size()) {
				HeightAugType ct = CrossingTime(sortedList[e.ui], sortedList[e.ui+1]);
				if (ct >= 0) Q.push(event(sortedList[e.ui],sortedList[e.ui+1],ct,e.ui,e.ui+1));
			}

			// Remove same events from the priority queue
			while(!Q.empty() && Q.top().v.get() == e.v.get() && Q.top().u.get() == e.u.get()) Q.pop();
			time_before_event = e.t+(Q.empty()?1:(Q.top().t-e.t)/2);

			/* Debug Code Compare with static contour tree */
/*		
			{
				typedef Vertex2D<double> EV;
				typedef typename std::reference_wrapper<const EV> VF;
				std::vector<Vertex2D<double> > P;
				if (!Q.empty())	assert(Q.top().t-e.t > 0);
				for (const auto p : points) {
					HeightType h = p.Height();
					P.push_back(EV(p.Coord().first, p.Coord().second, h(time_before_event)));
				}
				EdgeTerrain< EV > terrain(P, edges);
				MergeTree<EdgeTerrain< EV > > mt(terrain);
				SplitTree<EdgeTerrain< EV > > st(terrain);
				ContourTree< LinkCutTree<EV, Statistic, true>, EdgeTerrain< EV > > ct(mt, st);
				sort(P.begin(), P.end(), VertexRefLess<EV>());
				for (size_t i = 0 ; i < sortedList.size() ; ++i) {
					assert(sortedList[i].get().Coord() == P[i].Coord());
				}
				assert(Same(ct, *this));
			}
*/
			/* Debug Code Compare with both static contour tree and dynamic contour tree*/
			/*
			{
				typedef Vertex2D<double> EV;
				typedef typename std::reference_wrapper<const EV> VF;
				std::vector<Vertex2D<double> > P;
				for (const auto p : points) {
					HeightType h = p.Height();
					P.push_back(EV(p.Coord().first, p.Coord().second, h(time_before_event)));
				}
				EdgeTerrain< EV > terrain(P, edges);
				DynamicContourTree<STContourTree<VF>, LinkCutTree<VF, Statistic, true>, EdgeTerrain< EV > > dct(terrain);
				for (auto & p : terrain.PointList()) {
					if (p.Coord() == sortedList[e.ui].get().Coord()) {
						HeightInstType ht = sortedList[e.vi].get().Height()(time_before_event)+.5;
						std::cout << "Update " << p << " -> " << ht << std::endl;
						dct.ChangeHeight(p, ht);
					}
				}

				P.clear();
				for (const auto p : points) {
					HeightType h = p.Height();
					P.push_back(EV(p.Coord().first, p.Coord().second, h(time_before_event)));
				}

				EdgeTerrain< EV > terrain2(P, edges);
				MergeTree<EdgeTerrain< EV > > mt(terrain2);
				SplitTree<EdgeTerrain< EV > > st(terrain2);
				ContourTree< LinkCutTree<EV, Statistic, true>, EdgeTerrain< EV > > ct(mt, st);
				sort(P.begin(), P.end(), VertexRefLess<EV>());
				for (auto p : P) {
					std::cout << p << " " << p.Type() << " ";
				}
				std::cout << std::endl;
				for (size_t i = 0 ; i < sortedList.size() ; ++i) {
					std::cout << sortedList[i].get() << " ";
					assert(sortedList[i].get().Coord() == P[i].Coord());
				}
				assert(Same(ct, dct));
				assert(Same(ct, *this));
			}
			*/
		}
//		std::cerr << "Augment Event: " << n_aug << "\tLocal Event: " << n_local << "\tInterchange Event: " << n_interchange << std::endl;
//		std::cerr << "\t" << n_aug << "\t" << n_local << "\t" << n_interchange << std::endl;
		std::cerr << "\t" << n_aug << "\t" << DCTB::InterchangeEventCnt << "\t" << DCTB::LocalEventCnt << "\t" << DCTB::ShiftEventCnt << "\t" << DCTB::BirthEventCnt << "\t" << DCTB::DeathEventCnt << "\t" << DCTB::AuxEventCnt  << std::endl;
	}

	/* Debug Code */
	template <class CT, class DCT>
		bool Same(const CT & ct, const DCT & dct) {
			typedef typename CT::VertexRef VertexRef;
			size_t root_cnt = 0;
			for (const auto p : ct.VertexToNode) {
				VertexRef v = ct.Parent(p.first);
				if (v.get() == p.first.get()) {
					root_cnt++;
					continue;
				}   
				bool exist = false;
				for (const auto n2 : dct.VertexToNode) {
					if (p.first.get().Coord() == n2.first.get().Coord()) {
						for (const auto n3 : dct.Neighbors(n2.second))
							if (v.get().Coord() == n3->key.get().Coord()) exist = true;
						break;
					}
				}
				if (!exist) return false;
			} 
			return root_cnt == 1;
		}

private:
	std::vector<VertexRef> sortedList;
	std::priority_queue<event> Q;

	/* vectors for debug */
	/*
	std::vector< Vertex > points;
	std::vector< std::pair<size_t, size_t> > edges;
	*/

};


#endif /* __DYNAMIC_CONTOUR_TREE_H__ */

