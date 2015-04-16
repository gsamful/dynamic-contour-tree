#ifndef __VERTEX_H__
#define __VERTEX_H__

#include <vector>
#include <list>
#include <cassert>
#include <numeric>
#include <algorithm>


enum VertexType {
	REGULAR = 0,
	SADDLE = 1,
	MINIMUM = 2,
	MAXIMUM = 3,
	MONKEY = 4,
};

template <class Vertex>
struct VertexPair {
	std::reference_wrapper<Vertex> u;
	std::reference_wrapper<Vertex> v;
	VertexPair (Vertex &u, Vertex &v) : u(u), v(v) {}
};

template <class Vertex>
struct VertexRefEqual{
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	bool operator() (const VertexRef &v1, const VertexRef &v2) const {
		return v1.get() == v2.get();
	}
};

template <class Vertex>
struct VertexRefLess{
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	bool operator() (const VertexRef &v1, const VertexRef &v2) const {
		if (v1.get().Height() == v2.get().Height()) return v1.get() < v2.get();
		return v1.get().Height() < v2.get().Height();
	}
};


template <class T>
class Vertex {
	public:
	typedef T HeightType;
	typedef typename std::reference_wrapper<Vertex> VertexRef;
	typedef std::vector< VertexRef > VertexList;

	virtual const VertexList Neighbors() const = 0;
	virtual const VertexType Type() const = 0;
	virtual void Update(const HeightType & h) = 0;
	virtual HeightType Height() = 0; 
};


template <class HT>
class NaiveVertex /*: public Vertex<HT> */ {
	public:
	typedef HT HeightType;
	typedef NaiveVertex<HT> Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef typename std::list<VertexRef> VertexList;
	typedef VertexRefLess<Vertex> LessVertexRef;

	private: 
	HeightType height;
	VertexType type;

	VertexList neighbors;

	enum Dir {
		UP = -1,
		DOWN = 1,
		UNKNOWN = 0,
	};

	public:
	// This assumes that there is no vertices with the same height
	struct Hash {
		size_t operator()(const Vertex &v) const {
			return std::hash<size_t>()(v.Height());
		}
	};

	NaiveVertex(const HeightType &h) : 
		height(h), type(REGULAR) {}

	void TypeCheck() {
		bool prev_dir = LessVertexRef()(*this, neighbors.back());
		int cnt = 0;
		for (const Vertex& neighbor : neighbors) {
			bool dir = LessVertexRef()(*this, neighbor);
			if (dir != prev_dir) cnt++;
			prev_dir = dir;
		}
		if (cnt > 3) type = SADDLE;
		else if (cnt > 1) type = REGULAR;
		else if (prev_dir) type = MINIMUM;
		else type = MAXIMUM;
	}

	void ChangeHeight(const HeightType &h) {
		height = h;
	}

	void ChangeType(const VertexType t) {
		type = t;
	}

	bool IsConnected(const Vertex &v) {
		for (auto neighbor : neighbors)
			if (neighbor.get() == v) return true;
		return false;
	}

	void Connect(const Vertex &v) {
		if (!IsConnected(v)) neighbors.push_back(v);
	}
	
	void ConnectBack(const Vertex &v) {
		if (!IsConnected(v)) neighbors.push_back(v);
	}

	void ConnectFront(const Vertex &v) {
		if (!IsConnected(v)) neighbors.push_front(v);
	}
	
	VertexRef Successor(VertexRef v) const {
		for (auto neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor)
			if (neighbor->get() == v.get()) {
				++neighbor;
				if (neighbor != neighbors.end()) return *neighbor;
				return *(neighbors.begin());
			}
		assert(false);
	}

	VertexRef Predecessor(VertexRef v) const {
		for (auto neighbor = neighbors.rbegin(); neighbor != neighbors.rend(); ++neighbor)
			if (neighbor->get() == v.get()) {
				++neighbor;
				if (neighbor != neighbors.rend()) return *neighbor;
				return *(neighbors.rbegin());
			}
		assert(false);
	}

	const VertexList & Neighbors() const {
		return neighbors;
	}

	const VertexType Type() const {
		return type;
	}

	HeightType Height() const {
		return height;
	}

	// This assumes there is no vertices with same height
	bool operator == (const Vertex &v) const {
		return Height() == v.Height();
	}
	bool operator != (const Vertex &v) const {
		return Height() != v.Height();
	}
	bool operator < (const Vertex &v) const {
		return height < v.height;
	}
	bool operator > (const Vertex &v) const {
		return height > v.height;
	}
};

// This class is copy-constructable only when neighbor is not formed yet..
template <class CT, class HT = CT>
class Vertex2D {
	public:
	typedef HT HeightType;
	typedef CT CoordType;
	typedef Vertex2D<CT, HT> Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef typename std::vector<VertexRef> VertexList;
	typedef VertexRefLess<Vertex> LessVertexRef;

	protected:
	CoordType x;
	CoordType y;
	HeightType height;
	VertexType type;

	VertexList neighbors;

	enum Dir {
		UP = -1,
		DOWN = 1,
		UNKNOWN = 0,
	};

	struct CircleOrder {
		const Vertex & o;

		CircleOrder(const Vertex & v) : o(v) {}

		bool operator()(const Vertex& a, const Vertex& b) const {
			if (a.x - o.x >= 0 && b.x - o.x < 0)
				return true;
			if (a.x - o.x < 0 && b.x - o.x >= 0)
				return false;
			if (a.x - o.x == 0 && b.x - o.x == 0) {
				if (a.y - o.y >= 0 || b.y - o.y >= 0)
					return a.y > b.y;
				return b.y > a.y;
			}

			CT det = (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
			if (det < 0)
				return true;
			if (det > 0)
				return false;

			// points a and b are on the same line from the center
			// check which point is closer to the center
			CT d1 = (a.x - o.x) * (a.x - o.x) + (a.y - o.y) * (a.y - o.y);
			CT d2 = (b.x - o.x) * (b.x - o.x) + (b.y - o.y) * (b.y - o.y);
			return d1 > d2;
		}
	};

	public:
	struct Hash {
		size_t operator()(const Vertex &v) const {
			std::pair<CT, CT> co = v.Coord();
			return std::hash<CT>()(co.first) ^ std::hash<CT>()(co.second);
		}
	};

	void TypeCheck() {
		bool prev_dir = LessVertexRef()(*this, neighbors.back());
		int cnt = 0;
		for (const Vertex& neighbor : neighbors) {
			bool dir = LessVertexRef()(*this, neighbor);
			if (dir != prev_dir) cnt++;
			prev_dir = dir;
		}
		if (cnt > 3) type = SADDLE;
		else if (cnt > 1) type = REGULAR;
		else if (prev_dir) type = MINIMUM;
		else type = MAXIMUM;
	}

	Vertex2D(const CoordType &x, const CoordType &y, const HeightType &h) : 
		x(x), y(y), height(h), type(REGULAR) {}

	void ChangeHeight(const HeightType &h) {
		height = h;
	}

	void ChangeType(const VertexType t) {
		type = t;
	}

	bool IsConnected(const Vertex &v) const {
		for (auto neighbor : neighbors)
			if (neighbor.get() == v) return true;
		return false;
	}

	void Connect(const Vertex &v) {
		if (!IsConnected(v)) neighbors.push_back(v);
	}

	// This assume that all edges are sorted in Circle order
	void ConnectBoundary(const Vertex &inf, const Vertex &prev, const Vertex &next) {
		if (!IsConnected(inf)) {
			for (size_t i = 0 ; i < neighbors.size() ; ++i) {
				if ((neighbors[i].get() == prev && neighbors[(i+1)%neighbors.size()].get() == next) ||
					(neighbors[i].get() == next && neighbors[(i+1)%neighbors.size()].get() == prev)) {
					std::rotate(neighbors.begin(), neighbors.begin() + i + 1, neighbors.end());
					Connect(inf);
					break;
				}
			}
		}
	}

	void SortNeighbor() {
		std::sort(neighbors.begin(), neighbors.end(), CircleOrder(*this));
	}

	VertexRef Successor(VertexRef v) const {
		for (auto neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor)
			if (neighbor->get() == v.get()) {
				++neighbor;
				if (neighbor != neighbors.end()) return *neighbor;
				return *(neighbors.begin());
			}
		assert(false);
	}

	VertexRef Predecessor(VertexRef v) const {
		for (auto neighbor = neighbors.rbegin(); neighbor != neighbors.rend(); ++neighbor)
			if (neighbor->get() == v.get()) {
				++neighbor;
				if (neighbor != neighbors.rend()) return *neighbor;
				return *(neighbors.rbegin());
			}
		assert(false);
	}

	bool IsNeighbor(VertexRef v) const {
		for (auto neighbor : neighbors) 
			if (neighbor.get() == v.get()) return true;
		return false;
	}

	const VertexList & Neighbors() const {
		return neighbors;
	}

	const VertexType Type() const {
		return type;
	}

	HeightType Height() const {
		return height;
	}

	std::pair<CoordType, CoordType> Coord() const {
		return std::make_pair(x, y);
	}

	// This assumes there is no vertices with same height
	bool operator == (const Vertex &v) const {
		return Coord() == v.Coord();
	}

	bool operator != (const Vertex &v) const {
		return !(*this == v);
	}

	bool operator < (const Vertex &v) const {
		return Coord() < v.Coord();
	}

	bool operator > (const Vertex &v) const {
		return Coord() > v.Coord();
	}
	
	struct PositionComp{
		bool operator ()(const Vertex &v1, const Vertex &v2) const {
			return v1.Coord() < v2.Coord();
		}
	};
	
	static Vertex2D InfVertex(){
		return Vertex2D(std::numeric_limits<CT>::lowest(), std::numeric_limits<CT>::lowest(), std::numeric_limits<HT>::lowest());
	}
};

template <class T>
struct LinearFunction{
	typedef T InstType;
	typedef T AugType;

	T a, b;
	LinearFunction(T a, T b): a(a), b(b) {}

	InstType operator()(AugType t = AugType()) const {
		return a * t + b;
	}

	LinearFunction operator - (const LinearFunction & f) const {
		return LinearFunction(a - f.a, b - f.b);
	}

	LinearFunction operator + (const LinearFunction & f) const {
		return LinearFunction(a + f.a, b + f.b);
	}

	AugType solve() const {
		if (a == 0) return std::numeric_limits<T>::lowest();
		return -b/(AugType)a;
	}

	bool operator < (const LinearFunction & f) const {
		if (a == f.a) return b < f.b;
		return a < f.a;
	}

	bool operator > (const LinearFunction & f) const {
		if (a == f.a) return b > f.b;
		return a > f.a;
	}

	bool operator == (const LinearFunction & f) const {
		return this->a == f.a && this->b == f.b;
	}
	bool operator != (const LinearFunction & f) const {
		return !(*this == f);
	}

	static LinearFunction lowest() {
		return LinearFunction(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest());
	}
};

namespace std {
template <class T>
struct numeric_limits< LinearFunction<T> >  {
	static LinearFunction<T> lowest() { return LinearFunction<T>::lowest(); }
};
};

template <class T>
std::ostream & operator<< (std::ostream &os, const LinearFunction<T> &f) {
	os << f.a << "t + " << f.b;
	return os;
}

template <class CT, class T>
struct VertexRefLess<Vertex2D<CT, LinearFunction<T> > >{
	typedef Vertex2D<CT, LinearFunction<T> > Vertex;
	typedef typename std::reference_wrapper<const Vertex> VertexRef;
	typedef typename LinearFunction<T>::AugType HeightAugType;
	typedef typename LinearFunction<T>::InstType HeightInstType;
	HeightAugType t;
	VertexRefLess() : t(HeightAugType()) {}
	bool operator() (const VertexRef &v1, const VertexRef &v2) const {
		if (v2.get().Height() == LinearFunction<T>::lowest()) return false;
		if (v1.get().Height() == LinearFunction<T>::lowest()) return true;
		HeightInstType v1h = v1.get().Height()(t), v2h = v2.get().Height()(t);
		if (v1h == v2h) return v1.get() < v2.get();
		return v1h < v2h;
	}
};

template <class Vertex>
std::ostream & operator<< (std::ostream &os, const VertexPair<Vertex> &v) {
	os << (Vertex &)v.u << " - " << (Vertex &)v.v;
	return os;
}

template <class HT>
std::ostream &operator<<(std::ostream &os, const NaiveVertex<HT> &v) {
	return os << v.Height();
}

template <class CT, class HT>
std::ostream &operator<<(std::ostream &os, const Vertex2D<CT, HT> &v) {
	return os << "(" << v.Coord().first << ", " << v.Coord().second << " [" << v.Height() << "])";
}


#endif /* __VERTEX_H__ */
