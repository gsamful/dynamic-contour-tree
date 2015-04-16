#include <iostream>
#include <vector>
#include <ctime>
#include <sys/time.h>
#include <cstdlib>
#include <cmath>
#include "terrain.h"
#include "vertex.h"
#include "contour_tree.h"
#include "dynamic_contour_tree.h"
#include "dynamic_contour_tree_interfaces.h"

using namespace std;

template <class CT, class DCT>
bool Same(CT & ct, DCT & dct) {
	typedef typename CT::VertexRef VertexRef;
	size_t root_cnt = 0;
	for (auto p : ct.VertexToNode) {
		VertexRef v = ct.Parent(p.first);
		if (v.get() == p.first.get()) {
			root_cnt++;
			continue;
		}
		bool exist = false;
		for (auto n2 : dct.Neighbors(p.first)) {
			if (v.get() == n2->key.get()) exist = true;
		}
		if (!exist) return false;
	}
	return root_cnt == 1;
}

typedef long double HeightType;
typedef HeightType CoordType;
typedef LinearFunction<HeightType> HT;
typedef Vertex2D<CoordType, HT> TV;
typedef EdgeTerrain<TV> TT;
typedef typename std::reference_wrapper<const TV> TimeVertexRef;

std::vector< pair<size_t, size_t> > set_grid(size_t sqrn) {
	size_t n = sqrn * sqrn;
	vector< pair<size_t, size_t> > edges(0);
	for (size_t i = 0 ; i < n; i++) {
		size_t t;
		size_t r = i/sqrn, c = i%sqrn;
		if (i % 2 == 0) { // diagonal
			t = i - (sqrn+1);
			if (t < n && t/sqrn+1 == r) edges.push_back(make_pair(i, t));
			t = i - (sqrn-1);
			if (t < n && t/sqrn != r) edges.push_back(make_pair(i, t));
			t = i + (sqrn-1);
			if (t < n && t/sqrn != r) edges.push_back(make_pair(i, t));
			t = i + (sqrn+1);
			if (t < n && t/sqrn-1 == r) edges.push_back(make_pair(i, t));
		}
		t = i - sqrn;
		if (t < n) edges.push_back(make_pair(i, t));
		t = i - 1;
		if (t < n && t/sqrn == r) edges.push_back(make_pair(i, t));
		t = i + 1;
		if (t < n && t/sqrn == r) edges.push_back(make_pair(i, t));
		t = i + sqrn;
		if (t < n) edges.push_back(make_pair(i, t));
	}
	return edges;
}

// 0 is the center
std::vector< pair<size_t, size_t> > set_ring(size_t n) {
	vector< pair<size_t, size_t> > edges(0);
	edges.push_back(make_pair(0,1));
	edges.push_back(make_pair(1,n-1));
	for (size_t i = 2; i < n ; ++i) {
		edges.push_back(make_pair(0, i));
		edges.push_back(make_pair(i, i-1));
	}
	return edges;
}

template <typename T>
void ExecuteAndVarify(const std::vector<T> _points, const std::vector<pair<size_t, size_t > > edges) {
	std::vector<T> points = _points;
	for (auto &point : points) {
		HT ct = point.Height();
		point.ChangeHeight(HT(ct.b, ct.a));
	}
	TT bt(points, edges);
	MergeTree<TT> mt(bt);
	SplitTree<TT> st(bt);
	ContourTree< LinkCutTree<TimeVertexRef, Statistic, true>, TT > ct(mt,st);
	for (auto &point : points) {
		HT ct = point.Height();
		point.ChangeHeight(HT(ct.b, ct.a));
	}
	TT t(points, edges);
	TimeVaryingDynamicContourTree < STContourTree<TimeVertexRef>, LinkCutTree<TimeVertexRef, Statistic, true>, TT > dct(t);
	dct.run();
	Same(ct, dct);

}

template <typename T>
void GenerateTopDownLinearFunction(std::vector<T> & a, std::vector<T> & b, size_t n) {
	a.resize(n), b.resize(n);
	for (size_t i = 0 ; i < n ; ++i) a[i] = (rand()%(n*n*81)), b[i] = rand()%(n*37);
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());
	reverse(b.begin(), b.end());
	for (size_t i = 1 ; i < n ; ++i) {
		while (b[i-1] <= b[i])	b[i]-=1;
	}
	while (a[n-1] <= a[n-2]) a[n-2] = (rand()%(n*n*81));
	for (size_t i = n-3 ; i < n ; --i) {
		// To reduce the precision error, only compare the number with one row above
		// This function works only for sqrt(n) * sqrt(n) grid
		// For general case, target should be n-1 always
		size_t target = min(n-1, (size_t)(i+2 * sqrt(n))+1);
		a[i] = a[target] - (b[i] - b[target]) * (a[i+2] - a[i+1]) / (long double) (b[i+1] - b[i+2]) - 1;
	}
}
 
template <typename T>
void GenerateBottomUpLinearFunction(std::vector<T> & a, std::vector<T> & b, size_t n) {
	a.resize(n), b.resize(n);
	for (size_t i = 0 ; i < n ; ++i) a[i] = (rand()%(n*n*81)), b[i] = rand()%(n*37);
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());
	reverse(b.begin(), b.end());
	for (size_t i = 1 ; i < n ; ++i) {
		while (b[i-1] <= b[i])	b[i]-=1;
	}
	for (size_t i = 2 ; i < n ; ++i) {
		a[i] = a[0] - (b[i]-b[0])*(a[i-2]-a[i-1])/(long double) (b[i-1] - b[i-2]) + 1;
	}
}

template <typename T>
void GenerateNeighborChangeLinearFunction(std::vector<T> & a, std::vector<T> & b, size_t n) {
	a.resize(n), b.resize(n);
	for (size_t i = 0 ; i < n ; ++i) a[i] = (rand()%(n*n*81)), b[i] = rand()%(n*37);
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());
	reverse(b.begin(), b.end());
	for (size_t i = 1 ; i < n ; ++i) {
		while (b[i-1] <= b[i])	b[i]-=1;
	}
	for (size_t i = 2 ; i < n ; ++i) {
		a[i] = (a[i-1]-a[i-2]) * (b[i-1]-b[i]) / (long double)(b[i-2]-b[i-1]) + a[i-1] + 1;
	}
}
 
void TimeVaryingGrid(size_t sqrn) {
	size_t n = sqrn*sqrn;
	vector<HeightType> a(n, 0), b(n, 0);	
	vector< TV > points;
	vector< pair<size_t, size_t> > edges;
	GenerateTopDownLinearFunction(a, b, n);
	for (size_t i = 0 ; i < n ; ++i) {
		size_t r = i/sqrn, c = i%sqrn;
		HT h(0,0);
		if (!((i/sqrn)%2)) h = HT(a[i],b[i]);
		else h = HT(a[r*sqrn + (sqrn-c-1)], b[r*sqrn + (sqrn-c-1)]);
		points.push_back(TV(r, c, h));
	}
	std::cerr << n;
	edges = set_grid(sqrn);
	ExecuteAndVarify(points, edges);
}

void TimeVaryingGridRandom(size_t sqrn) {
	size_t n = sqrn*sqrn;
	vector<HeightType> a(n, 0), b(n, 0);	
	vector< TV > points;
	vector< pair<size_t, size_t> > edges;
	assert(points.empty() && edges.empty());
	for (size_t i = 0 ; i < n ; ++i) {
		size_t r = i/sqrn, c = i%sqrn;
		a[i] = (rand()%(n*n*81)), b[i] = rand()%(n*n*37);
		HT h(a[i],b[i]);
		points.push_back(TV(r, c, h));
	}
	edges = set_grid(sqrn);
	std::cerr << n;
	ExecuteAndVarify(points, edges);
}

void TimeVaryingCircle(size_t n, double r = 1000) {
	CoordType x, y;
	double theta = 2 * acos(-1) / (n-1);
	vector<HeightType> a(n, 0), b(n, 0);	
	vector< TV > points;
	vector< pair<size_t, size_t> > edges;
	x = 0, y = 0;
	GenerateNeighborChangeLinearFunction(a, b, n);
	a.resize(n), b.resize(n);
	a[n-2] = (b[0]-b[n-2]) * (a[n-3] - a[n-4]) / (long double)(b[n-4]-b[n-3]) + a[0] + (rand()%7+1);
	a[n-1] = (b[0]-b[n-1]) * (a[n-2] - a[n-3]) / (long double)(b[n-3]-b[n-2]) + a[0] + (rand()%7+1);

	points.push_back(TV(x, y, HT(a[n-2], b[n-2])));
	points.push_back(TV(x, y = r, HT(a[n-1], b[n-1])));
	for (size_t i = n-3 ; i < n ; --i) {
		CoordType nx = cos(theta)*x - sin(theta)*y;
		CoordType ny = sin(theta)*x + cos(theta)*y;
		points.push_back(TV(nx, ny, HT(a[i],b[i])));
		x = nx, y = ny;
	}
	edges = set_ring(n);
	std::cerr << n;
	ExecuteAndVarify(points, edges);
}


int main() {
	srand(time(0));
//	for (size_t i = 3 ; i <= 60 ; i+=2) TimeVaryingGrid(i);
//	for (size_t i = 3 ; i <= 60 ; i+=2) TimeVaryingCircle(i*i);
	for (size_t i = 3 ; i <= 60 ; i+=2) TimeVaryingGridRandom(i);
}
