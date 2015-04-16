#ifndef __UNION_FIND_H__
#define __UNION_FIND_H__

#include <map>
#include <cassert>

template <class T, typename Hash=std::hash<T>, typename Equal=std::equal_to<T> >
class UnionFind {
	private:
	typedef typename std::unordered_map<T, T, Hash, Equal> Map;
	Map uf;

	public:
	T Root(T i) {
		// Path halving
		while (!Equal()(uf.at(i),i)) {
			T p = uf.at(i);
			uf.at(i) = uf.at(p);
			i = p;
		}
		return i;
	}
	void Add(T i) {
		if (uf.count(i) == 0) 
			uf.insert(std::make_pair(i,i));
	}
	void Empty() {
		return uf.empty();
	}
	void Join(T f, T t) {
		uf.at(Root(f)) = t;
	}
};

#endif /* __UNION_FIND_H__ */
