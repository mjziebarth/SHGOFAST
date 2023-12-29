/*
 * A C++ translation of the Vertex classes from SHGO.
 * Original files within the SHGO library:
 *    shgo/_shgo_lib/_vertex.py
 *    shgo/_shgo_lib/_complex.py
 * by Stefan Endres.
 * Translation and additions by Malte J. Ziebarth (mjz.science@fmvkb.de)
 *
 * Copyright (C) 2017 Stefan Endres,
 *               2023 Malte J. Ziebarth
 *
 * This code is licensed under the MIT license (see LICENSE).
 * SPDX-License-Identifier: MIT
 */


#include <compare>
#include <vector>
#include <memory>
#include <functional>
#include <limits>
#include <cmath>

/*
 * Default to boost::container::flat_set since iteration over the set seems
 * to be the bottle neck - at least for low dimensions.
 * Defining STD_SET_ONLY will compile using std::set instead.
 */
#ifdef STD_SET_ONLY
#include <set>
#else
#include <boost/container/flat_set.hpp>
#endif

#ifndef SHGOPATCH_VERTEX_HPP
#define SHGOPATCH_VERTEX_HPP

namespace shgopatch {

template<typename value>
class Vertex {
public:
	typedef value value_t;

	struct VertexKey {
		long hash = 0;
		std::vector<double> x;

		std::strong_ordering operator<=>(const VertexKey& other) const
		{
			std::strong_ordering res = hash <=> other.hash;
			if (res == std::strong_ordering::equal){
				std::partial_ordering resd = x <=> other.x;
				if (resd == std::partial_ordering::unordered)
					throw std::runtime_error("Could not compare two "
						                     "coordinates.");
				else if (resd == std::partial_ordering::less)
					res = std::strong_ordering::less;
				else if (resd == std::partial_ordering::greater)
					res = std::strong_ordering::greater;
			}
			return res;
		}
		
		bool operator==(const VertexKey& other) const {
			return hash == other.hash && x == other.x;
		}
	};

	Vertex() : f(std::numeric_limits<value_t>::quiet_NaN())
	{}

	Vertex(long hash, std::vector<double> x)
	   : _key({.hash=hash, .x=x}),
	     f(std::numeric_limits<value_t>::quiet_NaN())
	{
		for (double xi : _key.x)
			if (std::isinf(xi) || std::isnan(xi))
				throw std::runtime_error("Cannot handle NaN or infinite "
				                         "coordinates.");
	}
	
	~Vertex() {
		for (auto nb : nn){
			nb.v->nn.erase(*this);
		}
	}
	
	Vertex& operator=(const Vertex& other) = delete;
	
	Vertex& operator=(Vertex&& other) {
		/* Delete all links to this object: */
		for (auto nb : nn){
			nb.v->nn.erase(*this);
		}

		/* Update to values of `other`: */
		std::swap(_key, other._key);
		f = other.f;
		std::swap(nn, other.nn);
		_max = other._max;
		check_max = true;
		_min = other._min;
		check_min = true;

		/* Set all links to `other` to this one: */
		for (auto nb : nn){
			auto it = nb.v->nn.find(*this);
			if (it != nb.v->nn.end())
				it->v = this;
		}
	}
	
	const std::vector<double>& x() const {
		return _key.x;
	}

	size_t hash() const {
		return _key.hash;
	}
	
	const VertexKey& key() const {
		return _key;
	}

	void set_fval(value_t fval)
	{
		f = fval;
		check_max = true;
		check_min = true;
	}

	double get_fval() const
	{
		return f;
	}

	void connect(Vertex& other)
	{
		
		if (other == *this)
			return;
		auto it = nn.lower_bound(other);
		if (it == nn.cend() || *(it->v) != other){
			/* Vertex is not in nearest neighbors. */
			nn.emplace_hint(it, other);
			other.nn.emplace(*this);

			/* Need to recheck: */
			check_max = true;
			check_min = true;
			other.check_max = true;
			other.check_min = true;
		}
	}
	
	void disconnect(Vertex& other)
	{
		auto it = nn.find(other);
		if (it != nn.cend()){
			/* Delete the links: */
			nn.erase(it);
			other.nn.erase(*this);
			
			/* Need to recheck: */
			check_max = true;
			check_min = true;
			other.check_max = true;
			other.check_min = true;
		}
	}

	bool minimiser() const {
		if (check_min){
			_min = true;
			for (auto v : nn){
				if (f >= v.v->f)
					_min = false;
			}
			check_min = false;
		}
		return _min;
	}

	bool maximiser() const {
		if (check_max){
			_max = true;
			for (auto v : nn){
				if (f <= v.v->f)
					_max = false;
			}
			check_max = false;
		}
		return _max;
	}

	bool operator==(const Vertex& other) const {
		return _key == other._key;
	}

	std::vector<VertexKey> star() const {
		std::vector<VertexKey> s(nn.size()+1);
		auto it = s.begin();
		for (auto v : nn){
			*it = v.v->_key;
			++it;
		}
		*it = _key;
		return s;
	}

	std::vector<VertexKey> nearest_neighbors() const {
		std::vector<VertexKey> ret(nn.size());
		auto it = ret.begin();
		for (auto v : nn){
			*it = v.v->_key;
			++it;
		}
		return ret;
	}

	std::strong_ordering operator<=>(const Vertex& other) const {
		return _key <=> other._key;
	}

	/*
	 * Internals:
	 */
	struct vertex_ptr_t {
		long hash;
		Vertex* v = nullptr;

		vertex_ptr_t(Vertex& v) : hash(v.hash()), v(&v)
		{
		}

		std::strong_ordering operator<=>(const vertex_ptr_t& other) const {
			std::strong_ordering res = hash <=> other.hash;
			if (res == std::strong_ordering::equal){
				std::partial_ordering resd = v->_key.x <=> other.v->_key.x;
				if (resd == std::partial_ordering::unordered)
					throw std::runtime_error("Could not compare two "
					                         "coordinates.");
				else if (resd == std::partial_ordering::less)
					res = std::strong_ordering::less;
				else if (resd == std::partial_ordering::greater)
					res = std::strong_ordering::greater;
			}
			return res;
		}

		std::strong_ordering operator<=>(const Vertex& other) const {
			std::strong_ordering res = hash <=> other.hash();
			if (res == std::strong_ordering::equal){
				std::partial_ordering resd = v->_key.x <=> other._key.x;
				if (resd == std::partial_ordering::unordered)
					throw std::runtime_error("Could not compare two "
					                         "coordinates.");
				else if (resd == std::partial_ordering::less)
					res = std::strong_ordering::less;
				else if (resd == std::partial_ordering::greater)
					res = std::strong_ordering::greater;
			}
			return res;
		}
	};

	#ifdef STD_SET_ONLY
	typedef std::set<vertex_ptr_t> vptr_set_t;
	#else
	typedef boost::container::flat_set<vertex_ptr_t> vptr_set_t;
	#endif

	const vptr_set_t& get_nn() const {
		return nn;
	}

private:

	VertexKey _key;
	value f;
	vptr_set_t nn;
	mutable bool _max = true;
	mutable bool check_max = true;
	mutable bool _min = true;
	mutable bool check_min = true;

};


/*
 * Specializations:
 */
typedef Vertex<double> ScalarVertex;

struct compare_vertex_shared_ptr_t {
	typedef std::shared_ptr<ScalarVertex> vsp_t;

	bool operator()(const vsp_t& v0, const vsp_t& v1) const {
		if (!v0 || !v1)
			throw std::runtime_error("Nullpointer comparison.");
		return *v0 < *v1;
	}
};

/*
 * The method "proc_minimiser" is, in pure Python mode, with negligible
 * function evaluation time, and at large number of `iter`, the bottleneck
 * of the SHGO routine.
 * By implementng the graph logic in C++ and storing the relevant data
 * structures in a C++ container, `proc_minimiser_set` can be 
 */
#ifdef STD_SET_ONLY
typedef std::set<std::shared_ptr<ScalarVertex>,compare_vertex_shared_ptr_t>
        ScalarVertexPtrSet;
#else
typedef boost::container::flat_set<std::shared_ptr<ScalarVertex>,
                                   compare_vertex_shared_ptr_t>
        ScalarVertexPtrSet;
#endif

static void proc_minimiser_set(ScalarVertexPtrSet& svps)
{
	for (std::shared_ptr<ScalarVertex>& vptr : svps)
	{
		if (vptr){
			vptr->minimiser();
			vptr->maximiser();
		}
	}
}

static std::vector<ScalarVertex::VertexKey>
merge_neighborhoods(const ScalarVertex& v0, const ScalarVertex& v1)
{
	/* Build up the merged neighborhood in a set container: */
	#ifdef STD_SET_ONLY
	std::set<ScalarVertex::VertexKey> nn;
	#else
	boost::container::flat_set<ScalarVertex::VertexKey> nn;
	#endif

	/* Use v0 as a base: */
	for (const ScalarVertex::vertex_ptr_t& vptr : v0.get_nn()){
		if (vptr.v){
			nn.insert(vptr.v->key());
		}
	}

	/* Insert v1: */
	for (const ScalarVertex::vertex_ptr_t& vptr : v1.get_nn()){
		if (!vptr.v)
			continue;
		auto it = nn.lower_bound(vptr.v->key());
		if (it == nn.end() || *it != vptr.v->key())
			nn.insert(it, vptr.v->key());
	}

	/* Assemble the corresponding vector: */
	return std::vector<ScalarVertex::VertexKey>(nn.cbegin(), nn.cend());
}


}


#endif