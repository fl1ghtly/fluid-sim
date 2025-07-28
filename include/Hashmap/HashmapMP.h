#pragma once
#include <oneapi/tbb/concurrent_hash_map.h>
#include "Hashmap.h"

template <typename Key, typename Value>
class HashmapMP : public Hashmap<Key, Value> {
	public: 
		void clear() override {
			m.clear();
		};

		void insert(Key k, Value v) override {
			typename tbb::concurrent_hash_map<Key, Value>::accessor a;
			m.insert(a, {k, v});
		};

		Value& get(Key k) override {
			static Value empty;
			typename tbb::concurrent_hash_map<Key, Value>::accessor a;
			if (m.find(a, k)) return a->second;
			return empty;
		}
		
		Value get(Key k) const override {
			typename tbb::concurrent_hash_map<Key, Value>::const_accessor ca;
			if (m.find(ca, k)) return ca->second;
			return Value();
		}

		bool find(Key k) override {
			typename tbb::concurrent_hash_map<Key, Value>::const_accessor ca;
			return m.find(ca, k);
		};
	private:
		tbb::concurrent_hash_map<Key, Value> m;
};

template <typename Key, typename Value>
Hashmap<Key, Value>* createHashmap() {
	return new HashmapMP<Key, Value>();
}
