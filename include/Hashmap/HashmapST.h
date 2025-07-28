#pragma once
#include <unordered_map>
#include "Hashmap.h"

template <typename Key, typename Value>
class HashmapST : public Hashmap<Key, Value> {
	public: 
		void clear() override {
			m.clear();
		}

		void insert(Key k, Value v) override {
			m[k] = v;
		}

		Value& get(Key k) override {
			return m[k];
		}

		Value get(Key k) const override {
			auto it = m.find(k);
			return (it != m.end()) ? it->second : Value();
		}

		bool find(Key k) override {
			const auto iter = m.find(k);
			if (iter == m.end()) return false;
			return true;
		}
	
	private:
		std::unordered_map<Key, Value> m;
};

template <typename Key, typename Value>
Hashmap<Key, Value>* createHashmap() {
	return new HashmapST<Key, Value>();
}