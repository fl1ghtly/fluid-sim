#pragma once

template <typename Key, typename Value>
class Hashmap {
	public:
		virtual void clear() = 0;
		virtual void insert(Key k, Value v) = 0;
		virtual Value& get(Key k) = 0;
		virtual Value get(Key k) const = 0;
		virtual bool find(Key k) = 0;
};