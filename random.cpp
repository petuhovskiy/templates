#include <bits/stdc++.h>

using namespace std;


namespace Random {

	long long rdtsc() {
		long long tmp;
		asm("rdtsc": "=A"(tmp));
		return tmp;
	}
	
	void initialize() {
		srand(rdtsc());
	}

	inline long long mr() {
		long long r = (rand() << 15LL) ^ rand();
		return r;
	}

	long long random(long long limit) {
		return mr() % limit;
	}

};

int main() {
	

	return 0;
}