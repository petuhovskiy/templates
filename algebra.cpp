#include <bits/stdc++.h>

using namespace std;

namespace Algebra {

	/// Returns the amount of numbers coprime with n from range 1..n
	/// O(sqrt(N))
	template<typename T> T phi(T n) {
		T r = n;
		for (T i = 2; i * i <= n; i++) {
			if (n % i == 0) {
				r -= (r / i);
				while (n % i == 0) {
					n /= i;
				}
			}
		}
		if (n > 1) {
			r -= r / n;
		}
		return r;
	}

	/// Returns the amount of numbers coprime with n from range 1..n
	/// O(NlogN) for preprocessing, O(1) for a single call
	/// Not recommended to use with limit more than 10^6, n should be up to limit
	/// If you call this many times, keep limit non-increasing
	int numerous_phi(int limit, int n) {
		static vector<int> phi;
		static bool called = false;
		if (called) {
			return phi[n];
		} else {
			phi.resize(limit + 1);
			for (int i = limit; i; i--) {
				phi[i] = i;
			}
			for (int i = 2; i <= limit; i++) {
				if (phi[i] == i) {
					for (int j = i; j <= n; j += i) {
						phi[j] = phi[j] / i * (i - 1);
					}
				}
			}
		}
		called = true;
		return phi[n];
	}

	/// Calculates a^n 
	/// O(logN)
	template<typename U, typename V> U binpow(U a, V n) {
		U r = 1;
		while (n) {
			if (n & 1) {
				r = (r * a);
			}
			a = (a * a);
			n >>= 1;
		}
		return r;
	}

	/// Calculates (a^n) % mod
	/// O(logN)
	template<typename U, typename V> U binpow(U a, V n, U mod) {
		U r = 1;
		while (n) {
			if (n & 1) {
				r = (static_cast<long long>(r) * a) % mod;
			}
			a = (static_cast<long long>(a) * a) % mod;
			n >>= 1;
		}
		return r;
	}

	/// Calculates greatest common divisor of a and b
	/// O(log(max(a, b)))
	template<typename T> T gcd(T a, T b) {
		while (a && b) {
			if (a >= b) {
				a %= b;
			} else {
				b %= a;
			}
		}
		return a + b;
	}

	/// Returns all primes up to n 
	/// O(NlogN) for each call
	vector<int> primes_up_to(int r) {
		vector<bool> sieve;
		sieve.assign(r + 1, false);
		vector<int> res;
		for (int i = 2; i <= r; i++) {
			if (sieve[i]) {
				continue;
			}
			res.push_back(i);
			for (int j = i; j <= r; j += i) {
				sieve[j] = true;
			}
		}
		return res;
	}

	/// Returns greatest common divisor of a and b, calculates x and y
	/// such that a*x + b*y = gcd
	template<typename T> T extended_gcd(T a, T b, T &x, T&y) {
		if (a == 0) {
			x = 0;
			y = 1;
			return b;
		}

		T xx, yy;
		T d = extended_gcd(b % a, a, xx, yy);
		x = yy - (b / a) * xx;
		y = xx;
		return d;
	}

	/// Calculate n'th Fibonacci number
	/// Probably O(logN), though works incredibly fast
	/// Result stored in x
	void calculate_fibonacci(long long n, long long &x, long long &y) {
		if (n == 0) {
			x = 0;
			y = 1;
			return;
		}

		if (n & 1) {
			calculate_fibonacci(n - 1, y, x);
			y = (y + x);
		} else {
			long long a, b;
			calculate_fibonacci(n >> 1, a, b);
			y = (a * a + b * b);
			x = (a * b + a * (b - a));
		}
	}

	/// Calculate n'th Fibonacci number modulo mod
	/// Probably O(logN), though works incredibly fast
	/// Result stored in x
	void calculate_fibonacci(long long n, long long &x, long long &y, long long mod) {
		if (n == 0) {
			x = 0;
			y = 1;
			return;
		}

		if (n & 1) {
			calculate_fibonacci(n - 1, y, x, mod);
			y = (y + x) % mod;
		} else {
			long long a, b;
			calculate_fibonacci(n >> 1, a, b, mod);
			y = (a * a + b * b) % mod;
			x = (a * b + a * (b - a + mod)) % mod;
		}
	}

	/// Returns gray code
	template<typename T> T gray_code(T n) {
		return n ^ (n >> 1);
	}

	/// Finds n, such that n's gray code is gray
	template<typename T> T inverse_gray_code(T gray) {
		T n = 0;
		for (; gray; gray >>= 1) {
			n ^= gray;
		}
		return n;
	}	

	/// Find any solution for ax + by = c
	template<typename T> bool diofantite2_any_solution(T a, T b, T c, T &x, T &y, T &g) {
		g = extended_gcd(abs(a), abs(b), x, y);
		if (c % g != 0) {
			return false;
		}
		x *= c / g;
		y *= c / g;
		if (a < 0) {
			x *= -1;
		}
		if (b < 0) {
			y *= -1;
		}
		return true;
	}

	/// Returns least prime divisor of n
	/// O(NlogN) for preprocessing, O(1) for a single call
	/// Not recommended to use with limit more than 10^6, n should be up to limit
	/// If you call this many times, keep limit non-increasing
	int least_prime_divisor(int limit, int n) {
		static vector<int> lp;
		static bool called = false;
		if (called) {
			return lp[n];
		} else {
			vector<int> pr;
			lp.resize(limit + 1);
			for (int i = 2; i <= limit; i++) {
				if (lp[i] == 0) {
					lp[i] = i;
					pr.push_back(i);
				}
				for (int j = 0; j < (int)pr.size() && pr[j] <= lp[i] && i * pr[j] <= limit; j++) {
					lp[i * pr[j]] = pr[j];
				}
			}
		}
		called = true;
		return lp[n];
	}
	
};

int main() {

	// Testing & usage example of phi
	assert(Algebra::phi(1) == 1);
	assert(Algebra::phi(2) == 1);
	assert(Algebra::phi(3) == 2);
	assert(Algebra::phi(4) == 2);
	assert(Algebra::phi(5) == 4);

	// Testing & usage example of numerous_phi
	assert(Algebra::numerous_phi(5, 5) == 4);
	assert(Algebra::numerous_phi(5, 4) == 2);
	assert(Algebra::numerous_phi(5, 3) == 2);
	assert(Algebra::numerous_phi(5, 2) == 1);
	assert(Algebra::numerous_phi(5, 1) == 1);

	// Testing & usage example of binpow
	assert(Algebra::binpow(2, 5) == 32);
	assert(Algebra::binpow(3, 11) == 177147);
	assert(Algebra::binpow(5, 7) == 78125);

	// Testing & usage example of binpow with modulo
	assert(Algebra::binpow(2, 5, 31) == 1);
	assert(Algebra::binpow(3, 11, 177147) == 0);
	assert(Algebra::binpow(2, 12344, 2) == 0);

	// Testing & usage example of gcd
	assert(Algebra::gcd(2, 6) == 2);
	assert(Algebra::gcd(33, 11) == 11);
	assert(Algebra::gcd(354, 2) == 2);
	assert(Algebra::gcd(0, 5) == 5);

	// Testing & usage example of primes_up_to
	vector<int> a = Algebra::primes_up_to(10);
	assert(a == ((vector<int>){2, 3, 5, 7}));
	
	// Testing & usage example of extended_gcd
	int x, y;
	assert(Algebra::extended_gcd(12, 30, x, y) == 6);
	assert(x * 12 + y * 30 == 6);

	// Testing & usage example of calculate_fibonacci
	long long u, v;
	Algebra::calculate_fibonacci(5, u, v);
	assert(u == 5);

	// Testing & usage example of calculate_fibonacci with modulo
	Algebra::calculate_fibonacci(5, u, v, 3);
	assert(u == 2);

	// Testing & usage example of gray_code
	assert(Algebra::gray_code(5) == 7);
	assert(Algebra::gray_code(1) == 1);

	// Testing & usage example of inverse_gray_code
	assert(Algebra::inverse_gray_code(7) == 5);
	assert(Algebra::inverse_gray_code(1) == 1);

	// Testing & usage example of diofantite2_any_solution
	long long q, b, c, d, e, f;
	assert(Algebra::diofantite2_any_solution(q, b, c, d, e, f));
	assert(d * q + e * b == c);

	// Testing & usage example of least_prime_divisor
	assert(Algebra::least_prime_divisor(10, 2) == 2);
	assert(Algebra::least_prime_divisor(10, 6) == 2);
	assert(Algebra::least_prime_divisor(10, 7) == 7);

	return 0;
}
