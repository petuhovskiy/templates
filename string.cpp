#include <bits/stdc++.h>

using namespace std;

namespace String {

	vector<int> z_function(string s) {
		int n = (int)s.size();
		vector<int> z(n);
		for (int i = 1, l = 0, r = 0; i < n; i++) {
			if (i <= r) {
				z[i] = min(r - i + 1, z[i - l]);
			}
			while (i + z[i] < n && s[z[i]] == s[i + z[i]]) {
				z[i]++;
			}
			if (i + z[i] - 1 > r) {
				l = i;
				r = i + z[i] - 1;
			}
		}
		return z;
	}

	vector<int> prefix_function(string s) {
		int n = (int)s.size();
		vector<int> pi(n);
		for (int i = 1; i < n; i++) {
			int j = pi[i - 1];
			while (j > 0 && s[i] != s[j]) {
				j = pi[j - 1];
			}
			if (s[i] == s[j]) {
				j++;
			}
			pi[i] = j;
		}
		return pi;
	}

	pair<vector<int>, vector<int> > manacher(string v) {
		int n = (int)v.size();
		vector<int> p1(n), p2(n);
		int l = 0, r = -1;
		for (int i = 0; i < n; i++) {
		    int k = (i > r ? 1 : min(p1[l + r - i], r - i + 1));
		    while (i + k < n && i - k >= 0 && v[i + k] == v[i - k]) k++;
		    p1[i] = k--;
		    if (i + k > r) {
		        l = i - k;
		        r = i + k;
		    }
		}
		l = 0;
		r = -1;
		for (int i = 0; i < n; i++) {
		    int k = (i >= r ? 0 : min(p2[l + r - i - 1], r - i));
		    while (i - k >= 0 && i + k + 1 < n && v[i - k] == v[i + k + 1]) k++;
		    p2[i] = k--;
		    if (i + k + 1 > r) {
		        l = i - k;
		        r = i + k + 1;
		    }
		}
		return make_pair(p1, p2);
	}

	string min_cyclic_shift(string s) {
		s += s;
		int n = (int)s.size();
		int i = 0, ans = 0;
		while (i < n / 2) {
			ans = i;
			int j = i + 1, k = i;
			while (j < n && s[k] <= s[j]) {
				if (s[k] < s[j]) {
					k = i;
				} else {
					k++;
				}
				j++;
			}
			while (i <= k) {
				i += j - k;
			}
		}
		return s.substr(ans, n / 2);
	}

	vector<int> suffix_array(string s) {
		s += '#';
		int n = (int)s.size();
		vector<int> p(n), cnt(300), c(n);
		fill(cnt.begin(), cnt.end(), 0);
		for (int i = 0; i < n; i++) {
			cnt[s[i]]++;
		}
		for (int i = 1; i < 300; i++) {
			cnt[i] += cnt[i - 1];
		}
		for (int i = 0; i < n; i++) {
			p[--cnt[s[i]]] = i;
		}
		c[p[0]] = 0;
		int classes = 1;
		for (int i = 1; i < n; i++) {
			if (s[p[i]] != s[p[i - 1]]) {
				classes++;
			}
			c[p[i]] = classes - 1;
		}
		vector<int> pn(n), cn(n);
		for (int h = 0; (1 << h) < n; h++) {
			for (int i = 0; i < n; i++) {
				pn[i] = p[i] - (1 << h);
				if (pn[i] < 0) {
					pn[i] += n;
				}
			}
			fill(cnt.begin(), cnt.end(), 0);
			for (int i = 0; i < n; i++) {
				cnt[c[pn[i]]]++;
			}
			for (int i = 1; i < classes; i++) {
				cnt[i] += cnt[i - 1];
			}
			for (int i = n - 1; i >= 0; i--) {
				p[--cnt[c[pn[i]]]] = pn[i];
			}
			cn[p[0]] = 0;
			classes = 1;
			for (int i = 1; i < n; i++) {
				int mid1 = (p[i] + (1 << h)) % n;
				int mid2 = (p[i - 1] + (1 << h)) % n;
				if (c[p[i]] != c[p[i - 1]] || c[mid1] != c[mid2]) {
					classes++;
				}
				cn[p[i]] = classes - 1;
			}
			for (int i = 0; i < n; i++) {
				c[i] = cn[i];
			}
		}
		vector<pair<int, int> > sa;
		for (int i = 0; i < n; i++) {
			sa.push_back(make_pair(c[i], i));
		}
		sort(sa.begin(), sa.end());
		vector<int> id;
		for (int i = 0; i < n; i++) {
			id.push_back(sa[i].second);

		}
		return id;
	}

};

int main() {


	return 0;
}