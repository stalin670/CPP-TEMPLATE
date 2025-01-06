// x.cpp

#include<bits/stdc++.h>
using namespace std;

#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>

using namespace __gnu_pbds;

#define MOD 1000000007
#define MOD1 998244353
#define INF 1000000000000000000LL
#define endl "\n"
#define nline cout << endl
#define pb push_back
#define ppb pop_back
#define m_p make_pair
#define ff first
#define ss second
#define bg begin
#define lbd lower_bound
#define ubd upper_bound
#define pll pair<ll, ll>
#define PI 3.141592653589793238462
#define set_bits(x) __builtin_popcountll(x)
#define sz(x) ((int)(x).size())
#define all(x) (x).begin(), (x).end()
#define allR(x) (x).rbegin(), (x).rend()

// -------------------------------------------------------------------------------------------------------//

template <typename T> void read(T& t) { cin >> t; }
template <typename T, typename... Args> void read(T& t, Args&... args) { read(t); read(args...); }
template <typename T> void read(vector<T>& vec) { for (auto& element : vec) { cin >> element; } }
template <typename T> void print(const T& t) { cout << t << " "; }
template <typename T> void prints(const T& t) { cout << t; }
template <typename T, typename... Args> void print(const T& t, const Args&... args) { print(t); print(args...); }
template <typename T, typename... Args> void prints(const T& t, const Args&... args) { prints(t); print(args...); }
template <typename T> void print(const vector<T>& vec) { for (const auto& element : vec) { cout << element << " "; } cout << "\n"; }
template<class T> using oset = tree<T, null_type, less_equal<T>, rb_tree_tag, tree_order_statistics_node_update>;

// -------------------------------------------------------------------------------------------------------//

using ll = long long;
using ull = unsigned long long;
using lld = long double;

void solve() {
	ll n, m;
	read(n, m);
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif
	ll tc;
	tc = 1;
	// read(tc);
	for (ll i = 1; i <= tc; i++) {
		solve();
	}

	return 0;
}
