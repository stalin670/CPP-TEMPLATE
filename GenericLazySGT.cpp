// F.cpp

#include<bits/stdc++.h>
using namespace std;

#define MOD 1000000007
#define MOD1 998244353
#define INF 1e18
#define endl "\n"
#define pb push_back
#define ppb pop_back
#define mp make_pair
#define ff first
#define ss second
#define PI 3.141592653589793238462
#define set_bits(x) __builtin_popcountll(x)
#define sz(x) ((int)(x).size())
#define all(x) (x).begin(), (x).end()

using ll = long long;
using ull = unsigned long long;
using lld = long double;

template<typename Node, typename Update>
struct LazySGT {
	vector<Node> tree;
	vector<bool> lazy;
	vector<Update> updates;
	vector<ll> arr;
	ll n;
	ll s;

	LazySGT(ll a_len, vector<ll> &a) {
		arr = a;
		n = a_len;
		s = 1;
		while (s < 2 * n) {
			s = s << 1;
		}
		tree.resize(s); fill(all(tree), Node());
		lazy.resize(s); fill(all(lazy), false);
		updates.resize(s); fill(all(updates), Update());
		build(0, n - 1, 1);
	}

	void build(ll start, ll end, ll index) {
		if (start == end)   {
			tree[index] = Node(arr[start]);
			return;
		}
		ll mid = (start + end) / 2;
		build(start, mid, 2 * index);
		build(mid + 1, end, 2 * index + 1);
		tree[index].merge(tree[2 * index], tree[2 * index + 1]);
	}

	void pushdown(ll index, ll start, ll end) {
		if (lazy[index]) {
			ll mid = (start + end) / 2;
			apply(2 * index, start, mid, updates[index]);
			apply(2 * index + 1, mid + 1, end, updates[index]);
			updates[index] = Update();
			lazy[index] = 0;
		}
	}

	void apply(ll index, ll start, ll end, Update& u) {
		if (start != end) {
			lazy[index] = 1;
			updates[index].combine(u, start, end);
		}
		u.apply(tree[index], start, end);
	}

	void update(ll start, ll end, ll index, ll left, ll right, Update& u) {  // Never Change this
		if (start > right || end < left)
			return;
		if (start >= left && end <= right) {
			apply(index, start, end, u);
			return;
		}
		pushdown(index, start, end);
		ll mid = (start + end) / 2;
		update(start, mid, 2 * index, left, right, u);
		update(mid + 1, end, 2 * index + 1, left, right, u);
		tree[index].merge(tree[2 * index], tree[2 * index + 1]);
	}

	Node query(ll start, ll end, ll index, ll left, ll right) {
		if (start > right || end < left)
			return Node();
		if (start >= left && end <= right) {
			pushdown(index, start, end);
			return tree[index];
		}
		pushdown(index, start, end);
		ll mid = (start + end) / 2;
		Node l, r, ans;
		l = query(start, mid, 2 * index, left, right);
		r = query(mid + 1, end, 2 * index + 1, left, right);
		ans.merge(l, r);
		return ans;
	}

	void make_update(ll left, ll right, ll val) {
		Update new_update = Update(val);
		update(0, n - 1, 1, left, right, new_update);
	}

	Node make_query(ll left, ll right) {
		return query(0, n - 1, 1, left, right);
	}
};

struct Node1 {
	ll val;
	Node1() {
		val = 0;
	}
	Node1(ll p1) {
		val = p1;
	}
	void merge(Node1 &l, Node1 &r) {
		val = (l.val + r.val);
	}
};

struct Update1 {
	ll val;
	Update1() {
		val = 0;
	}
	Update1(ll val1) {
		val = val1;
	}
	void apply(Node1 &a, ll start, ll end) {
		a.val = (a.val * val) % MOD;
	}
	void combine(Update1& new_update, ll start, ll end) {
		val = (val * new_update.val) % MOD;
	}
};

void solve() {
	ll n, q;
	cin >> n >> q;
	vector<ll> arr(n);
	for (ll i = 0; i < n; i++) {
		cin >> arr[i];
	}

	LazySGT<Node1, Update1> seg1 = LazySGT<Node1, Update1>(n, arr);

	while (q--) {
		ll type;
		cin >> type;
		if (type == 1) {
			ll lb, rb, val;
			cin >> lb >> rb >> val;
			seg1.make_update(lb, rb - 1, val);
		}
		else {
			ll l, r;
			cin >> l >> r;
			Node1 qry = seg1.make_query(l, r - 1);
			cout << qry.val << endl;
		}
	}

}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif
	ll tc;
	tc = 1;
	// cin >> tc;
	for (ll i = 1; i <= tc; i++) {
		solve();
	}

	return 0;
}

