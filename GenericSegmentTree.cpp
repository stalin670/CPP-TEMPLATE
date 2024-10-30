// G.cpp

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
struct SegTree {
	vector<Node> tree;
	vector<ll> arr;
	int n;
	int s;

	SegTree(int a_len, vector<ll> &a) {
		arr = a;
		n = a_len;
		s = 1;
		while (s < 2 * n) {
			s = s << 1;
		}
		tree.resize(s);
		fill(all(tree), Node());
		build(0, n - 1, 1);
	}

	void build(int start, int end, int index)
	{
		if (start == end)	{
			tree[index] = Node(arr[start]);
			return;
		}
		int mid = (start + end) / 2;
		build(start, mid, 2 * index);
		build(mid + 1, end, 2 * index + 1);
		tree[index].merge(tree[2 * index], tree[2 * index + 1]);
	}

	void update(int start, int end, int index, int query_index, Update &u)
	{
		if (start == end) {
			u.apply(tree[index]);
			return;
		}
		int mid = (start + end) / 2;
		if (mid >= query_index)
			update(start, mid, 2 * index, query_index, u);
		else
			update(mid + 1, end, 2 * index + 1, query_index, u);
		tree[index].merge(tree[2 * index], tree[2 * index + 1]);
	}

	Node query(int start, int end, int index, int left, int right) {
		if (start > right || end < left)
			return Node();
		if (start >= left && end <= right)
			return tree[index];
		int mid = (start + end) / 2;
		Node l, r, ans;
		l = query(start, mid, 2 * index, left, right);
		r = query(mid + 1, end, 2 * index + 1, left, right);
		ans.merge(l, r);
		return ans;
	}

	void make_update(int index, ll val) {
		Update new_update = Update(val);
		update(0, n - 1, 1, index, new_update);
	}

	Node make_query(int left, int right) {
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
		val = l.val ^ r.val;
	}
};

struct Update1 {
	ll val;
	Update1(ll p1) {
		val = p1;
	}
	void apply(Node1 &a) {
		a.val = val;
	}
};


void solve() {
	ll n, q;
	cin >> n >> q;
	vector<ll> arr(n);
	for (ll i = 0; i < n; i++) {
		cin >> arr[i];
	}

	SegTree<Node1, Update1> seg1 = SegTree<Node1, Update1>(n, arr);

	for (ll i = 0; i < q; i++) {
		ll lb, rb;
		cin >> lb >> rb;
		lb--, rb--;
		Node1 ans = seg1.make_query(lb, rb);
		cout << ans.val << endl;
	}
}

int main() {
#ifndef ONLINE_JUDGE
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
#endif
	int tc;
	tc = 1;
	// cin >> tc;
	for (int i = 1; i <= tc; i++) {
		solve();
	}

	return 0;
}

