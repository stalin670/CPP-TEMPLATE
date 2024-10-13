//KNOWLEDGE BASE

//#pragma pack(1)
//TEMPLATE
// #include<ext/pb_ds/assoc_container.hpp>
// #include<ext/pb_ds/tree_policy.hpp>

// using namespace __gnu_pbds;

// #define pb(a) push_back(a);
// #define op(a) pop_back(a);
// #define ft front;
// #define bk back;
// #define first ff;
// #define second ss;
// #define ub upper_bound;
// #define lb lower_bound;
// #define vs vector<str>
// #define vs vector<str>
// #define forij(j, n) for(ll j = 0; j < n; j++)
// #define forijm(j, m) for(ll j = 0; j < m; j++)
// #define foriji(j, i, n) for(ll j = i + 1; j < n; j++)
// #define ford(i, n) for(ll i = n - 1; i >= 0; i--)
// #define forii(i, n) for(int i = 0; i < n; i++)
// #define fordi(i, n) for(int i = n - 1; i >= 0; i--)

// typedef tree<int, null_type,less<int>, rb_tree_tag,tree_order_statistics_node_update> pbds
// #define ordered_set tree<int, null_type,less<int>, rb_tree_tag,tree_order_statistics_node_update>

// using pi = pair<int, int>;
// using vi = vector<int , int>;
// using mi = map<int, int>;
// using si = set<int, int>;

//INFORMATION
// class compare {
//     public:
//         bool operator () (const ll &a, const ll &b) return a > b;
// }
//sort(a.begin, a.end(), greater<ll>());
//sort(a.begin(), a.end(), less<ll>());
//set<ll, compare> s;
//priortiy_queue<ll, vll, compare> pq; // MIN_HEAP;
// auto lambda_function_name = [](parameters)-> return_type {function_definition};
//[&] pass all external variables with reference
//[=] pass all external variables with value
//[a, &b] pass variable a with the value(pass by value) and b(pass by reference)

// map<pair<int, int>, int> mymap;
// mymap.insert(make_pair(make_pair(key1,key2), value));
//access mymap[make_pair(key1, key2)];
// How to print the map of pairs
// for(auto &it: mp) { // it always took key as a pair
//     auto key_pair = it.first;
//     cout << key_pair.first << " " << key_pair.second << " " << it.second << endl;
// }
// for(const auto& [key, value]: unordered_map) {
//   cout << key << " -  " << value << endl;
// }

// int logi(int a, int b) {
//     int t = 0;
//     i64 v = 1;
//     while (v < b) {
//         v *= a;
//         t++;
//     }
//     return t;
// }

// int llog(int a, int b) {
//     if (a <= b) {
//         int l = logi(a, b);
//         return (l == 0 ? 0 : std::__lg(2 * l - 1));
//     }
//     int l = logi(b, a + 1) - 1;
//     assert(l > 0);
//     return -std::__lg(l);
// }

// CODE SNIPPETS
// Binary Exponentiation using mod operation
// ll binaryexpo(ll base, ll x) {
// 	ll ans = 1;
// 	while(x) {
// 		if(x % 2) {
// 			ans = (ans * base) % mod;
// 			x = x - 1;
// 		}
// 		else {
// 			base = (base * base) % mod;
// 			x = x / 2;
// 		}
// 	}
// 	return ans;
// }

//to find only modular inverse value
// ll inv(int a) {
//     return a <= 1 ? a : mod - (ll)(mod / a) * inv(mod % a) % mod;
// }

// Combinatorics
// ll fact[1e6 + 4];
// ll modinv[1e6 + 4];
// void precompfact() {
// 	modinv[0] = 1;
// 	fact[0] = 1;
// 	for(ll i = 1; i <= 1e6; i++) {
// 		fact[i] = (fact[i - 1] * i) % mod;
// 		modinv[i] = binaryexpo(fact[i], mod - 2);
// 	}
// }

// ll ncr(ll n, ll r) {
// 	if(n < 0 || r < 0 || r > n) return 0;
// 	return(((fact[n] * modinv[r]) % mod) * modinv[n - r]) % mod;
// }


// full nCr
// int  factorialNumInverse[N + 1];

// // array to precompute inverse of 1! to N!
// int  naturalNumInverse[N + 1];

// // array to store factorial of first N numbers
// int  fact[N + 1];

// // Function to precompute inverse of numbers
// void InverseofNumber(int  p)
// {
//     naturalNumInverse[0] = naturalNumInverse[1] = 1;
//     for (int i = 2; i <= N; i++)
//         naturalNumInverse[i] = naturalNumInverse[p % i] * (p - p / i) % p;
// }
// // Function to precompute inverse of factorials
// void InverseofFactorial(int  p)
// {
//     factorialNumInverse[0] = factorialNumInverse[1] = 1;

//     // precompute inverse of natural numbers
//     for (int i = 2; i <= N; i++)
//         factorialNumInverse[i] = (naturalNumInverse[i] * factorialNumInverse[i - 1]) % p;
// }

// // Function to calculate factorial of 1 to N
// void factorial(int  p)
// {
//     fact[0] = 1;

//     // precompute factorials
//     for (int i = 1; i <= N; i++) {
//         fact[i] = (fact[i - 1] * i) % p;
//     }
// }

// // Function to return nCr % p in O(1) time
// int  Binomial(int  N, int  R, int  p) {
//     int  ans = ((fact[N] * factorialNumInverse[R])
//               % p * factorialNumInverse[N - R])
//              % p;
//     return ans;
// }

//Without Rank
// class DisjointSet {
// public:
//     vector<int> parent, size;
//     DisjointSet(int n) {
//         parent.resize(n+1);
//         size.resize(n+1);
//         for(int i = 0;i<=n;i++) {
//             parent[i] = i;
//             size[i] = 1;
//         }
//     }

//     int findUPar(int node) {
//         if(node == parent[node])
//             return node;
//         return parent[node] = findUPar(parent[node]);
//     }

//     void unionBySize(int u, int v) {
//         int ulp_u = findUPar(u);
//         int ulp_v = findUPar(v);
//         if(ulp_u == ulp_v) return;
//         if(size[ulp_u] < size[ulp_v]) {
//             parent[ulp_u] = ulp_v;
//             size[ulp_v] += size[ulp_u];
//         }
//         else {
//             parent[ulp_v] = ulp_u;
//             size[ulp_u] += size[ulp_v];
//         }
//     }
// };

//With Rank
// class DisjointSet {
//     vector<int> rank, parent, size;
// public:
//     DisjointSet(int n) {
//         rank.resize(n+1, 0);
//         parent.resize(n+1);
//         size.resize(n+1);
//         for(int i = 0;i<=n;i++) {
//             parent[i] = i;
//             size[i] = 1;
//         }
//     }

//     int findUPar(int node) {
//         if(node == parent[node])
//             return node;
//         return parent[node] = findUPar(parent[node]);
//     }

//     void unionByRank(int u, int v) {
//         int ulp_u = findUPar(u);
//         int ulp_v = findUPar(v);
//         if(ulp_u == ulp_v) return;
//         if(rank[ulp_u] < rank[ulp_v]) {
//             parent[ulp_u] = ulp_v;
//         }
//         else if(rank[ulp_v] < rank[ulp_u]) {
//             parent[ulp_v] = ulp_u;
//         }
//         else {
//             parent[ulp_v] = ulp_u;
//             rank[ulp_u]++;
//         }
//     }

//     void unionBySize(int u, int v) {
//         int ulp_u = findUPar(u);
//         int ulp_v = findUPar(v);
//         if(ulp_u == ulp_v) return;
//         if(size[ulp_u] < size[ulp_v]) {
//             parent[ulp_u] = ulp_v;
//             size[ulp_v] += size[ulp_u];
//         }
//         else {
//             parent[ulp_v] = ulp_u;
//             size[ulp_u] += size[ulp_v];
//         }
//     }
// };

// ll get_hash(str s) {
//     ll h = 0;
//     for(char c: s) h = (h * 31 + (c - 'a' + 1)) % mod;
//     return h;
// }

// ll power(ll x, ll y) {
//     ll res = 1;
//     x = x % mod;
//     while(y > 0) {
//         if(y & 1)
//             res = (res * x) % mod;
//         y = y >> 1;
//         x = (x * x) % mod;
//     }
//     return res;
// }


//Sieve of Eratosthenes
// vll is_prime(10000001, 1);
// for(ll i = 2; i <= 10000000; i++) {
//     if(is_prime[i] == 1) {
//         for(ll j = i * i; j <= 10000000; j += i)
//             is_prime[j] = 0;
//     }
// }

// vector<int> is_prime(100001, 1);
// void SieveOfEratosthenes(int n)
// {
//     is_prime[0] = is_prime[1] = 0;
//     for (int i = 2; i * i <= 100000; i++) {
//         if (is_prime[i]) {
//             for (int j = i * i; j <= 100000; j += i)
//                 is_prime[j] = 1;
//         }
//     }
// }

//SPF
// ll maxi = 1e6;
// vll spf(maxi + 1, 0);
// for(ll i = 2; i <= maxi; i++) spf[i] = i;
// for(ll i = 2; i * i <= maxi; i++) {
// 	if(sieve[i] != i) continue;
// 	for(ll j = i * i; j <= maxi; j += i) {
// 		if(spf[j] == j) spf[j] = i;
// 	}
// }

// bool isPrime(ll n) {
//     if(n <= 1) return false;
//     for(ll i = 2; i * i <= n; i++) {
//         if(n % i == 0) return false;
//     }
//     return true;
// }

// vll prifac[100001];
// for(ll i = 1; i <= 100000; i++) {
//     for(ll j = i; j <= 100000; j += i) {
//         prifac[j].push_back(i);
//     }
// }


//Print the all only unique prime factors of n
//If you want to store instead of print store in some data Structure
// ll primeFactors(ll n) {
//     ll c = 2;
//     ll pre = -1;
//     while(n > 1) {
//         if(n % c == 0) {
//             if(pre != c) {
//                 cout << c << endl;
//             }
//             n /= c;
//             pre = c;
//         }
//         else
//             c++;
//     }
//     cout << endl;
//     return 1;
// }

//Find all only unique prime factors of n in optimized way
//and stored in vector and return it
// vll wheel(ll n) {
//     vll factorization;
//     for (ll d : {2, 3, 5}) {
//         if(n % d == 0)
//             factorization.push_back(d);
//         while(n % d == 0) {
//             n /= d;
//         }
//     }
//     static array<int, 8> increments = {4, 2, 4, 2, 4, 6, 2, 6};
//     ll i = 0;
//     for (ll d = 7; d * d <= n; d += increments[i++]) {
//         if(n % d == 0)
//             factorization.push_back(d);
//         while (n % d == 0) {
//             n /= d;
//         }
//         if (i == 8)
//             i = 0;
//     }
//     if (n > 1)
//         factorization.push_back(n);
//     return factorization;
// }



//Find the all only prime factors of n with repitions
//and they will stored in vector v
// void primeFactors(ll n, vll &v) {
//     while (n % 2 == 0) {
//         v.push_back(2);
//         n = n/2;
//     }

//     for (ll i = 3; i <= sqrt(n); i = i + 2) {
//         while (n % i == 0) {
//             v.push_back(1);
//             n = n/i;
//         }
//     }

//     if (n > 2)
//         v.push_back(n);
// }

//Print all divisors of n and store in map
// void printDivisors(ll n, mll &mp) {
//     for(ll i = 1; i <= sqrt(n); i++) {
//         if(n % i == 0) {
//             if (n / i == i)
//                 mp[i]++;
//             else {
//                 mp[i]++;
//                 mp[n / i]++;
//             }
//         }
//     }
// }

// bool balance_str(string expr) {
//     stack<char> temp;
//     for (ll i = 0; i < expr.length(); i++) {
//         if (temp.empty()) temp.push(expr[i]);
//         else if ((temp.top() == '(' && expr[i] == ')') || (temp.top() == '{' && expr[i] == '}') || (temp.top() == '[' && expr[i] == ']')) temp.pop();
//         else temp.push(expr[i]);
//     }
//     if (temp.empty()) return true;
//     return false;
// }

//Range Minimum Query use for min, max, gcd (if there is not update in array) // Time Complexity - O(N(log(N) + Q))
// int n;
// cin >> n;
// bin_log[1] = 0;
// for(int i = 2; i <= n; i++) {
// 	bin_log[i] = bin_log[i/2]+1;
// }
// for(int i = 0; i < n; i++) {
// 	cin >> a[i];
// 	m[i][0] = a[i];
// }
// // 2) preprocessing, O(N*log(N))
// for(int k = 1; k < LOG; k++) {
// 	for(int i = 0; i + (1 << k) - 1 < n; i++) {
// 		m[i][k] = min(m[i][k-1], m[i+(1<<(k-1))][k-1]);
// 	}
// }
// 3) answer queries
// int q;
// cin >> q;
// for(int i = 0; i < q; i++) {
// 	int L, R;
// 	cin >> L >> R;
// 	cout << query(L, R) << "\n";
// }

// const int MAX_N = 100'005;
// const int LOG = 17;
// int a[MAX_N];
// int m[MAX_N][LOG]; // m[i][j] is minimum among a[i..i+2^j-1]
// int bin_log[MAX_N];

// int query(int L, int R) { // O(1)
// 	int length = R - L + 1;
// 	int k = bin_log[length];
// 	return min(m[L][k], m[R-(1<<k)+1][k]);
// }


//Range Minimum Query use for min, max, gcd (if there is not update in array)  with space optimization
// vector<Query> bucket[1 + lgmax];
// int lg(int number) {
//   return 31 - __builtin_clz(number);
// }
// int rmq[5 + nmax], sol[5 + nmax];

// int n, q;
//   std::cin >> n >> q;
//   for(int i = 1; i <= n; i++)
//     std::cin >> rmq[i];
//   for(int i = 1;i <= q; i++) {
//     int x, y;
//     std::cin >> x >> y;
//     bucket[lg(y - x + 1)].push_back({x, y, i});
//   }

//   for(int h = 0; h < lgmax; h++) {
//     for(int i = 0; i < bucket[h].size(); i++) {
//       int x = bucket[h][i].x;
//       int y = bucket[h][i].y;
//       int id = bucket[h][i].id;
//       sol[id] = std::min(rmq[x], rmq[y - (1 << h) + 1]);
//     }

//     for(int j = 1;j <= n - (1 << (h + 1)) + 1; j++)
//       rmq[j] = std::min(rmq[j], rmq[j + (1 << h)]);
//   }

//   for(int i = 1;i <= q; i++)
//     std::cout << sol[i] << '\n';

//generic segment Tree
// this is how you will create segment tree and call build function
// this is zero based indexing segmeent tree
// Lsegtree<ll, ll> st(n, 0, 0);
// st.build(a);
// template<class T, class U>
// struct Lsegtree {
// 	vector<T> st;
// 	vector<U> lazy;
// 	ll n;
// 	T identity_element;
// 	U identity_update;
//Identity element
//in case of no overlap, the value of identity_element
//should be such that the answer remains unchanged.
//Identity update
//such a value, that if we apply this update to a node,
// thenthe value stored in that node does not change
// 	Lsegtree(ll n, T identity_element, U identity_update) {
// 		this -> n = n;
// 		this -> identity_element = identity_element;
// 		this -> identity_update = identity_update;
// 		st.assign(4 * n, identity_element);
// 		lazy.assign(4 * n, identity_update);
// 	}
// 	T combine(T l, T r) {
// 		//first location where we need to change
// 		T ans = (l + r);
// 		return ans;
// 	}
// 	void buildUtil(ll v, ll tl, ll tr, vector<T> &a) {
// 		if(tl == tr) {
// 			st[v] = a[tl];
// 			return;
// 		}
// 		ll tm = (tl + tr) >> 1;
// 		buildUtil(2 * v + 1, tl, tm, a);
// 		buildUtil(2 * v + 2, tm + 1, tr, a);
// 		st[v] = combine(st[2 * v + 1], st[2 * v + 2]);
// 	}
//	void build(vector<T> a) {
// assert((ll)a.size() == n);
// buildUtil(0, 0, n - 1, a);
//	}
// 	T apply(T curr, U upd, ll tl, ll tr) {
//      second location we need to change
// 		T ans = (tr - tl + 1) * upd;
// 		return ans;
// 	}
// 	U combineUpdate(U old_upd, U new_upd, ll tl, ll tr) {
//      third location we need to change
// 		U ans = old_upd;
// 		ans = new_upd;
// 		return ans;
// 	}
// 	void push_down(ll v, ll tl, ll tr) {
// 		if(lazy[v] == identity_update) return;
// 		st[v] = apply(st[v], lazy[v], tl, tr);
// 		if(2 * v + 1 <= 4 * n) {
// 			ll tm = (tl + tr) >> 1;
// 			lazy[2 * v + 1] = combineUpdate(lazy[2 * v + 1], lazy[v], tl, tm);
// 			lazy[2 * v + 2] = combineUpdate(lazy[2 * v + 2], lazy[v], tm + 1, tr);
// 		}
// 		lazy[v] = identity_update;
// 	}
// 	T queryUtil(ll v, ll tl, ll tr, ll l, ll r) {
// 		push_down(v, tl, tr);
// 		if(l > r) return identity_element;
// 		if(tr < l or tl > r) return identity_element;
// 		if(l <= tl and r >= tr) return st[v];
// 		ll tm = (tl + tr) >> 1;
// 		return combine(queryUtil(2 * v + 1, tl , tm, l, r), queryUtil(2 * v + 2, tm + 1, tr, l, r));
// 	}
// T query(ll l, ll r) {
// return queryUtil(0, 0, n - 1, l, r);
// }
// 	void updateUtil(ll v, ll tl, ll tr, ll l, ll r, U upd) {
// 		push_down(v, tl, tr);
// 		if(tr < l or tl > r) return;
// 		if(tl >= l and tr <= r) {
// 			lazy[v] = combineUpdate(lazy[v], upd, tl, tr);
// 			push_down(v, tl, tr);
// 		}
// 		else {
// 			ll tm = (tl + tr) >> 1;
// 			updateUtil(2 * v + 1, tl, tm, l, r, upd);
// 			updateUtil(2 * v + 2, tm + 1, tr, l, r, upd);
// 			st[v] = combine(st[2 * v + 1], st[2 * v + 2]);
// 		}
// 	}
// void update(ll l, ll r, U upd) {
// updateUtil(0, 0, n - 1, l, r, upd);
// }
// };

//Segment Tree (all associative operations are allowed)
// class SGTree {
// 	vector<int> seg;
// public:
// 	SGTree(int n) {
// 		seg.resize(4 * n + 1);
// 	}

// 	void build(int ind, int low, int high, int arr[]) {
// 		if (low == high) {
// 			seg[ind] = arr[low];
// 			return;
// 		}

// 		int mid = (low + high) / 2;
// 		build(2 * ind + 1, low, mid, arr);
// 		build(2 * ind + 2, mid + 1, high, arr);
// 		seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
// 	}

// 	int query(int ind, int low, int high, int l, int r) {
// 		if (r < low || high < l) return INT_MAX;

// 		if (low >= l && high <= r) return seg[ind];

// 		int mid = (low + high) >> 1;
// 		int left = query(2 * ind + 1, low, mid, l, r);
// 		int right = query(2 * ind + 2, mid + 1, high, l, r);
// 		return min(left, right);
// 	}
// 	void update(int ind, int low, int high, int i, int val) {
// 		if (low == high) {
// 			seg[ind] = val;
// 			return;
// 		}

// 		int mid = (low + high) >> 1;
// 		if (i <= mid) update(2 * ind + 1, low, mid, i, val);
// 		else update(2 * ind + 2, mid + 1, high, i, val);
// 		seg[ind] = min(seg[2 * ind + 1], seg[2 * ind + 2]);
// 	}
// };


//lazy propagation in segment Tree
//this code if for range max queries and range increment updates
// int lazy[4 * n];
// int tr[4 * n];
// ll query(ll v, ll tl, ll tr, ll l, ll r) {
// 	if(lazy[v] > 0) {
// 		tr[v] += lazy[v];
//      in case range increment queries tr[v] += lazy[v] * (tr - tl + 1);
// 		if(tl != tr) {
// 			lazy[2 * v] += lazy[v];
// 			lazy[2 * v + 1] += lazy[v];
// 		}
// 		lazy[v] = 0;
// 	}

// 	if(l > r) return 0;
// 	if(tr < l || tl > r) return 0; //no overlap

// 	if(tl >= l && tr <= r) return tr[v] //full overlap
// 	else {
// 		ll tm = (tl + tr) / 2;
// 		max(query(2 * v, l, tm, l, r), query(2 * v + 1, tm + 1, tr, l, r));
// 	}
// }
// ll update(ll v, ll tl, ll tr, ll l, ll r, ll x) {
// 	if(lazy[v] > 0) {
// 		tr[v] += lazy[v];
//		in case range increment queries tr[v] += lazy[v] * (tr - tl + 1);
// 		if(tl != tr) {
// 			lazy[2 * v] += lazy[v];
// 			lazy[2 * v + 1] += lazy[v];
// 		}
// 		lazy[v] = 0;
// 	}

// 	if(tr < l or tl > r) return; //no overlap

// 	if(tl >= l and tr <= r) { //full overlap
// 		tr[v] += x;
//		in case range increment queries tr[v] += x * (tr - tl + 1);
// 		lazy[2 * v] += x;
// 		lazy[2 * v + 1] += x;
// 	}
// 	else {
// 		ll tem = (tl + tr) / 2;
// 		update(2 * v, l, tm, l, r, x);
// 		update(2 * v + 1, tm + 1, tr, l, r, x);
// 	}
// }




//wanted to check length of cycles in graph
// template<class Fun> class y_combinator_result {
//     Fun fun_;
// public:
//     template<class T> explicit y_combinator_result(T &&fun): fun_(std::forward<T>(fun)) {}
//     template<class ...Args> decltype(auto) operator()(Args &&...args) { return fun_(std::ref(*this), std::forward<Args>(args)...); }
// };
// template<class Fun> decltype(auto) y_combinator(Fun &&fun) { return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun)); }

// struct functional_graph {
//     int n = 0;
//     vector<int> F;
//     vector<bool> in_cycle;
//     vector<int> which_cycle;
//     vector<int> cycle_root;
//     vector<int> cycle_position;
//     vector<int> depth;
//     vector<vector<int>> children;
//     vector<vector<int>> cycles;
//     vector<int> jump, jump_dist;

//     functional_graph(const vector<int> &f = {}) {
//         if (!f.empty())
//             build(f);
//     }

//     void build(const vector<int> &f) {
//         F = f;
//         n = int(F.size());
//         in_cycle.assign(n, false);
//         which_cycle.assign(n, -1);
//         cycle_root.assign(n, -1);
//         cycle_position.assign(n, -1);
//         depth.assign(n, -1);
//         children.assign(n, {});
//         cycles.clear();

//         vector<int> seen(n, -1);

//         for (int i = 0; i < n; i++) {
//             if (seen[i] >= 0)
//                 continue;

//             int x = i;

//             do {
//                 seen[x] = i;
//                 x = F[x];
//             } while (seen[x] < 0);

//             if (seen[x] != i)
//                 continue;

//             vector<int> cycle;
//             int y = x;

//             do {
//                 in_cycle[y] = true;
//                 cycle_position[y] = int(cycle.size());
//                 which_cycle[y] = int(cycles.size());
//                 cycle.push_back(y);
//                 y = F[y];
//             } while (y != x);

//             cycles.push_back(cycle);
//         }

//         seen.assign(n, 0);

//         for (int i = 0; i < n; i++)
//             seen[F[i]] += !in_cycle[i];

//         for (int i = 0; i < n; i++)
//             children[i].reserve(seen[i]);

//         seen.clear();

//         for (int i = 0; i < n; i++)
//             if (!in_cycle[i])
//                 children[F[i]].push_back(i);

//         jump.assign(n, -1);
//         jump_dist.assign(n, -1);
//         int cyc_root = -1, which_cyc = -1;

//         auto dfs = y_combinator([&](auto self, int node, int parent) -> void {
//             depth[node] = parent < 0 ? 0 : depth[parent] + 1;
//             cycle_root[node] = cyc_root;
//             which_cycle[node] = which_cyc;
//             jump[node] = parent < 0 ? node : jump_dist[parent] == jump_dist[jump[parent]] ? jump[jump[parent]] : parent;
//             jump_dist[node] = depth[node] - depth[jump[node]];

//             for (int child : children[node])
//                 self(child, node);
//         });

//         for (int i = 0; i < n; i++)
//             if (in_cycle[i]) {
//                 cyc_root = i;
//                 which_cyc = which_cycle[i];
//                 dfs(i, -1);
//             }
//     }

//     int cycle_length(int node) const {
//         return int(cycles[which_cycle[node]].size());
//     }

//     int go_forward(int v, int64_t k) const {
//         if (k >= depth[v]) {
//             int root = cycle_root[v];
//             k -= depth[v];
//             int cyc = which_cycle[root];
//             int64_t position = (cycle_position[root] + k) % cycles[cyc].size();
//             return cycles[cyc][position];
//         }

//         while (k > 0)
//             if (jump_dist[v] <= k) {
//                 k -= jump_dist[v];
//                 v = jump[v];
//             } else {
//                 k--;
//                 v = F[v];
//             }

//         return v;
//     }
// };

//same copied above segment tree little bit modify
// template<class T, class U>
// struct Lsegtree {
//     vector<T> st;
//     vector<U> lazy;
//     ll n;
//     T identity_element;
//     U identity_update;

//     Lsegtree(ll n, T identity_element, U identity_update) {
//         this -> n = n;
//         this -> identity_element = identity_element;
//         this -> identity_update = identity_update;
//         st.assign(4 * n, identity_element);
//         lazy.assign(4 * n, identity_update);
//     }
//     T combine(T l, T r) {
//         //first location where we need to change
//         T ans = (l + r);
//         return ans;
//     }
//     void buildUtil(ll v, ll tl, ll tr, vector<T> &a) {
//         if(tl == tr) {
//             st[v] = a[tl];
//             return;
//         }
//         ll tm = (tl + tr) >> 1;
//         buildUtil(2 * v + 1, tl, tm, a);
//         buildUtil(2 * v + 2, tm + 1, tr, a);
//         st[v] = combine(st[2 * v + 1], st[2 * v + 2]);
//     }
//     void build(vector<T> a) {
//         assert((ll)a.size() == n);
//         buildUtil(0, 0, n - 1, a);
//     }
//     T apply(T curr, U upd, ll tl, ll tr) {
//         // second location we need to change
//         T ans = curr + (tr - tl + 1) * upd;
//         return ans;
//     }
//     U combineUpdate(U old_upd, U new_upd, ll tl, ll tr) {
//         // third location we need to change
//         U ans = old_upd + new_upd;
//         return ans;
//     }
//     void push_down(ll v, ll tl, ll tr) {
//         if(lazy[v] == identity_update) return;
//         st[v] = apply(st[v], lazy[v], tl, tr);
//         if(2 * v + 1 <= 4 * n) {
//             ll tm = (tl + tr) >> 1;
//             lazy[2 * v + 1] = combineUpdate(lazy[2 * v + 1], lazy[v], tl, tm);
//             lazy[2 * v + 2] = combineUpdate(lazy[2 * v + 2], lazy[v], tm + 1, tr);
//         }
//         lazy[v] = identity_update;
//     }
//     T queryUtil(ll v, ll tl, ll tr, ll l, ll r) {
//         push_down(v, tl, tr);
//         if(l > r) return identity_element;
//         if(tr < l or tl > r) return identity_element;
//         if(l <= tl and r >= tr) return st[v];
//         ll tm = (tl + tr) >> 1;
//         return combine(queryUtil(2 * v + 1, tl , tm, l, r), queryUtil(2 * v + 2, tm + 1, tr, l, r));
//     }
//     T query(ll l, ll r) {
//         return queryUtil(0, 0, n - 1, l, r);
//     }
//     void updateUtil(ll v, ll tl, ll tr, ll l, ll r, U upd) {
//         push_down(v, tl, tr);
//         if(tr < l or tl > r) return;
//         if(tl >= l and tr <= r) {
//             lazy[v] = combineUpdate(lazy[v], upd, tl, tr);
//             push_down(v, tl, tr);
//         }
//         else {
//             ll tm = (tl + tr) >> 1;
//             updateUtil(2 * v + 1, tl, tm, l, r, upd);
//             updateUtil(2 * v + 2, tm + 1, tr, l, r, upd);
//             st[v] = combine(st[2 * v + 1], st[2 * v + 2]);
//         }
//     }
//     void update(ll l, ll r, U upd) {
//         updateUtil(0, 0, n - 1, l, r, upd);
//     }
// };

//operation with mod
// int mod = 1e9 + 7;

// int add(int a, int b) {
//     return a + b >= mod ? a + b - mod : a + b;
// }
// int sub(int a, int b) {
//     return a >= b ? a - b : a + mod - b;
// }
// int mul(int a, int b) {
//     return 1ll * a * b % mod;
// }
// int inv(int x) {
//     return x == 1 ? 1 : sub(0, mul(mod / x, inv(mod % x)));
// }
