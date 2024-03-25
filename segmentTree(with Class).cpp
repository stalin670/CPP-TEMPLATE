class SGTree {
   public:
      vector<ll> segTree;
   public:
      SGTree(ll n) {
         segTree.resize(4 * n + 1);
      }

      void build_seg(ll ind, ll low, ll high, vi arr) {
         if(low == high) {
            segTree[ind] = arr[low];
            return;
         }

         ll mid = low + (high - low) / 2;
         build_seg(2 * ind + 1, low , mid, arr);
         build_seg(2 * ind + 2, mid + 1 , high, arr);

         segTree[ind] = segTree[2 * ind + 1] + segTree[2 * ind + 2];
      }

      ll query_seg(ll ind, ll low, ll high, ll l, ll r) {
         // No Overlap
         if(r < low or high < l)
            return 0;

         // complete Overlap
         if(low >= l and high <= r)
            return segTree[ind];

         ll mid = low + (high - low) / 2;
         ll left = query_seg(2 * ind + 1, low, mid, l, r);
         ll right = query_seg(2 * ind + 2, mid + 1, high, l, r);
         return left + right;
      }

      void update_seg(ll ind, ll low, ll high, ll point, ll val) {
         if(low == high) {
            segTree[ind] = val;
            return;
         }

         ll mid = low + (high - low) / 2;
         if(point <= mid)
            update_seg(2 * ind + 1, low, mid, point, val);
         else
            update_seg(2 * ind + 2, mid + 1, high, point, val);

         segTree[ind] = segTree[2 * ind + 1] + segTree[2 * ind + 2];         
      }

};


void solve(){
   ll n, q;
   cin >> n >> q;
   vi arr(n);
   for(ll i = 0; i < n; i++)
      cin >> arr[i];

   SGTree st(n);
   st.build_seg(0, 0, n - 1, arr);

   while(q--) {
      ll l, r;
      cin >> l >> r;
      l--;
      r--;
      cout << st.query_seg(0, 0, n - 1, l, r) << endl;
   }  

   return;
}
