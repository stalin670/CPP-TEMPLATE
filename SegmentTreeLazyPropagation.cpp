class SGTree {
   public:
      vector<ll> segTree, lazyTree;
   public:
      SGTree(ll n) {
         segTree.resize(4 * n + 1);
         lazyTree.resize(4 * n + 1);
      }

      void build_seg(ll ind, ll low, ll high, vi & arr) {
         if(low == high) {
            segTree[ind] = arr[low];
            return;
         }

         ll mid = low + (high - low) / 2;
         build_seg(2 * ind + 1, low , mid, arr);
         build_seg(2 * ind + 2, mid + 1 , high, arr);

         segTree[ind] = segTree[2 * ind + 1] + segTree[2 * ind + 2];
      }

      void update_lazy(ll ind, ll low, ll high) {
         if(low == high) {
            segTree[ind] += lazyTree[ind];
            lazyTree[ind] = 0;
         }
         else {
            segTree[ind] += (high - low + 1) * lazyTree[ind];
            lazyTree[2 * ind + 1] += lazyTree[ind];
            lazyTree[2 * ind + 2] += lazyTree[ind];
            lazyTree[ind] = 0;
         }
      }

      ll query_seg(ll ind, ll low, ll high, ll l, ll r) {
         
         update_lazy(ind, low, high);

         // No Overlap
         if(r < low or high < l or low > high)
            return 0;

         if(low == high)
            return segTree[ind];

         // complete Overlap
         if(low >= l and high <= r)
            return segTree[ind];

         ll mid = low + (high - low) / 2;
         if(l <= mid)
            return query_seg(2 * ind + 1, low, mid, l, r);
         else
            return query_seg(2 * ind + 2, mid + 1, high, l, r);
      }

      void update_seg(ll ind, ll low, ll high, ll l, ll r, ll val) {
         
         update_lazy(ind, low, high);

         // No overlap
         if(high < l or r < low)
            return;

         // complete overlap 
         // l low high r 
         if(low >= l and high <= r) {
            segTree[ind] += (high - low + 1) * val; 
            // if a leaf node, it will have childrens
            if(low != high) {
               lazyTree[2 * ind + 1] += val; 
               lazyTree[2 * ind + 2] += val; 
            }
            return; 
         }

         ll mid = low + (high - low) / 2;
         update_seg(2 * ind + 1, low, mid, l, r, val);
         update_seg(2 * ind + 2, mid + 1, high, l, r, val);

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
      ll type;
      cin >> type;
      if(type == 1) {
         ll a, b, val;
         cin >> a >> b >> val;
         a--;
         b--;
         st.update_seg(0, 0, n - 1, a, b, val);
      }
      else {
         ll l;
         cin >> l;
         l--;
         ll ans = st.query_seg(0, 0, n - 1, l, l);
         cout << ans << endl;
      }
   }  

   return;
}
