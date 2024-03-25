const ll maxs = 1e7;

ll segTree[maxs];
vi a(maxs);

// Update

void build_seg(ll ind, ll lb, ll rb) {

   if(lb == rb) {
      segTree[ind] = a[lb];
      return;
   }

   ll mb = lb + (rb - lb) / 2;
   build_seg(2 * ind + 1, lb , mb);
   build_seg(2 * ind + 2, mb + 1, rb);

   segTree[ind] = segTree[2 * ind + 1] + segTree[2 * ind + 2];
}

// Querying

ll segQuery(ll ind, ll lb, ll rb, ll l, ll r) {
   if(rb < l or r < lb)
      return 0;

   if(lb >= l and rb <= r)
      return segTree[ind];

   ll mb = lb + (rb - lb) / 2;
   ll left = segQuery(2 * ind + 1, lb , mb, l, r);
   ll right = segQuery(2 * ind + 2, mb + 1, rb, l, r);

   return left + right;
}

// Updating

void segUpdate(ll ind, ll lb, ll rb, ll point, ll val) {
   if(lb == rb) {
      segTree[ind] = val;
      return;
   }

   ll mb = lb + (rb - lb) / 2;
   if(point <= mb)
      segUpdate(2 * ind + 1, lb, mb, point, val);
   else 
      segUpdate(2 * ind + 2, mb + 1, rb, point, val);

   segTree[ind] = segTree[2 * ind + 1] + segTree[2 * ind + 2];

   return;
}
