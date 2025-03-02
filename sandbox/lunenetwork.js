// LICENSE: CC0
//
// To the extent possible under law, the person who associated CC0 with
// this project has waived all copyright and related or neighboring rights
// to this project.
// 
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

// Some motivation for this program:
//
// The Relative Neighborhood Graph (RNG) is calculated from a set of points
// in (Euclidean0 space by joining any two points that don't have another
// in the middle of their lune.
//
// If the distance between two points is $d$, the lune is the intersection
// of the two spheres of radius $d$.
//
// The motivation was to create a graph that "matched the human perception
// of the shape of the set".
//
// Here, we restrict ourselves to random 2d and 3d points in Euclidean space.
// That is, we assume a Poisson process for point distribution and try to
// create relatively run-time efficient algorithms to find the RNG from
// these set of points.
//
// I'm starting with 2d to work out some proofs of concept but hope to extend
// to 3d.
//
// The naive algorithm would be to start with a fully connected graph of all
// points, then remove edges when considering all triples of points.
// This algorithm has $O(n^3)$ and can be improved on.
//
// In general, the RNG is a sub-graph of the Delaunay Triangulation (DT) and,
// in general, one can't do better than the run-time efficiency of DT.
// DT has worst case $O(n^2)$ in 2d and 3d (I believe), so that's a hard
// upper bound.
//
// To construct the RNG from the DT, you can restrict the lune tests
// to only the neighbors implied by the DT. DT has an upper bound of
// degree 6 in 2d so this yields an O(n^2) RNG algorithm in 2d.
// It looks like DT in 3d has worst case O(n) degree but the average
// edge count is bounded by O(n^2), so in the aggregate, doing the same
// greedy algorithm will yield O(n^2) to extract the RNG.
//
// When the points are in generic (e.g. random) position, some optimizations
// can be employed as the average neighbor degree is effectively bounde,
// either for the RNG or the DT.
//
// The idea that we're trying to implement here is an expected linear
// time algorithm for the construction of the RNG assuming random point placement.
//
// The key insight is to look at "wedge" regions around a point $p$ and realize
// that any other point $q$ inside this wedge region will preclude points further
// away but still in the wedge.
//
// For 2d, this means creating wedges that octants:
//              
//        \ R2| R1/
//         \  |  /
//       R3 \ | / R0
//     ------ p -----
//       R4 / | \  R7
//         /  |  \
//        /R5 | R6\
//
// Consider choosing a point, $q$, somewhere in R1.
// The lune created by by $p$ and $q$ will fully the R1
// region that is closer than dist(p,q).
// So if $q$ is chosen to be the closest point to $p$ in R1,
// we know no other points in R1 further away need be considered.
//
// The closest points in each of the R0 to R7 regions are candidates
// and might not be part of the RNG but we know if we have a candidate
// list of $q_0$ to $q_7$ (one for each region, say), only they can
// be part of the RNG connecting to $p$, and no other.
//
// Since the wedge region count is finite and if we can construct 
// the set of closest points for each $p$ considered in $O(1)$ time, we have a linear
// algorithm.
// The "expected linear time" comes from our expectation of constructing
// the candidate $q$ points for each $p$ and would be violated if
// there were special structure or an adversarial choice policy.
//
// For 3d, I'm not sure what the best analogue is, but one way, I believe (I
// have to confirm), is that you can create the analogue of the wedges in each
// of the six directions (+-x, +-y, +-z), giving 48 regions (6*8) to consider.
//
// The idea is:
//
// * construct a grid of sqrt(n) x sqrt(n) size and bin points in their appropriate
//   coarse grid cell location
// * for each point $p$:
//   - use a spiral-out pattern of traversing the grid, anchored where $p$ is and
//     keep adding the closest point encountered in each of the wedge regions
//   - when each region has a point or the boundaries have been encountered, stop
//     the spiral-out search
//   - take the list of candidates and do a naive RNG construction to find which
//     each of the $q$s point $p$ should connect to.
//
// The grid construction is O(n) ($\sqrt(n) \mult \sqrt(n)$). I'm not adept enough
// to work out the probabilities but I think occupancy works out so that expected
// run times stay constant for individual edge construction.
//
// Some notes on implementation:
//
// * for each point $p$, we find the integral grid point as a starting point
// * we do a full catalogue of points in the base grid point and grid points
//   that are 1 away
// * if we haven't found all wedge occupancy, we sweep out the boundary, at a given
//   radius, to iterate through the grid cells to try to find candidate neighbors
//   - we have to make sure to take the 'overhang' grid cell locations as the wedge
//     might eat into boundary walls near the edges
// * when occupancy has been satisfied, do a naive calculation of the RNG
//
//

var njs = require("./numeric.js");
var srand = require("./seedrandom.js");

var fs = require("fs");

var _rnd = srand("lunenetwork");

var DEBUG_LEVEL = 0;

var _debug_stat = {
  "count": {},
  "val": {}
};

var prof_ctx = {
};


var N = 10;
var pnt = [];

function poisson_point(N, D) {
  D = ((typeof D === "undefined") ? 2 : D);
  let pnt = [];
  for (let i=0; i<N; i++) {
    let p = [];
    for (let j=0; j<D; j++) {
      //p.push( Math.random() );
      p.push( _rnd() );
    }
    pnt.push(p);
  }
  return pnt;
}

function print_point(pnt, disconnect) {
  for (let i=0; i<pnt.length; i++) {
    if (pnt[i].length == 2) {
      console.log(pnt[i][0], pnt[i][1]);
    }
    else {
      console.log(pnt[i][0], pnt[i][1], pnt[i][2]);
    }
    if (disconnect) { console.log("\n"); }
  }
  console.log("");
}

function in_lune(pnt_a, pnt_b, tst_c) {

  let dist_ca = njs.norm2( njs.sub(tst_c, pnt_a) );
  let dist_cb = njs.norm2( njs.sub(tst_c, pnt_b) );
  let dist_ab = njs.norm2( njs.sub(pnt_a, pnt_b) );

  //console.log("# pnt_a(", pnt_a, "), pnt_b(", pnt_b, "), tst_c(", tst_c, ")");
  //console.log("# dist_ca:", dist_ca, "dist_cb:", dist_cb, "dist_ab:", dist_ab);

  if ((dist_ca <= dist_ab) &&
      (dist_cb <= dist_ab)) {
    return true;
  }

  return false;
}

function _lerp(v0, v1, t) {
  return v0 + ((v1-v0)*t)
}

function lune_points( a, b, seg, connect ) {
  seg = ((typeof seg === "undefined") ? 8 : seg);
  connect = ((typeof connect === "undefined") ? true : connect);

  let ba = njs.sub(b,a);
  let theta = -Math.atan2(ba[1], ba[0]);
  let r = njs.norm2(ba);

  let _pi3 = Math.PI/3;

  let p = [];

  let bxy = [ [1,0], [0,1] ];

  let p0 = [];

  for (let idx=0; idx<seg; idx++) {
    let t = _lerp( theta - _pi3, theta + _pi3, idx / seg);
    let m = [ [ Math.cos(t), -Math.sin(t) ], [ Math.sin(t), Math.cos(t) ] ];
    let bt = njs.mul( r, njs.dot(m, bxy) );
    p.push( njs.add(a, bt[0]) );

    if (p0.length == 0) {
      p0 = njs.clone(p[0]);
    }
  }

  theta += Math.PI;

  for (let idx=0; idx<seg; idx++) {
    let t = _lerp( theta - _pi3, theta + _pi3, idx / seg);
    let m = [ [ Math.cos(t), -Math.sin(t) ], [ Math.sin(t), Math.cos(t) ] ];
    let bt = njs.mul( r, njs.dot(m, bxy) );
    p.push( njs.add(b, bt[0]) );

  }

  if (connect) {
    p.push(p0);
  }

  return p;
}


// 3D vector to idir
//
function v2idir(v) {
  let max_xyz = 0;
  let max_val = v[0];
  for (let xyz=0; xyz<3; xyz++) {
    if (Math.abs(v[xyz]) > max_val) {
      max_xyz = xyz;
      max_val = Math.abs(v[xyz]);
    }
  }

  if (v[max_xyz] < 0) { return (2*max_xyz)+1; }
  return 2*max_xyz;
}

// (3D) angle between p and q vectors
//
function v3theta(p,q) {
  let s = njs.norm2( cross3(p,q) );
  let c = njs.dot(p,q);
  return Math.atan2(s,c);
}

// plane equation:
//
// P(u) = N_p \dot (u - p)
//
function plane_f(u, Np, p) {
  return njs.dot( Np, njs.sub(u, p) );
}

// Np here is 4 vector plane, where first
// three elements are normal, fourth is distance
// to origin.
//
// Return plane equation, using derived point on
// plane from reduced representation.
//
// Np = [nx, ny, nz, nw]
//
// P(u) = [nx,ny,nz] \dot (u - [ nw*nx, nw*ny, nw*nz ] )
//
function _Pf(u, Np) {
  let Pp = [
    Np[3]*Np[0],
    Np[3]*Np[1],
    Np[3]*Np[2]
  ];

  let _u = [ u[0], u[1], u[2] ];
  let _Np = [ Np[0], Np[1], Np[2] ];

  return plane_f(_u, _Np, Pp);
}

// return "time" value of line to plane intersection
// line(t) = v0 + t v
// plane(u) = Np . ( u - p )
//
// -> Np . ( v0 + t v - p ) = 0
// -> t = ( (Np . p) - (Np . v0) ) / (Np . v)
//
function t_plane_line(Np, p, v0, v) {
  let _eps = 1/(1024*1024*1024);
  let _d = njs.dot(Np,v);
  if (Math.abs(_d) < _eps) { return NaN; }

  let t = (njs.dot(Np,p) - njs.dot(Np,v0)) / _d;
  return t;
}

// evaluate line
//
// V(t) = v0 + t v
//
function Vt( v0, v, t ) {
  return njs.add(v0, njs.mul(t, v));
}

// 3d cross product.
//
function cross3(p,q) {
  let c0 = ((p[1]*q[2]) - (p[2]*q[1])),
      c1 = ((p[2]*q[0]) - (p[0]*q[2])),
      c2 = ((p[0]*q[1]) - (p[1]*q[0]));

  return [c0,c1,c2];
}

// Three points to plane equation,
// returns 4x1 vector,
// where first three elemenst are normal,
// last is distance to plane from origin.
//
function p3toP(p0, p1, p2) {
  let _eps = 1/(1024*1024*1024);

  let p10 = njs.sub(p1,p0);
  let p20 = njs.sub(p2,p0);

  let Np = cross3(p10,p20);

  let Pk = -njs.dot(Np, p0);

  return [ Np[0], Np[1], Np[2], Pk ]
}


// p0 : point on plane
// u  : normal to plane
// ds : frustum scaling factor (default 1)
//
// returns:
//
// { 
//   idir       : if q-plane fully intersects the frustum vectors, holds idir that this happens
//                -1 if none found
//   idir_t     : four vector of time values (positive) that the intersection happens
//                default if idir < 0
//
//   frustum_idir : frustum q-point sits in
//   frustum_t    : 'time' values of q-plane intersection to each frustum vector
//   frustum_v    : 3d vectors of frustum vectors, origin centered ([idir][f_idx][xyz])
//   
//   // WIP
//   frame_t    : frame time (frame edge in frustum order) (source, dest)
//   frame_updated : 1 source/dest frame_t updated
// }
//
// So, here's what I think should happen:
//
// frustum_t[k] frustum_t[k+1] both in (0,1), wndow closed
// frame_updated 1 -> look at frame_t to see if window needs updating
//
function frustum3d_intersection(p_world, q_world, ds) {
  ds = ((typeof ds === "undefined") ? 1 : ds);

  let q = njs.sub(q_world, p_world);

  let s3 = 1/Math.sqrt(3);
  let s3ds = s3*ds;

  let L = ds;
  let _eps = (1.0 / (1024*1024*1024));
  let oppo = [1,0, 3,2, 5,4];

  let _res_t = [
    [-1,-1, 0, 0, 0, 0 ],
    [-1,-1, 0, 0, 0, 0 ],
    [ 0, 0,-1,-1, 0, 0 ],
    [ 0, 0,-1,-1, 0, 0 ],
    [ 0, 0, 0, 0,-1,-1 ],
    [ 0, 0, 0, 0,-1,-1 ]
  ];

  let _frustum_t = [
    [ 0, 0, 0, 0 ],
    [ 0, 0, 0, 0 ],
    [ 0, 0, 0, 0 ],
    [ 0, 0, 0, 0 ],
    [ 0, 0, 0, 0 ],
    [ 0, 0, 0, 0 ]
  ];

  let v_idir = [
    [1,0,0], [-1,0,0],
    [0,1,0], [0,-1,0],
    [0,0,1], [0,0,-1]
  ];

  // idir:
  // 0 : +x, 1 : -x
  // 2 : +y, 3 : -y
  // 4 : +z, 5 : -z
  //
  // ccw order
  //
  let frustum_v = [
    [ [ L, L, L ], [ L,-L, L ], [ L,-L,-L ], [ L, L,-L ] ],
    [ [-L, L, L ], [-L, L,-L ], [-L,-L,-L ], [-L,-L, L ] ],

    [ [ L, L, L ], [ L, L,-L ], [-L, L,-L ], [-L, L, L ] ],
    [ [ L,-L, L ], [-L,-L, L ], [-L,-L,-L ], [ L,-L,-L ] ],

    [ [ L, L, L ], [-L, L, L ], [-L,-L, L ], [ L,-L, L ] ],
    [ [ L, L,-L ], [ L,-L,-L ], [-L,-L,-L ], [-L, L,-L ] ]
  ];

  let v_norm = njs.norm2( [L,L,L] );

  let qn = njs.norm2(q);
  let q2 = njs.norm2Squared(q);
  if (q2 < _eps) { return; }

  let found_idir = -1;

  let n_q = njs.mul(1/qn, q);

  // res_t calculation
  // simple case of where each frustum vector intersects the q-plane
  //
  // n, normal to plane: q / |q|
  // plane(u) = n . (u - q)
  // v(t) = t . v_k  (point on frustum vector, $t \in \mathbb{R}$ parameter)
  // => n . ( t . v_k - q ) = 0
  // => t = ( q . n ) / (n . v _k)
  //      = ( q . (q / |q|) ) / ( (q / |q|) . v_k )
  //      = |q|^2 / (q . v_k)
  //  
  for (idir=0; idir<6; idir++) {
    let fv_count = 0;
    let fv_n = frustum_v[idir].length;

    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v = frustum_v[idir][f_idx];

      let qv = njs.dot(q,v);
      if (Math.abs(qv) < _eps) { continue; }

      let t = q2 / qv;
      _frustum_t[idir][f_idx] = t;
      if (t < 0) { continue; }
      fv_count++;
    }

    if (fv_count < fv_n) { continue; }

    // we've found a positive time intersection for each of the
    // four frustum vectors to the q-plane.
    // Remember the idir we've found it in

    found_idir = idir;

    // now fill out res_t with actual time values for the intersection
    //
    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v = frustum_v[idir][f_idx];

      for (let pn=-1; pn<2; pn+=2) {
        let v_nei = frustum_v[idir][(f_idx+pn+fv_n)%fv_n];
        let win_edge = njs.sub(v_nei, v)

        // plane(u) = n . (u - q)
        // w(t) = w_0 + t w_v
        // => n . ( w_0 + t w_v - q ) = 0
        // => t = [ (n . q) - (n . w_0) ] / (n . w_v)
        //      = [ ((q / |q|) . q) - ((q / |q|) . w_0) ] / ((q / |q|) . w_v)
        //      = [ |q|^2 - (q . w_0) ] / (q . w_v)

        let _d = njs.dot(q, v_nei);
        if (Math.abs(_d) < _eps) { continue; }

        let t_w = ( q2 - njs.dot(q,v) ) / _d;

        let edge_idir = v2idir(win_edge);

        if (njs.dot(n_q, njs.sub(v, q)) < 0) {
          edge_idir = oppo[edge_idir];
        }

        _res_t[idir][edge_idir] = t_w;

      }

    }

  }

  let frustum_idir=-1;
  for (idir = 0; idir < 6; idir++) {
    let part_count = 0;
    let _n = frustum_v[idir].length;
    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v_cur = frustum_v[idir][f_idx];
      let v_nxt = frustum_v[idir][(f_idx+1)%_n];
      let vv = cross3(v_cur, v_nxt);
      if (njs.dot(vv, q) >= 0) { part_count++; }
    }
    if (part_count == _n) { frustum_idir = idir; }
  }

  let frame_d = -1;
  let frame_sd_t = [ [-1,-1], [-1,-1], [-1,-1], [-1,-1] ];
  let frame_sd_updated = [ [0,0], [0,0], [0,0], [0,0] ];
  let frame_sd_side = [ 0,0,0,0 ];
  if (frustum_idir>=0) {
    let idir = frustum_idir;
    let _n = frustum_v[idir].length;
    let Nq = njs.mul( 1 / njs.norm2(q), q );
    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v_cur = frustum_v[idir][f_idx];
      let v_nxt = frustum_v[idir][(f_idx+1)%_n];
      let vv = njs.sub(v_nxt, v_cur);

      let t1 = t_plane_line( Nq, q, v_cur, vv );

      if ((t1 < 0) || (t1 > 1)) { continue; }

      frame_sd_side[f_idx] = plane_f(v_cur, Nq, q);

      if (plane_f(v_cur, Nq, q) > 0) {
        frame_sd_t[f_idx][0] = t1;
        frame_sd_updated[f_idx][0] = 1;
      }
      else {
        frame_sd_t[f_idx][1] = t1;
        frame_sd_updated[f_idx][1] = 1;
      }
    }
  }

  return {
    "idir": found_idir,
    "idir_t": _res_t,
    "frustum_t": _frustum_t,
    "frustum_v": frustum_v,
    "frustum_idir": frustum_idir,
    "frame_t" : frame_sd_t,
    "frame_updated": frame_sd_updated
  };

}



//---
//---
//---

function check_cmp(res, edge) {
  let resE = [];
  let n = res.P.length;

  for (let i=0; i<n; i++) {
    resE.push([]);
    for (let j=0; j<n; j++) {
      resE[i].push(0);
    }
  }

  for (let i=0; i<res.E.length; i++) {
    let a = res.E[i][0];
    let b = res.E[i][1];

    resE[a][b] = 1;
    resE[b][a] = 1;
  }

  let mismatch_count = 0;

  for (let i=0; i<n; i++) {
    for (let j=0; j<n; j++) {
      if (resE[i][j] != edge[i][j]) {

        if (DEBUG_LEVEL > 0) {
          console.log("# mismatch", i, j);
        }
          console.log("# mismatch", i, j, "(res:", resE[i][j], "edge:", edge[i][j], ")");

        mismatch_count++;
      }
    }
  }

  if (mismatch_count > 0) { return false; }
  return true;
}

//----
//----
//----

function print_E(pnt, edge) {

  for (let i=0; i<edge.length; i++) {
    let a = edge[i][0];
    let b = edge[i][1];

    let p = pnt[a];
    let q = pnt[b];

    if (p.length == 2) {
      console.log(p[0], p[1]);
      console.log(q[0], q[1]);
    }
    else {
      console.log(p[0], p[1], p[2]);
      console.log(q[0], q[1], q[2]);
    }
    console.log("\n");
  }
}

function print_E_naive(pnt, edge) {

  for (let i=0; i<edge.length; i++) {
    for (let j=0; j<edge.length; j++) {

      if (edge[i][j] == 0) { continue; }

      let p = pnt[i];
      let q = pnt[j];

      if (p.length == 2) {
        console.log(p[0], p[1]);
        console.log(q[0], q[1]);
      }
      else {
        console.log(p[0], p[1], p[2]);
        console.log(q[0], q[1], q[2]);
      }
      console.log("\n");
    }
  }
}

// small tests for validation
//
function small_2d_tests() {
  let _N = [10,20,30,40,50,100,200,300];

  let _debug_output = true;

  let passed = true;

  for (let ii=0; ii<_N.length; ii++) {
    N = _N[ii];
    pnt = poisson_point(N, 2);

    let res = gen_instance_2d_fence(N, pnt);
    let naive_res = naive_relnei_E(pnt);
    let Echeck = naive_res.A;

    let _cr = check_cmp(res, Echeck);
    let sfx = (_cr ? "ok" : "error");

    if (_debug_output) {
      console.log("# n:", N, ", got:", _cr, "(", sfx, ")");
    }

    if (!_cr) { passed = false; }
  }

  return passed;
}


//-----------
//-----------
//-----------
//-----------
//-----------
//-----------
//-----------
//-----------
//-----------

function naive_relnei_E(pnt) {
  let n = pnt.length;

  let A = [];

  for (let i=0; i<n; i++) {
    A.push([]);
    for (let j=0; j<n; j++) {
      A[i].push( (i==j) ? 0 : 1 );
    }
  }

  for (let i=0; i<pnt.length; i++) {
    for (let j=0; j<pnt.length; j++) {
      if (i==j) { continue; }
      for (let k=0; k<n; k++) {
        if ((i==k) || (j==k)) { continue; }
        if (in_lune(pnt[i], pnt[j], pnt[k])) {
          A[i][j] = 0;
          A[j][i] = 0;
        }
      }
    }

  }

  let E = [];

  for (let i=0; i<n; i++) {
    for (let j=0; j<n; j++) {
      if (A[i][j] > 0.5) { E.push([i,j]); }
    }
  }

  return { "P": pnt, "E": E, "A": A };
}

function oob(p, B) {
  B = ((typeof B === "undefined") ? [[0,0,0], [1,1,1]] : B);

  for (i=0; i<p.length; i++) {
    if (p[i] < B[0][i]) { return true; }
    if (p[i] >= B[1][i]) { return true; }
  }

  return false;
}


// Notes for future me:
//
// For 
//
// When calculating the sweep 
// the R1 octent wedge might intersect the topmost square at i_radius=1
// but subsequent intersections will all be on the rightmost face.
// So gathering points for comparison should use all of i_radius=1
// and then higher i_radius can use only the sides.
//
function grid_sweep_2d(ctx, pnt, i_radius, r_idx) {
  let pi8 = Math.PI / 8;
  let theta = r_idx * pi8;

  let wedge_ends = [
    [ [ 1, 0], [ 1, 1] ],
    [ [ 1, 1], [ 0, 1] ],

    [ [ 0, 1], [-1, 1] ],
    [ [-1, 1], [-1, 0] ],

    [ [-1, 0], [-1,-1] ],
    [ [-1,-1], [ 0,-1] ],

    [ [ 0,-1], [ 1,-1] ],
    [ [ 1,-1], [ 1, 0] ]
  ];

  let face_dir = [
    [ 0, 1],
    [-1, 0],
    [ 0,-1],
    [ 1, 0]
  ];

  let cell_size = ctx.grid_cell_size;

  let ipnt = [
    Math.floor(pnt[0] / cell_size[0]),
    Math.floor(pnt[1] / cell_size[1])
  ];

  let info = {
    "path": [],
    "p": [ pnt[0], pnt[1] ],
    "s": [], "S": [],
    "v": [], "V": [],
    "n": 0
  };

  let grid_bbox = [ [ 0, 0 ], [ ctx.grid_n, ctx.grid_n ] ];

  let q_dir = Math.floor(((r_idx+1)%8)/2);
  let vdir = face_dir[q_dir];

  let s_wedge = njs.mul(i_radius, wedge_ends[r_idx][0]);
  let e_wedge = njs.mul(i_radius, wedge_ends[r_idx][1]);

  let n = Math.abs( (e_wedge[0] - s_wedge[0]) + (e_wedge[1] - s_wedge[1]) );
  n = i_radius;

  let s_ipnt = [ ipnt[0] + s_wedge[0], ipnt[1] + s_wedge[1] ];
  let cur_ipnt = [0,0];

  info.s = s_ipnt;
  info.S = [ s_ipnt[0]*cell_size[0], s_ipnt[1]*cell_size[1] ];
  info.v = vdir;
  info.V = [vdir[0]*cell_size[0], vdir[1]*cell_size[1]];
  info.n = n;

  if (i_radius == 0) {
    info.path.push([s_ipnt[0], s_ipnt[1]]);
    return info;
  }

  let grid_path = [];

  if ((r_idx%2)==1) {
    let cdir = face_dir[ (q_dir+3)%4 ];
    let _p = [ s_ipnt[0] - cdir[0], s_ipnt[1] - cdir[1] ];
    if (!oob(_p, grid_bbox)) {
      grid_path.push( _p );
    }
  }

  for (let i=0; i<=n; i++) {
    let _p = [ s_ipnt[0] + vdir[0]*i, s_ipnt[1] + vdir[1]*i ];
    cur_ipnt = _p;
    if (oob(_p, grid_bbox)) { continue; }
    grid_path.push( _p );
  }

  if ((r_idx%2)==0) {
    let cdir = face_dir[ (q_dir+1)%4 ];
    let _p = [ cur_ipnt[0] + cdir[0], cur_ipnt[1] + cdir[1] ];
    if (!oob(_p, grid_bbox)) {
      grid_path.push( [ cur_ipnt[0] + cdir[0], cur_ipnt[1] + cdir[1] ] );
    }
  }


  info.path = grid_path;
  return info;
}

function debug_sweep() {

  let ctx = {
    "grid_cell_size": [1/8, 1/8],
    "bbox": [[0,0], [1,1]],
    "n": 100
  };

  let grid_s = Math.sqrt(ctx.n);
  let grid_n = Math.ceil(grid_s);

  ctx["grid_s"] = grid_s;
  ctx["grid_n"] = grid_n;
  ctx["ds"] = 1 / grid_n;
  ctx["grid_cell_size"] = [ 1/grid_n, 1/grid_n ];


  let pnt = [ Math.random(), Math.random() ];
  pnt = [.9,0];

  console.log("# anchor point:", pnt, "grid_n:", grid_n, "ds:", ctx.ds);

  for (let ir=0; ir<4; ir++) {
    for (let r_idx=0; r_idx<8; r_idx++) {
      let info = grid_sweep_2d(ctx, pnt, ir, r_idx);

      let f = 1/4;
      let dxy = [ f*Math.random(), f*Math.random() ];

      for (let i=0; i<info.path.length; i++) {
        console.log(info.path[i][0] + dxy[0], info.path[i][1] + dxy[1]);
      }
      console.log("");

    }
  }
}

//debug_sweep();
//process.exit();


function octant_index_2d(p,q) {

  let octant_lookup = [ 4,5,6,7, 0,1,2,3 ];
  let pi4 = Math.PI / 4;

  let dy = q[1] - p[1];
  let dx = q[0] - p[0];

  let theta = Math.atan2(dy,dx) + Math.PI;
  let idx = Math.floor(theta / pi4);
  idx %= 8;

  return octant_lookup[idx];
}

function prof_s(ctx, name) {
  if (!(name in ctx)) {
    ctx[name] = { "s": 0, "e": 0, "c": 0, "t":0 };
  }
  ctx[name].s = Date.now();
  return ctx[name].s
}

function prof_e(ctx, name) {
  ctx[name].e = Date.now();
  ctx[name].t += (ctx[name].e - ctx[name].s);
  ctx[name].c++;
  return ctx[name].e;
}

function prof_avg(ctx, name) {
  if (!(name in ctx)) { return 0; }
  if (ctx[name].c == 0) { return 0; }
  return ctx[name].t / (1000*ctx[name].c);
}

function prof_reset(ctx) {
  for (let name in ctx) {
    ctx[name].s = 0;
    ctx[name].e = 0;
    ctx[name].t = 0;
    ctx[name].c = 0;
  }
  return ctx;
}

function prof_print(ctx) {
  for (let name in ctx) {
    console.log("#", name, (ctx[name].t / ctx[name].c) / 1000, "s", "(", ctx[name].t, "/", ctx[name].c, ")");
  }
}

function perf_experiment() {

  let NREP = 10;
  for (let n=1000; n<20001; n+=1000) {


    for (let rep=0; rep<NREP; rep++) {
      prof_s(prof_ctx, "tot");
      gen_instance_2d(n, [[0,0],[1,1]]);
      prof_e(prof_ctx, "tot");
    }

    console.log("#n:", n);

    console.log(n, prof_avg(prof_ctx, "tot"));

    //prof_print(prof_ctx);
    prof_reset(prof_ctx);
    //console.log("\n#---\n");
  }
  process.exit();

}

function alloc_info_3d(n, B, pnts) {
  pnts = ((typeof _point === "undefined") ? [] : pnts);

  let _eps = 1 / (1024*1024*1024);

  let info = {
    "dim": 3,
    "start": [0,0,0],
    "size": [1,1,1],
    "point": [],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1,-1],
    "bbox": [[0,0,0], [1,1,1]],
    "grid": [],
    "edge": [],
    "P": [],
    "E": []
  };

  let grid_s = Math.cbrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  let grid_size  = [ 1, 1, 1 ];
  let grid_start = [ 0, 0, 0 ];
  let grid_cell_size = [ ds, ds, ds ];


  let s3 = Math.sqrt(3)/2;


  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_cell_size[2] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;

  info.size = grid_size;
  info.start = grid_start;

  // initialize points, creating random ones if ncessary
  //
  for (let i=0; i<n; i++) {
    if (i < pnts.length) {
      info.point.push(pnts[i]);
    }
    else {
      let pnt = [
        Math.random()*grid_size[0] + grid_start[0],
        Math.random()*grid_size[1] + grid_start[1],
        Math.random()*grid_size[2] + grid_start[2]
    ];
      info.point.push(pnt);
    }
    info.point_grid_bp.push([-1,-1,-1]);
    //info.edge.push([]);
  }


  info.P = info.point;

  // init grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
      for (let k=0; k<grid_n; k++) {
        info.grid[i][j].push([]);
      }
    }
  }


  // setup lll grid binning
  //
  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.point[i][0]*grid_n);
    let iy = Math.floor(info.point[i][1]*grid_n);
    let iz = Math.floor(info.point[i][2]*grid_n);
    info.grid[iz][iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy,iz];
  }


  return info;
}

// WIP!!!
//
function lune_network_3d_shrinking_fence(n, B, _point) {
  _point = ((typeof _point === "undefined") ? [] : _point);

  let _eps = 1 / (1024*1024*1024);

  let info = {
    "dim": 3,
    "start": [0,0,0],
    "size": [1,1,1],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1,-1],
    "bbox": [[0,0,0], [1,1,1]],
    "grid": [],
    "P": [],
    "E": []
  };

  let s3 = Math.sqrt(3)/2;

  let v_idir = [
    [1,0,0], [-1,0,0],
    [0,1,0], [0,-1,0],
    [0,0,1], [0,0,-1]
  ];

  // 0 : +x, 1 : -x
  // 2 : +y, 3 : -y
  // 4 : +z, 5 : -z
  //
  let frustum_v = [
    [ [ s3, s3, s3 ], [ s3,-s3, s3 ], [ s3,-s3,-s3 ], [ s3, s3,-s3 ] ],
    [ [-s3, s3, s3 ], [-s3, s3,-s3 ], [-s3,-s3,-s3 ], [-s3,-s3, s3 ] ],

    [ [ s3, s3, s3 ], [ s3, s3,-s3 ], [-s3, s3,-s3 ], [-s3, s3, s3 ] ],
    [ [ s3,-s3, s3 ], [-s3,-s3, s3 ], [-s3,-s3,-s3 ], [ s3,-s3,-s3 ] ],

    [ [ s3, s3, s3 ], [-s3, s3, s3 ], [-s3,-s3, s3 ], [ s3,-s3, s3 ] ],
    [ [ s3, s3,-s3 ], [ s3,-s3,-s3 ], [-s3,-s3,-s3 ], [-s3, s3,-s3 ] ]

  ];


  //DEBUG
  //DEBUG
  let debug_frustum = false;
  if (debug_frustum) {
    for (let idir=0; idir<frustum_v.length; idir++) {
      let fr = frustum_v[idir];
      for (let v_idx=0; v_idx<fr.length; v_idx++) {
        console.log(fr[v_idx][0],fr[v_idx][1],fr[v_idx][2]);
      }
      console.log(fr[0][0],fr[0][1],fr[0][2]);
      console.log("\n");
    }
  }
  //DEBUG
  //DEBUG

  let grid_s = Math.cbrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_cell_size[2] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;

  // init grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
      for (let k=0; k<grid_n; k++) {
        info.grid[i][j].push([]);
      }
    }
  }

  let grid_size  = [ 1, 1, 1 ];
  let grid_start = [ 0,0, 0 ];
  let grid_cell_size = [ ds, ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  // initialize points, creating random ones if ncessary
  //
  for (let i=0; i<n; i++) {
    if (i < _point.length) {
      info.P.push(_point[i]);
    }
    else {
      let pnt = [
        Math.random()*grid_size[0] + grid_start[0],
        Math.random()*grid_size[1] + grid_start[1],
        Math.random()*grid_size[2] + grid_start[2]
      ];
      info.P.push(pnt);
    }
    info.point_grid_bp.push([-1,-1,-1]);
    //info.edge.push([]);
  }

  let P = info.P;
  let E = [];

  // setup lll grid binning
  //
  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.P[i][0]*grid_n);
    let iy = Math.floor(info.P[i][1]*grid_n);
    let iz = Math.floor(info.P[i][2]*grid_n);
    info.grid[iz][iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy,iz];
  }

  //DEBUG
  //DEBUG
  //print grid
  let _debug_grid = false;
  if (_debug_grid) {
    console.log("#grid (grid_n:", grid_n, "ds:", ds,  ")");
    for (let iz=0; iz<info.grid.length; iz++) {
      for (let iy=0; iy<info.grid[iz].length; iy++) {
        for (let ix=0; ix<info.grid[iz][iy].length; ix++) {
          console.log( ix*ds, iy*ds, iz*ds );
          console.log( (ix+1)*ds, iy*ds, iz*ds );
          console.log("\n");

          console.log( ix*ds, iy*ds, iz*ds );
          console.log( ix*ds, (iy+1)*ds, iz*ds );
          console.log("\n");

          console.log( ix*ds, iy*ds, iz*ds );
          console.log( ix*ds, iy*ds, (iz+1)*ds );
          console.log("\n");

        }
      }
    }
  }
  //DEBUG
  //DEBUG

  let _debug = false;

  for (let p_idx = 0; p_idx < info.P.length; p_idx++) {

    //FAILURE DEBUGGING
    //if ((p_idx == 14) || (p_idx==15)){ _debug = true; }
    //else { _debug = false; }

    //let p_idx = 0;
    let p = info.P[p_idx];

    let Wp = [ p[0]*grid_n, p[1]*grid_n, p[2]*grid_n ];
    let ip = Wp.map( Math.floor );

    let p_fence = [
      grid_n-ip[0], ip[0],
      grid_n-ip[1], ip[1],
      grid_n-ip[2], ip[2]
    ];

    let p_near_idir = 1;
    let l0 = Wp[0] - ip[0];
    for (let xyz=0; xyz<3; xyz++) {
      let _l = Wp[xyz] - ip[xyz];
      if (_l < l0) {
        p_near_idir = 2*xyz + 1;
        l0 = _l;
      }
      _l = 1 - (Wp[xyz] - ip[xyz]);
      if (_l < l0) {
        p_near_idir = 2*xyz + 0;
        l0 = _l;
      }
    }
    l0 *= ds;
    let t0 = l0*Math.sqrt(3);

    if (_debug) {
      console.log("# P[", p_idx, "]:", p, "ip:", ip, "Wp:", Wp, "l0:", l0, "(near_idir:", p_near_idir, ")");

      console.log(p[0], p[1], p[2]);
      console.log(p[0] + l0*v_idir[p_near_idir][0], p[1] + l0*v_idir[p_near_idir][1], p[2] + l0*v_idir[p_near_idir][2]);
      console.log("\n");
    }

    let p_f_v = [];
    for (let idx=0; idx<frustum_v.length; idx++) {
      p_f_v.push([]);
      for (let ii=0; ii<frustum_v[idx].length; ii++) {
        //p_f_v[idx].push( njs.add( njs.mul(ds, frustum_v[idx][ii]), p ) );
        p_f_v[idx].push( njs.mul(ds, frustum_v[idx][ii]) );
      }
    }

    if (_debug) {
      //DEBUG
      //DEBUG
      console.log("# local p[", p_idx, "], fence:", p_fence, " frustum:");
      for (let idx=0; idx<p_f_v.length; idx++) {
        for (let ii=0; ii<p_f_v[idx].length; ii++) {
          console.log(p[0], p[1], p[2]);
          //console.log(p_f_v[idx][ii][0], p_f_v[idx][ii][1], p_f_v[idx][ii][2]);
          let _vt = njs.add(p, p_f_v[idx][ii]);
          console.log( _vt[0], _vt[1], _vt[2] );
          console.log("\n");

          console.log(p[0], p[1], p[2]);
          //console.log(p_f_v[idx][ii][0], p_f_v[idx][ii][1], p_f_v[idx][ii][2]);
          _vt = njs.add(p, njs.mul(grid_n, p_f_v[idx][ii]));
          console.log( _vt[0], _vt[1], _vt[2] );
          console.log("\n");
        }
      }

      console.log("# l0 test");
      for (let ii=0; ii<4; ii++) {
        console.log(p[0], p[1], p[2]);
        let _vt = njs.add( p, njs.mul( t0, frustum_v[p_near_idir][ii] ) );
        console.log(_vt[0], _vt[1], _vt[2]);
        console.log("\n");
      }
      //DEBUG
      //DEBUG
    }

    //  points in fence (index)
    //
    let pif_list = [];

    for (let ir = 0; ir < grid_n; ir++) {

      if (_debug) {
        console.log("# ir:", ir, "fence:", p_fence.join(","));
      }


      let fenced_in = true;
      for (let ii=0; ii<p_fence.length; ii++) {
        if (p_fence[ii] >= ir) { fenced_in = false; break; }
      }
      if (fenced_in) {

        if (_debug) {
          console.log("#fenced in (fence:", p_fence, "), breaking");
        }

        break;
      }

      let sweep = grid_sweep_perim_3d(info, info.P[p_idx], ir);

      let sweep_q_idx = [];
      for (let path_idx = 0; path_idx < sweep.path.length; path_idx++) {
        let ixyz = sweep.path[path_idx];
        let grid_bin = info.grid[ixyz[2]][ixyz[1]][ixyz[0]];
        for (let bin_idx = 0; bin_idx < grid_bin.length; bin_idx++) {
          sweep_q_idx.push( grid_bin[bin_idx] );
        }

        //DEBUG
        if (_debug) {
          console.log("# grid_bin", ixyz, ":", grid_bin.join(","));
        }
        //DEBUG

      }

      for (let nei_idx=0; nei_idx < sweep_q_idx.length; nei_idx++) {
        let q_idx = sweep_q_idx[nei_idx];
        if (q_idx == p_idx) { continue; }

        if (_debug) {
          console.log("# adding q_idx:", q_idx, "to pif_list (", pif_list.join(","), ") (p_idx:", p_idx, ")");
        }

        pif_list.push( q_idx );

        let q = info.P[q_idx];

        let dqp = njs.sub(q,p);
        let qp2 = njs.norm2Squared(dqp);

        let t_frustum = [];
        for (let idir=0; idir<p_f_v.length; idir++) {
          t_frustum.push([]);

          let pos_count = 0,
              min_tI = grid_n,
              max_tI = -1;

          let _debug_vidx = -1, _debug_t = -1;

          for (let ii=0; ii<p_f_v[idir].length; ii++) {

            let v = p_f_v[idir][ii];

            let qp_v = njs.dot(dqp,v);

            if ( Math.abs(qp_v) < _eps) {

              if (_debug) {
                console.log("#skipping p_frustum[", idir, "][", ii, "]: ( (q-p).v =", qp_v, ")");
              }

              continue;
            }

            let t = qp2 / qp_v;
            if (t < 0) { continue; }

            let tI = Math.floor(t - t0);

            if (_debug) {
              console.log("##>> F[", idir, "][", ii, "] p:", p, "q:", q, "t:", t, "t0:", t0, "tI:", tI);
            }

            pos_count++;
            if (tI < min_tI) {
              min_tI = tI;
            }

            if (tI > max_tI) {
              max_tI = tI;
              _debug_vidx = ii;
              _debug_t = t;
            }

          }

          if (pos_count == 4) {
            //if (p_fence[idir] > min_tI) {
              //p_fence[idir] = min_tI;

            if (p_fence[idir] > max_tI) {
              p_fence[idir] = max_tI;

              //DEBUG
              //DEBUG
              if (_debug) {
                if (_debug_vidx >= 0) {
                  console.log(q[0], q[1], q[2]);
                  let _v = njs.add(p, njs.mul(_debug_t, p_f_v[idir][_debug_vidx])) ;
                  console.log( _v[0], _v[1], _v[2] );
                  console.log("\n");
                }

                console.log("## UPDATING FENCE (p_idx:", p_idx, "ir:", ir, ")>> pos_count:", pos_count,
                  "idir:", idir, "fence now:", p_fence, "from q[", q_idx, "]:", q);
              }
              //DEBUG
              //DEBUG

            }
          }

        }

      } // END for nei_idx in sweep

    } // END for ir


    // naive relative neighborhood graph detection by considering
    // all points in fence
    //

    if (_debug) {
      console.log("# naive rng for p_idx:", p_idx, "on pif_list:", pif_list.join(","));
    }

    for (let i = 0; i < pif_list.length; i++) {
      let q_idx = pif_list[i];

      let _found = true;
      for (let j = 0; j < pif_list.length; j++) {
        if (i==j) { continue; }

        let u_idx = pif_list[j];

        if (_debug) {
          console.log("# u:", u_idx, "in lune(", p_idx, q_idx,")?", in_lune(P[p_idx], P[q_idx], P[u_idx]));
        }

        if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {

          if (_debug) {
            console.log("# u_idx:", u_idx, "in lune of (", p_idx, ",", q_idx, ")");
          }

          _found = false;
          break;
        }
      }
      if (_found) {

        if (_debug) {
          console.log("# adding (", p_idx, q_idx, ")");
        }

        E.push([p_idx, q_idx]);

      }
      else {
        if (_debug) {
          console.log("# skipping (", p_idx, q_idx, ")");
        }
      }

    }

  }

  info.E = E;

  return { "P": P, "E": E, "info": info };
}




//
function gen_instance_3d_fence(n, B, _point) {
  _point = ((typeof _point === "undefined") ? [] : _point);

  let _eps = 1 / (1024*1024*1024);

  let info = {
    "dim": 3,
    "start": [0,0,0],
    "size": [1,1,1],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1,-1],
    "bbox": [[0,0,0], [1,1,1]],
    "grid": [],
    "P": [],
    "E": []
  };

  let s3 = 1/Math.sqrt(3);

  let v_idir = [
    [1,0,0], [-1,0,0],
    [0,1,0], [0,-1,0],
    [0,0,1], [0,0,-1]
  ];

  // 0 : +x, 1 : -x
  // 2 : +y, 3 : -y
  // 4 : +z, 5 : -z
  //
  let frustum_v = [
    [ [ s3, s3, s3 ], [ s3,-s3, s3 ], [ s3,-s3,-s3 ], [ s3, s3,-s3 ] ],
    [ [-s3, s3, s3 ], [-s3, s3,-s3 ], [-s3,-s3,-s3 ], [-s3,-s3, s3 ] ],

    [ [ s3, s3, s3 ], [ s3, s3,-s3 ], [-s3, s3,-s3 ], [-s3, s3, s3 ] ],
    [ [ s3,-s3, s3 ], [-s3,-s3, s3 ], [-s3,-s3,-s3 ], [ s3,-s3,-s3 ] ],

    [ [ s3, s3, s3 ], [-s3, s3, s3 ], [-s3,-s3, s3 ], [ s3,-s3, s3 ] ],
    [ [ s3, s3,-s3 ], [ s3,-s3,-s3 ], [-s3,-s3,-s3 ], [-s3, s3,-s3 ] ]

  ];


  //DEBUG
  //DEBUG
  let debug_frustum = false;
  if (debug_frustum) {
    for (let idir=0; idir<frustum_v.length; idir++) {
      let fr = frustum_v[idir];
      for (let v_idx=0; v_idx<fr.length; v_idx++) {
        console.log(fr[v_idx][0],fr[v_idx][1],fr[v_idx][2]);
      }
      console.log(fr[0][0],fr[0][1],fr[0][2]);
      console.log("\n");
    }
  }
  //DEBUG
  //DEBUG

  let grid_s = Math.cbrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_cell_size[2] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;

  // init grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
      for (let k=0; k<grid_n; k++) {
        info.grid[i][j].push([]);
      }
    }
  }

  let grid_size  = [ 1, 1, 1 ];
  let grid_start = [ 0,0, 0 ];
  let grid_cell_size = [ ds, ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  // initialize points, creating random ones if ncessary
  //
  for (let i=0; i<n; i++) {
    if (i < _point.length) {
      info.P.push(_point[i]);
    }
    else {
      let pnt = [
        Math.random()*grid_size[0] + grid_start[0],
        Math.random()*grid_size[1] + grid_start[1],
        Math.random()*grid_size[2] + grid_start[2]
      ];
      info.P.push(pnt);
    }
    info.point_grid_bp.push([-1,-1,-1]);
    //info.edge.push([]);
  }

  let P = info.P;
  let E = [];

  // setup lll grid binning
  //
  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.P[i][0]*grid_n);
    let iy = Math.floor(info.P[i][1]*grid_n);
    let iz = Math.floor(info.P[i][2]*grid_n);
    info.grid[iz][iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy,iz];
  }

  //DEBUG
  //DEBUG
  //print grid
  let _debug_grid = false;
  if (_debug_grid) {
    console.log("#grid (grid_n:", grid_n, "ds:", ds,  ")");
    for (let iz=0; iz<info.grid.length; iz++) {
      for (let iy=0; iy<info.grid[iz].length; iy++) {
        for (let ix=0; ix<info.grid[iz][iy].length; ix++) {
          console.log( ix*ds, iy*ds, iz*ds );
          console.log( (ix+1)*ds, iy*ds, iz*ds );
          console.log("\n");

          console.log( ix*ds, iy*ds, iz*ds );
          console.log( ix*ds, (iy+1)*ds, iz*ds );
          console.log("\n");

          console.log( ix*ds, iy*ds, iz*ds );
          console.log( ix*ds, iy*ds, (iz+1)*ds );
          console.log("\n");

        }
      }
    }
  }
  //DEBUG
  //DEBUG

  let _debug = false;

  for (let p_idx = 0; p_idx < info.P.length; p_idx++) {
    let p = info.P[p_idx];

    let Wp = [ p[0]*grid_n, p[1]*grid_n, p[2]*grid_n ];
    let ip = Wp.map( Math.floor );

    let p_fence = [
      grid_n-ip[0], ip[0],
      grid_n-ip[1], ip[1],
      grid_n-ip[2], ip[2]
    ];

    let p_near_idir = 1;
    let l0 = Wp[0] - ip[0];
    for (let xyz=0; xyz<3; xyz++) {
      let _l = Wp[xyz] - ip[xyz];
      if (_l < l0) {
        p_near_idir = 2*xyz + 1;
        l0 = _l;
      }
      _l = 1 - (Wp[xyz] - ip[xyz]);
      if (_l < l0) {
        p_near_idir = 2*xyz + 0;
        l0 = _l;
      }
    }
    l0 *= ds;
    let t0 = l0*Math.sqrt(3);

    if (_debug) {
      console.log("# P[", p_idx, "]:", p, "ip:", ip, "Wp:", Wp, "l0:", l0, "(near_idir:", p_near_idir, ")");

      console.log(p[0], p[1], p[2]);
      console.log(p[0] + l0*v_idir[p_near_idir][0], p[1] + l0*v_idir[p_near_idir][1], p[2] + l0*v_idir[p_near_idir][2]);
      console.log("\n");
    }

    let p_f_v = [];
    for (let idx=0; idx<frustum_v.length; idx++) {
      p_f_v.push([]);
      for (let ii=0; ii<frustum_v[idx].length; ii++) {
        //p_f_v[idx].push( njs.add( njs.mul(ds, frustum_v[idx][ii]), p ) );
        p_f_v[idx].push( njs.mul(ds, frustum_v[idx][ii]) );
      }
    }

    if (_debug) {
      //DEBUG
      //DEBUG
      console.log("# local p[", p_idx, "], fence:", p_fence, " frustum:");
      for (let idx=0; idx<p_f_v.length; idx++) {
        for (let ii=0; ii<p_f_v[idx].length; ii++) {
          console.log(p[0], p[1], p[2]);
          //console.log(p_f_v[idx][ii][0], p_f_v[idx][ii][1], p_f_v[idx][ii][2]);
          let _vt = njs.add(p, p_f_v[idx][ii]);
          console.log( _vt[0], _vt[1], _vt[2] );
          console.log("\n");

          console.log(p[0], p[1], p[2]);
          //console.log(p_f_v[idx][ii][0], p_f_v[idx][ii][1], p_f_v[idx][ii][2]);
          _vt = njs.add(p, njs.mul(grid_n, p_f_v[idx][ii]));
          console.log( _vt[0], _vt[1], _vt[2] );
          console.log("\n");
        }
      }

      console.log("# l0 test");
      for (let ii=0; ii<4; ii++) {
        console.log(p[0], p[1], p[2]);
        let _vt = njs.add( p, njs.mul( t0, frustum_v[p_near_idir][ii] ) );
        console.log(_vt[0], _vt[1], _vt[2]);
        console.log("\n");
      }
      //DEBUG
      //DEBUG
    }

    //  points in fence (index)
    //
    let pif_list = [];

    for (let ir = 0; ir < grid_n; ir++) {

      if (_debug) {
        console.log("# ir:", ir, "fence:", p_fence.join(","));
      }


      let fenced_in = true;
      for (let ii=0; ii<p_fence.length; ii++) {
        if (p_fence[ii] >= ir) { fenced_in = false; break; }
      }
      if (fenced_in) {

        if (_debug) {
          console.log("#fenced in (fence:", p_fence, "), breaking");
        }

        break;
      }

      let sweep = grid_sweep_perim_3d(info, info.P[p_idx], ir);

      let sweep_q_idx = [];
      for (let path_idx = 0; path_idx < sweep.path.length; path_idx++) {
        let ixyz = sweep.path[path_idx];
        let grid_bin = info.grid[ixyz[2]][ixyz[1]][ixyz[0]];
        for (let bin_idx = 0; bin_idx < grid_bin.length; bin_idx++) {
          sweep_q_idx.push( grid_bin[bin_idx] );
        }

        //DEBUG
        if (_debug) {
          console.log("# grid_bin", ixyz, ":", grid_bin.join(","));
        }
        //DEBUG

      }

      for (let nei_idx=0; nei_idx < sweep_q_idx.length; nei_idx++) {
        let q_idx = sweep_q_idx[nei_idx];
        if (q_idx == p_idx) { continue; }

        if (_debug) {
          console.log("# adding q_idx:", q_idx, "to pif_list (", pif_list.join(","), ") (p_idx:", p_idx, ")");
        }

        pif_list.push( q_idx );

        let q = info.P[q_idx];

        let dqp = njs.sub(q,p);
        let qp2 = njs.norm2Squared(dqp);

        let t_frustum = [];
        for (let idir=0; idir<p_f_v.length; idir++) {
          t_frustum.push([]);

          let pos_count = 0,
              min_tI = grid_n,
              max_tI = -1;

          let _debug_vidx = -1, _debug_t = -1;

          for (let ii=0; ii<p_f_v[idir].length; ii++) {

            let v = p_f_v[idir][ii];

            let qp_v = njs.dot(dqp,v);

            if ( Math.abs(qp_v) < _eps) {

              if (_debug) {
                console.log("#skipping p_frustum[", idir, "][", ii, "]: ( (q-p).v =", qp_v, ")");
              }

              continue;
            }

            let t = qp2 / qp_v;
            if (t < 0) { continue; }

            let tI = Math.floor(t - t0);

            if (_debug) {
              console.log("##>> F[", idir, "][", ii, "] p:", p, "q:", q, "t:", t, "t0:", t0, "tI:", tI);
            }

            pos_count++;
            if (tI < min_tI) {
              min_tI = tI;
            }

            if (tI > max_tI) {
              max_tI = tI;
              _debug_vidx = ii;
              _debug_t = t;
            }

          }

          if (pos_count == 4) {
            //if (p_fence[idir] > min_tI) {
              //p_fence[idir] = min_tI;

            if (p_fence[idir] > max_tI) {
              p_fence[idir] = max_tI;

              //DEBUG
              //DEBUG
              if (_debug) {
                if (_debug_vidx >= 0) {
                  console.log(q[0], q[1], q[2]);
                  let _v = njs.add(p, njs.mul(_debug_t, p_f_v[idir][_debug_vidx])) ;
                  console.log( _v[0], _v[1], _v[2] );
                  console.log("\n");
                }

                console.log("## UPDATING FENCE (p_idx:", p_idx, "ir:", ir, ")>> pos_count:", pos_count,
                  "idir:", idir, "fence now:", p_fence, "from q[", q_idx, "]:", q);
              }
              //DEBUG
              //DEBUG

            }
          }

        }

      } // END for nei_idx in sweep

    } // END for ir


    // naive relative neighborhood graph detection by considering
    // all points in fence
    //

    if (_debug) {
      console.log("# naive rng for p_idx:", p_idx, "on pif_list:", pif_list.join(","));
    }

    for (let i = 0; i < pif_list.length; i++) {
      let q_idx = pif_list[i];

      let _found = true;
      for (let j = 0; j < pif_list.length; j++) {
        if (i==j) { continue; }

        let u_idx = pif_list[j];

        if (_debug) {
          console.log("# u:", u_idx, "in lune(", p_idx, q_idx,")?", in_lune(P[p_idx], P[q_idx], P[u_idx]));
        }

        if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {

          if (_debug) {
            console.log("# u_idx:", u_idx, "in lune of (", p_idx, ",", q_idx, ")");
          }

          _found = false;
          break;
        }
      }
      if (_found) {

        if (_debug) {
          console.log("# adding (", p_idx, q_idx, ")");
        }

        E.push([p_idx, q_idx]);

      }
      else {
        if (_debug) {
          console.log("# skipping (", p_idx, q_idx, ")");
        }
      }

    }

  }

  info.E = E;

  return { "P": P, "E": E, "info": info };
}

//perf_experiment();

// this looks like it scales as O(n^2) and I'm not so
// clear as to why. The original algorithm by Katajainen and
// Nevalamen needed to take special consideration for
// points near the grid boundaries, so that's might be
// what's going on.
//
// Some suspicious behavior:
//
// Looking at how many comparisons are done for points
// in grid cell locations, it looks like points near
// [grid_n-1,*] and even more so near [grid_n-1,grid_n-1]
// take more time than others.
//
// I don't know why there would be an assymetry here.
// If points near the edge have issues, it should be
// symmetric.
//
// So this means there's something I really don't understand
// or there's a bug somewhere.
//
// update:
// I was using grid_s to scale points and I've switched it to grid_n.
// This mitigates the issue but there still looks to be an ever
// so slight bias as now in the positive x direction and the
// near 0 direction...
//
// This might be on account of the ordering used to resolve the wedges.
//
// The increasing size perf experiments were jagged, with quadratic
// increases, then jumps down, then creeping back up.
// This is most likely due to artifacts from choosing grid_s as opposed
// to grid_n, where the gap from grid_s and grid_n was causing points
// on the periphery to take more time. When the point count got
// closer to a perfect square, the gap would close, causing the
// decrease in run time.
//
// I still don't understand why points on the edges aren't taking
// more time.
//
//
//
// ok, well, this has a bug...does not properly calculate rng
//
function gen_instance_2d(n, B, _point) {

  _point = ((typeof _point === "undefined") ? [] : _point);

  let info = {
    "dim": 2,
    "start": [0,0],
    "size": [1,1],
    "point": [],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1],
    "bbox": [[0,0], [1,1]],
    "grid": [],
    "edge": [],
    "P": [],
    "E": []
  };

  let grid_s = Math.sqrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;


  //PROFILING
  prof_s(prof_ctx, "init_grid");
  //PROFILING


  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
    }
  }

  let grid_size  = [ 1, 1 ];
  let grid_start = [ 0,0 ];
  let grid_cell_size = [ ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  for (let i=0; i<n; i++) {
    if (i < _point.length) {
      info.P.push(_point[i]);
    }
    else {
      let pnt = [ Math.random()*grid_size[0] + grid_start[0], Math.random()*grid_size[1] + grid_start[1] ];
      info.P.push(pnt);
    }
    info.point_grid_bp.push([-1,-1]);
  }

  //PROFILING
  prof_e(prof_ctx, "init_grid");
  //PROFILING

  //PROFILING
  prof_s(prof_ctx, "setup");
  //PROFILING


  for (let i=0; i<n; i++) {
    //let ix = Math.floor(info.P[i][0]*grid_s);
    //let iy = Math.floor(info.P[i][1]*grid_s);
    let ix = Math.floor(info.P[i][0]*grid_n);
    let iy = Math.floor(info.P[i][1]*grid_n);
    info.grid[iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy];
  }

  //PROFILING
  prof_e(prof_ctx, "setup");
  //PROFILING


  let P = info.P;
  let G = info.grid;

  let E = [];

  console.log("# ds:", ds, "grid_s:", grid_s);


  //PROFILING
  prof_s(prof_ctx, "main_loop");
  //PROFILING

  for (let p_idx=0; p_idx < P.length; p_idx++) {
    let pnt = P[p_idx];

    let wedge_resolved  = [ 0, 0,  0, 0,  0, 0,  0, 0];
    let wedge_nei       = [-2,-2, -2,-2, -2,-2, -2,-2];

    let resolved_count = 0;

    let i_anch_x = Math.floor(pnt[0] / ds);
    let i_anch_y = Math.floor(pnt[1] / ds);

    //STAT
    //
    /*
    let _min_dist = i_anch_x;
    if ( (grid_n-i_anch_x)  < _min_dist ) { _min_dist = (grid_n-i_anch_x); }
    if ( (i_anch_y)         < _min_dist ) { _min_dist = (i_anch_y); }
    if ( (grid_n-i_anch_y)  < _min_dist ) { _min_dist = (grid_n-i_anch_y); }
    if (!(_min_dist in _debug_stat.count)) {
      _debug_stat.count[_min_dist] = 0;
      _debug_stat.val[_min_dist] = 0;
    }
    _debug_stat.count[_min_dist]++;
    */
    let key = i_anch_x.toString() + " " + i_anch_y.toString();
    if (!(key in _debug_stat.count)) {
      _debug_stat.count[key] = 0;
      _debug_stat.val[key] = 0;
    }
    _debug_stat.count[key]++;
    //
    //STAT

    //PROFILING
    prof_s(prof_ctx, "main_loop.0");
    //PROFILING


    for (let g_idx=0; g_idx < G[i_anch_y][i_anch_x].length; g_idx++) {
      let q_idx = G[i_anch_y][i_anch_x][g_idx];
      if (q_idx == p_idx) { continue; }

      let r_idx = octant_index_2d( P[p_idx], P[q_idx] );
      if (wedge_nei[r_idx] < 0) {
        wedge_nei[r_idx] = q_idx;
        resolved_count++;
        continue;
      }

      let d_prv = njs.norm2( njs.sub(P[p_idx], P[wedge_nei[r_idx]]) );
      let d_cur = njs.norm2( njs.sub(P[p_idx], P[q_idx]) );

      if (d_cur < d_prv) {
        wedge_nei[r_idx] = q_idx;
      }
    }

    //PROFILING
    prof_e(prof_ctx, "main_loop.0");
    //PROFILING

    //PROFILING
    prof_s(prof_ctx, "main_loop.1");
    //PROFILING

    if (resolved_count < 8) {

      // just easier to enumerate all points in initial 3x3 grid region
      //
      for (let idx_y=-1; idx_y<2; idx_y++) {
        for (let idx_x=-1; idx_x<2; idx_x++) {

          let ix = i_anch_x + idx_x;
          let iy = i_anch_y + idx_y;

          if ((ix < 0) || (ix >= grid_n) ||
              (iy < 0) || (iy >= grid_n)) { continue; }

          for (let g_idx=0; g_idx < G[iy][ix].length; g_idx++) {
            let q_idx = G[iy][ix][g_idx];
            if (q_idx == p_idx) { continue; }

            let r_idx = octant_index_2d( P[p_idx], P[q_idx] );
            if (wedge_nei[r_idx] < 0) {
              wedge_nei[r_idx] = q_idx;
              resolved_count++;
              continue;
            }

            let d_prv = njs.norm2( njs.sub(P[p_idx], P[wedge_nei[r_idx]]) );
            let d_cur = njs.norm2( njs.sub(P[p_idx], P[q_idx]) );

            if (d_cur < d_prv) {
              wedge_nei[r_idx] = q_idx;
            }
          }

        }

      }

    }

    //PROFILING
    prof_e(prof_ctx, "main_loop.1");
    //PROFILING

    //console.log("# p_idx:", p_idx, P[p_idx], "wedge_nei:", wedge_nei);

    //PROFILING
    prof_s(prof_ctx, "main_loop.2");
    //PROFILING

    if (resolved_count < 8) {

      for (let r_idx=0; r_idx<8; r_idx++) {

        if (resolved_count==8) { break; }

        if (wedge_nei[r_idx] == -1) { continue; }
        if (wedge_nei[r_idx] >= 0) { continue; }

        for (let ir=2; ir<grid_n; ir++) {
          if (resolved_count==8) { break; }

          let wedge_info = grid_sweep_2d(info, P[p_idx], ir, r_idx);

          //STAT
          //
          //_debug_stat.val[_min_dist]++;
          _debug_stat.val[key]++;
          //
          //STAT

          //console.log(">> ir:", ir, "wedge_info:", wedge_info);

          if (wedge_info.path.length == 0) {
            if (wedge_nei[r_idx] == -2) {
              wedge_nei[r_idx] = -1;
              resolved_count++;
            }
            break;
          }

          for (let w_idx=0; w_idx < wedge_info.path.length; w_idx++) {
            let grid_point = wedge_info.path[w_idx];
            let ix = grid_point[0];
            let iy = grid_point[1];

            //console.log("#>>>>", ix, iy, "grid_n:", grid_n);

            for (let g_idx=0; g_idx<G[iy][ix].length; g_idx++) {
              let q_idx = G[iy][ix][g_idx];
              let qr_idx = octant_index_2d( P[p_idx], P[q_idx] );
              if (qr_idx != r_idx) { continue; }

              if (wedge_nei[r_idx] < 0) {
                wedge_nei[r_idx] = q_idx;
                resolved_count++;
                continue;
              }

              let d_prv = njs.norm2( njs.sub(P[p_idx], P[wedge_nei[r_idx]]) );
              let d_cur = njs.norm2( njs.sub(P[p_idx], P[q_idx]) );

              if (d_cur < d_prv) {
                wedge_nei[r_idx] = q_idx;
              }

            }

          }

          // we've reached the end and haven't found any points in this wedge
          //
          if ((ir == (grid_n-1)) &&
              (wedge_nei[r_idx] == -2)) {
            resolved_count++;
            wedge_nei[r_idx] = -1;
          }

        }
      }

    }

    //PROFILING
    prof_e(prof_ctx, "main_loop.2");
    //PROFILING


    // DEBUG PRINT
    //
    let _debug_octant = false;
    if (_debug_octant) {
      console.log("# p_idx:", p_idx, "pnt:", pnt, "resolved:", resolved_count);

      let _f = 1/64;
      let _dxy = [ _f*Math.random(), _f*Math.random() ];

      for (let r_idx=0; r_idx<8; r_idx++) {
        let q_idx = wedge_nei[r_idx];
        if (q_idx == -2) {
          console.log("### ERROR!!!", "p_idx:", p_idx, P[p_idx], "q_idx:", q_idx, "r_idx:", r_idx);
          return;
        }
        if (q_idx == -1) { continue; }
        console.log(P[p_idx][0] + _dxy[0], P[p_idx][1] + _dxy[1]);
        console.log(P[q_idx][0] + _dxy[0], P[q_idx][1] + _dxy[1]);
        console.log("");
      }

    }

    //PROFILING
    prof_s(prof_ctx, "main_loop.lune_test");
    //PROFILING

    // now that we have fully resolved occupancy of the wedge regions,
    // we go through and test for rng connections
    //
    let nei_pnt = [],
        nei_idx = [];
    for (let r_idx=0; r_idx<8; r_idx++) {
      if (wedge_nei[r_idx] < 0) { continue; }

      let q_idx = wedge_nei[r_idx];

      nei_idx.push(q_idx);
      nei_pnt.push( P[q_idx] );
    }

    for (let i=0; i<nei_pnt.length; i++) {
      let _found = true;
      for (let j=0; j<nei_pnt.length; j++) {
        if (i==j) { continue; }
        if (in_lune(P[p_idx], nei_pnt[i], nei_pnt[j])) {
          _found = false;
          break;
        }
      }
      if (_found) {
        E.push([p_idx, nei_idx[i]]);
      }
    }

    //PROFILING
    prof_e(prof_ctx, "main_loop.lune_test");
    //PROFILING



  }

  //PROFILING
  prof_e(prof_ctx, "main_loop");
  //PROFILING


  let _debug_edge = false;
  if (_debug_edge) {
    for (let e_idx=0; e_idx<E.length; e_idx++) {
      let p_idx = E[e_idx][0];
      let q_idx = E[e_idx][1];

      let p = P[p_idx];
      let q = P[q_idx];

      let _f = 1/64;
      _f = 0;
      let _dxy = [_f*Math.random(), _f*Math.random()];

      console.log(p[0] + _dxy[0], p[1] + _dxy[1]);
      console.log(q[0] + _dxy[0], q[1] + _dxy[1]);
      console.log("");
    }
  }

  //info.P = info.point;
  info.E = E;

  let _debug_grid = false;
  if (_debug_grid) {
    for (let iy=0; iy<grid_n; iy++) {

      for (let ix=0; ix<grid_n; ix++) {

        console.log( ix*ds, iy*ds);
        console.log( (ix+1)*ds, iy*ds);
        console.log("");
        console.log( ix*ds, iy*ds);
        console.log( ix*ds, (iy+1)*ds);
        console.log("");

        for (let ii=0; ii<info.grid[iy][ix].length; ii++) {
          let idx = info.grid[iy][ix][ii];
          let p = info.P[idx];
          console.log("#", ii, ix, iy, "-> [", p[0],p[1] ,"] {", idx, "}");
        }
      }
    }
  }

  return info;
}

function grid_sweep_perim_3d(ctx, pnt, ir) {
  let info = {
    "path": []
  };


  let grid_n = ctx.grid_n;
  let cell_size = ctx.grid_cell_size;
  let cell_offset = [0,0,0];

  let _grid_bbox = [[0,0,0], [grid_n, grid_n, grid_n]];

  let ipnt = [
    Math.floor(pnt[0] / cell_size[0]),
    Math.floor(pnt[1] / cell_size[1]),
    Math.floor(pnt[2] / cell_size[2])
  ];

  let mxyz = [ ipnt[0] - ir, ipnt[1] - ir, ipnt[2] - ir ];
  let Mxyz = [ ipnt[0] + ir+1, ipnt[1] + ir+1, ipnt[2] + ir+1 ];

  //console.log("#", mxyz, Mxyz);

  for (let ix=mxyz[0]; ix<Mxyz[0]; ix++) {
    for (let iy=mxyz[1]; iy<Mxyz[1]; iy++) {
      if (!oob([ix,iy,mxyz[2]], _grid_bbox)) {

        //console.log("## xy-:", ix, iy, mxyz[2]);

        info.path.push([ix,iy,mxyz[2]]);
      }
      if (mxyz[2] == (Mxyz[2]-1)) { continue; }
      if (!oob([ix,iy,Mxyz[2]-1], _grid_bbox)) {

        //console.log("## xy+:", ix, iy, Mxyz[2]-1);

        info.path.push([ix,iy,Mxyz[2]-1]);
      }
    }
  }

  //console.log("#\n#");

  for (let iy=mxyz[1]; iy<Mxyz[1]; iy++) {
    for (let iz=(mxyz[2]+1); iz<(Mxyz[2]-1); iz++) {
      if (!oob([mxyz[0], iy, iz], _grid_bbox)) {

        //console.log("## yz-:", mxyz[0], iy, iz);

        info.path.push([mxyz[0], iy, iz]);
      }
      if (mxyz[0] == (Mxyz[0]-1)) { continue; }
      if (!oob([Mxyz[0]-1, iy, iz], _grid_bbox)) {

        //console.log("## yz+:", Mxyz[0]-1, iy, iz);

        info.path.push([Mxyz[0]-1, iy, iz]);
      }
    }
  }

  //console.log("#\n#");

  for (let ix=(mxyz[0]+1); ix<(Mxyz[0]-1); ix++) {
    for (let iz=(mxyz[2]+1); iz<(Mxyz[2]-1); iz++) {
      if (!oob([ix, mxyz[1], iz], _grid_bbox)) {

        //console.log("## xz-:", ix, mxyz[1], iz);

        info.path.push([ix, mxyz[1], iz]);
      }
      if (mxyz[1] == (Mxyz[1]-1)) { continue; }
      if (!oob([ix, Mxyz[1]-1, iz], _grid_bbox)) {

        //console.log("## xz+:", ix, Mxyz[1]-1, iz);

        info.path.push([ix, Mxyz[1]-1, iz]);
      }
    }
  }

  //console.log("#\n#");

  return info;
}

//DEBUG
//DEBUG
function test_grid_sweep_perim_3d() {
  let B = [[0,0,0],[1,1,1]];

  let info = alloc_info_3d(150, B);

  let sweep = {};

  //print_point(info.P, 1);

  /*
  sweep = grid_sweep_perim_3d(info, info.P[0], 0);
  console.log("#P[0]:", info.P[0]);
  //console.log("#sweep path ir:0", sweep.path);

  for (let ii=0; ii<sweep.path.length; ii++) {
    console.log( sweep.path[ii][0], sweep.path[ii][1], sweep.path[ii][2]);
  }
  console.log("\n");
  */

  sweep = grid_sweep_perim_3d(info, info.P[0], 1);
  console.log("#P[0]:", info.P[0]);
  //console.log("#sweep path ir:1", sweep.path);
  for (let ii=0; ii<sweep.path.length; ii++) {
    console.log( sweep.path[ii][0], sweep.path[ii][1], sweep.path[ii][2]);
  }
  console.log("\n");

  return;

  sweep = grid_sweep_perim_3d(info, info.P[0], 2);
  console.log("#P[0]:", info.P[0]);
  //console.log("#sweep path ir:2", sweep.path);

  for (let ii=0; ii<sweep.path.length; ii++) {
    console.log( sweep.path[ii][0], sweep.path[ii][1], sweep.path[ii][2]);
  }
  console.log("\n");


  process.exit();
}
//test_grid_sweep_perim_3d();
//process.exit();
//DEBUG
//DEBUG


function grid_sweep_perim_2d(ctx, pnt, ir) {
  let face_dir = [
    [ 0, 1],
    [-1, 0],
    [ 0,-1],
    [ 1, 0]
  ];

  let cell_size = ctx.grid_cell_size;

  let cell_offset = [0,0];

  let ipnt = [
    Math.floor(pnt[0] / cell_size[0]),
    Math.floor(pnt[1] / cell_size[1])
  ];

  // fence is r, u, l, d
  //   [p0 , v] (p0 + v(t))
  //
  // in world coordinates
  //
  let info = {
    "path": [],
    "p": [ pnt[0], pnt[1] ],
    "i_p": [ ipnt[0], ipnt[1] ],
    "fence" : [
      [[0,0], [0,0]],
      [[0,0], [0,0]],
      [[0,0], [0,0]],
      [[0,0], [0,0]]
    ],
    "perim_bbox": [],
    "n": 0
  };

  let grid_bbox = [ [ 0, 0 ], [ ctx.grid_n, ctx.grid_n ] ];


  let perim_bbox = [
    [ ipnt[0] - ir, ipnt[1] - ir ],
    [ ipnt[0] + ir + 1, ipnt[1] + ir + 1]
  ];

  let virt_bbox = [
    [ perim_bbox[0][0], perim_bbox[0][1] ],
    [ perim_bbox[1][0], perim_bbox[1][1] ]
  ];

  if (perim_bbox[0][0] < grid_bbox[0][0]) { perim_bbox[0][0] = grid_bbox[0][0]; }
  if (perim_bbox[0][1] < grid_bbox[0][1]) { perim_bbox[0][1] = grid_bbox[0][1]; }

  if (perim_bbox[1][0] > grid_bbox[1][0]) { perim_bbox[1][0] = grid_bbox[1][0]; }
  if (perim_bbox[1][1] > grid_bbox[1][1]) { perim_bbox[1][1] = grid_bbox[1][1]; }

  info.perim_bbox = perim_bbox;

  // right fence, lower right start, move up
  //
  info.fence[0][0][0] = perim_bbox[1][0]*cell_size[0] + cell_offset[0];
  info.fence[0][0][1] = perim_bbox[0][1]*cell_size[1] + cell_offset[1];

  info.fence[0][1][0] = 0;
  info.fence[0][1][1] = (perim_bbox[1][1] - perim_bbox[0][1])*cell_size[1] ;

  // top fence, upper right, move left
  //
  info.fence[1][0][0] = perim_bbox[1][0]*cell_size[0] + cell_offset[0];
  info.fence[1][0][1] = perim_bbox[1][1]*cell_size[1] + cell_offset[1];

  info.fence[1][1][0] = (perim_bbox[0][0] - perim_bbox[1][0])*cell_size[0];
  info.fence[1][1][1] = 0;

  // left fence, upper left, move down
  //
  info.fence[2][0][0] = perim_bbox[0][0]*cell_size[0] + cell_offset[0];
  info.fence[2][0][1] = perim_bbox[1][1]*cell_size[1] + cell_offset[1];

  info.fence[2][1][0] = 0;
  info.fence[2][1][1] = (perim_bbox[0][1] - perim_bbox[1][1])*cell_size[1];

  // bottom fence, lower left, move right
  //
  info.fence[3][0][0] = perim_bbox[0][0]*cell_size[0] + cell_offset[0];
  info.fence[3][0][1] = perim_bbox[0][1]*cell_size[1] + cell_offset[1];

  info.fence[3][1][0] = (perim_bbox[1][0] - perim_bbox[0][0])*cell_size[0];
  info.fence[3][1][1] = 0;


  //console.log("#bbox: [[", perim_bbox[0][0], perim_bbox[0][1], "],[", perim_bbox[1][0], perim_bbox[1][1], "]], ipnt:", ipnt);

  let use_fence_perim = false;

  if (use_fence_perim) {

    let dx = perim_bbox[1][0] - perim_bbox[0][0],
        dy = perim_bbox[1][1] - perim_bbox[0][1];

    let ix = perim_bbox[0][0],
        iy = perim_bbox[0][1];
    for (ix=perim_bbox[0][0]; ix<perim_bbox[1][0]; ix++) {
      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }
      info.path.push([ix,iy]);
    }

    ix = perim_bbox[1][0]-1;
    for (iy=perim_bbox[0][1]; iy<perim_bbox[1][1]; iy++) {
      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }
      info.path.push([ix,iy]);
    }

    iy = perim_bbox[1][1]-1;
    for (ix=(perim_bbox[1][0]-1); ix>=perim_bbox[0][0]; ix--) {
      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }
      info.path.push([ix,iy]);
    }

    ix = perim_bbox[0][0];
    for (iy=(perim_bbox[1][1]-1); iy>=perim_bbox[0][1]; iy--) {
      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }
      info.path.push([ix,iy]);
    }
  }
  else {

    let dx = virt_bbox[1][0] - virt_bbox[0][0],
        dy = virt_bbox[1][1] - virt_bbox[0][1];

    let ix = virt_bbox[0][0],
        iy = virt_bbox[0][1];
    for (ix=virt_bbox[0][0]; ix<virt_bbox[1][0]; ix++) {
      if ((ix < grid_bbox[0][0]) || (ix >= grid_bbox[1][0]) ||
          (iy < grid_bbox[0][1]) || (iy >= grid_bbox[1][1])) {
        continue;
      }

      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }

      info.path.push([ix,iy]);
    }

    ix = virt_bbox[1][0]-1;
    for (iy=virt_bbox[0][1]; iy<virt_bbox[1][1]; iy++) {
      if ((ix < grid_bbox[0][0]) || (ix >= grid_bbox[1][0]) ||
          (iy < grid_bbox[0][1]) || (iy >= grid_bbox[1][1])) {
        continue;
      }

      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }

      info.path.push([ix,iy]);
    }

    iy = virt_bbox[1][1]-1;
    for (ix=(virt_bbox[1][0]-1); ix>virt_bbox[0][0]; ix--) {
      if ((ix < grid_bbox[0][0]) || (ix >= grid_bbox[1][0]) ||
          (iy < grid_bbox[0][1]) || (iy >= grid_bbox[1][1])) {
        continue;
      }

      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }

      info.path.push([ix,iy]);
    }

    ix = virt_bbox[0][0];
    for (iy=(virt_bbox[1][1]-1); iy>virt_bbox[0][1]; iy--) {
      if ((ix < grid_bbox[0][0]) || (ix >= grid_bbox[1][0]) ||
          (iy < grid_bbox[0][1]) || (iy >= grid_bbox[1][1])) {
        continue;
      }

      let idx = info.path.length-1;
      if ((idx >= 0) &&
          (info.path[idx][0] == ix) &&
          (info.path[idx][1] == iy)) { continue; }

      info.path.push([ix,iy]);
    }

  }

  return info;
}

function test_grid_sweep_perim_2d() {
  let n = 107;
  let m = Math.ceil(Math.sqrt(n));
  let ctx = {
    "grid_cell_size": [ 1/m, 1/m ],
    "grid_n": m,
    "n" : n
  };

  let p = [ Math.random(), Math.random() ];
  //p = [ 0.7564873100990699, 0.6250385838518433 ];
  //p = [ 0.99, 0.99 ];
  //p = [ 0, 0.99];
  //p = [ .99, 0];
  //p = [0,0];

  let ir = 1;
  let info = grid_sweep_perim_2d(ctx, p, ir);

  console.log(info.i_p[0], info.i_p[1], "\n\n");

  console.log("# n:", ctx.n, "grid_n:", ctx.grid_n, "p:", p, "ir:", ir);
  for (let i=0; i<info.path.length; i++) {
    console.log(info.path[i][0], info.path[i][1]);
  }

  console.log("\n\n");
  console.log(p[0], p[1]);
  console.log("\n");
  for (let i=0; i<info.fence.length; i++) {
    let dxy = [ Math.random()/64, Math.random()/64 ];
    dxy = [0,0];
    let pv = info.fence[i];
    console.log(pv[0][0] + dxy[0], pv[0][1] + dxy[1] );
    console.log(pv[0][0] + pv[1][0] + dxy[0], pv[0][1] + pv[1][1] + dxy[1] );
    console.log("\n");
  }

}

function print_fence(fence) {
  for (let i=0; i<fence.length; i++) {
    let dxy = [ Math.random()/64, Math.random()/64 ];
    dxy = [0,0];
    let pv = fence[i];
    console.log(pv[0][0] + dxy[0], pv[0][1] + dxy[1] );
    console.log(pv[0][0] + pv[1][0] + dxy[0], pv[0][1] + pv[1][1] + dxy[1] );
    console.log("\n");
  }

}

// find bounding box of lune created by p,q
function lune_bbox(p,q) {

  let r = njs.sub(p,q);

  for (let knockout_xyz=2; knockout_xyz>=0; knockout--) {

    let u = [0,0];
    let v = [0,0];

    let idx = 0;
    for (let xyz=0; xyz<3; xyz++) {
      if (xyz == knockout_xyz) { continue; }
      if (xyz >= p.length) { continue; }
      if (xyz >= q.length) { continue; }

      u[idx] = p[xyz];
      v[idx] = q[xyz];
      idx++;
    }

    let dvu = njs.sub(v,u);

    let theta = Math.atan2(dvu[1], dvu[0]);


  }


}

//
//
function lune_network_shrinking_fence_2d(n, _point) {
  _point = ((typeof _point === "undefined") ? [] : _point);

  let _eps = 1 / (1024*1024*1024);

  let info = {
    "dim": 2,
    "start": [0,0],
    "size": [1,1],
    "point": [],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1],
    "bbox": [[0,0], [1,1]],
    "grid": [],
    "edge": []
  };

  let pi4 = Math.PI/4;
  let s2 = Math.sqrt(2)/2;

  let grid_s = Math.sqrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;
  let s2ds = s2*ds;

  let v_idir = [
    [1,0], [-1,0],
    [0,1], [0,-1]
  ];

  let frustum_v = [
    [ [ s2ds, s2ds ], [ s2ds, -s2ds ] ],
    [ [-s2ds, s2ds ], [-s2ds, -s2ds ] ],

    [ [ s2ds, s2ds ], [-s2ds,  s2ds ] ],
    [ [ s2ds,-s2ds ], [-s2ds, -s2ds ] ]
  ];

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;

  if (DEBUG_LEVEL > 0) {
    console.log("# grid_n:", info.grid_n, "cell_size:", info.grid_cell_size, "grid_s:", info.grid_s, "ds:", ds);
  }

  // alloc grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
    }
  }

  let grid_size  = [ 1, 1 ];
  let grid_start = [ 0,0 ];
  let grid_cell_size = [ ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  // alloc and create random points
  //
  for (let i=0; i<n; i++) {
    let pnt = [0,0];
    if ( i < _point.length ) { pnt = _point[i]; }
    else {
      pnt = [ _rnd()*grid_size[0] + grid_start[0], _rnd()*grid_size[1] + grid_start[1] ];
    }
    info.point.push(pnt);
    info.point_grid_bp.push([-1,-1]);
  }

  // push points into grid, linear linked list/array for dups
  //
  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.point[i][0]*grid_n);
    let iy = Math.floor(info.point[i][1]*grid_n);
    info.grid[iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy];
  }


  if (DEBUG_LEVEL > 1) {
    print_point(info.point, true);
  }

  let P = info.point;
  let E = [];

  let octant2idir = [ 0, 2,2, 1,1, 3,3, 0 ];

  for (let p_idx = 0; p_idx < info.point.length; p_idx++) {

    let p = info.point[p_idx];

    let Wp = [ grid_n*p[0], grid_n*p[1] ];
    let ip = Wp.map( Math.floor );
    let p_fence = [
      grid_n - ip[0], ip[0],
      grid_n - ip[1], ip[1]
    ];

    // +x [+y,-y], -x [+y,-y],
    // +y [+x,-x], -y [+x,-x]
    //
    let p_fence_win = [
      [ 0, 1 ], [ 0, 1 ],
      [ 0, 1 ], [ 0, 1 ]
    ];

    // points in fence
    //
    let pif_list = [];


    if (DEBUG_LEVEL > 1) {
      console.log("# p:", p, "(idx:", p_idx, ")");
      console.log("# init p_fence_idx:", p_fence_idx, "(ip:", ip, ")");
      console.log("# p (point[0])");
      console.log(p[0],p[1], "\n");
    }


    // calculate l0, the distance to learest grid bin
    // and t0, the time it takes for a point along
    // the frustum line to hit the initial grid boundary
    //

    let l0 = Wp[0] - ip[0];
    for (let xy=0; xy<2; xy++) {
      let _l = Wp[xy] - ip[xy];
      if (_l < l0) { l0 = _l; }

      _l = 1.0 - (Wp[xy] - ip[xy]);
      if (_l < l0) { l0 = _l; }
    }
    l0 *= ds;
    let t0 = l0*Math.sqrt(2);

    for (let ir=0; ir<info.grid_n; ir++) {

      if (DEBUG_LEVEL > 1) {
        console.log("#ir:", ir);
      }

      // fence within current fence radius
      //
      let fenced_in = true;
      for (let idir=0; idir<4; idir++) {
        if (ir < p_fence_idx[idir]) { fenced_in = false; break; }
      }

      if (fenced_in) {
        if (DEBUG_LEVEL > 2) {
          console.log("#end search (0) (ir", ir, ", fence:", p_fence_idx , ", fence_wind:", p_fence_win, ")");
        }
        break;
      }

      // window fence has completely closed
      //
      fenced_in = true;
      for (let idir=0; idir<4; idir++) {
        if (p_fence_win[idir][0] < p_fence_win[idir][1]) {
          fenced_in = false;
          break;
        }
      }

      if (fenced_in) {
        if (DEBUG_LEVEL > 2) {
          console.log("#end search (1) (ir", ir, ", fence:", p_fence_idx ,", fence_win:", p_fence_win, ")");
        }
        break;
      }

      let sweep = grid_sweep_perim_2d(info, p, ir);

      if (DEBUG_LEVEL > 2) {
        console.log("# sweep.path:", sweep.path.join(";"));
      }

      let sweep_q_idx = [];
      for (let path_idx=0; path_idx < sweep.path.length; path_idx++) {

        let ix = sweep.path[path_idx][0];
        let iy = sweep.path[path_idx][1];
        let grid_bin = info.grid[iy][ix];
        for (let bin_idx=0; bin_idx < grid_bin.length; bin_idx++) {
          sweep_q_idx.push( grid_bin[bin_idx] );
        }

        for (let nei_idx=0; nei_idx < sweep_q_idx.length; nei_idx++) {

          let q_idx = sweep_q_idx[nei_idx];
          if (q_idx == p_idx) { continue; }

          pif_list.push( q_idx );

          let q = info.P[q_idx];

          let sqp = njs.sub(q,p);
          let qp2 = njs.norm2Squared(dqp);

          for (let idir=0; idir<frustum_v.length; idir++) {

            let pos_count = 0,
                min_tI = grid_n,
                max_tI = -1;

            // count the number of frustum intersections,
            //
            for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
              let v = frustum_v[idir][f_idx];
              let qp_v = njs.dot(sqp,v);

              // perpendicular, ignore
              //
              if ( Math.abs(qp_v) < _eps ) { continue; }

              let t = qp2 / qp_v;

              // in the other direction, skip
              //
              if (t < 0) { continue; }

              let tI = Math.floor(t - t0);

              pos_count++;
              if (tI < min_tI) { min_tI = tI; }
              if (tI > max_tI) { max_tI = tI; }

            }

            // we're in the frustum and we have a new secure point
            // for it
            //
            if ((pos_count == 2) &&
                (p_fence[idir] > max_tI)) {
              p_fence[idir] = max_tI;
            }

          }

        }

      }
    }

    if (DEBUG_LEVEL > 2) {
      console.log("#pif_list.length:", pif_list.length);
    }

    // do a naive RNG centered on p
    //
    for (let i=0; i<pif_list.length; i++) {

      let q_idx = pif_list[i];

      let _found = true;
      for (let j=0; j<pif_list.length; j++) {
        if (i==j) { continue; }

        let u_idx = pif_list[j];

        if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {
          _found = false;
          break;
        }
      }
      if (_found) {
        E.push([p_idx, q_idx]);
      }
    }


  // for p_idx
  }

  return { "P": P, "E": E };
}

function gen_instance_2d_fence(n, _point) {
  _point = ((typeof _point === "undefined") ? [] : _point);

  let _eps = 1 / (1024*1024*1024);

  let info = {
    "dim": 2,
    "start": [0,0],
    "size": [1,1],
    "point": [],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1],
    "bbox": [[0,0], [1,1]],
    "grid": [],
    "edge": [],

    "stat": {
      "q_count": 0,
      "grid_count": 0,
      "shell_count": 0
    }
  };

  let grid_s = Math.sqrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;

  if (DEBUG_LEVEL > 0) {
    console.log("# grid_n:", info.grid_n, "cell_size:", info.grid_cell_size, "grid_s:", info.grid_s, "ds:", ds);
  }


  // alloc grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
    }
  }

  let grid_size  = [ 1, 1 ];
  let grid_start = [ 0,0 ];
  let grid_cell_size = [ ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  // alloc and create random points
  //
  for (let i=0; i<n; i++) {
    let pnt = [0,0];
    if ( i < _point.length ) { pnt = _point[i]; }
    else {
      pnt = [ _rnd()*grid_size[0] + grid_start[0], _rnd()*grid_size[1] + grid_start[1] ];
    }
    info.point.push(pnt);
    info.point_grid_bp.push([-1,-1]);
  }

  // push points into grid, linear linked list/array for dups
  //
  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.point[i][0]*grid_n);
    let iy = Math.floor(info.point[i][1]*grid_n);
    info.grid[iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy];
  }

  // octant to opposite frustum vector
  //
  //  1,4            2,7
  //   ._______________.
  //   |\      |      /|
  //   |  \  2 | 1  /  |
  //   |    \  |  /    |
  //   |  3   \|/   0  |
  //   |-------|-------|
  //   |  4   /|\   7  |    ^
  //   |    /  |  \    |    |
  //   |  /  5 | 6  \  |   +y
  //   |/      |      \|    |
  //   .---------------.    |
  //  3,6             0,5
  //
  //  ---+x--->
  //
  //
  // v_lookup vectors normalized
  //

  let pi4 = Math.PI/4;
  let s2 = Math.sqrt(2)/2;
  let v_lookup = [
    [ s2, -s2 ],
    [-s2,  s2 ],

    [ s2,  s2 ],
    [-s2, -s2 ],

    [-s2,  s2 ],
    [ s2, -s2 ],

    [-s2, -s2 ],
    [ s2,  s2 ]
  ];

  if (DEBUG_LEVEL > 1) {
    print_point(info.point, true);
  }

  let P = info.point;
  let E = [];


  for (let p_idx = 0; p_idx < info.point.length; p_idx++) {

    let p = info.point[p_idx];

    if (DEBUG_LEVEL > 1) {
      console.log("# p:", p, "(idx:", p_idx, ")");
    }

    // points in fence
    //
    let pif_list = [];

    //let p_fence_idx = [ grid_n, grid_n, grid_n, grid_n ];

    let p_grid = [ Math.floor(p[0]/ds), Math.floor(p[1]/ds) ];
    let p_fence_idx = [
      grid_n - p_grid[0], grid_n - p_grid[1],
      p_grid[0], p_grid[1]
    ];

    if (DEBUG_LEVEL > 1) {
      console.log("# init p_fence_idx:", p_fence_idx, "(p_grid:", p_grid, ")");
    }

    let octant2quadrent = [ 0, 1,1, 2,2, 3,3, 0 ];

    if (DEBUG_LEVEL > 1) {
      console.log("# p (point[0])");
      console.log(p[0],p[1], "\n");
    }

    // we're trying to find the intersection of the perpendicular line
    // from p to q as it intersects the edges of the enclosing fence
    // aroudn p
    //
    // l0 min dist between point and grid cell boundary p sits in
    //

    let gpi0 = grid_sweep_perim_2d(info, p, 0);
    let l0 = Math.abs(p[0] - (ds*gpi0.perim_bbox[0][0]));

    let _l = Math.abs(p[0] - (ds*gpi0.perim_bbox[1][0]));
    if ( _l < l0) { l0 = _l; }

    _l = Math.abs(p[1] - (ds*gpi0.perim_bbox[0][1]));
    if ( _l < l0) { l0 = _l; }

    _l = Math.abs(p[1] - (ds*gpi0.perim_bbox[1][1]));
    if ( _l < l0) { l0 = _l; }

    if (DEBUG_LEVEL > 2) {
      console.log("# l0:", l0);
      console.log("# debug print fence:");
      for (let ir=0; ir<5; ir++) {
        let _gpi = grid_sweep_perim_2d(info, p, ir);
        print_fence(_gpi.fence);
      }
    }


    for (let ir=0; ir<info.grid_n; ir++) {

      if (DEBUG_LEVEL > 1) {
        console.log("#ir:", ir, "current fence:", p_fence_idx);
      }

      let grid_perim_info = grid_sweep_perim_2d(info, p, ir);

      if (DEBUG_LEVEL > 2) {
        console.log("#path:", grid_perim_info.path.join(";"));
      }

      let end_search = true;
      for (let i=0; i<4; i++) {
        if (ir <= p_fence_idx[i]) { end_search = false; break; }
      }
      if (end_search) {

        info.stat.shell_count += ir;

        if (DEBUG_LEVEL > 2) {
          console.log("#end search (ir", ir, ", fence:", p_fence_idx ,")");
        }
        break;
      }

      for (let grid_perim_idx=0; grid_perim_idx < grid_perim_info.path.length; grid_perim_idx++) {

        info.stat.grid_count ++;

        let ix = grid_perim_info.path[grid_perim_idx][0];
        let iy = grid_perim_info.path[grid_perim_idx][1];

        let gpi = grid_perim_info;

        if (DEBUG_LEVEL > 2) {
          console.log("# ixy:", ix, iy);
          console.log("# fence:");
          print_fence(grid_perim_info.fence);
          console.log("\n");
        }

        for (let bin_idx=0; bin_idx < info.grid[iy][ix].length; bin_idx++) {

          let q_idx = info.grid[iy][ix][bin_idx];
          if (q_idx == p_idx) { continue; }

          info.stat.q_count++;

          pif_list.push(q_idx);

          let q = info.point[q_idx];

          let Rk = octant_index_2d(p,q);

          let v = v_lookup[Rk];

          // l0 represents the initial size of the enclosing fence around p
          //

          // u = rot_90_ccw(p-q)
          // v = frustum vector opposite to octant occupied by q
          //
          // w_a(t_a) = q + u t_a
          // w_b(t_b) = p + v t_b
          //
          // -> q_x + u_x t_a = p_x + v_x t_b
          //    q_y + u_y t_a = p_y + v_y t_b
          //
          // -> t_b = [ u_y (p_x - q_x) - u_x (p_y - q_y) ] / [ v_y u_x - v_x u_y ]
          //

          let _u = njs.sub(p,q);
          let u = [0,0];

          let _a = Math.atan2(_u[1], _u[0]);
          let _a_idx = Math.floor((_a + Math.PI) / pi4);
          if ((_a_idx % 2) == 0) {
            u = njs.mul( 1/njs.norm2(_u), [-_u[1],  _u[0] ] );
          }
          else {
            u = njs.mul( 1/njs.norm2(_u), [ _u[1], -_u[0] ] );
          }

          let _denom = Math.abs((v[1]*u[0]) - (v[0]*u[1]));
          if (_denom < _eps) { continue; }

          /*
          let t0 = ((v[0]*(q[1] - p[1])) - (v[1]*(q[0] - p[0]))) / ((v[1]*u[0]) - (v[0]*u[1]));

          let sq = njs.add(q, njs.mul(t0, u));
          let t1 = njs.dot(v, njs.sub(sq, p));
          //let tI = Math.ceil( (Math.sqrt(2)*t1) - l0 ) + 1;
          let p1 = njs.sub(p, njs.mul( Math.sqrt(2)*l0, v ));
          let tI = Math.ceil(njs.dot(v, njs.mul( 1/(ds*Math.sqrt(2)), njs.sub( sq, p1 ) )));
          */

          let tp = ( (u[1]*(p[0]-q[0])) - (u[0]*(p[1]-q[1])) ) / ((u[0]*v[1]) - (u[1]*v[0]));
          let tI = Math.abs(Math.ceil( (tp/(Math.sqrt(2)*ds)) - l0 ));

          if (DEBUG_LEVEL > 1) {

            console.log("#v[", Rk, "]:", v);


            let sq = njs.add(p, njs.mul(tp, v));
            console.log("# q (point[", q_idx, "])");
            console.log(q[0], q[1]);
            console.log(sq[0], sq[1],  "\n");

            console.log("# p diag");
            console.log(p[0], p[1]);
            console.log(sq[0], sq[1], "\n");

            console.log("# tp:", tp, "tI:", tI);
          }

          let quadrent_idx = octant2quadrent[Rk];

          if (DEBUG_LEVEL > 1) {
            console.log("# tI:", tI, ", Rk:", Rk, ", quad:", quadrent_idx, ", q_fence_idx:", p_fence_idx, "(q[", q_idx, "]:", q,")");
          }

          if (tI < p_fence_idx[quadrent_idx]) {
            p_fence_idx[quadrent_idx] = tI;

            if (DEBUG_LEVEL > 1) {
              console.log("# adding tI (p_fence_idx now:", p_fence_idx, ")");
            }

          }

        }
      }
    }

    if (DEBUG_LEVEL > 2) {
      console.log("#pif_list.length:", pif_list.length);
    }

    for (let i=0; i<pif_list.length; i++) {

      let q_idx = pif_list[i];

      let _found = true;
      for (let j=0; j<pif_list.length; j++) {
        if (i==j) { continue; }

        let u_idx = pif_list[j];

        if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {
          _found = false;
          break;
        }
      }
      if (_found) {
        E.push([p_idx, q_idx]);
      }
    }


  // for p_idx
  }

  return { "P": P, "E": E };
}

function check_answer(ctx) {

  let info = ctx.info;

  let grid_n = info.grid_n;

  let P = info.P;
  let E = info.E;
  let A = info.A;

  let cell_size = info.grid_cell_size;

  let v_idir = [
    [1,0,0], [-1,0,0],
    [0,1,0], [0,-1,0],
    [0,0,1], [0,0,-1]
  ];

  let p_nei = [];
  for (let i=0; i<P.length; i++) {
    p_nei.push([]);
  }

  for (let e_idx=0; e_idx < E.length; e_idx++) {
    let p_idx = E[e_idx][0];
    let q_idx = E[e_idx][0];
    p_nei[p_idx].push(q_idx);
    p_nei[q_idx].push(p_idx);
  }

  let bbox = [ [0,0,0], [1,1,1] ];
  let grid_bbox = [[0,0,0], [grid_n, grid_n, grid_n]];
  for (let p_idx=0; p_idx < p_nei.length; p_idx++) {

    for (let nei_idx=0; nei_idx < p_nei[p_idx].length; nei_idx++) {
      let q_idx = p_nei[p_idx][nei_idx];

      let r = njs.norm2( P[p_idx], P[q_idx] );

      if (nei_idx == 0) {
        for (let xyz=0; xyz<3; xyz++) {
          bbox[0][xyz] = P[p_idx][xyz] - r;
          bbox[1][xyz] = P[p_idx][xyz] + r;
        }
      }

      for (let xyz=0; xyz<3; xyz++) {
        if ((P[q_idx][xyz]-r) < bbox[0][xyz]) {
          bbox[0][xyz] = P[q_idx][xyz]-r;
        }
        if ((P[q_idx][xyz]+r) > bbox[1][xyz]) {
          bbox[1][xyz] = P[q_idx][xyz]+r;
        }
      }

    }

    let isxyz = [
      Math.floor( bbox[0][0] / cell_size[0] ),
      Math.floor( bbox[0][1] / cell_size[1] ),
      Math.floor( bbox[0][2] / cell_size[2] )
    ];

    let iexyz = [
      Math.ceil( bbox[1][0] / cell_size[0] ),
      Math.ceil( bbox[1][1] / cell_size[1] ),
      Math.ceil( bbox[1][2] / cell_size[2] )
    ];


    let pnt_idx = [];
    for (let iz=isxyz[2]; iz<iexyz[2]; iz++) {
      for (let iy=isxyz[1]; iy<iexyz[1]; iy++) {
        for (let ix=isxyz[0]; ix<iexyz[0]; ix++) {
          if (oob([ix,iy,iz], grid_bbox)) { continue; }
          let grid_bin = info.grid[iz][iy][ix];

          for (let i=0; i<grid_bin.length; i++) {
            pnt_idx.push(grid_bin[i]);
          }

        }
      }
    }

    for (let i=0; i<pnt_idx.length; i++) {
      let q_idx = pnt_idx[i];
      if (p_idx == q_idx) { continue; }

      let pq_connect = 1;
      for (let j=0; j<pnt_idx.length; j++) {
        let u_idx = pnt_idx[j];
        if (u_idx == p_idx) { continue; }
        if (u_idx == q_idx) { continue; }

        if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {
          pq_connect = 0;
          break;
        }

      }

      if ( A[p_idx][q_idx] != pq_connect ) {
        console.log("!!!ERROR: p_idx:", p_idx, "q_idx:", q_idx, "A:", A[p_idx][q_idx], "pq_connect:", pq_connect);
        return 0;
      }

    }

  }

  return 1;
}

function print_edge(P, E) {
  for (let i=0; i<P.length; i++) {
    console.log("#", i);
    console.log(P[i][0], P[i][1], "\n");
  }

  for (let e_idx=0; e_idx<E.length; e_idx++) {
    let p = P[E[e_idx][0]];
    let q = P[E[e_idx][1]];

    console.log(p[0], p[1]);
    console.log(q[0], q[1]);
    console.log("");
  }

}


function __cruft() {
  let res0 = gen_instance_2d_fence(10000);
  print_edge(res0.P, res0.E);
  process.exit();

  var prof_ctx = {};
  for (let n=10000; n<100000; n+=1000) {
    console.log(n);
    prof_s(prof_ctx, "n" + n.toString());
    let res0 = gen_instance_2d_fence(n);

    prof_e(prof_ctx, "n" + n.toString());

    prof_print(prof_ctx);
    
    //print_edge(res0.P, res0.E);
  }


  process.exit();

  let res1 = gen_instance_2d(res0.P.length, undefined, res0.P);
  print_edge(res1.P, res1.E);

  let _lp = lune_points( res0.P[157], res0.P[33] );
  print_point(_lp);

  process.exit();

}


//test_grid_sweep_perim_2d();
//process.exit();


function fence_secure(fence, ir) {
  for (let i=0; i<fence.length; i++) {
    if (fence[i] < 0) { return false; }
    if (fence[i] > ir) { return false; }
  }
  return true;
}

function _gen_instance_2d_fence(n, B) {
  let info = {
    "dim": 2,
    "start": [0,0],
    "size": [1,1],
    "point": [],
    "point_grid_bp": [],
    "grid_cell_size": [-1,-1],
    "bbox": [[0,0], [1,1]],
    "grid": [],
    "edge": []
  };

  let grid_s = Math.sqrt(n);
  let grid_n = Math.ceil(grid_s);

  let ds = 1 / grid_n;

  info.grid_cell_size[0] = ds;
  info.grid_cell_size[1] = ds;
  info.grid_n = grid_n;
  info.grid_s = grid_s;


  //PROFILING
  prof_s(prof_ctx, "init_grid");
  //PROFILING


  // alloc grid
  //
  for (let i=0; i<grid_n; i++) {
    info.grid.push([]);
    for (let j=0; j<grid_n; j++) {
      info.grid[i].push([]);
    }
  }

  let grid_size  = [ 1, 1 ];
  let grid_start = [ 0,0 ];
  let grid_cell_size = [ ds, ds ];

  info.start = grid_start;
  info.size = grid_size;

  // alloc and create random points
  //
  for (let i=0; i<n; i++) {
    let pnt = [ Math.random()*grid_size[0] + grid_start[0], Math.random()*grid_size[1] + grid_start[1] ];
    info.point.push(pnt);
    info.point_grid_bp.push([-1,-1]);
    //info.edge.push([]);
  }

  //PROFILING
  prof_e(prof_ctx, "init_grid");
  //PROFILING

  //PROFILING
  prof_s(prof_ctx, "setup");
  //PROFILING


  // push points into grid, linear linked list/array for dups
  //
  for (let i=0; i<n; i++) {
    //let ix = Math.floor(info.point[i][0]*grid_s);
    //let iy = Math.floor(info.point[i][1]*grid_s);
    let ix = Math.floor(info.point[i][0]*grid_n);
    let iy = Math.floor(info.point[i][1]*grid_n);
    info.grid[iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy];
  }

  //PROFILING
  prof_e(prof_ctx, "setup");
  //PROFILING


  let P = info.point;
  let G = info.grid;

  let E = [];

  let octant_dir_map = [ 0, 1, 1, 2, 2, 3, 3, 0 ];

  // 0 - right
  // 1 - left
  // 2 - up
  // 3 - down
  //

  // octant_index_2d(p,q)

  for (let p_idx=0; p_idx < P.length; p_idx++) {
    let pnt = P[p_idx];

    let i_anch_x = Math.floor(pnt[0] / ds);
    let i_anch_y = Math.floor(pnt[1] / ds);

    let fence_q = [-1,-1,-1,-1];

    let iR = 0;
    while (iR < grid_n) {

      if (fence_secure(i_r, fence_q)) { break; }

      for (let g_idx=0; g_idx < G[i_anch_y][i_anch_x].length; g_idx++) {
        let q_idx = G[i_anch_y][i_anch_x][g_idx];
        if (q_idx == p_idx) { continue; }

      }

      iR++;
    }

  }

}

/*
let a = lune_points([0,0], [0,1], 32, 1);
print_point(a);
console.log("\n");

a = lune_points([0,0], [Math.sqrt(2), -1], 32, 1);
print_point(a);
console.log("\n");

a = lune_points([0,0], [1, -Math.sqrt(2)/2], 32, 1);
print_point(a);
console.log("\n");

process.exit();
*/

// trying to figure out why there's a polynomial time blowup
//

//info = gen_instance_2d(3000, [[0,0],[1,1]]);
//process.exit();

function cell_stat_experiment() {

  for (let it=0; it<10; it++) {
    console.log("#it:", it);
    let info = gen_instance_2d(10000, [[0,0],[1,1]]);
  }

  for (let dist in _debug_stat.count) {
    console.log("#", dist, _debug_stat.val[dist], _debug_stat.count[dist] );
    console.log( dist, _debug_stat.val[dist] / _debug_stat.count[dist]);
  }
  process.exit();
}

function __main() {
  let _single_run = false;
  if (_single_run) {
    let info = gen_instance_2d(1000, [[0,0],[1,1]]);
    print_point(info.point, 1);
    process.exit();

  }
  else {

    for (let t=100; t<5000; t+=100) {

      let s_t = new Date().getTime();

      let info = gen_instance_2d(t, [[0,0],[1,1]]);
      //print_point(info.point, 1);

      let e_t = new Date().getTime();

      console.log("##", t, (e_t - s_t)/1000);
    }
    process.exit();

  }

}

// CRUF!!!

function _xxx() {
  let anchor_pnt = [ Math.random(), Math.random() ];

  for (let ir=0; ir<5; ir++) {
    for (let r_idx=0; r_idx<8; r_idx++) {
      console.log("###", ir, r_idx);
      let gp = grid_sweep_2d([0.05,0.05], anchor_pnt, ir, r_idx);


      for (let ii=0; ii<gp.n; ii++) {
        //console.log(gp.s[0] + gp.v[0]*ii, gp.s[1] + gp.v[1]*ii);
        console.log(gp.S[0] + gp.V[0]*ii, gp.S[1] + gp.V[1]*ii);
      }

      console.log("");

      //console.log(gp);
    }
    console.log("\n");
  }

  process.exit();

  let naive_res = naive_relnei_E(pnt);
  let A = naive_res.A;

  let print_lune = true;
  let print_graph = true;

  //print_point(pnt);
  if (print_lune) {
    for (let i=0; i<N; i++) {
      for (let j=(i+1); j<N; j++) {
        if (A[i][j] < 0.5) { continue; }
        let a = pnt[i];
        let b = pnt[j];
        let lune_pnt = lune_points(a,b);
        print_point(lune_pnt);
      }
    }
  }

  if (print_graph) {
    for (let i=0; i<N; i++) {
      for (let j=(i+1); j<N; j++) {
        if (A[i][j] < 0.5) { continue; }
        let a = pnt[i];
        let b = pnt[j];
        console.log(a[0], a[1]);
        console.log(b[0], b[1]);
        console.log("");
      }
    }
  }

}

function failing0() {
  let P = [
    [0.44584097480654905,0.2630677410107003,0.978063734762359],
    [0.8792005129276691,0.9569835189557732,0.4981144569316643],

    // !!
    [0.7302476517208656,0.4085662385644728,0.8950140702477447],
    // !!

    [0.4021124780579683,0.3027373978303754,0.6800684457807135],

    // !!
    [0.40536140905436313,0.6930924860922895,0.8439412700390785],
    // !!

    [0.8046564969848393,0.02887287542893153,0.9445998333212928],
    [0.4472281866578808,0.7399787104118357,0.940831312948197],
    [0.4024911152400458,0.9972875928707101,0.6140172162073474],
    [0.45340850104715535,0.48780232400880913,0.2571819484664662],
    [0.1309191022868963,0.29401818315842815,0.024937234825944187]
  ];

  let res_fence = gen_instance_3d_fence(P.length, [[0,0,0],[1,1,1]], P);
  let res_naive = naive_relnei_E(P);

  print_point(res_fence.P, 1);
  print_E(res_fence.P, res_fence.E);

  let _cmp_res = check_cmp(res_fence, res_naive.A);

  if (!_cmp_res) {
    console.log("# mismatch!!");
  }
}

function failing1() {
  let P = [

[0.48236153119999786,0.5829562578210701,0.5894499262133782],
[0.40226423198824285,0.18571662125302568,0.8676381195370179],
[0.33215985928550007,0.7072200784349136,0.8788649781571626],
[0.48895317751943884,0.3456950558866528,0.7701554499785753],
[0.2902892853701606,0.17521479079563637,0.9474716241453544],
[0.2801319609391316,0.9000756739056867,0.10790588342455201],
[0.7257716226289987,0.9872984729307688,0.6558161569627607],
[0.20155954725594033,0.05141859488332581,0.4873397088251342],

  //!!
[0.38641533071261625,0.19886654954131305,0.1412455118739861],
  //!!

[0.9430324013584659,0.5791526363231108,0.8927393680820273],
[0.4569842977471846,0.9361777812913996,0.8644125468436177],
[0.044553879420293645,0.643883319861713,0.9199383885412326],
[0.3995710006357906,0.3124296141199602,0.7713934061210598],
[0.052040577538081334,0.47600458330531326,0.4871463395672549],

//!!
[0.37520407169511805,0.7350236931845794,0.1331874245080567],
[0.7484014367367805,0.3342644552089469,0.12295008151606544],
//!!

[0.8238535262366712,0.9172865900773475,0.629972129987517],
[0.26076817231777033,0.798646842066488,0.1639779261157275],
[0.1922848467550962,0.6117787464338532,0.9838191584318533],
[0.9281548838754002,0.8630286046201179,0.16849134805951754]
];


  let res_fence = gen_instance_3d_fence(P.length, [[0,0,0],[1,1,1]], P);
  let res_naive = naive_relnei_E(P);

  print_point(res_fence.P, 1);
  print_E(res_fence.P, res_fence.E);

  let _cmp_res = check_cmp(res_fence, res_naive.A);

  let p_idx = 14;
  let q_idx = 15;
  let u_idx = 8;

  let info = res_fence.info;

  let cell_size = info.grid_cell_size;
  console.log("#ds:", cell_size);

  console.log("# u_idx:", u_idx, "in lune(", p_idx, q_idx,"):", in_lune(P[p_idx], P[q_idx], P[u_idx]));
  console.log("# u_idx:", u_idx, "in lune(", q_idx, p_idx,"):", in_lune(P[q_idx], P[p_idx], P[u_idx]));

  let p_ixyz = [
    Math.floor( P[p_idx][0] / cell_size[0]),
    Math.floor( P[p_idx][1] / cell_size[1]),
    Math.floor( P[p_idx][2] / cell_size[2])
  ]

  let q_ixyz = [
    Math.floor( P[q_idx][0] / cell_size[0]),
    Math.floor( P[q_idx][1] / cell_size[1]),
    Math.floor( P[q_idx][2] / cell_size[2])
  ]

  let u_ixyz = [
    Math.floor( P[u_idx][0] / cell_size[0]),
    Math.floor( P[u_idx][1] / cell_size[1]),
    Math.floor( P[u_idx][2] / cell_size[2])
  ];
  console.log("# p_idx:", p_idx, "P[", p_idx,"]:", P[p_idx], "ixyz:", p_ixyz);
  console.log("# q_idx:", q_idx, "P[", q_idx,"]:", P[q_idx], "ixyz:", q_ixyz);
  console.log("# u_idx:", u_idx, "P[", u_idx,"]:", P[u_idx], "ixyz:", u_ixyz);

  for (let u_idx=0; u_idx < P.length; u_idx++) {
    if ((u_idx == p_idx) || (u_idx == q_idx)){ continue; }
    if (in_lune(P[p_idx], P[q_idx], P[u_idx])) {
      console.log("##!!!! u_idx:", u_idx, "in lune(", p_idx, q_idx, ")");
    }
  }

  //let u_idx = 15;
  //P = res_fence.P;
  //console.log("#0: u_idx", u_idx, "in lune(", p_idx, q_idx, ")?", in_lune(P[p_idx], P[q_idx], P[u_idx]));
  //P = res_naive.P;
  //console.log("#1: u_idx", u_idx, "in lune(", p_idx, q_idx, ")?", in_lune(P[p_idx], P[q_idx], P[u_idx]));

  //console.log("#???", res_naive.A[14][15], res_naive.A[15][14]);

  if (!_cmp_res) {
    console.log("# mismatch!!");
  }
  else {
    console.log("# ok, match");
  }
}

function debug_bad_2d(fn) {
  let data = fs.readFileSync( fn ).toString().split("\n");

  let _pnt = [];

  for (let i=0; i<data.length; i++) {
    if (data[i].length == 0) { continue; }
    if (data[i][0] == '#') { continue; }
    let tok = data[i].split(" ");
    if (tok.length != 2) { continue; }
    let x = parseFloat(tok[0]);
    let y = parseFloat(tok[1]);
    _pnt.push([x,y]);
  }

  console.log("##...", _pnt.length);


  let info = gen_instance_2d_fence(_pnt.length, _pnt);
  let naive_res = naive_relnei_E(_pnt);

  let _cmp_res = check_cmp(info, naive_res.A);

  console.log("# _cmp:", _cmp_res);

}


/*
  //console.log("# DEBUG BAD 2d 100");
  //debug_bad_2d("bad.100");
  console.log("# DEBUG BAD 2d 400");
  debug_bad_2d("bad_2d.400");
  process.exit();
*/

function main_test2d() {
  let Ntest = [
    10,10,10,20,20,20,30,30,30,
    /*
    10,10,10,20,20,20,30,30,30,
    10,10,10,20,20,20,30,30,30,
    10,10,10,20,20,20,30,30,30,
    10,10,10,20,20,20,30,30,30,
    10,10,10,20,20,20,30,30,30,
    */
    50,50,
    100,100,
    200,200,
    300,400,
    500, 600,
    700, 800,
    900, 1000
  ];

  console.log("# main_test_2d");

  for (let n_idx=0; n_idx < Ntest.length; n_idx++) {
    let N = Ntest[n_idx];
    let _pnts = poisson_point(N, 2);

    console.log("## test N:", N, "(", n_idx, "/", Ntest.length, ")");


    console.log("# fence start");

    let info = gen_instance_2d_fence(N, _pnts);

    console.log("# fence done");

    //print_point(info.P, 1);
    //print_E(info.P, info.E);

    console.log("# naive start");

    let naive_res = naive_relnei_E(info.P);

    console.log("# naive done");

    //print_E(naive_res.P, naive_res.E);

    let _cmp_res = check_cmp(info, naive_res.A);

    if (!_cmp_res) {
      console.log("#check failed, points:");
      print_point(info.P);
      break;
    }

    console.log("#", N, "{", n_idx,"}", _cmp_res);

  }


}

function _avgAdeg(A) {
  let tot = 0;
  for (let i=0; i<A.length; i++) {
    for (let j=0; j<A[i].length; j++) {
      if (A[i][j] > 0) { tot++; }
    }
  }

  return tot / A.length;
}

function main_test3d() {

  //failing0();
  //failing1();
  //return;

  let Ntest = [
    10,10,10,20,20,20,30,30,30,
    50,50,
    100,100,
    200,200,
    300,400,
    500, 600,
    700, 800,
    900, 1000
  ];

  //Ntest = [ 1000, 2000, 3000 ];

  /*
  Ntest = [];
  for (let i=0; i<100; i++) {
    Ntest.push(20);
  }
  */

  for (let n_idx=0; n_idx < Ntest.length; n_idx++) {
    let N = Ntest[n_idx];
    let _pnts = poisson_point(N, 3);


    console.log("# fence start");

    let info = gen_instance_3d_fence(N, [[0,0,0],[1,1,1]], _pnts);

    console.log("# fence done");

    //print_point(info.P, 1);
    //print_E(info.P, info.E);

    console.log("# naive start");

    let naive_res = naive_relnei_E(info.P);

    console.log("# naive done");

    console.log("## naive avg deg:", _avgAdeg(naive_res.A));

    //print_E(naive_res.P, naive_res.E);

    let _cmp_res = check_cmp(info, naive_res.A);

    if (!_cmp_res) {
      console.log("#check failed, points:");
      print_point(info.P);
      break;
    }

    console.log("#", N, "{", n_idx,"}", _cmp_res);

  }


  process.exit();


  let _n = 100000;
  let _pnt = poisson_point(_n, 2);
  let res = gen_instance_2d_fence(_n, _pnt);

  process.exit();
}

function create_A(info) {
  let n = info.P.length;

  let A = [];
  for (let i=0; i<n; i++) {
    A.push([]);
    for (let j=0; j<n; j++) {
      A[i].push(0);
    }
  }

  let E = info.E;

  for (let i=0; i<E.length; i++) {
    A[ E[i][0] ][ E[i][1] ] = 1;
    A[ E[i][1] ][ E[i][0] ] = 1;
  }

  return A;
}

function _main() {
  let N = 500;
  let info = gen_instance_3d_fence(N, [[0,0,0],[1,1,1]]);

  info["A"] = create_A(info);
  info["info"]["A"] = info["A"];

  //print_point(info.P, 1);
  //print_E(info.P, info.E);

  check_answer(info);
}

function naive_avg_deg_3d() {
  let M = [ 20, 50, 100, 500 ];
  let rep = 10;
  for (let m_idx=0; m_idx<M.length; m_idx++) {
    let N = M[m_idx];
    let avg_sum = 0;
    for (let it=0; it<rep; it++) {
      let _pnt = poisson_point(N, 3);
      let res = naive_relnei_E(_pnt);
      avg_sum += _avgAdeg(res.A);
    }

    console.log(N, avg_sum / rep);
  }
}

function fence_avg_deg_3d() {
  let M = [ 20, 50, 100, 500, 1000, 2000, 4000, 5000, 10000 ];
  let rep = 10;
  for (let m_idx=0; m_idx<M.length; m_idx++) {
    let N = M[m_idx];
    let avg_sum = 0;
    for (let it=0; it<rep; it++) {
      let _pnt = poisson_point(N, 3);
      let info = gen_instance_3d_fence(N, [[0,0,0],[1,1,1]], _pnt);
      avg_sum += info.E.length / _pnt.length;
    }

    console.log(N, avg_sum / rep);
  }
}

function fence_dist_plot_3d() {
  let N = 10000;
  let info = gen_instance_3d_fence(N);

  let avg_dist = 0;
  for (let i=0; i<info.E.length; i++) {
    let p_idx = info.E[i][0];
    let q_idx = info.E[i][1];

    let p = info.P[p_idx];
    let q = info.P[q_idx];

    console.log(njs.norm2(njs.sub(p,q)), 1);
  }
}

function main() {
  let N = 7000;
  let info = gen_instance_3d_fence(N, [[0,0,0],[1,1,1]]);

  print_E(info.P, info.E);
}



main_test2d();
//main();
//main_test3d();

//fence_avg_deg_3d();
//fence_dist_plot_3d();

