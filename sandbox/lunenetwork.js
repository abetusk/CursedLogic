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

var N = 10;
var pnt = [];

function poisson_point(N, D) {
  D = ((typeof D === "undefined") ? 2 : D);
  let pnt = [];
  for (let i=0; i<N; i++) {
    let p = [];
    for (let j=0; j<D; j++) {
      p.push( Math.random() );
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
    if (disconnect) { console.log(""); }
  }
  console.log("");
}

function in_lune(pnt_a, pnt_b, tst_c) {

  let dist_ca = njs.norm2( njs.sub(tst_c, pnt_a) );
  let dist_cb = njs.norm2( njs.sub(tst_c, pnt_b) );
  let dist_ab = njs.norm2( njs.sub(pnt_a, pnt_b) );

  console.log("# pnt_a(", pnt_a, "), pnt_b(", pnt_b, "), tst_c(", tst_c, ")");
  console.log("# dist_ca:", dist_ca, "dist_cb:", dist_cb, "dist_ab:", dist_ab);

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

N = 10;
pnt = poisson_point(N, 2);

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

  return A;
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

function gen_instance_2d(n, B) {

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
    let pnt = [ Math.random()*grid_size[0] + grid_start[0], Math.random()*grid_size[1] + grid_start[1] ];
    info.point.push(pnt);
    info.point_grid_bp.push([-1,-1]);
    info.edge.push([]);
  }

  for (let i=0; i<n; i++) {
    let ix = Math.floor(info.point[i][0]*grid_s);
    let iy = Math.floor(info.point[i][1]*grid_s);
    info.grid[iy][ix].push(i);
    info.point_grid_bp[i] = [ix,iy];
  }

  let P = info.point;
  let G = info.grid;

  //console.log("# ds:", ds, "grid_s:", grid_s);

  for (let p_idx=0; p_idx < P.length; p_idx++) {
    let pnt = P[p_idx];

    let wedge_resolved  = [ 0, 0,  0, 0,  0, 0,  0, 0];
    let wedge_nei       = [-2,-2, -2,-2, -2,-2, -2,-2];

    let resolved_count = 0;

    let i_anch_x = Math.floor(pnt[0] / ds);
    let i_anch_y = Math.floor(pnt[1] / ds);

    // just easier to enumerate all points in initial 3x3 grid region
    //
    for (let idy=-1; idy<2; idy++) {
      for (let idx=-1; idx<2; idx++) {

        let ix = i_anch_x + idx;
        let iy = i_anch_y + idy;

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

    //console.log("# p_idx:", p_idx, P[p_idx], "wedge_nei:", wedge_nei);

    for (let r_idx=0; r_idx<8; r_idx++) {

      if (wedge_nei[r_idx] == -1) { continue; }

      for (let ir=2; ir<grid_n; ir++) {
        let wedge_info = grid_sweep_2d(info, P[p_idx], ir, r_idx);

        //console.log(">> ir:", ir, "wedge_info:", wedge_info);

        if (wedge_info.path.length == 0) {
          wedge_nei[r_idx] = -1;
          resolved_count++;
          break;
        }

        for (let w_idx=0; w_idx < wedge_info.path.length; w_idx++) {
          let grid_point = wedge_info.path[w_idx];
          let ix = grid_point[0];
          let iy = grid_point[1];

          console.log("#>>>>", ix, iy, "grid_n:", grid_n);

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
          wedge_nei[r_idx] = -1;
        }

      }
    }

    console.log("# p_idx:", p_idx, "pnt:", pnt, "resolved:", resolved_count);

    for (let r_idx=0; r_idx<8; r_idx++) {
      let q_idx = wedge_nei[r_idx];
      if (q_idx == -2) {
        console.log("### ERROR!!!", "p_idx:", p_idx, P[p_idx], "q_idx:", q_idx, "r_idx:", r_idx);
        return;
      }
      if (q_idx == -1) { continue; }
      console.log(P[p_idx][0], P[p_idx][1]);
      console.log(P[q_idx][0], P[q_idx][1]);
      console.log("");
    }

  }

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
        let p = info.point[idx];
        console.log("#", ii, ix, iy, "-> [", p[0],p[1] ,"] {", idx, "}");
      }
    }
  }

  return info;
}


let info = gen_instance_2d(10);
print_point(info.point, 1);
process.exit();


// CRUF!!!

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

let A = naive_relnei_E(pnt);

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


