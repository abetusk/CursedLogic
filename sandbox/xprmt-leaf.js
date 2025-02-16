//
// To the extent possible under law, the person who associated CC0 with
// this code has waived all copyright and related or neighboring rights
// to this code.
//
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

// export for commonJS or browser
//
var qt;
if ((typeof module !== 'undefined') &&
    (typeof module.exports !== 'undefined')) {
  qt = require("./js-quadtree.js");
  var halton = require("./halton.js");
}
else {
  qt = QT;
  //halton = halton;
}



let N = 2500;
N = 5500;
let S = 10000;
let rho =  S / Math.sqrt(N);

let VERBOSE = 0;

function sanity_vein() {

  let vein = g_info.vein;

  let leaf_v = [];
  let leaf_map = {};
  let leaf_count=0;

  for (let ii=0; ii<vein.length; ii++) {
    let v = vein[ii];
    if (v.state == 0) {
      leaf_v.push(v);
      leaf_map[ v.idx ] = 1;
      leaf_count++;
    }
  }

  for (let ii=0; ii<vein.length; ii++) {
    let v = vein[ii];
    if (v.p_idx < 0) { continue; }

    if (v.p_idx in leaf_map) {
      console.log("v", v, " has parent ", vein[v.p_idx], " but parent is leaf?");
      return -1;
    }
  }

  console.log(leaf_count);

  return 0;
}

function new_vein(info, x, y, dx, dy, p_idx, c_idx) {
  p_idx = ((typeof p_idx === "undefined") ? -1 : p_idx);
  c_idx = ((typeof c_idx === "undefined") ? -1 : c_idx);
  let vein_info = {
    "x" : x,
    "y" : y,
    "p_idx" : p_idx,
    "c_idx" : c_idx,
    "type"  : "vein",
    "state" : 0,
    "a_idx" : -1,
    "d"     : -1,
    "idx"   : info.vein_id,
    "auxin_count" : 0,
    "dir_count"   : 0,
    "orig_dx" : dx,
    "orig_dy" : dy,
    "dx"      : dx,
    "dy"      : dy
  };

  info.vein_id++;

  return vein_info;
}


let g_info = {
  "n": N,
  "bb": {"x":-S/2, "y":-S/2, "w":S, "h": S},
  "dx": S, "dy": S,

  "ds": 125,
  "max_iter" : 1000,

  "rho": rho,

  //"r_kill": 0.5*rho,
  //"r_bias": 1.0*rho,

  //"r_kill": 1.125*rho,
  "r_kill": 1.125*rho,

  //"r_bias": 1.25*rho,
  //"r_bias": 2.75*rho,
  "r_bias": 2.5*rho,

  "eps": 1/(1024*1024),
  "auxin" : [],
  "vein": [],
  "wf_idx": [],

  "parent_vein_weight": 0.0,
  "parent_vein_denom": 1,

  "auxin_tree" : {},
  "vein_tree": {},

  "disp_scale" : 400/S,
  "disp_dxy" : [ S / 2, S / 2 ],
  "disp_pnt" : []


};

function _transform(pnt, dx, dy, s) {
  for (let ii=0; ii<pnt.length; ii++) {
    pnt[ii].x = (pnt[ii].x + dx)*s;
    pnt[ii].y = (pnt[ii].y + dy)*s;
  }
  return pnt;
}

function _transforml(line, dx, dy, s) {
  for (let ii=0; ii<line.length; ii++) {
    line[ii].x = (line[ii].x + dx)*s;
    line[ii].y = (line[ii].y + dy)*s;

    line[ii].dx *= s;
    line[ii].dy *= s;
  }
  return line;
}

function _transformL(line, dx, dy, s) {
  for (let ii=0; ii<line.length; ii++) {
    line[ii].x0 = (line[ii].x0 + dx)*s;
    line[ii].y0 = (line[ii].y0 + dy)*s;

    line[ii].x1 = (line[ii].x1 + dx)*s;
    line[ii].y1 = (line[ii].y1 + dy)*s;
  }
  return line;
}

function space_col_init() {
  //let rand_alg = "halton"; // "random"
  let rand_alg = "random";

  g_info.auxin = [];

  let tree = new qt.QuadTree( new qt.Box( g_info.bb.x, 
                                          g_info.bb.y,
                                          g_info.bb.w,
                                          g_info.bb.h ) );
  g_info.auxin_tree = tree;

  let halton_cx = 0.5,
      halton_cy = 0.5;

  let seq = [];
  if (rand_alg == "halton") {
    seq = halton.seq2d(2,3,g_info.n);
  }
  else if (rand_alg == "random") {
    seq = [];
    for (let ii=0; ii<g_info.n; ii++) {
      seq.push([ Math.random(), Math.random() ]);
    }
  }

  for (let ii=0; ii<seq.length; ii++) {

    seq[ii][0] -= halton_cx;
    seq[ii][1] -= halton_cy;

    let _x = Math.floor(g_info.dx*seq[ii][0]);
    let _y = Math.floor(g_info.dy*seq[ii][1]);
    let pnt_info = {
      "x" : _x,
      "y" : _y,
      "idx": ii,
      "type": "auxin",
      "state": 0,
      "d": -1,
      "v_list": [],
      "v_idx": -1
    };
    g_info.auxin.push( new qt.Point(_x,_y, pnt_info) );

  }

  tree.insert(g_info.auxin);

}

// The algorihm iterates through five phases:
//
// * auxin kill       - kill auxin nodes that are too close to vein nodes
// * auxin init       - sweep nearby auxin nodes from vein nodes and initialize them
// * auxin queue      - associate the appropriate vein nodes to each auxin node within
//                      influence distance
// * auxin pull       - gather forces for each influenced vein node
// * vein wavefront   - create new vein nodes from vein nodes that have at least
//                      one auxin node influencing it, keeping vein nodes already in
//                      the wavefront if they are still within influence radius of
//                      an auxin node
//
// Initially the vein wavefront list has the root node.
// The above steps are iterated until the vein wave front list is cleared or a maximum
// iteration count is exceeded.
//
// Auxin nodes are kept in a quadtree for fast neighbor querying.
//
// g_info.vein holds the vein nodes.
//
//
function open_space_col_start() {
  let tree = g_info.auxin_tree;

  // root point
  //
  g_info.vein_id = 0;
  let root_vein = new_vein(g_info, 0, g_info.bb.y, 0, g_info.ds);
  g_info.vein.push(root_vein);

  let query_kill = new qt.Circle(0, 0, g_info.r_kill);
  let query_bias = new qt.Circle(0, 0, g_info.r_bias);

  if (VERBOSE > 0) {
    console.log("# r_kill:", g_info.r_kill, ", r_bias:", g_info.r_bias);
  }

  let vein_wf = [ g_info.vein[0] ];

  let _time_s = Date.now();

  let iter = -1;
  while ( (iter < g_info.max_iter) &&
          (vein_wf.length > 0) ) {
    iter++;

    if (VERBOSE > 0) {
      console.log("### iter:", iter, ", vein_wf.length:", vein_wf.length);
    }

    let auxin_queue = [];
    let nxt_vein_wf = [];

    // AUXIN KILL
    //
    // init current vein wavefront and do a kill sweep
    //
    for (let ii=0; ii<vein_wf.length; ii++) {
      let v = vein_wf[ii];

      v.auxin_count = 0;
      v.dir_count = 0;

      query_kill.x = v.x;
      query_kill.y = v.y;
      let killpnt = tree.query(query_kill);
      for (let _k=0; _k<killpnt.length; _k++) {

        if (VERBOSE > 0) {
          console.log("# v:", v.x, v.y, "kill", killpnt[_k].data.x, killpnt[_k].data.y);
        }

        tree.remove( killpnt[_k] );
      }
    }

    // AUXIN INIT
    //
    // init auxin sweep
    //
    for (let ii=0; ii<vein_wf.length; ii++) {
      let v = vein_wf[ii];
      query_bias.x = v.x;
      query_bias.y = v.y;
      let querypnt = tree.query(query_bias);
      for (let _k=0; _k<querypnt.length; _k++) {
        let a = querypnt[_k];
        a.data.v_idx = -1;

        if (VERBOSE > 0) {
          console.log("# v:", v.x, v.y, ", auxin init:", a.data.x, a.data.y);
        }
      }
    }

    // AUXIN QUEUE
    //
    // distance query
    //
    for (let ii=0; ii<vein_wf.length; ii++) {
      let v = vein_wf[ii];

      query_bias.x = v.x;
      query_bias.y = v.y;
      let querypnt = tree.query(query_bias);
      for (let _k=0; _k<querypnt.length; _k++) {
        let a = querypnt[_k];

        v.auxin_count++;

        let _dx = a.data.x - v.x;
        let _dy = a.data.y - v.y;
        let _dist = Math.sqrt(  _dx*_dx + _dy*_dy );

        if (a.data.v_idx < 0) {
          auxin_queue.push(a);
        }
        else { }

        if ((a.data.v_idx < 0) ||
            (_dist < a.data.d)) {
          a.data.d = _dist;
          a.data.v_idx = v.idx;

          if (VERBOSE > 0) {
            console.log("#> updating auxin (", a.data.x, a.data.y, "),  v_idx:", a.data.v_idx, ", dist:", a.data.d );
          }
          //console.log(a, v);
        }
      }
    }

    // AUXIN PULL
    //
    // go through all auxin nodes that have influence over vein nodes
    // and update their thralled vein nodes
    //
    for (let ii=0; ii<auxin_queue.length; ii++) {

      let a = auxin_queue[ii];
      let v_idx = a.data.v_idx;
      let v = g_info.vein[v_idx];

      let _dx = a.data.x - v.x;
      let _dy = a.data.y - v.y;
      let _dist = Math.sqrt( _dx*_dx + _dy*_dy );

      let ds = g_info.ds;

      v.dx += ds*_dx / _dist;
      v.dy += ds*_dy / _dist;
      v.dir_count++;

    }

    // VEIN WAVEFRONT
    //
    // remove vein elements that have no auxin nodes within bias distance
    //
    for (let ii=0; ii<vein_wf.length; ii++) {
      let v = vein_wf[ii];
      if (v.auxin_count > 0) { nxt_vein_wf.push( v ); }
      if (v.dir_count > 0) {

        let _den = g_info.parent_vein_denom;
        let dx = v.dx / (v.dir_count+_den);
        let dy = v.dy / (v.dir_count+_den);

        //let dst_v = new_vein(g_info, v.x + dx, v.y + dy, dx, dy, v.idx);

        let new_dx = 0, new_dy = 0;
        let _w = g_info.parent_vein_weight;
        if (v.p_idx >= 0) {
          new_dx = _w*(v.x - g_info.vein[v.p_idx].x);
          new_dy = _w*(v.y - g_info.vein[v.p_idx].y);
        }
        let dst_v = new_vein(g_info, v.x + dx, v.y + dy, new_dx, new_dy, v.idx);
        v.state = 1;

        g_info.vein.push(dst_v);

        nxt_vein_wf.push(dst_v);

        v.dx = dx;
        v.dy = dy;
      }
    }
    vein_wf = nxt_vein_wf;

    if (VERBOSE > 0) {
      for (let ii=0; ii<vein_wf.length; ii++) {
        console.log("# vein_wf[", ii, "/", vein_wf.length, "] (", vein_wf[ii].x, vein_wf[ii].y, ")");
      }
      console.log("#\n#");
    }



  }

  let _time_e = Date.now();
  console.log(_time_s, _time_e, "==>", _time_e - _time_s, "ms to render");

  plot_auxin();
  plot_vein();
}

function plot_auxin() {
  let _auxin = g_info.auxin;
  let disp_pnt = [];
  for (let ii=0; ii<_auxin.length; ii++) {
    let c = "#f003";
    if (_auxin[ii].data.d < 0) { c = "#0f05"; }
    disp_pnt.push({
      "c": c,
      "x": _auxin[ii].x,
      "y": _auxin[ii].y
    });
  }
  _transform(disp_pnt, g_info.disp_dxy[0], g_info.disp_dxy[1], g_info.disp_scale);
  _plot( disp_pnt );
}

function plot_vein() {
  disp_line = [];
  for (let ii=0; ii<g_info.vein.length; ii++) {
    let p_idx = g_info.vein[ii].p_idx;
    if (p_idx < 0) { continue; }

    let v  = g_info.vein[ii];
    let pv = g_info.vein[p_idx];

    disp_line.push({
      "c": "#777",
      "x": pv.x,
      "y": pv.y,
      "dx": v.x - pv.x,
      "dy": v.y - pv.y
    });

  }

  _transforml(disp_line, g_info.disp_dxy[0], g_info.disp_dxy[1], g_info.disp_scale);
  _plotl(disp_line);
}

function plot_ok(line) {

  _transformL(line, g_info.disp_dxy[0], g_info.disp_dxy[1], g_info.disp_scale);
  _plotL(line);
}


function relative_neighborhood_graph() {
  let tree = g_info.auxin_tree;
  let auxin = g_info.auxin;
  let vein = g_info.vein;

  let bbox = new qt.Box( g_info.bb.x, 
                         g_info.bb.y,
                         g_info.bb.w,
                         g_info.bb.h );
  let vein_tree = new qt.QuadTree( bbox );
  g_info.vein_tree = vein_tree;

  g_info.vein_id = 0;
  let root_vein = new_vein(g_info, 0, g_info.bb.y, 0, g_info.ds);
  g_info.vein.push(root_vein);

  vein_tree.insert( new qt.Point(root_vein.x, root_vein.y, root_vein) );

  let query_kill = new qt.Circle(0, 0, g_info.r_kill);
  let query_bias = new qt.Circle(0, 0, g_info.r_bias);
  let query_nei  = new qt.Circle(0, 0, 2*g_info.rho);

  if (VERBOSE > 0) {
    console.log("# r_kill:", g_info.r_kill, ", r_bias:", g_info.r_bias);
  }

  let vein_wf = [ g_info.vein[0] ];

  let _time_s = Date.now();


  for (let it=0; it<2; it++) {

    let _line = [];

    let vein_visited = {};
    let vein_queue = [];

    for (let ii=0; ii<auxin.length; ii++) {
      let a = auxin[ii];

      a.v_list = [];

      query_nei.x = a.x;
      query_nei.y = a.y;
      let v_test = vein_tree.query( query_nei );

      for (let _i=0; _i<v_test.length; _i++) {
        let v0 = v_test[_i];
        let d_a_v0 = Math.sqrt( (v0.x - a.x)*(v0.x - a.x) + (v0.y - a.y)*(v0.y - a.y) );
        let rng_accept = true;
        for (let _j=(_i+1); _j<v_test.length; _j++) {
          let v1 = v_test[_j];
          let d_a_v1 = Math.sqrt( (v1.x - a.x)*(v1.x - a.x) + (v1.y - a.y)*(v1.y - a.y) );
          let d_v0_v1 = Math.sqrt( (v1.x - v0.x)*(v1.x - v0.x) + (v1.y - v0.y)*(v1.y - v0.y) );

          let d = d_a_v1;
          if (d_v0_v1 > d)  { d = d_v0_v1; }

          if (d >= d_a_v0) { rng_accept = false; break; }
        }

        if (rng_accept) {
          console.log(ii, "a", a, "--rng--", v0);

          _line.push( {"c": "#877", "x0": a.x, "y0": a.y, "x1": v0.x, "y1": v0.y } );

          let dav = [ a.x - v0.x, a.y - v0.y ];
          let _d = Math.sqrt( dav.x * dav.x + dav.y * dav.y );
          if (_d < (1/1024)) { _d = 0; }
          else { _d = 1/_d; }

          v0.dir_count++;
          v0.x += _d*dav.x;
          v0.y += _d*dav.y;

          a.v_list.push( v0.data.idx );

          if ( !(v0.data.idx in vein_visited) ) {
            vein_visited[v0.data.idx] = 1;
            vein_queue.push(v0);
          }

        }
      }

    }

    for (let ii=0; ii<vein_queue.length; ii++) {
      let v = vein_queue[ii];

      v.x /= v.dir_count;
      v.y /= v.dir_count;

      let nv = new_vein(g_info, v.x, v.y, v.idx, g_info.ds);


    }

    let x = tree.query( query_nei );
    console.log(">>", x);

    plot_auxin();
    plot_ok(_line);

  }

}

var __iter_iter=0;
function __iter() {
  __iter_iter++;

  space_col_init();
  relative_neighborhood_graph();

}


