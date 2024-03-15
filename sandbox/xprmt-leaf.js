//
// To the extent possible under law, the person who associated CC0 with
// this code has waived all copyright and related or neighboring rights
// to this code.
//
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//


var qt = require("./js-quadtree.js");
var halton = require("./halton.js");

let N = 500;
let S = 10000;
let rho =  S / Math.sqrt(N);

let VERBOSE = 0;

console.log("## N:", N, "S:", S, "rho:", rho);

let g_info = {
  "n": N,
  "bb": {"x":-S/2, "y":-S/2, "w":S, "h": S},
  "dx": S, "dy": S,
  "ds": 50,
  //"r_kill": 0.5*rho,
  //"r_bias": 1.0*rho,
  "r_kill": 0.2*rho,
  //"r_bias": 1.25*rho,
  "r_bias": 7.0*rho,
  "eps": 1/(1024*1024),
  "auxin" : [],
  "vein": [],
  "wf_idx": []
};


let halton_cx = 0.5,
    halton_cy = 0.5;
let tree = new qt.QuadTree( new qt.Box( g_info.bb.x, 
                                        g_info.bb.y,
                                        g_info.bb.w,
                                        g_info.bb.h ) );

let seq = halton.seq2d(2,3,g_info.n);

seq = [];
for (let ii=0; ii< g_info.n; ii++) {
  seq.push( [ Math.random(), Math.random() ] );
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
    "v_idx": -1
  };
  g_info.auxin.push( new qt.Point(_x,_y, pnt_info) );
}

// root point
//

tree.insert(g_info.auxin);

function new_vein(info, x, y, dx, dy) {
  let vein_info = {
    "x" : x,
    "y" : y,
    "type"  : "vein",
    "state" : 0,
    "a_idx" : -1,
    "d"     : -1,
    "idx"    : info.vein_id,
    "auxin_count": 0,
    "dir_count" : 0,
    "orig_dx": dx,
    "orig_dy": dy,
    "dx"    : dx,
    "dy"    : dy
  };

  info.vein_id++;

  return vein_info;
}


g_info.vein_id = 0;
/*
let vein_info = {
  "type"  : "vein",
  "state" : 0,
  "a_idx" : -1,
  "d"     : -1,
  "idx"    : g_info.vein_id,
  "auxin_count": 0,
  "dir_count" : 0,
  "orig_dx": 0,
  "orig_dy": 1,
  "dx"    : 0,
  "dy"    : 1
};
*/
let root_vein = new_vein(g_info, 0, g_info.bb.y, 0, g_info.ds);

g_info.vein.push(root_vein);

//g_info.vein.push( new qt.Point(0, g_info.bb.y, vein_info) );
//g_info.ven_id++;

//tree.insert( g_info.vein );

//let wf = [];
//wf.push( g_info.vein[0] );

let query_kill = new qt.Circle(0, 0, g_info.r_kill);
let query_bias = new qt.Circle(0, 0, g_info.r_bias);

console.log("# r_kill:", g_info.r_kill, ", r_bias:", g_info.r_bias);

let vein_wf = [ g_info.vein[0] ];

let iter = -1;
while (vein_wf.length > 0) {

  iter++;

  if (VERBOSE > 0) {
    console.log("### iter:", iter, ", vein_wf.length:", vein_wf.length);
  }
  //console.log(vein_wf);




  let auxin_queue = [];
  let nxt_vein_wf = [];

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

  // distance query
  //
  for (let ii=0; ii<vein_wf.length; ii++) {
    let v = vein_wf[ii];

    //console.log("vein_wf[", ii, "]:", v);

    query_bias.x = v.x;
    query_bias.y = v.y;
    let querypnt = tree.query(query_bias);
    for (let _k=0; _k<querypnt.length; _k++) {
      let a = querypnt[_k];

      v.auxin_count++;

      let _dx = a.data.x - v.x;
      let _dy = a.data.y - v.y;
      let _dist = Math.sqrt(  _dx*_dx + _dy*_dy );

      //console.log("## dist:", _dx, _dy, _dist);

      if (a.data.v_idx < 0) {
        //console.log(">>> pushing auxin", a.data.idx, a);
        auxin_queue.push(a);

        //console.log("pushed auxin:", a.data.idx);
      }
      else {
        //console.log(" >> auxin", a.data.idx, "already there...");
      }

      if ((a.data.v_idx < 0) ||
          (_dist < a.data.d)) {
        a.data.d = _dist;
        a.data.v_idx = v.idx;

        //console.log(" ....auxin:", a.data.idx, " connected to vein:", v.idx);

        if (VERBOSE > 0) {
          console.log("#> updating auxin (", a.data.x, a.data.y, "),  v_idx:", a.data.v_idx, ", dist:", a.data.d );
        }
        //console.log(a, v);
      }
    }
  }

  //DEBUG
  if (VERBOSE > 0) {
    for (let ii=0; ii<vein_wf.length; ii++) {
      console.log("# vein_wf[", ii, "]: (", vein_wf[ii].x, vein_wf[ii].y, "), auxin_count:", vein_wf[ii].auxin_count );
    }
  }

  //DEBUG
  //for (let ii=0; ii<auxin_queue.length; ii++) {
  //  let a = auxin_queue[ii];
  //  console.log("# auxin_queue[", ii, "/", auxin_queue.length, "]:", a.data.x, a.data.y, ", v_idx:", a.data.v_idx);
  //}

  // go through all auxin nodes that have influence over vein nodes
  // and update their thralled vein nodes
  //
  for (let ii=0; ii<auxin_queue.length; ii++) {

    let a = auxin_queue[ii];
    let v_idx = a.data.v_idx;
    let v = g_info.vein[v_idx];

    //console.log(">>>v_idx:", v_idx);

    let _dx = a.data.x - v.x;
    let _dy = a.data.y - v.y;
    let _dist = Math.sqrt( _dx*_dx + _dy*_dy );

    let ds = g_info.ds;

    //v.orig_dx = v.dx;
    //v.orig_dy = v.dy;
    v.dx += ds*_dx / _dist;
    v.dy += ds*_dy / _dist;
    v.dir_count++;

    //DEBUG
    if (VERBOSE > 0) {
      console.log("# auxin_queue[", ii, "/", auxin_queue.length, "]:", a.data.x, a.data.y, ", v_idx:", a.data.v_idx,
      ", v(", v.x, v.y, ") + (", v.dx, v.dy, ") c:", v.dir_count);
    }


  }

  // remove vein elements that have no auxin nodes within bias distance
  //
  for (let ii=0; ii<vein_wf.length; ii++) {
    let v = vein_wf[ii];
    if (v.auxin_count > 0) { nxt_vein_wf.push( v ); }
    if (v.dir_count > 0) {

      let dx = v.dx / (v.dir_count+1);
      let dy = v.dy / (v.dir_count+1);
      let dst_v = new_vein(g_info, v.x + dx, v.y + dy, dx, dy);

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

for (let ii=0; ii<g_info.auxin.length; ii++) {
  console.log(g_info.auxin[ii].x, g_info.auxin[ii].y);
  console.log("\n");
}
console.log("\n");

for (let ii=0; ii<g_info.vein.length; ii++) {
  let v = g_info.vein[ii];
  console.log(v.x, v.y);
  console.log(v.x + v.dx, v.y + v.dy );
  console.log();
}
//console.log(g_info.vein);

process.exit();


let cpnt = new qt.Circle(0, 0, g_biaso.r_inf);
cpnt.x = 0.0;
cpnt.y = 0.0;
let xx = tree.query(cpnt);

let tpnt = new qt.Circle( g_info.bb.x + g_info.dx/2, g_info.bb.y , g_info.r_bias);
console.log(tpnt, tree.query(tpnt));

//for (let ii=0; ii<xx.length; ii++) {
//  console.log(xx[ii].x, xx[ii].y, xx[ii]);
//}

//let idx = Math.floor( N * Math.random() );
//let wf = [ seq[idx][0], seq[idx][1] ];

//console.log(">>", tree.objects.length, tree.nodes.length);


//console.log(xx);
//console.log("##");

//for (let ii=0; ii<xx.length; ii++) { console.log(xx[ii].x, xx[ii].y); }



