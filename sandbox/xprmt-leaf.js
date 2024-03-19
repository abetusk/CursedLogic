//
// To the extent possible under law, the person who associated CC0 with
// this code has waived all copyright and related or neighboring rights
// to this code.
//
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

//export for commonJS or browser
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



let N = 500;
let S = 10000;
let rho =  S / Math.sqrt(N);

let VERBOSE = 0;

function new_vein(info, x, y, dx, dy, p_idx) {
  p_idx = ((typeof p_idx === "undefined") ? -1 : p_idx);
  let vein_info = {
    "x" : x,
    "y" : y,
    "p_idx": p_idx,
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


//console.log("## N:", N, "S:", S, "rho:", rho);

let g_info = {
  "n": N,
  "bb": {"x":-S/2, "y":-S/2, "w":S, "h": S},
  "dx": S, "dy": S,
  "ds": 200,
  //"r_kill": 0.5*rho,
  //"r_bias": 1.0*rho,
  "r_kill": 1.125*rho,
  //"r_bias": 1.25*rho,
  "r_bias": 2.75*rho,
  "eps": 1/(1024*1024),
  "auxin" : [],
  "vein": [],
  "wf_idx": [],

  "max_iter" : 80,

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

function space_col_init() {
  let tree = new qt.QuadTree( new qt.Box( g_info.bb.x, 
                                          g_info.bb.y,
                                          g_info.bb.w,
                                          g_info.bb.h ) );
  g_info.auxin_tree = tree;

  let halton_cx = 0.5,
      halton_cy = 0.5;


  let seq = halton.seq2d(2,3,g_info.n);

  //seq = [];
  //for (let ii=0; ii< g_info.n; ii++) {
  //  seq.push( [ Math.random(), Math.random() ] );
  //}

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

  tree.insert(g_info.auxin);

}

function space_col_start() {
  __x();
}

function __x() {

  let tree = g_info.auxin_tree;

  // root point
  //
  g_info.vein_id = 0;
  let root_vein = new_vein(g_info, 0, g_info.bb.y, 0, g_info.ds);
  g_info.vein.push(root_vein);

  let query_kill = new qt.Circle(0, 0, g_info.r_kill);
  let query_bias = new qt.Circle(0, 0, g_info.r_bias);

  console.log("# r_kill:", g_info.r_kill, ", r_bias:", g_info.r_bias);

  let vein_wf = [ g_info.vein[0] ];

  // ---disp----
  /*
  let _tp = g_info.auxin_tree.getAllPoints();
  disp_pnt = [];
  for (let ii=0; ii<_tp.length; ii++) {
    disp_pnt.push({
      "c": "#ff0",
      "x": (_tp[ii].x + g_info.disp_dxy[0])*g_info.disp_scale,
      "y": (_tp[ii].y + g_info.disp_dxy[1])*g_info.disp_scale
    });
  }
  _plot( disp_pnt );
  */
  // ---disp----

  let iter = -1;
  while ( (iter < g_info.max_iter) &&
          (vein_wf.length > 0) ) {
  //for (let _iter=0; _iter<60; _iter++) {
  //while (vein_wf.length > 0) {

    //  ---DISPLAY---
    //  ---DISPLAY---
    //  ---DISPLAY---
    disp_pnt = [];

    for (let ii=0; ii<g_info.vein.length; ii++) {
      let c = "#707";
      if (g_info.vein[ii].auxin_count > 0) { c = "#f0f"; }
      disp_pnt.push({
        "c": c,
        "x": g_info.vein[ii].x,
        "y": g_info.vein[ii].y
      });
    }

    //  ---DISPLAY---
    //  ---DISPLAY---
    //  ---DISPLAY---


    iter++;

    //console.log(">>>", iter, vein_wf);

    if (VERBOSE > 0) {
      console.log("### iter:", iter, ", vein_wf.length:", vein_wf.length);
    }
    //console.log(vein_wf);


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

        disp_pnt.push({
          "c": "#f00",
          "x": killpnt[_k].x,
          "y": killpnt[_k].y
        });

        tree.remove( killpnt[_k] );

        //console.log("KILL", killpnt[_k]);
      }
    }

    // ---disp----
    let _tp = g_info.auxin_tree.getAllPoints();
    for (let ii=0; ii<_tp.length; ii++) {
      disp_pnt.push({
        "c": "#070",
        "x": _tp[ii].x,
        "y": _tp[ii].y
      });
    }
    _transform(disp_pnt, g_info.disp_dxy[0], g_info.disp_dxy[1], g_info.disp_scale);
    //_plot( disp_pnt );
    // ---disp----


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
        let dst_v = new_vein(g_info, v.x + dx, v.y + dy, dx, dy, v.idx);

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

  //_plot( disp_pnt );

  console.log("done");

  /*
  disp_pnt = [];
  for (let ii=0; ii<g_info.vein.length; ii++) {
    disp_pnt.push({
      "c": "f0f",
      "x": g_info.vein[ii].x,
      "y": g_info.vein[ii].y
    });
  }
  _transform(disp_pnt, g_info.disp_dxy[0], g_info.disp_dxy[1], g_info.disp_scale);
  _plot(disp_pnt);
  */

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

  /*
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
  */


  return;

  /*
  let cpnt = new qt.Circle(0, 0, g_biaso.r_inf);
  cpnt.x = 0.0;
  cpnt.y = 0.0;
  let xx = tree.query(cpnt);

  let tpnt = new qt.Circle( g_info.bb.x + g_info.dx/2, g_info.bb.y , g_info.r_bias);
  console.log(tpnt, tree.query(tpnt));
  */
}

//for (let ii=0; ii<xx.length; ii++) {
//  console.log(xx[ii].x, xx[ii].y, xx[ii]);
//}

//let idx = Math.floor( N * Math.random() );
//let wf = [ seq[idx][0], seq[idx][1] ];

//console.log(">>", tree.objects.length, tree.nodes.length);


//console.log(xx);
//console.log("##");

//for (let ii=0; ii<xx.length; ii++) { console.log(xx[ii].x, xx[ii].y); }

space_col_init();


