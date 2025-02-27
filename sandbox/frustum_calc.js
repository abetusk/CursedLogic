// LICENSE: CC0
//
// To the extent possible under law, the person who associated CC0 with
// this project has waived all copyright and related or neighboring rights
// to this project.
// 
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

var njs = require("./numeric.js");

function print_e(p,q) {
  console.log(p[0], p[1], p[2]);
  console.log(q[0], q[1], q[2], "\n\n")
}

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

function v3theta(p,q) {
  let s = njs.norm2( cross3(p,q) );
  let c = njs.dot(p,q);
  return Math.atan2(s,c);
}

function plane_f(u, Np, p) {
  return njs.dot( Np, njs.sub(u, p) );
}

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

function t_plane_line(Np, p, v0, v) {
  let _eps = 1/(1024*1024*1024);
  let _d = njs.dot(Np,v);
  if (Math.abs(_d) < _eps) { return NaN; }

  let t = (njs.dot(Np,p) - njs.dot(Np,v0)) / _d;
  return t;
}

function Vt( v0, v, t ) {
  return njs.add(v0, njs.mul(t, v));
}

function cross3(p,q) {
  let c0 = ((p[1]*q[2]) - (p[2]*q[1])),
      c1 = ((p[2]*q[0]) - (p[0]*q[2])),
      c2 = ((p[0]*q[1]) - (p[1]*q[0]));

  return [c0,c1,c2];
}

function p3toP(p0, p1, p2) {
  let _eps = 1/(1024*1024*1024);

  let p10 = njs.sub(p1,p0);
  let p20 = njs.sub(p2,p0);

  let Np = cross3(p10,p20);

  let Pk = -njs.dot(Np, p0);

  return [ Np[0], Np[1], Np[2], Pk ]
}

function debug_shell(q, frustum_v) {

  let Nq = njs.mul( 1/njs.norm2(q), q );
  let p0 = q;

  let idir_descr = [ "+x", "-x", "+y", "-y", "+z", "-z" ];

  for (let D=1; D<3.1; D+=1) {
    for (let idir=0; idir<6; idir++) {

      for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
        let v_cur = njs.mul(D, frustum_v[idir][f_idx]);
        console.log(v_cur[0], v_cur[1], v_cur[2]);
      }
      console.log("\n\n");

    }
  }

  for (let idir=0; idir<6; idir++) {
    let _n = frustum_v[idir].length;
    for (let D=1; D<3.1; D+=1) {
      for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
        let v_cur = njs.mul(D, frustum_v[idir][f_idx]);
        let v_nxt = njs.mul(D, frustum_v[idir][(f_idx+1)%_n]);

        let vnc = njs.sub(v_nxt, v_cur);

        let t = t_plane_line(Nq, p0, v_cur, vnc);

        let _t_idir = v2idir(vnc);

        console.log("#D:", D, "idir:", idir, "f_idx:", f_idx, "t:", t, "(f_dir:", idir_descr[_t_idir], ")");

        if ((t > 0) && (t < 1)) {
          let u = Vt(v_cur, vnc, t);
          console.log(u[0], u[1], u[2]);
        }
      }

      console.log("\n\n");
    }
  }



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
//   frame_d    : D values (e.g. [2,4])
//   frame_v    : frame vectors (for D = frame_d[idx])
//   frame_t    : frame time (q-plane intersect to frustum frame for scale D)
// }
//
function frustum3d_intersection(q, ds) {
  ds = ((typeof ds === "undefined") ? 1 : ds);
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


  //DEBUG
  let debug_frustum = true;
  if (debug_frustum) {

    debug_shell(q, frustum_v);

    console.log(q[0], q[1], q[2], "\n");
    for (let idir=0; idir<6; idir++) {
      for (let f_idx=0; f_idx<frustum_v[idir].length; f_idx++) {
        print_e([0,0,0], frustum_v[idir][f_idx]);
      }
    }

    console.log("\n\n");

    for (let idir=0; idir<6; idir++) {
      for (let f_idx=0; f_idx<frustum_v[idir].length; f_idx++) {
        let v = frustum_v[idir][f_idx];
        console.log(v[0], v[1], v[2]);
      }
      console.log("\n\n");
    }


    /*
    let qq = njs.dot( n_q, q );
    for (let it=0; it<1000; it++) {
      let ux = 2*(Math.random());
      let uy = 2*(Math.random()-0.5);
      let uz = (qq - n_q[0]*ux - n_q[1]*uy) / n_q[2];

      console.log(ux, uy, uz, "\n");
    }
    */
  }
  //DEBUG

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

    // frustum idir check
    //
    // wip
    //
    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v = frustum_v[idir][f_idx];
    }
    //
    // frustum idir check


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

    found_idir = idir;

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

    // frustum idir check
    //
    // wip
    //

    let part_count = 0;
    let _n = frustum_v[idir].length;
    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v_cur = frustum_v[idir][f_idx];
      let v_nxt = frustum_v[idir][(f_idx+1)%_n];
      let vv = cross3(v_cur, v_nxt);
      if (njs.dot(vv, q) >= 0) { part_count++; }

      console.log("## f_idx:", f_idx, "v_cur:", v_cur, "v_nxt:", v_nxt, "vv:", vv, "vv.q:", njs.dot(vv,q), "(part_count now", part_count, ")");

    }

    console.log("## part_count:", part_count, "(/", _n, ")");

    if (part_count == _n) { frustum_idir = idir; }

    //
    // frustum idir check

  }

  let frame_d = [ 2, 4];
  let frame_t = [ [ 0,0,0,0 ], [0,0,0,0] ];
  let frame_v = [ [], [] ];
  if (frustum_idir>=0) {

    let idir = frustum_idir;
    let _n = frustum_v[idir].length;

    let Nq = njs.mul( 1 / njs.norm2(q), q );

    let D = 2;


    for (let f_idx=0; f_idx < frustum_v[idir].length; f_idx++) {
      let v_cur = frustum_v[idir][f_idx];
      let v_nxt = frustum_v[idir][(f_idx+1)%_n];
      let vv = njs.sub(v_nxt, v_cur);

      let vd_cur = njs.mul(D, v_cur);
      let vd_nxt = njs.mul(D, v_nxt);
      let vvd = njs.sub(vd_nxt, vd_cur);

      let t1 = t_plane_line( Nq, q, v_cur, vv );
      let td = t_plane_line( Nq, q, vd_cur, vvd );

      //t1 /= njs.norm2(vv);
      //td /= njs.norm2(vvd);

      frame_t[0][f_idx] = t1;
      frame_t[1][f_idx] = td;

      frame_d[0] = njs.norm2(vv);
      frame_d[1] = njs.norm2(vvd);

      frame_v[0].push( [ v_cur, vv ] );
      frame_v[1].push( [ vd_cur, vvd ] );
    }

  }


  return {
    "idir": found_idir,
    "idir_t": _res_t,
    "frustum_t": _frustum_t,
    "frustum_v": frustum_v,

    "frustum_idir": frustum_idir,
    "frame_v": frame_v,
    "frame_d": frame_d,
    "frame_t" : frame_t
  };

}

function __analyze(info, q) {

  let Nq = njs.mul( 1/njs.norm2(q), q );
  let Pq = q;

  let idir_descr = [ "+x", "-x", "+y", "-y", "+z", "-z" ];

  let frame_v = info.frame_v;
  let famee_d = info.frame_d;
  let frame_t = info.frame_t;

  let idir = info.frustum_idir;

  let okvec = [0,0,0,0];

  //console.log("point in idir:", info.frustum_idir);

  for (let i=0; i<frame_t[0].length; i++) {
    let _f = frame_t[1][i] / frame_t[0][i];

    let frame_idir = v2idir(frame_v[0][i][1]);
    //console.log(" frame_idir:", frame_idir, idir_descr[frame_idir], ">>>", (_f > 1) ? "+++" : "---" );

    //console.log(frame_v[0][i][0]);
    let xxx = plane_f( frame_v[0][i][0], Nq, Pq );

    console.log("# xxx:", xxx, "_f:", _f);

    if ((xxx > 0) && ( _f > 1)) {
      okvec[i] = 1;
      //console.log("  away (ok)");
    }
    else if ((xxx < 0) && (_f < 1)) {
      okvec[i] = 1;
      //console.log("  closer (ok)");
    }
    else {
      okvec[i] = 0;
      //console.log("  no");
    }
  }

  return okvec;
}

function _rnd3C() {
  return [
    2*(Math.random()-0.5),
    2*(Math.random()-0.5),
    2*(Math.random()-0.5)
  ];
}

function comment_stringify(data) {
  let lines = JSON.stringify(data, undefined, 2).split("\n");
  for (let i=0; i<lines.length; i++) {
    console.log("#", lines[i]);
  }
}

function investigate_q_point() {
  let q;
  q =  [ 0.4608165114850644, 0.21948347420131942, 0.24588673712113795 ];
  q =  [ 0.99, 0.01, 0.97 ];

  //console.log("#", q);

  console.log(0,0,0);
  console.log(q[0], q[1], q[2], "\n\n");

  let res = frustum3d_intersection(q);


  comment_stringify(res);

  let __a = __analyze(res, q);

  console.log("#a:", __a);

  return;

  /*
  console.log(res);
  console.log( JSON.stringify( res.frame_v[0] ) );
  console.log( JSON.stringify( res.frame_v[1] ) );
  */

  let v0 = [ 1, 1, 1];
  let v1 = [ 1,-1, 1];
  let v2 = [-1, 1, 1];
  let v3 = [ 1, 1,-1];

  let u0 = [-1,-1, 1];

  let v01 = njs.sub(v1,v0);
  let v02 = njs.sub(v2,v0);
  let v03 = njs.sub(v3,v0);

  let Nv0 = njs.mul( 1/njs.norm2(v0), v0 );
  let Nv1 = njs.mul( 1/njs.norm2(v1), v1 );
  let Nv2 = njs.mul( 1/njs.norm2(v2), v2 );
  let Nv3 = njs.mul( 1/njs.norm2(v3), v3 );

  let Nu0 = njs.mul( 1/njs.norm2(u0), u0 );

  let Nv01 = njs.mul( 1/njs.norm2(v01), v01 );
  let Nv02 = njs.mul( 1/njs.norm2(v02), v02 );
  let Nv03 = njs.mul( 1/njs.norm2(v03), v03 );

  let _s3 = Math.sqrt(3);

  console.log("Nv0:", Nv0);
  let _t = njs.mul( -1/_s3, njs.add( Nv01, njs.add(Nv02, Nv03) ) );
  console.log(" _t:", _t);

  console.log(">>> ang(Nv01, Nv02):", v3theta(Nv01, Nv02));
  console.log(">>> ang(Nv01, Nv03):", v3theta(Nv01, Nv03));
  console.log(">>> ang(Nv02, Nv03):", v3theta(Nv02, Nv03));

  let a = v3theta(Nv0, Nv01);
  console.log(">>> ang(Nv0, Nv01):", a, Math.PI - a );

  console.log(">>> ang(Nv0, Nv02):", v3theta(Nv0, Nv02));
  console.log(">>> ang(Nv0, Nv03):", v3theta(Nv0, Nv03));

  a = v3theta(Nv0, Nu0);
  console.log(">>> ang(Nv0, Nu0):", a, Math.PI - a );

  console.log(">>> ang(Nv0, Nv1)", v3theta(Nv0, Nv1));
  console.log(">>> ang(Nv0, Nv3)", v3theta(Nv0, Nv3));

  __analyze(res, q);

  return;

}

function experiment_x() {

  for (let it=0; it < 10000; it++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    let _three_cut=0;

    for (let idir=0; idir<6; idir++) {
      let _count = 0;
      for (let ii=0; ii<4; ii++) {
        if (res.frustum_t[idir][ii] > 0) { _count++; }
      }

      if (_count==3) {
        _three_cut ++;
      }
    }

    if (_three_cut > 0) {
      let ok = __analyze(res, q);
      if ((ok[0] == 0) ||
          (ok[1] == 0) ||
          (ok[2] == 0) ||
          (ok[3] == 0)) {
        console.log(q[0], q[1], q[2]);
      }
    }

  }

}

function full_cut_square_region() {
  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    if (res.idir >= 0) {

      if ( (res.frustum_t[res.idir][0] < 1) &&
           (res.frustum_t[res.idir][1] < 1) &&
           (res.frustum_t[res.idir][2] < 1) &&
           (res.frustum_t[res.idir][3] < 1) ) {
        console.log(q[0], q[1], q[2]);
        c++;
      }
    }
  }

  console.log("#", c / N);

}

function four_cut_region() {
  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    let _four_cut = 0;

    for (let idir=0; idir<6; idir++) {
      let _count = 0;
      for (let ii=0; ii<4; ii++) {
        if (res.frustum_t[idir][ii] > 0) { _count++; }
      }

      if (_count==4) {
        _four_cut ++;
      }
    }

    if (_four_cut > 0) {
      console.log(q[0], q[1], q[2]);
      c++;
    }

  }

  console.log("#", c / N);


}

function three_cut_region() {
  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    let _three_cut = 0;

    for (let idir=0; idir<6; idir++) {
      let _count = 0;
      for (let ii=0; ii<4; ii++) {
        if (res.frustum_t[idir][ii] > 0) { _count++; }
      }

      if (_count==3) {
        _three_cut ++;
        //console.log(q[0], q[1], q[2]);
        //c++;
      }
    }

    if (_three_cut > 0) {
      console.log(q[0], q[1], q[2]);
      c++;
    }

  }

  console.log("#", c / N);

}

function mult_three_cut_region() {
  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    let _three_cut = 0;

    let idir_cut = [0,0,0,0,0,0];

    for (let idir=0; idir<6; idir++) {
      let _count = 0;
      for (let ii=0; ii<4; ii++) {
        if (res.frustum_t[idir][ii] > 0) { _count++; }
      }

      if (_count==3) {
        _three_cut ++;
        idir_cut[idir] = 1;
        //console.log(q[0], q[1], q[2]);
        //c++;
      }
    }

    if (_three_cut == 3) {

      if (idir_cut[0] && idir_cut[2] && idir_cut[4]) {
        console.log(q[0], q[1], q[2]);
        c++;
      }
    }

  }

  console.log("#", c / N);

}

// should be zer0?
// unless we get a point exactly on the frustum planes...
//
function only_two_cut_region() {
  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    let _two_count = 0;
    for (let idir=0; idir<6; idir++) {
      let _count = 0;
      for (let ii=0; ii<4; ii++) {
        if (res.frustum_t[idir][ii] > 0) { _count++; }
      }

      if (_count==2) { _two_count++; }
    }

    if (_two_count == 6) {
      console.log(q[0], q[1], q[2]);
      c++;
    }
  }

  console.log("#", c / N);

}

function main() {

  //experiment_x();
  //return;

  investigate_q_point();
  return;

  //let q = _rnd3C();

  //full_cut_square_region();
  //return;

  //four_cut_region();
  //return;

  //three_cut_region();
  //return;

  mult_three_cut_region();
  return;

  only_two_cut_region();
  return;


  let N = 100000;
  let c= 0;
  for (let i=0; i<N; i++) {
    let q = _rnd3C();
    let res = frustum3d_intersection(q);

    console.log(res);
    return;

    if (res.idir >= 0) {

      if ( (res.frustum_t[res.idir][0] < 1) &&
           (res.frustum_t[res.idir][1] < 1) &&
           (res.frustum_t[res.idir][2] < 1) &&
           (res.frustum_t[res.idir][3] < 1) ) {
        console.log(q[0], q[1], q[2]);
        c++;
      }
    }
  }

  console.log("#", c / N);

  //console.log(res);


}

main();

