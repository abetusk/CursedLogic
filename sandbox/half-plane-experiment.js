
var njs = require("./numeric.js");

var _faces = [
  [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1]],
  [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]],
  [[0, 0, 0], [0, 0, 1], [1, 0, 0], [1, 0, 1]],
  [[0, 1, 0], [0, 1, 1], [1, 1, 0], [1, 1, 1]],
  [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]],
  [[0, 0, 1], [0, 1, 1], [1, 0, 1], [1, 1, 1]]
];

// 0 : -x
// 1 : +x
// 2 : -y
// 3 : +y
// 4 : -z
// 5 : +z
//
var faces = [
  [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0] ],
  [[1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0] ],
  [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0] ],
  [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0] ],
  [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0] ],
  [[0, 0, 1], [0, 1, 1], [1, 1, 1], [1, 0, 1] ]
];

function rand_normal() {
  let pi  = 2 * Math.PI * Math.random();
  let R   = Math.sqrt(-2 * Math.log(Math.random()));
  let x   = R * Math.cos(pi);
  let y   = R * Math.sin(pi);

  return [ x, y ];
}

function rand_sphere() {
  let sqrt3 = Math.cbrt(Math.random());
  let u01 = rand_normal();
  let u23 = rand_normal();

  let v = [ u01[0], u01[1], u23[0] ];
  let d = njs.norm2(v);

  return njs.mul(sqrt3, njs.mul(1/d, v));
}

function _debug_randsphere() {
  for (let i=0; i<1000; i++) {
    let v = rand_sphere();
    console.log(v[0], v[1], v[2]);
  }
  process.exit();
}

function _debug_randsphere1() {
  for (let i=0; i<1000; i++) {
    let v = rand_sphere();
    let vn = njs.mul( 1/njs.norm2(v), v);

    console.log(vn[0], vn[1], vn[2]);
  }
  process.exit();
}



function line_plane_intersect(v0, vn, p0, pn) {
  let denom = njs.dot( vn, pn );
  let _eps = 1/(1024*1024);

  if (denom < _eps) { return { "T":"parallel", "t": 0, "p":[0,0,0] } }

  let t = njs.dot( njs.sub(p0,v0), pn ) / denom;

  return {"T": "x", "t": t, "p": njs.add(v0, njs.mul(t, vn))};
}

function face_plane_intersect(face_line, p0, pn) {

  for (let i=0; i<face_line.length; i++) {
    let info = line_plane_intersect(face_line[i][0], face_line[i][1], p0, pn);
    if (info.T == "parallel") { continue; }

    if ((info.t > 0) && (info.t < 1)) { return true; }
  }

  return false;
}

function rand_plane_point(p0, pn) {
  let x = (Math.random() - 0.5)*2;
  let y = (Math.random() - 0.5)*2;

  let z = (njs.dot(p0, pn) - (x*pn[0]) - (y*pn[1])) / pn[2];

  return [x,y,z];
}

let face_lines = [ [], [], [], [], [], [] ];

for (let face_idx=0; face_idx < faces.length; face_idx++) {
  for (let i=0; i<faces[face_idx].length; i++) {
    let i_nxt = (i+1)%(faces[face_idx].length);
    let midpoint = [0.5, 0.5, 0.5];
    face_lines[face_idx].push( [ njs.sub(faces[face_idx][i], midpoint), njs.sub(faces[face_idx][i_nxt], faces[face_idx][i]) ] );
  }
}

function experiment0() {
  let stat_count  = [0,0, 0,0, 0,0];
  let stat_tot    = [0,0, 0,0, 0,0];

  for (let it=0; it<100000; it++) {
    let pn = [ Math.random() - 0.5, Math.random() - 0.5, Math.random() - 0.5 ];
    let _n = njs.norm2(pn);
    pn = njs.mul( 1/_n, pn );

    //console.log("#", pn);

    for (let face_idx=0; face_idx < 6; face_idx++) {

      let r = face_plane_intersect(face_lines[face_idx], [0,0,0], pn);

      if (!r) { stat_count[face_idx]++; }
      stat_tot[face_idx]++;

    }

  }

  for (let i=0; i<6; i++) {
    console.log("#", i, stat_count[i], stat_tot[i]);
    console.log(i, stat_count[i] / stat_tot[i]);
  }
}

function experiment1() {
  let stat_count  = [];
  let stat_tot    = 0;

  let n_it = 100000;

  let T = 100;

  for (let i=0; i<=T; i++) { stat_count.push(0); }

  for (let it=0; it<n_it; it++) {

    let face_hit = [ 0,0, 0,0, 0,0 ];

    let t = 0;
    while (t<T) {
      //let pn = [ Math.random() - 0.5, Math.random() - 0.5, Math.random() - 0.5 ];
      let pn = rand_sphere();
      let _n = njs.norm2(pn);
      pn = njs.mul( 1/_n, pn );

      for (let face_idx=0; face_idx < 6; face_idx++) {

        let r = face_plane_intersect(face_lines[face_idx], [0,0,0], pn);

        if (!r) { face_hit[face_idx] = 1; }
      }

      let s = 0;
      for (let i=0; i<face_hit.length; i++) { s += face_hit[i]; }
      if (s == 6) { break; }

      t++;
    }

    stat_count[t]++;
    stat_tot++;

  }

  for (let t=0; t<stat_count.length; t++) {
    if (stat_count[t] > 0) {
      console.log(t, stat_count[t] / stat_tot);
    }

  }

}

function experiment2() {

  let n_it = 1000000;
  let ds = 1/32;

  for (let d=0; d<0.25; d+=ds) {
    let stat_count  = [];
    let stat_tot    = 0;


    let T = 100;

    for (let i=0; i<=T; i++) { stat_count.push(0); }

    for (let it=0; it<n_it; it++) {

      let face_hit = [ 0,0, 0,0, 0,0 ];

      let t = 0;
      while (t<T) {
        let p0 = njs.mul(d, rand_sphere());
        let pn = rand_sphere();
        let _n = njs.norm2(pn);
        pn = njs.mul( 1/_n, pn );

        for (let face_idx=0; face_idx < 6; face_idx++) {
          let r = face_plane_intersect(face_lines[face_idx], p0, pn);
          if (!r) { face_hit[face_idx] = 1; }
        }

        let s = 0;
        for (let i=0; i<face_hit.length; i++) { s += face_hit[i]; }
        if (s == 6) { break; }

        t++;
      }

      stat_count[t]++;
      stat_tot++;

    }

    let E = 0;
    for (let t=0; t<stat_count.length; t++) {
      if (stat_count[t] > 0) {
        E += t*(stat_count[t] / stat_tot)
        console.log("#", d, t, stat_count[t] / stat_tot);
        console.log(t, stat_count[t] / stat_tot);
      }
    }
    console.log("# E:", E);
    console.log("\n");

  }
}

function experiment3(d) {
  d = ((typeof d === "undefined") ? 0.25 : d);

  let n_it = 100000;
  let ds = 1/32;

  let stat_count  = [];
  let stat_tot    = 0;


  let T = 100;

  for (let i=0; i<=T; i++) { stat_count.push(0); }

  for (let it=0; it<n_it; it++) {

    let face_hit = [ 0,0, 0,0, 0,0 ];

    let t = 0;
    while (t<T) {
      let p0 = njs.mul(d, rand_sphere());
      let _n = njs.norm2(p0);
      let pn = njs.mul( 1/_n, p0 );

      for (let face_idx=0; face_idx < 6; face_idx++) {
        let r = face_plane_intersect(face_lines[face_idx], p0, pn);
        if (!r) { face_hit[face_idx] = 1; }
      }

      let s = 0;
      for (let i=0; i<face_hit.length; i++) { s += face_hit[i]; }
      if (s == 6) { break; }

      t++;
    }

    stat_count[t]++;
    stat_tot++;

  }

  let E = 0;
  for (let t=0; t<stat_count.length; t++) {
    if (stat_count[t] > 0) {
      E += t*(stat_count[t] / stat_tot)
      console.log("#", d, t, stat_count[t] / stat_tot);
      console.log(t, stat_count[t] / stat_tot);
    }
  }
  console.log("# E:", E);
  console.log("\n");
}

//experiment3();
//process.exit();

function _cruft() {
  for (let i=0; i<100; i++) {
    let p = rand_plane_point([0,0,0], pn);
    console.log(p[0], p[1], p[2]);
    console.log("\n");
  }

  for (let face_idx = 0; face_idx < 6; face_idx++) {
    let r = face_plane_intersect(face_lines[face_idx], [0,0,0], pn);
    console.log("#face:", face_idx, "plane_n:", pn, "res:", r);
  }

  process.exit();

  let curface = 0;
  for (let i=0; i<face_lines[curface].length; i++) {

    let f = 1/32;
    let dxy = [ f*Math.random(), f*Math.random(), f*Math.random() ];

    let v0 = njs.add( njs.clone(face_lines[curface][i][0]), dxy);
    let vn = njs.add( njs.clone(face_lines[curface][i][1]), dxy);

    let vt = njs.add( v0, njs.mul(1, vn) );

    console.log( v0[0], v0[1], v0[2] );
    console.log( vt[0], vt[1], vt[2] );
    console.log("\n");
  }

  process.exit();
}


function _cruft1() {
  let perm = [0,1,3,2];

  for (let i=0; i<faces.length; i++) {
    for (let j=0; j<faces[i].length; j++) {
      let k = perm[j];
      console.log(faces[i][k][0],faces[i][k][1],faces[i][k][2]);
    }
    console.log(faces[i][0][0],faces[i][0][1],faces[i][0][2]);
    console.log("\n");
  }
}
