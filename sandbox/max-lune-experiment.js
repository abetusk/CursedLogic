

var njs = require("./numeric.js");

function rnd_vec(dim) {
  dim = ((typeof dim === "undefined") ? 3: dim);
  let u = [];
  for (let i=0; i<dim; i++) {
    u.push(Math.random());
  }
  return u;
}

function cross3(u,v) {

  let M = [
    [    0, -u[2],  u[1] ],
    [ u[2],     0, -u[0] ],
    [-u[1],  u[0],     0 ]
  ];

  return njs.dot(M, v);
}

// euler rotation or olinde rodrigues
// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
//
function rodrigues(v0, _vr, theta) {
  let c = Math.cos(theta);
  let s = Math.sin(theta);

  let v_r = njs.mul( 1 / njs.norm2(_vr), _vr );

  return njs.add(
    njs.mul(c, v0),
    njs.add(
      njs.mul( s, cross3(v_r,v0)),
      njs.mul( (1-c) * njs.dot(v_r, v0), v_r )
    )
  );
}


let v_r = [ Math.random(), Math.random(), Math.random() ];
let v0 = [ Math.random(), Math.random(), Math.random() ];

//v_r = njs.mul( 1/ njs.norm2(v_r), v_r );

console.log(0,0,0);
console.log(v0[0], v0[1], v0[2], "\n");

console.log(0,0,0);
console.log(v_r[0], v_r[1], v_r[2], "\n");

for (let idx = 0; idx <= 100; idx++ ) {
  let theta = 2*Math.PI*idx / 100;
  theta *= 1/25;
  let u = rodrigues(v0,v_r,theta);
  console.log( u[0], u[1], u[2] );
}

process.exit();


let R = 3;

let p = [ Math.random(), Math.random(), Math.random() ];
p = njs.mul( R/njs.norm2(p), p );

let p2 = njs.norm2Squared(p);
let _d = ((2*p[2]*p[2]) - p2) + (R*R*p2);

_d = (p2*p2) - (4*p[2]*R);

console.log("#", _d);

let Mp = 0;

if (_d > 0) {

  /*
  Mp = [
    ( ((2*p[2]*p[2]) - njs.norm2Squared(p)) + Math.sqrt(_d) ) / (2*p[2]),
    ( ((2*p[2]*p[2]) - njs.norm2Squared(p)) - Math.sqrt(_d) ) / (2*p[2])
  ];
  */

  Mp = [
    (p2 + Math.sqrt(_d)) / (2*p[2]*p[2]),
    (p2 + Math.sqrt(_d)) / (2*p[2]*p[2]),
    (p2 - Math.sqrt(_d)) / (2*p[2]*p[2])
  ];

  Mp = [
    p[2] * R / Math.sqrt( p2 - (p[2]*p[2]) ),
    p[2] * R / Math.sqrt( p2 - (p[2]*p[2]) )
  ];

  console.log("# Mp:", Mp);
}

for (let ii=0; ii<300; ii++) {
  let _u = njs.sub( rnd_vec(), [0.5, 0.5, 0.5] );
  let u = njs.mul( R / njs.norm2(_u), _u );

  console.log(u[0], u[1], u[2], "\n");
  console.log(u[0] + p[0], u[1] + p[1], u[2] + p[2], "\n");

  let xyz0 = [ Math.random()*4*R + p[0], Math.random()*4*R + p[1], Mp[0] ];
  let xyz1 = [ Math.random()*4*R + p[0], Math.random()*4*R + p[1], Mp[1] ];

  console.log(xyz0[0], xyz0[1], xyz0[2], "\n");
  console.log(xyz1[0], xyz1[1], xyz1[2], "\n");
}




