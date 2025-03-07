//import * as THREE from './three.js';
import * as THREE from 'three';

import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import * as GeometryUtils from 'three/addons/utils/GeometryUtils.js';

import { FontLoader } from 'three/addons/loaders/FontLoader.js';


var _numeric = await import("./njs.mjs");
var njs = _numeric.Numeric;

//------
//------
//------

// plane equation:
//
// P(u) = N_p \dot (u - p)
//
function plane_f(u, Np, p) {
  return njs.dot( Np, njs.sub(u, p) );
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





// p0 : point on plane
// u  : normal to plane
// box_r: frustum scaling factor (default 1)
//        radius of box (distance of side of frustum box to p)
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
function frustum3d_intersection(p_world, q_world, box_r) {
  box_r = ((typeof box_r === "undefined") ? 1 : box_r);

  let q = njs.sub(q_world, p_world);

  let s3 = 1/Math.sqrt(3);
  let s3br = s3*box_r;

  let L = box_r;
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
  for (let idir=0; idir<6; idir++) {
    let fv_count = 0;
    let fv_n = frustum_v[idir].length;

    // test to see if there are the four frustum vectors that
    // have positive 'time' intersection to q-plane
    //
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
    //

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

  // now calculate the 'time' component of intersection
  // to the frame for the frustum the q point sits fully
  // inside
  //
  // First calculate which frustum idir the q point sits in
  // by seeing if q is comletely enclosed by four planes
  // making up the frustum in this idir.
  //
  let frustum_idir=-1;
  for (let idir = 0; idir < 6; idir++) {
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

  // once the frustum idir of where q is sitting in is
  // determined, find the times for the q-plane to each
  // of the frustum edge frame intersections.
  // `frame_sd_t` holds the source/dest times for the frame
  // edge, where which source or destination is filled
  // depends on which side the source of the frustum frame
  // edge (v_cur) sits on the q-plane.
  //
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

function simple_line(u,v,c) {
  c = ((typeof c === "undefined") ? [Math.random(), Math.random(), Math.random() ] : c );

  let matLine = new LineMaterial( {
    color: 0xffffff,
    linewidth: 2,
    vertexColors: true,
    dashed: false,
    alphaToCoverage: true
  } );


  let points = [];
  let colors = [];
  points.push( u[0], u[1], u[2] );
  points.push( v[0], v[1], v[2] );
  colors.push( c[0], c[1], c[2] );
  colors.push( c[0], c[1], c[2] );

  let g = new LineGeometry();
  g.setPositions( points );
  g.setColors( colors );

  let l = new Line2( g, matLine );
  l.computeLineDistances();
  l.scale.set( 1, 1, 1 );
  scene.add( l );


}

function setup_viz() {
  let v_idir = [
    [1,0,0], [-1,0,0],
    [0,1,0], [0,-1,0],
    [0,0,1], [0,0,-1]
  ];


  let p = [ 0.15087696063772743, 0.6433609372904464, 0.14637496046413653 ];
  let q = [ 0.17501859380897827, 0.6267461187478763, 0.4539898881747171 ];


  let grid_n = 3;
  let ds = 1 / grid_n;

  let Wp = [ p[0]*grid_n, p[1]*grid_n, p[2]*grid_n ];
  let ip = Wp.map( Math.floor );

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

  let ir = 1;
  let frustum_box_r = l0 + (ds*ir);

  let fi_info = frustum3d_intersection(p, q, frustum_box_r);

  console.log(fi_info);

  let color = new THREE.Color();
  let matLine = new LineMaterial( {

    color: 0xffffff,
    linewidth: 2, // in world units with size attenuation, pixels otherwise
    vertexColors: true,

    dashed: false,
    alphaToCoverage: true,

  } );

  let xyz_legend = [
    [ [-0.5, -0.5, -0.5], [-0.25, -0.5, -0.5], [0.5, 0, 0] ],
    [ [-0.5, -0.5, -0.5], [-0.5, -0.25, -0.5], [0, 0.5, 0] ],
    [ [-0.5, -0.5, -0.5], [-0.5, -0.5, -0.25], [0, 0, 0.5] ]
  ];

  for (let xyz=0; xyz<3; xyz++) {
    simple_line( xyz_legend[xyz][0], xyz_legend[xyz][1], xyz_legend[xyz][2] );
  }

  simple_line( p, q, [1,1,1] );
  simple_line( p, njs.add(p, njs.mul(l0, v_idir[p_near_idir])), [1,0,1] );

  console.log(p_near_idir);

  let _shift = 1/512;

  for (let idir=0; idir<6; idir++) {

    color.setHSL( Math.random(), 1.0, 0.5, THREE.SRGBColorSpace );
    for (let v_idx=0; v_idx<fi_info.frustum_v[idir].length; v_idx++) {
      let pp = njs.add(p, njs.mul(_shift, v_idir[idir]));
      let tv = njs.add(pp, fi_info.frustum_v[idir][v_idx]);
      simple_line( pp, tv, [color.r, color.g, color.b ] );
    }

    let frame_points = [];
    let frame_colors = [];

    color.setHSL( Math.random(), 1.0, 0.5, THREE.SRGBColorSpace );
    for (let v_idx=0; v_idx<fi_info.frustum_v[idir].length; v_idx++) {
      let pp = njs.add(p, njs.mul(_shift, v_idir[idir]));
      let tv = njs.add(pp, fi_info.frustum_v[idir][v_idx]);
      frame_points.push( tv[0], tv[1], tv[2] );
      frame_colors.push(color.r, color.g, color.b);
    }
    let pp = njs.add(p, njs.mul(1/128, v_idir[idir]));
    let tv = njs.add(pp, fi_info.frustum_v[idir][0]);
    frame_points.push( tv[0], tv[1], tv[2] );
    frame_colors.push(color.r, color.g, color.b);

    let fg = new LineGeometry();
    fg.setPositions( frame_points );
    fg.setColors( frame_colors );

    let fl = new Line2( fg, matLine );
    fl.computeLineDistances();
    fl.scale.set( 1, 1, 1 );
    scene.add( fl );

  }


}

//------
//------
//------

export var g_data = {
  "njs": njs,
  "camera": "ortho"
};

let camera_persp,
    camera_ortho;

let line, renderer, scene;
let controls_ortho,
    controls_persp;
let line1;
let matLine, matLineBasic, matLineDashed;

// viewport
let insetWidth;
let insetHeight;


const pointer = new THREE.Vector2();
const frustumSize = 2;

init();

function simple_text(txt, pos, scale, color) {
  pos = ((typeof pos === "undefined") ? [0,0,0] : pos );
  scale = ((typeof scale === "undefined") ? 0.125 : scale );
  color = ((typeof color === "undefined") ? 0x006699 : color );

  let font = g_data["font"];

  const matDark = new THREE.LineBasicMaterial( {
    color: color,
    side: THREE.DoubleSide
  });

  const matLite = new THREE.MeshBasicMaterial( {
    color: color,
    transparent: true,
    opacity: 0.4,
    side: THREE.DoubleSide
  });

  const message = txt;
  const shapes = font.generateShapes( message, scale );
  const geometry = new THREE.ShapeGeometry( shapes );
  geometry.computeBoundingBox();
  geometry.translate( pos[0], pos[1], pos[2] );

  const text = new THREE.Mesh( geometry, matLite );
  scene.add( text );
}

function init() {

  const loader = new FontLoader();
  loader.load( 'fonts/helvetiker_regular.typeface.json', function ( font ) {
    g_data["font"] = font;
    simple_text('z', [-0.65, -0.55, -0.2 ], 0.125 );
    simple_text('x', [-0.2, -0.55, -0.5 ], 0.125 );
    simple_text('y', [-0.5, -0.2, -0.5 ], 0.125 );
  });

  renderer = new THREE.WebGLRenderer( { antialias: true } );
  renderer.setPixelRatio( window.devicePixelRatio );
  renderer.setSize( window.innerWidth, window.innerHeight );
  renderer.setClearColor( 0xff0000, 0.0 );
  //renderer.setClearColor( 0xffffffff, 0.0 );
  renderer.setAnimationLoop( animate );
  document.body.appendChild( renderer.domElement );

  scene = new THREE.Scene();

  camera_persp = new THREE.PerspectiveCamera( 40, window.innerWidth / window.innerHeight, 0.1, 10 );
  camera_persp.position.set( - 1, 0, 2 );

  const aspect = window.innerWidth / window.innerHeight;
  camera_ortho = new THREE.OrthographicCamera( frustumSize * aspect / - 2, frustumSize * aspect / 2, frustumSize / 2, frustumSize / - 2, 0.1, 10 );
  camera_ortho.position.set( - 1, 0, 2 );

  controls_ortho = new OrbitControls( camera_ortho, renderer.domElement );
  controls_ortho.enableDamping = true;
  controls_ortho.minDistance = 0.1;
  controls_ortho.maxDistance = 500;

  controls_persp = new OrbitControls( camera_persp, renderer.domElement );
  controls_persp.enableDamping = true;
  controls_persp.minDistance = 0.1;
  controls_persp.maxDistance = 500;


  setup_viz();
  window.addEventListener( 'resize', onWindowResize );
  onWindowResize();
}


function onWindowResize() {
  camera_ortho.aspect = window.innerWidth / window.innerHeight;
  camera_ortho.updateProjectionMatrix();

  camera_persp.aspect = window.innerWidth / window.innerHeight;
  camera_persp.updateProjectionMatrix();

  renderer.setSize( window.innerWidth, window.innerHeight );

  insetWidth = window.innerHeight / 4;
  insetHeight = window.innerHeight / 4;
}


function animate() {

  renderer.setClearColor( 0x222222, 1 );
  renderer.setViewport( 0, 0, window.innerWidth, window.innerHeight );
  controls_ortho.update();
  controls_persp.update();

  if (g_data.camera == "ortho") {
    renderer.render( scene, camera_ortho );
  }
  else {
    renderer.render( scene, camera_persp );
  }

  renderer.setClearColor( 0x222222, 1 );
  renderer.clearDepth();
  renderer.setScissorTest( true );
  renderer.setScissor( 20, 20, insetWidth, insetHeight );
  renderer.setViewport( 20, 20, insetWidth, insetHeight );
  renderer.setScissorTest( false );

}


