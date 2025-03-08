//import * as THREE from './three.js';
import * as THREE from 'three';

import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

import { Line2 } from 'three/addons/lines/Line2.js';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineGeometry } from 'three/addons/lines/LineGeometry.js';
import * as GeometryUtils from 'three/addons/utils/GeometryUtils.js';

import { FontLoader } from 'three/addons/loaders/FontLoader.js';

var legend_anchor = [ -0.1, -0.1, -0.1 ],
    legend_dxyz = [ 0.1, 0.1, 0.1 ];

var _numeric = await import("./njs.mjs");
var njs = _numeric.Numeric;

var IDIR_COLOR = [
  [ 1,0,0 ], [ 1,0,0 ],
  [ 0,1,0 ], [ 0,1,0 ],
  [ 0,0,1 ], [ 0,0,1 ]
];

function simple_line(u,v,c, lw) {
  c = ((typeof c === "undefined") ? [Math.random(), Math.random(), Math.random() ] : c );
  lw = ((typeof lw === "undefined") ? 2 : lw);

  let matLine = new LineMaterial( {
    color: 0xffffff,
    linewidth: lw,
    vertexColors: true,
    dashed: false,
    alphaToCoverage: true
  } );


  let points = [];
  let colors = [];
  points.push( u[0], u[1], u[2] );
  points.push( v[0], v[1], v[2] );
  colors.push( c[0], c[1], c[2] );
  if (c.length >= 6) {
    colors.push( c[3], c[4], c[5] );
  }
  else {
    colors.push( c[0], c[1], c[2] );
  }

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

  let fi_info = lunet.frustum3d_intersection(p, q, frustum_box_r);

  // VIZ
  // VIZ
  // VIZ

  let idir_descr = [ "+x", "-x", "+y", "-y", "+z", "-z" ];

  console.log(fi_info);
  console.log("frustum_idir:", fi_info.frustum_idir, "(", idir_descr[fi_info.frustum_idir], ")");

  // frame vector viz
  //

  let fidir = fi_info.frustum_idir;
  let fv = fi_info.frustum_v[fidir];
  let frame_t = fi_info.frame_t;
  let frame_u = fi_info.frame_updated;
  for (let idx=0; idx < frame_t.length; idx++) {

    let vs = njs.add(p, fv[idx]);
    let ve = njs.add(p, fv[(idx+1)%frame_t.length]);

    let dv = njs.sub(ve,vs);
    let du = njs.sub(vs,ve);

    let cv = IDIR_COLOR[lunet.v2idir(dv)];
    let cu = IDIR_COLOR[lunet.v2idir(du)];

    if (frame_u[idx][0] == 1) {
      let t = frame_t[idx][0];
      simple_line( vs, njs.add(vs, njs.mul(t, dv)), [0.0,0,0, cv[0],cv[1],cv[2]], 1.5);
    }

    if (frame_u[idx][1] == 1) {
      let t = 1-frame_t[idx][1];
      simple_line( ve, njs.add(ve, njs.mul(t, du)), [0.0,0,0, cu[0],cu[1],cu[2]], 1.5);


    }

  }

  //
  // frame vector viz



  let ldx = njs.add( legend_anchor, [ legend_dxyz[0], 0, 0 ] ),
      ldy = njs.add( legend_anchor, [ 0, legend_dxyz[1], 0 ] ),
      ldz = njs.add( legend_anchor, [ 0, 0, legend_dxyz[2] ] );

  let xyz_legend = [
    [ legend_anchor, ldx, [0.5, 0, 0] ],
    [ legend_anchor, ldy, [0, 0.5, 0] ],
    [ legend_anchor, ldz, [0, 0, 0.5] ]
  ];

  for (let xyz=0; xyz<3; xyz++) {
    simple_line( xyz_legend[xyz][0], xyz_legend[xyz][1], xyz_legend[xyz][2] );
  }

  // base grid
  //
  for (let z=0; z<=grid_n; z++) {
    for (let y=0; y<=grid_n; y++) {
      for (let x=0; x<=grid_n; x++) {
        let cur = [ x/grid_n, y/grid_n, z/grid_n ];

        let nxt = [ (x+1)/grid_n, y/grid_n, z/grid_n ];
        if ((nxt[0] <= 1) && (nxt[1] <= 1) && (nxt[2] <= 1)) {
          simple_line( cur, nxt, [0.2, 0.2, 0.2], 0.5);
        }

        nxt = [ (x)/grid_n, (y+1)/grid_n, z/grid_n ];
        if ((nxt[0] <= 1) && (nxt[1] <= 1) && (nxt[2] <= 1)) {
          simple_line( cur, nxt, [0.2, 0.2, 0.2], 0.5);
        }

        nxt = [ (x)/grid_n, (y)/grid_n, (z+1)/grid_n ];
        if ((nxt[0] <= 1) && (nxt[1] <= 1) && (nxt[2] <= 1)) {
          simple_line( cur, nxt, [0.2, 0.2, 0.2], 0.5);
        }


      }
    }
  }

  simple_line( p, q, [1,1,1] );
  simple_line( p, njs.add(p, njs.mul(l0, v_idir[p_near_idir])), [1,0,1] );

  let _shift = 1/512;
  _shift = 0;

  for (let idir=0; idir<6; idir++) {

    for (let v_idx=0; v_idx<fi_info.frustum_v[idir].length; v_idx++) {
      let pp = njs.add(p, njs.mul(_shift, v_idir[idir]));
      let tv = njs.add(pp, fi_info.frustum_v[idir][v_idx]);
      simple_line( pp, tv, [0.75, 0.75, 0.75 ], 0.5 );
    }

    let prv = njs.add(p, fi_info.frustum_v[idir][ fi_info.frustum_v[idir].length-1 ]);
    for (let v_idx=0; v_idx<fi_info.frustum_v[idir].length; v_idx++) {
      let tv = njs.add(p, fi_info.frustum_v[idir][v_idx]);
      simple_line(prv, tv, [0.35, 0.2, 0.2], 0.45 );
      prv = tv;
    }

  }


  //----

  // perpendicular q-plane viz
  //

  let qt = njs.sub( q, p );
  let Nqt = njs.mul( 1/njs.norm2(qt), qt );

  let thetax = Math.abs(lunet.v3theta( [1,0,0], Nqt ));
  let thetay = Math.abs(lunet.v3theta( [0,1,0], Nqt ));
  let thetaz = Math.abs(lunet.v3theta( [0,0,1], Nqt ));

  let qt_plane_pnt = [0,0,0];

  if ((thetax <= thetay) && (thetax <= thetaz)) {
    qt_plane_pnt[0] = -Nqt[1] / Nqt[0];
    qt_plane_pnt[1] = 1;
  }
  else if ((thetay <= thetax) && (thetay <= thetaz)) {
    qt_plane_pnt[1] = -Nqt[0] / Nqt[1];
    qt_plane_pnt[0] = 1;
  }
  else {
    qt_plane_pnt[2] = -Nqt[1] / Nqt[2];
    qt_plane_pnt[1] = 1;
  }

  qt_plane_pnt = njs.mul( 1/njs.norm2(qt_plane_pnt), qt_plane_pnt );

  let plane_r = frustum_box_r * Math.sqrt(3);

  let n_seg = 32;
  for (let seg=0; seg<n_seg; seg++) {
    let theta_cur = 2*Math.PI * (seg+0) / n_seg;
    let theta_nxt = 2*Math.PI * (seg+1) / n_seg;
    let qr0 = lunet.rodrigues( qt_plane_pnt, Nqt, theta_cur );
    let qr1 = lunet.rodrigues( qt_plane_pnt, Nqt, theta_nxt );

    let qrt = njs.add( q, njs.mul( plane_r, qr0 ) );

    let _u0 = njs.add(q, njs.mul(plane_r, qr0));
    let _u1 = njs.add(q, njs.mul(plane_r, qr1));

    simple_line(q, qrt, [0,0,0, 0.1, 0.2, 0.3], 0.3);
    simple_line(_u0, _u1, [0.1, 0.2, 0.3], 0.3);
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

    let ldx = njs.add(legend_anchor, [legend_dxyz[0], 0, 0 ]),
        ldy = njs.add(legend_anchor, [0, legend_dxyz[1], 0 ]),
        ldz = njs.add(legend_anchor, [0, 0, legend_dxyz[1] ]);

    g_data["font"] = font;
    //simple_text('z', [-0.65, -0.55, -0.2 ], 0.035 );
    simple_text('z', ldz, 0.035 );
    simple_text('y', ldy, 0.035 );
    simple_text('x', ldx, 0.035 );
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


