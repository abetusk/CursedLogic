// LICENSE: CC0
//
//



var njs = require("./numeric.js");
var lunet = require("./lunenetwork.js");
var getopt = require("posix-getopt");


function main(argv) {

  let _eps = 1 / (1024*1024*1024);

  let p = [0,0,0],
      q = [0,0,0],
      ir = 0,
      ds = 0;

  let p_str = "",
      q_str = "",
      ir_str = "",
      ds_str = "";

  let opt = {};
  let parser = new getopt.BasicParser("p:(p)q:(q)d:(ds)r:(ir)", argv);
  while ((opt = parser.getopt()) !== undefined) {
    switch(opt.option) {
      case 'h':
        console.log("help");
        break;
      case 'p':
        p_str = opt.optarg; 
        break;
      case 'q':
        q_str = opt.optarg; 
        break;
      case 'r':
        ir_str = opt.optarg;
        break;
      case 'd':
        ds_str = opt.optarg;
        break;
      default:
        console.log("unknown option", opt);
        break;
    }
  }

  let tok = p_str.split(",");
  if (tok.length == 3) {
    p[0] = parseFloat(tok[0]);
    p[1] = parseFloat(tok[1]);
    p[2] = parseFloat(tok[2]);
  }

  tok = q_str.split(",");
  if (tok.length == 3) {
    q[0] = parseFloat(tok[0]);
    q[1] = parseFloat(tok[1]);
    q[2] = parseFloat(tok[2]);
  }

  if (ds_str.length > 0) {
    ds = parseFloat(ds_str);
  }

  if (ir_str.length > 0) {
    ir = parseInt(ir_str);
  }

  if ((Math.abs(ds) < _eps) ||
      (njs.norm2( njs.sub(p,q) ) < _eps)) {
    console.log("error, ds < _eps (", ds, "<", _eps, ") or p ~= q (", p, "~=", q, ")");
    return;
  }

  // main functionality...
  //

  let v_idir = [
    [1,0,0], [-1, 0, 0],
    [0,1,0], [ 0,-1, 0],
    [0,0,1], [ 0, 0,-1]
  ];


  let grid_n = Math.floor( (1/ds) + _eps );

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


  let fi_info = lunet.frustum3d_intersection(p, q, ds);

  /*
  let p_f_v = [];
  for (let idx=0; idx<fi_info.frustum_v.length; idx++) {
    p_f_v.push([]);
    for (let ii=0; ii<fi_info.frustum_v[idx].length; ii++) {
      p_f_v[idx].push( njs.mul(ds, fi_info.frustum_v[idx][ii]) );
    }
  }

  for (let idx=0; idx<p_f_v.length; idx++) {
    for (let ii=0; ii<p_f_v[idx].length; ii++) {
      //let _vt = njs.add(p, p_f_v[idx][ii]);
      let _vt = njs.add(p, njs.mul( t0, p_f_v[idx][ii]) );
      console.log( _vt[0], _vt[1], _vt[2] );
    }
    //let _vt = njs.add(p, p_f_v[idx][0]);
    let _vt = njs.add(p, njs.mul( t0, p_f_v[idx][0]) );
    console.log( _vt[0], _vt[1], _vt[2] );
    console.log("\n");
  }

  for (let idx=0; idx<p_f_v.length; idx++) {
    for (let ii=0; ii<p_f_v[idx].length; ii++) {
      console.log(p[0], p[1], p[2]);
      let _vt = njs.add(p, p_f_v[idx][ii]);
      console.log( _vt[0], _vt[1], _vt[2] );
      console.log("\n");

      console.log(p[0], p[1], p[2]);
      _vt = njs.add(p, njs.mul(grid_n, p_f_v[idx][ii]));
      console.log( _vt[0], _vt[1], _vt[2] );
      console.log("\n");
    }
  }

  */


  console.log("#grid:");
  for (let x=0; x<=grid_n; x++) {
    for (let y=0; y<=grid_n; y++) {
      for (let z=0; z<=grid_n; z++) {

        for (let pidir=0; pidir<6; pidir+=2) {
          let dxyz = njs.add( [x,y,z], v_idir[pidir] );

          if ((dxyz[0] > grid_n) ||
              (dxyz[1] > grid_n) ||
              (dxyz[2] > grid_n)) { continue; }

          console.log( x*ds, y*ds, z*ds );
          console.log( ds*dxyz[0], ds*dxyz[1], ds*dxyz[2] );
          console.log("\n");
        }
      }
    }
  }

  console.log("#p:");
  console.log(p[0], p[1], p[2], "\n");

  console.log("#q:");
  console.log(p[0], p[1], p[2]);
  console.log(q[0], q[1], q[2], "\n");


  console.log("#Wp", Wp, "ip", ip, "l0:", l0, "t0:", t0, "p_near_idir:", p_near_idir);

  console.log("# p to nearest grid edge straight vector (v_idir dir)");

  console.log(p[0], p[1], p[2]);
  console.log(p[0] + l0*v_idir[p_near_idir][0], p[1] + l0*v_idir[p_near_idir][1], p[2] + l0*v_idir[p_near_idir][2]);
  console.log("\n");


  console.log("#p-frustum:");

  for (let _tr=0; _tr<2; _tr++) {
    console.log("###tr:", _tr);
    for (let idir=0; idir<6; idir++) {
      console.log("##idir:", idir, "(tr:", _tr, ")");
      for (let i=0; i<fi_info.frustum_v[idir].length; i++) {

        //let v = njs.add( p, njs.mul( (t0*Math.sqrt(3)) + ds*_tr*Math.sqrt(3), fi_info.frustum_v[idir][i] ) );
        let v = [ 0,0,0];
        if (_tr==0) {
          v = njs.add( p, njs.mul( t0*Math.sqrt(3), fi_info.frustum_v[idir][i] ) );
        }

        //WIP!!!
        else {
          let f0  = t0*Math.sqrt(3);
          v = njs.add( p, njs.mul( f0 + Math.sqrt(3)*ds*_tr, fi_info.frustum_v[idir][i] ) );
        }

        console.log(p[0], p[1], p[2]);
        console.log(v[0], v[1], v[2]);
        console.log("\n");
      }
    }
  }

  /*
  console.log("# l0 test");
  for (let ii=0; ii<4; ii++) {
    console.log(p[0], p[1], p[2]);
    let _vt = njs.add( p, njs.mul( t0, fi_info.frustum_v[p_near_idir][ii] ) );
    console.log(_vt[0], _vt[1], _vt[2]);
    console.log("\n");
  }


  for (let idir=0; idir < 6; idir++) {
    let fv = fi_info.frustum_v[idir];
    for (let ii=0; ii<fv.length; ii++) {
      console.log(fv[ii][0], fv[ii][1], fv[ii][2]);
    }
    console.log(fv[0][0], fv[0][1], fv[0][2], "\n");
  }
  */


}



main(process.argv);


