
let njs = numeric;

var g_ctx = {

  "t_prv": 0,
  "q_update": false,

  "focus": "",
  "shiftkey": false,

  "g_x": {
    "size": [300,300],
    "two" : {},
    "state" : {},
    "data_ref": 50,
    "two_data": [],
    "cursor": {},
    "transform": [],
    "two_data_aux": []
  },

  "g_y": {
    "size": [300,300],
    "two" : {},
    "state" : {},
    "data_ref": 200,
    "two_data": [],
    "cursor": {},
    "transform": [],
    "two_data_aux": []
  },

  "f": {
    "size": [300,300],
    "two" : {},
    "state" : {},
    "data": [],
    "data_ref": 0,
    "two_data": [],

    "two_data_growth": [],

    "cursor": {},
    "transform": [],
    "two_data_aux": []
  }

};


function viz2world(two_data, transform) {
  let viz_pnt = [];
  for (let i=0; i<two_data.length; i++) {
    viz_pnt.push( [ two_data[i].position.x, two_data[i].position.y, 1 ] );
  }
  return njs.transpose( njs.dot( transform, njs.transpose( viz_pnt ) ) );
}

function _print_point( pnt ) {
  for (let i=0; i<pnt.length; i++) {
    console.log( pnt[i][0], pnt[i][1], pnt[i][2] );
  }
}

function afx(x, a) {
  let n = a.length;

  if (x <= a[0][0]) { return a[0][1]; }
  if (x >= a[n-1][0]) { return a[n-1][1]; }

  let iter = 0;
  let max_iter = 30;

  let l_x = a[0][0];
  let l_y = a[0][1];

  let r_x = a[n-1][0];
  let r_y = a[n-1][0];

  let l_idx = 0;
  let r_idx = n-1;
  let m_idx = 0;

  // I'm almost positive ther'es a bug here.
  // I can never remember the boundary conditions.
  // Investigate further
  //
  while ( (r_idx-l_idx) > 1 ) {
    m_idx = l_idx + Math.floor( (r_idx-l_idx)/2 );
    if ( x <= a[m_idx][0] ) {
      r_idx = m_idx;
    }
    else {
      l_idx = m_idx;
    }

    //console.log(l_idx, m_idx, r_idx, "(x:", x, ", a[", l_idx, "]:", a[l_idx], "a[", m_idx, "]", a[m_idx], "a[", r_idx, "]", a[r_idx], ")");

    iter++;
    if (iter >= max_iter) {
      console.log("!!!!", l_x, m_x, r_x, x );
      break;
    }

  }

  if (a[m_idx][0] < x) {

    //console.log("cp0");

    l_x = a[m_idx][0];
    l_y = a[m_idx][1];

    if ( (m_idx+1) < n ) {
      r_x = a[m_idx+1][0];
      r_y = a[m_idx+1][1];
    }
  }
  else {
    //console.log("cp1");

    if ( (m_idx-1) >= 0 ) {
      l_x = a[m_idx-1][0];
      l_y = a[m_idx-1][1];
    }

    r_x = a[m_idx][0];
    r_y = a[m_idx][1];
  }


  let p = (x - l_x) / (r_x - l_x);

  //console.log( "p:", p, "r_x - l_x:", r_x - l_x, "x - l_x", x - l_x, "r_y-l_y:", r_y-l_y, "*p:", (r_y-l_y)*p);

  return l_y + ((r_y - l_y)*p);

}


function _normalize_g(a) {
  let S = 0.0;
  for (let i=1; i<a.length; i++) {
    S += a[i]
  }
}

function grow_f() {

  let f_pnt = viz2world( g_ctx.f.two_data, g_ctx.f.transform );
  let g_x_pnt = viz2world( g_ctx.g_x.two_data, g_ctx.g_x.transform );
  let g_y_pnt = viz2world( g_ctx.g_y.two_data, g_ctx.g_y.transform );

  let W = g_ctx.f.size[0];
  let H = g_ctx.f.size[1];
  let w2 = W/2;

  g_x_pnt.sort( function(a,b) { return (a[0] < b[0]) ? -1 : ( (a[0] == b[0]) ? 0 : 1 ) ; } )
  g_y_pnt.sort( function(a,b) { return (a[0] < b[0]) ? -1 : ( (a[0] == b[0]) ? 0 : 1 ) ; } )

  _normalize_g(g_x_pnt);
  _normalize_g(g_y_pnt);


  let ds = 1/6;

  if ( g_ctx.f.two_data_growth.length > 0 ) {
    g_ctx.f.two.remove( g_ctx.f.two_data_growth );
    g_ctx.f.two_data_growth = [];
    g_ctx.f.two.update();
  }

  let boundary_point = [];
  let cur_boundary = f_pnt;

  boundary_point.push(cur_boundary);

  for (let it=0; it<10; it++) {

    let nxt_boundary = [];

    let _gpnt = [];
    for (let i=0; i<cur_boundary.length; i++) {
      let x = cur_boundary[i][0];
      let y = cur_boundary[i][1];

      // experimenting...
      let dx = ds*afx( y, g_y_pnt );
      let dy = ds*afx( y, g_x_pnt );

      let x1 = 0;
      let y1 = 0;

      let blech = false;
      if (blech) {
        if (x < 0) { dx *= -1; }
        let x1 = x + dx;
        let y1 = y + dy;
      }
      else {
        if (x < 0) {
          x1 = -x*(1+dx/20);
        }
        else {
          x1 = x*(1+dx/20);
        }
        y1 = y*(1+dy/20);
      }

      nxt_boundary.push( [x1, y1, 1 ] );

      let vzx = x1+w2;
      let vzy = H-y1;

      let _c = g_ctx.f.two.makeCircle( vzx, vzy, 0.5 );
      _c.fill = "rgba(255,0,0,0.2)";
      _c.stroke = "rgba(255,0,0,0.2)";

      g_ctx.f.two_data_growth.push(_c);
    }

    boundary_point.push( nxt_boundary );
    cur_boundary = nxt_boundary;

  }

  g_ctx.f.two.update();
}

function _mouseenter(id, ev) {
  g_ctx[id].state = "me";
  g_ctx.focus = id;
}

function _mousedown(id, ev) {
  g_ctx[id].state = "md";
  g_ctx.focus = id;

  if (g_ctx.focus == "g_x") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.9)';
    g_ctx[id].two.update();
  }

  if (g_ctx.focus == "g_y") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.9)';
    g_ctx[id].two.update();
  }


}

function _mouseout(id, ev) {
  g_ctx[id].state = "mo";
  g_ctx.focus = '';

  if (g_ctx.focus == "g_x") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';
    g_ctx[id].two.update();
  }

  if (g_ctx.focus == "g_y") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';
    g_ctx[id].two.update();
  }


}

function _mouseup(id, ev) {
  g_ctx[id].state = "mu";

  if (g_ctx.focus == "g_x") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';
    g_ctx[id].two.update();
  }

  if (g_ctx.focus == "g_y") {
    g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';
    g_ctx[id].two.update();
  }

}

function _mousemove(id, ev) {

  if (g_ctx.focus == "g_x") {

    if ((g_ctx[id].state == 'mu') ||
        (g_ctx[id].state == 'me')) {
      let _x = ev.offsetX;
      let _y = ev.offsetY;

      g_ctx[id].cursor.position.x = _x;
      g_ctx[id].cursor.position.y = _y;

      g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';

      let t = Date.now();
      if ( (t - g_ctx.t_prv) > 10 ) {
        g_ctx.t_prv = t;
        g_ctx[id].two.update();
      }

    }
    else if (g_ctx[id].state == 'md') {

      let _x = ev.offsetX;
      let _y = ev.offsetY;

      g_ctx[id].cursor.position.x = _x;
      g_ctx[id].cursor.position.y = _y;

      g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.9)';

      let rx = g_ctx[id].cursor.width/2;
      let ry = g_ctx[id].cursor.height/2;

      let pnt = g_ctx[id].two_data;
      for (let ii=0; ii<pnt.length; ii++) {
        let px = pnt[ii].position.x;
        let py = pnt[ii].position.y;

        let dx = px-_x;
        let dy = py-_y;

        let _d2 = (dx*dx/(rx*rx)) + (dy*dy/(ry*ry));

        if ( _d2 < 1 ) {

          let _disc = 1 - ((dy*dy)/(ry*ry));
          if (_disc < 0) { continue; }

          let dx0 = rx * Math.sqrt(_disc);
          let dx1 = -dx0;

          let tx0 = [ _x + dx0, py ];
          let tx1 = [ _x + dx1, py ];

          let jx0 = Math.abs( px - tx0[0] );
          let jx1 = Math.abs( px - tx1[0] );

          if (ev.shiftKey) {
            if (tx1[0] < g_ctx[id].data_ref) {
              pnt[ii].position.x = g_ctx[id].data_ref;
            }
            else {
              pnt[ii].position.x = tx1[0];
            }
          }

          else if ((jx1 < jx0) &&
                   (tx1[0] > g_ctx[id].data_ref)) {
            pnt[ii].position.x = tx1[0];
          }
          else {
            pnt[ii].position.x = tx0[0];
          }

        }

      }

      let t = Date.now();
      if ( (t - g_ctx.t_prv) > 10 ) {
        g_ctx.t_prv = t;
        g_ctx[id].two.update();
      }

    }

  }

  if (g_ctx.focus == "g_y") {

    let _x = ev.offsetX;
    let _y = ev.offsetY;

    g_ctx[id].cursor.position.x = _x;
    g_ctx[id].cursor.position.y = _y;

    if ((g_ctx[id].state == 'mu') ||
        (g_ctx[id].state == 'me')) {
      g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.1)';

    }
    else if (g_ctx[id].state == 'md') {
      g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.9)';

      let rx = g_ctx[id].cursor.width/2;
      let ry = g_ctx[id].cursor.height/2;

      let pnt = g_ctx[id].two_data;
      for (let ii=0; ii<pnt.length; ii++) {
        let px = pnt[ii].position.x;
        let py = pnt[ii].position.y;

        let dx = px-_x;
        let dy = py-_y;

        let _d2 = (dx*dx/(rx*rx)) + (dy*dy/(ry*ry));

        if ( _d2 < 1 ) {

          let _disc = 1 - ((dx*dx)/(rx*rx));
          if (_disc < 0) { continue; }

          let dy0 = ry * Math.sqrt(_disc);
          let dy1 = -dy0;

          let t0 = [ _x, _y + dy0 ];
          let t1 = [ _x, _y + dy1 ];

          let jy0 = Math.abs( py - t0[1] );
          let jy1 = Math.abs( py - t1[1] );

          if (ev.shiftKey) {
            if (t0[1] > g_ctx[id].data_ref) {
              pnt[ii].position.y = g_ctx[id].data_ref;
            }
            else {
              pnt[ii].position.y = t0[1];
            }
          }

          else if ((jy0 < jy1) &&
                   (t0[1] < g_ctx[id].data_ref)) {
            pnt[ii].position.y = t0[1];
          }
          else {
            pnt[ii].position.y = t1[1];
          }

        }

      }

    }

    let t = Date.now();
    if ( (t - g_ctx.t_prv) > 10 ) {
      g_ctx.t_prv = t;
      g_ctx[id].two.update();
    }

    /*
    else if (g_ctx[id].state == 'md') {

      let _x = ev.offsetX;
      let _y = ev.offsetY;

      g_ctx[id].cursor.position.x = _x;
      g_ctx[id].cursor.position.y = _y;

      g_ctx[id].cursor.stroke = 'rgba(0,0,0,0.9)';

      let rx = g_ctx[id].cursor.width/2;
      let ry = g_ctx[id].cursor.height/2;

      let r = rx;

      let t = Date.now();
      if ( (t - g_ctx.t_prv) > 10 ) {
        g_ctx.t_prv = t;
        g_ctx[id].two.update();
      }

    }
    */

  }



}

function _mousewheel(ev) {

  if (g_ctx.focus == "g_x") {


    if (ev.shiftKey) {
      let rx = g_ctx[ g_ctx.focus ].cursor.width/2;
      let ry = g_ctx[ g_ctx.focus ].cursor.height/2;
      if (ev.wheelDelta < 0) {
        rx--;
        if (ev.buttons == 0) { ry--; }
        if ((rx > 5) &&
            (ry > 5)) {
          g_ctx[ g_ctx.focus ].cursor.width = 2*rx;
          g_ctx[ g_ctx.focus ].cursor.height = 2*ry;
          g_ctx[ g_ctx.focus ].two.update();
        }
      }
      else if (ev.wheelDelta > 0) {
        rx++;
        if (ev.buttons == 0) { ry++; }
        if ((rx < 100) &&
            (ry < 100)) {
          g_ctx[ g_ctx.focus ].cursor.width = 2*rx;
          g_ctx[ g_ctx.focus ].cursor.height = 2*ry;
          g_ctx[ g_ctx.focus ].two.update();
        }

      }
    }

  }

  if (g_ctx.focus == "g_y") {
    let ctx = g_ctx[ g_ctx.focus ];

    if (ev.shiftKey) {
      let rx = ctx.cursor.width/2;
      let ry = ctx.cursor.height/2;
      if (ev.wheelDelta < 0) {
        ry--;
        if (ev.buttons == 0) { rx--; }
        if ((rx > 5) &&
            (ry > 5)) {
          ctx.cursor.width = 2*rx;
          ctx.cursor.height = 2*ry;
          ctx.two.update();
        }
      }
      else if (ev.wheelDelta > 0) {
        ry++;
        if (ev.buttons == 0) { rx++; }
        if ((rx < 100) &&
            (ry < 100)) {
          ctx.cursor.width = 2*rx;
          ctx.cursor.height = 2*ry;
          ctx.two.update();
        }

      }
    }

    ctx.two.update();
  }
}

function setup_g_x_graph() {
  let N = 100;

  let ctx = g_ctx.g_x;

  let H = ctx.size[1];
  let W = ctx.size[0];

  let _lx = ctx.data_ref;
  let _lh = ctx.size[1];
  let _ref_line = ctx.two.makeLine( _lx, 0, _lx, _lh );
  _ref_line.stroke = "rgba(0,0,0,0.3)";
  ctx.two_data_aux.push( _ref_line );

  let stride = 30;
  for (let y = (ctx.size[1]-stride); y >= 0; y -= stride) {
    let _guide_line = ctx.two.makeLine( 0, y, ctx.size[0], y );
    _guide_line.stroke = "rgba(32,32,32,0.15)";
    ctx.two_data_aux.push( _guide_line );
  }

  for (let i=0; i<N; i++) {
    let _y = i * ctx.size[1] / N;
    let _x = ctx.data_ref;
    let _c = ctx.two.makeCircle( _x, _y, 0.5 );
    ctx.two_data.push(_c);
  }
  ctx.two.update();

  ctx.transform = [
    [ 0,-1,  H ],
    [ 1, 0, -ctx.data_ref  ],
    [ 0, 0, 1 ]
  ];

  //console.log( njs.dot( ctx.transform, [ctx.data_ref, 20, 1] ) );
  //console.log( njs.dot( ctx.transform, [ctx.data_ref+5, 25, 1] ) );
  //console.log( njs.dot( ctx.transform, [ctx.data_ref, H, 1] ) );
}

function setup_g_y_graph() {

  let N = 100;

  let ctx = g_ctx.g_y;

  let W = ctx.size[0];
  let H = ctx.size[1];

  let w2 = W/2;

  let _ly = ctx.data_ref;
  let _lw = ctx.size[0];
  let _ref_line = ctx.two.makeLine( 0, _ly, _lw, _ly );
  _ref_line.stroke = "rgba(0,0,0,0.3)";
  ctx.two_data_aux.push( _ref_line );

  let stride = 30;
  for (let x = 0; x < ctx.size[0]; x += stride) {
    let _guide_line = ctx.two.makeLine( x, 0, x, ctx.size[1] );
    _guide_line.stroke = "rgba(32,32,32,0.15)";
    ctx.two_data_aux.push( _guide_line );
  }

  for (let i=0; i<N; i++) {
    let _x = i * ctx.size[0] / N;
    let _y = ctx.data_ref;
    let _c = ctx.two.makeCircle( _x, _y, 0.5 );
    ctx.two_data.push(_c);
  }
  ctx.two.update();


  ctx.transform = [
    //[ 1, 0,  -w2 ],
    [ 1, 0,  0 ],
    [ 0, -1, +ctx.data_ref],
    [ 0, 0, 1 ]
  ];

  let _debug = false;
  if (_debug) {
    console.log( njs.dot( ctx.transform, [0, ctx.data_ref, 1] ) );
    console.log( njs.dot( ctx.transform, [w2, ctx.data_ref, 1] ) );
    console.log( njs.dot( ctx.transform, [w2, ctx.data_ref - 30, 1] ) );
    console.log( njs.dot( ctx.transform, [W, ctx.data_ref, 1] ) );
  }
}

function setup_f_graph() {

  let ctx = g_ctx.f;

  let H = ctx.size[1];
  let W = ctx.size[0];

  let w2 = W/2;

  let seedling = [];

  let h0 = 50;
  let w0 = 25;
  let r = w0;

  let n_l = 400;
  let n_r = 400;
  let n_t = 501;

  for (let i=0; i<n_l; i++) {
    let y = H - (h0*i/n_l);
    let x = w2 - w0;
    let _c = ctx.two.makeCircle( x,y, 0.5 );
    ctx.two_data.push(_c);
  }

  let cxy = {"x": w2, "y": H-h0 };
  for (let i=0; i<n_t; i++) {
    let theta = Math.PI - (i*Math.PI/(n_t-1));
    let x = Math.cos(-theta)*r + cxy.x;
    let y = Math.sin(-theta)*r + cxy.y;
    let _c = ctx.two.makeCircle( x,y, 0.5 );
    ctx.two_data.push(_c);
  }

  for (let i=(n_r-1); i>=0; i--) {
    let y = H - (h0*i/n_r);
    let x = w2 + w0;
    let _c = ctx.two.makeCircle( x,y, 0.5 );
    ctx.two_data.push(_c);
  }

  ctx.transform = [
    [ 1, 0, -w2 ],
    [ 0,-1,  H ],
    [ 0, 0, 1]
  ];

  // for context
  //
  let _debug = false;
  if (_debug) {
    console.log( njs.dot( ctx.transform, [w2-w0, H, 1] ));
    console.log( njs.dot( ctx.transform, [w2-w0, H-h0, 1] ));
    console.log( njs.dot( ctx.transform, [w2, H-h0-w0, 1] ));
    console.log( njs.dot( ctx.transform, [w2+w0, H-h0, 1] ));
    console.log( njs.dot( ctx.transform, [w2+w0, H, 1] ));
  }

  ctx.two.update();
}

function init() {

  document.addEventListener("wheel", function(ev) {
    _mousewheel(ev);
  });

  let div_container_id = [ "ui_g_y", "ui_g_x", "ui_f" ];
  let key = [ "g_y", "g_x", "f" ];

  let _debug = [ "#c22", "#2c2", "#22c" ];

  for (let ii=0; ii<div_container_id.length; ii++) {

    let w = g_ctx[key[ii]].size[0];
    let h = g_ctx[key[ii]].size[1];

    let two = new Two({
      "width": w,
      "height": h,
      "type": Two.Types.canvas
    });
    g_ctx[ key[ii] ].two = two;

    //two.renderer.domElement.style.background = 'rgb(255,255,255)';

    let ele = document.getElementById( div_container_id[ii] );

    ele.style.userSelect = 'none';

    two.appendTo(ele);
    let canvas = ele.getElementsByTagName('canvas')[0];

    canvas.style.border = "1px solid rgba(0,0,0,0.5)";

    canvas.addEventListener('mousedown', function(ev) {
      _mousedown(key[ii], ev);
    });

    canvas.addEventListener('mousemove', function(ev) {
      _mousemove(key[ii], ev);
    });

    canvas.addEventListener('mouseout', function(ev) {
      _mouseout(key[ii], ev);
    });

    canvas.addEventListener('mouseup', function(ev) {
      _mouseup(key[ii], ev);
    });

    canvas.addEventListener('mouseenter', function(ev) {
      _mouseenter(key[ii], ev);
    });

    //g_ctx[key[ii]].cursor = g_ctx[key[ii]].two.makeCircle(50,50,20);
    g_ctx[key[ii]].cursor = g_ctx[key[ii]].two.makeEllipse(50,50,20, 20);

    if (key[ii] == "f") {
      g_ctx[key[ii]].cursor.stroke = "rgba(0,0,0,0.0)";
    }
    else {
      g_ctx[key[ii]].cursor.stroke = "rgba(0,0,0,0.1)";
    }

    two.update();
  }

  setup_g_x_graph();
  setup_g_y_graph();
  setup_f_graph();

}
