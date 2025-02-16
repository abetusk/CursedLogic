

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
    "two_data_aux": []
  },

  "g_y": {
    "size": [300,300],
    "two" : {},
    "state" : {},
    "data_ref": 200,
    "two_data": [],
    "cursor": {},
    "two_data_aux": []
  },

  "f": {
    "size": [300,300],
    "two" : {},
    "state" : {},
    "data": [],
    "data_ref": 0,
    "two_data": [],
    "cursor": {},
    "two_data_aux": []
  }

};

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

          if ((jx1 < jx0) &&
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

          if ((jy0 < jy1) &&
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

  let _lx = g_ctx.g_x.data_ref;
  let _lh = g_ctx.g_x.size[1];
  let _ref_line = g_ctx.g_x.two.makeLine( _lx, 0, _lx, _lh );
  g_ctx.g_x.two_data_aux.push( _ref_line );

  for (let i=0; i<N; i++) {
    let _y = i * g_ctx.g_x.size[1] / N;
    let _x = g_ctx.g_x.data_ref;
    let _c = g_ctx.g_x.two.makeCircle( _x, _y, 1 );
    g_ctx.g_x.two_data.push(_c);
  }
  g_ctx.g_x.two.update();
}

function setup_g_y_graph() {

  let N = 100;

  let _ly = g_ctx.g_y.data_ref;
  let _lw = g_ctx.g_y.size[0];
  let _ref_line = g_ctx.g_y.two.makeLine( 0, _ly, _lw, _ly );
  g_ctx.g_y.two_data_aux.push( _ref_line );

  for (let i=0; i<N; i++) {
    let _x = i * g_ctx.g_y.size[0] / N;
    let _y = g_ctx.g_y.data_ref;
    let _c = g_ctx.g_y.two.makeCircle( _x, _y, 1 );
    g_ctx.g_y.two_data.push(_c);
  }
  g_ctx.g_y.two.update();
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
    g_ctx[key[ii]].cursor.stroke = "rgba(0,0,0,0.1)";

    two.update();
  }

  setup_g_x_graph();
  setup_g_y_graph();

}
