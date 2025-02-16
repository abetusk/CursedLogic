

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
    "data_ref": 50,
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
}

function _mouseout(id, ev) {
  g_ctx[id].state = "mo";
  g_ctx.focus = '';
}

function _mouseup(id, ev) {
  g_ctx[id].state = "mu";
}

function _mousemove(id, ev) {

  if (g_ctx[id].state == 'md') {

    let _x = ev.offsetX;
    let _y = ev.offsetY;

    g_ctx[id].cursor.position.x = _x;
    g_ctx[id].cursor.position.y = _y;

    let r = g_ctx[id].cursor.radius;

    //console.log(_x, _y);

    let pnt = g_ctx[id].two_data;
    for (let ii=0; ii<pnt.length; ii++) {
      let px = pnt[ii].position.x;
      let py = pnt[ii].position.y;

      let dx = px-_x;
      let dy = py-_y;

      if ( Math.sqrt(dx*dx + dy*dy) < r ) {

        let dx0 = Math.sqrt( (r*r) - (dy*dy) );
        let dx1 = -dx0;

        let theta = Math.atan2(-(py-_y), px-_x);

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
  //g_ctx[id].state = "mu":
}

function _mousewheel(ev) {
  console.log("?", ev.shiftKey, g_ctx.focus, (ev.wheelDelta < 0) ? -1 : 1);

  if (g_ctx.focus == "g_x") {
    if (ev.shiftKey) {
      let r = g_ctx[ g_ctx.focus ].cursor.radius;
      if (ev.wheelDelta < 0) {
        r--;
        if (r > 5) {
          g_ctx[ g_ctx.focus ].cursor.radius = r;
          g_ctx[ g_ctx.focus ].two.update();
        }
      }
      else if (ev.wheelDelta > 0) {
        r++;
        if (r < 100) {
          g_ctx[ g_ctx.focus ].cursor.radius = r;
          g_ctx[ g_ctx.focus ].two.update();
        }

      }
    }
  }
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

    two.renderer.domElement.style.background = 'rgb(255,255,255)';

    let ele = document.getElementById( div_container_id[ii] );

    ele.style.userSelect = 'none';

    two.appendTo(ele);
    let canvas = ele.getElementsByTagName('canvas')[0];

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

    g_ctx[key[ii]].cursor = g_ctx[key[ii]].two.makeCircle(50,50,20);

    two.update();
  }

  // setup g_x (growth in x direction
  //

  let N = 100;

  let _lx = g_ctx.g_x.data_ref;
  let _lh = g_ctx.g_x.size[1];
  let _ref_line = g_ctx.g_x.two.makeLine( _lx, 0, _lx, _lh );
  g_ctx.g_x.two_data_aux.push( _ref_line );



  for (let i=0; i<N; i++) {
    let _y = i * g_ctx.g_x.size[1] / N;
    let _x = g_ctx.g_x.data_ref;

    //g_ctx.g_x.data.push({ "x": _x, "y": _y });

    let _c = g_ctx.g_x.two.makeCircle( _x, _y, 1 );
    g_ctx.g_x.two_data.push(_c);
  }
  g_ctx.g_x.two.update();



}
