// LICENSE: CC0
//
// To the extent possible under law, the person who associated CC0 with
// this file has waived all copyright and related or neighboring rights
// to this file.
//
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//

var g_plotter_info = {
  "pnt":[],
  "two": null,
  "data": []
};

//---
//---
//---

function update_two() {


}

function init_two() {
  g_plotter_info.two = new Two({
    "type": Two.Types.svg,
    "fullscreen": true,
    "autostart": true
  }).appendTo(document.body);

  g_plotter_info.two.bind('update', update_two);

  g_plotter_info.pnt.push( new Two.Points( new Two.Vector(1,2) ) );
  g_plotter_info.pnt.push( new Two.Points( new Two.Vector(100,50) ) );
  g_plotter_info.pnt.push( new Two.Points( new Two.Vector(23,90) ) );
  g_plotter_info.pnt.push( new Two.Points( new Two.Vector(200,50) ) );

  for (let ii=0; ii<g_plotter_info.pnt.length; ii++) {
    g_plotter_info.pnt[ii].size = 3;
    g_plotter_info.pnt[ii].fill = '#777';
    g_plotter_info.two.add( g_plotter_info.pnt[ii] );
  }

  //g_plotter_info.two.makePoints(1,2,5,9);
}

function _plot(pnt, fx, fy) {
  fx = ((typeof fx === "undefined") ? 1.0 : fx);
  fy = ((typeof fy === "undefined") ? 1.0 : fy);

  let two = g_plotter_info.two;
  let _pnt = g_plotter_info.pnt;

  for (let ii=0; ii<_pnt.length; ii++) {
    two.remove( _pnt[ii] );
  }

  g_plotter_info.pnt = [];
  for (let ii=0; ii<pnt.length; ii++) {
    g_plotter_info.pnt.push( new Two.Points( new Two.Vector( fx*pnt[ii].x, fy*pnt[ii].y ) ) );
    g_plotter_info.pnt[ii].size = 3;
    g_plotter_info.pnt[ii].fill = '#777';
    two.add( g_plotter_info.pnt[ii] );
  }

}

//---
//---
//---


function init() {

  console.log("...init...");

  //init_ws();

  init_two();
}


/*
var WEBSOCKET_URL = "ws://localhost";
var WEBSOCKET_PORT = "3001";

//---
//---
//---


// main entry point for new data.
//
function ws_message(ev) {
  let _data = [];

  let z = ev.data.split(",");
  for (let ii=0; ii<z.length; ii++) {
    tok = z[ii].split(" ");
    if (tok.length != 2) { continue; }
    _data.push( { "x": parseFloat(tok[0]), "y": parseFloat(tok[1]) } );
  }

  g_plotter_info.data = _data;

  //console.log("ws_message:", _data);
}

function init_ws() {
  var ws = new WebSocket(WEBSOCKET_URL + ":" + WEBSOCKET_PORT);
  ws.onmessage = ws_message;
}
*/


