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
  "line":[],
  "Line":[],
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
    g_plotter_info.pnt[ii].stroke = '#777';
    g_plotter_info.pnt[ii].fill = '#777';
    if ("c" in pnt[ii]) { g_plotter_info.pnt[ii].fill = pnt[ii].c; }
    if ("c" in pnt[ii]) { g_plotter_info.pnt[ii].stroke= pnt[ii].c; }
    two.add( g_plotter_info.pnt[ii] );
  }

}

function _plotl(line, fx, fy) {
  fx = ((typeof fx === "undefined") ? 1.0 : fx);
  fy = ((typeof fy === "undefined") ? 1.0 : fy);

  let two = g_plotter_info.two;
  let _line = g_plotter_info.line;

  for (let ii=0; ii<_line.length; ii++) {
    two.remove( _line[ii] );
  }

  g_plotter_info.line = [];
  for (let ii=0; ii<line.length; ii++) {
    g_plotter_info.line.push( new Two.Line(line[ii].x, line[ii].y, line[ii].x+line[ii].dx, line[ii].y+line[ii].dy) );
    two.add( g_plotter_info.line[ii] );
  }

}

function _plotL(line, fx, fy) {
  fx = ((typeof fx === "undefined") ? 1.0 : fx);
  fy = ((typeof fy === "undefined") ? 1.0 : fy);

  let two = g_plotter_info.two;
  let _line = g_plotter_info.Line;

  for (let ii=0; ii<_line.length; ii++) {
    two.remove( _line[ii] );
  }

  g_plotter_info.line = [];
  for (let ii=0; ii<line.length; ii++) {
    g_plotter_info.Line.push( new Two.Line(line[ii].x0, line[ii].y0, line[ii].x1, line[ii].y1) );
    two.add( g_plotter_info.Line[ii] );
  }

}

//---
//---
//---


