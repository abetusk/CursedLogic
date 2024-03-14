var qt = require("./js-quadtree.js");
var halton = require("./halton.js");

let N = 300;
let g_info = {
  "n": 300,
  "bb": {"x":0, "y":0, "w":1000, "h": 10000 },
  "dx": 10000, "dy": 10000,
  "r_kill": -1,
  "r_inf": -1,
  "eps": -1,
  "p" : []
};

let rho = g_info.dx / Math.sqrt(g_info.n);

g_info.r_kill = 0.5*rho;
g_info.r_inf = 1*rho;
g_info.eps = 1/(1024*1024);

let halton_cx = 0.5,
    halton_cy = 0.5;
let tree = new qt.QuadTree( new qt.Box( g_info.bb.x, 
                                        g_info.bb.y,
                                        g_info.bb.w,
                                        g_info.bb.h ) );

let seq = halton.seq2d(2,3,g_info.n);

for (let ii=0; ii<seq.length; ii++) {
  g_info.p.push( new qt.Point(Math.floor(g_info.bb.w*seq[ii][0]), Math.floor(g_info.bb.h*seq[ii][1])) );
  //tree.insert(new qt.Point(seq[ii][0]-cx, seq[ii][1]-cy));
}

tree.insert(g_info.p);

//for (let ii=0; ii<g_info.p.length; ii++) {
  //let x = g_info.p[ii].x;
  //let y = g_info.p[ii].y;

  //tree.insert({"x":x, "y":y, "width": r, "height": r});
  //tree.insert(g_info.p[ii]);
  //console.log(g_info.p[ii].x, g_info.p[ii].y);

  //console.log(g_info.p[ii]);
//}

//process.exit();

//let xx = tree.retrieve({"x":-2, "y":-2, "width":10, "height":10});

let cpnt = new qt.Circle(0, 0, g_info.r_inf);
cpnt.x = 0.0;
cpnt.y = 0.0;
let xx = tree.query(cpnt);

//console.log(tree.query( new qt.Box(-10,-10, 10,10) ) );

for (let ii=0; ii<xx.length; ii++) {
  console.log(xx[ii].x, xx[ii].y);
}

let idx = Math.floor( N * Math.random() );
let wf = [ seq[idx][0], seq[idx][1] ];

//console.log(">>", tree.objects.length, tree.nodes.length);


//console.log(xx);
//console.log("##");

//for (let ii=0; ii<xx.length; ii++) { console.log(xx[ii].x, xx[ii].y); }



