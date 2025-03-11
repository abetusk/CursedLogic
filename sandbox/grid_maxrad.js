// To the extent possible under law, the person who associated CC0 with
// this project has waived all copyright and related or neighboring rights
// to this project.
// 
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.


// Helper functions to decorate grid elements with maximum radius information.
// For computing the relative neighborhood graph with an expanding fence and/or
// shrinking window, it could be easy to put the point distribution in a geometry
// that is ill conditioned to the algorithm.
// For example, if there were a thin sliver of points just above the $\pi/4$ line,
// this would eat into the top right fence/window but would have difficulty securing
// or closing the left hand side.
//
// If there are grid points that completely enclose the fence corner, we can use that
// to update the fence and window.
//
// If the fence is large enough, there might be a comparitively small grid cell that
// encloses the fence endpoint.
// We want to further update the window by taking the longest contiguous region of the
// grid that fully encloses the fence side.
// We could do this by walking the grid until we hit the first point occupied grid cell,
// but my worry is that this will become quadratic.
//
// Instead, we can pre-process the grid by decorating the grid cells with the maximum
// square/cube that it can expand to before it hits a grid cell with at least one point
// in it.
// The decoration is pretty much a dynamic programming solution, with two passes, one
// for each diagonal direction, filling out the maximum lower right/upper left square/cube
// that can be created.
// The last pass uses that information to infor maximum square/cube size of the
// grid point to the nearest point occupied grid cell.
//
//

function _min3(a,b,c) {
  if ((a <= b) && (a <= c)) { return a; }
  if ((b <= a) && (b <= c)) { return b; }
  return c;
}

function grid_max_rad_2d(G) {
  let Ny = G.length;
  let Nx = G[0].length;

  let G_rad = [];

  // init,
  // b < 0  :  uninitialized
  // b = 0  :  has at least one point in grid location
  //
  for (let y=0; y<Ny; y++) {
    G_rad.push([]);
    for (let x=0; x<Nx; x++) {
      let b = ((G[y][x].length > 0) ? 0 : -1);
      G_rad[y].push({"b": b, "B": b, "r": -1});
    }
  }

  // dynamic programming ot fill in max bottomr right
  // largest square
  //
  for (let y=0; y<Ny; y++) {
    for (let x=0; x<Nx; x++) {

      if (G_rad[y][x].b >= 0) { continue; }

      let min_v = Nx*Ny;
      for (let dy=-1; dy<1; dy++) {
        for (let dx=-1; dx<1; dx++) {
          if ((dx==0) && (dy==0)) { continue; }
          if (((x+dx) < 0) || ((y+dy) < 0)) {
            min_v = 0;
            break;
          }
          if (min_v > G_rad[y+dy][x+dx].b) {
            min_v = G_rad[y+dy][x+dx].b;
          }
        }
        if (min_v <= 0) { break; }
      }
      G_rad[y][x].b = min_v+1;
    }

  }

  for (let y=(Ny-1); y>=0; y--) {
    for (let x=(Nx-1); x>=0; x--) {

      if (G_rad[y][x].B >= 0) { continue; }

      min_v = Nx*Ny;
      for (let dy=0; dy<2; dy++) {
        for (let dx=0; dx<2; dx++) {
          if ((dx==0) && (dy==0)) { continue; }
          if (((x+dx) >= Nx) || ((y+dy) >= Ny)) {
            min_v = 0;
            break;
          }
          if (min_v > G_rad[y+dy][x+dx].B) {
            min_v = G_rad[y+dy][x+dx].B;
          }
        }
        if (min_v <= 0) { break; }
      }
      G_rad[y][x].B = min_v+1;

    }
  }


  for (let y=0; y<Ny; y++) {
    for (let x=0; x<Nx; x++) {

      let b = G_rad[y][x].b;
      if (b>0) {

        let dy = -Math.floor(b/2);
        let dx = -Math.floor(b/2);
        if ((b%2)==0) {
          G_rad[y+dy][x+dx].r = ((b-1)/2);

          dy++;
          dx++;
          G_rad[y+dy][x+dx].r = ((b-1)/2);
        }
        else {
          G_rad[y+dy][x+dx].r = (b/2);
        }

      }


      let B = G_rad[y][x].B;
      if (B>0) {
        let dy = Math.floor(B/2);
        let dx = Math.floor(B/2);
        if ((B%2) == 0) {
          G_rad[y+dy][x+dx].r = ((B-1)/2);

          dy--;
          dx--;
          G_rad[y+dy][x+dx].r = ((B-1)/2);
        }
        else {
          G_rad[y+dy][x+dx].r = (B/2);
        }

      }

    }
  }


  return G_rad;
}

function grid_max_rad_3d(G) {

  let Nz = G.length;
  let Ny = G[0].length;
  let Nx = G[0][0].length;

  let G_rad = [];

  // init,
  // b < 0  :  uninitialized
  // b = 0  :  has at least one point in grid location
  //
  for (let z=0; z<Nz; z++) {
    G_rad.push([]);
    for (let y=0; y<Ny; y++) {
      G_rad[z].push([]);
      for (let x=0; x<Nx; x++) {
        let b = ((G[z][y][x].length > 0) ? 0 : -1);
        G_rad[z][y].push({"b": b, "B": b, "r": -1});
      }
    }
  }

  // dynamic programming ot fill in max bottomr right
  // largest square
  //
  for (let z=0; z<Nz; z++) {
    for (let y=0; y<Ny; y++) {
      for (let x=0; x<Nx; x++) {

        if (G_rad[z][y][x].b >= 0) { continue; }

        let min_v = Nz*Nx*Ny;
        for (let dz=-1; dz<1; dz++) {
          for (let dy=-1; dy<1; dy++) {
            for (let dx=-1; dx<1; dx++) {
              if ((dz==0) && (dx==0) && (dy==0)) { continue; }
              if (((x+dx) < 0) ||
                  ((y+dy) < 0) ||
                  ((z+dz) < 0)) {
                min_v = 0;
                break;
              }
              if (min_v > G_rad[z+dz][y+dy][x+dx].b) {
                min_v = G_rad[z+dz][y+dy][x+dx].b;
              }
            }
            if (min_v <= 0) { break; }
          }
          if (min_v <= 0) { break; }
        }
        G_rad[z][y][x].b = min_v+1;
      }

    }
  }

  for (let z=(Nz-1); z>=0; z--) {
    for (let y=(Ny-1); y>=0; y--) {
      for (let x=(Nx-1); x>=0; x--) {

        if (G_rad[z][y][x].B >= 0) { continue; }

        min_v = Nz*Nx*Ny;
        for (let dz=0; dz<2; dz++) {
          for (let dy=0; dy<2; dy++) {
            for (let dx=0; dx<2; dx++) {
              if ((dz==0) && (dx==0) && (dy==0)) { continue; }
              if (((x+dx) >= Nx) || ((y+dy) >= Ny) || ((z+dz) >= Nz)) {
                min_v = 0;
                break;
              }
              if (min_v > G_rad[z+dz][y+dy][x+dx].B) {
                min_v = G_rad[z+dz][y+dy][x+dx].B;
              }
            }
            if (min_v <= 0) { break; }
          }
          if (min_v <= 0) { break; }
        }
        G_rad[z][y][x].B = min_v+1;

      }
    }
  }


  for (let z=0; z<Nz; z++) {
    for (let y=0; y<Ny; y++) {
      for (let x=0; x<Nx; x++) {

        let b = G_rad[z][y][x].b;
        if (b>0) {

          let dz = -Math.floor(b/2);
          let dy = -Math.floor(b/2);
          let dx = -Math.floor(b/2);
          if ((b%2)==0) {
            G_rad[z+dz][y+dy][x+dx].r = ((b-1)/2);

            dz++;
            dy++;
            dx++;
            G_rad[z+dz][y+dy][x+dx].r = ((b-1)/2);
          }
          else {
            G_rad[z+dz][y+dy][x+dx].r = (b/2);
          }

        }


        let B = G_rad[z][y][x].B;
        if (B>0) {
          let dz = Math.floor(B/2);
          let dy = Math.floor(B/2);
          let dx = Math.floor(B/2);
          if ((B%2) == 0) {
            G_rad[z+dz][y+dy][x+dx].r = ((B-1)/2);

            dz--;
            dy--;
            dx--;
            G_rad[z+dz][y+dy][x+dx].r = ((B-1)/2);
          }
          else {
            G_rad[z+dz][y+dy][x+dx].r = (B/2);
          }

        }

      }
    }
  }

  return G_rad;
}

function printG3(G, gr) {
  let Nz = G.length;
  let Ny = G[0].length;
  let Nx = G[0][0].length;

  for (let z=0; z<Nz; z++) {
    for (let y=0; y<Ny; y++) {
      for (let x=0; x<Nx; x++) {
        console.log("[",x,y,z,"]:", JSON.stringify(G[z][y][x]));
      }
    }
  }
  console.log();

  for (let z=0; z<Nz; z++) {
    console.log("\n#z:", z);
    for (let y=0; y<Ny; y++) {
      let _s = [];
      for (let x=0; x<Nx; x++) {
        _s.push( (G[z][y][x].length > 0) ? 'x' : '.' );
      }
      console.log(_s.join(""));
    }
  }
  console.log();

  for (let z=0; z<Nz; z++) {
    console.log("#z:", z);
    for (let y=0; y<Ny; y++) {
      let _s = [];
      for (let x=0; x<Nx; x++) {
        _s.push( (G[z][y][x].length > 0) ? '  x' : '  .' );
      }
      console.log(_s.join(""));
    }
  }
  console.log();

  let space = "    ";
  for (let z=0; z<Ny; z++) {
    console.log("#z:", z);
    for (let y=0; y<Ny; y++) {
      let _s = [];
      for (let x=0; x<Nx; x++) {
        let b = gr[z][y][x].b;
        let bs = b.toString();
        _s.push( space.slice(bs.length) + bs );
      }
      console.log(_s.join(""));
    }
  }
  console.log("");

  space = "    ";
  for (let z=0; z<Nz; z++) {
    console.log("#z:", z);
    for (let y=0; y<Ny; y++) {
      let _s = [];
      for (let x=0; x<Nx; x++) {
        let B = gr[z][y][x].B;
        let Bs = B.toString();
        _s.push( space.slice(Bs.length) + Bs );
      }
      console.log(_s.join(""));
    }
  }
  console.log("");

  space = "    ";
  for (let z=0; z<Nz; z++) {
    console.log("#z:", z);
    for (let y=0; y<Ny; y++) {
      let _s = [];
      for (let x=0; x<Nx; x++) {
        let r = gr[z][y][x].r;
        let rs = r.toString();
        if (r<0) { rs = "."; }
        _s.push( space.slice(rs.length) + rs );
      }
      console.log(_s.join(""));
    }
  }
  console.log();


}


function printG2(G, gr) {
  let Ny = G.length;
  let Nx = G[0].length;

  for (let y=0; y<Ny; y++) {
    for (let x=0; x<Nx; x++) {
      console.log("[",x,y,"]:", JSON.stringify(G[y][x]));
    }
  }
  console.log();

  for (let y=0; y<Ny; y++) {
    let _s = [];
    for (let x=0; x<Nx; x++) {
      _s.push( (G[y][x].length > 0) ? 'x' : '.' );
    }
    console.log(_s.join(""));
  }
  console.log();

  for (let y=0; y<Ny; y++) {
    let _s = [];
    for (let x=0; x<Nx; x++) {
      _s.push( (G[y][x].length > 0) ? '  x' : '  .' );
    }
    console.log(_s.join(""));
  }
  console.log();

  let space = "    ";
  for (let y=0; y<Ny; y++) {
    let _s = [];
    for (let x=0; x<Nx; x++) {
      let b = gr[y][x].b;
      let bs = b.toString();
      _s.push( space.slice(bs.length) + bs );
    }
    console.log(_s.join(""));
  }
  console.log("");

  space = "    ";
  for (let y=0; y<Ny; y++) {
    let _s = [];
    for (let x=0; x<Nx; x++) {
      let B = gr[y][x].B;
      let Bs = B.toString();
      _s.push( space.slice(Bs.length) + Bs );
    }
    console.log(_s.join(""));
  }
  console.log("");

  space = "    ";
  for (let y=0; y<Ny; y++) {
    let _s = [];
    for (let x=0; x<Nx; x++) {
      let r = gr[y][x].r;
      let rs = r.toString();
      if (r<0) { rs = "."; }
      _s.push( space.slice(rs.length) + rs );
    }
    console.log(_s.join(""));
  }


}

function _rand_G3() {
  let N = 8*8*8;
  let Nx = Math.ceil(Math.cbrt(N)),
      Ny = Math.ceil(Math.cbrt(N)),
      Nz = Math.ceil(Math.cbrt(N));

  let G = [];

  N = 10;

  for (let z=0; z<Nz; z++) {
    G.push([]);
    for (let y=0; y<Ny; y++) {
      G[z].push([]);
      for (let x=0; x<Nx; x++) {
        G[z][y].push([]);
      }
    }
  }

  for (let i=0; i<N; i++) {
    let x = Math.random();
    let y = Math.random();
    let z = Math.random();

    let ix = Math.floor(x*Nx);
    let iy = Math.floor(y*Ny);
    let iz = Math.floor(z*Nz);

    G[iz][iy][ix].push([x,y,z]);
  }

  let gr = grid_max_rad_3d(G);

  printG3(G, gr);

}

function _rand_G2() {

  let N = 900;
  let Nx = Math.ceil(Math.sqrt(N)),
      Ny = Math.ceil(Math.sqrt(N));

  N = 100;

  let G = [];

  for (let y=0; y<Ny; y++) {
    G.push([]);
    for (let x=0; x<Nx; x++) {
      G[y].push([]);
    }
  }

  for (let i=0; i<N; i++) {
    let x = Math.random();
    let y = Math.random();

    let ix = Math.floor(x*Nx);
    let iy = Math.floor(y*Ny);

    G[iy][ix].push([x,y]);
  }


  let gr = grid_max_rad_2d(G);

  printG2(G, gr);

}

//_rand_G2();
_rand_G3();
