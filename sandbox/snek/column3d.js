// To the extent possible under law, the person who associated CC0 with
// this project has waived all copyright and related or neighboring rights
// to this project.
// 
// You should have received a copy of the CC0 legalcode along with this
// work.  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
//


function experiment0() {
  let rho = Math.PI/6,
      n_strip = 10,
      H = 10,
      ds = 1/32;

  rho = Math.PI / 5;
  let R = 1;

  let n = Math.ceil( H / ds );

  let strip = {
    "L": [],
    "R": []
  };

  for (let strip_idx=0; strip_idx < n_strip; strip_idx++) {

    for (let sign_idx=0; sign_idx < 2; sign_idx++) {

      let strip_a = [];

      let s = ( (sign_idx == 0) ? 1 : -1 );

      for (let i=0; i<n; i++) {
        let gamma = (strip_idx / n_strip) * Math.PI * 2;
        let z = i*ds;
        let x = R*Math.cos(s * rho * z + gamma),
            y = R*Math.sin(s * rho * z + gamma);
        strip_a.push( [x,y,z] );
      }

      if (sign_idx == 0) { strip.L.push(strip_a); }
      else               { strip.R.push(strip_a); }
    }


  }

  return strip;
}

function experiment1() {
  let n_strip = 10,
      H = 10,
      dz = 1/2;

  let R = 1;

  let n = Math.ceil( H / dz );

  let strip = {
    "L": [],
    "R": []
  };

  let point_ring = [];

  for (let i=0; i<n; i++) {
    let z = dz*i;


    let point_lvl = [];

    for (let j=0; j<n_strip; j++) {
      let t = (i%2)*(Math.PI/n_strip);
      let theta = (2*Math.PI*j/n_strip) + t;

      let x = R*Math.cos( theta );
      let y = R*Math.sin( theta );

      point_lvl.push( [x,y,z] );
    }

    point_ring.push(point_lvl);
  }

  for (let i=1; i<point_ring.length; i++) {

    let m = point_ring[i].length;

    for (let j=0; j<point_ring[i].length; j++) {
      console.log( point_ring[i+0][j][0], point_ring[i+0][j][1], point_ring[i+0][j][2] );
      console.log( point_ring[i-1][j][0], point_ring[i-1][j][1], point_ring[i-1][j][2] );
      console.log("\n");

      let jj = (j+1)%m;

      if ((i%2)==0) {
        jj = (j+m-1)%m;
      }

      console.log( point_ring[i+0][j][0], point_ring[i+0][j][1], point_ring[i+0][j][2] );
      console.log( point_ring[i-1][jj][0], point_ring[i-1][jj][1], point_ring[i-1][jj][2] );
      console.log("\n");
    }

  }

}

function print_strip(strip) {
  for (let strip_idx=0; strip_idx < strip.L.length; strip_idx++) {
    for (let i=0; i<strip.L[strip_idx].length; i++) {
      console.log( strip.L[strip_idx][i][0], strip.L[strip_idx][i][1], strip.L[strip_idx][i][2]);
    }
    console.log("\n");
  }

  for (let strip_idx=0; strip_idx < strip.R.length; strip_idx++) {
    for (let i=0; i<strip.R[strip_idx].length; i++) {
      console.log( strip.R[strip_idx][i][0], strip.R[strip_idx][i][1], strip.R[strip_idx][i][2]);
    }
    console.log("\n");
  }

}

experiment1();

//print_strip(experiment0());
