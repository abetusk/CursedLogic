
var jimp = require("jimp").Jimp;

function fxy(pxl, x, y, w, h, t) {
  t = ((typeof t == "undefined") ? [128,128,128,-1] : t);

  let n = pxl.length;

  if ((x < 0) || (x >= w) ||
      (y < 0) || (y >= h)) {
    return 0;
  }

  let p = 4*(y*w + x);

  for (let i=0; i<4; i++) {
    if ((p+i) >= n) { continue; }
    if (t[i] < 0) { continue; }
    if ( pxl[p+i] < t[i] ) return 0;
  }

  return 1;
}

async function main() {

  let img = await jimp.read("data/leaf_sil0.png");

  let w = img.width;
  let h = img.height;

  let img_data = img.bitmap.data;

  /*
  let map = {};
  for (let i=0; i<img_data.length; i+=4) {
    let pxl = img_data[i+0].toString() + ":" +
              img_data[i+1].toString() + ":" +
              img_data[i+2].toString() + ":" +
              img_data[i+3].toString();
    map[pxl] = 1;
  }

  for (let key in map) {
    console.log(key);
  }
  */

  /*
  for (let y=0; y<h; y++) {
    let line = [];
    for (let x=0; x<w; x++) {

      line.push( (fxy(img_data, x,y, w,h)  < 0.5) ? '.' : 'x' );
    }
    console.log(line.join(""));
  }
  */

  let lvl_map = {};


  lvl_map[0] = {};

  lvl_map[0]['.'] = { "count": -1, "parent": {}, "child": {} };
  lvl_map[0]['x'] = { "count": -1, "parent": {}, "child": {} };
  let uniq_count = 0;

  for (let r_win=1; r_win < 10; r_win++) {

    lvl_map[r_win] = {};

    for (let y=0; y<h; y++) {
      for (let x=0; x<w; x++) {


        let key_a = [];
        for (let wy=(y-r_win); wy < (y+r_win+1); wy++) {
          for (let wx = (x-r_win); wx < (x+r_win+1); wx++) {
            key_a.push( (fxy(img_data, wx,wy, w,h) < 0.5) ? '.' : 'x' );
          }
        }

        let key = key_a.join("");
        if (!(key in lvl_map[r_win])) {
          lvl_map[r_win][key] = { "count": 0, "parent": {}, "child": {} };
          uniq_count++;
        }

        lvl_map[r_win][key].count++;

      }

    }

    for (let key in lvl_map[r_win]) {
      let sub_key_a = [];
      for (let y=1; y<(2*r_win); y++) {
        for (let x=1; x<(2*r_win); x++) {
          sub_key_a.push( key[(y*(2*r_win+1)) + x] );
        }
      }

      let sub_key = sub_key_a.join("");

      lvl_map[r_win][key].parent[sub_key] = 1;

      //console.log(key, sub_key, lvl_map[(r_win-1)][sub_key]);
      lvl_map[(r_win-1)][sub_key].child[key] = 1;

    }

    console.log(r_win);

  }



  console.log("timing test...");

  let _s = Date.now();

  let c = 0;
  for (let r_win=1; r_win < 10; r_win++) {
    for (let key in lvl_map[r_win]) {
      c++;
    }
  }

  let _e = Date.now();

  console.log("s:", _s, "e:", _e, "dt:", _e - _s);
  console.log("c:", c);

}

main();
