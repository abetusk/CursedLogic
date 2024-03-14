
function dig2num(dig, p) {
  let b = 1/p,
    v = 0;

  for (let idx=0; idx<dig.length; idx++) {
    v += dig[idx]*b;
    b *= 1/p;
  }

  return v;
}

function num2dig(num, p) {
  let dig = [];
  num = Math.floor(num);

  if (num==0) { return [0]; }

  while (num > 0) {
    dig.push( num % p );
    num = Math.floor(num / p);
  }
  return dig;
}

function dig_incr(dig, p) {
  let carry = 1,
      idx=0;

  while (carry>0) {

    carry=0;
    if (dig.length == idx) { dig.push(0); }

    dig[idx]++;
    if (dig[idx] >= p) {
      dig[idx] -= p;
      carry=1;
    }

    idx++;
  }

  return dig;
}

function halton_sequence2d(p,q,n,skip_init,take_every) {
  n = ((typeof n === "undefined") ? 1 : n);
  skip_init = ((typeof skip_init === "undefined") ? 0 : skip_init);
  take_every = ((typeof take_every === "undefined") ? 1 : take_every);

  let dig_p = [0], dig_q = [0];

  let seq = [];

  while (seq.length != n) {

    if ((seq.length >= skip_init) &&
        ((seq.length % take_every)==0)) {
      seq.push( [ dig2num(dig_p, p), dig2num(dig_q, q) ] );
    }

    dig_incr(dig_p, p);
    dig_incr(dig_q, q);
  }

  return seq;
}

function ___xxx() {
  let d = [0,0,0,0];
  for (let i=0; i<20; i++) {
    dig_incr(d, 3);
    console.log(d, dig2num(d,3));
  }

  let seq = halton_sequence2d(2,3,1000);
  for (let ii=0; ii<seq.length; ii++) {
    console.log(seq[ii][0], seq[ii][1]);
  }

  for (let ii=0; ii<20; ii++) {
    console.log( ii, num2dig(ii, 3));
  }

}

var halton = {
  "dig2num": dig2num,
  "num2dig": num2dig,
  "dig_incr": dig_incr,
  "seq2d": halton_sequence2d
};


//export for commonJS or browser
//
if ((typeof module !== 'undefined') &&
    (typeof module.exports !== 'undefined')) {
  module.exports = halton;
}
else {
  window.halton = halton;
}


