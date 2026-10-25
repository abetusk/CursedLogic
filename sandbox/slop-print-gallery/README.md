# Print Gallery transform

Single-page WebGL2 app (`index.html`, no build, no server) that applies
Escher's "Print Gallery" logarithmic transform to an image.

References (downloaded to `refs/`, not tracked):

- J. Leys, "The mathematics behind the Droste effect",
  https://www.josleys.com/article_show.php?id=82
- B. de Smit, H. W. Lenstra Jr., "The Mathematical Structure of Escher's
  Print Gallery", Notices AMS 50(4), 2003,
  https://www.ams.org/notices/200304/fea-escher.pdf

## Modes

- **Rectangle**: the picture repeats inside a (scaled, rotated) inner frame.
- **Spiral**: for pictures with a painted spiral. Pick the centre c, a point
  A on the arm, and the matching point A' one step further in; the step
  s = (A' - c)/(A - c) (plus "extra turns") maps the arm onto itself. The
  region is any polygon that every ray from c crosses once. The log-polar
  view shows a true logarithmic spiral as a straight line. The (n1, n2)
  list under Spiral shows each continuous choice's arm pitch in the
  output (0 deg = rings) and zoom per repeat.
