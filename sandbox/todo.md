Todo
===

* viz
  - highlight idir frustum when picking points
  - highlight updated frame edges
  - consider displaying all q points in fence
* check/confirm frame edge calculation done correctly

done
---

* make naive method and compare on 10,100,1000, if possible
* 3d rng
* check solution (from fence 2d)
  - for each point
    + find angle between connected points
    + rotate to the +- pi/2 axis aligned point to find
      max/min point relative to each point
    + construct bounding box
    + do naive calculation relative to bounding box


