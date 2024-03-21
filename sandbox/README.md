Space Colonization Algorithm
===

These are some notes on the space colonization algorithm.

---

### Open Venetian Pattern

The basis of the space colonization algorithm is to start with
a list of auxin nodes in some space, 2D say, that act as attractors
for vein nodes or line segments.

Each auxin node influences at most one vein node within a bias distance,
with an auxin node being destroyed if it gets within a kill distance of
any vein.
New vein node are created in the direction of the auxin nodes,
linked to the parent vein node.

That is, auxin nodes link to a vein node that are within the influence
zone larger than the kill distance but smaller than the bias distance,
creating a new vein node whose direction is determined by the biasing
auxin nodes.

Auxin nodes are so named after the plant hormone that causes cell
growth.


I find it useful to break down the main loop of the algorithm into five parts:

* `KILL AUXIN` - kill auxin nodes within `killDist` of any vein node
* `INIT AUXIN` - initialize all auxin nodes within `biasDist` of any vein node
* `QUEUE AUXIN` - link all auxin nodes to vein nodes that are within `biasDist` and queue them for processing
* `PULL AUXIN` - go through the auxin queue, tallying the influence on each of the linked vein nodes
* `CREATE VEIN` - create new vein nodes

As an optimization, a "wave front" vein node list can be created to process only vein nodes
that fall within a bias distance of any auxin nodes.


The pseudo code is as follows:

```
killDist, biasDist, ds
auxin = RandomPoints(N)
vein, veinWF = [], []

veinWF.push( createRootNode() )

while |veinWF| > 0 {

  // kill
  foreach a in auxin within killDist of each v in veinWF
    remove a from auxin

  // init
  foreach v in veinWF that has a within biasDist
    initialize v
    initialize a

  // queue
  auxinQueue = []
  foreach v in veinWF that has a within biasDist {
    if a.curDist undefined {
      auxinQueue.push( a )
      a.curDist = |v-a|
      a.vein = v
    }
    if |v-a| < a.curDist {
      a.curDist = |v-a|
      a.vein = v
    }
  }
  
  // pull
  foreach a in auxinQueue {
    a.vein.nextPos += (a.vein.curPos-a) / |a.vein.curPos-a|
    a.vein.auxinCount++
  }

  // create
  nextVeinWF = []
  foreach v in veinWF {
    if v.auxinCount == 0 { continue }
    vNew = createVein(a.vein, (ds) * (a.vein.NextPos / v.auxinCount) )
    vein.push(vNew)
    nextVeinWF.push(vNew)
  }

  veinWF = nextVeinWF

}
```

One alteration that can be made is to give the vein an initial direction,
from its parent say, in which case the ` ... / v.auxinCount` normalization would need to
have 1 added to it.


Some points that initially caused me confusion:

* an auxin node can only have at maximum one vein node that it influences,
  taking the closest vein node that is within the bias distance
* vein nodes can have multiple auxin nodes influencing it
* `ds` is the step size
  - too small, the algorithm will run too long
  - too large, joints will be too long and other instability can occur if it
    overshoots an auxin node that's attracting it that can't be killed because
    the vein nodes are bouncing all around it
* aside from some pathological cases, auxin nodes that are within a bias distance
  of some vein should eventually be killed
  - auxin nodes will create new vein nodes that trundle along in a line until they
    either converge in on their location and are killed or vein branches are created
    as the main trunk shoots past them, with the branch then converging in to kill them
  - one pathological case can occur when a vein node is directly in the middle of two
    auxin nodes, with the inability to create a vein node that is closer to one of
    the two auxin nodes
  


