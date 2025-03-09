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
  

### Closed Venetian Pattern

Here, auxin nodes now influence all vein nodes within their relative neighborhood graph.

The relative neighborhood graph, $\text{RNG}(V)$,
of an embedded graph $(V,E)$ in some 3D Euclidean space, say,
is defined as:

$$
\begin{array}{ll}
\textbf{B} _ {0} (x, r) & = \\{ y : |x-y| < r \\} \\
\Lambda _ {p,q} & = \textbf{B} _ {0} (p, |p-q|) \ \cap \ \textbf{B} _ {0} (q, |p-q|) \\
\text{RNG}(V) & = \\{ (p,q) : (p,q) \in E , \ \Lambda _ {p,q} \ \cap V = \emptyset\\}
\end{array}
$$

$\Lambda _ {p,q}$ is called a *lune*.

For any two points, $p$ and $q$, if there are no other points in their lune, they are part
of the relative neighbor graph.

An alternative, equivalent, definition is:

$$
\text{RNG}(V) = \\{ (p,q) : (p,q) \in E, |p-q| \le \max  _ {v \in V / \{p,q\} } ( |p-v|, |q-v| ) \\}
$$

The Closed Venetian Pattern algorithm considers all vein nodes within the relative neighborhood
graph of each auxin node.
That is, instead of each auxin node linking to at most one vein node, they now link up to multiple.

Once a vein node enters the kill zone of an auxin node, the auxin node is not immediately removed
but is earmarked for future removal when all descendent vein nodes have converged into the kill zone.


#### Further Notes

Some more verbose thoughts:

* For each auxin node, $a$, calculate the RNG to the vein nodes
* Add new vein nodes if the source vein node isn't within the kill distance of the auxin
  node $a$
  - Earmark newly added vein nodes that are RNG connected to auxin node $a$ and within the
    kill distance of $a$
* Remove auxin node $a$ if all RNG connected vein nodes are within the kill distance

So the API for the RNG is:

* add points (sets up the vein nodes, puts them in grids, etc)
* add single auxin point and calculate RNG(a)
* remove single auxin point $a$

References
---

* [medium](https://medium.com/@jason.webb/space-colonization-algorithm-in-javascript-6f683b743dc5)
* [Relative Neighborhoods and Their Relatives by Jaromczyk and Toussaint](https://cgm.cs.mcgill.ca/~godfried/publications/proximity.survey.pdf)
* [Lune](https://en.wikipedia.org/wiki/Lune_%28geometry%29)
