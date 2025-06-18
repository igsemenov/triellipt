<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.trimesh

Triangle mesh object.

## TriMesh

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>TriMesh</b>(points=<span>None</span>, triangs=<span>None</span>)</pre>

Triangle mesh.

<b>Attributes</b>

<p><span class="vardef"><code>points</code> : <em>flat-complex-array</em></span></p>

<dl><dd>
  Mesh points in a complex plane.
</dd></dl>

<p><span class="vardef"><code>triangs</code> : <em>3-column-int-table</em></span></p>

<dl><dd>
  Triangles vertices in CCW order.
</dd></dl>

<b>Properties</b>

Name       | Description
-----------|--------------------------------------
`triu`     | Generator of triplot arguments.
`points2d` | Mesh points as two float rows.
`centrs2d` | Triangle centers as two float rows.

### submesh()

<pre class="py-sign">TriMesh.<b>submesh</b>(<em>self</em>, *trinums)</pre>

Extracts a submesh.

<b>Parameters</b>

<p><span class="vardef"><code>trinums</code> : <em>*int</em></span></p>

<dl><dd>
  Numbers of triangles to include.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  New mesh object.
</dd></dl>

### deltriangs()

<pre class="py-sign">TriMesh.<b>deltriangs</b>(<em>self</em>, *trinums)</pre>

Removes triangles from the mesh.

<b>Parameters</b>

<p><span class="vardef"><code>trinums</code> : <em>*int</em></span></p>

<dl><dd>
  Numbers of triangles to delete.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  New mesh object.
</dd></dl>

### delghosts()

<pre class="py-sign">TriMesh.<b>delghosts</b>(<em>self</em>)</pre>

Removes ghost points from the mesh.

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  New mesh.
</dd></dl>

<b>Notes</b>

Related methods:

- `.hasghosts()` shows if there are any ghosts
- `.getghosts()` returns ghost numbers, if any

### getvoids()

<pre class="py-sign">TriMesh.<b>getvoids</b>(<em>self</em>)</pre>

Finds empty triangles (voids).

<b>Returns</b>

<p><span class="vardef"><em>flat-int-array</em></span></p>

<dl><dd>
  Numbers of empty triangles.
</dd></dl>

<b>Notes</b>

Related methods:

- `.hasvoids()` shows if there are any voids
- `.delvoids()` deletes voids from the mesh

### alignnodes()

<pre class="py-sign">TriMesh.<b>alignnodes</b>(<em>self</em>, *anchors)</pre>

Performs the edge-core ordering of the mesh nodes.

<b>Parameters</b>

<p><span class="vardef"><code>anchors</code> : <em>*(float, float)</em></span></p>

<dl><dd>
  Points to synchronize the mesh boundary.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh | None</em></span></p>

<dl><dd>
  New mesh, or <em>None</em> if the mesh boundary can not be fetched.
</dd></dl>

### renumed()

<pre class="py-sign">TriMesh.<b>renumed</b>(<em>self</em>, permuter)</pre>

Renumbers the mesh nodes.

<b>Parameters</b>

<p><span class="vardef"><code>permuter</code> : <em>flat-int-array</em></span></p>

<dl><dd>
  Permutation of mesh nodes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  Mesh with the nodes renumbered.
</dd></dl>

### shuffled()

<pre class="py-sign">TriMesh.<b>shuffled</b>(<em>self</em>, permuter)</pre>

Shuffles the mesh triangles.

<b>Parameters</b>

<p><span class="vardef"><code>permuter</code> : <em>flat-int-array</em></span></p>

<dl><dd>
  Permutation of mesh triangles.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  Mesh with the triangles permuted.
</dd></dl>

### meshedge()

<pre class="py-sign">TriMesh.<b>meshedge</b>(<em>self</em>)</pre>

Extracts the mesh edge.

<b>Returns</b>

<p><span class="vardef"><em>MeshEdge</em></span></p>

<dl><dd>
  Mesh edge object.
</dd></dl>

### edgesmap()

<pre class="py-sign">TriMesh.<b>edgesmap</b>(<em>self</em>)</pre>

Maps inner mesh edges.

<b>Returns</b>

<p><span class="vardef"><em>EdgesMap</em></span></p>

<dl><dd>
  Map of inner mesh edges.
</dd></dl>

### nodesmap()

<pre class="py-sign">TriMesh.<b>nodesmap</b>(<em>self</em>)</pre>

Maps nodes to hosting triangles.

<b>Returns</b>

<p><span class="vardef"><em>NodesMap</em></span></p>

<dl><dd>
  Nodes-to-triangles map.
</dd></dl>

### supertriu()

<pre class="py-sign">TriMesh.<b>supertriu</b>(<em>self</em>)</pre>

Creates a super triangulation.

<b>Returns</b>

<p><span class="vardef"><em>SuperTriu</em></span></p>

<dl><dd>
  Super triangulation object.
</dd></dl>

### reduced()

<pre class="py-sign">TriMesh.<b>reduced</b>(<em>self</em>, shrink=<span>None</span>, detach=<span>False</span>, seed=<span>None</span>)</pre>

Tries to compress the mesh.

<b>Parameters</b>

<p><span class="vardef"><code>shrink</code> : <em>int = None</em></span></p>

<dl><dd>
  Controls shrinking of super-triangulations (i).
</dd></dl>

<p><span class="vardef"><code>detach</code> : <em>bool = False</em></span></p>

<dl><dd>
  Runs the edge detachment before compression, if <em>True</em>.
</dd></dl>

<p><span class="vardef"><code>seed</code> : <em>(float, float) = None</em></span></p>

<dl><dd>
  Seed point to start reduction.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  New mesh.
</dd></dl>

<b>Notes</b>

(i) Number of shrinking steps after one compression event.

### split()

<pre class="py-sign">TriMesh.<b>split</b>(<em>self</em>)</pre>

Splits the mesh into homogeneous parts.

<b>Returns</b>

<p><span class="vardef"><em>list</em></span></p>

<dl><dd>
  List of homogeneous submeshes.
</dd></dl>

## MeshEdge

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>MeshEdge</b>(mesh=<span>None</span>, data=<span>None</span>)</pre>

Mesh edge object.

<b>Properties</b>

Primary data:

Name       | Description
-----------|-----------------------------------
`trinums`  | Host triangles (HT) across edges.
`locnums`  | Local numbers of edges in HTs

Nodes numbers:

Name       | Description
-----------|--------------------------------
`nodnums1` | Start nodes across edges.
`nodnums2` | End nodes across edges.
`nodnums3` | Apexes across edges.

Tables:

Name      | Description
----------|-------------------------------------------
`edges2d` | Edges as two rows of nodes numbers.
`nodes2d` | Start nodes positions as two float rows.

Unique values:

Name             | Description
-----------------|-----------------------------------
`nodnums_unique` | Unique numbers of edge points.
`trinums_unique` | Unique numbers of edge triangles.

### getloops()

<pre class="py-sign">MeshEdge.<b>getloops</b>(<em>self</em>)</pre>

Splits the mesh edge into loops.

<b>Returns</b>

<p><span class="vardef"><em>list</em></span></p>

<dl><dd>
  A list of <code>EdgeLoop</code> objects.
</dd></dl>

## EdgeLoop

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>EdgeLoop</b>(mesh=<span>None</span>, data=<span>None</span>)</pre>

Loop on the mesh edge.

- Inherits `MeshEdge` properties.
- Nodes are oriented in CCW order.

### synctoedge()

<pre class="py-sign">EdgeLoop.<b>synctoedge</b>(<em>self</em>, edgeind)</pre>

Synchronizes to the specified segment.

<b>Parameters</b>

<p><span class="vardef"><code>edgeind</code> : <em>int</em></span></p>

<dl><dd>
  Index of the segment to synchronize to.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>EdgeLoop</em></span></p>

<dl><dd>
  New loop.
</dd></dl>

### synctonode()

<pre class="py-sign">EdgeLoop.<b>synctonode</b>(<em>self</em>, nodenum)</pre>

Synchronizes to the specified node.

<b>Parameters</b>

<p><span class="vardef"><code>nodenum</code> : <em>int</em></span></p>

<dl><dd>
  Global node number to synchronize to.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>EdgeLoop</em></span></p>

<dl><dd>
  New loop.
</dd></dl>

## EdgesMap

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>EdgesMap</b>(mesh=<span>None</span>, data=<span>None</span>)</pre>

Map of inner mesh edges.

<b>Properties</b>

Name       | Description
-----------|---------------------------
`trinums1` | Host triangle one (T1).
`trinums2` | Host triangle two (T2).
`locnums1` | Local edge number in T1.
`locnums2` | Local edge number in T2.

### getspec()

<pre class="py-sign">EdgesMap.<b>getspec</b>(<em>self</em>)</pre>

Classifies triangles based on edges pairing.

<b>Returns</b>

<p><span class="vardef"><em>dict</em></span></p>

<dl><dd>
  Triangles numbers for each category (i).
</dd></dl>

<b>Notes</b>

(i) Triangles categories:

- "heads" — 1 pair
- "links" — 2 pairs
- "cores" — 3 pairs
- "spots" — 0 pairs

## NodesMap

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>NodesMap</b>(mesh=<span>None</span>, data=<span>None</span>)</pre>

Nodes-to-triangles map.

<b>Properties</b>

Nodes data:

Name      | Description
----------|------------------------
`nodnums` | Global numbers.
`trinums` | Host triangles (HTs).
`locnums` | Local numbers in HTs.

Neighbors local numbers:

Name       | Description
-----------|----------------------
`locnums1` | Next CCW node.
`locnums2` | Next-next CCW node.

Neighbors global numbers:

Name       | Description
-----------|----------------------
`nodnums1` | Next CCW node.
`nodnums2` | Next-next CCW node.

### atrank()

<pre class="py-sign">NodesMap.<b>atrank</b>(<em>self</em>, rank)</pre>

Extracts data for nodes with the specified rank.

<b>Parameters</b>

<p><span class="vardef"><code>rank</code> : <em>int</em></span></p>

<dl><dd>
  Rank of the nodes to extract.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>tuple</em></span></p>

<dl><dd>
  Triplet <code>(nodnums, trinums, locnums)</code> for the extracted nodes.
</dd></dl>

### atnode()

<pre class="py-sign">NodesMap.<b>atnode</b>(<em>self</em>, nodenum)</pre>

Extracts the map of a single node.

<b>Parameters</b>

<p><span class="vardef"><code>nodenum</code> : <em>int</em></span></p>

<dl><dd>
  Number of the target node.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>NodesMap</em></span></p>

<dl><dd>
  The resulting single-node map.
</dd></dl>

## SuperTriu

<pre class="py-sign"><b><em>class</em></b> triellipt.trimesh.<b>SuperTriu</b>(mesh=<span>None</span>, data=<span>None</span>)</pre>

Super triangulation.

<b>Properties</b>

- Hosts are background triangles with three neighbors.
- Super-triangles are made from apexes of hosts neighbors.

General:

Name       | Description
-----------|----------------------------
`trinums`  | Numbers of host triangles.
`kermesh`  | Submesh of host triangles.
`supmesh`  | Mesh from super-vertices.

Neighbors:

Name       | Description
-----------|-------------------
`trinums1` | 1st CCW neighbor. 
`trinums2` | 2nd CCW neighbor.
`trinums3` | 3rd CCW neighbor.

Super-vertices:

Name       | Description
-----------|-------------------
`nodnums1` | 1st CCW vertex. 
`nodnums2` | 2nd CCW vertex.
`nodnums3` | 3rd CCW vertex.

<b>Notes</b>

Non-standard nodes pairing for `supmesh-kermesh` transition:

- (0, 1) → 1
- (1, 2) → 2
- (2, 0) → 0

### strip()

<pre class="py-sign">SuperTriu.<b>strip</b>(<em>self</em>)</pre>

Remove links from a super-triangulation.

<b>Notes</b>

See `EdgesMap.getspec()` for links definition.

### smooth()

<pre class="py-sign">SuperTriu.<b>smooth</b>(<em>self</em>, iterate=<span>True</span>)</pre>

Removes heads and spots from a super-triangulation.

<b>Parameters</b>

<p><span class="vardef"><code>iterate</code> : <em>bool = True</em></span></p>

<dl><dd>
  Runs smoothing until possible, if <em>True</em>.
</dd></dl>

<b>Notes</b>

See `EdgesMap.getspec()` for heads and spots definition.

### detach()

<pre class="py-sign">SuperTriu.<b>detach</b>(<em>self</em>)</pre>

Removes super-triangles touching the background mesh edge.

### reduce()

<pre class="py-sign">SuperTriu.<b>reduce</b>(<em>self</em>, seed=<span>None</span>, iterate=<span>True</span>)</pre>

Extracts a compact super-triangulation, if possible.

<b>Parameters</b>

<p><span class="vardef"><code>seed</code> : <em>(float, float) = None</em></span></p>

<dl><dd>
  Seed point to start reduction (b).
</dd></dl>

<p><span class="vardef"><code>iterate</code> : <em>bool = True</em></span></p>

<dl><dd>
  Triggers cleaning and retry in case of failure (a).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>SuperTriu | None</em></span></p>

<dl><dd>
  Compact super-triangulation or <em>None</em>, if failed.
</dd></dl>

<b>Notes</b>

(a) Cleaning is a strip-and-smooth action.

(b) First super-triangle is used by default.

## merge_mesh()

<pre class="py-sign">triellipt.trimesh.<b>merge_mesh</b>(omesh, imesh)</pre>

Merges a mesh from an outer part and an inner parts.

<b>Parameters</b>

<p><span class="vardef"><code>omesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Outer mesh.
</dd></dl>

<p><span class="vardef"><code>imesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Inner mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh | None</em></span></p>

<dl><dd>
  New mesh or <em>None</em>, if failed.
</dd></dl>

<b>Notes</b>

Parts must be from the same point set.