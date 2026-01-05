<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.amr

Mesh refinement tools.

## getunit()

<pre class="py-sign">triellipt.amr.<b>getunit</b>(mesh)</pre>

Creates an AMR unit.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Input triangle mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>AMRUnit</em></span></p>

<dl><dd>
  AMR unit with the mesh.
</dd></dl>

## AMRUnit

<pre class="py-sign"><b><em>class</em></b> triellipt.amr.<b>AMRUnit</b>(mesh)</pre>

Mesh refinement unit.

<b>Attributes</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Current triangle mesh.
</dd></dl>

<p><span class="vardef"><code>data</code> : <em>dict</em></span></p>

<dl><dd>
  Data defined on the unit.
</dd></dl>

<b>Properties</b>

Name        | Description
------------|----------------------------------
`refiner`   | Data-refiner after refinement.
`collector` | Data-collector after coarsening.

### refine()

<pre class="py-sign">AMRUnit.<b>refine</b>(<em>self</em>, trinums=<span>None</span>)</pre>

Performs a static mesh refinement.

<b>Parameters</b>

<p><span class="vardef"><code>trinums</code> : <em>Iterable = None</em></span></p>

<dl><dd>
  Numbers of triangles to refine, if <em>None</em> takes all triangles.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>AMRUnit</em></span></p>

<dl><dd>
  Unit with the refined mesh.
</dd></dl>

<b>Notes</b>

- The `data-refiner` is included in the mesh metadata.
- The void ears are not refined to keep the mesh 1-irregular.

### coarsen()

<pre class="py-sign">AMRUnit.<b>coarsen</b>(<em>self</em>, trinums_cores)</pre>

Performs a static mesh coarsening.

<b>Parameters</b>

<p><span class="vardef"><code>trinums_cores</code> : <em>Iterable</em></span></p>

<dl><dd>
  Numbers of the super-triangle-cores to coarsen.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>AMRUnit</em></span></p>

<dl><dd>
  Unit with the coarsened mesh.
</dd></dl>

<b>Notes</b>

- The `data-collector` is included in the mesh metadata.

### find_node()

<pre class="py-sign">AMRUnit.<b>find_node</b>(<em>self</em>, anchor)</pre>

Finds the neighborhood of an anchor point.

<b>Parameters</b>

<p><span class="vardef"><code>anchor</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Coordinates of an anchor point.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-int-array</em></span></p>

<dl><dd>
  Numbers of triangles near the anchor point.
</dd></dl>

### find_subset()

<pre class="py-sign">AMRUnit.<b>find_subset</b>(<em>self</em>, count, anchor, remove_heads=<span>False</span>)</pre>

Finds a convex subset of a mesh.

<b>Parameters</b>

<p><span class="vardef"><code>count</code> : <em>int</em></span></p>

<dl><dd>
  Seed number of triangles in a subset.
</dd></dl>

<p><span class="vardef"><code>anchor</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Anchor point for a starting triangle.
</dd></dl>

<p><span class="vardef"><code>remove_heads</code> : <em>bool = False</em></span></p>

<dl><dd>
  Removes single-paired triangles, if <em>True</em>.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-int-array</em></span></p>

<dl><dd>
  Numbers of triangles in a subset.
</dd></dl>

### find_masked()

<pre class="py-sign">AMRUnit.<b>find_masked</b>(<em>self</em>, mask)</pre>

Finds triangles by a mask-function.

<b>Parameters</b>

<p><span class="vardef"><code>mask</code> : <em>function</em></span></p>

<dl><dd>
  Boolean mask <code>(x, y) </code> on the triangles centroids.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-int-array</em></span></p>

<dl><dd>
  Numbers of the found triangles.
</dd></dl>

### front_coarse()

<pre class="py-sign">AMRUnit.<b>front_coarse</b>(<em>self</em>)</pre>

Finds a front of coarse triangles.

### front_fine()

<pre class="py-sign">AMRUnit.<b>front_fine</b>(<em>self</em>)</pre>

Finds a front of fine triangles.

### from_func()

<pre class="py-sign">AMRUnit.<b>from_func</b>(<em>self</em>, func)</pre>

Creates data from a function.

<b>Parameters</b>

<p><span class="vardef"><code>func</code> : <em>Callable</em></span></p>

<dl><dd>
  Function to create data.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Data defined on mesh nodes.
</dd></dl>

### constrain()

<pre class="py-sign">AMRUnit.<b>constrain</b>(<em>self</em>, data)</pre>

Constrains data on hanging nodes.

<b>Parameters</b>

<p><span class="vardef"><code>data</code> : <em>float-flat-array</em></span></p>

<dl><dd>
  Input data defined on mesh nodes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>float-flat-array</em></span></p>

<dl><dd>
  Constrained data.
</dd></dl>

### getinterp()

<pre class="py-sign">AMRUnit.<b>getinterp</b>(<em>self</em>, xnodes, ynodes)</pre>

Creates an interpolator on a mesh.

<b>Parameters</b>

<p><span class="vardef"><code>xnodes</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  x-coordinates of the interpolation nodes.
</dd></dl>

<p><span class="vardef"><code>ynodes</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  y-coordinates of the interpolation nodes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriInterp</em></span></p>

<dl><dd>
  Callable interpolator.
</dd></dl>

<b>Notes</b>

`TriInterp` object has the following attributes:

- `xnodes` contains interpolation x-nodes
- `xnodes` contains interpolation y-nodes

`TriInterp()` takes nodes-data and returns interpolated one.

## TriFront

<pre class="py-sign"><b><em>class</em></b> triellipt.amr.<b>TriFront</b>(unit=<span>None</span>, data=<span>None</span>)</pre>

Front of triangles.

<b>Properties</b>

 Name       | Description
------------|----------------------------------------
`trinums`   | Indices of the front-facing triangles.
`voidnums`  | Indices of void triangles in the front.

### atrank()

<pre class="py-sign">TriFront.<b>atrank</b>(<em>self</em>, rank)</pre>

Selects the front with the specified rank.

### angles()

<pre class="py-sign">TriFront.<b>angles</b>(<em>self</em>)</pre>

Computes the orientation angles of the front.

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Angles between the voids and the front centroids.
</dd></dl>

### scales()

<pre class="py-sign">TriFront.<b>scales</b>(<em>self</em>)</pre>

Computes the normalized front scales.

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Normalized distances between the front centroids and voids.
</dd></dl>

### filter_by_mask()

<pre class="py-sign">TriFront.<b>filter_by_mask</b>(<em>self</em>, mask)</pre>

Filters the front by the mask.

<b>Parameters</b>

<p><span class="vardef"><code>mask</code> : <em>function</em></span></p>

<dl><dd>
  Boolean mask <code>(x, y) </code> on the triangles centroids.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriFront</em></span></p>

<dl><dd>
  New front.
</dd></dl>

### filter_by_angle()

<pre class="py-sign">TriFront.<b>filter_by_angle</b>(<em>self</em>, angmin, angmax)</pre>

Filters the front by the orientation angles.

### filter_by_scale()

<pre class="py-sign">TriFront.<b>filter_by_scale</b>(<em>self</em>, minval, maxval)</pre>

Filters the front by the scales.

## join_meshes()

<pre class="py-sign">triellipt.amr.<b>join_meshes</b>(mesh1, mesh2, tol=<span>None</span>)</pre>

Join the meshes along a shared boundary, if available.

<b>Parameters</b>

<p><span class="vardef"><code>mesh1</code> : <em>TriMesh</em></span></p>

<dl><dd>
  1-st input mesh.
</dd></dl>

<p><span class="vardef"><code>mesh2</code> : <em>TriMesh</em></span></p>

<dl><dd>
  2-nd input mesh.
</dd></dl>

<p><span class="vardef"><code>tol</code> : <em>int = None</em></span></p>

<dl><dd>
  Absolute tolerance in decimal places for detecting nearby points.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh | None</em></span></p>

<dl><dd>
  New mesh or <em>None</em>, if failed.
</dd></dl>