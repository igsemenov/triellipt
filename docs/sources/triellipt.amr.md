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
  Numbers of triangles to refine, if <i>None</i> takes all triangles.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>AMRUnit</em></span></p>

<dl><dd>
  Unit with the refined mesh.
</dd></dl>

<b>Notes</b>

- The `data-refiner` is included in the mesh metadata.
- The voids ears are not refined to keep the mesh 1-irregular.

### coarsen()

<pre class="py-sign">AMRUnit.<b>coarsen</b>(<em>self</em>, trinums_cores)</pre>

Performs a static mesh coarsening.

<b>Parameters</b>

<p><span class="vardef"><code>trinums_cores</code> : <em>Iterable</em></span></p>

<dl><dd>
  Numbers of the super-triangles-cores to coarsen.
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
  Removes single-paired triangles, if <i>True</i>.
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

### constrain()

<pre class="py-sign">AMRUnit.<b>constrain</b>(<em>self</em>, data)</pre>

Constrains data on a mesh.

<b>Parameters</b>

<p><span class="vardef"><code>data</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Nodal data to be constrained.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Constained nodal data.
</dd></dl>