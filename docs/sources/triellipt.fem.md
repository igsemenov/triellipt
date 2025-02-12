<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.fem

Finite-element solver.

## getunit()

<pre class="py-sign">triellipt.fem.<b>getunit</b>(mesh, anchors=<span>None</span>)</pre>

Creates a FEM computing unit.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Input triangle mesh.
</dd></dl>

<p><span class="vardef"><code>anchors</code> : <em>tuple = None</em></span></p>

<dl><dd>
  Nodes numbers to synchronize the mesh boundary.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMUnit</em></span></p>

<dl><dd>
  FEM computing unit.
</dd></dl>

## FEMUnit

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMUnit</b>(mesh=<span>None</span>, meta=<span>None</span>)</pre>

FEM computing unit.

<b>Properties</b>

FEM operators as data-streams:

Name        | Description
------------|---------------------
`massmat`   | Mass-matrix
`massdiag`  | Mass-matrix lumped
`diff_1y`   | 1st-y derivative
`diff_1x`   | 1st-x derivative
`diff_2y`   | 2nd-y derivative
`diff_2x`   | 2nd-x derivative

FEM operators as FEM matrices:

Name           | Description
-------------- |----------------------------------------
`laplace_fem`  | Laplace operator
`massmat_fem`  | Mass-matrix (no constraints)
`massmat_amr`  | Mass-matrix (with constraints)
`massdiag_fem` | Mass-matrix lumped (no constraints)
`massdiag_amr` | Mass-matrix lumped (with constraints)

Others:

Name      | Description
----------|----------------------------------
`grad`    | Gradient operator.
`perm`    | Mesh-to-unit nodes permutation.
`loops`   | List of the mesh loops.

### fem_factory()

<pre class="py-sign">FEMUnit.<b>fem_factory</b>(<em>self</em>, add_constraints=<span>True</span>)</pre>

Creates a factory of FEM matrices.

<b>Parameters</b>

<p><span class="vardef"><code>add_constraints</code> : <em>bool = True</em></span></p>

<dl><dd>
  If <i>True</i>, forces constraints to be included in the matrix.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMFactory</em></span></p>

<dl><dd>
  Resulting factory of FEM matrices.
</dd></dl>

### new_vector()

<pre class="py-sign">FEMUnit.<b>new_vector</b>(<em>self</em>)</pre>

Returns a new FEM vector.

### getinterp()

<pre class="py-sign">FEMUnit.<b>getinterp</b>(<em>self</em>, xnodes, ynodes)</pre>

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
  Interpolator object.
</dd></dl>

### partition()

<pre class="py-sign">FEMUnit.<b>partition</b>(<em>self</em>, loopnum, splitang, parttspec=<span>None</span>)</pre>

Creates the FEM unit partition.

<b>Parameters</b>

<p><span class="vardef"><code>loopnum</code> : <em>int</em></span></p>

<dl><dd>
  Number of the mesh loop to partition.
</dd></dl>

<p><span class="vardef"><code>splitang</code> : <em>float</em></span></p>

<dl><dd>
  Threshold angle for the loop splitting.
</dd></dl>

<p><span class="vardef"><code>parttspec</code> : <em>list = None</em></span></p>

<dl><dd>
  List of <code>(edge1, edge2, operation)</code> triplets (a).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>Partition</em></span></p>

<dl><dd>
  Partition object.
</dd></dl>

<b>Notes</b>

(a) Color-operations on the edges:

- "l" — left-repainting
- "r" — right-repainting
- "s" — switching-of-colors

## FEMFactory

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMFactory</b>(unit=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

Factory of FEM matrices.

<b>Attributes</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Parent mesh.
</dd></dl>

<p><span class="vardef"><code>body</code> : <em>sparse-matrix</em></span></p>

<dl><dd>
  FEM-matrix pattern.
</dd></dl>

<p><span class="vardef"><code>meta</code> : <em>dict</em></span></p>

<dl><dd>
  Factory metadata.
</dd></dl>

### feed_data()

<pre class="py-sign">FEMFactory.<b>feed_data</b>(<em>self</em>, data)</pre>

Transmits data to the FEM matrix.

<b>Parameters</b>

<p><span class="vardef"><code>data</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Combination of local FEM operators (a).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Resulting FEM matrix.
</dd></dl>

<b>Notes</b>

(a) Data stream compatible with ij-stream of the FEM unit.

## MatrixFEM

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>MatrixFEM</b>(unit=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

Global FEM matrix.

### with_name()

<pre class="py-sign">MatrixFEM.<b>with_name</b>(<em>self</em>, name)</pre>

Prescribes the matrix name.

<b>Parameters</b>

<p><span class="vardef"><code>name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the matrix.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Copy of the matrix with the new name.
</dd></dl>

### dirichsplit()

<pre class="py-sign">MatrixFEM.<b>dirichsplit</b>(<em>self</em>)</pre>

Creates a Dirichlet (core-edge) partition.

### partitioned()

<pre class="py-sign">MatrixFEM.<b>partitioned</b>(<em>self</em>, meta)</pre>

Creates a partitioned matrix.

<b>Parameters</b>

<p><span class="vardef"><code>meta</code> : <em>dict</em></span></p>

<dl><dd>
  Partition meta-data (a).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Copy of the matrix with the partition defined.
</dd></dl>

<b>Notes</b>

(a) Partition meta:

- Keys are section IDs.
- Values are section indices.

### getblock()

<pre class="py-sign">MatrixFEM.<b>getblock</b>(<em>self</em>, row_key, col_key)</pre>

Extracts a block from a partitioned matrix.

<b>Parameters</b>

<p><span class="vardef"><code>row_key</code> : <em>str-or-int</em></span></p>

<dl><dd>
  ID of the vertical section.
</dd></dl>

<p><span class="vardef"><code>col_key</code> : <em>str-or-int</em></span></p>

<dl><dd>
  ID of the horizontal section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>csc-matrix</em></span></p>

<dl><dd>
  Block of a partitioned matrix.
</dd></dl>

## VectorFEM

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>VectorFEM</b>(unit=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

FEM vector.

### with_name()

<pre class="py-sign">VectorFEM.<b>with_name</b>(<em>self</em>, name)</pre>

Prescribes the vector name.

<b>Parameters</b>

<p><span class="vardef"><code>name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the vector.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the new name.
</dd></dl>

### with_body()

<pre class="py-sign">VectorFEM.<b>with_body</b>(<em>self</em>, value)</pre>

Defines the vector body.

<b>Parameters</b>

<p><span class="vardef"><code>value</code> : <em>scalar | flat-float-array</em></span></p>

<dl><dd>
  Data that defines the vector body.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the body updated.
</dd></dl>

### from_func()

<pre class="py-sign">VectorFEM.<b>from_func</b>(<em>self</em>, func)</pre>

Defines the vector body from a function.

<b>Parameters</b>

<p><span class="vardef"><code>func</code> : <em>Callable</em></span></p>

<dl><dd>
  Function with <code>(x, y)</code> call that returns the vector body.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the body updated.
</dd></dl>

### dirichsplit()

<pre class="py-sign">VectorFEM.<b>dirichsplit</b>(<em>self</em>)</pre>

Creates a Dirichlet (edge-core) partition.

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the partition defined.
</dd></dl>

### partitioned()

<pre class="py-sign">VectorFEM.<b>partitioned</b>(<em>self</em>, meta)</pre>

Creates a partitioned vector.

<b>Parameters</b>

<p><span class="vardef"><code>meta</code> : <em>dict</em></span></p>

<dl><dd>
  Partition meta data (a).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the partition defined.
</dd></dl>

<b>Notes</b>

(a) Same as for the FEM matrix.

### getsect()

<pre class="py-sign">VectorFEM.<b>getsect</b>(<em>self</em>, key)</pre>

Returns a copy of the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>key</code> : <em>str-or-int</em></span></p>

<dl><dd>
  ID of the vector section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Copy of the vector section.
</dd></dl>

### setsect()

<pre class="py-sign">VectorFEM.<b>setsect</b>(<em>self</em>, key, data)</pre>

Defines the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>key</code> : <em>str-or-int</em></span></p>

<dl><dd>
  ID of the vector section.
</dd></dl>

<p><span class="vardef"><code>data</code> : <em>scalar | flat-float-array</em></span></p>

<dl><dd>
  Data that defines the vector section.
</dd></dl>

### sectxy()

<pre class="py-sign">VectorFEM.<b>sectxy</b>(<em>self</em>, key)</pre>

Returns xy-coordinates of the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>key</code> : <em>str-or-int</em></span></p>

<dl><dd>
  ID of the vector section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>two-row-float-array</em></span></p>

<dl><dd>
  xy-coordinates of the vector section.
</dd></dl>

## mesh_metric()

<pre class="py-sign">triellipt.fem.<b>mesh_metric</b>(mesh)</pre>

Returns the mesh metric properties.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Triangular mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MeshMetric</em></span></p>

<dl><dd>
  Object with the metric properties of triangles.
</dd></dl>