<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.fem

FEM (P1) solver.

## getunit()

<pre class="py-sign">triellipt.fem.<b>getunit</b>(mesh, anchors=<span>None</span>, mode=<span>None</span>)</pre>

Creates a FEM computing unit.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Input triangle mesh.
</dd></dl>

<p><span class="vardef"><code>anchors</code> : <em>Iterable = None</em></span></p>

<dl><dd>
  Provides <code>(float, float)</code> points to synchronize the mesh boundary.
</dd></dl>

<p><span class="vardef"><code>mode</code> : <em>str = None</em></span></p>

<dl><dd>
  Solver mode — <i>&quot;fvm&quot;</i> or <i>&quot;fem&quot;</i> (default).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMUnit</em></span></p>

<dl><dd>
  FEM computing unit.
</dd></dl>

## getdtn()

<pre class="py-sign">triellipt.fem.<b>getdtn</b>(mesh, anchors, mode=<span>None</span>)</pre>

Creates a FEM-DtN unit.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Input triangle mesh.
</dd></dl>

<p><span class="vardef"><code>anchors</code> : <em>Iterable</em></span></p>

<dl><dd>
  Provides <code>(float, float)</code> points to partition the mesh boundary.
</dd></dl>

<p><span class="vardef"><code>mode</code> : <em>str = None</em></span></p>

<dl><dd>
  Solver mode — <i>&quot;fvm&quot;</i> or <i>&quot;fem&quot;</i> (default).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMDtN</em></span></p>

<dl><dd>
  FEM-DtN computing unit.
</dd></dl>

<b>Notes</b>

FEM-DtN units are only supported for simply-connected meshes.

## FEMUnit

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMUnit</b>(mesh=<span>None</span>, meta=<span>None</span>)</pre>

FEM computing unit.

<b>Properties</b>

Basic FEM operators as data-streams:

Name        | Description
------------|----------------------
`massmat`   | Mass-matrix
`massdig`   | Mass-matrix lumped
`diff_1y`   | 1st y-derivative (weak)
`diff_1x`   | 1st x-derivative (weak)
`grad_1y`   | 1st y-derivative (strong)
`grad_1x`   | 1st x-derivative (strong)
`diff_2y`   | 2nd y-derivative (weak)
`diff_2x`   | 2nd x-derivative (weak)
`diff_yx`   | 2nd yx-derivative (weak)
`diff_xy`   | 2nd xy-derivative (weak)

General properties:

Name      | Description
----------|-------------------------------
`grad`    | Gradient operator (a).
`geom`    | Geometric properties (b).
`perm`    | Mesh-to-unit permutation (c).
`base`    | Basic edge-core partition.
`loops`   | List of the mesh loops.
`partts`  | Map of the unit partitions.

<b>Notes</b>

- (a) `TriGrad` as returned by `mesh_grad()`
- (b) `MeshGeom` as returned by `mesh_geom()`
- (c) `DataPerm` object with the attributes:
    - `mesh` is the parent mesh
    - `perm` is the permutation from the parent mesh

### add_partition()

<pre class="py-sign">FEMUnit.<b>add_partition</b>(<em>self</em>, partt_spec)</pre>

Adds new partition to the unit.

<b>Parameters</b>

<p><span class="vardef"><code>partt_spec</code> : <em>dict</em></span></p>

<dl><dd>
  Partition specification.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>self</em></span></p>

<dl><dd>
  Unit with the partition added.
</dd></dl>

### get_partition()

<pre class="py-sign">FEMUnit.<b>get_partition</b>(<em>self</em>, partt_name)</pre>

Fetches the unit partition.

<b>Parameters</b>

<p><span class="vardef"><code>partt_name</code> : <em>str</em></span></p>

<dl><dd>
  Partition name.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMPartt</em></span></p>

<dl><dd>
  Desired unit partition.
</dd></dl>

### del_partition()

<pre class="py-sign">FEMUnit.<b>del_partition</b>(<em>self</em>, name) → <em>None</em></pre>

Deletes the specified partition from the unit.

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
  Callable interpolator.
</dd></dl>

<b>Notes</b>

`TriInterp` object has the following attributes:

- `xnodes` contains interpolation x-nodes
- `xnodes` contains interpolation y-nodes

`TriInterp()` takes nodes-data and returns interpolated one.

### massopr()

<pre class="py-sign">FEMUnit.<b>massopr</b>(<em>self</em>, is_lumped, add_constr)</pre>

Creates the mass operator from the base partition.

<b>Parameters</b>

<p><span class="vardef"><code>is_lumped</code> : <em>bool</em></span></p>

<dl><dd>
  Creates a lumped mass operator, if <em>True</em>.
</dd></dl>

<p><span class="vardef"><code>add_constr</code> : <em>bool</em></span></p>

<dl><dd>
  Adds constraints, if <em>True</em>.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Mass operator as a matrix.
</dd></dl>

### average()

<pre class="py-sign">FEMUnit.<b>average</b>(<em>self</em>, data)</pre>

Converts node-based data to triangle-averaged values.

<b>Parameters</b>

<p><span class="vardef"><code>data</code> : <em>float-flat-array</em></span></p>

<dl><dd>
  Node-based data.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>float-flat-array</em></span></p>

<dl><dd>
  Triangle-averaged values.
</dd></dl>

## FEMDtN

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMDtN</b>(unit=<span>None</span>)</pre>

FEM-DtN computing unit.

<b>Properties</b>

Name        | Description
------------|----------------------
`fem`       | Parent FEM unit.
`dtn`       | DtN partition.

### dirich_sides()

<pre class="py-sign">FEMDtN.<b>dirich_sides</b>(<em>self</em>) → <em>list</em></pre>

Returns a list of Dirichlet sides.

### switch_side()

<pre class="py-sign">FEMDtN.<b>switch_side</b>(<em>self</em>, key)</pre>

Switches the DtN state of the boundary part.

<b>Parameters</b>

<p><span class="vardef"><code>key</code> : <em>int</em></span></p>

<dl><dd>
  Index of the boundary part.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>self</em></span></p>

<dl><dd>
  Unit itself.
</dd></dl>

## FEMPartt

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMPartt</b>(unit, meta, edge)</pre>

FEM unit partition.

<b>Properties</b>

Name        | Description
------------|-------------------------
`core`      | Core section.
`edge`      | Map of edge sections.
`meta`      | Partition metadata.

### new_vector()

<pre class="py-sign">FEMPartt.<b>new_vector</b>(<em>self</em>)</pre>

Creates a new FEM vector.

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  New empty FEM vector.
</dd></dl>

### new_matrix()

<pre class="py-sign">FEMPartt.<b>new_matrix</b>(<em>self</em>, operator, add_constr=<span>False</span>)</pre>

Creates a new FEM matrix.

<b>Parameters</b>

<p><span class="vardef"><code>operator</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Linear combination of the basic FEM operators.
</dd></dl>

<p><span class="vardef"><code>add_constr</code> : <em>bool = False</em></span></p>

<dl><dd>
  Constraints are included in the matrix, if <em>True</em>.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Resulting FEM matrix.
</dd></dl>

### get_nodes()

<pre class="py-sign">FEMPartt.<b>get_nodes</b>(<em>self</em>, key)</pre>

Retrieves the points of the partition section.

<b>Parameters</b>

<p><span class="vardef"><code>key</code> : <em>int</em></span></p>

<dl><dd>
  Number of the partition section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>two-row-float-array</em></span></p>

<dl><dd>
  Points of the partition section stacked horizontally.
</dd></dl>

<b>Notes</b>

This method is also available as `nodes2d(key)`.

## MatrixFEM

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>MatrixFEM</b>(partt=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

Global FEM matrix.

### getblock()

<pre class="py-sign">MatrixFEM.<b>getblock</b>(<em>self</em>, row_id, col_id)</pre>

Extracts a block of a matrix.

<b>Parameters</b>

<p><span class="vardef"><code>row_id</code> : <em>int</em></span></p>

<dl><dd>
  ID of the vertical section.
</dd></dl>

<p><span class="vardef"><code>col_id</code> : <em>int</em></span></p>

<dl><dd>
  ID of the horizontal section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>csc-matrix</em></span></p>

<dl><dd>
  Matrix bock in CSC format.
</dd></dl>

## VectorFEM

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>VectorFEM</b>(partt=<span>None</span>, body=<span>None</span>)</pre>

FEM vector.

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

Defines the vector via a function on the mesh nodes.

<b>Parameters</b>

<p><span class="vardef"><code>func</code> : <em>Callable</em></span></p>

<dl><dd>
  Function <code>(x, y)</code> that returns the vector body.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the body updated.
</dd></dl>

### constrained()

<pre class="py-sign">VectorFEM.<b>constrained</b>(<em>self</em>)</pre>

Constrains the vector on the parent mesh.

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the body constrained.
</dd></dl>

### getsection()

<pre class="py-sign">VectorFEM.<b>getsection</b>(<em>self</em>, sec_id)</pre>

Returns a copy of the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>sec_id</code> : <em>int</em></span></p>

<dl><dd>
  ID of the vector section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Copy of the vector section.
</dd></dl>

### setsection()

<pre class="py-sign">VectorFEM.<b>setsection</b>(<em>self</em>, sec_id, data) → <em>None</em></pre>

Defines the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>sec_id</code> : <em>int</em></span></p>

<dl><dd>
  ID of the vector section.
</dd></dl>

<p><span class="vardef"><code>data</code> : <em>scalar | flat-float-array</em></span></p>

<dl><dd>
  Data that defines the vector section.
</dd></dl>

## mesh_grad()

<pre class="py-sign">triellipt.fem.<b>mesh_grad</b>(mesh)</pre>

Returns the mesh gradient operator.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Triangular mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriGrad</em></span></p>

<dl><dd>
  Gradient operator on the mesh.
</dd></dl>

<b>Notes</b>

- `TriGrad` operates on nodes-based data
- `TriGrad.diff_1x()` returns x-grad over triangles
- `TriGrad.diff_1y()` returns y-grad over triangles
- `TriGrad()` returns xy-grad over triangles

## mesh_geom()

<pre class="py-sign">triellipt.fem.<b>mesh_geom</b>(mesh)</pre>

Returns the mesh geometric properties.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Triangular mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MeshGeom</em></span></p>

<dl><dd>
  Object with the geometric properties of triangles.
</dd></dl>

<b>Notes</b>

`MeshGeom` object has the following attributes:

- `areas` is a flat-array of triangle areas
- `sides` is a 3col-table of triangle sides

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

<b>Notes</b>

`MeshMetric` object has the following attributes:

- `bcoeffs` is a 3col-table of b-coefficients.
- `ccoeffs` is a 3col-table of c-coefficients.
- `jacobis` is a flat-array of Jacobians.