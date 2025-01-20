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

FEM local operators:

Name        | Description
------------|---------------------
`massmat`   | Mass-matrix
`massdiag`  | Mass-matrix lumped
`laplace`   | Laplace operator
`diff_1y`   | 1st-y derivative
`diff_1x`   | 1st-x derivative
`diff_2y`   | 2nd-y derivative
`diff_2x`   | 2nd-x derivative

FEM global matrices:

Name           | Description
-------------- |---------------------
`massmat_fem`  | Mass-matrix
`massdiag_fem` | Mass-matrix lumped
`laplace_fem`  | Laplace operator

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

### makecoeff()

<pre class="py-sign">FEMUnit.<b>makecoeff</b>(<em>self</em>, mesh_data)</pre>

Generates a coefficient for local FEM operators.

<b>Parameters</b>

<p><span class="vardef"><code>mesh_data</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Coefficient defined over the mesh triangles.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Coefficient matching the sizes of local FEM operators.
</dd></dl>

### getinterp()

<pre class="py-sign">FEMUnit.<b>getinterp</b>(<em>self</em>, xnodes, ynodes)</pre>

Creates an interpolator on a mesh.

<b>Parameters</b>

<p><span class="vardef"><code>xnodes</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  X-coordinates of the interpolation nodes.
</dd></dl>

<p><span class="vardef"><code>ynodes</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Y-coordinates of the interpolation nodes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriInterp</em></span></p>

<dl><dd>
  Interpolator object.
</dd></dl>

## FEMFactory

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMFactory</b>(unit=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

Factory of FEM matrices.

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
  Partition meta data (name and sections) (i).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Copy of the matrix with the partition defined.
</dd></dl>

<b>Notes</b>

(i) Partition map:

- Keys are the names of the sections.
- Values are flat arrays with nodes numbers.

### getblock()

<pre class="py-sign">MatrixFEM.<b>getblock</b>(<em>self</em>, rowname, colname)</pre>

Extracts a block from a partitioned matrix.

<b>Parameters</b>

<p><span class="vardef"><code>rowname</code> : <em>str</em></span></p>

<dl><dd>
  Name of the vertical section.
</dd></dl>

<p><span class="vardef"><code>colname</code> : <em>str</em></span></p>

<dl><dd>
  Name of the horizontal section.
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
  Partition meta data (name and sections) (i).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the partition defined.
</dd></dl>

<b>Notes</b>

(i) Same as for the FEM matrix.

### getsection()

<pre class="py-sign">VectorFEM.<b>getsection</b>(<em>self</em>, name)</pre>

Returns a copy of the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the vector section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Copy of the vector section.
</dd></dl>

### setsection()

<pre class="py-sign">VectorFEM.<b>setsection</b>(<em>self</em>, name, data) → <em>None</em></pre>

Defines the vector section.

<b>Parameters</b>

<p><span class="vardef"><code>name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the vector section.
</dd></dl>

<p><span class="vardef"><code>data</code> : <em>scalar | flat-float-array</em></span></p>

<dl><dd>
  Data that defines the vector section.
</dd></dl>

### sectionxy()

<pre class="py-sign">VectorFEM.<b>sectionxy</b>(<em>self</em>, name)</pre>

Returns xy-coordinates of the vector section nodes.

<b>Parameters</b>

<p><span class="vardef"><code>name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the vector section.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>two-row-float-array</em></span></p>

<dl><dd>
  xy-coordinates of the vector section.
</dd></dl>

## mesh_metric()

<pre class="py-sign">triellipt.fem.<b>mesh_metric</b>(mesh) → <em>dict</em></pre>

Returns the mesh metric properties.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Triangular mesh.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>dict</em></span></p>

<dl><dd>
  Metric properties of triangles.
</dd></dl>