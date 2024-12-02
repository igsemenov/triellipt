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

### get_factory()

<pre class="py-sign">FEMUnit.<b>get_factory</b>(<em>self</em>, with_constraints=<span>True</span>)</pre>

Creates a factory of FEM matrices.

<b>Parameters</b>

<p><span class="vardef"><code>with_constraints</code> : <em>bool = True</em></span></p>

<dl><dd>
  If <i>False</i>, constraints are not included in the matrix.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>FEMFactory</em></span></p>

<dl><dd>
  Resulting factory of FEM matrices.
</dd></dl>

### get_interp()

<pre class="py-sign">FEMUnit.<b>get_interp</b>(<em>self</em>, xnodes, ynodes)</pre>

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

### new_vector()

<pre class="py-sign">FEMUnit.<b>new_vector</b>(<em>self</em>)</pre>

Returns a new FEM vector.

## FEMFactory

<pre class="py-sign"><b><em>class</em></b> triellipt.fem.<b>FEMFactory</b>(unit=<span>None</span>, body=<span>None</span>, meta=<span>None</span>)</pre>

Factory of FEM matrices.

### feed_data()

<pre class="py-sign">FEMFactory.<b>feed_data</b>(<em>self</em>, data)</pre>

Transmits data to the matrix.

<b>Parameters</b>

<p><span class="vardef"><code>data</code> : <em>flat-float-array</em></span></p>

<dl><dd>
  Data stream compatible with ij-stream of the FEM unit.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Resulting FEM matrix.
</dd></dl>

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
  Partition meta data (name and sections).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MatrixFEM</em></span></p>

<dl><dd>
  Copy of the matrix with the partition defined.
</dd></dl>

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
  Partition meta data (name and sections).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>VectorFEM</em></span></p>

<dl><dd>
  Copy of the vector with the partition defined.
</dd></dl>

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