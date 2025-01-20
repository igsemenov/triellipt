<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.amr

Mesh refinement tools.

## refine_mesh()

<pre class="py-sign">triellipt.amr.<b>refine_mesh</b>(mesh, trinums=<span>None</span>)</pre>

Performs static mesh refinement.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Triangle mesh to refine.
</dd></dl>

<p><span class="vardef"><code>trinums</code> : <em>Iterable[int] = None</em></span></p>

<dl><dd>
  Numbers of triangles to refine, if <i>None</i> takes all triangles.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  Refined mesh (a).
</dd></dl>

<b>Notes</b>

- (a) The data refiner is included in the mesh metadata.
- (b) The neighborhood of the voids is not refined.

## mesh_subset()

<pre class="py-sign">triellipt.amr.<b>mesh_subset</b>(mesh, count, anchor)</pre>

Finds a convex subset of a mesh.

<b>Parameters</b>

<p><span class="vardef"><code>mesh</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Input triangle mesh.
</dd></dl>

<p><span class="vardef"><code>count</code> : <em>int</em></span></p>

<dl><dd>
  Seed number of triangles in a subset.
</dd></dl>

<p><span class="vardef"><code>anchor</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Anchor point to find a starting triangle.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-int-array</em></span></p>

<dl><dd>
  Numbers of triangles in a subset.
</dd></dl>

## mass_collector()

<pre class="py-sign">triellipt.amr.<b>mass_collector</b>(mesh1, mesh2, mode)</pre>

Creates a mass collector for a master-slave pair of meshes.

- Only for meshes with no hanging nodes.
- No scaling by the area in a structured mode.

<b>Parameters</b>

<p><span class="vardef"><code>mesh1</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Master (parent) mesh.
</dd></dl>

<p><span class="vardef"><code>mesh2</code> : <em>TriMesh</em></span></p>

<dl><dd>
  Mesh obtained by the refinement of the parent mesh.
</dd></dl>

<p><span class="vardef"><code>mode</code> : <em>str</em></span></p>

<dl><dd>
  Defines the type of a collector (i).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MassCollector</em></span></p>

<dl><dd>
  Mass collector object.
</dd></dl>

<b>Notes</b>

(i) Collector types:

- "scaled" — for general unstructured meshes
- "structed" — for meshes with equal triangles