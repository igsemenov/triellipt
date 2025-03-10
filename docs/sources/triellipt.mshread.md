<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.mshread

Reader of Gmsh meshes.

## getreader()

<pre class="py-sign">triellipt.mshread.<b>getreader</b>(path)</pre>

Creates a reader of Gmsh meshes.

<b>Parameters</b>

<p><span class="vardef"><code>path</code> : <em>str</em></span></p>

<dl><dd>
  Path to the folder with Gmsh meshes.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>MSHReader</em></span></p>

<dl><dd>
  Reader of Gmsh meshes.
</dd></dl>

## MSHReader

<pre class="py-sign"><b><em>class</em></b> triellipt.mshread.<b>MSHReader</b>(root_path=<span>None</span>)</pre>

Reader of Gmsh meshes.

### listmeshes()

<pre class="py-sign">MSHReader.<b>listmeshes</b>(<em>self</em>)</pre>

Returns list of `.msh` files in the root directory.

### read_mesh()

<pre class="py-sign">MSHReader.<b>read_mesh</b>(<em>self</em>, file_name)</pre>

Reads a mesh from an `.msh` file.

<b>Parameters</b>

<p><span class="vardef"><code>file_name</code> : <em>str</em></span></p>

<dl><dd>
  Name of the <code>.msh</code> file.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  Mesh object.
</dd></dl>