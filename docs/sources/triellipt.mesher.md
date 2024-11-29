<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.mesher

Meshing tools.

## trigrid()

<pre class="py-sign">triellipt.mesher.<b>trigrid</b>(xsize, ysize, slopes)</pre>

Creates a triangle grid.

<b>Parameters</b>

<p><span class="vardef"><code>xsize</code> : <em>int</em></span></p>

<dl><dd>
  Number of nodes in x-direction.
</dd></dl>

<p><span class="vardef"><code>ysize</code> : <em>int</em></span></p>

<dl><dd>
  Number of nodes in y-direction.
</dd></dl>

<p><span class="vardef"><code>slopes</code> : <em>str</em></span></p>

<dl><dd>
  Controls the orientation of triangles (i).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>TriMesh</em></span></p>

<dl><dd>
  New mesh.
</dd></dl>

<b>Notes</b>

(i) Triangulation in terms of grid cell division:

- "west-slope" — by west diagonals
- "east-slope" — by east diagonals
- "west-snake" — snake-like, starting from the west slope
- "east-snake" — snake-like, starting from the east slope
- "cross-wise" — by both diagonals