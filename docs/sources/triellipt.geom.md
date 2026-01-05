<!--
{
  "webtitle": "Modules \u2014 triellipt documentation",
  "codeblocks": false
}
-->

# triellipt.geom

Geometry module.

## line()

<pre class="py-sign">triellipt.geom.<b>line</b>(startpoint, endpoint)</pre>

Creates a line between two points.

<b>Parameters</b>

<p><span class="vardef"><code>startpoint</code> : <em>complex</em></span></p>

<dl><dd>
  Start point.
</dd></dl>

<p><span class="vardef"><code>endpoint</code> : <em>complex</em></span></p>

<dl><dd>
  End point.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>Line</em></span></p>

<dl><dd>
  Line as a curve-like object.
</dd></dl>

## hyperb()

<pre class="py-sign">triellipt.geom.<b>hyperb</b>(axes, ksis, etas)</pre>

Creates a hyperbolic curve.

<b>Parameters</b>

<p><span class="vardef"><code>axes</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Hyperbolic parameters <code>(a, b)</code>.
</dd></dl>

<p><span class="vardef"><code>ksis</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Interval of <code>ksi</code> parameter (in [-1, 1]).
</dd></dl>

<p><span class="vardef"><code>etas</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Interval of <code>eta</code> parameter (≥0).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>Hyperb</em></span></p>

<dl><dd>
  Hyperbola as a curve-like object.
</dd></dl>

## elliparc()

<pre class="py-sign">triellipt.geom.<b>elliparc</b>(center, axes, phis, tilt=0)</pre>

Creates an elliptic arc.

<b>Parameters</b>

<p><span class="vardef"><code>center</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Ellipse center.
</dd></dl>

<p><span class="vardef"><code>axes</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Major and minor ellipse axes.
</dd></dl>

<p><span class="vardef"><code>phis</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Start and end arc angles.
</dd></dl>

<p><span class="vardef"><code>tilt</code> : <em>float = 0</em></span></p>

<dl><dd>
  Ellipse rotation angle around the center. 
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>EllipArc</em></span></p>

<dl><dd>
  Elliptic arc as a curve-like object.
</dd></dl>

## bezier2()

<pre class="py-sign">triellipt.geom.<b>bezier2</b>(point0, point1, point2)</pre>

Creates a quadratic Bezier curve.

<b>Parameters</b>

<p><span class="vardef"><code>point0</code> : <em>complex</em></span></p>

<dl><dd>
  1st control point (startpoint).
</dd></dl>

<p><span class="vardef"><code>point1</code> : <em>complex</em></span></p>

<dl><dd>
  2nd control point.
</dd></dl>

<p><span class="vardef"><code>point2</code> : <em>complex</em></span></p>

<dl><dd>
  3rd control point (endpoint).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>Bezier2</em></span></p>

<dl><dd>
  Curve-like object.
</dd></dl>

## bezier3()

<pre class="py-sign">triellipt.geom.<b>bezier3</b>(point0, point1, point2, point3)</pre>

Creates a cubic Bezier curve.

<b>Parameters</b>

<p><span class="vardef"><code>point0</code> : <em>complex</em></span></p>

<dl><dd>
  1st control point (startpoint).
</dd></dl>

<p><span class="vardef"><code>point1</code> : <em>complex</em></span></p>

<dl><dd>
  2nd control point.
</dd></dl>

<p><span class="vardef"><code>point2</code> : <em>complex</em></span></p>

<dl><dd>
  3rd control point.
</dd></dl>

<p><span class="vardef"><code>point3</code> : <em>complex</em></span></p>

<dl><dd>
  4th control point (endpoint).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>Bezier3</em></span></p>

<dl><dd>
  Curve-like object.
</dd></dl>

## makeloop()

<pre class="py-sign">triellipt.geom.<b>makeloop</b>(*curves)</pre>

Creates a loop of connected curves.

<b>Parameters</b>

<p><span class="vardef"><code>curves</code> : <em>*one-of-the-curves</em></span></p>

<dl><dd>
  Sequence of connected curves.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>CurvesLoop</em></span></p>

<dl><dd>
  Loop of connected curves.
</dd></dl>

## makerect()

<pre class="py-sign">triellipt.geom.<b>makerect</b>(corner, dims)</pre>

Creates a rectangle as a loop.

<b>Parameters</b>

<p><span class="vardef"><code>corner</code> : <em>(float, float)</em></span></p>

<dl><dd>
  South-west rectangle corner.
</dd></dl>

<p><span class="vardef"><code>dims</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Width and height of the rectangle.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>CurvesLoop</em></span></p>

<dl><dd>
  Loop of rectangle sides (south-east-north-west).
</dd></dl>

## makeellip()

<pre class="py-sign">triellipt.geom.<b>makeellip</b>(center, axes, tilt=0)</pre>

Creates a closed ellipse as a loop.

<b>Parameters</b>

<p><span class="vardef"><code>center</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Ellipse center.
</dd></dl>

<p><span class="vardef"><code>axes</code> : <em>(float, float)</em></span></p>

<dl><dd>
  Major and minor axes.
</dd></dl>

<p><span class="vardef"><code>tilt</code> : <em>float = 0</em></span></p>

<dl><dd>
  Ellipse rotation angle.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>CurvesLoop</em></span></p>

<dl><dd>
  Ellipse as a single-curve loop.
</dd></dl>

## makecycle()

<pre class="py-sign">triellipt.geom.<b>makecycle</b>(path)</pre>

Creates a cyclic path.

<b>Parameters</b>

<p><span class="vardef"><code>path</code> : <em>flat-complex-array</em></span></p>

<dl><dd>
  Input path.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>CycPath</em></span></p>

<dl><dd>
  Path closed to a cycle.
</dd></dl>

## Curve

<pre class="py-sign"><b><em>class</em></b> triellipt.geom.<b>Curve</b>()</pre>

Base class for parametric curves on [0, 1].

### getpath()

<pre class="py-sign">Curve.<b>getpath</b>(<em>self</em>, *args)</pre>

Picks points on the curve.

<b>Parameters</b>

<p><span class="vardef"><code>args</code> : <em>*float</em></span></p>

<dl><dd>
  Parameters in [0, 1].
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-complex-array</em></span></p>

<dl><dd>
  Points on the curve.
</dd></dl>

### linspace()

<pre class="py-sign">Curve.<b>linspace</b>(<em>self</em>, nparts)</pre>

Splits the curve uniformly in the parameter space.

<b>Parameters</b>

<p><span class="vardef"><code>nparts</code> : <em>int</em></span></p>

<dl><dd>
  Number of intervals in the parameter space.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-complex-array</em></span></p>

<dl><dd>
  The resulting polygonal path.
</dd></dl>

### partition()

<pre class="py-sign">Curve.<b>partition</b>(<em>self</em>, nparts, ratio=1)</pre>

Splits the curve into segments based on length.

<b>Parameters</b>

<p><span class="vardef"><code>nparts</code> : <em>int</em></span></p>

<dl><dd>
  Number of segments in the partition.
</dd></dl>

<p><span class="vardef"><code>ratio</code> : <em>float = 1</em></span></p>

<dl><dd>
  Ratio of segments lengths (last-to-first).
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>flat-complex-array</em></span></p>

<dl><dd>
  The resulting polygonal chain.
</dd></dl>

### length()

<pre class="py-sign">Curve.<b>length</b>(<em>self</em>, *, places=4, maxitr=10)</pre>

Estimates the curve length iteratively.

<b>Parameters</b>

<p><span class="vardef"><code>places</code> : <em>int = 4</em></span></p>

<dl><dd>
  Number of decimal places to resolve.
</dd></dl>

<p><span class="vardef"><code>maxitr</code> : <em>int = 10</em></span></p>

<dl><dd>
  Maximum number of iterations.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>float</em></span></p>

<dl><dd>
  Length estimate.
</dd></dl>

## CurvesLoop

<pre class="py-sign"><b><em>class</em></b> triellipt.geom.<b>CurvesLoop</b>(curves=<span>None</span>)</pre>

Loop of connected curves.

<b>Attributes</b>

<p><span class="vardef"><code>curves</code> : <em>tuple</em></span></p>

<dl><dd>
  Sequence of connected curves.
</dd></dl>

<b>Properties</b>

 Name         | Description
--------------|----------------------
`startpoints` | Curves start points.
`endpoints`   | Curves end points.

### discretize()

<pre class="py-sign">CurvesLoop.<b>discretize</b>(<em>self</em>, *params)</pre>

Discretizes the curves loop.

<b>Parameters</b>

<p><span class="vardef"><code>params</code> : <em>*(nparts, ratio)</em></span></p>

<dl><dd>
  Partition parameters for each curve, see <code>Curve().partition()</code>
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>PathsMap</em></span></p>

<dl><dd>
  Polygonal path colored at curves.
</dd></dl>

## PathMap

<pre class="py-sign"><b><em>class</em></b> triellipt.geom.<b>PathMap</b>(nodes=<span>None</span>)</pre>

Polygonal path with colored nodes.

<b>Properties</b>

 Name      | Description
-----------|------------------------------
`colors`   | Nodes colors.
`numbers`  | Nodes numbers.
`points`   | Nodes positions (complex).
`points2d` | Nodes positions (xy-rows).

### togeo()

<pre class="py-sign">PathMap.<b>togeo</b>(<em>self</em>, geopath, seeds) → <em>None</em></pre>

Dumps path to the geo file.

<b>Parameters</b>

<p><span class="vardef"><code>geopath</code> : <em>str</em></span></p>

<dl><dd>
  Absolute path to the <code>.geo</code> file.
</dd></dl>

<p><span class="vardef"><code>seeds</code> : <em>dict</em></span></p>

<dl><dd>
  Maps colours to the seed mesh sizes.
</dd></dl>

### atcolors()

<pre class="py-sign">PathMap.<b>atcolors</b>(<em>self</em>, *colors)</pre>

Fetches a subpath with the specified colors.

<b>Parameters</b>

<p><span class="vardef"><code>colors</code> : <em>*int</em></span></p>

<dl><dd>
  Colors in the subpath.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>PathMap</em></span></p>

<dl><dd>
  The resulting subpath.
</dd></dl>

### repaint()

<pre class="py-sign">PathMap.<b>repaint</b>(<em>self</em>, color, newcolor)</pre>

Changes the specified color to the new one.

### rshift()

<pre class="py-sign">PathMap.<b>rshift</b>(<em>self</em>, color1, color2)</pre>

Shifts the contact of two colors to the right by one node.

### lshift()

<pre class="py-sign">PathMap.<b>lshift</b>(<em>self</em>, color1, color2)</pre>

Shifts the contact of two colors to the left by one node.

## CycPath

<pre class="py-sign"><b><em>class</em></b> triellipt.geom.<b>CycPath</b>(nodes=<span>None</span>)</pre>

Cyclic polygonal path.

<b>Properties</b>

 Name      | Description
-----------|----------------------
`points`   | Nodes positions.
`numbers`  | Nodes numbers.

### angles()

<pre class="py-sign">CycPath.<b>angles</b>(<em>self</em>)</pre>

Provides rotation angles of edges.

<b>Returns</b>

<p><span class="vardef"><em>flat-float-array</em></span></p>

<dl><dd>
  Rotation angles of edges at each node.
</dd></dl>

### dissect()

<pre class="py-sign">CycPath.<b>dissect</b>(<em>self</em>, angle)</pre>

Splits the cycle based on rotation angle.

<b>Parameters</b>

<p><span class="vardef"><code>angle</code> : <em>float</em></span></p>

<dl><dd>
  Threshold angle for a node to become a corner.
</dd></dl>

<b>Returns</b>

<p><span class="vardef"><em>PathMap</em></span></p>

<dl><dd>
  Partition of a cycle.
</dd></dl>