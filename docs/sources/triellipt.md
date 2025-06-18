<!--
{
  "webtitle": "Modules — triellipt documentation",
  "doctitle": "triellipt — Modules"
}
-->

## Annotations

Module                   | Description           
------------------------ | ----------------------
<b>triellipt.fem</b>     | Finite-element solver.
<b>triellipt.amr</b>     | Mesh refinement tools.
<b>triellipt.geom</b>    | Geometry module.      
<b>triellipt.mesher</b>  | Mesh generation tools.
<b>triellipt.trimesh</b> | Triangle mesh object. 
<b>triellipt.mshread</b> | Reader of Gmsh meshes.

## Reference

### triellipt.fem

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.fem.md#triellipt.fem">triellipt.fem</a>
        <ul>
            <li><a href="triellipt.fem.md#getunit">getunit()</a></li>
            <li><a href="triellipt.fem.md#femunit">FEMUnit</a>
                <ul>
                    <li><a href="triellipt.fem.md#add_partition">add_partition()</a></li>
                    <li><a href="triellipt.fem.md#get_partition">get_partition()</a></li>
                    <li><a href="triellipt.fem.md#set_partition">set_partition()</a></li>
                    <li><a href="triellipt.fem.md#del_partition">del_partition()</a></li>
                    <li><a href="triellipt.fem.md#getinterp">getinterp()</a></li>
                    <li><a href="triellipt.fem.md#massopr">massopr()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.fem.md#fempartt">FEMPartt</a>
                <ul>
                    <li><a href="triellipt.fem.md#new_vector">new_vector()</a></li>
                    <li><a href="triellipt.fem.md#new_matrix">new_matrix()</a></li>
                    <li><a href="triellipt.fem.md#get_nodes">get_nodes()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.fem.md#matrixfem">MatrixFEM</a>
                <ul>
                    <li><a href="triellipt.fem.md#getblock">getblock()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.fem.md#vectorfem">VectorFEM</a>
                <ul>
                    <li><a href="triellipt.fem.md#with_body">with_body()</a></li>
                    <li><a href="triellipt.fem.md#from_func">from_func()</a></li>
                    <li><a href="triellipt.fem.md#constrained">constrained()</a></li>
                    <li><a href="triellipt.fem.md#getsection">getsection()</a></li>
                    <li><a href="triellipt.fem.md#setsection">setsection()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.fem.md#mesh_grad">mesh_grad()</a></li>
            <li><a href="triellipt.fem.md#mesh_geom">mesh_geom()</a></li>
            <li><a href="triellipt.fem.md#mesh_metric">mesh_metric()</a></li>
            <li><a href="triellipt.fem.md#gettransp">gettransp()</a></li>
            <li><a href="triellipt.fem.md#transpunit">TranspUnit</a>
                <ul>
                    <li><a href="triellipt.fem.md#constr">constr()</a></li>
                    <li><a href="triellipt.fem.md#transp">transp()</a></li>
                    <li><a href="triellipt.fem.md#source">source()</a></li>
                    <li><a href="triellipt.fem.md#newdata">newdata()</a></li>
                    <li><a href="triellipt.fem.md#newcoeff">newcoeff()</a></li>
                </ul>
            </li>
        </ul>
    </li>
</ul>
</p>

### triellipt.amr

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.amr.md#triellipt.amr">triellipt.amr</a>
        <ul>
            <li><a href="triellipt.amr.md#getunit">getunit()</a></li>
            <li><a href="triellipt.amr.md#amrunit">AMRUnit</a>
                <ul>
                    <li><a href="triellipt.amr.md#refine">refine()</a></li>
                    <li><a href="triellipt.amr.md#coarsen">coarsen()</a></li>
                    <li><a href="triellipt.amr.md#find_node">find_node()</a></li>
                    <li><a href="triellipt.amr.md#find_subset">find_subset()</a></li>
                    <li><a href="triellipt.amr.md#find_masked">find_masked()</a></li>
                    <li><a href="triellipt.amr.md#front_coarse">front_coarse()</a></li>
                    <li><a href="triellipt.amr.md#front_fine">front_fine()</a></li>
                    <li><a href="triellipt.amr.md#makedata">makedata()</a></li>
                    <li><a href="triellipt.amr.md#getinterp">getinterp()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.amr.md#trifront">TriFront</a>
                <ul>
                    <li><a href="triellipt.amr.md#atrank">atrank()</a></li>
                    <li><a href="triellipt.amr.md#angles">angles()</a></li>
                    <li><a href="triellipt.amr.md#scales">scales()</a></li>
                    <li><a href="triellipt.amr.md#filter_by_mask">filter_by_mask()</a></li>
                    <li><a href="triellipt.amr.md#filter_by_angle">filter_by_angle()</a></li>
                    <li><a href="triellipt.amr.md#filter_by_scale">filter_by_scale()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.amr.md#join_meshes">join_meshes()</a></li>
        </ul>
    </li>
</ul>
</p>

### triellipt.geom

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.geom.md#triellipt.geom">triellipt.geom</a>
        <ul>
            <li><a href="triellipt.geom.md#line">line()</a></li>
            <li><a href="triellipt.geom.md#elliparc">elliparc()</a></li>
            <li><a href="triellipt.geom.md#bezier2">bezier2()</a></li>
            <li><a href="triellipt.geom.md#bezier3">bezier3()</a></li>
            <li><a href="triellipt.geom.md#makeloop">makeloop()</a></li>
            <li><a href="triellipt.geom.md#makerect">makerect()</a></li>
            <li><a href="triellipt.geom.md#makeellip">makeellip()</a></li>
            <li><a href="triellipt.geom.md#makecycle">makecycle()</a></li>
            <li><a href="triellipt.geom.md#curve">Curve</a>
                <ul>
                    <li><a href="triellipt.geom.md#getpath">getpath()</a></li>
                    <li><a href="triellipt.geom.md#linspace">linspace()</a></li>
                    <li><a href="triellipt.geom.md#partition">partition()</a></li>
                    <li><a href="triellipt.geom.md#length">length()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.geom.md#curvesloop">CurvesLoop</a>
                <ul>
                    <li><a href="triellipt.geom.md#discretize">discretize()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.geom.md#pathmap">PathMap</a>
                <ul>
                    <li><a href="triellipt.geom.md#atcolors">atcolors()</a></li>
                    <li><a href="triellipt.geom.md#repaint">repaint()</a></li>
                    <li><a href="triellipt.geom.md#rshift">rshift()</a></li>
                    <li><a href="triellipt.geom.md#lshift">lshift()</a></li>
                    <li><a href="triellipt.geom.md#split">split()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.geom.md#cycpath">CycPath</a>
                <ul>
                    <li><a href="triellipt.geom.md#angles">angles()</a></li>
                    <li><a href="triellipt.geom.md#dissect">dissect()</a></li>
                    <li><a href="triellipt.geom.md#split-1">split()</a></li>
                </ul>
            </li>
        </ul>
    </li>
</ul>
</p>

### triellipt.mesher

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.mesher.md#triellipt.mesher">triellipt.mesher</a>
        <ul>
            <li><a href="triellipt.mesher.md#trigrid">trigrid()</a></li>
            <li><a href="triellipt.mesher.md#trilattice">trilattice()</a></li>
        </ul>
    </li>
</ul>
</p>

### triellipt.trimesh

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.trimesh.md#triellipt.trimesh">triellipt.trimesh</a>
        <ul>
            <li><a href="triellipt.trimesh.md#trimesh">TriMesh</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#submesh">submesh()</a></li>
                    <li><a href="triellipt.trimesh.md#deltriangs">deltriangs()</a></li>
                    <li><a href="triellipt.trimesh.md#delghosts">delghosts()</a></li>
                    <li><a href="triellipt.trimesh.md#getvoids">getvoids()</a></li>
                    <li><a href="triellipt.trimesh.md#alignnodes">alignnodes()</a></li>
                    <li><a href="triellipt.trimesh.md#renumed">renumed()</a></li>
                    <li><a href="triellipt.trimesh.md#shuffled">shuffled()</a></li>
                    <li><a href="triellipt.trimesh.md#meshedge">meshedge()</a></li>
                    <li><a href="triellipt.trimesh.md#edgesmap">edgesmap()</a></li>
                    <li><a href="triellipt.trimesh.md#nodesmap">nodesmap()</a></li>
                    <li><a href="triellipt.trimesh.md#supertriu">supertriu()</a></li>
                    <li><a href="triellipt.trimesh.md#reduced">reduced()</a></li>
                    <li><a href="triellipt.trimesh.md#split">split()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#meshedge-1">MeshEdge</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#getloops">getloops()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#edgeloop">EdgeLoop</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#synctoedge">synctoedge()</a></li>
                    <li><a href="triellipt.trimesh.md#synctonode">synctonode()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#edgesmap-1">EdgesMap</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#getspec">getspec()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#nodesmap-1">NodesMap</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#atrank">atrank()</a></li>
                    <li><a href="triellipt.trimesh.md#atnode">atnode()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#supertriu-1">SuperTriu</a>
                <ul>
                    <li><a href="triellipt.trimesh.md#strip">strip()</a></li>
                    <li><a href="triellipt.trimesh.md#smooth">smooth()</a></li>
                    <li><a href="triellipt.trimesh.md#detach">detach()</a></li>
                    <li><a href="triellipt.trimesh.md#reduce">reduce()</a></li>
                </ul>
            </li>
            <li><a href="triellipt.trimesh.md#merge_mesh">merge_mesh()</a></li>
        </ul>
    </li>
</ul>
</p>

### triellipt.mshread

<p>
<ul class="ref-list" id="mod-refs">
    <li><a href="triellipt.mshread.md#triellipt.mshread">triellipt.mshread</a>
        <ul>
            <li><a href="triellipt.mshread.md#getreader">getreader()</a></li>
            <li><a href="triellipt.mshread.md#mshreader">MSHReader</a>
                <ul>
                    <li><a href="triellipt.mshread.md#listmeshes">listmeshes()</a></li>
                    <li><a href="triellipt.mshread.md#read_mesh">read_mesh()</a></li>
                </ul>
            </li>
        </ul>
    </li>
</ul>
</p>