# Nearly Optimal Coverings of Orientation Space

Charles F. F. Karney (charles@karney.com)<br>
Version 1.0, March 6, 2006<br>
Version 1.1 (with minor revisions), October 12, 2015

## Introduction

<p>
Here we give various sets of orientations which cover orientation space
nearly optimally.  These are suitable for searching orientation space
and for integrating over orientation (together with the provided
weights).  The background to this work is given in Section 8 of
<blockquote>
Charles F. F. Karney,<br>
<a href="http://dx.doi.org/10.1016/j.jmgm.2006.04.002"><i>Quaternions
in molecular modeling</i></a>,<br>
J. Mol. Graph. Mod. <b>25</b>(5), 595&ndash;604 (Jan. 2007),<br>
Preprint:
<a href="http://arxiv.org/abs/physics/0506177">arXiv:physics/0506177</a>.
</blockquote>

<p>
Given a set of <i>N</i> orientations, we define its covering radius,
&alpha;, as the maximum amount by which an arbitrary orientation needs
to be rotated to align it with the closest member of the set.  The
coverage, <i>c</i>, is defined as
<blockquote>
    <i>c</i> = <i>N</i>(&alpha; &minus; sin &alpha;)/&pi;
</blockquote>

<p>
A set of orientations is <i>optimal</i> if
<ul>
<li> there are no other sets with the same number of elements with a
 smaller &alpha; (and <i>c</i>).
<li> all sets with a smaller number of elements have a larger &alpha;.
</ul>

<p>
For any set of <i>N</i> orientations, we can perform a Voronoi
tessellation of orientation space, associating with each member of the
set, <b>q</b><sub><i>i</i></sub>, all orientations for which
<b>q</b><sub><i>i</i></sub> is the closest orientation.  We define the
relative weight of <b>q</b><sub><i>i</i></sub> as
<blockquote>
    <i>w</i><sub><i>i</i></sub> =
    <i>N</i> (volume of Voronoi cell <i>i</i>&thinsp;) / (volume of orientation space)
</blockquote>
We can approximate an orientational average of <i>f</i>&thinsp;(<b>q</b>) with
<blockquote>
    &langle;<i>f</i>&thinsp;&rangle; =
    &sum;<sub><i>i</i></sub> <i>w</i><sub><i>i</i></sub>
    <i>f</i>&thinsp;(<b>q</b><sub><i>i</i></sub>) / <i>N</i>
</blockquote>
Assuming that the variation in <i>f</i> is bounded, we expect that, for
a given <i>N</i>, the error in this approximation to be minimized with
an optimal set of orientations.

<p>
Expressing the orientation as a unit quaternion or a pair of opposite
points on S<sup>3</sup>, we see that this problem is just a 4-dimensional
generalization of the &ldquo;spherical covering&rdquo; problem.  See
<blockquote>
    R. H. Hardin, N. J. A. Sloane, and W. D. Smith,<br>
    <a href="http://neilsloane.com/coverings/"><i>Spherical coverings</i></a>,
    (Feb. 1984),</a>
</blockquote>
with the additional constraint that the points come as opposite pairs.
The formula for <i>c</i> above involves the &ldquo;area&rdquo; of a
spherical cap on S<sup>3</sup> of radius &alpha;/2.
A quaternion [<i>q</i><sub>0</sub>, <i>q</i><sub>1</sub>,
<i>q</i><sub>2</sub>, <i>q</i><sub>3</sub>] represents the rotation give
by the matrix whose components are
<blockquote>
<table>
<tr align=center>
<td> 1 &minus; 2<i>q</i><sub>2</sub><sup>2</sup> &minus; 2<i>q</i><sub>3</sub><sup>2</sup>
<td> 2<i>q</i><sub>1</sub><i>q</i><sub>2</sub> &minus; 2<i>q</i><sub>0</sub><i>q</i><sub>3</sub>
<td> 2<i>q</i><sub>1</sub><i>q</i><sub>3</sub> + 2<i>q</i><sub>0</sub><i>q</i><sub>2</sub>
<tr align=center>
<td> 2<i>q</i><sub>2</sub><i>q</i><sub>1</sub> + 2<i>q</i><sub>0</sub><i>q</i><sub>3</sub>
<td> 1 &minus; 2<i>q</i><sub>3</sub><sup>2</sup> &minus; 2<i>q</i><sub>1</sub><sup>2</sup>
<td> 2<i>q</i><sub>2</sub><i>q</i><sub>3</sub> &minus; 2<i>q</i><sub>0</sub><i>q</i><sub>1</sub>
<tr align=center>
<td> 2<i>q</i><sub>3</sub><i>q</i><sub>1</sub> &minus; 2<i>q</i><sub>0</sub><i>q</i><sub>2</sub>
<td> 2<i>q</i><sub>3</sub><i>q</i><sub>2</sub> + 2<i>q</i><sub>0</sub><i>q</i><sub>1</sub>
<td> 1 &minus; 2<i>q</i><sub>1</sub><sup>2</sup> &minus; 2<i>q</i><sub>2</sub><sup>2</sup>
</table>
</blockquote>
where we have assumed that the quaternion is normalized, i.e.,
<i>q</i><sub>0</sub><sup>2</sup> + <i>q</i><sub>1</sub><sup>2</sup> +
<i>q</i><sub>2</sub><sup>2</sup> + <i>q</i><sub>3</sub><sup>2</sup> = 1.

<p>
The problem of determining good orientation sets for the purposes of
averaging is discussed in
<blockquote>
M. Ed&eacute;n and M. H. Levitt,<br>
<a href="http://dx.doi.org/10.1006/jmre.1998.1427"><i>Computation of orientational
averages in solid state NMR by Gaussian spherical quadrature</i></a>,<br>
J. Magn. Reson. <b>132</b>, 220&ndash;239 (1998).<br>
</blockquote>

## Table of orientation sets

<p>
Because determining optimal sets of points is a hard problem, we provide
here &ldquo;nearly&rdquo; optimal sets of points.  We begin by providing
a table of the orientation sets ranked by decreasing &alpha;.

<p>
<center>
<table>
<thead align=left>
<tr>
  <th>name
  <th align=center><i>N</i>
  <th align=center>&alpha; (&deg;)
  <th align=center><i>c</i>
  <th align=center>&delta;
  <th align=center>&sigma;
  <th>download
<tbody align=left>
<tr>
<td>c48u1
<td align=right> 24
<td align=right> 62.80
<td> 1.57514
<td> 0.70000
<td> 0.00
<td>
 <a href="data/c48u1.quat">quat</a>
 <a href="data/c48u1.grid">grid</a>
 <a href="data/c48u1.euler">euler</a>
<tr>
<td>c600v
<td align=right> 60
<td align=right> 44.48
<td> 1.44480
<td>
<td>
<td>
 <a href="data/c600v.quat">quat</a>
 <a href="data/c600v.euler">euler</a>
<tr>
<td>c48u9
<td align=right> 216
<td align=right> 38.45
<td> 3.38698
<td> 0.41422
<td> 0.00
<td>
 <a href="data/c48u9.quat">quat</a>
 <a href="data/c48u9.grid">grid</a>
 <a href="data/c48u9.euler">euler</a>
<tr>
<td>c48n9
<td align=right> 216
<td align=right> 36.47
<td> 2.89689
<td> 0.26091
<td> 7.00
<td>
 <a href="data/c48n9.quat">quat</a>
 <a href="data/c48n9.grid">grid</a>
 <a href="data/c48n9.euler">euler</a>
<tr>
<td>c600vc
<td align=right> 360
<td align=right> 27.78
<td> 2.15246
<td>
<td>
<td>
 <a href="data/c600vc.quat">quat</a>
 <a href="data/c600vc.euler">euler</a>
<tr>
<td>c600vec
<td align=right> 720
<td align=right> 22.25
<td> 2.22117
<td>
<td>
<td>
 <a href="data/c600vec.quat">quat</a>
 <a href="data/c600vec.euler">euler</a>
<tr>
<td>c48u27
<td align=right> 648
<td align=right> 20.83
<td> 1.64091
<td> 0.33582
<td> 0.00
<td>
 <a href="data/c48u27.quat">quat</a>
 <a href="data/c48u27.grid">grid</a>
 <a href="data/c48u27.euler">euler</a>
<tr>
<td>c48u83
<td align=right> 1992
<td align=right> 16.29
<td> 2.42065
<td> 0.25970
<td> 0.00
<td>
 <a href="data/c48u83.quat">quat</a>
 <a href="data/c48u83.grid">grid</a>
 <a href="data/c48u83.euler">euler</a>
<tr>
<td>c48u157
<td align=right> 3768
<td align=right> 14.49
<td> 3.22614
<td> 0.20710
<td> 0.00
<td>
 <a href="data/c48u157.quat">quat</a>
 <a href="data/c48u157.grid">grid</a>
 <a href="data/c48u157.euler">euler</a>
<tr>
<td>c48u181
<td align=right> 4344
<td align=right> 12.29
<td> 2.27013
<td> 0.19415
<td> 0.00
<td>
 <a href="data/c48u181.quat">quat</a>
 <a href="data/c48u181.grid">grid</a>
 <a href="data/c48u181.euler">euler</a>
<tr>
<td>c48u309
<td align=right> 7416
<td align=right> 10.07
<td> 2.13338
<td> 0.15846
<td> 0.00
<td>
 <a href="data/c48u309.quat">quat</a>
 <a href="data/c48u309.grid">grid</a>
 <a href="data/c48u309.euler">euler</a>
<tr>
<td>c48n309
<td align=right> 7416
<td align=right> 9.72
<td> 1.91567
<td> 0.15167
<td> 1.86
<td>
 <a href="data/c48n309.quat">quat</a>
 <a href="data/c48n309.grid">grid</a>
 <a href="data/c48n309.euler">euler</a>
<tr>
<td>c48u519
<td align=right> 12456
<td align=right> 9.05
<td> 2.60257
<td> 0.13807
<td> 0.00
<td>
 <a href="data/c48u519.quat">quat</a>
 <a href="data/c48u519.grid">grid</a>
<tr>
<td>c48u527
<td align=right> 12648
<td align=right> 8.43
<td> 2.13318
<td> 0.13229
<td> 0.00
<td>
 <a href="data/c48u527.quat">quat</a>
 <a href="data/c48u527.grid">grid</a>
<tr>
<td>c48n527
<td align=right> 12648
<td align=right> 8.17
<td> 1.94334
<td> 0.12599
<td> 1.86
<td>
 <a href="data/c48n527.quat">quat</a>
 <a href="data/c48n527.grid">grid</a>
<tr>
<td>c48u815
<td align=right> 19560
<td align=right> 7.40
<td> 2.23719
<td> 0.11607
<td> 0.00
<td>
 <a href="data/c48u815.quat">quat</a>
 <a href="data/c48u815.grid">grid</a>
<tr>
<td>c48u1153
<td align=right> 27672
<td align=right> 6.60
<td> 2.23735
<td> 0.10330
<td> 0.00
<td>
 <a href="data/c48u1153.quat">quat</a>
 <a href="data/c48u1153.grid">grid</a>
<tr>
<td>c48u1201
<td align=right> 28824
<td align=right> 6.48
<td> 2.20918
<td> 0.09999
<td> 0.00
<td>
 <a href="data/c48u1201.quat">quat</a>
 <a href="data/c48u1201.grid">grid</a>
<tr>
<td>c48u1641
<td align=right> 39384
<td align=right> 5.75
<td> 2.10646
<td> 0.08993
<td> 0.00
<td>
 <a href="data/c48u1641.quat">quat</a>
 <a href="data/c48u1641.grid">grid</a>
<tr>
<td>c48u2219
<td align=right> 53256
<td align=right> 5.27
<td> 2.20117
<td> 0.08249
<td> 0.00
<td>
 <a href="data/c48u2219.quat">quat</a>
 <a href="data/c48u2219.grid">grid</a>
<tr>
<td>c48u2867
<td align=right> 68808
<td align=right> 5.24
<td> 2.79649
<td> 0.07531
<td> 0.00
<td>
 <a href="data/c48u2867.quat">quat</a>
 <a href="data/c48u2867.grid">grid</a>
<tr>
<td>c48u2947
<td align=right> 70728
<td align=right> 4.71
<td> 2.07843
<td> 0.07359
<td> 0.00
<td>
 <a href="data/c48u2947.quat">quat</a>
 <a href="data/c48u2947.grid">grid</a>
<tr>
<td>c48u3733
<td align=right> 89592
<td align=right> 4.37
<td> 2.11197
<td> 0.06836
<td> 0.00
<td>
 <a href="data/c48u3733.quat">quat</a>
 <a href="data/c48u3733.grid">grid</a>
<tr>
<td>c48u4701
<td align=right> 112824
<td align=right> 4.22
<td> 2.39041
<td> 0.06372
<td> 0.00
<td>
 <a href="data/c48u4701.quat">quat</a>
 <a href="data/c48u4701.grid">grid</a>
<tr>
<td>c48u4749
<td align=right> 113976
<td align=right> 4.00
<td> 2.05300
<td> 0.06248
<td> 0.00
<td>
 <a href="data/c48u4749.quat">quat</a>
 <a href="data/c48u4749.grid">grid</a>
<tr>
<td>c48u5879
<td align=right> 141096
<td align=right> 3.74
<td> 2.07325
<td> 0.05837
<td> 0.00
<td>
 <a href="data/c48u5879.quat">quat</a>
 <a href="data/c48u5879.grid">grid</a>
<tr>
<td>c48u7111
<td align=right> 170664
<td align=right> 3.53
<td> 2.11481
<td> 0.05514
<td> 0.00
<td>
 <a href="data/c48u7111.quat">quat</a>
 <a href="data/c48u7111.grid">grid</a>
<tr>
<td>c48u8649
<td align=right> 207576
<td align=right> 3.26
<td> 2.02898
<td> 0.05094
<td> 0.00
<td>
 <a href="data/c48u8649.quat">quat</a>
 <a href="data/c48u8649.grid">grid</a>
</table>
</center>

<p>
The orientation sets can be downloaded by the links in the
&ldquo;download&rdquo; column.  The &ldquo;quat&rdquo; and
&ldquo;euler&rdquo; files are in the following format:
<ul>
<li>Any number of initial comment lines beginning with &ldquo;#&rdquo;.
<li>A line containing either &ldquo;format quaternion&rdquo; or
    &ldquo;format euler&rdquo;.
<li>A line containing: <i>N</i> &alpha;(&deg;) <i>c</i>.
<li><i>N</i> lines containing:
    <i>q</i><sub>0<i>i</i></sub>
    <i>q</i><sub>1<i>i</i></sub>
    <i>q</i><sub>2<i>i</i></sub>
    <i>q</i><sub>3<i>i</i></sub>
    <i>w</i><sub><i>i</i></sub> (for quaternions)
<li><i>or</i> <i>N</i> lines containing:
    <i>a</i><sub><i>i</i></sub>
    <i>b</i><sub><i>i</i></sub>
    <i>g</i><sub><i>i</i></sub>
    <i>w</i><sub><i>i</i></sub> (for Euler angles).
</ul>
The convention for Euler angles is that the rotation is given by
<i>R</i><sub><i>z</i></sub>(<i>a</i>)
<i>R</i><sub><i>y</i></sub>(<i>b</i>)
<i>R</i><sub><i>z</i></sub>(<i>g</i>).

<p>
The &ldquo;grid&rdquo; links provide a compact representation for those
orientation sets based on the 48-cell with &delta; and &sigma;
determining the grid spacing.  This is described below.
Only a subset of orientations sets in Euler angle format is provided
here.

<p>
The following orientation sets are non-optimal:
<ul>
<li>
c48u9 (beaten by c48n9),
<li>
c600vec (beaten by c48u27),
<li>
c48u309 (beaten by c48n309),
<li>
c48u527 (beaten by c48n527).
</ul>

<p>
The following orientation sets are sub-optimal with a substantially
thinner covering achieved by another set with somewhat more points:
<ul>
<li>
c48u157 (use c48u181 instead),
<li>
c48u519 (use c48n527 instead),
<li>
c48u2867 (use c48u2947 instead),
<li>
c48u4701 (use c48u4749 instead).
</ul>

## Sets based on regular and semi-regular polytopes

<p>
One strategy for evenly spacing points on S<sup>3</sup> is to place
the points using the vertices or cell centers of regular and
semi-regular polytopes.  The vertices of all the regular 4-dimensional
polytopes is given in
<blockquote>
<a href="http://paulbourke.net/geometry/hyperspace/">http://paulbourke.net/geometry/hyperspace/</a>
</blockquote>
c48u1 and c660v are two such sets.  The
points in c48u1 are placed at the centers of the cells of a
truncated-cubic tetracontoctachoron (48-cell), see
<blockquote>
<a href="https://en.wikipedia.org/wiki/Truncated_24-cells#Bitruncated_24-cell">https://en.wikipedia.org/wiki/Tetracontoctachoron</a>.
</blockquote>
These points are obtained by using the vertices of 2 24-cells in their
mutually dual positions.  Similarly the points in c600v are the
vertices of a 600-cell.  Both c48u1 and c600v probably <i>are</i>
optimal sets.

The set c600v can be extended by adding the centers of the cells of
the 600-cell (equivalent to the vertices of its dual, the 120-cell) to
give the set c600vc and by adding, in addition, the midpoints of the
edges of the 600-cell, to give the set c600vec.

## Sets based on gridding the 48-cell

<p>
In order to obtain larger sets we seek a systematic way to place
multiple points with the cells of a polytope.  The 48-cell is
convenient for this purpose.  The cells are all identical truncated
cubes and thus a body-center-cubic lattice lines up nicely with the
cells.  [A body-center-cubic lattice provides the thinnest covering of
R<sup>3</sup>; see R. P. Bambah, <i>On lattice coverings by spheres</i>,
Proc. Nat. Inst. Sci. India <b>20</b>, 25&ndash;52 (1954).]

<p>
Here is the procedure.  Each cell of the 48-cell is a truncated cube.
Define the primary cell as
<blockquote>
    <i>p</i><sub>0</sub> = 1,<br>
    <i>p</i><sub><i>i</i>&gt;0</sub> &lt; &radic;2 &minus; 1,<br>
    |<i>p</i><sub>1</sub>| +
    |<i>p</i><sub>2</sub>| +
    |<i>p</i><sub>3</sub>| &lt; 1.
</blockquote>
The other cells are generated from this by the application of the
rotational symmetry group of the cube.

<p>
Place a body-centered-cubic lattice, with lattice spacing &delta;, within
the primary cell (including only points lying within the cell).  Thus we
take
<blockquote>
   <i>p</i><sub>0</sub> = 1,<br>
   [<i>p</i><sub>1</sub>, <i>p</i><sub>2</sub>, <i>p</i><sub>3</sub>] =
    [<i>k</i>, <i>l</i>, <i>m</i>] &delta;/2
</blockquote>
where [<i>k</i>, <i>l</i>, <i>m</i>] are either all even or all odd integers
(to give a BCC lattice).  These points are then normalized with
<b>q</b> = <b>p</b>/|<b>p</b>| to place them on S<sup>3</sup>.

<p>
As &delta; is varied the number of points within the cell (<i>N</i>/24)
varies.  For a given <i>N</i>, pick the &delta; with the smallest
covering.  (To obtain the sets given here, we systematically varied
&delta; in steps of 0.00001.)  Discard any <i>N</i> for which there is a
smaller <i>N</i> with a smaller (or equal) &alpha;.

<p>
This procedure yields the sets c48uMMM where MMM =
<i>N</i><sub><i>c</i></sub> = <i>N</i>/24 is the
number of points per cell.

<p>
There are many ways in which we might imagine improving these sets.
One possibility is to use a non-uniform lattice spacing using
<blockquote>
   <i>p</i><sub>1</sub> = sinh(&sigma; <i>k</i> &delta;/2)/&sigma;,
</blockquote>
and similarly for <i>p</i><sub>2</sub> and <i>p</i><sub>3</sub>.
(The uniform lattice is recovered in the limit &sigma; &rarr; 0.)  The
increasing lattice spacing afforded by the sinh function counteracts
the bunching of points occurring when the lattice points are projected
onto S<sup>3</sup>.)  The two sets c48n309 and c48n527 are two examples
with reasonably thin coverage.  In the case of c48n9, &sigma; is used
merely to delay the entry of a new set of points into the primary cell.

<p>
One might also offset the lattice and remove or perturb the points near
the surface of the cells.  However because these strategies make the
search for good sets considerably more complex, the simple procedure
with the uniform lattice outlined above probably suffices for most
purposes.

<p>
Because of the regular way that the grids are obtained, we can define a
compact represtentation of the orientation set with a file in the
&ldquo;grid&rdquo; format.  This consists of
<ul>
<li>Any number of initial comment lines beginning with &ldquo;#&rdquo;.
<li>A line containing &ldquo;format grid&rdquo;.
<li>A line containing:
    &delta;
    &sigma;
    <i>N</i>
    <i>N</i><sub><i>c</i></sub>
    <i>N</i><sub><i>d</i></sub>
    &alpha;(&deg;)
    <i>c</i>.
<li><i>N</i><sub><i>d</i></sub> lines containing:
    <i>k</i><sub><i>i</i></sub>
    <i>l</i><sub><i>i</i></sub>
    <i>m</i><sub><i>i</i></sub>
    <i>w</i><sub><i>i</i></sub>
    <i>r</i><sub><i>i</i></sub>(&deg;)
    <i>M</i><sub><i>i</i></sub>.
</ul>
where <i>N</i><sub><i>d</i></sub> is the number of distinct entries,
<i>r</i><sub><i>i</i></sub> is the radius of the
<i>i</i>&thinsp;th Voronoi cell [thus &alpha; =
max<sub><i>i</i></sub>(<i>r</i><sub><i>i</i></sub>)], and
<i>M</i><sub><i>i</i></sub> is the multiplicity of the entry.  In the
file, we restrict <i>k</i><sub><i>i</i></sub> &ge;
<i>l</i><sub><i>i</i></sub> &ge; <i>m</i><sub><i>i</i></sub> &ge; 0.  For
each such [<i>k</i><sub><i>i</i></sub>, <i>l</i><sub><i>i</i></sub>,
<i>m</i><sub><i>i</i></sub>], we generate <i>M</i><sub><i>i</i></sub>
distinct permutations by changing the order and the signs of the
elements.

<p>
Code to produce the full orientation sets for the grid form is available
in <a href="ExpandSet.cpp">ExpandSet.cpp</a>.  After compiling this
code, you can generate a quaternion orientation set with, e.g.,
<blockquote>
   ./ExpandSet &lt; c48u527.grid &gt; c48u527.quat
</blockquote>
Supply the &ldquo;-e&rdquo; option to obtain the corresponding file of
Euler angles.

<p>
Additional denser orientation sets are provided below.  Because these
sets contain a large number of orientations (up to
25&times;10<sup>6</sup>), they are provided only the grid format.

## Further remarks

<p>
For each set, can obtain new sets by performing an arbitrary rotation of
R<sup>4</sup> via
<blockquote>
    <b>q</b><sub><i>i</i></sub>&prime; =
    <b>r</b> <b>q</b><sub><i>i</i></sub> <b>s</b>,
</blockquote>
where <b>r</b> and <b>s</b> are fixed (possibly random) unit
quaternions.  The pre- and post-multiplication allows all rotations of
R<sup>4</sup> to be accessed.

<p>
One way of estimating the error in the numerical quadrature is to repeat
the calculation several times with the same set of points but choosing
different random <b>r</b> and <b>s</b> each time.

<p>
(Note that original sets possess symmetry that if
<b>q</b><sub><i>i</i></sub> is a member of the set then so is the
inverse rotation <b>q</b><sub><i>i</i></sub>*.  The new sets
<b>q</b><sub><i>i</i></sub>&prime; do not have this property, in general.)

<p>
Estimated accuracy (ulp = units in last place):
<ul>
<li>&delta;: exact
  (search for &ldquo;optimum&rdquo; was with resolution 10<sup>&minus;5</sup>)
<li>&sigma;: exact
<li>&alpha;: 0.006&deg; (0.6 ulp)
<li><i>c</i>: 0.6&times;10<sup>&minus;5</sup> (0.6 ulp)
<li><i>w</i>: average 1.5&times;10<sup>&minus;6</sup> (1 ulp),
        maximum 4&times;10<sup>&minus;6</sup> (4 ulp)
        (last digit adjusted to give
         &sum;<sub><i>i</i></sub> <i>w</i><sub><i>i</i></sub> = <i>N</I>)
<li><i>r</i>: 0.015&deg; (1.5 ulp)
<li><i>q</i><sub>0</sub>, <i>q</i><sub>1</sub>,
    <i>q</i><sub>2</sub>, <i>q</i><sub>3</sub>:
    0.51&times;10<sup>&minus;9</sup> (0.51 ulp)
</ul>

## ZCW3 Orientation Sets

<p>
Ed&eacute;n and Levitt studied the ZCW3 orientation sets.
These are based on gridding the space of Euler angles.  These yield
less thin coverings of orientation space that the sets given above.
Here is the data

<p>
<center>
<table>
<thead align=left>
<tr>
  <th>name
  <th align=center><i>N</i>
  <th align=center>&alpha; (&deg;)
  <th align=center><i>c</i>
<tbody align=left>
<tr>
<td> ZCW3_50
<td align=right> 50
<td> 69.66
<td align=right> 4.426
<tr>
<td> ZCW3_100
<td align=right> 100
<td> 56.05
<td align=right> 4.735
<tr>
<td> ZCW3_144
<td align=right> 144
<td> 42.44
<td align=right> 3.021
<tr>
<td> ZCW3_200
<td align=right> 200
<td> 48.07
<td align=right> 6.050
<tr>
<td> ZCW3_300
<td align=right> 300
<td> 40.25
<td align=right> 5.384
<tr>
<td> ZCW3_538
<td align=right> 538
<td> 32.53
<td align=right> 5.142
<tr>
<td>ZCW3_1154
<td align=right> 1154
<td> 26.81
<td align=right> 6.203
<tr>
<td>ZCW3_3722
<td align=right> 3722
<td> 18.33
<td align=right> 6.436
<tr>
<td>ZCW3_6044
<td align=right> 6044
<td> 18.10
<td align=right> 10.051
</table>
</center>

## Denser orientation sets

<p>
The procedure used to obtain the orientation sets based on the 48-cell
can be continued to obtain denser orientation sets.  Here are the
results for uniform grids (&sigma; = 0):

<p>
<center>
<table>
<thead align=left>
<tr>
  <th>name
  <th align=center><i>N</i>
  <th align=center>&alpha; (&deg;)
  <th align=center><i>c</i>
  <th align=center>&delta;
  <th align=center>approx &delta;
  <th>download
<tbody align=left>
<tr>
<td>
<td align=right> 141096
<td> 3.735
<td> 2.07261
<td> 0.058364
<td> 2/34.2973
<td>
<tr>
<td>
<td align=right> 170664
<td> 3.529
<td> 2.11458
<td> 0.055138
<td> 2/36.2973
<td>
<tr>
<td>
<td align=right> 207576
<td> 3.260
<td> 2.02803
<td> 0.050932
<td> 2<i>s</i>/16.2657
<td>
<tr>
<td> c48u10305
<td align=right> 247320
<td> 3.102
<td> 2.08130
<td> 0.048456
<td> 2/41.2973
<td><a href="data/c48u10305.grid">grid</a>
<tr>
<td> c48u12083
<td align=right> 289992
<td> 3.096
<td> 2.42678
<td> 0.046023
<td>
<td><a href="data/c48u12083.grid">grid</a>
<tr>
<td> c48u12251
<td align=right> 294024
<td> 2.903
<td> 2.02950
<td> 0.045354
<td> 2<i>s</i>/18.2657
<td><a href="data/c48u12251.grid">grid</a>
<tr>
<td> c48u14251
<td align=right> 342024
<td> 2.767
<td> 2.04269
<td> 0.043215
<td> 2/46.2973
<td><a href="data/c48u14251.grid">grid</a>
<tr>
<td> c48u16533
<td align=right> 396792
<td> 2.655
<td> 2.09385
<td> 0.041421
<td> 2/48.2973
<td><a href="data/c48u16533.grid">grid</a>
<tr>
<td> c48u19181
<td align=right> 460344
<td> 2.497
<td> 2.02149
<td> 0.039000
<td> 2<i>s</i>/21.2450
<td><a href="data/c48u19181.grid">grid</a>
<tr>
<td> c48u21863
<td align=right> 524712
<td> 2.403
<td> 2.05419
<td> 0.037534
<td> 2/53.2973
<td><a href="data/c48u21863.grid">grid</a>
<tr>
<td> c48u25039
<td align=right> 600936
<td> 2.282
<td> 2.01458
<td> 0.035641
<td> 2<i>s</i>/23.2450
<td><a href="data/c48u25039.grid">grid</a>
<tr>
<td> c48u28329
<td align=right> 679896
<td> 2.197
<td> 2.03407
<td> 0.034313
<td> 2/58.2973
<td><a href="data/c48u28329.grid">grid</a>
<tr>
<td> c48u31793
<td align=right> 763032
<td> 2.162
<td> 2.17361
<td> 0.033137
<td>
<td><a href="data/c48u31793.grid">grid</a>
<tr>
<td> c48u32081
<td align=right> 769944
<td> 2.116
<td> 2.05852
<td> 0.032786
<td>
<td><a href="data/c48u32081.grid">grid</a>
<tr>
<td> c48u35851
<td align=right> 860424
<td> 2.024
<td> 2.01113
<td> 0.031601
<td> 2/63.2973
<td><a href="data/c48u35851.grid">grid</a>
<tr>
<td> c48u40003
<td align=right> 960072
<td> 1.962
<td> 2.04420
<td> 0.030633
<td> 2/65.2973
<td><a href="data/c48u40003.grid">grid</a>
<tr>
<td> c48u44709
<td align=right> 1073016
<td> 1.877
<td> 2.00081
<td> 0.029307
<td> 2<i>s</i>/28.2657
<td><a href="data/c48u44709.grid">grid</a>
<tr>
<td> c48u49397
<td align=right> 1185528
<td> 1.822
<td> 2.02304
<td> 0.028453
<td> 2/70.2973
<td><a href="data/c48u49397.grid">grid</a>
<tr>
<td> c48u54799
<td align=right> 1315176
<td> 1.753
<td> 1.99776
<td> 0.027370
<td> 2<i>s</i>/30.2657
<td><a href="data/c48u54799.grid">grid</a>
<tr>
<td> c48u60279
<td align=right> 1446696
<td> 1.701
<td> 2.00892
<td> 0.026563
<td> 2/75.2973
<td><a href="data/c48u60279.grid">grid</a>
<tr>
<td> c48u65985
<td align=right> 1583640
<td> 1.657
<td> 2.03291
<td> 0.025876
<td> 2/77.2973
<td><a href="data/c48u65985.grid">grid</a>
<tr>
<td> c48u72521
<td align=right> 1740504
<td> 1.596
<td> 1.99529
<td> 0.024918
<td> 2<i>s</i>/33.2450
<td><a href="data/c48u72521.grid">grid</a>
<tr>
<td> c48u79099
<td align=right> 1898376
<td> 1.557
<td> 2.01914
<td> 0.024303
<td> 2/82.2973
<td><a href="data/c48u79099.grid">grid</a>
<tr>
<td> c48u86451
<td align=right> 2074824
<td> 1.505
<td> 1.99648
<td> 0.023504
<td> 2<i>s</i>/35.2450
<td><a href="data/c48u86451.grid">grid</a>
<tr>
<td> c48u93701
<td align=right> 2248824
<td> 1.467
<td> 2.00411
<td> 0.022911
<td> 2/87.2973
<td><a href="data/c48u93701.grid">grid</a>
<tr>
<td> c48u101477
<td align=right> 2435448
<td> 1.447
<td> 2.07920
<td> 0.022389
<td>
<td><a href="data/c48u101477.grid">grid</a>
<tr>
<td> c48u101917
<td align=right> 2446008
<td> 1.444
<td> 2.07768
<td> 0.022222
<td>
<td><a href="data/c48u101917.grid">grid</a>
<tr>
<td> c48u110143
<td align=right> 2643432
<td> 1.388
<td> 1.99316
<td> 0.021669
<td> 2/92.2973
<td><a href="data/c48u110143.grid">grid</a>
<tr>
<td> c48u118647
<td align=right> 2847528
<td> 1.358
<td> 2.01352
<td> 0.021210
<td> 2/94.2973
<td><a href="data/c48u118647.grid">grid</a>
<tr>
<td> c48u128249
<td align=right> 3077976
<td> 1.318
<td> 1.98655
<td> 0.020574
<td> 2<i>s</i>/40.2657
<td><a href="data/c48u128249.grid">grid</a>
<tr>
<td> c48u137809
<td align=right> 3307416
<td> 1.290
<td> 2.00301
<td> 0.020142
<td> 2/99.2973
<td><a href="data/c48u137809.grid">grid</a>
<tr>
<td> c48u148395
<td align=right> 3561480
<td> 1.255
<td> 1.98744
<td> 0.019600
<td> 2<i>s</i>/42.2657
<td><a href="data/c48u148395.grid">grid</a>
<tr>
<td> c48u158763
<td align=right> 3810312
<td> 1.228
<td> 1.99130
<td> 0.019176
<td> 2/104.2973
<td><a href="data/c48u158763.grid">grid</a>
<tr>
<td> c48u169757
<td align=right> 4074168
<td> 1.205
<td> 2.01122
<td> 0.018815
<td> 2/106.2973
<td><a href="data/c48u169757.grid">grid</a>
<tr>
<td> c48u181909
<td align=right> 4365816
<td> 1.173
<td> 1.98631
<td> 0.018310
<td> 2<i>s</i>/45.2450
<td><a href="data/c48u181909.grid">grid</a>
<tr>
<td> c48u193767
<td align=right> 4650408
<td> 1.151
<td> 2.00013
<td> 0.017970
<td> 2/111.2973
<td><a href="data/c48u193767.grid">grid</a>
<tr>
<td> c48u207023
<td align=right> 4968552
<td> 1.123
<td> 1.98553
<td> 0.017535
<td> 2<i>s</i>/47.2450
<td><a href="data/c48u207023.grid">grid</a>
<tr>
<td> c48u220121
<td align=right> 5282904
<td> 1.102
<td> 1.99143
<td> 0.017197
<td> 2/116.2973
<td><a href="data/c48u220121.grid">grid</a>
<tr>
<td> c48u233569
<td align=right> 5605656
<td> 1.083
<td> 2.00765
<td> 0.016906
<td> 2/118.2973
<td><a href="data/c48u233569.grid">grid</a>
<tr>
<td> c48u248571
<td align=right> 5965704
<td> 1.056
<td> 1.98203
<td> 0.016488
<td> 2/121.2973
<td><a href="data/c48u248571.grid">grid</a>
<tr>
<td> c48u263339
<td align=right> 6320136
<td> 1.039
<td> 1.99944
<td> 0.016221
<td> 2/123.2973
<td><a href="data/c48u263339.grid">grid</a>
<tr>
<td> c48u279565
<td align=right> 6709560
<td> 1.015
<td> 1.98032
<td> 0.015850
<td> 2<i>s</i>/52.2657
<td><a href="data/c48u279565.grid">grid</a>
<tr>
<td> c48u295333
<td align=right> 7087992
<td> 0.999
<td> 1.99038
<td> 0.015589
<td> 2/128.2973
<td><a href="data/c48u295333.grid">grid</a>
<tr>
<td> c48u312831
<td align=right> 7507944
<td> 0.978
<td> 1.97997
<td> 0.015266
<td> 2<i>s</i>/54.2657
<td><a href="data/c48u312831.grid">grid</a>
<tr>
<td> c48u330023
<td align=right> 7920552
<td> 0.961
<td> 1.98309
<td> 0.015004
<td> 2/133.2973
<td><a href="data/c48u330023.grid">grid</a>
<tr>
<td> c48u347617
<td align=right> 8342808
<td> 0.947
<td> 1.99747
<td> 0.014782
<td> 2/135.2973
<td><a href="data/c48u347617.grid">grid</a>
<tr>
<td> c48u367113
<td align=right> 8810712
<td> 0.927
<td> 1.97956
<td> 0.014472
<td> 2<i>s</i>/57.2450
<td><a href="data/c48u367113.grid">grid</a>
<tr>
<td> c48u386211
<td align=right> 9269064
<td> 0.913
<td> 1.99027
<td> 0.014255
<td> 2/140.2973
<td><a href="data/c48u386211.grid">grid</a>
<tr>
<td> c48u407099
<td align=right> 9770376
<td> 0.896
<td> 1.98011
<td> 0.013983
<td> 2<i>s</i>/59.2450
<td><a href="data/c48u407099.grid">grid</a>
<tr>
<td> c48u427333
<td align=right> 10255992
<td> 0.882
<td> 1.98284
<td> 0.013765
<td> 2/145.2973
<td><a href="data/c48u427333.grid">grid</a>
<tr>
<td> c48u448437
<td align=right> 10762488
<td> 0.870
<td> 1.99711
<td> 0.013578
<td> 2/147.2973
<td><a href="data/c48u448437.grid">grid</a>
<tr>
<td> c48u471503
<td align=right> 11316072
<td> 0.852
<td> 1.97662
<td> 0.013307
<td> 2/150.2973
<td><a href="data/c48u471503.grid">grid</a>
<tr>
<td> c48u493799
<td align=right> 11851176
<td> 0.841
<td> 1.98949
<td> 0.013132
<td> 2/152.2973
<td><a href="data/c48u493799.grid">grid</a>
<tr>
<td> c48u518377
<td align=right> 12441048
<td> 0.826
<td> 1.97564
<td> 0.012891
<td> 2<i>s</i>/64.2657
<td><a href="data/c48u518377.grid">grid</a>
<tr>
<td> c48u542361
<td align=right> 13016664
<td> 0.814
<td> 1.98354
<td> 0.012715
<td> 2/157.2973
<td><a href="data/c48u542361.grid">grid</a>
<tr>
<td> c48u566819
<td align=right> 13603656
<td> 0.814
<td> 2.06674
<td> 0.012551
<td>
<td><a href="data/c48u566819.grid">grid</a>
<tr>
<td> c48u568499
<td align=right> 13643976
<td> 0.807
<td> 2.02337
<td> 0.012499
<td>
<td><a href="data/c48u568499.grid">grid</a>
<tr>
<td> c48u593755
<td align=right> 14250120
<td> 0.789
<td> 1.97681
<td> 0.012323
<td> 2/162.2973
<td><a href="data/c48u593755.grid">grid</a>
<tr>
<td> c48u619981
<td align=right> 14879544
<td> 0.780
<td> 1.98967
<td> 0.012173
<td> 2/164.2973
<td><a href="data/c48u619981.grid">grid</a>
<tr>
<td> c48u648549
<td align=right> 15565176
<td> 0.766
<td> 1.97599
<td> 0.011964
<td> 2<i>s</i>/69.2450
<td><a href="data/c48u648549.grid">grid</a>
<tr>
<td> c48u676103
<td align=right> 16226472
<td> 0.757
<td> 1.98293
<td> 0.011813
<td> 2/169.2973
<td><a href="data/c48u676103.grid">grid</a>
<tr>
<td> c48u706351
<td align=right> 16952424
<td> 0.747
<td> 1.98930
<td> 0.011627
<td> 2<i>s</i>/71.2450
<td><a href="data/c48u706351.grid">grid</a>
<tr>
<td> c48u735777
<td align=right> 17658648
<td> 0.735
<td> 1.97798
<td> 0.011475
<td> 2/174.2973
<td><a href="data/c48u735777.grid">grid</a>
<tr>
<td> c48u765729
<td align=right> 18377496
<td> 0.727
<td> 1.98881
<td> 0.011344
<td> 2/176.2973
<td><a href="data/c48u765729.grid">grid</a>
<tr>
<td> c48u798587
<td align=right> 19166088
<td> 0.715
<td> 1.97265
<td> 0.011155
<td> 2/179.2973
<td><a href="data/c48u798587.grid">grid</a>
<tr>
<td> c48u830491
<td align=right> 19931784
<td> 0.707
<td> 1.98390
<td> 0.011032
<td> 2/181.2973
<td><a href="data/c48u830491.grid">grid</a>
<tr>
<td> c48u865149
<td align=right> 20763576
<td> 0.696
<td> 1.97317
<td> 0.010863
<td> 2<i>s</i>/76.2657
<td><a href="data/c48u865149.grid">grid</a>
<tr>
<td> c48u898517
<td align=right> 21564408
<td> 0.688
<td> 1.97768
<td> 0.010735
<td> 2/186.2973
<td><a href="data/c48u898517.grid">grid</a>
<tr>
<td> c48u932999
<td align=right> 22391976
<td> 0.683
<td> 2.01538
<td> 0.010620
<td>
<td><a href="data/c48u932999.grid">grid</a>
<tr>
<td> c48u970447
<td align=right> 23290728
<td> 0.670
<td> 1.97320
<td> 0.010455
<td> 2/191.2973
<td><a href="data/c48u970447.grid">grid</a>
<tr>
<td>c48u1006449
<td align=right> 24154776
<td> 0.663
<td> 1.98364
<td> 0.010347
<td> 2/193.2973
<td><a href="data/c48u1006449.grid">grid</a>
<tr>
<td>c48u1045817
<td align=right> 25099608
<td> 0.653
<td> 1.97289
<td> 0.010197
<td> 2<i>s</i>/81.2450
<td><a href="data/c48u1045817.grid">grid</a>
<tr>
<td>c48u1083955
<td align=right> 26014920
<td> 0.646
<td> 1.97878
<td> 0.010086
<td> 2/198.2973
<td><a href="data/c48u1083955.grid">grid</a>
<tr>
<td>
<td align=right> 27960696
<td> 0.630
<td> 1.97374
<td> 0.009838
<td> 2/203.2973
<td>
<tr>
<td>
<td align=right> 28943544
<td> 0.624
<td> 1.98389
<td> 0.009742
<td> 2/205.2973
<td>
</table>
</center>

<p>
As before, <a href="ExpandSet.cpp">ExpandSet.cpp</a> can be used to
expand the grids into sets of quaternions or Euler angles.  The
computation of the weights in these high density grid files is carried
out by determining the volume of the Voronoi cells in axis-angle space.
Some corrections are applied to account for the fact that orientation
space is not flat and the resulting maximum error in the weights is
approximately (&delta;/12)<sup>2</sup>.

<p>
The column labeled &ldquo;approx &delta;&rdquo; illustrates that the
optimal grid spacing follow one of three patterns (here <i>s</i> =
&radic;2 &minus; 1).   This follows from the requirement that the
separation of the lattice planes at one set of faces (triangles or
octogons) be such that the maximum Voronoi radius of points near this
face match that of the point at the center of the cell.  At the same
time, the separation of the lattice planes at the other set of faces
must not result in larger Voronoi radii.
