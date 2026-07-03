package com.github.micycle1.circupack;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Queue;
import java.util.stream.Collectors;

import com.github.micycle1.circupack.linalg.AMG;
import com.github.micycle1.circupack.linalg.BiCGStabSolver;
import com.github.micycle1.circupack.linalg.Preconditioner;
import com.github.micycle1.circupack.linalg.BiCGStabSolver.SparseCSR;
import com.github.micycle1.circupack.triangulation.Triangulation;

/**
 * <p>
 * Core repacking engine for circle packings — a Java port of Gerald Orick's
 * GOPacker MATLAB code. Derives combinatorics (interior component, boundary
 * loop, orphans) from a {@link Triangulation}, then iteratively alternates
 * boundary placement, Tutte-style interior embedding (sparse linear solve) and
 * effective-radii updates until the packing converges.
 * </p>
 *
 * <p>
 * Usage:
 * </p>
 * <ol>
 * <li>Construct with a {@code Triangulation} (0-based indices; see its javadoc
 * for the flower/boundary conventions).</li>
 * <li>Optionally configure the target shape: {@link #setPolygonPack(int[])},
 * {@link #setRectanglePack()}, etc. The default is {@link Mode#MAX_PACK}
 * (packing the unit disc).</li>
 * <li>Run {@link #riffle(double)} to iterate until the maximum relative visual
 * error falls below the given threshold.</li>
 * <li>Retrieve results via {@link #getRadii()}, {@link #getCentersX()},
 * {@link #getCentersY()}.</li>
 * </ol>
 *
 * <p>
 * The class is not thread-safe; callers should synchronize externally if
 * needed.
 * </p>
 */
public class CircuPacker {

	/**
	 * Packing modes, mirroring GOPack: MAX_PACK packs the unit disc (boundary
	 * circles become horocycles); POLYGONAL packs so the boundary forms an n-gon
	 * with prescribed corner angles (side lengths emerge from the packing).
	 */
	public enum Mode {
		MAX_PACK, POLYGONAL
	}

	private static final int MAX_ITER = 200;

	// External triangulation (combinatorics)
	private final Triangulation tri;

	// Geometry flag: -1 hyp, 0 eucl, +1 sph; we operate in Euclidean throughout
	// (only convert on save or additional logic if needed).
	private int hes = 0;

	// Preferred normalizing vertices
	private int alpha = -1;
	private int gamma = -1;

	// All vertices: 0..n-1
	private int n;

	// combinatoric cached info
	// interior component containing alpha
	private int[] intVerts = new int[0];
	private int intCount = 0;

	// boundary loop CCW (closed in our storage with last == first)
	private int[] bdryListClosed = new int[0];
	private int bdryCount = 0;

	// orphan verts (neither in alpha-component interiors nor boundary loop)
	private int[] orphanVerts = new int[0];
	private int orphanCount = 0;

	// mapping convenience
	private boolean[] isBoundary; // boundary marking
	private boolean hasBoundary; // whether triangulation has a boundary

	// Radii and centers (Euclidean), both working (local) and final (best-known)
	private double[] radii; // final radii (after riffle)
	private double[] centersX; // final centers
	private double[] centersY;

	private double[] localRadii; // working radii during iteration
	private double[] localCentersX; // working centers
	private double[] localCentersY;

	// aims (angle-sum targets); by default interior 2*pi, boundary negative
	// (MAX_PACK)
	private double[] vAims;

	// layout/rim sets and indexing for linear system
	private int[] layoutVerts; // interior vertices whose centers are unknowns
	private int[] rimVerts; // boundary/rim vertices with fixed centers
	private int[] v2indx; // maps original vertex -> index in (layout+rims) space (1..m or m+1..m+n_?); 0
							// means not used
	private int[] indx2v; // inverse mapping: index -> original vertex
	private int layCount; // layout vertex count
	private int rimCount; // rim vertex count (not including closure repetition)

	// incircle radii per interior vertex (sequential ordering)
	private double[][] inRadii; // size layCount x deg(v) for interior vertices only
	private double[] conduct; // total conductances per interior vertex (size layCount)

	// Matrices and RHS data needed for layout centers
	// We build A (m x m) sparse, b = RHS (m), where m=layCount
	// A has -1 diagonal and off-diagonal entries are scaled weights among interior
	// neighbors
	// RHS sums contributions from rim neighbors.
	// We'll build A and b on each layout step from current localRadii/localCenters.
	// We don't store explicit "transition/rhs" topologies separately,
	// since we can compute them deterministically from neighbors lists each time.
	private Mode mode = Mode.MAX_PACK;

	// For polygonal mode
	private int[] corners = null; // CCW corners
	private List<int[]> sides = null; // each side: list of vertices including endpoints; CCW order

	// Monitoring
	private final List<Double> visErrMonitor = new ArrayList<>();

	/**
	 * Builds the engine from the given triangulation: classifies
	 * interior/boundary/orphan vertices, seeds working geometry, and defaults to
	 * {@link Mode#MAX_PACK} aims. Configure a polygonal target afterwards via
	 * {@link #setPolygonPack(int[])} or {@link #setRectanglePack()} if desired,
	 * then call {@link #riffle(double)}.
	 */
	public CircuPacker(Triangulation tri) {
		this.tri = Objects.requireNonNull(tri, "tri must not be null");

		n = tri.getVertexCount();
		if (n <= 0) {
			throw new IllegalStateException("Triangulation has no vertices");
		}

		// boundary presence
		List<Integer> bdryLoop = tri.getBoundaryLoop();
		hasBoundary = bdryLoop != null && !bdryLoop.isEmpty();

		isBoundary = new boolean[n];
		for (int v = 0; v < n; v++) {
			isBoundary[v] = tri.isBoundaryVertex(v);
		}

		alpha = pickAnInteriorVertexOrFallback();

		// Derive intVerts, bdryListClosed, and orphanVerts, similar to complex_count.m
		// logic
		deriveInteriorBoundaryOrphans();

		// Seed geometry: uniform radii, centers at origin (both are recomputed by
		// riffle; radii serve as the initial guess)
		radii = new double[n];
		Arrays.fill(radii, 0.5);
		centersX = new double[n];
		centersY = new double[n];

		localRadii = radii.clone();
		localCentersX = centersX.clone();
		localCentersY = centersY.clone();

		vAims = new double[n];
		setMaxPack();

		// Default layoutVerts = intVerts; rimVerts = boundary loop (closed)
		layoutVerts = intVerts.clone();
		layCount = layoutVerts.length;
		rimVerts = new int[bdryCount];
		System.arraycopy(bdryListClosed, 0, rimVerts, 0, bdryCount);
		rimCount = rimVerts.length;

		// indexing arrays
		buildIndexing();
	}

	/**
	 * Configure classic max packing (the default): interior angle sums aim at
	 * 2&pi;, boundary radii adjust GO-style so the boundary hugs the unit circle.
	 */
	public void setMaxPack() {
		mode = Mode.MAX_PACK;
		corners = null;
		sides = null;
		Arrays.fill(vAims, 2.0 * Math.PI);
		for (int v = 0; v < n; v++) {
			if (isBoundary[v]) {
				vAims[v] = -1.0;
			}
		}
	}

	/**
	 * Configure polygonal packing with {@code numSides} corners chosen evenly
	 * spaced along the boundary loop, all corner angles equal.
	 */
	public void setPolygonPack(int numSides) {
		setPolygonPack(chooseSpacedCorners(numSides), null);
	}

	/**
	 * Configure polygonal packing with the given corner vertices (boundary
	 * vertices, in any order — they are sorted CCW along the boundary anchored at
	 * the first). Corner angles default to those of a regular n-gon,
	 * &pi;(1&minus;2/n).
	 */
	public void setPolygonPack(int[] cornersIn) {
		setPolygonPack(cornersIn, null);
	}

	/**
	 * Configure polygonal packing: the boundary is packed into an n-gon whose
	 * corners are the given boundary vertices with the given interior angles (side
	 * lengths emerge from the packing).
	 *
	 * @param cornersIn    at least 3 distinct boundary vertices; sorted CCW along
	 *                     the boundary anchored at the first supplied corner
	 * @param cornerAngles interior angle per corner in (0, &pi;), summing to
	 *                     (n&minus;2)&pi;; null for equal angles
	 */
	public void setPolygonPack(int[] cornersIn, double[] cornerAngles) {
		if (hes > 0 || !hasBoundary || bdryCount < 3) {
			throw new IllegalStateException("Polygon mode requires a boundary with at least 3 vertices.");
		}
		if (cornersIn == null || cornersIn.length < 3) {
			throw new IllegalArgumentException("Polygon mode requires at least 3 corners");
		}

		// Validate boundary membership + uniqueness
		boolean[] seen = new boolean[n];
		for (int c : cornersIn) {
			if (c < 0 || c >= n || !isBoundary[c]) {
				throw new IllegalArgumentException("Corner " + c + " is not a boundary vertex");
			}
			if (seen[c]) {
				throw new IllegalArgumentException("Duplicate corner vertex: " + c);
			}
			seen[c] = true;
		}

		// Validate corner angles if given: turning angles (pi - angle) sum to 2*pi
		if (cornerAngles != null) {
			if (cornerAngles.length != cornersIn.length) {
				throw new IllegalArgumentException("Corner angles length mismatch");
			}
			double sum = 0.0;
			for (double a : cornerAngles) {
				if (!(a > 0.0 && a < Math.PI)) {
					throw new IllegalArgumentException("Each corner angle must be in (0, pi): " + a);
				}
				sum += a;
			}
			double expected = cornersIn.length * Math.PI - 2.0 * Math.PI;
			if (Math.abs(sum - expected) > 1e-2) {
				throw new IllegalArgumentException("Corner angles inconsistent with polygon turning angles; expected sum " + expected + " but got " + sum);
			}
		}

		// Reorder corners into CCW boundary order, preserving the first supplied
		// corner as anchor
		OrderedCorners ordered = orderCornersCCW(cornersIn.clone(), cornerAngles);
		this.corners = ordered.corners;
		this.sides = derivePolygonSidesFromCorners(this.corners);
		this.mode = Mode.POLYGONAL;

		// Aims: interior 2*pi; boundary pi; corners get their angles
		Arrays.fill(vAims, 2.0 * Math.PI);
		for (int w : rimVerts) {
			vAims[w] = Math.PI;
		}
		if (ordered.angles != null) {
			for (int i = 0; i < this.corners.length; i++) {
				vAims[this.corners[i]] = ordered.angles[i];
			}
		} else {
			int m = this.corners.length;
			double equalCorner = Math.PI * (1.0 - 2.0 / m);
			for (int c : this.corners) {
				vAims[c] = equalCorner;
			}
		}
	}

	/**
	 * Configure rectangular packing with 4 corners chosen evenly spaced along the
	 * boundary loop. See {@link #setRectanglePack(int[])}.
	 */
	public void setRectanglePack() {
		setRectanglePack(chooseSpacedCorners(4));
	}

	/**
	 * Configure rectangular packing: polygonal packing with 4 corners and all
	 * corner angles &pi;/2. The aspect ratio of the resulting rectangle emerges
	 * from the packing; query it with {@link #getAspect()}.
	 */
	public void setRectanglePack(int[] fourCorners) {
		if (fourCorners == null || fourCorners.length != 4) {
			throw new IllegalArgumentException("Rectangle packing requires exactly 4 corners");
		}
		double[] rightAngles = new double[4];
		Arrays.fill(rightAngles, Math.PI / 2.0);
		setPolygonPack(fourCorners, rightAngles);
	}

	/** The current packing mode. */
	public Mode getMode() {
		return mode;
	}

	/**
	 * The CCW-ordered corner vertices in use for polygonal packing, or null in
	 * MAX_PACK mode. Useful when corners were auto-chosen.
	 */
	public int[] getCorners() {
		return corners == null ? null : corners.clone();
	}

	/**
	 * Aspect ratio (top+bottom)/(left+right) of a rectangular packing, computed
	 * from the current centers. Requires POLYGONAL mode with exactly 4 corners;
	 * call after {@link #riffle(double)}.
	 */
	public double getAspect() {
		if (mode != Mode.POLYGONAL || corners == null || corners.length != 4) {
			throw new IllegalStateException("Aspect requires polygonal mode with 4 corners");
		}
		double top = hyp(corners[1], corners[0]);
		double rend = hyp(corners[3], corners[0]);
		double lend = hyp(corners[2], corners[1]);
		double bot = hyp(corners[3], corners[2]);
		return (top + bot) / (rend + lend);
	}

	/**
	 * Iterate boundary layout, interior embedding, and effective-radii updates
	 * until the maximum relative visual error falls below the threshold (or an
	 * iteration cap is hit). Returns the number of passes run.
	 *
	 * @param maxRelativeError in practice a lot visually lower!
	 */
	public int riffle(double maxRelativeError) {
		maxRelativeError = Math.max(1e-4, maxRelativeError);
		int pass = 0;
		double maxVis = Double.MAX_VALUE;

		while (maxVis > maxRelativeError && pass < MAX_ITER) {
			layoutBoundary();
			layoutCentersSolve();
			setEffectiveRadii();

			maxVis = updateVisErrorMonitor();
			pass++;
		}

		// final consistent geometry
		layoutBoundary();
		layoutCentersSolve();

		radii = localRadii.clone();
		centersX = localCentersX.clone();
		centersY = localCentersY.clone();
		return pass;
	}

	/****
	 * <p>
	 * Computes, for each interior vertex, the maximum "relative visual error" (RVE)
	 * across all incident edges and returns those maxima in an array.
	 * </p>
	 *
	 * <p>
	 * For a given interior vertex <code>v = layoutVerts[i]</code>, each neighbor
	 * <code>w</code> in v's flower is inspected and the edge RVE is computed as:
	 * </p>
	 *
	 * <pre>
	 * RVE(v,w) = | |z_v - z_w| - (r_v + r_w) | / r_v
	 * </pre>
	 *
	 * <ul>
	 * <li><code>z_v, z_w</code> are Euclidean centers
	 * (<code>localCentersX</code>/<code>localCentersY</code>).</li>
	 * <li><code>r_v, r_w</code> are Euclidean radii (<code>localRadii</code>).</li>
	 * <li><code>|z_v - z_w|</code> is the center-to-center distance (computed with
	 * <code>Math.hypot</code>).</li>
	 * </ul>
	 *
	 * <p>
	 * The per-vertex value is the maximum <code>RVE(v,w)</code> over all neighbors
	 * <code>w</code> of <code>v</code>. The returned array contains those
	 * per-vertex maxima in the same order as <code>layoutVerts</code>.
	 * </p>
	 *
	 * <h4>Units</h4>
	 * <p>
	 * Dimensionless (ratio). Interpretation examples:
	 * </p>
	 * <ul>
	 * <li><code>0.00</code> = perfect tangency on all incident edges (ideal).</li>
	 * <li><code>0.01</code> ≈ 1% of v's radius (very good).</li>
	 * <li><code>0.10</code> ≈ 10% of v's radius (noticeable).</li>
	 * <li><code>1.00</code> ≈ mismatch comparable to v's radius (poor).</li>
	 * </ul>
	 *
	 * <h4>Notes</h4>
	 * <ul>
	 * <li>Only vertices listed in <code>layoutVerts</code> (interior/layout
	 * vertices) are processed.</li>
	 * <li>Boundary-to-boundary edges are not directly included unless the boundary
	 * vertex is a neighbor of an interior vertex.</li>
	 * <li>To avoid division by zero when <code>r_v</code> is extremely small, the
	 * denominator uses <code>Math.max(1e-16, r_v)</code>. Very small radii can
	 * cause the ratio to spike.</li>
	 * </ul>
	 *
	 * @return an array of length <code>layCount</code> where element <code>i</code>
	 *         is the maximum relative visual error for vertex
	 *         <code>layoutVerts[i]</code> (dimensionless).
	 */
	public double[] visualErrors() {
		double[] errs = new double[layCount];
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double zvx = localCentersX[v], zvy = localCentersY[v];
			double rv = localRadii[v];
			List<Integer> flower = tri.getFlower(v);
			double maxErr = 0.0;
			for (int w : flower) {
				double dx = zvx - localCentersX[w];
				double dy = zvy - localCentersY[w];
				double cdiff = Math.sqrt(dx * dx + dy * dy);
				double rdiff = rv + localRadii[w];
				double me = Math.abs(cdiff - rdiff) / Math.max(1e-16, rv);
				if (me > maxErr) {
					maxErr = me;
				}
			}
			errs[i] = maxErr;
		}
		return errs;
	}

	public double[] angleSumErrors() {
		// diffs(k) = anglesum(v)-2*pi for interior layout vertices
		double target = -2.0 * Math.PI;
		double[] diffs = new double[layCount];
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double r = localRadii[v];
			List<Integer> flower = tri.getFlower(v);
			double diff = target;
			int m = flower.size();
			for (int j = 0; j < m; j++) {
				int w = flower.get(j);
				int u = flower.get((j + 1) % m);
				double cosang = MathUtil.cosAngle(r, localRadii[w], localRadii[u]);
				diff += Math.acos(Math.max(-1.0, Math.min(1.0, cosang)));
			}
			diffs[i] = diff;
		}
		return diffs;
	}

	/**
	 * <p>
	 * Remove orphan vertices (those not in the interior component and not on the
	 * rim).
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>Constructs a trimmed mapping that preserves the alpha-component interiors
	 * and the boundary, discarding orphan vertices and updating radii/centers
	 * arrays accordingly.</li>
	 * <li>Recomputes combinatorics and index mappings after removal.</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * Orphans complicate layout and may produce artifacts; pruned complexes are
	 * smaller and lead to simpler linear systems. This is an optional
	 * transformation that mutates internal arrays to operate on the reduced
	 * complex.
	 * </p>
	 */
	public void pruneOrphans() {
		// If no orphans, nothing to do
		if (orphanCount == 0) {
			return;
		}

		// Remove orphan vertices from triangulation view within the engine
		// We don't mutate Triangulation; we adjust our working sets and geometry
		// arrays.
		int newN = intCount + bdryCount; // keep alpha-component interiors and boundary
		boolean[] keep = new boolean[n];
		for (int v : intVerts) {
			keep[v] = true;
		}
		for (int i = 0; i < bdryCount; i++) {
			keep[bdryListClosed[i]] = true;
		}

		int[] old2new = new int[n];
		Arrays.fill(old2new, -1);
		int kk = 0;
		for (int v = 0; v < n; v++) {
			if (keep[v]) {
				old2new[v] = kk++;
			}
		}

		// Rebuild arrays
		double[] nr = new double[newN];
		double[] nx = new double[newN];
		double[] ny = new double[newN];
		for (int v = 0; v < n; v++) {
			if (keep[v]) {
				int nv = old2new[v];
				nr[nv] = localRadii[v];
				nx[nv] = localCentersX[v];
				ny[nv] = localCentersY[v];
			}
		}
		localRadii = nr;
		localCentersX = nx;
		localCentersY = ny;
		radii = nr.clone();
		centersX = nx.clone();
		centersY = ny.clone();

		// Rebuild interior/boundary indices
		int[] nInt = new int[intCount];
		for (int i = 0; i < intCount; i++) {
			nInt[i] = old2new[intVerts[i]];
		}
		int[] nBdryClosed = new int[bdryCount + 1];
		for (int i = 0; i < bdryCount + 1; i++) {
			nBdryClosed[i] = old2new[bdryListClosed[i]];
		}

		intVerts = nInt;
		intCount = nInt.length;
		bdryListClosed = nBdryClosed;
		bdryCount = nBdryClosed.length - 1;
		orphanVerts = new int[0];
		orphanCount = 0;

		// Layout sets
		layoutVerts = intVerts.clone();
		layCount = layoutVerts.length;
		rimVerts = new int[bdryCount];
		System.arraycopy(bdryListClosed, 0, rimVerts, 0, bdryCount);
		rimCount = rimVerts.length;

		buildIndexing(); // rebuild v2indx/indx2v
	}

	/** Final radii (valid after {@link #riffle(double)}). */
	public double[] getRadii() {
		return radii;
	}

	/** Final center x-coordinates (valid after {@link #riffle(double)}). */
	public double[] getCentersX() {
		return centersX;
	}

	/** Final center y-coordinates (valid after {@link #riffle(double)}). */
	public double[] getCentersY() {
		return centersY;
	}

	/** Max relative visual error recorded after each riffle pass. */
	public List<Double> getVisErrMonitor() {
		return Collections.unmodifiableList(visErrMonitor);
	}

	// Getters for testing
	public double[] getLocalRadii() {
		return localRadii;
	}

	public double[] getLocalCentersX() {
		return localCentersX;
	}

	public double[] getLocalCentersY() {
		return localCentersY;
	}

	public int[] getLayoutVerts() {
		return layoutVerts;
	}

	public int[] getRimVerts() {
		return rimVerts;
	}

	public int[] getIntVerts() {
		return intVerts;
	}

	public int[] getBdryLoopClosed() {
		return bdryListClosed;
	}

	private int[] chooseSpacedCorners(int sideN) {
		if (bdryCount < 3) {
			throw new IllegalStateException("Boundary must have at least 3 vertices");
		}
		sideN = Math.max(3, Math.min(sideN, bdryCount));

		int[] crn = new int[sideN];
		int step = Math.max(1, bdryCount / sideN);

		int start = 0;
		if (gamma >= 0) {
			for (int i = 0; i < bdryCount; i++) {
				if (bdryListClosed[i] == gamma) {
					start = i;
					break;
				}
			}
		}

		for (int e = 0; e < sideN; e++) {
			crn[e] = bdryListClosed[(start + e * step) % bdryCount];
		}
		return crn;
	}

	private int pickAnInteriorVertexOrFallback() {
		// Try provided interior vertices list
		List<Integer> ints = tri.getInteriorVertices();
		if (ints != null && !ints.isEmpty()) {
			return ints.get(0);
		}
		// Else pick any vertex whose flower is cyclic (interior)
		for (int v = 0; v < n; v++) {
			if (!isBoundary[v]) {
				return v;
			}
		}
		// fallback: just 0 if nothing else (we'll treat as sphere)
		return 0;
	}

	private void deriveInteriorBoundaryOrphans() {
		List<Integer> bdryLoop = tri.getBoundaryLoop();
		if (bdryLoop == null || bdryLoop.isEmpty()) {
			// Sphere case: no boundary
			hasBoundary = false;
			hes = 1; // spherical
			// Choose a "faux" boundary of 3 vertices: take far vertex a from alpha (BFS
			// depth),
			// then its first two neighbors as b and c
			int a = farVertex(Collections.singletonList(alpha));
			List<Integer> fl = tri.getFlower(a);
			int b = fl.get(1 % fl.size());
			int c = fl.get(0);
			bdryListClosed = new int[] { a, b, c, a };
			bdryCount = 3;

			// All others are interiors
			List<Integer> ints = new ArrayList<>();
			for (int v = 0; v < n; v++) {
				if (v != a && v != b && v != c) {
					ints.add(v);
				}
			}
			intVerts = ints.stream().mapToInt(Integer::intValue).toArray();
			intCount = intVerts.length;

			// orphans empty
			orphanVerts = new int[0];
			orphanCount = 0;
			gamma = a;
			return;
		}

		hasBoundary = true;
		hes = 0;

		// Use BFS to find interior component containing alpha
		boolean[] touchedInterior = new boolean[n];
		Queue<Integer> q = new ArrayDeque<>();
		if (!isBoundary[alpha]) {
			q.add(alpha);
			touchedInterior[alpha] = true;
		}

		while (!q.isEmpty()) {
			int v = q.poll();
			for (int w : tri.getFlower(v)) {
				if (!isBoundary[w] && !touchedInterior[w]) {
					touchedInterior[w] = true;
					q.add(w);
				}
			}
		}

		// Collect intVerts (alpha component)
		List<Integer> ints = new ArrayList<>();
		for (int v = 0; v < n; v++) {
			if (!isBoundary[v] && touchedInterior[v]) {
				ints.add(v);
			}
		}
		intVerts = ints.stream().mapToInt(Integer::intValue).toArray();
		intCount = intVerts.length;

		// Build bdry loop closed (from Triangulation)
		bdryCount = tri.getBoundaryLoop().size();
		bdryListClosed = new int[bdryCount + 1];
		for (int i = 0; i < bdryCount; i++) {
			bdryListClosed[i] = tri.getBoundaryLoop().get(i);
		}
		bdryListClosed[bdryCount] = bdryListClosed[0];

		// Orphans: not in intVerts and not in boundary loop
		boolean[] isInInt = new boolean[n];
		for (int v : intVerts) {
			isInInt[v] = true;
		}
		boolean[] isInB = new boolean[n];
		for (int i = 0; i < bdryCount; i++) {
			isInB[bdryListClosed[i]] = true;
		}

		List<Integer> orph = new ArrayList<>();
		for (int v = 0; v < n; v++) {
			if (!isInInt[v] && !isInB[v]) {
				orph.add(v);
			}
		}
		orphanVerts = orph.stream().mapToInt(Integer::intValue).toArray();
		orphanCount = orphanVerts.length;

		// pick gamma if not provided
		if (gamma < 0 || gamma >= n || !isBoundary[gamma]) {
			gamma = bdryListClosed[0];
		}
	}

	private int farVertex(List<Integer> seeds) {
		if (seeds == null || seeds.isEmpty()) {
			return (n > 0 ? 0 : -1);
		}
		int[] marks = new int[n]; // 0 = unvisited
		Queue<Integer> q = new ArrayDeque<>();
		for (int s : seeds) {
			q.add(s);
			marks[s] = 1;
		}
		int far = seeds.get(0);
		int depth = 1;
		while (!q.isEmpty() && depth < 10) {
			int sz = q.size();
			for (int i = 0; i < sz; i++) {
				int v = q.poll();
				far = v;
				for (int w : tri.getFlower(v)) {
					if (marks[w] == 0) {
						marks[w] = marks[v] + 1;
						q.add(w);
					}
				}
			}
			depth++;
		}
		return far;
	}

	/**
	 * <p>
	 * Build index mappings used for matrix assembly and neighbor lookups.
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>Creates {@code indx2v} (sequential indices: layout verts then rim
	 * verts).</li>
	 * <li>Creates {@code v2indx} mapping original vertex -> index in the combined
	 * list.</li>
	 * <li>Sets {@code layCount} and {@code rimCount} used by solver and update
	 * routines.</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * Sparse matrix assembly and RHS construction require a compact contiguous
	 * indexing; this method centralizes that translation so the rest of the code
	 * works with small dense/sparse structures indexed 0..m-1.
	 * </p>
	 */
	private void buildIndexing() {
		// indx2v: layout first, then rim vertices
		indx2v = new int[layCount + rimCount];
		v2indx = new int[n];
		Arrays.fill(v2indx, -1);
		for (int i = 0; i < layCount; i++) {
			indx2v[i] = layoutVerts[i];
			v2indx[layoutVerts[i]] = i;
		}
		for (int j = 0; j < rimCount; j++) {
			indx2v[layCount + j] = rimVerts[j];
			v2indx[rimVerts[j]] = layCount + j;
		}
	}

	/**
	 * <p>
	 * Place boundary centers according to the current packing mode.
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>Dispatches to {@code setHoroCenters}, {@code setPolyCenters}, or
	 * {@code setRectCenters} depending on the engine {@code mode} and polygon
	 * metadata.</li>
	 * <li>May also scale {@code localRadii} to match a normalized boundary (e.g.,
	 * unit circle).</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * Boundary placement is mode-specific: max-pack uses horocycles in the unit
	 * disc, polygonal mode places centers on a normalized n-gon. Centralizing the
	 * choice avoids duplicated logic and ensures consistent scaling before interior
	 * layout.
	 * </p>
	 */
	private void layoutBoundary() {
		if (mode == Mode.POLYGONAL) {
			setPolygonCenters();
		} else {
			setHoroCenters();
		}
	}

	/**
	 * <p>
	 * Lay out boundary circles as horocycles inside the unit circle for MAX_PACK.
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>If only three boundary vertices: place them symmetrically and set aims to
	 * zero.</li>
	 * <li>Otherwise: compute a scalar {@code R} by Newton iteration so that the
	 * angles formed by tangent circles sum to 2π, then scale radii by 1/R.</li>
	 * <li>Position boundary centers consecutively along the unit circle using the
	 * computed turning angles; the first center is placed on the positive imaginary
	 * axis.</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * This is the standard initialization for disc-style max packing: boundary
	 * circles are placed so that when interpreted as horocycles their tangency
	 * geometry matches the intended discrete turning angles, providing a stable
	 * starting layout for interior Tutte-style embedding.
	 * </p>
	 */
	private void setHoroCenters() {
		// Special case: only 3 boundary vertices => equilateral symmetric
		if (bdryCount <= 3) {
			double s3 = Math.sqrt(3);
			double brad = s3 / (2.0 + s3);
			int v1 = bdryListClosed[0];
			int v2 = bdryListClosed[1];
			int v3 = bdryListClosed[2];
			localRadii[v1] = brad;
			localRadii[v2] = brad;
			localRadii[v3] = brad;
			vAims[v1] = 0;
			vAims[v2] = 0;
			vAims[v3] = 0;
			// centers equally spaced on circle radius (1-brad)
			double R = 1.0 - brad;
			localCentersX[v1] = 0.0;
			localCentersY[v1] = +R;
			localCentersX[v2] = -R * Math.sqrt(3) / 2.0;
			localCentersY[v2] = -R * 0.5;
			localCentersX[v3] = +R * Math.sqrt(3) / 2.0;
			localCentersY[v3] = -R * 0.5;
			return;
		}

		// Newton iteration to solve for circle around unit circle
		// r[j] = localRadii[bdryList[j]]
		double[] r = new double[bdryCount + 1];
		double sum = 0.0, minrad = 0.0;
		for (int i = 0; i < bdryCount; i++) {
			r[i] = localRadii[bdryListClosed[i]];
			minrad = Math.max(minrad, r[i]);
			sum += r[i];
		}
		r[bdryCount] = r[0];
		double R = sum / Math.PI;
		if (R < 2.0 * minrad) {
			R = 3.0 * minrad;
		}

		int trys = 0;
		while (trys < 100) {
			trys++;
			double f = -2.0 * Math.PI;
			double fp = 0.0;
			for (int j = 0; j < bdryCount; j++) {
				double r1 = r[j], r2 = r[j + 1];
				double Rrr = R - r1 - r2;
				double RRrr = R * Rrr;
				double ab = r1 * r2;
				double val = (RRrr - ab) / (RRrr + ab);
				val = Math.max(-1.0, Math.min(1.0, val));
				f += Math.acos(val);
				// derivative approx (as in MATLAB)
				double denom = (RRrr + ab);
				double num = -1.0 * (R + Rrr) * Math.sqrt(Math.max(1e-16, ab / RRrr));
				fp += num / denom;
			}
			double newR = R - f / fp;
			newR = Math.max(R / 2.0, Math.min(2.0 * R, newR));
			if (Math.abs(newR - R) < 1e-5) {
				R = newR;
				break;
			}
			R = newR;
		}

		// Scale all radii by 1/R
		for (int v = 0; v < n; v++) {
			localRadii[v] /= R;
		}
		for (int j = 0; j <= bdryCount; j++) {
			r[j] /= R;
		}

		// Place boundary centers around unit circle, first on +i axis
		int v1 = bdryListClosed[0];
		localCentersX[v1] = 0.0;
		localCentersY[v1] = 1.0 - r[0];
		double arg = Math.PI / 2.0;
		double rPrev = r[0];
		for (int k = 1; k < bdryCount; k++) {
			double r2 = r[k];
			double RRrr = 1.0 - rPrev - r2;
			double ab = rPrev * r2;
			double val = (RRrr - ab) / (RRrr + ab);
			val = Math.max(-1.0, Math.min(1.0, val));
			double delta = Math.acos(val);
			arg += delta;
			double d = 1.0 - r2;
			int vk = bdryListClosed[k];
			localCentersX[vk] = d * Math.cos(arg);
			localCentersY[vk] = d * Math.sin(arg);
			rPrev = r2;
		}
	}

	/**
	 * <p>
	 * Place boundary centers and scale radii for polygonal packing (general n-gon).
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>Computes current side lengths from {@code localRadii}, derives target
	 * side lengths (triangles by law of sines, even-odd strategies otherwise), and
	 * scales radii.</li>
	 * <li>Computes edge directions from corner aims and lays out each side
	 * incrementally, then recenters and rescales so a normalization condition holds
	 * (first corner at 1+i or i).</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * Polygonal packing requires boundary centers to sit on a specific polygonal
	 * shape whose side lengths are consistent with circle diameters; this method
	 * enforces those proportions and normalization so interior embedding solves
	 * behave predictably.
	 * </p>
	 */
	private void setPolygonCenters() {
		if (corners == null || sides == null) {
			throw new IllegalStateException("Polygonal mode without corners; configure via setPolygonPack/setRectanglePack");
		}

		// if rectangle (4 corners) with all right angles, do rectangular layout
		if (corners.length == 4) {
			double bendErr = 0.0;
			for (int c : corners) {
				bendErr += Math.abs(vAims[c] - Math.PI / 2.0);
			}
			if (bendErr <= 1e-5) {
				setRectCenters();
				return;
			}
		}

		// general polygon: set side lengths from radii; compute target side lengths
		int numSides = corners.length;

		double[] sideLengths = new double[numSides];
		double fullLength = 0.0;
		double[] targetLength = new double[numSides];

		for (int i = 0; i < numSides; i++) {
			int[] side = sides.get(i);
			double L = localRadii[side[0]] + localRadii[side[side.length - 1]];
			for (int j = 1; j < side.length - 1; j++) {
				L += 2.0 * localRadii[side[j]];
			}
			sideLengths[i] = L;
			fullLength += L;
		}

		int halfn = numSides / 2;

		if (numSides == 3) {
			// triangle with corner aims provided in vAims at corners
			double opp1 = vAims[corners[2]];
			double opp2 = vAims[corners[0]];
			if (opp1 <= 0 || opp2 <= 0 || (opp1 + opp2) >= Math.PI) {
				throw new IllegalStateException("Invalid triangle corner aims: " + opp2 + ", " + vAims[corners[1]] + ", " + opp1);
			}
			double opp3 = Math.PI - (opp1 + opp2);
			targetLength[0] = 1.0;
			targetLength[1] = Math.sin(opp2) / Math.sin(opp1);
			targetLength[2] = Math.sin(opp3) / Math.sin(opp1);
			double sum = targetLength[0] + targetLength[1] + targetLength[2];
			double factor = sum / fullLength;
			scaleAllLocalBy(factor);
			for (int i = 0; i < numSides; i++) {
				sideLengths[i] *= factor;
			}
		} else if (2 * halfn == numSides) {
			// even: pair opposite sides, target total per pair equals average; set scale so
			// sum target ~ 6
			double factor = 6.0 / fullLength;
			scaleAllLocalBy(factor);
			for (int i = 0; i < numSides; i++) {
				sideLengths[i] *= factor;
			}
			for (int j = 0; j < halfn; j++) {
				targetLength[j] = (sideLengths[j] + sideLengths[halfn + j]) / 2.0;
				targetLength[halfn + j] = targetLength[j];
			}
		} else {
			// odd: make sides equal to 2*sin(pi/n)
			double spn = 2.0 * Math.sin(Math.PI / numSides);
			double factor = numSides * spn / fullLength;
			scaleAllLocalBy(factor);
			for (int i = 0; i < numSides; i++) {
				sideLengths[i] *= factor;
				targetLength[i] = spn;
			}
		}

		// Layout rules:
		// - Odd: first corner at (0,1), first edge direction is pi + (pi - aim1)/2
		// - Even: first corner at (1,1), first edge horizontal left (pi)
		double[] edgeArg = new double[numSides];
		if (2 * halfn != numSides) {
			edgeArg[0] = Math.PI + (Math.PI - vAims[corners[0]]) / 2.0;
		} else {
			edgeArg[0] = Math.PI;
		}
		for (int j = 1; j < numSides; j++) {
			edgeArg[j] = edgeArg[j - 1] + Math.PI - vAims[corners[j]];
		}

		double[] dirX = new double[numSides];
		double[] dirY = new double[numSides];
		for (int i = 0; i < numSides; i++) {
			dirX[i] = Math.cos(edgeArg[i]);
			dirY[i] = Math.sin(edgeArg[i]);
		}

		if (2 * halfn == numSides) {
			localCentersX[corners[0]] = 1.0;
			localCentersY[corners[0]] = 1.0;
		} else {
			localCentersX[corners[0]] = 0.0;
			localCentersY[corners[0]] = 1.0;
		}

		for (int k = 0; k < numSides; k++) {
			int[] side = sides.get(k);
			double factor = targetLength[k] / Math.max(1e-16, sideLengths[k]);
			int start = corners[k];
			localCentersX[side[0]] = localCentersX[start];
			localCentersY[side[0]] = localCentersY[start];

			double prev = localRadii[start];
			double spotX = localCentersX[start];
			double spotY = localCentersY[start];

			for (int i = 0; i < side.length - 1; i++) {
				int nxt = side[i + 1];
				double next = localRadii[nxt];
				spotX += factor * dirX[k] * (prev + next);
				spotY += factor * dirY[k] * (prev + next);
				localCentersX[nxt] = spotX;
				localCentersY[nxt] = spotY;
				prev = next;
			}
		}

		// center polygon: translate by centroid of corners
		double avgx = 0, avgy = 0;
		for (int c : corners) {
			avgx += localCentersX[c];
			avgy += localCentersY[c];
		}
		avgx /= corners.length;
		avgy /= corners.length;
		for (int v = 0; v < n; v++) {
			localCentersX[v] -= avgx;
			localCentersY[v] -= avgy;
		}

		// scale: even => make corner[0].x = 1; odd => make corner[0].y = 1
		double scal = (2 * halfn == numSides) ? Math.max(1e-16, localCentersX[corners[0]]) : Math.max(1e-16, localCentersY[corners[0]]);
		for (int v = 0; v < n; v++) {
			localCentersX[v] /= scal;
			localCentersY[v] /= scal;
			localRadii[v] /= scal;
		}
	}

	/**
	 * <p>
	 * Specialized layout for rectangular packings (four sides, right-angled
	 * corners).
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>Computes averaged top/bottom and left/right side lengths and an aspect
	 * ratio.</li>
	 * <li>Scales {@code localRadii} to fit a rectangle centered at origin with y =
	 * ±1 and x = ±aspect.</li>
	 * <li>Distributes boundary centers proportionally along each rectangle
	 * edge.</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * Rectangles are a common practical target; handling them explicitly yields
	 * better aspect preservation and more stable side-length matching than the
	 * general polygon routine.
	 * </p>
	 */
	private void setRectCenters() {
		// sides.size() == 4; corners count == 4
		double[] sideLen = new double[4];
		for (int i = 0; i < 4; i++) {
			int[] side = sides.get(i);
			double L = localRadii[side[0]] + localRadii[side[side.length - 1]];
			for (int j = 1; j < side.length - 1; j++) {
				L += 2.0 * localRadii[side[j]];
			}
			sideLen[i] = L;
		}
		double width = (sideLen[0] + sideLen[2]) / 2.0;
		double height = (sideLen[1] + sideLen[3]) / 2.0;
		double aspect = height / Math.max(1e-16, width);

		// scale radii so that total fits rectangle of width 2*aspect and height 2
		double factor = 2.0 * (aspect + 1.0) / Math.max(1e-16, (width + height));
		scaleAllLocalBy(factor);
		for (int i = 0; i < 4; i++) {
			sideLen[i] *= factor;
		}

		// rectangle corners: (1,aspect), (-1,aspect), (-1,-aspect), (1,-aspect)
		double[] crx = { 1.0, -1.0, -1.0, 1.0 };
		double[] cry = { aspect, aspect, -aspect, -aspect };
		double[] dirx = { -1.0, 0.0, 1.0, 0.0 };
		double[] diry = { 0.0, -1.0, 0.0, 1.0 };
		double[] slength = { 2.0, 2.0 * aspect, 2.0, 2.0 * aspect };

		for (int k = 0; k < 4; k++) {
			int[] side = sides.get(k);
			double sf = slength[k] / Math.max(1e-16, sideLen[k]);

			int start = corners[k];
			double spotX = crx[k];
			double spotY = cry[k];
			localCentersX[side[0]] = spotX;
			localCentersY[side[0]] = spotY;

			// place only the side's internal vertices: each corner stays exactly at
			// its rectangle corner (set when its own side starts)
			double prev = localRadii[start];
			for (int i = 0; i < side.length - 2; i++) {
				int nxt = side[i + 1];
				double next = localRadii[nxt];
				spotX += sf * dirx[k] * (prev + next);
				spotY += sf * diry[k] * (prev + next);
				localCentersX[nxt] = spotX;
				localCentersY[nxt] = spotY;
				prev = next;
			}
		}
	}

	private void scaleAllLocalBy(double factor) {
		for (int v = 0; v < n; v++) {
			localRadii[v] *= factor;
		}
	}

	private List<int[]> derivePolygonSidesFromCorners(int[] orderedCorners) {
		if (orderedCorners == null || orderedCorners.length < 3) {
			throw new IllegalArgumentException("Need at least 3 corners");
		}

		List<Integer> loop = Arrays.stream(bdryListClosed, 0, bdryCount).boxed().collect(Collectors.toList());
		List<int[]> result = new ArrayList<>(orderedCorners.length);

		for (int i = 0; i < orderedCorners.length; i++) {
			int a = orderedCorners[i];
			int b = orderedCorners[(i + 1) % orderedCorners.length];

			int ia = boundaryPosition(a);
			if (ia < 0) {
				throw new IllegalStateException("Corner " + a + " not found in boundary loop");
			}

			List<Integer> side = new ArrayList<>();
			int k = ia;
			side.add(loop.get(k));

			while (loop.get(k) != b) {
				k = (k + 1) % bdryCount;
				side.add(loop.get(k));
				if (side.size() > bdryCount + 1) {
					throw new IllegalStateException("Side extraction overflow");
				}
			}

			result.add(side.stream().mapToInt(Integer::intValue).toArray());
		}

		return result;
	}

	private void layoutCentersSolve() {
		updateVdata();

		int n = layCount;

		// First pass: count nnz per row (we always add diagonal)
		int[] rowCounts = new int[n];
		Arrays.fill(rowCounts, 1); // diagonal per row
		double[] bx = new double[n];
		double[] by = new double[n];

		for (int i = 0; i < n; i++) {
			int v = layoutVerts[i];
			double vrad = localRadii[v];
			var fl = tri.getFlower(v);
			int m = fl.size();
			if (m == 0) {
				continue;
			}

			double[] iR = inRadii[i];
			double invTot = 1.0 / Math.max(1e-16, conduct[i]);

			double prevR = iR[m - 1];
			for (int j = 0; j < m; j++) {
				int w = fl.get(j);
				double coeff = (prevR + iR[j]) * invTot / Math.max(1e-16, vrad + localRadii[w]);
				prevR = iR[j];

				int idxW = v2indx[w];
				if (0 <= idxW && idxW < n) {
					rowCounts[i]++; // internal neighbor contributes to A
				} else {
					// external neighbor contributes to RHS
					bx[i] -= coeff * localCentersX[w];
					by[i] -= coeff * localCentersY[w];
				}
			}
		}

		// Build CSR structure
		int[] rowPtr = new int[n + 1];
		rowPtr[0] = 0;
		for (int i = 0; i < n; i++) {
			rowPtr[i + 1] = rowPtr[i] + rowCounts[i];
		}
		int nnz = rowPtr[n];
		int[] colIdx = new int[nnz];
		double[] val = new double[nnz];
		double[] Minv = new double[n];

		// Second pass: fill values
		int[] next = new int[n];
		System.arraycopy(rowPtr, 0, next, 0, n);

		for (int i = 0; i < n; i++) {
			// Put diagonal first
			int pos = next[i]++;
			colIdx[pos] = i;
			val[pos] = -1.0; // your system has diagonal -1
			Minv[i] = -1.0 != 0.0 ? 1.0 / (-1.0) : 1.0; // Jacobi inverse

			int v = layoutVerts[i];
			double vrad = localRadii[v];
			var fl = tri.getFlower(v);
			int m = fl.size();
			if (m == 0) {
				continue;
			}

			double[] iR = inRadii[i];
			double invTot = 1.0 / Math.max(1e-16, conduct[i]);

			double prevR = iR[m - 1];
			for (int j = 0; j < m; j++) {
				int w = fl.get(j);
				double coeff = (prevR + iR[j]) * invTot / Math.max(1e-16, vrad + localRadii[w]);
				prevR = iR[j];

				int idxW = v2indx[w];
				if (0 <= idxW && idxW < n) {
					int p = next[i]++;
					colIdx[p] = idxW;
					val[p] = coeff;
				} else {
					// already accounted in bx/by in first pass
				}
			}
		}

		SparseCSR A = new SparseCSR(n, nnz, rowPtr, colIdx, val, Minv);

		Preconditioner precond;
//		precond = new SSOR(A.n, A.rowPtr, A.colIdx, A.val, 1.333);
//		precond = Jacobi.fromConstantDiagonal(A.n, -1.0); // fastest (~2x)
		precond = new AMG(A.n, A.rowPtr, A.colIdx, A.val);
//		precond = new ILU0(A.n, A.rowPtr, A.colIdx, A.val);

		double tol = 1e-7; // NOTE magic constant. seems sufficient...
		int maxIters = Math.max(1000, 10 * A.n);

		double[] solx = new double[A.n];
		double[] soly = new double[A.n];

		// writes to solx+soly
		BiCGStabSolver.solve(A, bx, solx, tol, maxIters, precond);
		BiCGStabSolver.solve(A, by, soly, tol, maxIters, precond);

		// Copy results back
		for (int i = 0; i < n; i++) {
			int v = layoutVerts[i];
			localCentersX[v] = solx[i];
			localCentersY[v] = soly[i];
		}
	}

	/**
	 * <p>
	 * Update per-layout-vertex incircle radii and total conductances.
	 * </p>
	 *
	 * <p>
	 * What:
	 * </p>
	 * <ul>
	 * <li>For each layout vertex v, compute an array
	 * {@code inRadii[v][j] = sqrt((rv*ru*rw)/(rv+ru+rw))} for each face (triple)
	 * adjacent to v.</li>
	 * <li>Compute {@code conduct[v]} as the sum of edge conductances (t1+t2)/(rv +
	 * r_w) over petals.</li>
	 * </ul>
	 *
	 * <p>
	 * Why:
	 * </p>
	 * <p>
	 * These values are the geometric weights used when assembling {@code A} and the
	 * RHS in {@code layoutCentersSolve()}. They reflect the edge-level tangency
	 * geometry and ensure the resulting embedding respects local circle
	 * configuration.
	 * </p>
	 */
	private void updateVdata() {
		// Compute inRadii and node conduct for each interior (layout) vertex
		inRadii = new double[layCount][];
		conduct = new double[layCount];

		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double vr = localRadii[v];
			List<Integer> fl = tri.getFlower(v);
			int m = fl.size();
			double[] data = new double[m];
			for (int j = 0; j < m; j++) {
				int wj = fl.get(j);
				int wjp1 = fl.get((j + 1) % m);
				double ur = localRadii[wjp1];
				double wr = localRadii[wj];
				// incircle radius: sqrt(vr*ur*wr / (vr+ur+wr))
				data[j] = Math.sqrt(Math.max(0.0, vr * ur * wr / Math.max(1e-16, (vr + ur + wr))));
			}
			inRadii[i] = data;

			// total conductance
			double sum = 0.0;
			for (int j = 0; j < m; j++) {
				int w = fl.get(j);
				double t1 = data[(j - 1 + m) % m];
				double t2 = data[j];
				sum += (t1 + t2) / Math.max(1e-16, (vr + localRadii[w]));
			}
			conduct[i] = sum;
		}
	}

	/**
	 * <p>
	 * Set {@code localRadii} to effective radii derived from current centers (port
	 * of {@code setEffective.m}).
	 * </p>
	 *
	 * <ul>
	 * <li>Positive aim (interiors; boundary vertices in polygonal mode): compute
	 * sector areas from neighboring centers and set r_eff = sqrt(area /
	 * (vAims(v)/2)).</li>
	 * <li>Negative aim (boundary in max packing): GO-style update r_eff =
	 * sqrt(2*area/angsum), averaged with the old radius to moderate
	 * oscillation.</li>
	 * <li>Aims near zero (e.g. the 3-boundary-vertex disc case): radius left
	 * unchanged.</li>
	 * </ul>
	 */
	private void setEffectiveRadii() {
		// ---------------- interior ----------------
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double targetArea = vAims[v] / 2.0;
			if (targetArea <= 1e-3) {
				continue;
			}

			double area = sectorAreaAt(v);
			if (!(area > 0.0) || !Double.isFinite(area)) {
				continue;
			}

			double newr = Math.sqrt(area / targetArea);
			if (newr > 0.0 && Double.isFinite(newr)) {
				localRadii[v] = newr;
			}
		}

		// ---------------- boundary ----------------
		for (int k = 0; k < bdryCount; k++) {
			int w = bdryListClosed[k];
			double targetArea = vAims[w] / 2.0;

			List<Integer> fl = tri.getFlower(w); // open flower
			int m = fl.size();
			if (m < 2) {
				continue;
			}

			double area = 0.0;
			double angsum = 0.0;
			double zx = localCentersX[w], zy = localCentersY[w];

			for (int j = 0; j < m - 1; j++) {
				int jr = fl.get(j);
				int jl = fl.get(j + 1);

				double ang = MathUtil.angleAtCorner(zx, zy, localCentersX[jr], localCentersY[jr], localCentersX[jl], localCentersY[jl]);

				double r = 0.5 * (hyp(jr, w) + hyp(jl, w) - hyp(jr, jl));
				if (r > 0.0 && ang > 0.0 && Double.isFinite(r) && Double.isFinite(ang)) {
					area += 0.5 * r * r * ang;
					angsum += ang;
				}
			}

			if (!(area > 0.0) || !(angsum > 0.0) || !Double.isFinite(area) || !Double.isFinite(angsum)) {
				continue;
			}

			// Positive aim (polygonal mode): match sector area to the target.
			// Negative aim (max packing): GO-style update, averaged to moderate.
			double newr;
			if (targetArea > 1e-3) {
				newr = Math.sqrt(area / targetArea);
				if (newr > 0.0 && Double.isFinite(newr)) {
					localRadii[w] = newr;
				}
			} else if (targetArea < -1e-3) {
				newr = Math.sqrt(2.0 * area / Math.max(1e-16, angsum));
				if (newr > 0.0 && Double.isFinite(newr)) {
					localRadii[w] = 0.5 * (newr + localRadii[w]);
				}
			}
		}
	}

	private double hyp(int a, int b) {
		double dx = localCentersX[a] - localCentersX[b];
		double dy = localCentersY[a] - localCentersY[b];
		return Math.sqrt(dx * dx + dy * dy);
	}

	private double sectorAreaAt(int v) {
		List<Integer> fl = tri.getFlower(v);
		int m = fl.size();
		double zx = localCentersX[v], zy = localCentersY[v];
		double area = 0.0;
		for (int j = 0; j < m; j++) {
			int jr = fl.get(j);
			int jl = fl.get((j + 1) % m);

			double xjr = localCentersX[jr], yjr = localCentersY[jr];
			double xjl = localCentersX[jl], yjl = localCentersY[jl];

			double dxr = xjr - zx, dyr = yjr - zy;
			double dxl = xjl - zx, dyl = yjl - zy;
			double dxrl = xjr - xjl, dyrl = yjr - yjl;

			double ang = MathUtil.angleAtCorner(zx, zy, xjr, yjr, xjl, yjl);

			double a = Math.sqrt(dxr * dxr + dyr * dyr);
			double b = Math.sqrt(dxl * dxl + dyl * dyl);
			double c = Math.sqrt(dxrl * dxrl + dyrl * dyrl);

			double r = 0.5 * (a + b - c);
			area += 0.5 * r * r * ang;
		}
		return area;
	}

	private double updateVisErrorMonitor() {
		double[] ve = visualErrors();
		double max = 0.0;
		for (double x : ve) {
			max = Math.max(max, x);
		}
		visErrMonitor.add(max);
		return max;
	}

	private static final class OrderedCorners {
		final int[] corners;
		final double[] angles;

		OrderedCorners(int[] corners, double[] angles) {
			this.corners = corners;
			this.angles = angles;
		}
	}

	private OrderedCorners orderCornersCCW(int[] rawCorners, double[] rawAngles) {
		final int m = rawCorners.length;

		int[] pos = new int[m];
		for (int i = 0; i < m; i++) {
			pos[i] = boundaryPosition(rawCorners[i]);
			if (pos[i] < 0) {
				throw new IllegalArgumentException("Corner " + rawCorners[i] + " is not on the boundary loop");
			}
		}

		Integer[] ord = new Integer[m];
		for (int i = 0; i < m; i++) {
			ord[i] = i;
		}
		Arrays.sort(ord, Comparator.comparingInt(i -> pos[i]));

		// Rotate so the first user-supplied corner remains first
		int anchor = 0;
		for (int k = 0; k < m; k++) {
			if (rawCorners[ord[k]] == rawCorners[0]) {
				anchor = k;
				break;
			}
		}

		int[] outCorners = new int[m];
		double[] outAngles = (rawAngles == null) ? null : new double[m];

		for (int t = 0; t < m; t++) {
			int src = ord[(anchor + t) % m];
			outCorners[t] = rawCorners[src];
			if (outAngles != null) {
				outAngles[t] = rawAngles[src];
			}
		}

		return new OrderedCorners(outCorners, outAngles);
	}

	private int boundaryPosition(int v) {
		for (int i = 0; i < bdryCount; i++) {
			if (bdryListClosed[i] == v) {
				return i;
			}
		}
		return -1;
	}
}