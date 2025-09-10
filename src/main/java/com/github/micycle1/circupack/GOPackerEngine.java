package com.github.micycle1.circupack;

import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.sparse.FillReducing;
import org.ejml.sparse.csc.CommonOps_DSCC;
import org.ejml.sparse.csc.factory.LinearSolverFactory_DSCC;
import org.ojalgo.matrix.decomposition.LU;
import org.ojalgo.matrix.decomposition.QR;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.R064Store;
import org.ojalgo.matrix.store.SparseStore;
import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.interfaces.linsol.LinearSolverSparse;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Java port of the MATLAB @GOPacker engine. Focus: Euclidean MAX_PACK and
 * POLYGONAL modes. Supports spherical faux boundary, effective radii,
 * Tutte-like embedding, boundary layout (horocycles in unit disk, rectangle,
 * general polygon), conductance/incircle computation, visual error, angle-sum
 * error, pruning orphans, etc.
 *
 * Works with the provided Triangulation interface (0-based indices).
 */
public class GOPackerEngine {

	// External triangulation (combinatorics and optional geometry)
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
	private boolean[] isInterior; // interior vs boundary (simple classification)
	private boolean[] isBoundary; // boundary marking
	private boolean hasBoundary; // whether triangulation has a boundary

	// Radii and centers (Euclidean), both working (local) and final (best-known)
	public double[] radii; // final radii (after riffle)
	public double[] centersX; // final centers
	public double[] centersY;

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
	private Triangulation.Mode mode = Triangulation.Mode.MAX_PACK;

	// For polygonal mode
	private int[] corners = null; // CCW corners
	private List<int[]> sides = null; // each side: list of vertices including endpoints; CCW order

	// Monitoring
	private final List<Double> visErrMonitor = new ArrayList<>();

	// Constructor
	public GOPackerEngine(Triangulation tri) {
		this.tri = Objects.requireNonNull(tri, "tri must not be null");
	}

	// Initialize combinatorics, defaults, radii/centers, aims, mode
	public void initialize() {
		n = tri.getVertexCount();
		if (n <= 0)
			throw new IllegalStateException("Triangulation has no vertices");

		// read mode hint
		mode = tri.getMode() != null ? tri.getMode() : Triangulation.Mode.MAX_PACK;

		// set alpha, gamma if provided
		alpha = tri.getAlpha();
		gamma = tri.getGamma();

		// boundary presence
		List<Integer> bdryLoop = tri.getBoundaryLoop();
		hasBoundary = bdryLoop != null && !bdryLoop.isEmpty();

		// classify interior/boundary naive
		isBoundary = new boolean[n];
		for (int v = 0; v < n; v++) {
			isBoundary[v] = tri.isBoundaryVertex(v);
		}
		isInterior = new boolean[n];
		for (int v = 0; v < n; v++) {
			isInterior[v] = !isBoundary[v];
		}

		// Find interior connected component containing alpha; if alpha invalid, pick a
		// default interior
		if (alpha < 0 || alpha >= n || isBoundary[alpha]) {
			alpha = pickAnInteriorVertexOrFallback();
		}

		// Derive intVerts, bdryListClosed, and orphanVerts, similar to complex_count.m
		// logic
		deriveInteriorBoundaryOrphans();

		// Copy provided geometry if available; else default
		if (tri.hasRadii()) {
			radii = tri.getRadii().clone();
		} else {
			radii = new double[n];
			Arrays.fill(radii, 0.5);
		}
		if (tri.hasCenters()) {
			double[] cx = tri.getCentersX();
			double[] cy = tri.getCentersY();
			centersX = cx != null ? cx.clone() : new double[n];
			centersY = cy != null ? cy.clone() : new double[n];
		} else {
			centersX = new double[n];
			centersY = new double[n];
		}

		localRadii = radii.clone();
		localCentersX = centersX.clone();
		localCentersY = centersY.clone();

		// aims
		if (tri.hasVAims()) {
			double[] a = tri.getVAims();
			if (a != null && a.length == n)
				vAims = a.clone();
		}
		if (vAims == null) {
			vAims = new double[n];
			Arrays.fill(vAims, 2.0 * Math.PI);
			for (int v = 0; v < n; v++) {
				if (isBoundary[v])
					vAims[v] = -1.0;
			}
		}

		// polygonal metadata if present
		if (mode == Triangulation.Mode.POLYGONAL || mode == Triangulation.Mode.FIXED_CORNERS) {
			List<Integer> cn = tri.getCorners();
			if (cn != null && !cn.isEmpty()) {
				corners = cn.stream().mapToInt(Integer::intValue).toArray();
			}
			List<List<Integer>> sd = tri.getSides();
			if (sd != null && !sd.isEmpty()) {
				sides = new ArrayList<>();
				for (List<Integer> s : sd) {
					sides.add(s.stream().mapToInt(Integer::intValue).toArray());
				}
			}
		}

		// Default layoutVerts = intVerts; rimVerts = boundary loop (closed)
		layoutVerts = intVerts.clone();
		layCount = layoutVerts.length;
		rimVerts = new int[bdryCount];
		System.arraycopy(bdryListClosed, 0, rimVerts, 0, bdryCount);
		rimCount = rimVerts.length;

		// indexing arrays
		buildIndexing();
	}

	// Public API methods

	public void setMode(Triangulation.Mode desiredMode, int[] cornersIn, List<int[]> sidesIn, double[] cornerAngles) {
		if (desiredMode == null)
			desiredMode = Triangulation.Mode.MAX_PACK;

		if (hes > 0 && desiredMode != Triangulation.Mode.MAX_PACK) {
			// Sphere must be in MAX_PACK
			this.mode = Triangulation.Mode.MAX_PACK;
		} else {
			this.mode = desiredMode;
		}

		if (this.mode == Triangulation.Mode.MAX_PACK) {
			// aims: interior 2*pi; boundary -1
			Arrays.fill(vAims, 2.0 * Math.PI);
			for (int v = 0; v < n; v++)
				if (isBoundary[v])
					vAims[v] = -1.0;
			return;
		}

		// POLYGONAL / FIXED_CORNERS
		// Corners/sides optional. If cornersIn provided, validate boundary membership.
		if (cornersIn != null && cornersIn.length >= 3) {
			for (int c : cornersIn) {
				if (!isBoundary[c]) {
					throw new IllegalArgumentException("Corner " + c + " is not a boundary vertex");
				}
			}
			this.corners = cornersIn.clone();
		} else if (this.corners == null || this.corners.length < 3) {
			this.corners = chooseRandomCorners(4); // fallback 4-gon
		}

		if (sidesIn != null && !sidesIn.isEmpty()) {
			this.sides = new ArrayList<>(sidesIn);
		} else {
			// derive sides by splitting boundary loop at corners
			this.sides = derivePolygonSidesFromCorners(this.corners);
		}

		// set default aims: boundary pi except corners get equal angles by default
		for (int w : rimVerts)
			vAims[w] = Math.PI;
		if (cornerAngles != null) {
			if (cornerAngles.length != this.corners.length) {
				throw new IllegalArgumentException("Corner angles length mismatch");
			}
			for (int i = 0; i < this.corners.length; i++)
				vAims[this.corners[i]] = cornerAngles[i];
		} else {
			int m = this.corners.length;
			for (int i = 0; i < m; i++)
				vAims[this.corners[i]] = Math.PI * (1.0 - 2.0 / m);
		}
	}

	public void debugDumpSystemRow0() {
		if (layCount == 0)
			return;
		int v = layoutVerts[0];
		double vrad = localRadii[v];
		var fl = tri.getFlower(v);
		int m = fl.size();
		double[] data = inRadii[0];
		double total = Math.max(1e-16, conduct[0]);

		double sumCoeff = 0;
		double sumBx = 0, sumBy = 0;
		for (int j = 0; j < m; j++) {
			int w = fl.get(j);
			double t1 = data[(j - 1 + m) % m];
			double t2 = data[j];
			double coeff = ((t1 + t2) / Math.max(1e-16, (vrad + localRadii[w]))) / total;
			int idxW = v2indx[w];
			String kind = (0 <= idxW && idxW < layCount) ? "int" : "rim";
			System.out.printf("j=%d w=%d coeff=%.8f kind=%s  zw=(%.6f,%.6f)\n", j, w, coeff, kind, localCentersX[w], localCentersY[w]);
			sumCoeff += coeff;
			if (!"int".equals(kind)) {
				sumBx += coeff * localCentersX[w];
				sumBy += coeff * localCentersY[w];
			}
		}
		System.out.printf("diag=-1.0, sumCoeff(all nbrs)=%.8f\n", sumCoeff);
		System.out.printf("Expected interior center from rims only: (%.8f, %.8f)\n", sumBx, sumBy);
	}

	private void debugBoundaryStats() {
		double min = 1e9, max = -1e9, sum = 0;
		for (int i = 0; i < bdryCount; i++) {
			double r = localRadii[bdryListClosed[i]];
			min = Math.min(min, r);
			max = Math.max(max, r);
			sum += r;
		}
		System.out.printf("bdry r: min=%.6g max=%.6g mean=%.6g%n", min, max, sum / bdryCount);
	}

	public int riffle(int passNum) {
		if (passNum <= 0) {
			layoutBoundary();
			layoutCentersSolve();
			updateVisErrorMonitor();
			return 0;
		}

		double cutval = 0.01;
		int pass = 0;
		double maxVis = 2 * cutval;

		while (pass < passNum && maxVis > cutval) {
			layoutBoundary(); // set boundary centers (and possibly scale radii)
			layoutCentersSolve(); // solve A * Z = rhs for interior centers
			setEffectiveRadii(); // update effective radii
			maxVis = updateVisErrorMonitor();
			pass++;
			System.out.print(".");
		}
		radii = localRadii.clone();
		centersX = localCentersX.clone();
		centersY = localCentersY.clone();
		return pass;
	}

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
				double cdiff = Math.hypot(dx, dy);
				double rdiff = rv + localRadii[w];
				double me = Math.abs(cdiff - rdiff) / Math.max(1e-16, rv);
				if (me > maxErr)
					maxErr = me;
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

	public void pruneOrphans() {
		// If no orphans, nothing to do
		if (orphanCount == 0)
			return;

		// Remove orphan vertices from triangulation view within the engine
		// We don't mutate Triangulation; we adjust our working sets and geometry
		// arrays.
		int newN = intCount + bdryCount; // keep alpha-component interiors and boundary
		boolean[] keep = new boolean[n];
		for (int v : intVerts)
			keep[v] = true;
		for (int i = 0; i < bdryCount; i++)
			keep[bdryListClosed[i]] = true;

		int[] old2new = new int[n];
		Arrays.fill(old2new, -1);
		int kk = 0;
		for (int v = 0; v < n; v++)
			if (keep[v])
				old2new[v] = kk++;

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
		for (int i = 0; i < intCount; i++)
			nInt[i] = old2new[intVerts[i]];
		int[] nBdryClosed = new int[bdryCount + 1];
		for (int i = 0; i < bdryCount + 1; i++)
			nBdryClosed[i] = old2new[bdryListClosed[i]];

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

	// Export back to Triangulation (optional)
	public void writeBackToTriangulation() {
		tri.setRadii(radii.clone());
		tri.setCenters(centersX.clone(), centersY.clone());
		tri.setVAims(vAims.clone());
		tri.setAlpha(alpha);
		if (gamma >= 0)
			tri.setGamma(gamma);
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

	// --------------- Internals ---------------

	private int[] chooseRandomCorners(int sideN) {
		if (bdryCount < 3)
			throw new IllegalStateException("Boundary must have at least 3 vertices");
		if (sideN < 3)
			sideN = 3;
		if (sideN > bdryCount)
			sideN = bdryCount;

		int[] crn = new int[sideN];
		int step = Math.max(1, bdryCount / sideN);
		int start = 0; // deterministic start at gamma could be used instead
		// If you prefer to start at gamma:
		// int start = 0;
		// for (int i = 0; i < bdryCount; i++) if (bdryListClosed[i] == gamma) { start =
		// i; break; }

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
				if (v != a && v != b && v != c)
					ints.add(v);
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
			if (!isBoundary[v] && touchedInterior[v])
				ints.add(v);
		}
		intVerts = ints.stream().mapToInt(Integer::intValue).toArray();
		intCount = intVerts.length;

		// Build bdry loop closed (from Triangulation)
		bdryCount = tri.getBoundaryLoop().size();
		bdryListClosed = new int[bdryCount + 1];
		for (int i = 0; i < bdryCount; i++)
			bdryListClosed[i] = tri.getBoundaryLoop().get(i);
		bdryListClosed[bdryCount] = bdryListClosed[0];

		// Orphans: not in intVerts and not in boundary loop
		boolean[] isInInt = new boolean[n];
		for (int v : intVerts)
			isInInt[v] = true;
		boolean[] isInB = new boolean[n];
		for (int i = 0; i < bdryCount; i++)
			isInB[bdryListClosed[i]] = true;

		List<Integer> orph = new ArrayList<>();
		for (int v = 0; v < n; v++) {
			if (!isInInt[v] && !isInB[v])
				orph.add(v);
		}
		orphanVerts = orph.stream().mapToInt(Integer::intValue).toArray();
		orphanCount = orphanVerts.length;

		// pick gamma if not provided
		if (gamma < 0 || gamma >= n || !isBoundary[gamma]) {
			gamma = bdryListClosed[0];
		}
	}

	private int farVertex(List<Integer> seeds) {
		if (seeds == null || seeds.isEmpty())
			return (n > 0 ? 0 : -1);
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

	private void layoutBoundary() {
		if (mode == Triangulation.Mode.POLYGONAL || mode == Triangulation.Mode.FIXED_CORNERS) {
			setPolygonCenters();
		} else {
			setHoroCenters();
		}
	}

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
		if (R < 2.0 * minrad)
			R = 3.0 * minrad;

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
		for (int v = 0; v < n; v++)
			localRadii[v] /= R;
		for (int j = 0; j <= bdryCount; j++)
			r[j] /= R;

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

	private void setPolygonCenters() {
		// if rectangle (4 corners) with all right angles, do rectangular layout
		if (corners != null && corners.length == 4) {
			double bendErr = 0.0;
			for (int c : corners)
				bendErr += Math.abs(vAims[c] - Math.PI / 2.0);
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
			for (int j = 1; j < side.length - 1; j++)
				L += 2.0 * localRadii[side[j]];
			sideLengths[i] = L;
			fullLength += L;
		}

		int halfn = numSides / 2;

		if (numSides == 3) {
			// triangle with corner aims provided in vAims at corners
			double opp1 = vAims[corners[2]];
			double opp2 = vAims[corners[0]];
			if (opp1 <= 0 || opp2 <= 0 || (opp1 + opp2) >= Math.PI)
				return;
			double opp3 = Math.PI - (opp1 + opp2);
			targetLength[0] = 1.0;
			targetLength[1] = Math.sin(opp2) / Math.sin(opp1);
			targetLength[2] = Math.sin(opp3) / Math.sin(opp1);
			double sum = targetLength[0] + targetLength[1] + targetLength[2];
			double factor = sum / fullLength;
			scaleAllLocalBy(factor);
			for (int i = 0; i < numSides; i++)
				sideLengths[i] *= factor;
		} else if (2 * halfn == numSides) {
			// even: pair opposite sides, target total per pair equals average; set scale so
			// sum target ~ 6
			double factor = 6.0 / fullLength;
			scaleAllLocalBy(factor);
			for (int i = 0; i < numSides; i++)
				sideLengths[i] *= factor;
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

	private void setRectCenters() {
		// sides.size() == 4; corners count == 4
		double[] sideLen = new double[4];
		for (int i = 0; i < 4; i++) {
			int[] side = sides.get(i);
			double L = localRadii[side[0]] + localRadii[side[side.length - 1]];
			for (int j = 1; j < side.length - 1; j++)
				L += 2.0 * localRadii[side[j]];
			sideLen[i] = L;
		}
		double width = (sideLen[0] + sideLen[2]) / 2.0;
		double height = (sideLen[1] + sideLen[3]) / 2.0;
		double aspect = height / Math.max(1e-16, width);

		// scale radii so that total fits rectangle of width 2*aspect and height 2
		double factor = 2.0 * (aspect + 1.0) / Math.max(1e-16, (width + height));
		scaleAllLocalBy(factor);
		for (int i = 0; i < 4; i++)
			sideLen[i] *= factor;

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

			double prev = localRadii[start];
			for (int i = 0; i < side.length - 1; i++) {
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
		for (int v = 0; v < n; v++)
			localRadii[v] *= factor;
	}

	private List<int[]> derivePolygonSidesFromCorners(int[] crn) {
		// Split bdry loop at corners; ensure CCW order starting at gamma if set.
		List<Integer> loop = Arrays.stream(bdryListClosed, 0, bdryCount).boxed().collect(Collectors.toList());
		// reorder corners in CCW boundary order
		List<Integer> cornerList = Arrays.stream(crn).boxed().collect(Collectors.toList());
		cornerList.sort(Comparator.comparingInt(loop::indexOf));
		int m = cornerList.size();

		List<int[]> result = new ArrayList<>();
		for (int i = 0; i < m; i++) {
			int a = cornerList.get(i);
			int b = cornerList.get((i + 1) % m);
			int ia = loop.indexOf(a);
			int ib = loop.indexOf(b);
			List<Integer> side = new ArrayList<>();
			side.add(a);
			int k = ia;
			while (true) {
				k = (k + 1) % bdryCount;
				side.add(loop.get(k));
				if (loop.get(k) == b)
					break;
				if (side.size() > bdryCount + 5)
					throw new IllegalStateException("Side extraction overflow");
			}
			result.add(side.stream().mapToInt(Integer::intValue).toArray());
		}
		return result;
	}

	private void layoutCentersSolve() {
		updateVdata();

		int estNnz = 0;
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			estNnz += tri.getFlower(v).size() + 1;
		}

		// Triplet builder for A
		DMatrixSparseTriplet Atr = new DMatrixSparseTriplet(layCount, layCount, estNnz);
		double[] bx = new double[layCount];
		double[] by = new double[layCount];

		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double vrad = localRadii[v];
			var fl = tri.getFlower(v);
			int m = fl.size();
			double[] iR = inRadii[i];
			double totalCond = Math.max(1e-16, conduct[i]);

			// diagonal
			Atr.addItem(i, i, -1.0);

			for (int j = 0; j < m; j++) {
				int w = fl.get(j);
				double t1 = iR[(j - 1 + m) % m];
				double t2 = iR[j];
				double coeff = ((t1 + t2) / Math.max(1e-16, (vrad + localRadii[w]))) / totalCond;

				int idxW = v2indx[w];
				if (0 <= idxW && idxW < layCount) {
					Atr.addItem(i, idxW, coeff);
				} else {
					bx[i] += -coeff * localCentersX[w];
					by[i] += -coeff * localCentersY[w];
				}
			}
		}

		DMatrixSparseCSC A = DConvertMatrixStruct.convert(Atr, (DMatrixSparseCSC) null);
		DMatrixRMaj rhsX = new DMatrixRMaj(bx.length, 1, true, bx);
		DMatrixRMaj rhsY = new DMatrixRMaj(by.length, 1, true, by);
		DMatrixRMaj solX = solveSparseSystem(A, rhsX);
		DMatrixRMaj solY = solveSparseSystem(A, rhsY);

		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			localCentersX[v] = solX.get(i);
			localCentersY[v] = solY.get(i);
		}
	}

	private DMatrixRMaj solveSparseSystem(DMatrixSparseCSC A, DMatrixRMaj b) {
		LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solver = LinearSolverFactory_DSCC.lu(FillReducing.NONE);
		if (!solver.setA(A)) {
			throw new RuntimeException("Linear solver failed to set A (possibly singular)");
		}
		DMatrixRMaj x = new DMatrixRMaj(A.numCols, 1);
		solver.solve(b, x);
		return x;
	}

	private void layoutCentersSolveFast() {

		// ------------------------------------------------------------------
		// 1) collect data ---------------------------------------------------
		updateVdata();

		/* ---------- estimate nnz and create once-only data structures ----- */
		int estNnz = 0;
		for (int k = 0; k < layCount; k++)
			estNnz += tri.getFlower(layoutVerts[k]).size() + 1;

		DMatrixSparseTriplet Atr = new DMatrixSparseTriplet(layCount, layCount, estNnz);
		double[] bx = new double[layCount];
		double[] by = new double[layCount];

		/* ---------- assemble A , bx , by ---------------------------------- */
		for (int i = 0; i < layCount; i++) {

			int v = layoutVerts[i];
			double vrad = localRadii[v];
			var fl = tri.getFlower(v); // neighbour list
			int m = fl.size();
			double[] iR = inRadii[i];
			double invTot = 1.0 / Math.max(1e-16, conduct[i]);

			/* diagonal entry (always –1) */
			Atr.addItem(i, i, -1.0);

			/* nothing to do if the flower is empty -------------------------------- */
			if (m == 0)
				continue;

			/* iterate flower, using prevR to avoid (j-1+m)%m */
			double prevR = iR[m - 1];
			for (int j = 0; j < m; j++) {

				int w = fl.get(j);
				double coeff = (prevR + iR[j]) * invTot / Math.max(1e-16, vrad + localRadii[w]);
				prevR = iR[j];

				int idxW = v2indx[w];
				if (0 <= idxW && idxW < layCount) { // internal neighbour
					Atr.addItem(i, idxW, coeff);
				} else { // already solved
					bx[i] -= coeff * localCentersX[w];
					by[i] -= coeff * localCentersY[w];
				}
			}
		}

		// ------------------------------------------------------------------
		// 2) convert to CSC only once and solve both RHS at once -------------
		DMatrixSparseCSC A = DConvertMatrixStruct.convert(Atr, (DMatrixSparseCSC) null);

		/* build RHS in correct row-major order -------------------------------- */
		DMatrixRMaj RHS = new DMatrixRMaj(layCount, 2);
		for (int r = 0, p = 0; r < layCount; r++, p += 2) {
			RHS.data[p] = bx[r]; // column 0 (x)
			RHS.data[p + 1] = by[r]; // column 1 (y)
		}

		/* one LU factorisation – two solves ----------------------------------- */
		LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solver = LinearSolverFactory_DSCC.lu(FillReducing.NONE);
		if (!solver.setA(A))
			throw new RuntimeException("Linear solver failed to set A (possibly singular)");

		DMatrixRMaj sol = new DMatrixRMaj(layCount, 2);
		solver.solve(RHS, sol);

		// ------------------------------------------------------------------
		// 3) copy results back ----------------------------------------------
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			localCentersX[v] = sol.get(i, 0);
			localCentersY[v] = sol.get(i, 1);
		}
	}

	private void layoutCentersSolveoJ() {
		updateVdata();

		final int n = layCount;

		// Estimate nnz: diagonal + one per neighbour; add some headroom
		int estNnz = n;
		for (int i = 0; i < n; i++) {
			estNnz += tri.getFlower(layoutVerts[i]).size();
		}
		estNnz += n;

		// Assemble in SparseStore (COO-ish), then convert to CSR
		final SparseStore<Double> Acoo = SparseStore.R064.make(n, n);
		final double[] bx = new double[n];
		final double[] by = new double[n];

		for (int i = 0; i < n; i++) {
			final int v = layoutVerts[i];
			final double vrad = localRadii[v];
			final var fl = tri.getFlower(v);
			final int m = fl.size();
			final double[] iR = inRadii[i];

			final double invCond = 1.0 / Math.max(1e-16, conduct[i]);

			// Diagonal
			Acoo.add(i, i, -1.0);

			for (int j = 0; j < m; j++) {
				final int w = fl.get(j);

				final double t1 = iR[(j == 0 ? m - 1 : j - 1)];
				final double t2 = iR[j];
				final double denom = Math.max(1e-16, vrad + localRadii[w]);
				final double coeff = (t1 + t2) * (1.0 / denom) * invCond;

				final int idxW = v2indx[w];
				if (0 <= idxW && idxW < n) {
					// Interior neighbour -> A contribution (add() accumulates)
					Acoo.add(i, idxW, coeff);
				} else {
					// Boundary neighbour -> RHS contribution
					bx[i] -= coeff * localCentersX[w];
					by[i] -= coeff * localCentersY[w];
				}
			}
		}

		final MatrixStore<Double> A = Acoo.toCSR(); // or toCSC(); both work with sparse LU

		// Two RHS columns
		final R064Store B = R064Store.FACTORY.make(n, 2);
		for (int i = 0; i < n; i++) {
			B.set(i, 0, bx[i]);
			B.set(i, 1, by[i]);
		}

		// Sparse LU once, solve both RHS at once
		final LU<Double> lu = LU.R064.make();
		if (!lu.decompose(A) || !lu.isSolvable()) {
			throw new IllegalStateException("Failed to solve (singular/ill-conditioned system)");
		}
		final MatrixStore<Double> X = lu.getSolution(B);

		// Write back
		for (int i = 0; i < n; i++) {
			final int v = layoutVerts[i];
			localCentersX[v] = X.doubleValue(i, 0);
			localCentersY[v] = X.doubleValue(i, 1);
		}
	}

	// Robust solver: LU first (fast), QR fallback (robust)
	private R064Store solveSparseSystemOjAlgoR064(final MatrixStore<Double> A, final R064Store b) {
		// Try LU
		final LU<Double> lu = LU.newSparseR064();

		if (lu.decompose(A) && lu.isSolvable()) {
			final MatrixStore<Double> x = lu.getSolution(b); // single-argument solve(rhs)
			return R064Store.FACTORY.copy(x);
		}
		// Fallback QR
//	    final QR<Double> qr = QR.R064.make();
//	    if (qr.decompose(A) && qr.isSolvable()) {
//	        final MatrixStore<Double> x = qr.getSolution(b);
//	        return R064Store.FACTORY.copy(x);
//	    }
		throw new IllegalStateException("ojAlgo v55+: failed to solve (singular/ill-conditioned system)");
	}

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

	private void setEffectiveRadii() {
		// interior adjustments: sector-area sum vs target area (aim/2)
		for (int i = 0; i < layCount; i++) {
			int v = layoutVerts[i];
			double targetArea = vAims[v] / 2.0;
			if (targetArea <= 1e-3)
				continue;
			double area = sectorAreaAt(v);
			localRadii[v] = Math.sqrt(Math.max(0.0, area / targetArea));
		}
		// boundary adjustments (bdry aims <= 0 typically)
		double EPS = 1e-12;
		for (int k = 0; k < bdryCount; k++) {
			int w = bdryListClosed[k];
			double targetArea = vAims[w] / 2.0;

			var fl = tri.getFlower(w); // open
			int m = fl.size();
			double area = 0.0, angsum = 0.0;
			double zx = localCentersX[w], zy = localCentersY[w];

			for (int j = 0; j < m - 1; j++) {
				int jr = fl.get(j), jl = fl.get(j + 1);
				double ang = MathUtil.angleAtCorner(zx, zy, localCentersX[jr], localCentersY[jr], localCentersX[jl], localCentersY[jl]);
				double r = 0.5 * (hyp(jr, w) + hyp(jl, w) - hyp(jr, jl));
				if (r > 0 && ang > 0) {
					area += 0.5 * r * r * ang;
					angsum += ang;
				}
			}

			if (targetArea > 1e-3) {
				localRadii[w] = Math.sqrt(Math.max(0.0, area / targetArea));
			} else if (targetArea < -1e-3) {
				double newr = Math.sqrt(Math.max(0.0, 2.0 * area / Math.max(1e-16, angsum)));
				localRadii[w] = 0.5 * (newr + localRadii[w]);
			}
		}
	}

	private double hyp(int a, int b) {
		return Math.hypot(localCentersX[a] - localCentersX[b], localCentersY[a] - localCentersY[b]);
	}

	private double sectorAreaAt(int v) {
		List<Integer> fl = tri.getFlower(v);
		int m = fl.size();
		double zx = localCentersX[v], zy = localCentersY[v];
		double area = 0.0;
		for (int j = 0; j < m; j++) {
			int jr = fl.get(j);
			int jl = fl.get((j + 1) % m);
			double ang = MathUtil.angleAtCorner(zx, zy, localCentersX[jr], localCentersY[jr], localCentersX[jl], localCentersY[jl]);
			double r = 0.5 * (Math.hypot(localCentersX[jr] - zx, localCentersY[jr] - zy) + Math.hypot(localCentersX[jl] - zx, localCentersY[jl] - zy)
					- Math.hypot(localCentersX[jr] - localCentersX[jl], localCentersY[jr] - localCentersY[jl]));
			area += 0.5 * r * r * ang;
		}
		return area;
	}

	private double updateVisErrorMonitor() {
		double[] ve = visualErrors();
		double max = 0.0;
		for (double x : ve)
			max = Math.max(max, x);
		visErrMonitor.add(max);
		return max;
	}
}