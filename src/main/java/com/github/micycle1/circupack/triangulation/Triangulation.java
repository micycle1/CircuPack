package com.github.micycle1.circupack.triangulation;

import java.util.List;

/**
 * Triangulation: a minimal combinatorial + optional geometric interface
 * required by Circupack.
 * <p>
 * Conventions:
 * <ul>
 * <li>Vertex indices are 0..n-1 (0-based).</li>
 * <li>For each vertex v, getFlower(v) returns neighbors in COUNTER-CLOCKWISE
 * order.
 * <ul>
 * <li>interior vertex: flower is cyclic (length = degree). Do NOT repeat first
 * at end.</li>
 * <li>boundary vertex: flower is open (first != last).</li>
 * </ul>
 * </li>
 * <li>getBoundaryLoop() returns CCW-ordered boundary vertices, each exactly
 * once; return null or an empty array to indicate a sphere (no boundary).</li>
 * </ul>
 * <p>
 * Non-MAX_PACK / optional methods default to null/false/do-nothing via default
 * implementations.
 */
public interface Triangulation {

	// ----------------------
	// Core combinatorics (must be implemented)
	// ----------------------

	/** Number of vertices (n). */
	int getVertexCount();

	/**
	 * The neighbor list ("flower") for vertex v in CCW order. For interior vertices
	 * the list is cyclic; for boundary vertices it's open.
	 *
	 * NOTE MATLAB code often stored flowers with the first neighbor repeated at the
	 * end for interior vertices (flower(1) == flower(end)). Our interface forbids
	 * that repetition
	 * <p>
	 * The GOPacker engine assumes boundary flowers are open, CCW, and specifically
	 * run from prevBoundary to nextBoundary through the interior. It then sums over
	 * consecutive pairs (j, j+1) and intentionally does not wrap (no pair from last
	 * to first), so we cover only the interior wedge. it must start at prevBoundary
	 * and end at nextBoundary.
	 * <p>
	 * For each boundary v, getFlower(v) should be exactly [prev, center, next] in
	 * CCW order.
	 */
	List<Integer> getFlower(int v);

	/** True if v is a boundary vertex. */
	boolean isBoundaryVertex(int v);

	/**
	 * CCW ordered boundary loop (each boundary vertex once). Null/empty =>
	 * sphere/no-boundary. should not repeat/close the loop. It should return a
	 * simple CCW-ordered list of the boundary vertices where each boundary vertex
	 * appears exactly once and the closure (edge between last and first) is
	 * implied, not explicitly repeated.
	 */
	List<Integer> getBoundaryLoop();

	/** Convenience: list of interior vertices. */
	List<Integer> getInteriorVertices();

	// ----------------------
	// Optional: faces (default: not provided)
	// ----------------------

	/** Face list (each int[3]); default null meaning "not provided". */
	default List<int[]> getFaces() {
		return null;
	}

	/** Face count; default -1 when faces not provided. */
	default int getFaceCount() {
		return -1;
	}

	// ----------------------
	// Optional geometric data (default: not provided)
	// ----------------------

	/** Euclidean centers X array length n or null if omitted. Default: null. */
	default double[] getCentersX() {
		return null;
	}

	/** Euclidean centers Y array length n or null if omitted. Default: null. */
	default double[] getCentersY() {
		return null;
	}

	/**
	 * Set centers: default does nothing (no-op). Implementers override to accept
	 * arrays.
	 */
	default void setCenters(double[] centersX, double[] centersY) {
		/* no-op by default */
	}

	/** True if centers present. Default false. */
	default boolean hasCenters() {
		return false;
	}

	/** Radii array length n or null if omitted. Default: null. */
	default double[] getRadii() {
		return null;
	}

	/** Set radii: default no-op. */
	default void setRadii(double[] radii) {
		/* no-op */
	}

	/**
	 * vAims: target angle-sum array length n or null if omitted. Default: null
	 * (packer should treat interior default = 2*pi and boundary negative marker).
	 */
	default double[] getVAims() {
		return null;
	}

	/** Set vAims: default no-op. */
	default void setVAims(double[] vAims) {
		/* no-op */
	}

	/** True if vAims present. Default false. */
	default boolean hasVAims() {
		return false;
	}

	// ----------------------
	// Normalization / special vertices (defaults)
	// ----------------------

	/** Preferred alpha index or -1 if not set. Default -1. */
	default int getAlpha() {
		return -1;
	}

	/** Set alpha: default no-op. */
	default void setAlpha(int alpha) {
		/* no-op */ }

	/** Preferred gamma index or -1 if not set. Default -1. */
	default int getGamma() {
		return -1;
	}

	/** Set gamma: default no-op. */
	default void setGamma(int gamma) {
		/* no-op */
	}

	// ----------------------
	// Polygon / rectangular support (defaults: not provided)
	// ----------------------

	/** Corner vertices for polygonal mode (CCW) or null (default). */
	default List<Integer> getCorners() {
		return null;
	}

	/** Set corners: default no-op. */
	default void setCorners(List<Integer> corners) {
		/* no-op */
	}

	/**
	 * Sides for polygonal mode: list of CCW lists of vertex indices (each includes
	 * endpoints). Default: null (not provided).
	 */
	default List<List<Integer>> getSides() {
		return null;
	}

	/** Set sides: default no-op. */
	default void setSides(List<List<Integer>> sides) {
		/* no-op */
	}

	/** Auxiliary vlist: default null. */
	default List<Integer> getVlist() {
		return null;
	}

	/** Set vlist: default no-op. */
	default void setVlist(List<Integer> vlist) {
		/* no-op */
	}

	// ----------------------
	// Mode / hints (default: MAX_PACK)
	// ----------------------

	enum Mode {
		MAX_PACK, POLYGONAL, FROZEN_BDRY, ORTH, FIXED_CORNERS
	}

	/** Preferred packer mode hint; default MAX_PACK. */
	default Mode getMode() {
		return Mode.MAX_PACK;
	}

	/** Set mode: default no-op. */
	default void setMode(Mode mode) {
		/* no-op */
	}
}