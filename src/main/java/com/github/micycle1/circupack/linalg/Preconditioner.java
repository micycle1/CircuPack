package com.github.micycle1.circupack.linalg;

/**
 * <p>
 * Collection of light-weight preconditioner implementations for CSR sparse
 * matrices used with iterative solvers (e.g. BiCGStab). Each implementation
 * implements Preconditioner.apply(r,z) which computes
 * <code>z = M^{-1} r</code>.
 * </p>
 */
public interface Preconditioner {
	void apply(double[] r, double[] z);
}