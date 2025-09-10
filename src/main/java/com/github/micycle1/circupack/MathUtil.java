package com.github.micycle1.circupack;

public final class MathUtil {

	private MathUtil() {
	}

	// Equivalent to cosAngle.m
	// Cosine of angle at circle with radius r in a triple of tangent circles with
	// radii r, r1, r2
	public static double cosAngle(double r, double r1, double r2) {
		double c = r1 * r2;
		double denom = r * r + r * (r1 + r2) + c;
		if (denom <= 1e-16)
			return 1.0;
		double val = 1.0 - 2.0 * c / denom;
		if (val > 1.0)
			val = 1.0;
		if (val < -1.0)
			val = -1.0;
		return val;
	}

	// Equivalent to cosCorner.m, but we return angle directly to avoid clip/acos
	// repeated.
	// Angle at (x1,y1) in triangle (x1,y1)-(x2,y2)-(x3,y3)
	public static double angleAtCorner(double x1, double y1, double x2, double y2, double x3, double y3) {
		double l2 = Math.hypot(x2 - x1, y2 - y1);
		double l3 = Math.hypot(x3 - x1, y3 - y1);
		double l23 = Math.hypot(x3 - x2, y3 - y2);
		double denom = 2.0 * l2 * l3;
		if (denom <= 1e-16)
			return 0.0;
		double cosang = (l2 * l2 + l3 * l3 - l23 * l23) / denom;
		cosang = Math.max(-1.0, Math.min(1.0, cosang));
		return Math.acos(cosang);
	}
}