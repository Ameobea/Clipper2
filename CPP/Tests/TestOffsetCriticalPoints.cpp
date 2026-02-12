#include <gtest/gtest.h>
#include "clipper2/clipper.offset.h"
#include <cmath>

using namespace Clipper2Lib;

static Path64 MakeRect(int64_t left, int64_t bottom, int64_t right, int64_t top)
{
	Path64 path;
	path.emplace_back(left, bottom);
	path.emplace_back(right, bottom);	
	path.emplace_back(right, top);
	path.emplace_back(left, top);
	return path;
}

static void ExpectValidCriticals(const std::vector<double>& tvals)
{
	double prev = -1.0;
	for (double t : tvals)
	{
		EXPECT_GE(t, 0.0);
		EXPECT_LE(t, 1.0);
		EXPECT_GE(t + 1e-12, prev);
		prev = t;
	}
}

static double ComputeOptimalScale(double maxAbsCoord, double absDelta)
{
	const double safeMaxCoord = 1.0e15;
	const double maxScale = 1.0e15;
	const double minScale = 1.0;

	double maxExpectedCoord = maxAbsCoord + absDelta;
	if (maxExpectedCoord < 1.0e-10) return maxScale;

	double scale = safeMaxCoord / maxExpectedCoord;
	if (scale > maxScale) scale = maxScale;
	if (scale < minScale) scale = minScale;
	return scale;
}

static size_t RunSuperellipseCriticalCount(double angle_threshold)
{
	ClipperOffset co;
	Paths64 outputs;

	const double delta = 0.101;
	const double arc_tol = 0.000019999999494757503;
	const double maxAbsCoord = 1.0;
	const double scale = ComputeOptimalScale(maxAbsCoord, std::abs(delta));

	co.SetSuperellipseExponent(2.5);
	co.ArcTolerance(arc_tol * scale);
	co.SetChebyshevSpacing(true);
	co.SetAngleThreshold(angle_threshold);

	Path64 path_a;
	path_a.emplace_back(static_cast<int64_t>(std::llround(0.0 * scale)), static_cast<int64_t>(std::llround(-1.0 * scale)));
	path_a.emplace_back(static_cast<int64_t>(std::llround(0.0 * scale)), static_cast<int64_t>(std::llround(1.0 * scale)));

	Path64 path_b;
	path_b.emplace_back(static_cast<int64_t>(std::llround(-1.0 * scale)), static_cast<int64_t>(std::llround(0.0 * scale)));
	path_b.emplace_back(static_cast<int64_t>(std::llround(1.0 * scale)), static_cast<int64_t>(std::llround(0.0 * scale)));

	co.AddPath(path_a, JoinType::Superellipse, EndType::Superellipse);
	co.AddPath(path_b, JoinType::Superellipse, EndType::Superellipse);

	co.Execute(delta * scale, outputs);
	return co.GetCriticalTValues().size();
}

TEST(Clipper2Tests, OffsetCriticalPointsBasicSquare)
{
	ClipperOffset co;
	Paths64 outputs;
	Path64 square = MakeRect(0, 0, 10, 10);
	co.AddPath(square, JoinType::Miter, EndType::Polygon);
	co.Execute(2.0, outputs);

	ASSERT_EQ(outputs.size(), 1u);
	EXPECT_EQ(outputs[0].size(), 4u);

	const auto& critical = co.GetCriticalTValues();
	ExpectValidCriticals(critical);
	EXPECT_EQ(critical.size(), outputs[0].size());
}

TEST(Clipper2Tests, OffsetCriticalPointsUnionAddsCorners)
{
	ClipperOffset co;
	Paths64 outputs;

	Path64 rect_a = MakeRect(0, 0, 4, 2);
	Path64 rect_b = MakeRect(1, -1, 3, 3);
	co.AddPath(rect_a, JoinType::Miter, EndType::Polygon);
	co.AddPath(rect_b, JoinType::Miter, EndType::Polygon);

	co.Execute(1.0, outputs);

	ASSERT_EQ(outputs.size(), 1u);

	const auto& critical = co.GetCriticalTValues();
	ExpectValidCriticals(critical);
	EXPECT_GT(critical.size(), 8u);
}

TEST(Clipper2Tests, OffsetCriticalPointsSuperellipseJoin)
{
	ClipperOffset co;
	Paths64 outputs;

	Path64 square = MakeRect(0, 0, 1000, 1000);
	co.SetSuperellipseExponent(2.5);
	co.ArcTolerance(0.01);
	// With both angle and segment-fraction detection disabled, expect 0 critical points.
	// The superellipse smooths corners, so there are no sharp angles to detect.
	co.SetAngleThreshold(10.0); // radians (large value disables angle-based detection)
	co.SetCriticalSegmentFraction(1.0); // disable segment-fraction detection
	co.AddPath(square, JoinType::Superellipse, EndType::Polygon);
	co.Execute(100.0, outputs);

	ASSERT_EQ(outputs.size(), 1u);
	EXPECT_GT(outputs[0].size(), 4u);

	const auto& critical = co.GetCriticalTValues();
	ExpectValidCriticals(critical);
	// With smoothed corners and all detection disabled, expect no critical points
	EXPECT_EQ(critical.size(), 0u);
}

TEST(Clipper2Tests, OffsetCriticalPointsSegmentFractionDetection)
{
	// Test that segment-fraction detection finds corners adjacent to long segments
	ClipperOffset co;
	Paths64 outputs;

	// Create a tall thin rectangle where the long sides are a large fraction of perimeter
	Path64 rect = MakeRect(0, 0, 100, 1000);
	co.SetAngleThreshold(10.0); // disable angle-based detection
	co.SetCriticalSegmentFraction(0.15); // enable segment-fraction detection (default)
	co.AddPath(rect, JoinType::Miter, EndType::Polygon);
	co.Execute(10.0, outputs);

	ASSERT_EQ(outputs.size(), 1u);
	EXPECT_EQ(outputs[0].size(), 4u);

	const auto& critical = co.GetCriticalTValues();
	ExpectValidCriticals(critical);
	// All 4 corners have adjacent long segments (each long side is ~45% of perimeter)
	EXPECT_EQ(critical.size(), 4u);
}

TEST(Clipper2Tests, OffsetCriticalPointsSuperellipseBevelRepro)
{
	const size_t count_with_default = RunSuperellipseCriticalCount(0.3); // radians
	EXPECT_LE(count_with_default, 10u);

	const size_t count_without_angle = RunSuperellipseCriticalCount(10.0); // radians (disable angle test)
	EXPECT_LE(count_without_angle, 4u);
}

// Helper to compute bounding box dimensions
static void GetBoundingBox(const Paths64& paths, int64_t& width, int64_t& height)
{
	int64_t minX = INT64_MAX, minY = INT64_MAX;
	int64_t maxX = INT64_MIN, maxY = INT64_MIN;
	for (const auto& path : paths)
	{
		for (const auto& pt : path)
		{
			if (pt.x < minX) minX = pt.x;
			if (pt.x > maxX) maxX = pt.x;
			if (pt.y < minY) minY = pt.y;
			if (pt.y > maxY) maxY = pt.y;
		}
	}
	width = maxX - minX;
	height = maxY - minY;
}

// Test that superellipse end caps on angled lines are correctly oriented.
// Previously, the superellipse radius was computed using global coordinate angles,
// which caused the shape to be distorted on non-axis-aligned lines.
// After the fix, the superellipse shape should be consistent regardless of line orientation.
TEST(Clipper2Tests, OffsetSuperellipseAngledPolyline)
{
	// Test the user's specific case: polyline (0,0) -> (0,10) -> (10,1)
	// with superellipse join and end types
	const double scale = 1e6; // Use high precision
	const double delta = 0.2;
	const double superellipseExp = 10.0;

	ClipperOffset co;
	co.SetSuperellipseExponent(superellipseExp);
	co.ArcTolerance(1e-8 * scale);
	co.SetAngleThreshold(0.3);

	Path64 path;
	path.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(0.0 * scale));
	path.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(10.0 * scale));
	path.emplace_back(static_cast<int64_t>(10.0 * scale), static_cast<int64_t>(1.0 * scale));

	co.AddPath(path, JoinType::Superellipse, EndType::Superellipse);

	Paths64 outputs;
	co.Execute(delta * scale, outputs);

	// The offset should produce exactly one output path
	ASSERT_EQ(outputs.size(), 1u);

	// The output should have reasonable number of points
	// (superellipse with high exponent should have many points for the curved sections)
	EXPECT_GT(outputs[0].size(), 10u);

	// Verify output is valid (all points are at reasonable distance from original path)
	// Just check that the path has valid coordinates
	for (const auto& pt : outputs[0])
	{
		EXPECT_GE(pt.x, -1.0 * scale);
		EXPECT_LE(pt.x, 11.0 * scale);
		EXPECT_GE(pt.y, -1.0 * scale);
		EXPECT_LE(pt.y, 11.0 * scale);
	}
}

// Test that superellipse end caps have consistent shape regardless of line angle.
// We compare horizontal and 45-degree angled lines - the key property is that
// a 45-degree line should produce equal width and height (symmetric bounding box).
TEST(Clipper2Tests, OffsetSuperellipseShapeConsistency)
{
	const double scale = 1e6;
	const double delta = 1.0;
	const double lineLength = 10.0;
	const double superellipseExp = 4.0;

	// Test 1: Horizontal line segment - verify reasonable output
	{
		ClipperOffset co1;
		co1.SetSuperellipseExponent(superellipseExp);
		co1.ArcTolerance(1e-6 * scale);

		Path64 horizontal;
		horizontal.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(0.0 * scale));
		horizontal.emplace_back(static_cast<int64_t>(lineLength * scale), static_cast<int64_t>(0.0 * scale));

		co1.AddPath(horizontal, JoinType::Superellipse, EndType::Superellipse);

		Paths64 outputs1;
		co1.Execute(delta * scale, outputs1);
		ASSERT_EQ(outputs1.size(), 1u);

		int64_t width1, height1;
		GetBoundingBox(outputs1, width1, height1);

		// For a horizontal line, height should be approximately 2*D
		// (the superellipse may bulge slightly due to n > 2, but should be close)
		double expectedHeight = 2 * delta * scale;
		EXPECT_NEAR(static_cast<double>(height1), expectedHeight, 0.5 * scale);

		// Width should be at least lineLength + 2*delta
		double minExpectedWidth = (lineLength + 2 * delta) * scale;
		EXPECT_GE(static_cast<double>(width1), minExpectedWidth - 0.1 * scale);
	}

	// Test 2: 45-degree angled line segment
	// The key test: for a 45-degree line, the bounding box should be symmetric (width == height)
	// This verifies the superellipse shape is correctly oriented with the line direction,
	// not aligned with global coordinate axes.
	{
		ClipperOffset co2;
		co2.SetSuperellipseExponent(superellipseExp);
		co2.ArcTolerance(1e-6 * scale);

		// 45-degree line
		double diag = lineLength / std::sqrt(2.0);
		Path64 diagonal;
		diagonal.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(0.0 * scale));
		diagonal.emplace_back(static_cast<int64_t>(diag * scale), static_cast<int64_t>(diag * scale));

		co2.AddPath(diagonal, JoinType::Superellipse, EndType::Superellipse);

		Paths64 outputs2;
		co2.Execute(delta * scale, outputs2);
		ASSERT_EQ(outputs2.size(), 1u);

		int64_t width2, height2;
		GetBoundingBox(outputs2, width2, height2);

		// KEY TEST: For a 45-degree line, width and height must be equal
		// This would fail with the old code because the superellipse shape
		// was computed in global coordinates, causing asymmetry.
		EXPECT_NEAR(static_cast<double>(width2), static_cast<double>(height2), 0.01 * scale);

		// Both dimensions should be reasonable
		EXPECT_GT(width2, static_cast<int64_t>(diag * scale));
		EXPECT_GT(height2, static_cast<int64_t>(diag * scale));
	}
}

// Helper to find the maximum distance from a point along a given direction
static double GetMaxDistanceAlongDirection(const Paths64& paths, const Point64& origin, double dir_x, double dir_y)
{
	double max_dist = -1e20;
	for (const auto& path : paths)
	{
		for (const auto& pt : path)
		{
			double dx = static_cast<double>(pt.x - origin.x);
			double dy = static_cast<double>(pt.y - origin.y);
			double dist = dx * dir_x + dy * dir_y;
			if (dist > max_dist) max_dist = dist;
		}
	}
	return max_dist;
}

// Test superellipse join behavior based on standard superellipse math:
// - n=1: diamond shape → straight bevel line (apex at bevel distance)
// - n=2: circle → round arc (apex at delta distance)
// - n>2: squircle → fills toward miter (apex between round and miter)
// - n→∞: approaches square → fills entire miter triangle
//
// The affine transformation maps a unit superellipse |x|^n + |y|^n = 1 onto
// the corner's parallelogram. At the apex (t=π/4), the distance from the
// original vertex along the bisector is determined by the weight w = (1/√2)^(2/n).
TEST(Clipper2Tests, OffsetSuperellipseJoinProgression)
{
	const double scale = 1e6;
	const double delta = 1.0;

	// Create an L-shaped path with a 90-degree corner
	// The corner is at (0, 0), with arms going to (-10, 0) and (0, -10)
	// This gives a convex 90-degree outer corner when offset outward
	auto makeCornerPath = [scale]() {
		Path64 path;
		path.emplace_back(static_cast<int64_t>(-10.0 * scale), static_cast<int64_t>(0.0 * scale));
		path.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(0.0 * scale));
		path.emplace_back(static_cast<int64_t>(0.0 * scale), static_cast<int64_t>(-10.0 * scale));
		return path;
	};

	// For a 90-degree corner:
	// - Miter distance = delta / cos(45°) = delta * sqrt(2) ≈ 1.414 * delta
	// - Round distance = delta
	// - Bevel midpoint = delta * cos(45°) = delta / sqrt(2) ≈ 0.707 * delta
	double miter_dist = delta * std::sqrt(2.0);
	double round_dist = delta;
	double bevel_dist = delta / std::sqrt(2.0);

	// Bisector direction for the 90-degree corner at origin
	// The corner goes from (-1, 0) to (0, -1) direction
	// Bisector points outward at 45 degrees: (1, 1) / sqrt(2)
	double bisector_x = 1.0 / std::sqrt(2.0);
	double bisector_y = 1.0 / std::sqrt(2.0);
	Point64 origin(0, 0);

	// Test n=1 (diamond shape → bevel-like)
	double apex_dist_n1;
	{
		ClipperOffset co;
		co.SetSuperellipseExponent(1.0);
		co.ArcTolerance(1e-6 * scale);
		co.AddPath(makeCornerPath(), JoinType::Superellipse, EndType::Butt);

		Paths64 outputs;
		co.Execute(delta * scale, outputs);
		ASSERT_EQ(outputs.size(), 1u);

		apex_dist_n1 = GetMaxDistanceAlongDirection(outputs, origin, bisector_x, bisector_y) / scale;
		// n=1 should be close to bevel distance
		EXPECT_NEAR(apex_dist_n1, bevel_dist, 0.05);
	}

	// Test n=2 (circle → round)
	double apex_dist_n2;
	{
		ClipperOffset co;
		co.SetSuperellipseExponent(2.0);
		co.ArcTolerance(1e-6 * scale);
		co.AddPath(makeCornerPath(), JoinType::Superellipse, EndType::Butt);

		Paths64 outputs;
		co.Execute(delta * scale, outputs);
		ASSERT_EQ(outputs.size(), 1u);

		apex_dist_n2 = GetMaxDistanceAlongDirection(outputs, origin, bisector_x, bisector_y) / scale;
		// n=2 should be close to round distance
		EXPECT_NEAR(apex_dist_n2, round_dist, 0.05);
	}

	// Test n=10 (squircle → between round and miter)
	double apex_dist_n10;
	{
		ClipperOffset co;
		co.SetSuperellipseExponent(10.0);
		co.ArcTolerance(1e-6 * scale);
		co.AddPath(makeCornerPath(), JoinType::Superellipse, EndType::Butt);

		Paths64 outputs;
		co.Execute(delta * scale, outputs);
		ASSERT_EQ(outputs.size(), 1u);

		apex_dist_n10 = GetMaxDistanceAlongDirection(outputs, origin, bisector_x, bisector_y) / scale;
		// n=10 should be between round and miter, closer to miter
		EXPECT_GT(apex_dist_n10, round_dist);
		EXPECT_LT(apex_dist_n10, miter_dist);
	}

	// Verify monotonic progression: as n increases, apex moves toward miter
	// n=1 < n=2 < n=10
	EXPECT_LT(apex_dist_n1, apex_dist_n2);
	EXPECT_LT(apex_dist_n2, apex_dist_n10);
}
