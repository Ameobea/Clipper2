/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Date      :  22 January 2025                                                 *
* Website   :  https://www.angusj.com                                          *
* Copyright :  Angus Johnson 2010-2025                                         *
* Purpose   :  Path Offset (Inflate/Shrink)                                    *
* License   :  https://www.boost.org/LICENSE_1_0.txt                           *
*******************************************************************************/

#ifndef CLIPPER_OFFSET_H_
#define CLIPPER_OFFSET_H_

#include "clipper.core.h"
#include "clipper.engine.h"
#include <optional>

namespace Clipper2Lib {

enum class JoinType {
	Square,       // Joins are 'squared' at exactly the offset distance (more complex code)
	Bevel,        // Similar to Square, but offset distance varies with angle (simple & faster)
	Round,        // Rounded joins at offset points
	Miter,        // Sharp pointed joins (with miter limit consideration)
	Superellipse, // Squircle/notch shape controlled by exponent (2=circle, >2=squircle, <1=star)
	Knob,         // The "long way round" - 270 degree arc instead of 90
	Step,         // 90-degree staircase corner - two perpendicular segments
	Spike         // Sharp pointed protrusion extending outward from corner
};

enum class EndType {
	Polygon,      // Offsets only one side of a closed path
	Joined,       // Offsets both sides of a path, with joined ends
	Butt,         // Offsets both sides of a path, with square blunt ends
	Square,       // Offsets both sides of a path, with square extended ends
	Round,        // Offsets both sides of a path, with round extended ends
	// Extended end types
	Superellipse, // Superellipse end cap (squircle shape, controlled by exponent)
	Triangle,     // Triangular/pointed end cap
	Arrow,        // Arrow-shaped end with swept-back barbs
	Teardrop      // Pinched/teardrop round end
};

typedef std::function<double(const Path64& path, const PathD& path_normals, size_t curr_idx, size_t prev_idx)> DeltaCallback64;

class ClipperOffset {
private:

	class Group {
	public:
		Paths64 paths_in;
    std::optional<size_t> lowest_path_idx{};
		bool is_reversed = false;
		JoinType join_type;
		EndType end_type;
		Group(const Paths64& _paths, JoinType _join_type, EndType _end_type);
	};

	int   error_code_ = 0;
	double delta_ = 0.0;
	double group_delta_ = 0.0;
	double temp_lim_ = 0.0;
	double steps_per_rad_ = 0.0;
	double step_sin_ = 0.0;
	double step_cos_ = 0.0;
	PathD norms;
	Path64 path_out;
	Paths64* solution = nullptr;
	PolyTree64* solution_tree = nullptr;
	std::vector<Group> groups_;
	JoinType join_type_ = JoinType::Bevel;
	EndType end_type_ = EndType::Polygon;

	double miter_limit_ = 0.0;
	double arc_tolerance_ = 0.0;
	bool preserve_collinear_ = false;
	bool reverse_solution_ = false;

	// Extended configuration for new join/end types
	int custom_step_count_ = 0;           // 0 = use ArcTolerance/auto, >0 = fixed segment count
	double superellipse_exp_ = 2.5;       // Exponent: 2=circle, >2=squircle, <1=star/pinch
	double end_extension_scale_ = 1.0;    // Multiplier for end cap tip extension
	double arrow_back_sweep_ = 0.0;       // 0.0=flat triangle, >0=swept back barbs
	double teardrop_pinch_ = 0.5;         // 0.0=round, 1.0=sharp point

	// Angle threshold for applying special join types (in radians)
	// Joins with angles below this threshold use fallback_join_type_ instead
	double join_angle_threshold_ = 0.0;   // 0.0 = no threshold (always apply)
	JoinType fallback_join_type_ = JoinType::Bevel; // Used when angle < threshold
	bool use_chebyshev_spacing_ = false;  // Use Chebyshev nodes for better curvature sampling
	// Critical point detection (post-union)
	std::vector<double> critical_t_values_;
	double critical_angle_threshold_ = 0.3; // radians, ~17 degrees
	// Fraction of total perimeter - segments longer than this make adjacent vertices critical
	// 0 = use default (0.15), 1 = disabled (only angle threshold matters)
	double critical_segment_fraction_ = 0.0;
	// Simplification epsilon for output paths (0 = use default based on delta)
	double simplify_epsilon_ = 0.0;

	void DetectCriticalPoints(const Paths64& paths, double abs_delta);

#ifdef USINGZ
	ZCallback64 zCallback64_ = nullptr;
	void ZCB(const Point64& bot1, const Point64& top1,
		const Point64& bot2, const Point64& top2, Point64& ip);
#endif
	DeltaCallback64 deltaCallback64_ = nullptr;
	size_t CalcSolutionCapacity();
	bool CheckReverseOrientation();
	void DoBevel(const Path64& path, size_t j, size_t k);
	void DoSquare(const Path64& path, size_t j, size_t k);
	void DoMiter(const Path64& path, size_t j, size_t k, double cos_a);
	void DoRound(const Path64& path, size_t j, size_t k, double angle);
	// New join/end type helpers
	void DoSuperellipse(const Path64& path, size_t j, size_t k, double angle);
	void DoSuperellipseJoin(const Path64& path, size_t j, size_t k, double angle);
	void DoKnob(const Path64& path, size_t j, size_t k, double angle);
	void DoTriangle(const Path64& path, size_t j, size_t k);
	void DoArrow(const Path64& path, size_t j, size_t k);
	void DoTeardrop(const Path64& path, size_t j, size_t k, double angle);
	void DoStep(const Path64& path, size_t j, size_t k);
	void DoSpike(const Path64& path, size_t j, size_t k);
	int CalcSteps(double angle) const; // Helper to calculate steps considering custom_step_count_
	void BuildNormals(const Path64& path);
	void OffsetPolygon(Group& group, const Path64& path);
	void OffsetOpenJoined(Group& group, const Path64& path);
	void OffsetOpenPath(Group& group, const Path64& path);
	void OffsetPoint(Group& group, const Path64& path, size_t j, size_t k);
	void DoGroupOffset(Group &group);
	void ExecuteInternal(double delta);
public:
	explicit ClipperOffset(double miter_limit = 2.0,
		double arc_tolerance = 0.0,
		bool preserve_collinear = false,
		bool reverse_solution = false) :
		miter_limit_(miter_limit), arc_tolerance_(arc_tolerance),
		preserve_collinear_(preserve_collinear),
		reverse_solution_(reverse_solution) { };

	~ClipperOffset() { Clear(); };

	int ErrorCode() const { return error_code_; };
	void AddPath(const Path64& path, JoinType jt_, EndType et_);
	void AddPaths(const Paths64& paths, JoinType jt_, EndType et_);
	void Clear() { groups_.clear(); norms.clear(); critical_t_values_.clear(); };
	
	void Execute(double delta, Paths64& sols_64);
	void Execute(double delta, PolyTree64& polytree);
	void Execute(DeltaCallback64 delta_cb, Paths64& paths);

	double MiterLimit() const { return miter_limit_; }
	void MiterLimit(double miter_limit) { miter_limit_ = miter_limit; }

	//ArcTolerance: needed for rounded offsets (See offset_triginometry2.svg)
	double ArcTolerance() const { return arc_tolerance_; }
	void ArcTolerance(double arc_tolerance) { arc_tolerance_ = arc_tolerance; }

	bool PreserveCollinear() const { return preserve_collinear_; }
	void PreserveCollinear(bool preserve_collinear){preserve_collinear_ = preserve_collinear;}
	
	bool ReverseSolution() const { return reverse_solution_; }
	void ReverseSolution(bool reverse_solution) {reverse_solution_ = reverse_solution;}

#ifdef USINGZ
	void SetZCallback(ZCallback64 cb) { zCallback64_ = cb; }
#endif
	void SetDeltaCallback(DeltaCallback64 cb) { deltaCallback64_ = cb; }

	// Extended configuration setters for new join/end types
	// StepCount: 0 = auto (based on ArcTolerance), >0 = fixed segment count
	void SetStepCount(int steps) { custom_step_count_ = steps; }
	int StepCount() const { return custom_step_count_; }

	// SuperellipseExponent: 2.0 = circle, >2 = squircle, <1 = star/pinch
	void SetSuperellipseExponent(double exponent) { superellipse_exp_ = exponent; }
	double SuperellipseExponent() const { return superellipse_exp_; }

	// EndCapParams: scale controls tip extension, sweep controls arrow barb angle
	void SetEndCapParams(double scale, double sweep) {
		end_extension_scale_ = scale;
		arrow_back_sweep_ = sweep;
	}
	double EndExtensionScale() const { return end_extension_scale_; }
	double ArrowBackSweep() const { return arrow_back_sweep_; }

	// TeardropPinch: 0.0 = round, 1.0 = sharp point
	void SetTeardropPinch(double pinch) { teardrop_pinch_ = pinch; }
	double TeardropPinch() const { return teardrop_pinch_; }

	// JoinAngleThreshold: minimum angle (radians) at which special joins are applied
	// Angles smaller than this use fallback_join_type instead
	void SetJoinAngleThreshold(double threshold, JoinType fallback = JoinType::Bevel) {
		join_angle_threshold_ = threshold;
		fallback_join_type_ = fallback;
	}
	double JoinAngleThreshold() const { return join_angle_threshold_; }
	JoinType FallbackJoinType() const { return fallback_join_type_; }

	// ChebyshevSpacing: use Chebyshev nodes for better curvature sampling
	void SetChebyshevSpacing(bool use_chebyshev) { use_chebyshev_spacing_ = use_chebyshev; }
	bool ChebyshevSpacing() const { return use_chebyshev_spacing_; }

	// Critical point detection settings/accessors
	// Angle threshold is in radians.
	void SetAngleThreshold(double radians) { critical_angle_threshold_ = radians; }
	// Segment fraction: vertices adjacent to segments longer than this fraction of total
	// perimeter are marked critical (even with small angles). 0 = default (0.15), 1 = disabled.
	void SetCriticalSegmentFraction(double fraction) { critical_segment_fraction_ = fraction; }
	double CriticalSegmentFraction() const { return critical_segment_fraction_; }
	// Simplification epsilon for output paths. 0 = default (auto-scaled based on delta).
	void SetSimplifyEpsilon(double epsilon) { simplify_epsilon_ = epsilon; }
	double SimplifyEpsilon() const { return simplify_epsilon_; }
	const std::vector<double>& GetCriticalTValues() const { return critical_t_values_; }

};

}
#endif /* CLIPPER_OFFSET_H_ */
