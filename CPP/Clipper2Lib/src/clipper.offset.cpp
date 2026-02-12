/*******************************************************************************
 * Author    :  Angus Johnson * Date      :  11 October 2025 * Website   :
 *https://www.angusj.com                                          * Copyright :
 *Angus Johnson 2010-2025                                         * Purpose   :
 *Path Offset (Inflate/Shrink)                                    * License   :
 *https://www.boost.org/LICENSE_1_0.txt                           *
 *******************************************************************************/

#include "clipper2/clipper.h"
#include "clipper2/clipper.offset.h"

namespace Clipper2Lib {

const double floating_point_tolerance = 1e-12;

// Clipper2 approximates arcs by using series of relatively short straight
// line segments. And logically, shorter line segments will produce better arc
// approximations. But very short segments can degrade performance, usually
// with little or no discernable improvement in curve quality. Very short
// segments can even detract from curve quality, due to the effects of integer
// rounding. Since there isn't an optimal number of line segments for any given
// arc radius (that perfectly balances curve approximation with performance),
// arc tolerance is user defined. Nevertheless, when the user doesn't define
// an arc tolerance (ie leaves alone the 0 default value), the calculated
// default arc tolerance (offset_radius / 500) generally produces good (smooth)
// arc approximations without producing excessively small segment lengths.
// See also: https://www.angusj.com/clipper2/Docs/Trigonometry.htm
const double arc_const = 0.002; // <-- 1/500

//------------------------------------------------------------------------------
// Miscellaneous methods
//------------------------------------------------------------------------------

void GetLowestClosedPathInfo(const Paths64 &paths, std::optional<size_t> &idx,
                             bool &is_neg_area) {
  idx.reset();
  Point64 botPt = Point64(INT64_MAX, INT64_MIN);
  for (size_t i = 0; i < paths.size(); ++i) {
    double a = MAX_DBL;
    for (const Point64 &pt : paths[i]) {
      if ((pt.y < botPt.y) || ((pt.y == botPt.y) && (pt.x >= botPt.x)))
        continue;
      if (a == MAX_DBL) {
        a = Area(paths[i]);
        if (a == 0)
          break; // invalid closed path, so break from inner loop
        is_neg_area = a < 0;
      }
      idx = i;
      botPt.x = pt.x;
      botPt.y = pt.y;
    }
  }
}

inline double Hypot(double x, double y) {
  // given that this is an internal function, and given the x and y parameters
  // will always be coordinate values (or the difference between coordinate
  // values), x and y should always be within INT64_MIN to INT64_MAX.
  // Consequently, there should be no risk that the following computation will
  // overflow see https://stackoverflow.com/a/32436148/359538
  return std::sqrt(x * x + y * y);
}

// Minimum turning angle (radians) below which we never mark a vertex as critical,
// even if it's adjacent to long segments. This prevents marking points on
// essentially-straight segments. ~2 degrees.
const double min_critical_angle = 0.035;

// Default segment fraction threshold when user passes 0
const double default_segment_fraction = 0.15;

void ClipperOffset::DetectCriticalPoints(const Paths64 &paths,
                                         double /*abs_delta*/) {
  critical_t_values_.clear();
  if (paths.empty())
    return;

  double angle_threshold = std::max(0.0, critical_angle_threshold_);

  // Resolve segment fraction: 0 = use default, 1 = disabled
  double seg_fraction = critical_segment_fraction_;
  if (seg_fraction <= 0.0)
    seg_fraction = default_segment_fraction;
  bool use_segment_check = (seg_fraction < 1.0);

  for (size_t p = 0; p < paths.size(); ++p) {
    const Path64 &path = paths[p];
    size_t path_size = path.size();
    if (path_size < 3)
      continue;

    // Compute segment lengths and total perimeter
    std::vector<double> segment_lengths;
    segment_lengths.reserve(path_size);

    double total_length = 0.0;
    for (size_t i = 0; i < path_size; ++i) {
      size_t next = (i + 1) % path_size;
      double dx = static_cast<double>(path[next].x - path[i].x);
      double dy = static_cast<double>(path[next].y - path[i].y);
      double len = std::sqrt(dx * dx + dy * dy);
      segment_lengths.push_back(len);
      total_length += len;
    }

    if (total_length < 1e-10)
      continue;

    // Compute cumulative lengths for t-value calculation
    std::vector<double> cumulative_lengths;
    cumulative_lengths.reserve(path_size);
    cumulative_lengths.push_back(0.0);
    for (size_t i = 1; i < path_size; ++i) {
      cumulative_lengths.push_back(cumulative_lengths[i - 1] +
                                   segment_lengths[i - 1]);
    }

    // Long segment threshold (absolute length)
    double long_segment_threshold = seg_fraction * total_length;

    for (size_t i = 0; i < path_size; ++i) {
      size_t prev = (i == 0) ? path_size - 1 : i - 1;
      size_t next = (i + 1) % path_size;

      // Compute turning angle at this vertex
      PointD v1(static_cast<double>(path[i].x - path[prev].x),
                static_cast<double>(path[i].y - path[prev].y));
      PointD v2(static_cast<double>(path[next].x - path[i].x),
                static_cast<double>(path[next].y - path[i].y));

      double len1 = std::sqrt(v1.x * v1.x + v1.y * v1.y);
      double len2 = std::sqrt(v2.x * v2.x + v2.y * v2.y);

      bool is_critical = false;

      if (len1 > 1e-10 && len2 > 1e-10) {
        double cos_angle = (v1.x * v2.x + v1.y * v2.y) / (len1 * len2);
        if (cos_angle > 1.0)
          cos_angle = 1.0;
        else if (cos_angle < -1.0)
          cos_angle = -1.0;
        double angle = std::acos(cos_angle);
        // angle is the turning angle between consecutive segments:
        // 0 = straight, PI = U-turn.

        // Primary check: angle exceeds user threshold
        if (angle > angle_threshold) {
          is_critical = true;
        }
        // Secondary check: angle exceeds minimum AND adjacent to long segment
        else if (use_segment_check && angle > min_critical_angle) {
          // Check if either adjacent segment is "long"
          // segment_lengths[prev] is the incoming segment (prev -> i)
          // segment_lengths[i] is the outgoing segment (i -> next)
          if (segment_lengths[prev] >= long_segment_threshold ||
              segment_lengths[i] >= long_segment_threshold) {
            is_critical = true;
          }
        }
      }

      if (is_critical)
        critical_t_values_.push_back(cumulative_lengths[i] / total_length);
    }
  }

  std::sort(critical_t_values_.begin(), critical_t_values_.end());
  critical_t_values_.erase(
      std::unique(critical_t_values_.begin(), critical_t_values_.end(),
                  [](double a, double b) { return std::abs(a - b) < 1e-9; }),
      critical_t_values_.end());
}

static PointD GetUnitNormal(const Point64 &pt1, const Point64 &pt2) {
  if (pt1 == pt2)
    return PointD(0.0, 0.0);
  double dx = static_cast<double>(pt2.x - pt1.x);
  double dy = static_cast<double>(pt2.y - pt1.y);
  double inverse_hypot = 1.0 / Hypot(dx, dy);
  dx *= inverse_hypot;
  dy *= inverse_hypot;
  return PointD(dy, -dx);
}

inline bool AlmostZero(double value, double epsilon = 0.001) {
  return std::fabs(value) < epsilon;
}

inline PointD NormalizeVector(const PointD &vec) {
  double h = Hypot(vec.x, vec.y);
  if (AlmostZero(h))
    return PointD(0, 0);
  double inverseHypot = 1 / h;
  return PointD(vec.x * inverseHypot, vec.y * inverseHypot);
}

inline PointD GetAvgUnitVector(const PointD &vec1, const PointD &vec2) {
  return NormalizeVector(PointD(vec1.x + vec2.x, vec1.y + vec2.y));
}

inline bool IsClosedPath(EndType et) {
  return et == EndType::Polygon || et == EndType::Joined;
}

static inline Point64 GetPerpendic(const Point64 &pt, const PointD &norm,
                                   double delta) {
#ifdef USINGZ
  return Point64(pt.x + norm.x * delta, pt.y + norm.y * delta, pt.z);
#else
  return Point64(pt.x + norm.x * delta, pt.y + norm.y * delta);
#endif
}

inline PointD GetPerpendicD(const Point64 &pt, const PointD &norm,
                            double delta) {
#ifdef USINGZ
  return PointD(pt.x + norm.x * delta, pt.y + norm.y * delta, pt.z);
#else
  return PointD(pt.x + norm.x * delta, pt.y + norm.y * delta);
#endif
}

inline void NegatePath(PathD &path) {
  for (PointD &pt : path) {
    pt.x = -pt.x;
    pt.y = -pt.y;
#ifdef USINGZ
    pt.z = pt.z;
#endif
  }
}

//------------------------------------------------------------------------------
// ClipperOffset::Group methods
//------------------------------------------------------------------------------

ClipperOffset::Group::Group(const Paths64 &_paths, JoinType _join_type,
                            EndType _end_type)
    : paths_in(_paths), join_type(_join_type), end_type(_end_type) {
  bool is_joined =
      (end_type == EndType::Polygon) || (end_type == EndType::Joined);
  for (Path64 &p : paths_in)
    StripDuplicates(p, is_joined);

  if (end_type == EndType::Polygon) {
    bool is_neg_area;
    GetLowestClosedPathInfo(paths_in, lowest_path_idx, is_neg_area);
    // the lowermost path must be an outer path, so if its orientation is
    // negative, then flag the whole group is 'reversed' (will negate delta
    // etc.) as this is much more efficient than reversing every path.
    is_reversed = lowest_path_idx.has_value() && is_neg_area;
  } else {
    lowest_path_idx.reset();
    is_reversed = false;
  }
}

//------------------------------------------------------------------------------
// ClipperOffset methods
//------------------------------------------------------------------------------

void ClipperOffset::AddPath(const Path64 &path, JoinType jt_, EndType et_) {
  groups_.emplace_back(Paths64(1, path), jt_, et_);
}

void ClipperOffset::AddPaths(const Paths64 &paths, JoinType jt_, EndType et_) {
  if (paths.size() == 0)
    return;
  groups_.emplace_back(paths, jt_, et_);
}

void ClipperOffset::BuildNormals(const Path64 &path) {
  norms.clear();
  norms.reserve(path.size());
  if (path.size() == 0)
    return;
  Path64::const_iterator path_iter, path_stop_iter = --path.cend();
  for (path_iter = path.cbegin(); path_iter != path_stop_iter; ++path_iter)
    norms.emplace_back(GetUnitNormal(*path_iter, *(path_iter + 1)));
  norms.emplace_back(GetUnitNormal(*path_stop_iter, *(path.cbegin())));
}

void ClipperOffset::DoBevel(const Path64 &path, size_t j, size_t k) {
  PointD pt1, pt2;
  if (j == k) {
    double abs_delta = std::abs(group_delta_);
#ifdef USINGZ
    pt1 = PointD(path[j].x - abs_delta * norms[j].x,
                 path[j].y - abs_delta * norms[j].y, path[j].z);
    pt2 = PointD(path[j].x + abs_delta * norms[j].x,
                 path[j].y + abs_delta * norms[j].y, path[j].z);
#else
    pt1 = PointD(path[j].x - abs_delta * norms[j].x,
                 path[j].y - abs_delta * norms[j].y);
    pt2 = PointD(path[j].x + abs_delta * norms[j].x,
                 path[j].y + abs_delta * norms[j].y);
#endif
  } else {
#ifdef USINGZ
    pt1 = PointD(path[j].x + group_delta_ * norms[k].x,
                 path[j].y + group_delta_ * norms[k].y, path[j].z);
    pt2 = PointD(path[j].x + group_delta_ * norms[j].x,
                 path[j].y + group_delta_ * norms[j].y, path[j].z);
#else
    pt1 = PointD(path[j].x + group_delta_ * norms[k].x,
                 path[j].y + group_delta_ * norms[k].y);
    pt2 = PointD(path[j].x + group_delta_ * norms[j].x,
                 path[j].y + group_delta_ * norms[j].y);
#endif
  }
  path_out.emplace_back(pt1);
  path_out.emplace_back(pt2);
}

void ClipperOffset::DoSquare(const Path64 &path, size_t j, size_t k) {
  PointD vec;
  if (j == k)
    vec = PointD(norms[j].y, -norms[j].x);
  else
    vec = GetAvgUnitVector(PointD(-norms[k].y, norms[k].x),
                           PointD(norms[j].y, -norms[j].x));

  double abs_delta = std::abs(group_delta_);

  // now offset the original vertex delta units along unit vector
  PointD ptQ = PointD(path[j]);
  ptQ = TranslatePoint(ptQ, abs_delta * vec.x, abs_delta * vec.y);
  // get perpendicular vertices
  PointD pt1 = TranslatePoint(ptQ, group_delta_ * vec.y, group_delta_ * -vec.x);
  PointD pt2 = TranslatePoint(ptQ, group_delta_ * -vec.y, group_delta_ * vec.x);
  // get 2 vertices along one edge offset
  PointD pt3 = GetPerpendicD(path[k], norms[k], group_delta_);
  if (j == k) {
    PointD pt4 =
        PointD(pt3.x + vec.x * group_delta_, pt3.y + vec.y * group_delta_);
    PointD pt = ptQ;
    GetLineIntersectPt(pt1, pt2, pt3, pt4, pt);
    // get the second intersect point through reflecion
    path_out.emplace_back(ReflectPoint(pt, ptQ));
    path_out.emplace_back(pt);
  } else {
    PointD pt4 = GetPerpendicD(path[j], norms[k], group_delta_);
    PointD pt = ptQ;
    GetLineIntersectPt(pt1, pt2, pt3, pt4, pt);
    path_out.emplace_back(pt);
    // get the second intersect point through reflecion
    path_out.emplace_back(ReflectPoint(pt, ptQ));
  }
}

void ClipperOffset::DoMiter(const Path64 &path, size_t j, size_t k,
                            double cos_a) {
  double q = group_delta_ / (cos_a + 1);
#ifdef USINGZ
  path_out.emplace_back(path[j].x + (norms[k].x + norms[j].x) * q,
                        path[j].y + (norms[k].y + norms[j].y) * q, path[j].z);
#else
  path_out.emplace_back(path[j].x + (norms[k].x + norms[j].x) * q,
                        path[j].y + (norms[k].y + norms[j].y) * q);
#endif
}

void ClipperOffset::DoRound(const Path64 &path, size_t j, size_t k,
                            double angle) {
  if (deltaCallback64_) {
    // when deltaCallback64_ is assigned, group_delta_ won't be constant,
    // so we'll need to do the following calculations for *every* vertex.
    double abs_delta = std::fabs(group_delta_);
    double arcTol = (arc_tolerance_ > floating_point_tolerance
                         ? std::min(abs_delta, arc_tolerance_)
                         : abs_delta * arc_const);
    double steps_per_360 =
        std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
    step_sin_ = std::sin(2 * PI / steps_per_360);
    step_cos_ = std::cos(2 * PI / steps_per_360);
    if (group_delta_ < 0.0)
      step_sin_ = -step_sin_;
    steps_per_rad_ = steps_per_360 / (2 * PI);
  }

  Point64 pt = path[j];
  PointD offsetVec =
      PointD(norms[k].x * group_delta_, norms[k].y * group_delta_);

  if (j == k)
    offsetVec.Negate();
#ifdef USINGZ
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y, pt.z);
#else
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y);
#endif
  int steps = CalcSteps(angle);
  for (int i = 1; i < steps; ++i) // ie 1 less than steps
  {
    offsetVec = PointD(offsetVec.x * step_cos_ - step_sin_ * offsetVec.y,
                       offsetVec.x * step_sin_ + offsetVec.y * step_cos_);
#ifdef USINGZ
    path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y, pt.z);
#else
    path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y);
#endif
  }
  path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_));
}

int ClipperOffset::CalcSteps(double angle) const {
  double abs_angle = std::abs(angle);
  if (custom_step_count_ > 0) {
    // Distribute user-defined steps proportionally over the angle
    return std::max(1, static_cast<int>(std::ceil(custom_step_count_ *
                                                  abs_angle / (2 * PI))));
  }
  // Use existing arc tolerance based calculation
  return static_cast<int>(std::ceil(steps_per_rad_ * abs_angle));
}

// Helper to compute interpolation parameter with optional Chebyshev spacing
static inline double GetInterpolationParam(int i, int steps,
                                           bool use_chebyshev) {
  if (use_chebyshev && steps > 1) {
    // Chebyshev nodes: t = 0.5 * (1 - cos(pi * i / steps))
    // This gives better spacing for high-curvature regions
    return 0.5 * (1.0 - std::cos(PI * static_cast<double>(i) /
                                 static_cast<double>(steps)));
  }
  // Linear (uniform) spacing
  return static_cast<double>(i) / static_cast<double>(steps);
}

void ClipperOffset::DoSuperellipse(const Path64 &path, size_t j, size_t k,
                                   double angle) {
  // Superellipse formula: |x|^n + |y|^n = 1
  // At angle θ, the radius is: r(θ) = 1 / (|cos θ|^n + |sin θ|^n)^(1/n)
  // n > 2: squircle (more square), n = 2: circle, n < 2: star/pinch
  //
  // IMPORTANT: The superellipse radius must be computed using LOCAL angles
  // (relative to the arc being traced), not global coordinate angles.
  // This ensures the superellipse shape is correctly oriented regardless
  // of the line's angle in the coordinate system.

  if (deltaCallback64_) {
    double abs_delta = std::fabs(group_delta_);
    double arcTol = (arc_tolerance_ > floating_point_tolerance
                         ? std::min(abs_delta, arc_tolerance_)
                         : abs_delta * arc_const);
    double steps_per_360 =
        std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
    steps_per_rad_ = steps_per_360 / (2 * PI);
  }

  Point64 pt = path[j];
  double n = superellipse_exp_;
  double abs_angle = std::fabs(angle);

  // Calculate the incoming normal direction (unit vector)
  PointD norm_in = (j == k) ? PointD(-norms[k].x, -norms[k].y) : norms[k];

  // Start point: offset along incoming normal
#ifdef USINGZ
  path_out.emplace_back(pt.x + norm_in.x * group_delta_,
                        pt.y + norm_in.y * group_delta_, pt.z);
#else
  path_out.emplace_back(pt.x + norm_in.x * group_delta_,
                        pt.y + norm_in.y * group_delta_);
#endif

  int steps = CalcSteps(angle);
  if (steps < 2)
    steps = 2;

  for (int i = 1; i < steps; ++i) {
    double t_param = GetInterpolationParam(i, steps, use_chebyshev_spacing_);

    // LOCAL angle for superellipse radius calculation:
    // Maps t from 0->1 to local angle from -abs_angle/2 to +abs_angle/2
    // At t=0.5 (middle of arc), local_angle = 0, which is the "tip" direction
    // where the superellipse should have r=1 (axis-aligned in superellipse
    // coords)
    double local_angle = (t_param - 0.5) * abs_angle;

    // Calculate superellipse radius using the LOCAL angle
    double abs_cos_local = std::fabs(std::cos(local_angle));
    double abs_sin_local = std::fabs(std::sin(local_angle));

    double r;
    if (abs_cos_local < floating_point_tolerance) {
      r = 1.0; // At ±90° local angle (perpendicular to bisector)
    } else if (abs_sin_local < floating_point_tolerance) {
      r = 1.0; // At 0° local angle (along bisector/tip)
    } else {
      r = 1.0 /
          std::pow(std::pow(abs_cos_local, n) + std::pow(abs_sin_local, n),
                   1.0 / n);
    }

    // Calculate actual direction by rotating norm_in through the arc
    double rotation_angle = t_param * angle;
    double cos_rot = std::cos(rotation_angle);
    double sin_rot = std::sin(rotation_angle);

    // Rotate the incoming normal to get the current direction
    PointD rotated_dir;
    rotated_dir.x = norm_in.x * cos_rot - norm_in.y * sin_rot;
    rotated_dir.y = norm_in.x * sin_rot + norm_in.y * cos_rot;

    // Apply superellipse radius and group_delta (which includes the offset
    // sign)
    double ox = rotated_dir.x * r * group_delta_;
    double oy = rotated_dir.y * r * group_delta_;

#ifdef USINGZ
    path_out.emplace_back(pt.x + ox, pt.y + oy, pt.z);
#else
    path_out.emplace_back(pt.x + ox, pt.y + oy);
#endif
  }

  // End point: offset along outgoing normal
  path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_));
}

void ClipperOffset::DoSuperellipseJoin(const Path64& path, size_t j, size_t k, double angle) {
  // 1. Setup Geometry
  // -----------------
  Point64 vert = path[j];
  PointD n_in = norms[k];  // Normal of incoming segment
  PointD n_out = norms[j]; // Normal of outgoing segment
  
  // Calculate the "Miter Point" (M) where the two offset lines would naturally intersect.
  // Standard vector algebra: M = V + (n_in + n_out) * (delta / (1 + dot(n_in, n_out)))
  double dot_prod = n_in.x * n_out.x + n_in.y * n_out.y;
  
  // Safety check: If lines are parallel (180 deg), dot_prod is -1.0. 
  // However, Clipper only calls Join for convex angles, so dot_prod > -1.0 usually.
  // We add a tiny epsilon to denom to prevent division by zero just in case.
  double denom = 1.0 + dot_prod;
  if (denom < 1e-6) denom = 1e-6; 

  double miter_scale = group_delta_ / denom;
  PointD vec_miter = PointD((n_in.x + n_out.x) * miter_scale, 
                            (n_in.y + n_out.y) * miter_scale);
  PointD M = PointD(vert.x + vec_miter.x, vert.y + vec_miter.y);

  // Calculate Start (P_in) and End (P_out) of the join
  PointD P_in  = PointD(vert.x + n_in.x * group_delta_, vert.y + n_in.y * group_delta_);
  PointD P_out = PointD(vert.x + n_out.x * group_delta_, vert.y + n_out.y * group_delta_);

  // 2. Define Affine Basis
  // ----------------------
  // We want to map a unit superellipse onto the parallelogram defined by the corner.
  // The parallelogram is defined by an "Anchor" (C) opposite to the Miter Tip (M).
  // C = P_in + P_out - M
  PointD C;
  C.x = P_in.x + P_out.x - M.x;
  C.y = P_in.y + P_out.y - M.y;

  // Basis Vector U: From Anchor to Start
  PointD U = PointD(P_in.x - C.x, P_in.y - C.y);
  // Basis Vector W: From Anchor to End
  PointD W = PointD(P_out.x - C.x, P_out.y - C.y);

  // 3. Generate Curve
  // -----------------
  // Note: We skip the first point (P_in) because Clipper expects the previous segment 
  // to have already added the final vertex. We add points 1..steps-1.
  // The final point (P_out) is usually added by the caller or the next segment logic,
  // but in DoRound, Clipper adds the final point explicitly. We will follow that pattern.

  path_out.emplace_back(P_in.x, P_in.y
#ifdef USINGZ
    , vert.z
#endif
  );

  int steps = CalcSteps(angle);
  // Ensure we have enough resolution for the curve
  if (steps < 2) steps = 2;

  double n = superellipse_exp_;
  
  // Pre-calculate 2/n for efficiency
  double exp_factor = 2.0 / n;

  for (int i = 1; i < steps; ++i) {
    // Parameter t goes from 0 to PI/2
    double t = (double)i / (double)steps * (PI * 0.5);

    double c_val = std::cos(t);
    double s_val = std::sin(t);

    // Apply Superellipse power (using absolute values for safety, 
    // though 0..PI/2 is positive)
    double weight_u = std::pow(std::abs(c_val), exp_factor);
    double weight_w = std::pow(std::abs(s_val), exp_factor);

    // Affine combination: Point = C + (weight_u * U) + (weight_w * W)
    double px = C.x + weight_u * U.x + weight_w * W.x;
    double py = C.y + weight_u * U.y + weight_w * W.y;

    path_out.emplace_back(px, py
#ifdef USINGZ
      , vert.z
#endif
    );
  }

  // Add the final point (P_out)
  path_out.emplace_back(P_out.x, P_out.y
#ifdef USINGZ
    , vert.z
#endif
  );
}

void ClipperOffset::DoKnob(const Path64 &path, size_t j, size_t k,
                           double angle) {
  // Knob: the "long way round" - use the reflex angle instead of the interior
  // angle If normal join would go 90 degrees, knob goes 270 degrees
  double reflex_angle;
  if (angle >= 0) {
    reflex_angle = -(2 * PI - angle);
  } else {
    reflex_angle = 2 * PI + angle;
  }

  // Use DoRound with the reflex angle, but we need to handle direction
  if (deltaCallback64_) {
    double abs_delta = std::fabs(group_delta_);
    double arcTol = (arc_tolerance_ > floating_point_tolerance
                         ? std::min(abs_delta, arc_tolerance_)
                         : abs_delta * arc_const);
    double steps_per_360 =
        std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
    step_sin_ = std::sin(2 * PI / steps_per_360);
    step_cos_ = std::cos(2 * PI / steps_per_360);
    if (group_delta_ < 0.0)
      step_sin_ = -step_sin_;
    steps_per_rad_ = steps_per_360 / (2 * PI);
  }

  Point64 pt = path[j];
  PointD offsetVec =
      PointD(norms[k].x * group_delta_, norms[k].y * group_delta_);

  if (j == k)
    offsetVec.Negate();

#ifdef USINGZ
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y, pt.z);
#else
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y);
#endif

  // For the knob, we go the "long way" around
  int steps = CalcSteps(reflex_angle);

  // Adjust step direction based on reflex angle sign
  double local_step_sin = step_sin_;
  double local_step_cos = step_cos_;
  if (reflex_angle < 0) {
    local_step_sin = -std::abs(step_sin_);
  } else {
    local_step_sin = std::abs(step_sin_);
  }
  if (group_delta_ < 0.0)
    local_step_sin = -local_step_sin;

  for (int i = 1; i < steps; ++i) {
    offsetVec =
        PointD(offsetVec.x * local_step_cos - local_step_sin * offsetVec.y,
               offsetVec.x * local_step_sin + offsetVec.y * local_step_cos);
#ifdef USINGZ
    path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y, pt.z);
#else
    path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y);
#endif
  }

  path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_));
}

void ClipperOffset::DoTriangle(const Path64 &path, size_t j, size_t k) {
  // Triangle end cap: pointed triangular tip
  // The tip extends outward by (delta * end_extension_scale_)

  Point64 pt = path[j];
  double abs_delta = std::fabs(group_delta_);

  // For end caps (j == k), we create a triangle pointing outward
  PointD norm = norms[j];
  if (j == k) {
    // At line endpoints, norm points perpendicular to the line
    // We need the direction along the line (tangent) for the tip
    // Tangent is 90 degrees from normal
    PointD tangent = PointD(-norm.y, norm.x);

    // Tip point extends along the tangent
    double tip_extension = abs_delta * end_extension_scale_;

    // Shoulder points are perpendicular to tangent (along original normal)
    PointD shoulder1 =
        PointD(pt.x + group_delta_ * norm.x, pt.y + group_delta_ * norm.y);
    PointD shoulder2 =
        PointD(pt.x - group_delta_ * norm.x, pt.y - group_delta_ * norm.y);
    PointD tip = PointD(pt.x + tip_extension * tangent.x,
                        pt.y + tip_extension * tangent.y);

#ifdef USINGZ
    path_out.emplace_back(shoulder1.x, shoulder1.y, pt.z);
    path_out.emplace_back(tip.x, tip.y, pt.z);
    path_out.emplace_back(shoulder2.x, shoulder2.y, pt.z);
#else
    path_out.emplace_back(shoulder1.x, shoulder1.y);
    path_out.emplace_back(tip.x, tip.y);
    path_out.emplace_back(shoulder2.x, shoulder2.y);
#endif
  } else {
    // For joins (not end caps), create a triangular join
    PointD pt1 = PointD(pt.x + group_delta_ * norms[k].x,
                        pt.y + group_delta_ * norms[k].y);
    PointD pt2 = PointD(pt.x + group_delta_ * norms[j].x,
                        pt.y + group_delta_ * norms[j].y);

    // Calculate tip at the angle bisector
    PointD bisector = GetAvgUnitVector(norms[k], norms[j]);
    double tip_dist = abs_delta * end_extension_scale_;
    PointD tip =
        PointD(pt.x + tip_dist * bisector.x, pt.y + tip_dist * bisector.y);

#ifdef USINGZ
    path_out.emplace_back(pt1.x, pt1.y, pt.z);
    path_out.emplace_back(tip.x, tip.y, pt.z);
    path_out.emplace_back(pt2.x, pt2.y, pt.z);
#else
    path_out.emplace_back(pt1.x, pt1.y);
    path_out.emplace_back(tip.x, tip.y);
    path_out.emplace_back(pt2.x, pt2.y);
#endif
  }
}

void ClipperOffset::DoArrow(const Path64 &path, size_t j, size_t k) {
  // Arrow end cap: V-shaped with swept-back barbs
  // Similar to triangle but the shoulders are swept back

  Point64 pt = path[j];
  double abs_delta = std::fabs(group_delta_);

  PointD norm = norms[j];
  if (j == k) {
    // At line endpoints
    PointD tangent = PointD(-norm.y, norm.x);

    double tip_extension = abs_delta * end_extension_scale_;
    double barb_sweep = abs_delta * arrow_back_sweep_;

    // Shoulder points are swept back from perpendicular
    PointD shoulder1 =
        PointD(pt.x + group_delta_ * norm.x - barb_sweep * tangent.x,
               pt.y + group_delta_ * norm.y - barb_sweep * tangent.y);
    PointD shoulder2 =
        PointD(pt.x - group_delta_ * norm.x - barb_sweep * tangent.x,
               pt.y - group_delta_ * norm.y - barb_sweep * tangent.y);
    PointD tip = PointD(pt.x + tip_extension * tangent.x,
                        pt.y + tip_extension * tangent.y);

#ifdef USINGZ
    path_out.emplace_back(shoulder1.x, shoulder1.y, pt.z);
    path_out.emplace_back(tip.x, tip.y, pt.z);
    path_out.emplace_back(shoulder2.x, shoulder2.y, pt.z);
#else
    path_out.emplace_back(shoulder1.x, shoulder1.y);
    path_out.emplace_back(tip.x, tip.y);
    path_out.emplace_back(shoulder2.x, shoulder2.y);
#endif
  } else {
    // For joins, use similar logic to triangle but with swept barbs
    PointD bisector = GetAvgUnitVector(norms[k], norms[j]);
    PointD tangent = PointD(-bisector.y, bisector.x);

    double tip_dist = abs_delta * end_extension_scale_;
    double barb_sweep = abs_delta * arrow_back_sweep_;

    PointD pt1 =
        PointD(pt.x + group_delta_ * norms[k].x - barb_sweep * tangent.x,
               pt.y + group_delta_ * norms[k].y - barb_sweep * tangent.y);
    PointD pt2 =
        PointD(pt.x + group_delta_ * norms[j].x + barb_sweep * tangent.x,
               pt.y + group_delta_ * norms[j].y + barb_sweep * tangent.y);
    PointD tip =
        PointD(pt.x + tip_dist * bisector.x, pt.y + tip_dist * bisector.y);

#ifdef USINGZ
    path_out.emplace_back(pt1.x, pt1.y, pt.z);
    path_out.emplace_back(tip.x, tip.y, pt.z);
    path_out.emplace_back(pt2.x, pt2.y, pt.z);
#else
    path_out.emplace_back(pt1.x, pt1.y);
    path_out.emplace_back(tip.x, tip.y);
    path_out.emplace_back(pt2.x, pt2.y);
#endif
  }
}

void ClipperOffset::DoTeardrop(const Path64 &path, size_t j, size_t k,
                               double angle) {
  // Teardrop: rounded end that pinches to a point
  // Combines round with a pointed tip based on teardrop_pinch_

  if (deltaCallback64_) {
    double abs_delta = std::fabs(group_delta_);
    double arcTol = (arc_tolerance_ > floating_point_tolerance
                         ? std::min(abs_delta, arc_tolerance_)
                         : abs_delta * arc_const);
    double steps_per_360 =
        std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
    step_sin_ = std::sin(2 * PI / steps_per_360);
    step_cos_ = std::cos(2 * PI / steps_per_360);
    if (group_delta_ < 0.0)
      step_sin_ = -step_sin_;
    steps_per_rad_ = steps_per_360 / (2 * PI);
  }

  Point64 pt = path[j];
  double pinch =
      std::max(0.0, std::min(1.0, teardrop_pinch_)); // Clamp to [0, 1]

  PointD offsetVec =
      PointD(norms[k].x * group_delta_, norms[k].y * group_delta_);
  if (j == k)
    offsetVec.Negate();

#ifdef USINGZ
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y, pt.z);
#else
  path_out.emplace_back(pt.x + offsetVec.x, pt.y + offsetVec.y);
#endif

  int steps = CalcSteps(angle);
  double abs_angle = std::abs(angle);

  // For teardrop, we modify the radius as we go around
  // At the middle (peak), the radius is reduced based on pinch
  // pinch=0: full radius (circle), pinch=1: radius goes to 0 (sharp point)
  for (int i = 1; i < steps; ++i) {
    double t = GetInterpolationParam(i, steps, use_chebyshev_spacing_);
    double current_angle = t * abs_angle;

    // Calculate pinch factor: peaks at the middle of the arc (t=0.5)
    // Using sin curve: maximum pinch at t=0.5
    double pinch_factor = std::sin(t * PI);
    double radius_scale = 1.0 - (pinch * pinch_factor);

    // Calculate position using rotation from start
    double cos_a = std::cos(current_angle);
    double sin_a = std::sin(current_angle);
    if (angle < 0)
      sin_a = -sin_a;

    // Rotate the starting offset vector
    PointD start_vec =
        PointD(norms[k].x * group_delta_, norms[k].y * group_delta_);
    if (j == k)
      start_vec.Negate();

    PointD rotated = PointD(start_vec.x * cos_a - start_vec.y * sin_a,
                            start_vec.x * sin_a + start_vec.y * cos_a);

    // Apply pinch (scale toward center)
    PointD pinched = PointD(rotated.x * radius_scale, rotated.y * radius_scale);

#ifdef USINGZ
    path_out.emplace_back(pt.x + pinched.x, pt.y + pinched.y, pt.z);
#else
    path_out.emplace_back(pt.x + pinched.x, pt.y + pinched.y);
#endif
  }

  path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_));
}

void ClipperOffset::DoStep(const Path64 &path, size_t j, size_t k) {
  // Step: 90-degree staircase corner
  // Creates two perpendicular segments instead of a diagonal bevel
  // Result looks like: __|  or |__ depending on turn direction

  Point64 pt = path[j];

  // Start point: offset along incoming normal
  PointD pt1 = PointD(pt.x + group_delta_ * norms[k].x,
                      pt.y + group_delta_ * norms[k].y);

  // End point: offset along outgoing normal
  PointD pt2 = PointD(pt.x + group_delta_ * norms[j].x,
                      pt.y + group_delta_ * norms[j].y);

  // Middle point: step corner
  // We go perpendicular from pt1 toward the bisector direction
  // The step is at the intersection of:
  // - A line through pt1, perpendicular to the incoming edge (parallel to
  // norms[k])
  // - A line through pt2, perpendicular to the outgoing edge (parallel to
  // norms[j]) For a proper step, we pick the corner that makes 90-degree turns

  // The step point is at (pt1.x, pt2.y) or (pt2.x, pt1.y) depending on geometry
  // But for arbitrary angles, we need to compute properly:
  // Step corner is where we extend from pt1 along the outgoing edge direction
  // or from pt2 along the incoming edge direction

  // Incoming edge direction (tangent): perpendicular to norms[k]
  PointD tangent_in = PointD(-norms[k].y, norms[k].x);
  // Outgoing edge direction (tangent): perpendicular to norms[j]
  PointD tangent_out = PointD(-norms[j].y, norms[j].x);

  // Step point: from pt1, move along tangent_in until we're aligned with pt2's
  // perpendicular This is essentially finding where a line from pt1 along
  // tangent_in intersects a line from pt2 along tangent_out (going backwards)
  PointD step_pt;
  if (!GetLineIntersectPt(
          pt1, PointD(pt1.x + tangent_in.x, pt1.y + tangent_in.y), pt2,
          PointD(pt2.x - tangent_out.x, pt2.y - tangent_out.y), step_pt)) {
    // Fallback if lines are parallel (shouldn't happen for typical corners)
    step_pt = PointD((pt1.x + pt2.x) / 2.0, (pt1.y + pt2.y) / 2.0);
  }

#ifdef USINGZ
  path_out.emplace_back(pt1.x, pt1.y, pt.z);
  path_out.emplace_back(step_pt.x, step_pt.y, pt.z);
  path_out.emplace_back(pt2.x, pt2.y, pt.z);
#else
  path_out.emplace_back(pt1.x, pt1.y);
  path_out.emplace_back(step_pt.x, step_pt.y);
  path_out.emplace_back(pt2.x, pt2.y);
#endif
}

void ClipperOffset::DoSpike(const Path64 &path, size_t j, size_t k) {
  // Spike: sharp pointed protrusion at corner
  // Like a thorn or star point extending outward along the angle bisector

  Point64 pt = path[j];
  double abs_delta = std::fabs(group_delta_);

  // Calculate the angle bisector direction
  PointD bisector = GetAvgUnitVector(norms[k], norms[j]);

  // Spike tip extends along bisector by (delta * end_extension_scale_)
  // Default scale of 1.0 means spike tip is at same distance as offset
  // Scale > 1.0 means longer spikes
  double spike_length = abs_delta * end_extension_scale_;

  PointD pt1 = PointD(pt.x + group_delta_ * norms[k].x,
                      pt.y + group_delta_ * norms[k].y);
  PointD pt2 = PointD(pt.x + group_delta_ * norms[j].x,
                      pt.y + group_delta_ * norms[j].y);

  // Spike tip
  PointD tip = PointD(pt.x + spike_length * bisector.x,
                      pt.y + spike_length * bisector.y);

#ifdef USINGZ
  path_out.emplace_back(pt1.x, pt1.y, pt.z);
  path_out.emplace_back(tip.x, tip.y, pt.z);
  path_out.emplace_back(pt2.x, pt2.y, pt.z);
#else
  path_out.emplace_back(pt1.x, pt1.y);
  path_out.emplace_back(tip.x, tip.y);
  path_out.emplace_back(pt2.x, pt2.y);
#endif
}

void ClipperOffset::OffsetPoint(Group &group, const Path64 &path, size_t j,
                                size_t k) {
  // Let A = change in angle where edges join
  // A == 0: ie no change in angle (flat join)
  // A == PI: edges 'spike'
  // sin(A) < 0: right turning
  // cos(A) < 0: change in angle is more than 90 degree

  if (path[j] == path[k])
    return;

  double sin_a = CrossProduct(norms[j], norms[k]);
  double cos_a = DotProduct(norms[j], norms[k]);
  if (sin_a > 1.0)
    sin_a = 1.0;
  else if (sin_a < -1.0)
    sin_a = -1.0;

  if (deltaCallback64_) {
    group_delta_ = deltaCallback64_(path, norms, j, k);
    if (group.is_reversed)
      group_delta_ = -group_delta_;
  }
  if (std::fabs(group_delta_) <= floating_point_tolerance) {
    path_out.emplace_back(path[j]);
    return;
  }

  if (cos_a > -0.999 &&
      (sin_a * group_delta_ < 0)) // test for concavity first (#593)
  {
    // is concave
    // by far the simplest way to construct concave joins, especially those
    // joining very short segments, is to insert 3 points that produce negative
    // regions. These regions will be removed later by the finishing union
    // operation. This is also the best way to ensure that path reversals (ie
    // over-shrunk paths) are removed.
#ifdef USINGZ
    path_out.emplace_back(GetPerpendic(path[j], norms[k], group_delta_),
                          path[j].z);
    path_out.emplace_back(path[j]); // (#405, #873, #916)
    path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_),
                          path[j].z);
#else
    path_out.emplace_back(GetPerpendic(path[j], norms[k], group_delta_));
    path_out.emplace_back(path[j]); // (#405, #873, #916)
    path_out.emplace_back(GetPerpendic(path[j], norms[j], group_delta_));
#endif
  } else if (cos_a > 0.999 && join_type_ != JoinType::Round &&
             join_type_ != JoinType::Superellipse &&
             join_type_ != JoinType::Knob && join_type_ != JoinType::Step &&
             join_type_ != JoinType::Spike) {
    // almost straight - less than 2.5 degree (#424, #482, #526 & #724)
    DoMiter(path, j, k, cos_a);
  } else {
    // Calculate the actual angle for threshold checking
    double angle = std::atan2(sin_a, cos_a);
    double abs_angle = std::abs(angle);

    // Check if angle is below threshold for special join types
    // If so, use fallback join type instead
    JoinType effective_join = join_type_;
    if (join_angle_threshold_ > 0.0 && abs_angle < join_angle_threshold_) {
      // Angle is too small for special join, use fallback
      if (join_type_ == JoinType::Superellipse ||
          join_type_ == JoinType::Knob || join_type_ == JoinType::Step ||
          join_type_ == JoinType::Spike) {
        effective_join = fallback_join_type_;
      }
    }

    switch (effective_join) {
    case JoinType::Miter:
      // miter unless the angle is sufficiently acute to exceed ML
      if (cos_a > temp_lim_ - 1)
        DoMiter(path, j, k, cos_a);
      else
        DoSquare(path, j, k);
      break;
    case JoinType::Round:
      DoRound(path, j, k, angle);
      break;
    case JoinType::Bevel:
      DoBevel(path, j, k);
      break;
    case JoinType::Superellipse:
      DoSuperellipseJoin(path, j, k, angle);
      break;
    case JoinType::Knob:
      DoKnob(path, j, k, angle);
      break;
    case JoinType::Step:
      DoStep(path, j, k);
      break;
    case JoinType::Spike:
      DoSpike(path, j, k);
      break;
    default: // JoinType::Square
      DoSquare(path, j, k);
      break;
    }
  }
}

void ClipperOffset::OffsetPolygon(Group &group, const Path64 &path) {
  path_out.clear();
  for (Path64::size_type j = 0, k = path.size() - 1; j < path.size();
       k = j, ++j)
    OffsetPoint(group, path, j, k);
  solution->emplace_back(path_out);
}

void ClipperOffset::OffsetOpenJoined(Group &group, const Path64 &path) {
  OffsetPolygon(group, path);
  Path64 reverse_path(path);
  std::reverse(reverse_path.begin(), reverse_path.end());

  // rebuild normals
  std::reverse(norms.begin(), norms.end());
  norms.emplace_back(norms[0]);
  norms.erase(norms.begin());
  NegatePath(norms);

  OffsetPolygon(group, reverse_path);
}

void ClipperOffset::OffsetOpenPath(Group &group, const Path64 &path) {
  // do the line start cap
  if (deltaCallback64_)
    group_delta_ = deltaCallback64_(path, norms, 0, 0);

  if (std::fabs(group_delta_) <= floating_point_tolerance)
    path_out.emplace_back(path[0]);
  else {
    switch (end_type_) {
    case EndType::Butt:
      DoBevel(path, 0, 0);
      break;
    case EndType::Round:
      DoRound(path, 0, 0, PI);
      break;
    case EndType::Superellipse:
      DoSuperellipse(path, 0, 0, PI);
      break;
    case EndType::Triangle:
      DoTriangle(path, 0, 0);
      break;
    case EndType::Arrow:
      DoArrow(path, 0, 0);
      break;
    case EndType::Teardrop:
      DoTeardrop(path, 0, 0, PI);
      break;
    default:
      DoSquare(path, 0, 0);
      break;
    }
  }

  size_t highI = path.size() - 1;
  // offset the left side going forward
  for (Path64::size_type j = 1, k = 0; j < highI; k = j, ++j)
    OffsetPoint(group, path, j, k);

  // reverse normals
  for (size_t i = highI; i > 0; --i)
    norms[i] = PointD(-norms[i - 1].x, -norms[i - 1].y);
  norms[0] = norms[highI];

  // do the line end cap
  if (deltaCallback64_)
    group_delta_ = deltaCallback64_(path, norms, highI, highI);

  if (std::fabs(group_delta_) <= floating_point_tolerance)
    path_out.emplace_back(path[highI]);
  else {
    switch (end_type_) {
    case EndType::Butt:
      DoBevel(path, highI, highI);
      break;
    case EndType::Round:
      DoRound(path, highI, highI, PI);
      break;
    case EndType::Superellipse:
      DoSuperellipse(path, highI, highI, PI);
      break;
    case EndType::Triangle:
      DoTriangle(path, highI, highI);
      break;
    case EndType::Arrow:
      DoArrow(path, highI, highI);
      break;
    case EndType::Teardrop:
      DoTeardrop(path, highI, highI, PI);
      break;
    default:
      DoSquare(path, highI, highI);
      break;
    }
  }

  for (size_t j = highI - 1, k = highI; j > 0; k = j, --j)
    OffsetPoint(group, path, j, k);
  solution->emplace_back(path_out);
}

void ClipperOffset::DoGroupOffset(Group &group) {
  if (group.end_type == EndType::Polygon) {
    // a straight path (2 points) can now also be 'polygon' offset
    // where the ends will be treated as (180 deg.) joins
    if (!group.lowest_path_idx.has_value())
      delta_ = std::abs(delta_);
    group_delta_ = (group.is_reversed) ? -delta_ : delta_;
  } else
    group_delta_ = std::abs(delta_); // *0.5;

  double abs_delta = std::fabs(group_delta_);
  join_type_ = group.join_type;
  end_type_ = group.end_type;

  // Check if we need arc step calculations for this group
  bool needs_arc_steps = group.join_type == JoinType::Round ||
                         group.join_type == JoinType::Superellipse ||
                         group.join_type == JoinType::Knob ||
                         group.end_type == EndType::Round ||
                         group.end_type == EndType::Superellipse ||
                         group.end_type == EndType::Teardrop;

  if (needs_arc_steps) {
    // calculate the number of steps required to approximate a circle
    // (see https://www.angusj.com/clipper2/Docs/Trigonometry.htm)
    // arcTol - when arc_tolerance_ is undefined (0) then curve imprecision
    // will be relative to the size of the offset (delta). Obviously very
    // large offsets will almost always require much less precision.
    double arcTol = (arc_tolerance_ > floating_point_tolerance)
                        ? std::min(abs_delta, arc_tolerance_)
                        : abs_delta * arc_const;

    double steps_per_360 =
        std::min(PI / std::acos(1 - arcTol / abs_delta), abs_delta * PI);
    step_sin_ = std::sin(2 * PI / steps_per_360);
    step_cos_ = std::cos(2 * PI / steps_per_360);
    if (group_delta_ < 0.0)
      step_sin_ = -step_sin_;
    steps_per_rad_ = steps_per_360 / (2 * PI);
  }

  // double min_area = PI * Sqr(group_delta_);
  Paths64::const_iterator path_in_it = group.paths_in.cbegin();
  for (; path_in_it != group.paths_in.cend(); ++path_in_it) {
    Path64::size_type pathLen = path_in_it->size();
    path_out.clear();

    if (pathLen == 1) // single point
    {
      if (deltaCallback64_) {
        group_delta_ = deltaCallback64_(*path_in_it, norms, 0, 0);
        if (group.is_reversed)
          group_delta_ = -group_delta_;
        abs_delta = std::fabs(group_delta_);
      }

      if (group_delta_ < 1)
        continue;
      const Point64 &pt = (*path_in_it)[0];
      // single vertex so build a circle, superellipse, or square ...
      if (group.join_type == JoinType::Round ||
          group.join_type == JoinType::Knob) {
        double radius = abs_delta;
        size_t steps =
            steps_per_rad_ > 0
                ? static_cast<size_t>(std::ceil(steps_per_rad_ * 2 * PI))
                : 0; // #617
        path_out = Ellipse(pt, radius, radius, steps);
#ifdef USINGZ
        for (auto &p : path_out)
          p.z = pt.z;
#endif
      } else if (group.join_type == JoinType::Superellipse) {
        // Generate superellipse shape for single point
        double radius = abs_delta;
        size_t steps =
            steps_per_rad_ > 0
                ? static_cast<size_t>(std::ceil(steps_per_rad_ * 2 * PI))
                : 0;
        if (steps < 8)
          steps = 8;
        double power = 2.0 / superellipse_exp_;
        path_out.clear();
        for (size_t i = 0; i < steps; ++i) {
          double t =
              2 * PI * static_cast<double>(i) / static_cast<double>(steps);
          double cos_t = std::cos(t);
          double sin_t = std::sin(t);
          // Superellipse: x = sgn(cos)*|cos|^(2/n), y = sgn(sin)*|sin|^(2/n)
          double sx = (cos_t >= 0 ? 1 : -1) * std::pow(std::abs(cos_t), power);
          double sy = (sin_t >= 0 ? 1 : -1) * std::pow(std::abs(sin_t), power);
#ifdef USINGZ
          path_out.emplace_back(pt.x + radius * sx, pt.y + radius * sy, pt.z);
#else
          path_out.emplace_back(pt.x + radius * sx, pt.y + radius * sy);
#endif
        }
      } else {
        int d = (int)std::ceil(abs_delta);
        Rect64 r = Rect64(pt.x - d, pt.y - d, pt.x + d, pt.y + d);
        path_out = r.AsPath();
#ifdef USINGZ
        for (auto &p : path_out)
          p.z = pt.z;
#endif
      }

      solution->emplace_back(path_out);
      continue;
    } // end of offsetting a single point

    if ((pathLen == 2) && (group.end_type == EndType::Joined)) {
      // Map join types to corresponding end types for 2-point paths
      switch (group.join_type) {
      case JoinType::Round:
      case JoinType::Knob:
        end_type_ = EndType::Round;
        break;
      case JoinType::Superellipse:
        end_type_ = EndType::Superellipse;
        break;
      default:
        end_type_ = EndType::Square;
        break;
      }
    }

    BuildNormals(*path_in_it);
    if (end_type_ == EndType::Polygon)
      OffsetPolygon(group, *path_in_it);
    else if (end_type_ == EndType::Joined)
      OffsetOpenJoined(group, *path_in_it);
    else
      OffsetOpenPath(group, *path_in_it);
  }
}

#ifdef USINGZ
void ClipperOffset::ZCB(const Point64 &bot1, const Point64 &top1,
                        const Point64 &bot2, const Point64 &top2, Point64 &ip) {
  if (bot1.z && ((bot1.z == bot2.z) || (bot1.z == top2.z)))
    ip.z = bot1.z;
  else if (bot2.z && (bot2.z == top1.z))
    ip.z = bot2.z;
  else if (top1.z && (top1.z == top2.z))
    ip.z = top1.z;
  else if (zCallback64_)
    zCallback64_(bot1, top1, bot2, top2, ip);
}
#endif

size_t ClipperOffset::CalcSolutionCapacity() {
  size_t result = 0;
  for (const Group &g : groups_)
    result += (g.end_type == EndType::Joined) ? g.paths_in.size() * 2
                                              : g.paths_in.size();
  return result;
}

bool ClipperOffset::CheckReverseOrientation() {
  // nb: this assumes there's consistency in orientation between groups
  bool is_reversed_orientation = false;
  for (const Group &g : groups_)
    if (g.end_type == EndType::Polygon) {
      is_reversed_orientation = g.is_reversed;
      break;
    }
  return is_reversed_orientation;
}

void ClipperOffset::ExecuteInternal(double delta) {
  error_code_ = 0;
  critical_t_values_.clear();
  if (groups_.size() == 0)
    return;
  solution->reserve(CalcSolutionCapacity());

  if (std::abs(delta) < 0.5) // ie: offset is insignificant
  {
    Paths64::size_type sol_size = 0;
    for (const Group &group : groups_)
      sol_size += group.paths_in.size();
    solution->reserve(sol_size);
    for (const Group &group : groups_)
      copy(group.paths_in.begin(), group.paths_in.end(),
           back_inserter(*solution));
  } else {

    temp_lim_ = (miter_limit_ <= 1) ? 2.0 : 2.0 / (miter_limit_ * miter_limit_);

    delta_ = delta;
    std::vector<Group>::iterator git;
    for (git = groups_.begin(); git != groups_.end(); ++git) {
      DoGroupOffset(*git);
      if (!error_code_)
        continue; // all OK
      solution->clear();
    }
  }

  if (!solution->size())
    return;

  bool paths_reversed = CheckReverseOrientation();
  // clean up self-intersections ...
  Clipper64 c;
  c.PreserveCollinear(preserve_collinear_);
  // the solution should retain the orientation of the input
  c.ReverseSolution(reverse_solution_ != paths_reversed);
#ifdef USINGZ
  auto fp = std::bind(&ClipperOffset::ZCB, this, std::placeholders::_1,
                      std::placeholders::_2, std::placeholders::_3,
                      std::placeholders::_4, std::placeholders::_5);
  c.SetZCallback(fp);
#endif
  c.AddSubject(*solution);
  if (solution_tree) {
    if (paths_reversed)
      c.Execute(ClipType::Union, FillRule::Negative, *solution_tree);
    else
      c.Execute(ClipType::Union, FillRule::Positive, *solution_tree);
  } else {
    if (paths_reversed)
      c.Execute(ClipType::Union, FillRule::Negative, *solution);
    else
      c.Execute(ClipType::Union, FillRule::Positive, *solution);
  }

  // Simplify output paths to collapse nearly-collinear segments before
  // critical point detection. Uses a very strict tolerance to preserve
  // fine topology while only removing truly collinear points from the
  // union operation.
  double simplify_eps = simplify_epsilon_;
  if (simplify_eps <= 0.0) {
    // Default: use a very small absolute value to only collapse
    // truly collinear points without removing any meaningful geometry.
    // This is intentionally very conservative - it should only remove
    // points that lie exactly (within floating point tolerance) on a line.
    simplify_eps = 0.5; // Half a unit in integer coordinate space
  }

  Paths64 paths_for_critical;
  if (solution_tree) {
    paths_for_critical =
        SimplifyPaths(PolyTreeToPaths64(*solution_tree), simplify_eps, true);
  } else {
    *solution = SimplifyPaths(*solution, simplify_eps, true);
    paths_for_critical = *solution;
  }

  DetectCriticalPoints(paths_for_critical, std::abs(delta));
}

void ClipperOffset::Execute(double delta, Paths64 &paths64) {
  paths64.clear();
  solution = &paths64;
  solution_tree = nullptr;
  ExecuteInternal(delta);
}

void ClipperOffset::Execute(double delta, PolyTree64 &polytree) {
  polytree.Clear();
  solution_tree = &polytree;
  solution = new Paths64();
  ExecuteInternal(delta);
  delete solution;
  solution = nullptr;
}

void ClipperOffset::Execute(DeltaCallback64 delta_cb, Paths64 &paths) {
  deltaCallback64_ = delta_cb;
  Execute(1.0, paths);
}

} // namespace Clipper2Lib
