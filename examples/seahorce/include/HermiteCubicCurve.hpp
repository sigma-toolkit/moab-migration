/*
* Downloaded from
* 	https://signalsmith-audio.co.uk/writing/2021/monotonic-smooth-interpolation/hermite-cubic-curve.h
*/
#ifndef SIGNALSMITH_MONOTONIC_CURVE_H
#define SIGNALSMITH_MONOTONIC_CURVE_H

#include <vector>
#include <array>
#include <limits>
#include <cmath>

template<typename Value>
class HermiteCubicSegment {
	Value xScale;
	// Cubic coefficients
	Value a, b, c, d;

public:
	Value x0, x1;

	HermiteCubicSegment(double x0, double x1, double y0, double y1, double g0, double g1) : x0(x0), x1(x1) {
		xScale = 1/(x1 - x0);

		double m0 = g0*(x1 - x0);
		double m1 = g1*(x1 - x0);
		a = 2*(y0 - y1) + m0 + m1;
		b = 3*(y1 - y0) - 2*m0 - m1;
		c = m0;
		d = y0;
	}
	double at(Value x) const {
		Value t = (x - x0)*xScale;
		return d + t*(c + t*(b + t*a));
	}
};

template<typename Value>
class HermiteCubicCurve {
	struct Point {
		Value x, y;
		Value step;
		bool peak;
		Value lineGradient, curveGradient;
	};
	Point earliest, latest;
	// used to collect points before assembling segments
	std::vector<Point> points;

public:
	enum class GradientMode {quadraticX, quadraticY, dualQuadratic, fritschCarlson};
	GradientMode mode = GradientMode::dualQuadratic;

	std::vector<HermiteCubicSegment<Value>> segments;

	bool monotonic = true;
	Value edges = 0; // By default, enforce flat edges - you can still pass the first/last point in twice to get an angled corner

	void clear() {
		points.resize(0);
		Point point{0, 0, 0, false, 0, 0};
		earliest = latest = point;
	}
	void add(Value x, Value y) {
		Point point{x, y, 0, false, 0, 0};
		if (!points.empty()) {
			Point &prev = points.back();
			prev.step = x - points.back().x;
			if (prev.step != 0) {
				prev.lineGradient = (y - prev.y)/prev.step;
			}
			// Check for peaks
			if (points.size() > 1) {
				Point &prev2 = points[points.size() - 2];
				if ((prev2.lineGradient < 0 && prev.lineGradient > 0) || (prev2.lineGradient > 0 && prev.lineGradient < 0)) {
					prev.peak = true;
				}
			}
		}
		points.push_back(point);
	}
	void finish() {
		if (points.empty()) return;
		earliest = points[0];
		latest = points.back();

		// Calculate curve gradients
		for (size_t i = 0; i < points.size(); ++i) {
			Point &current = points[i];
			if (monotonic && current.peak) {
				current.curveGradient = 0;
				continue;
			}
			if (i > 0 && points[i - 1].step != 0) { // We have a previous point
				Point &prev = points[i - 1];
				if (i + 1 < points.size() && points[i].step != 0) { // We have a next point
					// Different gradient modes
					if (mode == GradientMode::quadraticX) {
						current.curveGradient = (current.lineGradient*prev.step + prev.lineGradient*current.step)/(prev.step + current.step);
					} else if (mode == GradientMode::quadraticY) {
						Value weightPrev = std::abs(current.y - prev.y);
						Value weightNext = std::abs(points[i + 1].y - current.y);
						if (weightPrev + weightNext != 0) {
							current.curveGradient = (current.lineGradient*weightPrev + prev.lineGradient*weightNext)/(weightPrev + weightNext);
						} else {
							current.curveGradient = (current.lineGradient*prev.step + prev.lineGradient*current.step)/(prev.step + current.step);
						}
					} else if (mode == GradientMode::dualQuadratic) {
						Value weightPrev = std::sqrt(std::abs(prev.step*(current.y - prev.y)));
						Value weightNext = std::sqrt(std::abs(current.step*(points[i + 1].y - current.y)));
						if (weightPrev + weightNext > 0) {
							current.curveGradient = (current.lineGradient*weightPrev + prev.lineGradient*weightNext)/(weightPrev + weightNext);
						} else {
							current.curveGradient = (current.lineGradient*prev.step + prev.lineGradient*current.step)/(prev.step + current.step);
						}
					} else {
						current.curveGradient = (current.lineGradient + prev.lineGradient)/2;
					}
					if (monotonic && std::abs(current.curveGradient) > std::abs(current.lineGradient)*3) {
						current.curveGradient = current.lineGradient*3;
					}
				} else if (i > 1 && points[i - 2].step != 0) { // We have a previous-but-one point
					Point &prev2 = points[i - 2];
					Value prev2Line = prev.peak ? 0 : prev2.lineGradient;
					current.curveGradient = prev2Line + (prev.lineGradient - prev2Line)*1.5;
					// Don't cross 0
					if ((current.curveGradient > 0 && prev.lineGradient < 0) || (current.curveGradient < 0 && prev.lineGradient > 0)) {
						current.curveGradient = 0;
					}
				} else {
					current.curveGradient = prev.lineGradient;
				}
				if (monotonic && std::abs(current.curveGradient) > std::abs(prev.lineGradient)*3) {
					current.curveGradient = prev.lineGradient*3;
				}
			} else if (i + 2 < points.size() && points[i + 1].step != 0) { // We have a next-but-one point
				Point &next2 = points[i + 1];
				Value next2Line = next2.peak ? 0 : next2.lineGradient;
				current.curveGradient = next2Line + (current.lineGradient - next2Line)*1.5;
				// Don't cross 0
				if ((current.curveGradient > 0 && current.lineGradient < 0) || (current.curveGradient < 0 && current.lineGradient > 0)) {
					current.curveGradient = 0;
				}
			} else if (i + 1 < points.size() && points[i].step != 0) {
				current.curveGradient = current.lineGradient;
			}
		}
		if (!points.empty()) {
			points[0].curveGradient *= edges;
			points.back().curveGradient *= edges;
		}
		// Create segments
		segments.resize(0, segments[0]/* not actually used */);
		segments.reserve(points.size() - 1);
		for (size_t i = 1; i < points.size(); ++i) {
			Point &prev = points[i - 1];
			Point &next = points[i];
			if (prev.x != next.x) {
				segments.push_back({prev.x, next.x, prev.y, next.y, prev.curveGradient, next.curveGradient});
			}
		}
		points.resize(0);
	}

	Value at(Value x) const {
		if (x <= earliest.x) return earliest.y;
		if (x >= latest.x) return latest.y;
		for (int i = (int)segments.size() - 1; i >= 0; --i) {
			auto &segment = segments[i];
			if (x >= segment.x0) return segment.at(x);
		}
		return earliest.y;
	}
};

#endif