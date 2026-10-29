#include "ActuatorEffectivenessDeltaRay.hpp"

#include <lib/mathlib/mathlib.h>
#include <px4_platform_common/log.h>
#include <cmath>
#include <cfloat>

using namespace matrix;

ActuatorEffectivenessDeltaRay::ActuatorEffectivenessDeltaRay(ModuleParams *parent)
    : ModuleParams(parent),
      _rotors(this, ActuatorEffectivenessRotors::AxisConfiguration::Configurable, true),
      _tilts(this)
{
}

bool ActuatorEffectivenessDeltaRay::validGeometry() const
{
    const auto &g = _rotors.geometry();

    if (g.num_rotors != N || _tilts.count() != N) { return false; }

    for (int i = 0; i < N; ++i) {
        const auto &r = g.rotors[i];
        const int k = r.tilt_index;

        // This example requires motor i to be paired with servo i.
        if (k != i || !PX4_ISFINITE(r.thrust_coef) || r.thrust_coef <= FLT_EPSILON ||
            !PX4_ISFINITE(r.moment_ratio)) { return false; }

        for (int j = 0; j < 3; ++j) {
            if (!PX4_ISFINITE(r.position(j))) { return false; }
        }

        const auto &t = _tilts.config(k);
        if (!PX4_ISFINITE(t.min_angle) || !PX4_ISFINITE(t.max_angle) ||
            t.min_angle < -PI_HALF - 1e-4f || t.max_angle > PI_HALF + 1e-4f ||
            t.max_angle - t.min_angle <= FLT_EPSILON ||
            t.tilt_direction != ActuatorEffectivenessTilts::TiltDirection::TowardsFront) {
            return false;
        }
    }

    return true;
}

bool ActuatorEffectivenessDeltaRay::getEffectivenessMatrix(
    Configuration &configuration, EffectivenessUpdateReason reason)
{
    if (reason == EffectivenessUpdateReason::NO_EXTERNAL_UPDATE) { return false; }

    _ready = false;
    if (!validGeometry()) {
        PX4_ERR("DeltaRay: need 2 fixed EDFs + 2 independent XZ nozzles");
        return false;
    }

    // Register physical motor/servo types in PX4's normal publishing order.
    // Initially built columns are overwritten with virtual-force columns below.
    if (!_rotors.addActuators(configuration)) { return false; }
    if (!_tilts.addActuators(configuration)) { return false; }

    // Matrix order = [h0,h1,v0,v1], with positive X reaction to rearward exhaust.
    // Positive alpha rotates the exhaust DOWN (FRD +Z), hence force UP (-Z).
    // MOTOR slots stand for h, SERVO slots stand for v *only during pseudoinverse*.
    // The actual physical commands replace all four slots in updateSetpoint().
    auto &B = configuration.effectiveness_matrices[0];
    B.setZero();
    const auto &g = _rotors.geometry();

    for (int i = 0; i < N; ++i) {
        const auto &r = g.rotors[i];
        const Vector3f fh{r.thrust_coef, 0.f, 0.f};
        const Vector3f fv{0.f, 0.f, -r.thrust_coef};
        const Vector3f mh = r.position.cross(fh) - r.moment_ratio * fh;
        const Vector3f mv = r.position.cross(fv) - r.moment_ratio * fv;

        for (int axis = 0; axis < 3; ++axis) {
            B(axis, i) = mh(axis);
            B(axis + 3, i) = fh(axis);
            B(axis, N + i) = mv(axis);
            B(axis + 3, N + i) = fv(axis);
        }
    }

    // No per-axis stop masks: a tiltable rotor cannot be labeled as 'upward-only'.
    _motor_direction_bitmasks = {};
    _ready = true;
    return true;
}

void ActuatorEffectivenessDeltaRay::projectSector(
    float &h, float &v, float amin, float amax, float radius) const
{
    if (!PX4_ISFINITE(h) || !PX4_ISFINITE(v) || radius <= 0.f || amin > amax) {
        h = 0.f; v = 0.f; return;
    }

    const float h0 = h, v0 = v;
    const float norm = hypotf(h0, v0);
    if (norm < 1e-7f) { h = 0.f; v = 0.f; return; }

    const float angle = atan2f(v0, h0);
    if (angle >= amin && angle <= amax) {
        const float scale = math::min(1.f, radius / norm);
        h *= scale; v *= scale; return;
    }

    // Orthogonal projection to the two sector edges or the origin.
    float best_dist2 = h0 * h0 + v0 * v0;
    float best_h = 0.f, best_v = 0.f;
    const float edges[2] = {amin, amax};

    for (float a : edges) {
        const float c = cosf(a), s = sinf(a);
        const float t = math::constrain(h0 * c + v0 * s, 0.f, radius);
        const float hc = t * c, vc = t * s;
        const float d2 = (hc - h0) * (hc - h0) + (vc - v0) * (vc - v0);

        if (d2 < best_dist2) {
            best_dist2 = d2; best_h = hc; best_v = vc;
        }
    }

    h = best_h; v = best_v;
}

void ActuatorEffectivenessDeltaRay::updateSetpoint(
    const Vector<float, NUM_AXES> &control_sp, int matrix_index,
    ActuatorVector &actuator_sp, const ActuatorVector &actuator_min,
    const ActuatorVector &actuator_max)
{
    if (matrix_index != 0) { return; }

    _desired = control_sp;
    _desired.nanToZero();
    _have_request = true;

    // Critical: save UNCLIPPED virtual q, then replace it with physical outputs.
    const float h[N] = {actuator_sp(0), actuator_sp(1)};
    const float v[N] = {actuator_sp(2), actuator_sp(3)};
    actuator_sp.setZero();
    if (!_ready) { return; }

    const hrt_abstime now = hrt_absolute_time();
    // First update has no previous output: allow full movement to target.
    const float dt = _last_time == 0 ? 1.f :
                     math::constrain((now - _last_time) * 1e-6f, 0.f, 0.05f);
    _last_time = now;

    for (int i = 0; i < N; ++i) {
        const auto &t = _tilts.config(i);
        const float range = t.max_angle - t.min_angle;
        // Physical motor and servo setpoint limits.
        const float max_motor = math::constrain(actuator_max(i), 0.f, 1.f);
        const float min_servo = math::constrain(actuator_min(N + i), -1.f, 1.f);
        const float max_servo = math::constrain(actuator_max(N + i), -1.f, 1.f);
        const float low_a = t.min_angle + 0.5f * (min_servo + 1.f) * range;
        const float high_a = t.min_angle + 0.5f * (max_servo + 1.f) * range;
        const float amin = math::max(low_a, _alpha[i] - MAX_TILT_RATE * dt);
        const float amax = math::min(high_a, _alpha[i] + MAX_TILT_RATE * dt);

        float hi = h[i], vi = v[i];
        projectSector(hi, vi, amin, amax, max_motor * max_motor);

        // q represents thrust fraction, not motor speed. Bridge->GZ speed is
        // linear; GZ thrust is proportional to speed squared.
        const float requested_T = hypotf(hi, vi);
        const float final_T = math::constrain(requested_T,
                                              math::max(0.f, _thrust[i] - MAX_THRUST_RATE * dt),
                                              math::min(max_motor * max_motor, _thrust[i] + MAX_THRUST_RATE * dt));
        const float final_a = (requested_T > 1e-6f) ? atan2f(vi, hi) :
                              math::constrain(_alpha[i], amin, amax);

        _thrust[i] = final_T;
        _alpha[i] = final_a;
        actuator_sp(i) = sqrtf(final_T);
        actuator_sp(N + i) = math::constrain(
            -1.f + 2.f * (final_a - t.min_angle) / range, -1.f, 1.f);
    }
}

matrix::Vector<float, ActuatorEffectiveness::NUM_AXES>
ActuatorEffectivenessDeltaRay::actualWrench() const
{
    Vector<float, NUM_AXES> wrench{};
    const auto &g = _rotors.geometry();

    for (int i = 0; i < N; ++i) {
        const float T = _thrust[i] * g.rotors[i].thrust_coef;
        const Vector3f force{T * cosf(_alpha[i]), 0.f, -T * sinf(_alpha[i])};
        const Vector3f moment = g.rotors[i].position.cross(force) -
                                g.rotors[i].moment_ratio * force;

        for (int j = 0; j < 3; ++j) {
            wrench(j) += moment(j);
            wrench(j + 3) += force(j);
        }
    }

    return wrench;
}

void ActuatorEffectivenessDeltaRay::getUnallocatedControl(
    int matrix_index, control_allocator_status_s &status)
{
    if (matrix_index != 0 || !_have_request) { return; }
    const Vector<float, NUM_AXES> residual = _desired - actualWrench();
    for (int i = 0; i < 3; ++i) {
        status.unallocated_torque[i] = residual(i);
        status.unallocated_thrust[i] = residual(i + 3);
    }
}
