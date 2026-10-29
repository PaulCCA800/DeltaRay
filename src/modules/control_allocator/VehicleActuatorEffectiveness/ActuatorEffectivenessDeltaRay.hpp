#pragma once

#include "control_allocation/actuator_effectiveness/ActuatorEffectiveness.hpp"
#include "ActuatorEffectivenessRotors.hpp"
#include "ActuatorEffectivenessTilts.hpp"
#include <px4_platform_common/module_params.h>
#include <drivers/drv_hrt.h>

// Two FIXED EDF motors with independently deflected exhaust nozzles.
// NO elevon, aileron, elevator or flap actuation: 2 motors + 2 nozzle servos only.
// alpha=0: exhaust -X and reaction thrust +X.
// alpha=+pi/2: exhaust +Z (down FRD), aircraft reaction force -Z (up).
// Physical motor command: sqrt(normalized thrust), for GZ velocity->force^2.
class ActuatorEffectivenessDeltaRay final : public ModuleParams, public ActuatorEffectiveness
{
public:
    explicit ActuatorEffectivenessDeltaRay(ModuleParams *parent);

    bool getEffectivenessMatrix(Configuration &configuration,
                                EffectivenessUpdateReason reason) override;

    void getDesiredAllocationMethod(AllocationMethod out[MAX_NUM_MATRICES]) const override
    {
        out[0] = AllocationMethod::PSEUDO_INVERSE;
    }

    void getNormalizeRPY(bool out[MAX_NUM_MATRICES]) const override { out[0] = false; }
    const char *name() const override { return "DeltaRay Virtual Forces"; }

    void updateSetpoint(const matrix::Vector<float, NUM_AXES> &control_sp, int matrix_index,
                        ActuatorVector &actuator_sp, const ActuatorVector &actuator_min,
                        const ActuatorVector &actuator_max) override;

    void getUnallocatedControl(int matrix_index, control_allocator_status_s &status) override;

private:
    static constexpr int N = 2;
    static constexpr float MAX_TILT_RATE = 2.0f;   // rad/s, tune to hardware
    static constexpr float MAX_THRUST_RATE = 3.0f; // fraction of full thrust per second
    static constexpr float PI_HALF = 1.5707963267948966f;

    void projectSector(float &h, float &v, float amin, float amax, float radius) const;
    bool validGeometry() const;
    matrix::Vector<float, NUM_AXES> actualWrench() const;

    ActuatorEffectivenessRotors _rotors;
    ActuatorEffectivenessTilts _tilts;
    float _thrust[N] {};
    float _alpha[N] {};
    matrix::Vector<float, NUM_AXES> _desired {};
    hrt_abstime _last_time{0};
    bool _ready{false};
    bool _have_request{false};
};
