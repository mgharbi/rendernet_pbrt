#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PATH_KPCN_H
#define PBRT_INTEGRATORS_PATH_KPCN_H

// integrators/pathkpcn.h*
#include "pbrt.h"
#include "integrator.h"

// PathKPCNIntegrator Declarations
class PathKPCNIntegrator : public SurfaceIntegrator {
public:
    // PathKPCNIntegrator Public Methods
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    Spectrum RecordedLi(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena, SampleRecord *sw, Camera *camera) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    PathKPCNIntegrator(int md) { maxDepth = md; }
private:
    // PathKPCNIntegrator Private Data
    int maxDepth;
#define SAMPLE_DEPTH 3
    LightSampleOffsets lightSampleOffsets[SAMPLE_DEPTH];
    int lightNumOffset[SAMPLE_DEPTH];
    BSDFSampleOffsets bsdfSampleOffsets[SAMPLE_DEPTH];
    BSDFSampleOffsets pathSampleOffsets[SAMPLE_DEPTH];
};


PathKPCNIntegrator *CreatePathKPCNSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PATH_KPCN_H
