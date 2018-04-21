#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_PATH_RENDERNET_H
#define PBRT_INTEGRATORS_PATH_RENDERNET_H

// integrators/pathrendernet.h*
#include "pbrt.h"
#include "integrator.h"
#include "samplerecord.h"

// PathRendernetIntegrator Declarations
class PathRendernetIntegrator : public SurfaceIntegrator {
public:
    // PathRendernetIntegrator Public Methods
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    RadianceQueryRecord RecordedLi(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena, SampleRecord *sw, Camera *camera) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    PathRendernetIntegrator(int md) { maxDepth_ = md; }
    virtual int maxDepth() {return maxDepth_;};
private:
    int maxDepth_;
    // PathRendernetIntegrator Private Data
#define SAMPLE_DEPTH 3
    LightSampleOffsets lightSampleOffsets[SAMPLE_DEPTH];
    int lightNumOffset[SAMPLE_DEPTH];
    BSDFSampleOffsets bsdfSampleOffsets[SAMPLE_DEPTH];
    BSDFSampleOffsets pathSampleOffsets[SAMPLE_DEPTH];
};


PathRendernetIntegrator *CreatePathRendernetSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_PATH_RENDERNET_H
