#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_RENDERERS_RENDERNETRENDERER_H
#define PBRT_RENDERERS_RENDERNETRENDERER_H

// renderers/rendernetrenderer.h*
#include "pbrt.h"
#include "renderer.h"
#include "parallel.h"
#include "samplerecord.h"


// RendernetRenderer Declarations
class RendernetRenderer : public Renderer {
public:
    // RendernetRenderer Public Methods
    RendernetRenderer(Sampler *s, Camera *c, SurfaceIntegrator *si,
                    VolumeIntegrator *vi, int tileSize, int recordedSamples, bool useCameraSpaceNormals);
    ~RendernetRenderer();
    void Render(const Scene *scene);
    Spectrum Li(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL) const;
    RadianceQueryRecord RecordedLi(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena,
        Intersection *isect = NULL, Spectrum *T = NULL, SampleRecord *sr = NULL) const;
    Spectrum Transmittance(const Scene *scene, const RayDifferential &ray,
        const Sample *sample, RNG &rng, MemoryArena &arena) const;
    int tileSize;
    int recordedSamples;
    int maxDepth;
    bool useCameraSpaceNormals;
    bool is_kpcn() const;
private:
    // RendernetRenderer Private Data
    Sampler *sampler;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
};



// RendernetRendererTask Declarations
class RendernetRendererTask : public Task {
public:
    // RendernetRendererTask Public Methods
    RendernetRendererTask(const Scene *sc, RendernetRenderer *ren, Camera *c,
                        ProgressReporter &pr, Sampler *ms, Sample *sam, 
                        int tn, int tc)
      : reporter(pr), max_radiance(0.f)
    {
        scene = sc; renderer = ren; camera = c; mainSampler = ms;
        origSample = sam; taskNum = tn; taskCount = tc;
    }
    void Run();

    char fname[64];
    float max_radiance;
private:
    // RendernetRendererTask Private Data
    const Scene *scene;
    const RendernetRenderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    int taskNum, taskCount;

};



#endif // PBRT_RENDERERS_RENDERNETRENDERER_H
