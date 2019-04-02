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
    RendernetRenderer(
        Sampler *s, Sampler *s2, Sampler *rs, Camera *c, SurfaceIntegrator *si,
        VolumeIntegrator *vi, int tileSize, 
        int recordedSamples, bool useCameraSpaceNormals);
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
    // Sampler *sampler2;
    Sampler *sampler_recorded;
    Camera *camera;
    SurfaceIntegrator *surfaceIntegrator;
    VolumeIntegrator *volumeIntegrator;
};



// RendernetRendererTask Declarations
class RendernetRendererTask : public Task {
public:
    // RendernetRendererTask Public Methods
    RendernetRendererTask(const Scene *sc, RendernetRenderer *ren, Camera *c,
                        ProgressReporter &pr, 
                        Sampler *ms, Sampler* rs,
                        Sample *sam, Sample *rsam,
                        int tn, int tc)
      : reporter(pr)
    {
        scene = sc; renderer = ren; camera = c; 
        mainSampler = ms; 
        // mainSampler2 = ms2; 
        recordedSampler = rs;
        origSample = sam; 
        // origSample2 = sam2; 
        recordedOrigSample = rsam; 
        taskNum = tn; taskCount = tc;
    }
    void Run();

    char fname[64];
private:
    // RendernetRendererTask Private Data
    const Scene *scene;
    const RendernetRenderer *renderer;
    Camera *camera;
    Sampler *mainSampler;
    // Sampler *mainSampler2;
    Sampler *recordedSampler;
    ProgressReporter &reporter;
    Sample *origSample;
    // Sample *origSample2;
    Sample *recordedOrigSample;
    int taskNum, taskCount;

};



#endif // PBRT_RENDERERS_RENDERNETRENDERER_H
