
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// renderers/samplerrenderer.cpp*
#include "stdafx.h"
#include "renderers/rendernetrenderer.h"
#include "integrators/pathrendernet.h"
#include "scene.h"
#include "film.h"
#include "volume.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "intersection.h"
#include "montecarlo.h"
#include "cameras/perspective.h"

static uint32_t hash(char *key, uint32_t len)
{
    uint32_t hash = 0, i;
    for (hash=0, i=0; i<len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
} 

// RendernetRendererTask Definitions
void RendernetRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _RendernetRendererTask_
    Sampler *sampler = mainSampler->GetSubSampler(taskNum, taskCount);
    if (!sampler)
    {
        reporter.Update();
        PBRT_FINISHED_RENDERTASK(taskNum);
        return;
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(taskNum);

    // Allocate space for samples and intersections
    int maxSamples = sampler->MaximumSampleCount() + renderer->recordedSamples;
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];

    Point sceneCenter;
    float sceneRadius;
    scene->WorldBound().BoundingSphere(&sceneCenter, &sceneRadius);

    // TODO: use dynamic cast for camera
    SampleRecord *sr = new SampleRecord(
        sampler->xPixelStart,
        sampler->yPixelStart,
        renderer->tileSize, 
        renderer->recordedSamples,
        sampler->samplesPerPixel,
        renderer->maxDepth,
        camera->film->xResolution, camera->film->yResolution,
        sceneRadius, ((ProjectiveCamera*) camera)->focalDistance,
        ((ProjectiveCamera*) camera)->lensRadius,
        ((PerspectiveCamera*) camera)->fov,
        renderer->useCameraSpaceNormals
        );

    // Get samples from _Sampler_ and update image
    int sampleCount;
    int pixel_id = 0;
    while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        int n_added = 0;
        int n_added_lowspp = 0;
        Spectrum radiance = 0.0f;
        Spectrum variance = 0.0f;
        Spectrum lowspp_radiance = 0.0f;
        Spectrum lowspp_variance = 0.0f;

        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
            float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
            PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);
            if(rayWeight != 1.0f) {
              printf("weight %.2f\n", rayWeight);
            }

            // Evaluate radiance along camera ray
            PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
            if (rayWeight > 0.f) {
              if(i < renderer->recordedSamples) {
                Ls[i] = rayWeight * renderer->RecordedLi(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i], sr);
              }
              else {
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
              }
            }
            else {
                Error("Ray has zero weight: not handled by sample saver!");
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }

            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned "
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);

            // First "recordedSamples" go to .bin, the rest to ground truth
            if(i < renderer->recordedSamples) {
              int pix_x = pixel_id % renderer->tileSize + sr->tile_x;
              int pix_y = pixel_id / renderer->tileSize + sr->tile_y;
              sr->pixel_x.push_back((float) pix_x);
              sr->pixel_y.push_back((float) pix_y);
              sr->subpixel_x.push_back(samples[i].imageX-(float)pix_x);
              sr->subpixel_y.push_back(samples[i].imageY-(float)pix_y);
              
              float lensU, lensV;
              ConcentricSampleDisk(samples[i].lensU, samples[i].lensV, &lensU, &lensV);
              lensU *= sr->aperture_radius;
              lensV *= sr->aperture_radius;
              sr->lens_u.push_back(lensU);
              sr->lens_v.push_back(lensV);

              sr->time.push_back(samples[i].time);

              Spectrum current = lowspp_radiance;
              Spectrum delta = (Ls[i] - current);
              current += delta/(n_added_lowspp+1);
              lowspp_variance += delta*(Ls[i] - current);
              lowspp_radiance += (Ls[i]-lowspp_radiance)/(n_added_lowspp+1);
              n_added_lowspp += 1;
            } else {  // skipping "recordedSamples" from the groundtruth ensure i/o are decorrelated
              Spectrum current = radiance;
              Spectrum delta = (Ls[i] - current);
              current += delta/(n_added+1);
              variance += delta*(Ls[i] - current);
              radiance += (Ls[i]-radiance)/(n_added+1);
              // Update radiance mean and std
              n_added += 1;
            }
            
        }

        // Add rendered pixel to record
        sr->ground_truth.push_back(radiance.ToRGBSpectrum());
        sr->ground_truth_variance.push_back(variance.ToRGBSpectrum());
        sr->lowspp.push_back(lowspp_radiance.ToRGBSpectrum());
        sr->lowspp_variance.push_back(lowspp_variance.ToRGBSpectrum());

        // Report sample results to _Sampler_, add contributions to image
        if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        {
            for (int i = 0; i < sampleCount; ++i)
            {
                PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
                camera->film->AddSample(samples[i], Ls[i]);
                PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
            }
        }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();

        // Increment pixel counter
        pixel_id += 1;
    }
    
    // Write sample data
    sprintf(fname, "%04d_%04d.bin", sampler->xPixelStart, sampler->yPixelStart);
    sr->save(fname);

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(sampler->xPixelStart,
        sampler->yPixelStart, sampler->xPixelEnd+1, sampler->yPixelEnd+1);
    delete sr;
    delete sampler;
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// RendernetRenderer Method Definitions
RendernetRenderer::RendernetRenderer(Sampler *s, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 int tSz, int recSamples, bool useCamSpaceNrm) {
    sampler = s;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    tileSize = tSz;
    recordedSamples = recSamples;
    useCameraSpaceNormals = useCamSpaceNrm;

    maxDepth = surfaceIntegrator->maxDepth();
    printf("sI's maxdepth %d\n", maxDepth);
}


RendernetRenderer::~RendernetRenderer() {
    delete sampler;
    delete camera;
    delete surfaceIntegrator;
    delete volumeIntegrator;
}


void RendernetRenderer::Render(const Scene *scene) {
    PBRT_FINISHED_PARSING();
    // Allow integrators to do preprocessing for the scene
    PBRT_STARTED_PREPROCESSING();
    surfaceIntegrator->Preprocess(scene, camera, this);
    volumeIntegrator->Preprocess(scene, camera, this);
    PBRT_FINISHED_PREPROCESSING();
    PBRT_STARTED_RENDERING();
    // Allocate and initialize _sample_
    Sample *sample = new Sample(sampler, surfaceIntegrator,
                                volumeIntegrator, scene);

    // Create and launch _RendernetRendererTask_s for rendering image

    // Compute number of _RendernetRendererTask_s to create for rendering
    
    int xstart, xend, ystart, yend;
    camera->film->GetPixelExtent(&xstart, &xend, &ystart, &yend);
    int xRes = xend-xstart;
    int yRes = yend-ystart;
    if(xRes % tileSize != 0) {
      Error("tile size does not divide xRes");
      return;
    }
    if(yRes % tileSize != 0) {
      Error("tile size does not divide yRes");
      return;
    }
    int nTasks = xRes*yRes / (tileSize*tileSize);

    int maxSamples = sampler->MaximumSampleCount();
    if(recordedSamples > maxSamples) {
      Error("You asked to record more samples (%d) than"
          " what is used for ground-truth (%d)", recordedSamples, maxSamples);
      return;
    }

    printf("Resolution %dx%d, %d tiles with size %d. Saving %d samples (ground truth at %d)\n",
        xRes, yRes, nTasks, tileSize, recordedSamples, maxSamples);

    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) {
      renderTasks.push_back(new RendernetRendererTask(scene, this, camera,
                                                    reporter, sampler, sample, 
                                                    nTasks-1-i, nTasks));
    }
    EnqueueTasks(renderTasks);
    WaitForAllTasks();

    for (uint32_t i = 0; i < renderTasks.size(); ++i) {
      delete renderTasks[i];
    }

    reporter.Done();
    PBRT_FINISHED_RENDERING();
    // Clean up after rendering and store final image
    delete sample;
    camera->film->WriteImage();
}

Spectrum RendernetRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
  return RecordedLi(scene, ray, sample, rng, arena, isect, T, NULL);
}


Spectrum RendernetRenderer::RecordedLi(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T, SampleRecord * sr) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    Spectrum Li = 0.f;
    if (scene->Intersect(ray, isect)) {
        Li = surfaceIntegrator->RecordedLi(scene, this, ray, *isect, sample,
                                   rng, arena, sr, camera);
    }
    else {
        // Handle ray that doesn't intersect any geometry
        Spectrum contrib;
        for (uint32_t i = 0; i < scene->lights.size(); ++i) {
           Li += scene->lights[i]->Le(ray);
        }

        if(sr) {
          Transform tx;
          camera->CameraToWorld.Interpolate(sample->time, &tx);
          Spectrum zero = 0.;
          Normal default_n;
          sr->normal_at_first.push_back(default_n);
          sr->depth_at_first.push_back(-1.0f);  // no intersection
          sr->normal.push_back(default_n);
          sr->depth.push_back(-1.0f);  // no intersection
          sr->visibility.push_back(0.0f);
          sr->albedo.push_back(Spectrum(0.0f));
          sr->albedo_at_first.push_back(Spectrum(0.0f));

          sr->radiance_diffuse.push_back(zero);
          sr->radiance_diffuse_indirect.push_back(zero);
          sr->radiance_specular.push_back(Li); // We only have the scene lights/envmap contributions
          std::vector<float> p(4*sr->maxDepth);
          sr->probabilities.push_back(p);
          std::vector<float> ld(2*sr->maxDepth);
          sr->light_directions.push_back(ld);
          std::vector<uint16_t> bt(sr->maxDepth);
          sr->bounce_type.push_back(bt);
        }
    }
    // NOTE(mgharbi): volume not accounted for, for now
    Spectrum Lvi = volumeIntegrator->Li(scene, this, ray, sample, rng,
                                        T, arena);
    // return *T * Li + Lvi;
    
    return *T * Li;
}


Spectrum RendernetRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


