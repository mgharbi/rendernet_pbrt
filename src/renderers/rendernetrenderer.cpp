
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
#include "integrators/pathkpcn.h"
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
#include "time.h"

#include <typeinfo>

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

bool RendernetRenderer::is_kpcn() const {
  PathKPCNIntegrator* p = dynamic_cast<PathKPCNIntegrator*>(surfaceIntegrator);
  if(p) {
    return true;
  }
  return false;
}

// RendernetRendererTask Definitions
void RendernetRendererTask::Run() {
    PBRT_STARTED_RENDERTASK(taskNum);
    // Get sub-_Sampler_ for _RendernetRendererTask_
    //
    Sampler* samplers[3] = {
      mainSampler->GetSubSampler(taskNum, taskCount), 
      mainSampler2->GetSubSampler(taskNum, taskCount), 
      recordedSampler->GetSubSampler(taskNum, taskCount)
    };

    Sample* origSamples[3] = {
      origSample,
      origSample2,
      recordedOrigSample,
    };

    for (int i = 0; i < 3; ++i) {
      if (!samplers[i]) {
          reporter.Update();
          PBRT_FINISHED_RENDERTASK(taskNum);
          return;
      }
    }

    // Declare local variables used for rendering loop
    MemoryArena arena;
    RNG rng(time(NULL)); // RNG rng(taskNum);

    Point sceneCenter;
    float sceneRadius;
    scene->WorldBound().BoundingSphere(&sceneCenter, &sceneRadius);

    PerspectiveCamera *pcam = dynamic_cast<PerspectiveCamera*>(camera);
    if (!pcam) {
      Error("Rendernet only supports ProjectiveCamera\n");
      return;
    }

    int xstart, xend, ystart, yend;
    camera->film->GetPixelExtent(&xstart, &xend, &ystart, &yend);
    printf("film: %d %d\n", xend-xstart, yend-ystart);
      

    SampleRecord *sr = new SampleRecord(
        samplers[0]->xPixelStart,
        samplers[0]->yPixelStart,
        renderer->tileSize, 
        samplers[2]->samplesPerPixel,  // Saved samples
        samplers[0]->samplesPerPixel,  // Image reference
        renderer->maxDepth,
        xend-xstart,
        yend-xstart,
        // camera->film->xResolution, 
        // camera->film->yResolution,
        sceneRadius, pcam->focalDistance,
        pcam->lensRadius, pcam->fov,
        renderer->useCameraSpaceNormals
        );

    if(renderer->is_kpcn()) {
      sr->set_kpcn();
    }

    // Allocate space for samples and intersections
    for (int sampler_idx = 0; sampler_idx < 3; ++sampler_idx) {
      Sampler *sampler = samplers[sampler_idx];

      int maxSamples = sampler->MaximumSampleCount();
      // TODO: isolate origsamples per sampler
      Sample *samples = origSamples[sampler_idx]->Duplicate(maxSamples);
      RayDifferential *rays = new RayDifferential[maxSamples];
      Spectrum *Ts = new Spectrum[maxSamples];
      Intersection *isects = new Intersection[maxSamples];

      // Get samples from _Sampler_ and update image
      int sampleCount;
      int pixel_id = 0;
      while ((sampleCount = sampler->GetMoreSamples(samples, rng)) > 0) {
        RadianceQueryRecord rq_rec;

        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
          // Find camera ray for _sample[i]_
          // PBRT_STARTED_GENERATING_CAMERA_RAY(&samples[i]);
          float rayWeight = camera->GenerateRayDifferential(samples[i], &rays[i]);
          rays[i].ScaleDifferentials(1.f / sqrtf(sampler->samplesPerPixel));
          // PBRT_FINISHED_GENERATING_CAMERA_RAY(&samples[i], &rays[i], rayWeight);

          // Evaluate radiance along camera ray
          // PBRT_STARTED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i]);
          if (sampler_idx == 2) {
            // we save the sample, and ignore the buffer data
            RadianceQueryRecord ret = renderer->RecordedLi(scene, rays[i], &samples[i], rng,
                arena, &isects[i], &Ts[i], sr); 

            // Record sample data
            int pix_x = pixel_id % renderer->tileSize + sr->tile_x;
            int pix_y = pixel_id / renderer->tileSize + sr->tile_y;
            float lensU, lensV;
            ConcentricSampleDisk(samples[i].lensU, samples[i].lensV, &lensU, &lensV);
            lensU *= sr->aperture_radius;
            lensV *= sr->aperture_radius;
            sr->pixel_x.push_back((float) pix_x);
            sr->pixel_y.push_back((float) pix_y);
            sr->subpixel_x.push_back(samples[i].imageX-(float)pix_x);
            sr->subpixel_y.push_back(samples[i].imageY-(float)pix_y);
            sr->lens_u.push_back(lensU);
            sr->lens_v.push_back(lensV);
            sr->time.push_back(samples[i].time);
          } else {
            RadianceQueryRecord ret = renderer->RecordedLi(scene, rays[i], &samples[i], rng,
                arena, &isects[i], &Ts[i], NULL);
            rq_rec.add(ret, rayWeight);
          }
          // PBRT_FINISHED_CAMERA_RAY_INTEGRATION(&rays[i], &samples[i], &Ls[i]);
          
        } // spp loop

        // We're constructing an image
        if( sampler_idx < 2 ) {
          // Add pixel data to .bin record
          sr->add_image_sample(rq_rec, sampler_idx);
        }

        // // Report sample results to _Sampler_, add contributions to image
        // if (sampler->ReportResults(samples, rays, Ls, isects, sampleCount))
        // {
        //     for (int i = 0; i < sampleCount; ++i)
        //     {
        //         PBRT_STARTED_ADDING_IMAGE_SAMPLE(&samples[i], &rays[i], &Ls[i], &Ts[i]);
        //         camera->film->AddSample(samples[i], Ls[i]);
        //         PBRT_FINISHED_ADDING_IMAGE_SAMPLE();
        //     }
        // }

        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();

        // Increment pixel counter
        pixel_id += 1;
      }

      delete[] samples;
      delete[] rays;
      delete[] Ts;
      delete[] isects;
    } // Loop over samplers
    
    // Write sample data
    sprintf(fname, "%04d_%04d.bin", samplers[0]->xPixelStart, samplers[0]->yPixelStart);
    sr->save(fname);

    // Clean up after _SamplerRendererTask_ is done with its image region
    camera->film->UpdateDisplay(samplers[0]->xPixelStart,
        samplers[0]->yPixelStart, samplers[0]->xPixelEnd+1, samplers[0]->yPixelEnd+1);
    delete sr;
    for (int i = 0; i < 3; ++i) {
      delete samplers[i];
    }
    // delete[] Ls;
    reporter.Update();
    PBRT_FINISHED_RENDERTASK(taskNum);
}



// RendernetRenderer Method Definitions
RendernetRenderer::RendernetRenderer(Sampler *s, Sampler *s2, Sampler *rs, Camera *c,
                                 SurfaceIntegrator *si, VolumeIntegrator *vi,
                                 int tSz, int recSamples, bool useCamSpaceNrm) {
    sampler = s;
    sampler2 = s2;
    sampler_recorded = rs;
    camera = c;
    surfaceIntegrator = si;
    volumeIntegrator = vi;
    tileSize = tSz;
    recordedSamples = recSamples;
    useCameraSpaceNormals = useCamSpaceNrm;

    maxDepth = surfaceIntegrator->maxDepth();
    if(maxDepth != 5) {
      Error("Rendernet's sampler structure only supports path length 5.\n");
    }
}


RendernetRenderer::~RendernetRenderer() {
    delete sampler;
    delete sampler2;
    delete sampler_recorded;
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
    Sample *sample2 = new Sample(sampler2, surfaceIntegrator,
                                volumeIntegrator, scene);
    Sample *rsample = new Sample(sampler_recorded, surfaceIntegrator,
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

    printf("Resolution %dx%d, %d tiles with size %d. References with %d and %d samples. Input with %d samples)\n",
        xRes, yRes, nTasks, tileSize, sampler->MaximumSampleCount(), sampler2->MaximumSampleCount(), sampler_recorded->MaximumSampleCount());

    ProgressReporter reporter(nTasks, "Rendering");
    vector<Task *> renderTasks;
    for (int i = 0; i < nTasks; ++i) {
      renderTasks.push_back(
          new RendernetRendererTask(
            scene, this, camera, reporter, sampler, sampler2,
            sampler_recorded, sample, sample2, rsample, nTasks-1-i, nTasks));
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
    delete sample2;
    delete rsample;
    camera->film->WriteImage();
}

Spectrum RendernetRenderer::Li(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T) const {
  // return RecordedLi(scene, ray, sample, rng, arena, isect, T, NULL).L;
  throw;
  return Spectrum(0.0f);
}


RadianceQueryRecord RendernetRenderer::RecordedLi(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena, Intersection *isect, Spectrum *T, SampleRecord * sr) const {
    Assert(ray.time == sample->time);
    Assert(!ray.HasNaNs());
    // Allocate local variables for _isect_ and _T_ if needed
    Spectrum localT;
    if (!T) T = &localT;
    Intersection localIsect;
    if (!isect) isect = &localIsect;
    RadianceQueryRecord rq_rec;
    if (scene->Intersect(ray, isect)) {
        rq_rec = surfaceIntegrator->RecordedLi(
            scene, this, ray, *isect, sample, rng, arena, sr, camera);
    } else {
        // Handle ray that doesn't intersect any geometry
        Spectrum Li;
        for (uint32_t i = 0; i < scene->lights.size(); ++i) {
           Li += scene->lights[i]->Le(ray);
        }

        // No intersection
        Normal default_n;
        rq_rec = RadianceQueryRecord(
            Li, Spectrum(0.0f), Spectrum(0.0f), default_n, -1.0f, false, false);

        if(sr) {
          // Transform tx;
          // camera->CameraToWorld.Interpolate(sample->time, &tx);
          Spectrum zero = 0.;
          sr->radiance_diffuse.push_back(zero);
          sr->radiance_diffuse_indirect.push_back(zero);
          sr->radiance_specular.push_back(Li); // We only have the scene lights/envmap contributions

          sr->albedo.push_back(Spectrum(0.0f));
          sr->albedo_at_first.push_back(Spectrum(0.0f));
          sr->normal_at_first.push_back(default_n);
          sr->normal.push_back(default_n);

          sr->depth_at_first.push_back(-1.0f);  // no intersection
          sr->depth.push_back(-1.0f);  // no intersection

          sr->visibility.push_back(0.0f);
          sr->hasHit.push_back(0.0f);

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
    
    // TODO: multiply radiance by *T if using transmissive media
    return rq_rec;
}


Spectrum RendernetRenderer::Transmittance(const Scene *scene,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    return volumeIntegrator->Transmittance(scene, this, ray, sample,
                                           rng, arena);
}


