// integrators/pathkpcn.cpp*
#include "stdafx.h"
#include "integrators/pathkpcn.h"
#include "core/camera.h"
#include "scene.h"
#include "intersection.h"
#include "paramset.h"

// PathKPCNIntegrator Method Definitions
void PathKPCNIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                    const Scene *scene) {
    for (int i = 0; i < SAMPLE_DEPTH; ++i) {
        lightSampleOffsets[i] = LightSampleOffsets(1, sample);
        lightNumOffset[i] = sample->Add1D(1);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(1, sample);
        pathSampleOffsets[i] = BSDFSampleOffsets(1, sample);
    }
}

Spectrum PathKPCNIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
  return RecordedLi(scene, renderer, r, isect, sample, rng, arena, NULL, NULL);
}


Spectrum PathKPCNIntegrator::RecordedLi(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena, SampleRecord *sr, Camera* camera) const {
    // Declare common path integration variables
    Spectrum pathThroughput = 1., L = 0.;
    Spectrum pathThroughputDiffuse = 1., Ldiffuse = 0.;
    RayDifferential ray(r);
    bool specularBounce = true;
    bool foundRough = false;
    Intersection localIsect;
    const Intersection *isectp = &isect;

    Spectrum runningAlbedo = 1.f;
    Spectrum visibility = 0.f;
    Spectrum visibilityRecord = 0.0f;

    bool recordedOutputValues = false;
    float hitDistance = 0.0f;

    // Default values
    Normal nrm;
    float depth = 0.0f;
    Spectrum albedo = 0.0f;

    for (int bounces = 0; ; ++bounces) {
        // Possibly add emitted light at path vertex
        // NOTE: tungsten samples this at all bounces
        if (bounces == 0 || specularBounce) {
          Spectrum contrib = isectp->Le(-ray.d);
          L += contrib*pathThroughput;
          if (foundRough) {
            Ldiffuse += contrib*pathThroughputDiffuse;
          }
        }

        // Sample illumination from lights to find path contribution
        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;

        bool bsdf_has_diffuse = false;
        if (bsdf) {
          bsdf_has_diffuse =
              (bsdf->NumComponents(BxDFType(BSDF_DIFFUSE|BSDF_REFLECTION)) > 0) ||
              (bsdf->NumComponents(BxDFType(BSDF_DIFFUSE|BSDF_TRANSMISSION)) > 0);
        }

        Vector depth_vector = p-ray.o;
        hitDistance += depth_vector.Length();

        Vector wo = -ray.d;
        Spectrum contrib(0.0f);
        Spectrum diffuse_lighting(0.0f);
        LightQueryRecord qr;
        if (bounces < SAMPLE_DEPTH) {
            contrib = UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng,
                     lightNumOffset[bounces], &lightSampleOffsets[bounces],
                     &bsdfSampleOffsets[bounces], &qr);
        } else {
            contrib = UniformSampleOneLight(scene, renderer, arena, p, n, wo,
                     isectp->rayEpsilon, ray.time, bsdf, sample, rng, 
                     -1, NULL, NULL, &qr);
        }

        // Binarize visibility to match KPCN's Tungsten output:
        if (qr.visibility != 0.0f) {
          qr.visibility = 1.0f;
        }

        if (!foundRough && bsdf_has_diffuse) {
          L += contrib*pathThroughput;
          Ldiffuse += qr.diffuse_lighting*pathThroughputDiffuse;
        } else {
          L += contrib*pathThroughput;
          if (foundRough) {
            Ldiffuse += contrib*pathThroughputDiffuse;
          }
        }
        
        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        BSDFSample outgoingBSDFSample;
        if (bounces < SAMPLE_DEPTH)
            outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[bounces],
                                            0);
        else
            outgoingBSDFSample = BSDFSample(rng);
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                    BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.) {
          depth = hitDistance;
          nrm = n;
          albedo = 0.0f;
          visibility = qr.visibility;
          recordedOutputValues = true;

          break;
        } 

        Spectrum bsdfWeight =  f * AbsDot(wi, n) / pdf;
        pathThroughput *= bsdfWeight;

        specularBounce = (flags & BSDF_SPECULAR) != 0;

        // If the brdf has a diffuse component, we found our first
        // diffuse interaction. The path is no longer purely specular.
        if (!foundRough && bsdf_has_diffuse > 0) {
          // TODO BSDF_TRANSMISSION too?
          Spectrum bsdfWeightDiffuse = specularBounce ? Spectrum(0.0f) : 
            bsdf->f(wo, wi, BxDFType(BSDF_DIFFUSE|BSDF_REFLECTION)) * AbsDot(wi, n) / pdf;
          pathThroughputDiffuse *= bsdfWeightDiffuse;
          foundRough = true;
        } else {
          pathThroughputDiffuse *= bsdfWeight;
        }

        // record value at first rough 
        if (sr && !recordedOutputValues && foundRough) {
          depth = hitDistance;
          nrm = n;
          albedo = runningAlbedo*bsdfWeight;
          visibility = qr.visibility;
          recordedOutputValues = true;

        } else if (!foundRough) {
          runningAlbedo *= bsdfWeight;
        }  // record values

        // Scatter
        ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);

        // Possibly terminate the path (russian roulette)
        if (bounces > 3) {
            float continueProbability = min(.5f, pathThroughput.y());
            if (rng.RandomFloat() > continueProbability)
                break;
            pathThroughput /= continueProbability;
        }
        if (bounces == maxDepth_)
            break;

        // Find next vertex of path
        if (!scene->Intersect(ray, &localIsect)) {
            if (specularBounce) {
                for (uint32_t i = 0; i < scene->lights.size(); ++i) {
                  Spectrum contrib = scene->lights[i]->Le(ray);
                  L += contrib*pathThroughput;
                  if (foundRough) {
                    Ldiffuse += contrib*pathThroughputDiffuse;
                  }
                }
            }
            break;
        }
        Spectrum transmittance = renderer->Transmittance(scene, ray, NULL, rng, arena);
        pathThroughput *= transmittance;
        pathThroughputDiffuse *= transmittance;
        isectp = &localIsect;
    }  // bounces loop

    if (sr) {
      sr->normal.push_back(nrm);
      sr->normal_at_first.push_back(nrm);
      sr->depth.push_back(depth);
      sr->depth_at_first.push_back(depth);
      sr->visibility.push_back(visibility.y());
      sr->albedo.push_back(albedo);
      sr->albedo_at_first.push_back(albedo);

      sr->radiance_diffuse.push_back(Ldiffuse);
      sr->radiance_diffuse_indirect.push_back(Spectrum(0.0f));
      sr->radiance_specular.push_back(L-Ldiffuse);

      // Fill our features with dummies
      std::vector<float> probabilities(4*(maxDepth_+1), 0.0f);
      std::vector<float> light_directions(2*(maxDepth_+1), 0.0f);
      std::vector<uint16_t> bounce_type((maxDepth_+1), 0);
      sr->probabilities.push_back(probabilities);
      sr->light_directions.push_back(light_directions);
      sr->bounce_type.push_back(bounce_type);
    }

    return L;
}


PathKPCNIntegrator *CreatePathKPCNSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    printf("Max depth %d\n", maxDepth);
    return new PathKPCNIntegrator(maxDepth);
}
