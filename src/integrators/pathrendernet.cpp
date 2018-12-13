// integrators/pathrendernet.cpp*
#include "stdafx.h"
#include "integrators/pathrendernet.h"
#include "core/camera.h"
#include "scene.h"
#include "intersection.h"
#include "paramset.h"
#include "core/samplerecord.h"

// PathRendernetIntegrator Method Definitions
void PathRendernetIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
                                    const Scene *scene) {
    for (int i = 0; i < SAMPLE_DEPTH; ++i) {
        lightSampleOffsets[i] = LightSampleOffsets(1, sample);
        lightNumOffset[i] = sample->Add1D(1);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(1, sample);
        pathSampleOffsets[i] = BSDFSampleOffsets(1, sample);
    }
}

Spectrum PathRendernetIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {
  // return RecordedLi(scene, renderer, r, isect, sample, rng, arena, NULL, NULL).L;
  throw;
  return Spectrum(0.0f); // TODO: fixme
}

RadianceQueryRecord PathRendernetIntegrator::RecordedLi(const Scene *scene, const Renderer *renderer,
        const RayDifferential &r, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena, SampleRecord *sr, Camera* camera) const {
    // Declare common path integration variables
    Spectrum pathThroughput = 1., L = 0.;
    Spectrum pathThroughputDiffuse = 1., Ldiffuse = 0., Ldiffuse_indirect = 0.;
    RayDifferential ray(r);
    bool specularBounce = false;
    bool foundRough = false;
    bool foundNonSpecular = false;
    Intersection localIsect;
    const Intersection *isectp = &isect;

    // Default values
    bool isLightVisible = false;
    bool recordedOutputValues = false;
    float hitDistance = 0.0f;
    Normal nrm;
    Normal nrm_at_first;
    float depth = 0.0f;
    float depth_at_first = 0.0f;
    Spectrum albedo = 0.f;
    Spectrum albedo_at_first = 0.f;

    std::vector<float> probabilities(4*(maxDepth_+1), 0.0f);
    std::vector<float> light_directions(2*(maxDepth_+1), 0.0f);
    std::vector<uint16_t> bounce_type((maxDepth_+1), 0);

    for (int bounces = 0; ; ++bounces) {
        // Possibly add emitted light at path vertex
        if (bounces == 0 || specularBounce) {
          Spectrum contrib = isectp->Le(-ray.d);
          L += contrib*pathThroughput;
          if (foundRough) {
            Ldiffuse += contrib*pathThroughputDiffuse;
            if (bounces > 0)  {
              Ldiffuse_indirect += contrib*pathThroughputDiffuse;
            }
          }
        }

        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;

        bool bsdf_has_diffuse = false;
        bsdf_has_diffuse =
            (bsdf->NumComponents(BxDFType(BSDF_DIFFUSE|BSDF_REFLECTION)) > 0);
        bool bsdf_has_nonspecular = false;
        bsdf_has_nonspecular = bsdf_has_diffuse ||
            (bsdf->NumComponents(BxDFType(BSDF_GLOSSY|BSDF_REFLECTION)) > 0) ||
            (bsdf->NumComponents(BxDFType(BSDF_GLOSSY|BSDF_TRANSMISSION)) > 0);

        Vector depth_vector = p-ray.o;
        hitDistance += depth_vector.Length();

        // Sample illumination from lights to find path contribution
        Vector wo = -ray.d;
        Spectrum contrib(0.0f);
        Spectrum diffuse_lighting(0.0f);
        Transform tx;
        camera->CameraToWorld.Interpolate(sample->time, &tx);
        tx = Inverse(tx);
        LightQueryRecord qr(tx);
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
        
        L += contrib*pathThroughput;
        if (!foundRough && bsdf_has_diffuse) {
          Ldiffuse += qr.diffuse_lighting*pathThroughputDiffuse;
          if (bounces > 0)  {
            Ldiffuse_indirect += qr.diffuse_lighting*pathThroughputDiffuse;
          }
        } else if (foundRough) {
          Ldiffuse += contrib*pathThroughputDiffuse;
          if (bounces > 0)  {
            Ldiffuse_indirect += contrib*pathThroughputDiffuse;
          }
        }
        std::copy(qr.pdfs, qr.pdfs+4, probabilities.begin()+4*bounces);
        light_directions[2*bounces + 0] = qr.theta;
        light_directions[2*bounces + 1] = qr.phi;

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
        bounce_type[bounces] = flags;
        Spectrum currAlbedo = bsdf->K();

        // If the brdf has a diffuse component, we found our first
        // diffuse interaction. The path is no longer purely specular.
        bool isFirstRough = false;
        if (!foundRough && bsdf_has_diffuse) {
          foundRough = true;
          isFirstRough = true;
        } 

        bool isFirstNonSpecular = false;
        if (!foundNonSpecular && bsdf_has_nonspecular) {
          foundNonSpecular = true;
          isFirstNonSpecular = true;
        } 

        // Save depth, normal, albedo at first  (should it be first non-specular?)
        // if (isFirstNonSpecular) {
        if (bounces == 0) {
          Normal ssn(n);
          if (Dot(ssn, ray.d) < 0) { //face forward
            ssn.x *= -1.0f;
            ssn.y *= -1.0f;
            ssn.z *= -1.0f;
          }

          // Camera-space normals
          Transform tx;
          camera->CameraToWorld.Interpolate(sample->time, &tx);
          nrm_at_first = Inverse(tx)(ssn);

          depth_at_first = hitDistance;
          albedo_at_first = currAlbedo;

          isLightVisible = isLightVisible || qr.isLightVisible;
        }

        // record value at first rough 
        if (sr && !recordedOutputValues && foundRough) {
          recordedOutputValues = true;
          depth = hitDistance;
          albedo = currAlbedo;
          Normal ssn(n);
          if (Dot(ssn, ray.d) < 0) { //face forward
            ssn.x *= -1.0f;
            ssn.y *= -1.0f;
            ssn.z *= -1.0f;
          }

          // Camera-space normals
          Transform tx;
          camera->CameraToWorld.Interpolate(sample->time, &tx);
          nrm = Inverse(tx)(ssn);
        } 

        if (f.IsBlack() || pdf == 0.) { // Stop propagation
          break;
        } 

        Spectrum bsdfWeight =  f * AbsDot(wi, n) / pdf;
        pathThroughput *= bsdfWeight;

        specularBounce = (flags & BSDF_SPECULAR) != 0;

        if (bsdfWeight.HasNaNs() || isinf(bsdfWeight.y())) {
          Error("Not-a-number in bsdfweight");
        }

        if(isFirstRough) {
          Spectrum bsdfWeightDiffuse = specularBounce ? Spectrum(0.0f) : 
            bsdf->f(wo, wi, BxDFType(BSDF_DIFFUSE|BSDF_REFLECTION|BSDF_GLOSSY)) * AbsDot(wi, n) / pdf;
          pathThroughputDiffuse *= bsdfWeightDiffuse;
        } else {
          pathThroughputDiffuse *= bsdfWeight;
        }
        
        // Scatter
        ray = RayDifferential(p, wi, ray, isectp->rayEpsilon);

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
                    if (bounces > 0)  {
                      Ldiffuse_indirect += contrib*pathThroughputDiffuse;
                    }
                  }
                }
            }
            break;
        }
        Spectrum transmittance = renderer->Transmittance(scene, ray, NULL, rng, arena);
        pathThroughput *= transmittance;
        pathThroughputDiffuse *= transmittance;
        isectp = &localIsect;
    } // bounces loop

    if (sr) {
      if (nrm_at_first.HasNaNs()) {
        Error("normal first has nans");
        // nrm_at_first.x = 0.0f;
        // nrm_at_first.y = 0.0f;
        // nrm_at_first.z = 0.0f;
      }
      if (nrm.HasNaNs()) {
        Error("normal has nans");
        // nrm.x = 0.0f;
        // nrm.y = 0.0f;
        // nrm.z = 0.0f;
      }
      if (albedo.HasNaNs()) {
        Error("albedo has nans");
        // albedo = Spectrum(0.f);
      }
      if (albedo_at_first.HasNaNs()) {
        Error("albedo at first has nans");
        // albedo_at_first = Spectrum(0.f);
      }
      if (albedo.y() > 101.0f || albedo_at_first.y() > 101.0f) {
          Error("albedo is too high");
          // throw;
      }

      if (Ldiffuse.HasNaNs()) {
          Error("diffuse has nan");
        // Ldiffuse = Spectrum(0.f);
      }
      if (Ldiffuse_indirect.HasNaNs()) {
          Error("diffuse indirect has nan");
        // Ldiffuse_indirect = Spectrum(0.f);
      }
      if (L.HasNaNs()) {
          Error("L  has nan");
        // L = Spectrum(0.f);
      }

      sr->normal_at_first.push_back(nrm_at_first);
      sr->normal.push_back(nrm);
      sr->depth_at_first.push_back(depth_at_first);
      sr->depth.push_back(depth);
      sr->visibility.push_back(isLightVisible ? 1.0 : 0.0);
      sr->hasHit.push_back(1.0);
      sr->albedo.push_back(albedo);
      sr->albedo_at_first.push_back(albedo_at_first);

      sr->radiance_diffuse.push_back(Ldiffuse);
      sr->radiance_diffuse_indirect.push_back(Ldiffuse_indirect);
      sr->radiance_specular.push_back(L - Ldiffuse);
      sr->probabilities.push_back(probabilities);
      sr->light_directions.push_back(light_directions);
      sr->bounce_type.push_back(bounce_type);
    }

    // TODO(mgharbi): nrm vs nrm at first for image data?
    return RadianceQueryRecord(
        L, Ldiffuse, albedo_at_first, nrm_at_first, depth_at_first, isLightVisible, true);
}


PathRendernetIntegrator *CreatePathRendernetSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new PathRendernetIntegrator(maxDepth);
}
