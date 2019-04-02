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
    Spectrum pathThroughputDiffuse = 1., Ldiffuse = 0. ;
      // , Ldiffuse_indirect = 0.;
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
            // if (bounces > 0)  {
            //   Ldiffuse_indirect += contrib*pathThroughputDiffuse;
            // }
          }
        }

        BSDF *bsdf = isectp->GetBSDF(ray, arena);
        const Point &p = bsdf->dgShading.p;
        const Normal &n = bsdf->dgShading.nn;

        // Characterize the brdf
        bool bsdf_has_diffuse =
            (bsdf->NumComponents(BxDFType(BSDF_DIFFUSE|BSDF_REFLECTION)) > 0);
        bool bsdf_has_nonspecular = bsdf_has_diffuse ||
            (bsdf->NumComponents(BxDFType(BSDF_GLOSSY|BSDF_REFLECTION)) > 0) ||
            (bsdf->NumComponents(BxDFType(BSDF_GLOSSY|BSDF_TRANSMISSION)) > 0);

        // Accumulate path length
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
          // if (bounces > 0)  {
          //   Ldiffuse_indirect += qr.diffuse_lighting*pathThroughputDiffuse;
          // }
        } else if (foundRough) {
          Ldiffuse += contrib*pathThroughputDiffuse;
          // if (bounces > 0)  {
          //   Ldiffuse_indirect += contrib*pathThroughputDiffuse;
          // }
        }

        // Store the lighting directions
        std::copy(qr.pdfs, qr.pdfs+4, probabilities.begin()+4*bounces);
        light_directions[2*bounces + 0] = qr.theta;
        light_directions[2*bounces + 1] = qr.phi;

        // Sample BSDF to get new path direction

        // Get _outgoingBSDFSample_ for sampling new path direction
        BSDFSample outgoingBSDFSample;
        if (bounces < SAMPLE_DEPTH) {
          outgoingBSDFSample = BSDFSample(sample, pathSampleOffsets[bounces], 0);
        }
        else {
          outgoingBSDFSample = BSDFSample(rng);
        }
        Vector wi;
        float pdf;
        BxDFType flags;
        Spectrum f = bsdf->Sample_f(wo, &wi, outgoingBSDFSample, &pdf,
                                    BSDF_ALL, &flags);
        bounce_type[bounces] = flags;
        Spectrum currAlbedo = bsdf->K();

        // If the brdf has a diffuse component and we have not found the first
        // rough bounce: this is it. we found our first
        // diffuse interaction. The path is no longer purely specular.
        bool isFirstRough = false;
        if (!foundRough && bsdf_has_diffuse) {
          foundRough = true;
          isFirstRough = true;
        } 

        bool isFirstNonSpecular = false;
        // If the brdf has a nonspecular component and we have not found any
        // yet: this it is.
        if (!foundNonSpecular && bsdf_has_nonspecular) {
          foundNonSpecular = true;
          isFirstNonSpecular = true;
        } 

        // Record depth, normal, albedo, visibility at first bounce
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

          // Flag whether light is directly visible (i.e. at first bounce)
          isLightVisible = isLightVisible || qr.isLightVisible;
        }

        // record value at first nonspecular 
        if (!recordedOutputValues && isFirstNonSpecular) {
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

        // After the first rough bounce, the path is no longer purely specular,
        // we accumulate radiance in the diffuse component
        // TODO(mgharbi): this looks odd, should it be at all diffuse bounce? or after the first
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
                    // if (bounces > 0)  {
                    //   Ldiffuse_indirect += contrib*pathThroughputDiffuse;
                    // }
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

    // Check fro NaNs
    if (nrm_at_first.HasNaNs()) {
      Error("normal first has nans");
    }
    if (nrm.HasNaNs()) {
      Error("normal has nans");
    }
    if (albedo.HasNaNs()) {
      Error("albedo has nans");
    }
    if (albedo_at_first.HasNaNs()) {
      Error("albedo at first has nans");
    }
    if (albedo.y() > 101.0f || albedo_at_first.y() > 101.0f) {
        Error("albedo is too high");
    }
    if (Ldiffuse.HasNaNs()) {
        Error("diffuse has nan");
    }
    // if (Ldiffuse_indirect.HasNaNs()) {
    //     Error("diffuse indirect has nan");
    // }
    if (L.HasNaNs()) {
        Error("L  has nan");
    }

    if (sr) {
      // Store decomposed radiance
      sr->radiance_diffuse.push_back(Ldiffuse);
      // sr->radiance_diffuse_indirect.push_back(Ldiffuse_indirect);
      sr->radiance_specular.push_back(L - Ldiffuse);

      // Store features at first bounce
      sr->normal_at_first.push_back(nrm_at_first);
      sr->depth_at_first.push_back(depth_at_first);
      sr->albedo_at_first.push_back(albedo_at_first);

      // Store other features
      sr->normal.push_back(nrm);
      sr->depth.push_back(depth);
      sr->visibility.push_back(isLightVisible ? 1.0 : 0.0);
      sr->hasHit.push_back(1.0);
      sr->albedo.push_back(albedo);
      sr->probabilities.push_back(probabilities);
      sr->light_directions.push_back(light_directions);
      sr->bounce_type.push_back(bounce_type);
    }

    return RadianceQueryRecord(
        L, Ldiffuse, albedo, nrm, depth, isLightVisible, true);
}


PathRendernetIntegrator *CreatePathRendernetSurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new PathRendernetIntegrator(maxDepth);
}
