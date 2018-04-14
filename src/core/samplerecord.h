#ifndef SAMPLERECORD_H_WVTEISKD
#define SAMPLERECORD_H_WVTEISKD

#include "pbrt.h"
#include "core/spectrum.h"
#include "core/diffgeom.h"
#include "core/transform.h"
#include <vector>
#include <sstream>

using std::vector;

class LightQueryRecord {

public:
  LightQueryRecord(Transform t) 
    : visibility(0.0f), diffuse_lighting(0.0f), world2cam(t), theta(0.f), phi(0.f)
  { memset(pdfs, 0, sizeof(float)*4); };
  Spectrum visibility;
  Spectrum diffuse_lighting;
  Transform world2cam;

  // Spherical coordinates of light_sample
  float theta;
  float phi;
  // Direct/BSDF light sampling pdf
  float pdfs[4]; // (light_pdf, bsdf_pdf)_light, (light_pdf, bsdf_pdf)_bsdf

  void set_angles(Vector wi);
};

class SampleRecord {
public:
  static int version;
  static int sample_features;
  static int pixel_features;

  SampleRecord(
      int x, int y, int tilesize, int sample_count, int spp, int maxDepth,
      int image_width, int image_height, float scene_radius,
      float focal_distance, float aperture_radius, float fov,
      bool useCameraSpaceNormals);
  void save(const char* fname);

  int tile_x;
  int tile_y;
  int tileSize;
  int sample_count;
  int spp;
  int maxDepth;

  int image_width;
  int image_height;
  float scene_radius;
  float focal_distance;
  float aperture_radius;
  float fov;

  bool useCameraSpaceNormals;

  // prefix
  vector<float> pixel_x;
  vector<float> pixel_y;
  vector<float> subpixel_x;
  vector<float> subpixel_y;
  vector<float> lens_u;
  vector<float> lens_v;
  vector<float> time;

  // data
  vector<RGBSpectrum> radiance_diffuse;
  vector<RGBSpectrum> radiance_diffuse_indirect;
  vector<RGBSpectrum> radiance_specular;
  vector<Normal> normal_at_first;
  vector<float> depth_at_first;
  vector<Normal> normal;
  vector<float> depth;
  vector<float> visibility;
  vector<RGBSpectrum> albedo;
  vector<RGBSpectrum> albedo_at_first;
  vector<vector<float> > probabilities;
  vector<vector<float> > light_directions;
  vector<vector<uint16_t> > bounce_type;  // see core/reflection.h

  // suffix
  vector<RGBSpectrum> ground_truth;
  vector<RGBSpectrum> ground_truth_variance;
  vector<RGBSpectrum> lowspp;
  vector<RGBSpectrum> lowspp_variance;

private:
  void check_sizes();

  void normalize_distances();
  void normalize_probabilities();

  void write_rgb_buffer(vector<RGBSpectrum> &src, std::ostream &f, float clamp = -1.0f);

  int write_compressed(std::stringstream &fi, std::ostream &f);

  void write_sample_buffer(int sidx, vector<float> &src, std::ostream &f);
  void write_rgb_sample_buffer(int sidx, vector<RGBSpectrum> &src, std::ostream &f, float clamp = -1.0f);
  void write_normal_sample_buffer(int sidx, vector<Normal> &src, std::ostream &f);
  void write_float_path_data(int sidx, int count, vector<vector<float> > &src, std::ostream &f);
  void write_bt_sample_buffer(int sidx, vector<vector<uint16_t> > &src, std::ostream &f);

};

#endif /* end of include guard: SAMPLERECORD_H_WVTEISKD */
