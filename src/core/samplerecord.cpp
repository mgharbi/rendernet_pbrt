#include "samplerecord.h"
#include "pbrt.h"
#include <fstream>
#include <iostream>
#include <lz4frame.h>

int SampleRecord::version = 20180330;
int SampleRecord::sample_features = 
  5   // dx, dy, u, v, t
  + 3*2 // rgb*(diffuse + specular)
  + 3   // normals_at_first
  + 3   // normals
  + 1   // depth_at_first
  + 1   // depth
  + 1   // visibility
  + 3   // albedo_at_first
  + 3;  // albedo
int SampleRecord::pixel_features = 3*4; // rgb * (gt + std + lowspp + lowspp_std)

void LightQueryRecord::set_angles(Vector wi) {
    // spherical coordinates of light direction
    // in camera space
    wi = world2cam(wi);

    float nrm = sqrt(wi.x*wi.x + wi.y*wi.y);
    
    if(nrm == 0.0f) {
      theta = 0.0f;
    } else {
      theta = atan2(wi.y, wi.x);
    }

    if (nrm == 0.0f && wi.z == 0) {
      phi = 0.0f;
    } else {
      phi = atan2(nrm, wi.z);
    }
    theta /= M_PI;
    phi /= M_PI;
}

SampleRecord::SampleRecord(
    int x, int y, int tilesize, int sample_count, int spp, int maxDepth,
    int image_width, int image_height, float scene_radius,
    float focal_distance, float aperture_radius, float fov, bool useCameraSpaceNormals)
  : tile_x(x), tile_y(y), tileSize(tilesize), sample_count(sample_count), spp(spp), maxDepth(maxDepth+1),  // We store path 0...maxDepth
  image_width(image_width), image_height(image_height), scene_radius(scene_radius),
  focal_distance(focal_distance), aperture_radius(aperture_radius), fov(fov),
  useCameraSpaceNormals(useCameraSpaceNormals), is_kpcn(false)
{
  // prefix
  pixel_x.reserve(tileSize*tileSize*sample_count);
  pixel_y.reserve(tileSize*tileSize*sample_count);
  subpixel_x.reserve(tileSize*tileSize*sample_count);
  subpixel_y.reserve(tileSize*tileSize*sample_count);
  lens_u.reserve(tileSize*tileSize*sample_count);
  lens_v.reserve(tileSize*tileSize*sample_count);
  time.reserve(tileSize*tileSize*sample_count);

  // sample data
  radiance_diffuse.reserve(tileSize*tileSize*sample_count);
  radiance_diffuse_indirect.reserve(tileSize*tileSize*sample_count);
  radiance_specular.reserve(tileSize*tileSize*sample_count);
  normal_at_first.reserve(tileSize*tileSize*sample_count);
  depth_at_first.reserve(tileSize*tileSize*sample_count);
  normal.reserve(tileSize*tileSize*sample_count);
  depth.reserve(tileSize*tileSize*sample_count);
  visibility.reserve(tileSize*tileSize*sample_count);
  albedo.reserve(tileSize*tileSize*sample_count);
  albedo_at_first.reserve(tileSize*tileSize*sample_count);
  probabilities.reserve(tileSize*tileSize*sample_count);
  light_directions.reserve(tileSize*tileSize*sample_count);
  bounce_type.reserve(tileSize*tileSize*sample_count);

  // suffix
  ground_truth.reserve(tileSize*tileSize);
  ground_truth_variance.reserve(tileSize*tileSize);
  lowspp.reserve(tileSize*tileSize);
  lowspp_variance.reserve(tileSize*tileSize);
}

void SampleRecord::check_sizes() {
  if (pixel_x.size() != tileSize*tileSize*sample_count)
    Error("incorrect pixel_x");
  if (pixel_y.size() != tileSize*tileSize*sample_count)
    Error("incorrect pixel_y");
  if (subpixel_x.size() != tileSize*tileSize*sample_count)
    Error("incorrect subpixel_x");
  if (subpixel_y.size() != tileSize*tileSize*sample_count)
    Error("incorrect subpixel_y");
  if (lens_u.size() != tileSize*tileSize*sample_count)
    Error("incorrect lens_u");
  if (lens_v.size() != tileSize*tileSize*sample_count)
    Error("incorrect lens_v");
  if (time.size() != tileSize*tileSize*sample_count)
    Error("incorrect time");

  if (radiance_diffuse.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_diffuse");
  if (radiance_diffuse_indirect.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_diffuse_indirect");
  if (radiance_specular.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_specular");
  if (normal_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect normal_at_first (got %d expected %d)", normal_at_first.size(), tileSize*tileSize*sample_count);
  if (depth_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect depth  at first (got %d expected %d)", depth_at_first.size(), tileSize*tileSize*sample_count);
  if (normal.size() != tileSize*tileSize*sample_count)
    Error("incorrect normal");
  if (depth.size() != tileSize*tileSize*sample_count)
    Error("incorrect depth");
  if (visibility.size() != tileSize*tileSize*sample_count)
    Error("incorrect visibility");
  if (albedo.size() != tileSize*tileSize*sample_count)
    Error("incorrect albedo");
  if (albedo_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect albedo at first (got %d expected %d)", albedo_at_first.size(), tileSize*tileSize*sample_count);
  if (probabilities.size() != tileSize*tileSize*sample_count)
    Error("incorrect probabilities at first (got %d expected %d)", probabilities.size(), tileSize*tileSize*sample_count);

  if (ground_truth.size() != tileSize*tileSize)
    Error("incorrect ground_truth");
  if (ground_truth_variance.size() != tileSize*tileSize)
    Error("incorrect ground_truth_variance");
  if (lowspp.size() != tileSize*tileSize)
    Error("incorrect lowspp");
  if (lowspp_variance.size() != tileSize*tileSize)
    Error("incorrect lowspp_variance");
}

void SampleRecord::write_rgb_buffer(
    vector<RGBSpectrum> &src, std::ostream &f, float clamp) {
  int npixels = tileSize*tileSize;
  float* tmp = new float[npixels*3];
  float rgb[3] = {0};
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id) {
    src[pixel_id].ToRGB(rgb);
    for(int i = 0; i < 3; ++i)  {
      if(clamp > 0.0f) {
        rgb[i] = min(rgb[i], clamp);
      }
      tmp[i*npixels + pixel_id] = rgb[i];
    }
  }

  f.write((char*) tmp, 3*npixels*sizeof(float));
  delete[] tmp;
}

void SampleRecord::write_sample_buffer(int sample_id, vector<float> &src, std::ostream &f) {
  int npixels = tileSize*tileSize;
  float* tmp = new float[npixels];
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id)
  {
    int src_idx = sample_id + sample_count*pixel_id;
    int tgt_idx = pixel_id;
    tmp[tgt_idx] = src[src_idx];
  }

  f.write((char*) tmp, npixels*sizeof(float));
  delete[] tmp;
}

void SampleRecord::write_rgb_sample_buffer(
    int sample_id, vector<RGBSpectrum> &src, std::ostream &f, float clamp) {
  int npixels = tileSize*tileSize;
  float* tmp = new float[npixels*3];
  float rgb[3] = {0};
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id)
  {
    int src_idx = sample_id + sample_count*pixel_id;
    src[src_idx].ToRGB(rgb);
    for(int i = 0; i < 3; ++i)  {
      if(clamp > 0.0f) {
        rgb[i] = min(rgb[i], clamp);
      }
      tmp[i*npixels + pixel_id] = rgb[i];
    }
  }

  f.write((char*) tmp, 3*npixels*sizeof(float));
  delete[] tmp;
}

void SampleRecord::write_normal_sample_buffer(int sample_id, vector<Normal> &src, std::ostream &f) {
  int npixels = tileSize*tileSize;
  float* tmp = new float[npixels*3];
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id)
  {
    int src_idx = sample_id + sample_count*pixel_id;
    Normal n = src[src_idx];
    tmp[0*npixels + pixel_id] = n.x;
    tmp[1*npixels + pixel_id] = n.y;
    tmp[2*npixels + pixel_id] = n.z;
  }

  f.write((char*) tmp, 3*npixels*sizeof(float));
  delete[] tmp;
}

void SampleRecord::write_float_path_data(int sample_id, int count, vector<vector<float> > &src, std::ostream &f) {
  int npixels = tileSize*tileSize;
  int n_proba = count*maxDepth;
  float* tmp = new float[npixels*n_proba];
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id)
  {
    int src_idx = sample_id + sample_count*pixel_id;
    vector<float> p = src[src_idx];
    for(int i = 0; i < n_proba; ++ i) {
      int tgt_idx = i*npixels + pixel_id;
      tmp[tgt_idx] = p[i];
    }
  }

  f.write((char*) tmp, n_proba*npixels*sizeof(float));
  delete[] tmp;
}

void SampleRecord::write_bt_sample_buffer(int sample_id, vector<vector<uint16_t> > &src, std::ostream &f) {
  int npixels = tileSize*tileSize;
  int n = maxDepth;
  uint16_t* tmp = new uint16_t[npixels*n];
  // convert to contiguous pixels
  for(int pixel_id = 0; pixel_id < npixels; ++pixel_id)
  {
    int src_idx = sample_id + sample_count*pixel_id;
    vector<uint16_t> p = src[src_idx];
    for(int i = 0; i < n; ++ i) {
      int tgt_idx = i*npixels + pixel_id;
      tmp[tgt_idx] = p[i];
    }
  }

  f.write((char*) tmp, n*npixels*sizeof(uint16_t));
  delete[] tmp;
}

void SampleRecord::normalize_distances() {
  int npixels = tileSize*tileSize;
  float normalizer = scene_radius > 0.0f ? 1.0f/(10.0f*scene_radius) : 1.0f;
  for(int i = 0; i < npixels*sample_count; ++i) {
    depth[i] *= normalizer; 
    depth_at_first[i] *= normalizer; 
    lens_u[i] *= normalizer; 
    lens_v[i] *= normalizer; 
  }
  focal_distance *= normalizer;
  // If aperture is rescaled, u, v should be as well
  aperture_radius *= normalizer;

}

void SampleRecord::normalize_probabilities() {
  const int npixels = tileSize*tileSize;
  const int n_proba = 4*maxDepth;
  const float eps = 1e-8f;
  const float nrm = 30.0f;
  for(int i = 0; i < npixels*sample_count; ++i)
  for(int j = 0; j < n_proba; ++j) {
    float p = probabilities[i][j];
    probabilities[i][j] = log(max(p, 0.0f) + eps) / nrm;
  }
}

int SampleRecord::write_compressed(std::stringstream &fi, std::ostream &f) {
  int compsize = LZ4F_compressFrameBound(fi.tellp(), NULL);
  if(compsize == 0) {
    throw "could not compress";
  }
  const std::string tmp = fi.str();
  const char* cstr = tmp.c_str();
  char *dst = new char[compsize];

  // int nbytes = LZ4_compress_limitedOutput(cstr, dst, fi.tellp(), compsize);
  int nbytes = LZ4F_compressFrame(dst, compsize, cstr, fi.tellp(), NULL);
  if(nbytes == 0) {
    throw "could not compress";
  }

  f.write((char*)&nbytes, sizeof(int));
  f.write(dst, nbytes);
  delete[] dst;

  return nbytes;
}

void SampleRecord::save(const char* fname) {
  std::ofstream f(fname, std::ios::binary);

  check_sizes();

  // Keep distances unnormalized for NFOR and KPCN
  if(!is_kpcn) {
    normalize_distances();
    normalize_probabilities();
  }

  if (sample_count <= 0 || sample_count > spp)
    Error("saved samples should be higher than 0 and lower than spp %d, got %d.",
        spp, sample_count);

  // Write metadata header
  {
    f.write((char*)&version, sizeof(int));
    f.write((char*)&tileSize, sizeof(int));
    f.write((char*)&image_width, sizeof(int));
    f.write((char*)&image_height, sizeof(int));
    f.write((char*)&sample_count, sizeof(int));
    f.write((char*)&sample_features, sizeof(int));
    f.write((char*)&pixel_features, sizeof(int));
    f.write((char*)&maxDepth, sizeof(int));
  }

  // Write globals
  {
    f.write((char*)&focal_distance, sizeof(float));
    f.write((char*)&aperture_radius, sizeof(float));
    float norm_fov = fov / 100.0f;
    f.write((char*)&norm_fov, sizeof(float));
    f.write((char*)&scene_radius, sizeof(float));
  }
  
  // Write tile coordinates 
  {
    f.write((char*)&tile_x, sizeof(int));
    f.write((char*)&tile_y, sizeof(int));
  }

  std::stringstream sstream;

  // Write pixel data
  write_rgb_buffer(ground_truth, sstream);
  write_rgb_buffer(ground_truth_variance, sstream);
  write_rgb_buffer(lowspp, sstream);
  write_rgb_buffer(lowspp_variance, sstream);
  int nb = write_compressed(sstream, f);

  // Write sample data
  for(int sample_id = 0; sample_id < sample_count; ++sample_id) {
    sstream.seekp(0);
    // write_sample_buffer(sample_id, pixel_x, sstream);
    // write_sample_buffer(sample_id, pixel_y, sstream);
    write_sample_buffer(sample_id, subpixel_x, sstream);
    write_sample_buffer(sample_id, subpixel_y, sstream);
    write_sample_buffer(sample_id, lens_u, sstream);
    write_sample_buffer(sample_id, lens_v, sstream);
    write_sample_buffer(sample_id, time, sstream);
    write_rgb_sample_buffer(sample_id, radiance_diffuse, sstream);
    write_rgb_sample_buffer(sample_id, radiance_specular, sstream);
    // write_rgb_sample_buffer(sample_id, radiance_diffuse_indirect, sstream);
    write_normal_sample_buffer(sample_id, normal_at_first, sstream);
    write_normal_sample_buffer(sample_id, normal, sstream);
    write_sample_buffer(sample_id, depth_at_first, sstream);
    write_sample_buffer(sample_id, depth, sstream);
    write_sample_buffer(sample_id, visibility, sstream);
    write_rgb_sample_buffer(sample_id, albedo_at_first, sstream);
    write_rgb_sample_buffer(sample_id, albedo, sstream);
    // write_float_path_data(sample_id, 4, probabilities, sstream);
    write_float_path_data(sample_id, 2, light_directions, sstream);
    write_bt_sample_buffer(sample_id, bounce_type, sstream);

    nb = write_compressed(sstream, f);
  }
  
  f.close();
}
