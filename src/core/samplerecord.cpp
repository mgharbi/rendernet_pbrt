#include "samplerecord.h"
#include "pbrt.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <lz4frame.h>

int SampleRecord::version = 20180626;
int SampleRecord::buffer_channels = 15;
int SampleRecord::sample_features = 
  5   // dx, dy, u, v, t
  + 3*2 // rgb*(diffuse + specular)
  + 3  // normals_at_first
  // + 3   // normals
  + 1   // depth_at_first
  // + 1   // depth
  + 1   // visibility
  + 1   // hit
  // + 3   // albedo_at_first
  + 3;  // albedo
int SampleRecord::pixel_features = SampleRecord::buffer_channels*2*2; // 2 references with variance

RadianceQueryRecord::RadianceQueryRecord() {
  count = 0;
  buffer = std::vector<float>(SampleRecord::buffer_channels, 0.0f);
  var_buffer = std::vector<float>(SampleRecord::buffer_channels, 0.0f);
}

RadianceQueryRecord::RadianceQueryRecord(
      Spectrum radiance, Spectrum diffuse, Spectrum albedo, Normal nrm, 
      float depth, bool visibility, bool hasHit) {
  count = 1;
  buffer = std::vector<float>(SampleRecord::buffer_channels, 0.0f);
  var_buffer = std::vector<float>(SampleRecord::buffer_channels, 0.0f);

  radiance = check_radiance(radiance);
  diffuse = check_radiance(diffuse);
  albedo = check_radiance(albedo);

  float rgb_d[3];
  diffuse.ToRGB(rgb_d);
  float rgb_s[3];
  Spectrum specular = radiance-diffuse;
  specular.ToRGB(rgb_s);  // specular
  float rgb_a[3];
  albedo.ToRGB(rgb_a);

  std::copy(rgb_d, rgb_d+3, buffer.begin());
  std::copy(rgb_s, rgb_s+3, buffer.begin() + 3);
  std::copy(rgb_a, rgb_a+3, buffer.begin() + 6);
  buffer[9]  = nrm.x;
  buffer[10] = nrm.y;
  buffer[11] = nrm.z;
  buffer[12] = depth;
  buffer[13] = visibility ? 1.0f : 0.0f;
  buffer[14] = hasHit ? 1.0f : 0.0f;
}

void RadianceQueryRecord::add(const RadianceQueryRecord &other, float rayWeight) {
  if (rayWeight != 1.0f) {
    Error("RayWeight should be 1.0. Not handled by sample saver");
    throw;
  }

  count += 1;
  for (int i = 0; i < (int)buffer.size(); ++i) {
    float mean = buffer[i];
    float delta = other.buffer[i] - mean;
    mean += delta / count;
    buffer[i] = mean;
    float delta2 = other.buffer[i] - mean;
    var_buffer[i] += delta*delta2;
  }
}

Spectrum RadianceQueryRecord::check_radiance(Spectrum &r) {
  if (r.HasNaNs()) {
    Error("Not-a-number radiance value returned "
        "for image sample.  Setting to black.");
    return Spectrum(0.0f);
  } else if (r.y() < 0) {
    Error("Negative luminance value, %f, returned "
        "for image sample.  Setting to black.", r.y());
    return Spectrum(0.0f);
  } else if (isinf(r.y())) {
    Error("Infinite luminance value returned "
        "for image sample.  Setting to black.");
    return Spectrum(0.0f);
  }
  return r;
}

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
  : tile_x(x), tile_y(y), tileSize(tilesize), sample_count(sample_count),
  spp(spp), maxDepth(maxDepth+1),  // We store path 0...maxDepth
  image_width(image_width), image_height(image_height),
  scene_radius(scene_radius), focal_distance(focal_distance),
  aperture_radius(aperture_radius), fov(fov),
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
  hasHit.reserve(tileSize*tileSize*sample_count);
  albedo.reserve(tileSize*tileSize*sample_count);
  albedo_at_first.reserve(tileSize*tileSize*sample_count);
  probabilities.reserve(tileSize*tileSize*sample_count);
  light_directions.reserve(tileSize*tileSize*sample_count);
  bounce_type.reserve(tileSize*tileSize*sample_count);

  // suffix
  image_data.reserve(buffer_channels*2*2); // 2 reference images (+ variance)
  for (int i = 0; i < 2*2*buffer_channels; ++i) {
    image_data.push_back(std::vector<float>());
    image_data[i].reserve(tileSize*tileSize);
  }
}

void SampleRecord::check_sizes() {
  if ((int)pixel_x.size() != tileSize*tileSize*sample_count)
    Error("incorrect pixel_x");
  if ((int)pixel_y.size() != tileSize*tileSize*sample_count)
    Error("incorrect pixel_y");
  if ((int)subpixel_x.size() != tileSize*tileSize*sample_count)
    Error("incorrect subpixel_x");
  if ((int)subpixel_y.size() != tileSize*tileSize*sample_count)
    Error("incorrect subpixel_y");
  if ((int)lens_u.size() != tileSize*tileSize*sample_count)
    Error("incorrect lens_u");
  if ((int)lens_v.size() != tileSize*tileSize*sample_count)
    Error("incorrect lens_v");
  if ((int)time.size() != tileSize*tileSize*sample_count)
    Error("incorrect time");

  if ((int)radiance_diffuse.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_diffuse");
  if ((int)radiance_diffuse_indirect.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_diffuse_indirect");
  if ((int)radiance_specular.size() != tileSize*tileSize*sample_count)
    Error("incorrect radiance_specular");
  if ((int)normal_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect normal_at_first (got %lu expected %d)",
        normal_at_first.size(), tileSize*tileSize*sample_count);
  if ((int)depth_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect depth  at first (got %lu expected %d)",
        depth_at_first.size(), tileSize*tileSize*sample_count);
  if ((int)normal.size() != tileSize*tileSize*sample_count)
    Error("incorrect normal");
  if ((int)depth.size() != tileSize*tileSize*sample_count)
    Error("incorrect depth");
  if ((int)visibility.size() != tileSize*tileSize*sample_count)
    Error("incorrect visibility");
  if ((int)hasHit.size() != tileSize*tileSize*sample_count)
    Error("incorrect hasHit");
  if ((int)albedo.size() != tileSize*tileSize*sample_count)
    Error("incorrect albedo");
  if ((int)albedo_at_first.size() != tileSize*tileSize*sample_count)
    Error("incorrect albedo at first (got %lu expected %d)",
        albedo_at_first.size(), tileSize*tileSize*sample_count);
  if ((int)probabilities.size() != tileSize*tileSize*sample_count)
    Error("incorrect probabilities at first (got %lu expected %d)",
        probabilities.size(), tileSize*tileSize*sample_count);

  for (int i = 0; i < buffer_channels ; ++i) {
    if (image_data[i].size() != tileSize*tileSize)
      Error("incorrect image data");
  }
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

void SampleRecord::write_sample_buffer(
    int sample_id, vector<float> &src, std::ostream &f) {
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

void SampleRecord::write_normal_sample_buffer(
    int sample_id, vector<Normal> &src, std::ostream &f) {
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

void SampleRecord::write_float_path_data
(int sample_id, int count, vector<vector<float> > &src, std::ostream &f) {
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

void SampleRecord::write_bt_sample_buffer(
    int sample_id, vector<vector<uint16_t> > &src, std::ostream &f) {
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

  int nbytes = LZ4F_compressFrame(dst, compsize, cstr, fi.tellp(), NULL);
  if(nbytes == 0) {
    throw "could not compress";
  }

  f.write((char*)&nbytes, sizeof(int));
  f.write(dst, nbytes);
  delete[] dst;

  return nbytes;
}

void SampleRecord::add_image_sample(const RadianceQueryRecord &r, int sampler_idx) {
  // TODO: normalize normals to be unit vectors...
  for (int i = 0; i < buffer_channels; ++i) {
    image_data[sampler_idx*buffer_channels*2 + i].push_back(r.buffer[i]);
    float var = 0.f;
    if (r.count > 1) {
      var = r.var_buffer[i] / (r.count-1);
    }

    // Monte-Carlo variance estimate: 1/n * var
    if (r.count > 0) {
      var /= r.count;
    }

    image_data[sampler_idx*buffer_channels*2 + i + buffer_channels].push_back(var);
  }
}

bool SampleRecord::has_nans() {
  for(int s = 0; s < tileSize*tileSize*sample_count; ++s) {
    if ( isnan(pixel_x[s]) )
      return true;
    if ( isnan(pixel_y[s]) )
      return true;
    if ( isnan(lens_u[s]) )
      return true;
    if ( isnan(lens_v[s]) )
      return true;
    if ( isnan(time[s]) )
      return true;
    if ( radiance_diffuse[s].HasNaNs() )
      return true;
    if ( radiance_diffuse_indirect[s].HasNaNs() )
      return true;
    if ( radiance_specular[s].HasNaNs() )
      return true;
    if ( isnan(normal[s].x) || isnan(normal[s].y) || isnan(normal[s].z))
      return true;
    if (isnan(normal_at_first[s].x) || isnan(normal_at_first[s].y) ||
        isnan(normal_at_first[s].z))
      return true;
    if ( isnan(depth[s]) )
      return true;
    if ( isnan(depth_at_first[s]) )
      return true;
    if ( isnan(visibility[s]) )
      return true;
    if ( isnan(hasHit[s]) )
      return true;
    if ( albedo[s].HasNaNs() )
      return true;
    if ( albedo_at_first[s].HasNaNs() )
      return true;
    for(int pidx=0; pidx < 4*maxDepth; ++pidx) {
      if ( isnan(probabilities[s][pidx]) )
        return true;
    }
    for(int pidx=0; pidx < 2*maxDepth; ++pidx) {
      if ( isnan(light_directions[s][pidx]) )
        return true;
    }
  }
  for(int p = 0; p < tileSize*tileSize; ++p) {
    for (int i = 0; i < buffer_channels ; ++i) {
      if ( isnan(image_data[i][p]) )
        return true;
    }
  }
  return false;
}

void SampleRecord::save(const char* fname) {
  check_sizes();
  if (has_nans() ){
    Error("NaNs in sample record, skipping save.");
    return;
  }

  std::ofstream f(fname, std::ios::binary);

  // Keep distances unnormalized for NFOR and KPCN
  if(!is_kpcn) {
    normalize_distances();
    normalize_probabilities();
  }

  if (sample_count <= 0 || sample_count > spp) {
    Error("saved samples should be higher than 0 and lower than spp %d, got %d.",
        spp, sample_count);
    return;
  }

  // Write metadata header
  {
    f.write((char*)&version, sizeof(int));
    f.write((char*)&tileSize, sizeof(int));
    f.write((char*)&image_width, sizeof(int));
    f.write((char*)&image_height, sizeof(int));
    f.write((char*)&sample_count, sizeof(int));
    f.write((char*)&spp, sizeof(int));
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

  {
    // Write pixel data
    std::stringstream sstream;
    std::ostream_iterator<float> it(sstream);
    for (int i = 0; i < (int)image_data.size(); ++i) {
      sstream.write((char*)image_data[i].data(), tileSize*tileSize*sizeof(float));
    }
    int nb = write_compressed(sstream, f);
  }

  // Write sample data
  for(int sample_id = 0; sample_id < sample_count; ++sample_id) {
    std::stringstream sstream;
    write_sample_buffer(sample_id, subpixel_x, sstream);
    write_sample_buffer(sample_id, subpixel_y, sstream);
    write_sample_buffer(sample_id, lens_u, sstream);
    write_sample_buffer(sample_id, lens_v, sstream);
    write_sample_buffer(sample_id, time, sstream);
    write_rgb_sample_buffer(sample_id, radiance_diffuse, sstream);
    write_rgb_sample_buffer(sample_id, radiance_specular, sstream);
    write_normal_sample_buffer(sample_id, normal_at_first, sstream);
    // write_normal_sample_buffer(sample_id, normal, sstream);
    write_sample_buffer(sample_id, depth_at_first, sstream);
    // write_sample_buffer(sample_id, depth, sstream);
    write_sample_buffer(sample_id, visibility, sstream);
    write_sample_buffer(sample_id, hasHit, sstream);
    write_rgb_sample_buffer(sample_id, albedo_at_first, sstream);
    // write_rgb_sample_buffer(sample_id, albedo, sstream);
    // write_float_path_data(sample_id, 4, probabilities, sstream);
    // write_float_path_data(sample_id, 2, light_directions, sstream);
    // write_bt_sample_buffer(sample_id, bounce_type, sstream);

    write_compressed(sstream, f);
  }
  
  f.close();
}
