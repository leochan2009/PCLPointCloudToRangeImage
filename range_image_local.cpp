/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010-2012, Willow Garage, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

//#include <cstddef>
#include <iostream>
#include <cmath>
#include <set>
//#include <pcl/common/eigen.h>
//#include <pcl/common/transformation_from_correspondences.h>
#include "range_image_local.h"
//#include "range_image_local.hpp"  // Definitions of templated and inline functions


bool RangeImageLocal::debug = false;
int RangeImageLocal::max_no_of_threads = 1;
const int RangeImageLocal::lookup_table_size = 20001;
std::vector<float> RangeImageLocal::asin_lookup_table;
std::vector<float> RangeImageLocal::atan_lookup_table;
std::vector<float> RangeImageLocal::cos_lookup_table;

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::createLookupTables ()
{
  if (!asin_lookup_table.empty ())
    return;
  
  asin_lookup_table.resize (lookup_table_size);
  for (int i=0; i<lookup_table_size; ++i) {
    float value = static_cast<float> (i-(lookup_table_size-1)/2)/static_cast<float> ((lookup_table_size-1)/2);
    asin_lookup_table[i] = asinf (value);
  }
  
  atan_lookup_table.resize (lookup_table_size);
  for (int i=0; i<lookup_table_size; ++i)
  {
    float value = static_cast<float> (i-(lookup_table_size-1)/2)/static_cast<float> ((lookup_table_size-1)/2);
    atan_lookup_table[i] = atanf (value);
  }
  
  cos_lookup_table.resize (lookup_table_size);
  
  for (int i = 0; i < lookup_table_size; ++i)
  {
    float value = static_cast<float> (i) * 2.0f * static_cast<float> (M_PI) / static_cast<float> (lookup_table_size-1);
    cos_lookup_table[i] = cosf (value);
  }
}



/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getCoordinateFrameTransformation (RangeImageLocal::CoordinateFrame coordinate_frame,
                                              Eigen::Affine3f& transformation)
{
  switch (coordinate_frame)
  {
    case LASER_FRAME:
      transformation (0,0)= 0.0f; transformation (0,1)= 0.0f; transformation (0,2)=1.0f; transformation (0,3)=0.0f;
      transformation (1,0)=-1.0f; transformation (1,1)= 0.0f; transformation (1,2)=0.0f; transformation (1,3)=0.0f;
      transformation (2,0)= 0.0f; transformation (2,1)=-1.0f; transformation (2,2)=0.0f; transformation (2,3)=0.0f;
      transformation (3,0)= 0.0f; transformation (3,1)= 0.0f; transformation (3,2)=0.0f; transformation (3,3)=1.0f;
      break;
    case CAMERA_FRAME:
    default:
      transformation.setIdentity ();
      break;
  }
}

/////////////////////////////////////////////////////////////////////////
RangeImageLocal::RangeImageLocal () :
RangeImageLocal::BaseClass (),
to_range_image_system_ (Eigen::Affine3f::Identity ()),
to_world_system_ (Eigen::Affine3f::Identity ()),
angular_resolution_x_ (0), angular_resolution_y_ (0),
angular_resolution_x_reciprocal_ (0), angular_resolution_y_reciprocal_ (0),
image_offset_x_ (0), image_offset_y_ (0),
unobserved_point ()
{
  createLookupTables ();
  reset ();
  unobserved_point.x = unobserved_point.y = unobserved_point.z = std::numeric_limits<float>::quiet_NaN ();
  unobserved_point.range = -std::numeric_limits<float>::infinity ();
}

/////////////////////////////////////////////////////////////////////////
RangeImageLocal::~RangeImageLocal ()
{
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::reset ()
{
  is_dense = true;
  width = height = 0;
  points.clear ();
  to_range_image_system_.setIdentity ();
  to_world_system_.setIdentity ();
  setAngularResolution (deg2radLocal (0.5f));
  image_offset_x_ = image_offset_y_ = 0;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::createEmpty (float angular_resolution, const Eigen::Affine3f& sensor_pose,
                         RangeImageLocal::CoordinateFrame coordinate_frame, float angle_width, float angle_height)
{
  createEmpty (angular_resolution, angular_resolution, sensor_pose, coordinate_frame, angle_width, angle_height);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::createEmpty (float angular_resolution_x, float angular_resolution_y, const Eigen::Affine3f& sensor_pose,
                         RangeImageLocal::CoordinateFrame coordinate_frame, float angle_width, float angle_height)
{
  setAngularResolution(angular_resolution_x, angular_resolution_y);
  width  = static_cast<uint32_t> (pcl_lrint (floor (angle_width * angular_resolution_x_reciprocal_)));
  height = static_cast<uint32_t> (pcl_lrint (floor (angle_height * angular_resolution_y_reciprocal_)));
  
  int full_width  = static_cast<int> (pcl_lrint (floor (deg2radLocal (360.0f)*angular_resolution_x_reciprocal_))),
  full_height = static_cast<int> (pcl_lrint (floor (deg2radLocal (180.0f)*angular_resolution_y_reciprocal_)));
  image_offset_x_ = (full_width-width)/2;
  image_offset_y_ = (full_height-height)/2;
  is_dense = false;
  getCoordinateFrameTransformation (coordinate_frame, to_world_system_);
  to_world_system_ = sensor_pose * to_world_system_;
  to_range_image_system_ = to_world_system_.inverse (Eigen::Isometry);
  unsigned int size = width*height;
  points.clear ();
  points.resize (size, unobserved_point);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::cropImage (int borderSize, int top, int right, int bottom, int left) {
  //MEASURE_FUNCTION_TIME;
  
  bool topIsDone=true, rightIsDone=true, bottomIsDone=true, leftIsDone=true;
  if (top < 0) {
    top=-1;
    topIsDone=false;
  }
  if (right < 0)
  {
    right= static_cast<int> (width);
    rightIsDone=false;
  }
  if (bottom < 0)
  {
    bottom= static_cast<int> (height);
    bottomIsDone=false;
  }
  if (left < 0)
  {
    left = -1;
    leftIsDone = false;
  }
  
  // Find top border
  while (!topIsDone && top<=bottom)
  {
    ++top;
    int lineStart = top*width;
    int min_x=std::max(0, left), max_x=std::min(static_cast<int> (width)-1, right);
    for (int x=min_x; x<=max_x && !topIsDone; ++x)
      if (pcl_isfinite (points[lineStart + x].range))
        topIsDone = true;
  }
  // Check if range image is empty
  if (top >= static_cast<int> (height))
  {
    points.clear ();
    width = height = 0;
    return;
  }
  // Find right border
  while (!rightIsDone)
  {
    --right;
    int min_y=std::max(0, top), max_y=std::min(static_cast<int> (height)-1, bottom);
    for (int y=min_y; y<=max_y && !rightIsDone; ++y)
      if (pcl_isfinite (points[y*width + right].range))
        rightIsDone = true;
  }
  // Find bottom border
  while (!bottomIsDone)
  {
    --bottom;
    int lineStart = bottom*width;
    int min_x=std::max(0, left), max_x=std::min(static_cast<int> (width)-1, right);
    for (int x=min_x; x<=max_x && !bottomIsDone; ++x)
      if (pcl_isfinite (points[lineStart + x].range))
        bottomIsDone = true;
  }
  // Find left border
  while (!leftIsDone)
  {
    ++left;
    int min_y=std::max(0, top), max_y=std::min(static_cast<int> (height)-1, bottom);
    for (int y=min_y; y<=max_y && !leftIsDone; ++y)
      if (pcl_isfinite (points[y*width + left].range))
        leftIsDone = true;
  }
  left-=borderSize; top-=borderSize; right+=borderSize; bottom+=borderSize;
  
  // Create copy without copying the old points - vector::swap only copies a few pointers, not the content
  PointCloud<PointWithRange>::VectorType tmpPoints;
  points.swap (tmpPoints);
  RangeImageLocal oldRangeImageLocal = *this;
  tmpPoints.swap (oldRangeImageLocal.points);
  
  width = right-left+1; height = bottom-top+1;
  image_offset_x_ = left+oldRangeImageLocal.image_offset_x_;
  image_offset_y_ = top+oldRangeImageLocal.image_offset_y_;
  points.resize (width*height);
  
  //std::cout << oldRangeImageLocal.width<<"x"<<oldRangeImageLocal.height<<" -> "<<width<<"x"<<height<<"\n";
  
  // Copy points
  for (int y=0, oldY=top; y< static_cast<int> (height); ++y,++oldY)
  {
    for (int x=0, oldX=left; x< static_cast<int> (width); ++x,++oldX)
    {
      PointWithRange& currentPoint = points[y*width + x];
      if (oldX<0 || oldX>= static_cast<int> (oldRangeImageLocal.width) || oldY<0 || oldY>= static_cast<int> (oldRangeImageLocal.height))
      {
        currentPoint = unobserved_point;
        continue;
      }
      currentPoint = oldRangeImageLocal.points[oldY*oldRangeImageLocal.width + oldX];
    }
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::recalculate3DPointPositions ()
{
  for (int y = 0; y < static_cast<int> (height); ++y)
  {
    for (int x = 0; x < static_cast<int> (width); ++x)
    {
      PointWithRange& point = points[y*width + x];
      if (!pcl_isinf (point.range))
        calculate3DPoint (static_cast<float> (x), static_cast<float> (y), point.range, point);
    }
  }
}

/////////////////////////////////////////////////////////////////////////
float*
RangeImageLocal::getRangesArray () const
{
  int arraySize = width * height;
  float* ranges = new float[arraySize];
  for (int i=0; i<arraySize; ++i)
    ranges[i] = points[i].range;
  return ranges;
}


/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getIntegralImage (float*& integral_image, int*& valid_points_num_image) const
{
  integral_image = new float[width*height];
  float* integral_image_ptr = integral_image;
  valid_points_num_image = new int[width*height];
  int* valid_points_num_image_ptr = valid_points_num_image;
  for (int y = 0; y < static_cast<int> (height); ++y)
  {
    for (int x = 0; x < static_cast<int> (width); ++x)
    {
      float& integral_pixel = * (integral_image_ptr++);
      integral_pixel = getPoint (x, y).range;
      int& valid_points_num = * (valid_points_num_image_ptr++);
      valid_points_num = 1;
      if (pcl_isinf (integral_pixel))
      {
        integral_pixel = 0.0f;
        valid_points_num = 0;
      }
      float left_value=0, top_left_value=0, top_value=0;
      int left_valid_points=0, top_left_valid_points=0, top_valid_points=0;
      if (x>0)
      {
        left_value = integral_image[y*width+x-1];
        left_valid_points = valid_points_num_image[y*width+x-1];
        if (y>0)
        {
          top_left_value = integral_image[ (y-1)*width+x-1];
          top_left_valid_points = valid_points_num_image[ (y-1)*width+x-1];
        }
      }
      if (y>0)
      {
        top_value = integral_image[ (y-1)*width+x];
        top_valid_points = valid_points_num_image[ (y-1)*width+x];
      }
      
      integral_pixel += left_value + top_value - top_left_value;
      valid_points_num += left_valid_points + top_valid_points - top_left_valid_points;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::setUnseenToMaxRange ()
{
  for (unsigned int i=0; i<points.size (); ++i)
    if (pcl_isinf (points[i].range))
      points[i].range = std::numeric_limits<float>::infinity ();
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getHalfImage (RangeImageLocal& half_image) const
{
  half_image.setAngularResolution (2.0f*angular_resolution_x_, 2.0f*angular_resolution_y_);
  half_image.image_offset_x_ = image_offset_x_/2;
  half_image.image_offset_y_ = image_offset_y_/2;
  half_image.width  = width/2;
  half_image.height = height/2;
  half_image.is_dense = is_dense;
  half_image.points.clear ();
  half_image.points.resize (half_image.width*half_image.height);
  
  int src_start_x = 2*half_image.image_offset_x_ - image_offset_x_,
  src_start_y = 2*half_image.image_offset_y_ - image_offset_y_;
  
  for (int dst_y=0; dst_y < static_cast<int> (half_image.height); ++dst_y)
  {
    for (int dst_x=0; dst_x < static_cast<int> (half_image.width); ++dst_x)
    {
      PointWithRange& dst_point = half_image.getPoint (dst_x, dst_y);
      dst_point=unobserved_point;
      int src_x_min = src_start_x + 2*dst_x,
      src_x_max = src_x_min + 1,
      src_y_min = src_start_y + 2*dst_y,
      src_y_max = src_y_min + 1;
      for (int src_x=src_x_min; src_x<=src_x_max; ++src_x)
      {
        for (int src_y=src_y_min; src_y<=src_y_max; ++src_y)
        {
          if (!isObserved (src_x, src_y))
            continue;
          const PointWithRange& src_point = getPoint (src_x, src_y);
          if (pcl_isfinite (dst_point.range) && src_point.range > dst_point.range)
            continue;
          dst_point = src_point;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getSubImage (int sub_image_image_offset_x, int sub_image_image_offset_y, int sub_image_width,
                         int sub_image_height, int combine_pixels, RangeImageLocal& sub_image) const
{
  sub_image.setAngularResolution (static_cast<float> (combine_pixels)*angular_resolution_x_,
                                  static_cast<float> (combine_pixels)*angular_resolution_y_);
  sub_image.image_offset_x_ = sub_image_image_offset_x;
  sub_image.image_offset_y_ = sub_image_image_offset_y;
  sub_image.width = sub_image_width;
  sub_image.height = sub_image_height;
  sub_image.is_dense = is_dense;
  sub_image.points.clear ();
  sub_image.points.resize (sub_image.width*sub_image.height);
  
  int src_start_x = combine_pixels*sub_image.image_offset_x_ - image_offset_x_,
  src_start_y = combine_pixels*sub_image.image_offset_y_ - image_offset_y_;
  
  for (int dst_y=0; dst_y < static_cast<int> (sub_image.height); ++dst_y)
  {
    for (int dst_x=0; dst_x < static_cast<int> (sub_image.width); ++dst_x)
    {
      PointWithRange& dst_point = sub_image.getPoint (dst_x, dst_y);
      dst_point=unobserved_point;
      int src_x_min = src_start_x + combine_pixels*dst_x,
      src_x_max = src_x_min + combine_pixels-1,
      src_y_min = src_start_y + combine_pixels*dst_y,
      src_y_max = src_y_min + combine_pixels-1;
      for (int src_x=src_x_min; src_x<=src_x_max; ++src_x)
      {
        for (int src_y=src_y_min; src_y<=src_y_max; ++src_y)
        {
          if (!isInImage (src_x, src_y))
            continue;
          const PointWithRange& src_point = getPoint (src_x, src_y);
          if (pcl_isfinite (dst_point.range) && src_point.range > dst_point.range)
            continue;
          dst_point = src_point;
        }
      }
    }
  }
}


/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getMinMaxRanges (float& min_range, float& max_range) const
{
  min_range = std::numeric_limits<float>::infinity ();
  max_range = -std::numeric_limits<float>::infinity ();
  for (unsigned int i=0; i<points.size (); ++i)
  {
    float range = points[i].range;
    if (!pcl_isfinite (range))
      continue;
    min_range = (std::min) (min_range, range);
    max_range = (std::max) (max_range, range);
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::change3dPointsToLocalCoordinateFrame ()
{
  to_world_system_.setIdentity ();
  to_range_image_system_.setIdentity ();
  recalculate3DPointPositions ();
}

/////////////////////////////////////////////////////////////////////////
float*
RangeImageLocal::getInterpolatedSurfaceProjection (const Eigen::Affine3f& pose, int pixel_size, float world_size) const
{
  float max_dist = 0.5f*world_size,
  cell_size = world_size/float (pixel_size);
  float world2cell_factor = 1.0f/cell_size,
  world2cell_offset = 0.5f*float (pixel_size)-0.5f;
  float cell2world_factor = cell_size,
  cell2world_offset = -max_dist + 0.5f*cell_size;
  Eigen::Affine3f inverse_pose = pose.inverse (Eigen::Isometry);
  
  int no_of_pixels = pixel_size*pixel_size;
  float* surface_patch = new float[no_of_pixels];
  SET_ARRAY (surface_patch, -std::numeric_limits<float>::infinity (), no_of_pixels);
  
  Eigen::Vector3f position = inverse_pose.translation ();
  int middle_x, middle_y;
  getImagePoint (position, middle_x, middle_y);
  int min_search_radius = 2;
  bool still_in_range = true;
  
  for (int radius=0;  still_in_range;  ++radius)
  {
    int x=middle_x-radius-1, y=middle_y-radius;  // Top left - 1
    still_in_range = radius<min_search_radius;
    for (int i=0; i<8*radius || (radius==0&&i==0); ++i)
    {
      if (i<=2*radius) ++x; else if (i<=4*radius) ++y; else if (i<=6*radius) --x; else --y;
      
      Eigen::Vector3f point1, point2, point3;
      if (!isValid (x,y) || !isValid (x+1,y+1))
        continue;
      getPoint (x, y, point1);
      point1 = pose*point1;
      if (fabs (point1[2]) > max_dist)
        continue;
      
      getPoint (x+1, y+1, point2);
      point2 = pose*point2;
      if (fabs (point2[2]) > max_dist)
        continue;
      
      for (int triangle_idx=0; triangle_idx<=1; ++triangle_idx)
      {
        if (triangle_idx==0)  // First triangle
        {
          if (!isValid (x,y+1))
            continue;
          getPoint (x, y+1, point3);
        }
        else  // Second triangle
        {
          if (!isValid (x+1,y))
            continue;
          getPoint (x+1, y, point3);
        }
        point3 = pose*point3;
        if (fabs (point3[2]) > max_dist)
          continue;
        
        // Are all the points either left, right, on top or below the bottom of the surface patch?
        if ( (point1[0] < -max_dist  &&  point2[0] < -max_dist  &&  point3[0] < -max_dist) ||
            (point1[0] >  max_dist  &&  point2[0] >  max_dist  &&  point3[0] >  max_dist) ||
            (point1[1] < -max_dist  &&  point2[1] < -max_dist  &&  point3[1] < -max_dist) ||
            (point1[1] >  max_dist  &&  point2[1] >  max_dist  &&  point3[1] >  max_dist))
        {
          continue;
        }
        
        still_in_range = true;
        
        // Now we have a valid triangle (point1, point2, point3) in our new patch
        float cell1_x = world2cell_factor*point1[0] + world2cell_offset,
        cell1_y = world2cell_factor*point1[1] + world2cell_offset,
        cell1_z = point1[2],
        cell2_x = world2cell_factor*point2[0] + world2cell_offset,
        cell2_y = world2cell_factor*point2[1] + world2cell_offset,
        cell2_z = point2[2],
        cell3_x = world2cell_factor*point3[0] + world2cell_offset,
        cell3_y = world2cell_factor*point3[1] + world2cell_offset,
        cell3_z = point3[2];
        
        int min_cell_x = (std::max) (0, int (pcl_lrint (ceil ( (std::min) (cell1_x, (std::min) (cell2_x, cell3_x)))))),
        max_cell_x = (std::min) (pixel_size-1, int (pcl_lrint (floor ( (std::max) (cell1_x,
                                                                                   (std::max) (cell2_x, cell3_x)))))),
        min_cell_y = (std::max) (0, int (pcl_lrint (ceil ( (std::min) (cell1_y, (std::min) (cell2_y, cell3_y)))))),
        max_cell_y = (std::min) (pixel_size-1, int (pcl_lrint (floor ( (std::max) (cell1_y,
                                                                                   (std::max) (cell2_y, cell3_y))))));
        if (max_cell_x<min_cell_x || max_cell_y<min_cell_y)
          continue;
        
        // We will now do the following:
        //   For each cell in the rectangle defined by the four values above,
        //   we test if it is in the original triangle (in 2D).
        //   If this is the case, we calculate the actual point on the 3D triangle, thereby interpolating the result
        //   See http://www.blackpawn.com/texts/pointinpoly/default.html
        Eigen::Vector2f cell1 (cell1_x, cell1_y),
        cell2 (cell2_x, cell2_y),
        cell3 (cell3_x, cell3_y),
        v0 = cell3 - cell1,
        v1 = cell2 - cell1;
        float dot00 = v0.dot (v0),
        dot01 = v0.dot (v1),
        dot11 = v1.dot (v1),
        invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
        
        for (int cell_x=min_cell_x; cell_x<=max_cell_x; ++cell_x)
        {
          for (int cell_y=min_cell_y; cell_y<=max_cell_y; ++cell_y)
          {
            Eigen::Vector2f current_cell (cell_x, cell_y),
            v2 = current_cell - cell1;
            float dot02 = v0.dot (v2),
            dot12 = v1.dot (v2),
            u = (dot11 * dot02 - dot01 * dot12) * invDenom,
            v = (dot00 * dot12 - dot01 * dot02) * invDenom;
            bool point_in_triangle = (u > -0.01) && (v >= -0.01) && (u + v <= 1.01);
            
            if (!point_in_triangle)
              continue;
            
            float new_value = cell1_z + u* (cell3_z-cell1_z) + v* (cell2_z-cell1_z);
            
            float& value = surface_patch[cell_y*pixel_size + cell_x];
            if (pcl_isinf (value))
              value = new_value;
            else
              value = (std::min) (value, new_value);
          }
        }
      }
    }
  }
  
  // Now find out, if there should be max ranges
  for (int cell_y=0; cell_y<pixel_size; ++cell_y)
  {
    for (int cell_x=0; cell_x<pixel_size; ++cell_x)
    {
      int index= cell_y*pixel_size + cell_x;
      float& value = surface_patch[index];
      if (!pcl_isinf (value))
        continue;
      
      // Go through immediate neighbors
      bool is_background = false;
      for (int cell2_y=cell_y-1; cell2_y<=cell_y+1&&!is_background; ++cell2_y)
      {
        for (int cell2_x=cell_x-1; cell2_x<=cell_x+1; ++cell2_x)
        {
          if (cell2_x<0||cell2_x>=pixel_size||cell2_y<0||cell2_y>=pixel_size || (cell2_x==cell_x && cell2_y==cell_y))
            continue;
          float neighbor_value = surface_patch[cell2_y*pixel_size + cell2_x];
          if (pcl_isfinite (neighbor_value))
          {
            float cell_pos_x = static_cast<float> (cell_x) + 0.6f * static_cast<float> (cell_x - cell2_x),
            cell_pos_y = static_cast<float> (cell_y) + 0.6f * static_cast<float> (cell_y - cell2_y);
            Eigen::Vector3f fake_point (cell2world_factor* (cell_pos_x)+cell2world_offset,
                                        cell2world_factor*cell_pos_y+cell2world_offset, neighbor_value);
            fake_point = inverse_pose*fake_point;
            float range_difference = getRangeDifference (fake_point);
            if (range_difference > max_dist)
            {
              value = std::numeric_limits<float>::infinity ();
              is_background = true;
              break;
            }
          }
        }
      }
      if (is_background)
      {
        // Set all -INFINITY neighbors to INFINITY
        for (int cell2_y=cell_y-1; cell2_y<=cell_y+1; ++cell2_y)
        {
          for (int cell2_x=cell_x-1; cell2_x<=cell_x+1; ++cell2_x)
          {
            if (cell2_x<0||cell2_x>=pixel_size||cell2_y<0||cell2_y>=pixel_size || (cell2_x==cell_x && cell2_y==cell_y))
              continue;
            int index2 = cell2_y*pixel_size + cell2_x;
            float& neighbor_value = surface_patch[index2];
            if (pcl_isinf (neighbor_value) && neighbor_value<0)
              neighbor_value = std::numeric_limits<float>::infinity ();
          }
        }
      }
    }
  }
  
  return surface_patch;
}

/////////////////////////////////////////////////////////////////////////
float* RangeImageLocal::getInterpolatedSurfaceProjection (const Eigen::Vector3f& point, int pixel_size,
                                                     float world_size) const
{
  Eigen::Affine3f pose = getTransformationToViewerCoordinateFrame (point);
  return (getInterpolatedSurfaceProjection (pose, pixel_size, world_size));
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getNormalBasedUprightTransformation (const Eigen::Vector3f& point, float max_dist,
                                                 Eigen::Affine3f& transformation) const
{
  int x, y;
  getImagePoint (point, x, y);
  
  Eigen::Vector3f neighbor;
  pcl::VectorAverage3f vector_average;
  float max_dist_squared=max_dist*max_dist, max_dist_reciprocal=1.0f/max_dist;
  
  bool still_in_range = true;
  for (int radius=1;  still_in_range;  ++radius)
  {
    int x2=x-radius-1, y2=y-radius;  // Top left - 1
    still_in_range = false;
    for (int i=0; i<8*radius; ++i)
    {
      if (i<=2*radius) ++x2; else if (i<=4*radius) ++y2; else if (i<=6*radius) --x2; else --y2;
      if (!isValid (x2, y2))
      {
        continue;
      }
      getPoint (x2, y2, neighbor);
      float distance_squared = (neighbor-point).squaredNorm ();
      if (distance_squared > max_dist_squared)
      {
        continue;
      }
      still_in_range = true;
      float distance = sqrtf (distance_squared),
      weight = distance*max_dist_reciprocal;
      vector_average.add (neighbor, weight);
    }
  }
  
  Eigen::Vector3f normal, point_on_plane;
  if (vector_average.getNoOfSamples () > 10)
  {
    Eigen::Vector3f eigen_values, eigen_vector2, eigen_vector3;
    vector_average.doPCA (eigen_values, normal, eigen_vector2, eigen_vector3);
    if (normal.dot ( (vector_average.getMean ()-getSensorPos ()).normalized ()) < 0.0f)
      normal *= -1.0f;
    point_on_plane = (normal.dot (vector_average.getMean ()) - normal.dot (point))*normal + point;
  }
  else
  {
    if (!getNormalForClosestNeighbors (x, y, 2, point, 15, normal, &point_on_plane, 1))
      return false;
  }
  pcl::getTransformationFromTwoUnitVectorsAndOrigin (Eigen::Vector3f (0.0f, 1.0f, 0.0f),
                                                normal, point_on_plane, transformation);
  
  return (true);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getSurfaceAngleChangeImages (int radius, float*& angle_change_image_x, float*& angle_change_image_y) const
{
  int size = width*height;
  angle_change_image_x = new float[size];
  angle_change_image_y = new float[size];
  for (int y=0; y<int (height); ++y)
  {
    for (int x=0; x<int (width); ++x)
    {
      int index = y*width+x;
      getSurfaceAngleChange (x, y, radius, angle_change_image_x[index], angle_change_image_y[index]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getAcutenessValueImages (int pixel_distance, float*& acuteness_value_image_x,
                                     float*& acuteness_value_image_y) const
{
  int size = width*height;
  acuteness_value_image_x = new float[size];
  acuteness_value_image_y = new float[size];
  for (int y=0; y<int (height); ++y)
  {
    for (int x=0; x<int (width); ++x)
    {
      int index = y*width+x;
      acuteness_value_image_x[index] = getAcutenessValue (x, y, x+pixel_distance, y);
      acuteness_value_image_y[index] = getAcutenessValue (x, y, x, y+pixel_distance);
    }
  }
}

/////////////////////////////////////////////////////////////////////////
float*
RangeImageLocal::getImpactAngleImageBasedOnLocalNormals (int radius) const
{
  int size = width*height;
  float* impact_angle_image = new float[size];
  for (int y=0; y<int (height); ++y)
  {
    for (int x=0; x<int (width); ++x)
    {
      impact_angle_image[y*width+x] = getImpactAngleBasedOnLocalNormal (x, y, radius);
    }
  }
  return impact_angle_image;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getRangeImageWithSmoothedSurface (int radius, RangeImageLocal& smoothed_range_image) const
{
  int step_size = (std::max) (1, radius/2);
  int no_of_nearest_neighbors = static_cast<int> (pow (static_cast<double> (radius / step_size + 1), 2.0));
  
  smoothed_range_image = *this;
  Eigen::Vector3f sensor_pos = getSensorPos ();
  for (int y=0; y<int (height); ++y)
  {
    for (int x=0; x<int (width); ++x)
    {
      PointWithRange& point = smoothed_range_image.getPoint (x, y);
      if (pcl_isinf (point.range))
        continue;
      Eigen::Vector3f normal, mean, eigen_values;
      float used_squared_max_distance;
      getSurfaceInformation (x, y, radius, point.getVector3fMap (), no_of_nearest_neighbors,
                             step_size, used_squared_max_distance,
                             normal, mean, eigen_values);
      
      Eigen::Vector3f viewing_direction = (point.getVector3fMap ()-sensor_pos).normalized ();
      float new_range = normal.dot (mean-sensor_pos) / normal.dot (viewing_direction);
      point.range = new_range;
      calculate3DPoint (static_cast<float> (x), static_cast<float> (y), point.range, point);
      
      const PointWithRange& original_point = getPoint (x, y);
      float distance_squared = squaredEuclideanDistance (original_point, point);
      if (distance_squared > used_squared_max_distance)
        point = original_point;
    }
  }
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::extractFarRanges (const sensor_msgs::PointCloud2& point_cloud_data,
                              PointCloud<PointWithViewpoint>& far_ranges)
{
  int x_idx = -1, y_idx = -1, z_idx = -1,
  vp_x_idx = -1, vp_y_idx = -1, vp_z_idx = -1, distance_idx = -1;
  for (int d = 0; d < static_cast<int> (point_cloud_data.fields.size ()); ++d)
  {
    if (point_cloud_data.fields[d].name == "x") x_idx = d;
    if (point_cloud_data.fields[d].name == "y") y_idx = d;
    if (point_cloud_data.fields[d].name == "z") z_idx = d;
    if (point_cloud_data.fields[d].name == "vp_x") vp_x_idx = d;
    if (point_cloud_data.fields[d].name == "vp_y") vp_y_idx = d;
    if (point_cloud_data.fields[d].name == "vp_z") vp_z_idx = d;
    if (point_cloud_data.fields[d].name == "distance") distance_idx = d;
  }
  
  if (x_idx<0 || y_idx<0 || z_idx<0 || vp_x_idx<0 || vp_y_idx<0 || vp_z_idx<0 || distance_idx<0)
  {
    return;
  }
  
  int point_step = point_cloud_data.point_step;
  const unsigned char* data = &point_cloud_data.data[0];
  int x_offset = point_cloud_data.fields[x_idx].offset,
  y_offset = point_cloud_data.fields[y_idx].offset,
  z_offset = point_cloud_data.fields[z_idx].offset,
  vp_x_offset = point_cloud_data.fields[vp_x_idx].offset,
  vp_y_offset = point_cloud_data.fields[vp_y_idx].offset,
  vp_z_offset = point_cloud_data.fields[vp_z_idx].offset,
  distance_offset = point_cloud_data.fields[distance_idx].offset;
  
  for (size_t point_idx = 0; point_idx < point_cloud_data.width*point_cloud_data.height; ++point_idx)
  {
    float x = *reinterpret_cast<const float*> (data+x_offset),
    y = *reinterpret_cast<const float*> (data+y_offset),
    z = *reinterpret_cast<const float*> (data+z_offset),
    vp_x = *reinterpret_cast<const float*> (data+vp_x_offset),
    vp_y = *reinterpret_cast<const float*> (data+vp_y_offset),
    vp_z = *reinterpret_cast<const float*> (data+vp_z_offset),
    distance = *reinterpret_cast<const float*> (data+distance_offset);
    data+=point_step;
    
    if (!pcl_isfinite (x) && pcl_isfinite (distance))
    {
      PointWithViewpoint point;
      point.x=distance; point.y=y; point.z=z;
      point.vp_x=vp_x; point.vp_y=vp_y; point.vp_z=vp_z;
      far_ranges.points.push_back (point);
    }
  }
  far_ranges.width= static_cast<uint32_t> (far_ranges.points.size ());  far_ranges.height = 1;
  far_ranges.is_dense = false;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getOverlap (const RangeImageLocal& other_range_image, const Eigen::Affine3f& relative_transformation,
                        int search_radius, float max_distance, int pixel_step) const
{
  int hits_counter=0, valid_points_counter=0;
  
  float max_distance_squared = max_distance*max_distance;
  
# pragma omp parallel for num_threads (max_no_of_threads) default (shared) schedule (dynamic, 1) \
reduction (+ : valid_points_counter) reduction (+ : hits_counter)
  for (int other_y=0; other_y<int (other_range_image.height); other_y+=pixel_step)
  {
    for (int other_x=0; other_x<int (other_range_image.width); other_x+=pixel_step)
    {
      const PointWithRange& point = other_range_image.getPoint (other_x, other_y);
      if (!pcl_isfinite (point.range))
        continue;
      ++valid_points_counter;
      Eigen::Vector3f transformed_point = relative_transformation * point.getVector3fMap ();
      int x,y;
      getImagePoint (transformed_point, x, y);
      float closest_distance = max_distance_squared;
      Eigen::Vector3f closest_point (0.0f, 0.0f, 0.0f);
      bool found_neighbor = false;
      for (int y2=y-pixel_step*search_radius; y2<=y+pixel_step*search_radius; y2+=pixel_step)
      {
        for (int x2=x-pixel_step*search_radius; x2<=x+pixel_step*search_radius; x2+=pixel_step)
        {
          const PointWithRange& neighbor = getPoint (x2, y2);
          if (!pcl_isfinite (neighbor.range))
            continue;
          float distance = (transformed_point-neighbor.getVector3fMap ()).squaredNorm ();
          if (distance < closest_distance)
          {
            closest_distance = distance;
            closest_point = neighbor.getVector3fMap ();
            found_neighbor = true;
          }
        }
      }
      
      if (found_neighbor)
      {
        ++hits_counter;
      }
    }
  }
  return static_cast<float> (hits_counter)/static_cast<float> (valid_points_counter);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getBlurredImageUsingIntegralImage (int blur_radius, float* integral_image, int* valid_points_num_image,
                                               RangeImageLocal& blurred_image) const
{
  this->copyTo(blurred_image);
  for (int y=0; y<int (height); ++y)
  {
    for (int x=0; x<int (width); ++x)
    {
      const PointWithRange& old_point = getPoint (x, y);
      PointWithRange& new_point = blurred_image.getPoint (x, y);
      if (!pcl_isfinite (old_point.range))
        continue;
      
      int top= (std::max) (-1, y-blur_radius-1), right = (std::min) (static_cast<int> (width)-1, x+blur_radius), bottom =
      (std::min) (static_cast<int> (height)-1, y+blur_radius), left= (std::max) (-1, x-blur_radius-1);
      
      float top_left_value=0, top_right_value=0,
      bottom_right_value=integral_image[bottom*width+right], bottom_left_value=0;
      int top_left_valid_points=0, top_right_valid_points=0,
      bottom_right_valid_points=valid_points_num_image[bottom*width+right], bottom_left_valid_points=0;
      if (left>=0)
      {
        bottom_left_value = integral_image[bottom*width+left];
        bottom_left_valid_points = valid_points_num_image[bottom*width+left];
        if (top>=0)
        {
          top_left_value = integral_image[top*width+left];
          top_left_valid_points = valid_points_num_image[top*width+left];
        }
      }
      if (top>=0)
      {
        top_right_value = integral_image[top*width+right];
        top_right_valid_points = valid_points_num_image[top*width+right];
      }
      int valid_points_num = bottom_right_valid_points + top_left_valid_points - bottom_left_valid_points -
      top_right_valid_points;
      new_point.range = (bottom_right_value + top_left_value - bottom_left_value - top_right_value) /
      static_cast<float> (valid_points_num);
    }
  }
  blurred_image.recalculate3DPointPositions ();
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getBlurredImage (int blur_radius, RangeImageLocal& blurred_image) const
{
  //MEASURE_FUNCTION_TIME;
  
  if (blur_radius > 1)  // For a high blur radius it's faster to use integral images
  {
    float* integral_image;
    int* valid_points_num_image;
    getIntegralImage (integral_image, valid_points_num_image);
    getBlurredImageUsingIntegralImage (blur_radius, integral_image, valid_points_num_image, blurred_image);
    delete[] integral_image;
    delete[] valid_points_num_image;
    return;
  }
  
  this->copyTo(blurred_image);
  
  if (blur_radius==0)
    return;
  
  for (int y=0; y < static_cast<int> (height); ++y)
  {
    for (int x=0; x < static_cast<int> (width); ++x)
    {
      PointWithRange& new_point = blurred_image.getPoint (x, y);
      const PointWithRange& original_point = getPoint (x, y);
      if (!pcl_isfinite (original_point.range))
        continue;
      
      new_point.range = 0.0f;
      float weight_sum = 0.0f;
      for (int y2=y-blur_radius; y2<y+blur_radius; ++y2)
      {
        for (int x2=x-blur_radius; x2<x+blur_radius; ++x2)
        {
          if (!isValid (x2,y2))
            continue;
          new_point.range += getPoint (x2, y2).range;
          weight_sum += 1.0f;
        }
      }
      new_point.range /= weight_sum;
    }
  }
  blurred_image.recalculate3DPointPositions ();
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::copyTo (RangeImageLocal& other) const
{
  other = *this;
}



/////////////////////////////////////////////////////////////////////////
inline float
RangeImageLocal::asinLookUp (float value)
{
  return (asin_lookup_table[
                            static_cast<int> (
                                              static_cast<float> (pcl_lrintf ( (static_cast<float> (lookup_table_size) / 2.0f) * value)) +
                                              static_cast<float> (lookup_table_size) / 2.0f)]);
}

/////////////////////////////////////////////////////////////////////////
inline float
RangeImageLocal::atan2LookUp (float y, float x)
{
  float ret;
  if (fabsf (x) < fabsf (y))
  {
    ret = atan_lookup_table[
                            static_cast<int> (
                                              static_cast<float> (pcl_lrintf ( (static_cast<float> (lookup_table_size) / 2.0f) * (x / y))) +
                                              static_cast<float> (lookup_table_size) / 2.0f)];
    ret = static_cast<float> (x*y > 0 ? M_PI/2-ret : -M_PI/2-ret);
  }
  else
    ret = atan_lookup_table[
                            static_cast<int> (
                                              static_cast<float> (pcl_lrintf ( (static_cast<float> (lookup_table_size) / 2.0f) * (y / x))) +
                                              static_cast<float> (lookup_table_size)/2.0f)];
  if (x < 0)
    ret = static_cast<float> (y < 0 ? ret-M_PI : ret+M_PI);
  
  return (ret);
}

/////////////////////////////////////////////////////////////////////////
inline float
RangeImageLocal::cosLookUp (float value)
{
  int cell_idx = static_cast<int> (pcl_lrintf ( (static_cast<float> (lookup_table_size)-1) * fabsf (value) / (2.0f * M_PI)));
  return (cos_lookup_table[cell_idx]);
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudType> void
RangeImageLocal::createFromPointCloud (const PointCloudType& point_cloud, float angular_resolution,
                                       float max_angle_width, float max_angle_height,
                                       const Eigen::Affine3f& sensor_pose, RangeImageLocal::CoordinateFrame coordinate_frame,
                                       float noise_level, float min_range, int border_size)
{
  createFromPointCloud (point_cloud, angular_resolution, angular_resolution, max_angle_width, max_angle_height,
                        sensor_pose, coordinate_frame, noise_level, min_range, border_size);
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudType> void
RangeImageLocal::createFromPointCloud (const PointCloudType& point_cloud,
                                       float angular_resolution_x, float angular_resolution_y,
                                       float max_angle_width, float max_angle_height,
                                       const Eigen::Affine3f& sensor_pose, RangeImageLocal::CoordinateFrame coordinate_frame,
                                       float noise_level, float min_range, int border_size)
{
  setAngularResolution (angular_resolution_x, angular_resolution_y);
  
  width  = static_cast<uint32_t> (pcl_lrint (floor (max_angle_width*angular_resolution_x_reciprocal_)));
  height = static_cast<uint32_t> (pcl_lrint (floor (max_angle_height*angular_resolution_y_reciprocal_)));
  
  int full_width  = static_cast<int> (pcl_lrint (floor (deg2radLocal (360.0f)*angular_resolution_x_reciprocal_))),
  full_height = static_cast<int> (pcl_lrint (floor (deg2radLocal (180.0f)*angular_resolution_y_reciprocal_)));
  image_offset_x_ = (full_width -static_cast<int> (width) )/2;
  image_offset_y_ = (full_height-static_cast<int> (height))/2;
  is_dense = false;
  
  getCoordinateFrameTransformation (coordinate_frame, to_world_system_);
  to_world_system_ = sensor_pose * to_world_system_;
  
  to_range_image_system_ = to_world_system_.inverse (Eigen::Isometry);
  //std::cout << "to_world_system_ is\n"<<to_world_system_<<"\nand to_range_image_system_ is\n"<<to_range_image_system_<<"\n\n";
  
  unsigned int size = width*height;
  points.clear ();
  points.resize (size, unobserved_point);
  
  int top=height, right=-1, bottom=-1, left=width;
  doZBuffer (point_cloud, noise_level, min_range, top, right, bottom, left);
  
  cropImage (border_size, top, right, bottom, left);
  
  recalculate3DPointPositions ();
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudType> void
RangeImageLocal::createFromPointCloudWithKnownSize (const PointCloudType& point_cloud, float angular_resolution,
                                                    const Eigen::Vector3f& point_cloud_center, float point_cloud_radius,
                                                    const Eigen::Affine3f& sensor_pose, RangeImageLocal::CoordinateFrame coordinate_frame,
                                                    float noise_level, float min_range, int border_size)
{
  createFromPointCloudWithKnownSize (point_cloud, angular_resolution, angular_resolution, point_cloud_center, point_cloud_radius,
                                     sensor_pose, coordinate_frame, noise_level, min_range, border_size);
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudType> void
RangeImageLocal::createFromPointCloudWithKnownSize (const PointCloudType& point_cloud,
                                                    float angular_resolution_x, float angular_resolution_y,
                                                    const Eigen::Vector3f& point_cloud_center, float point_cloud_radius,
                                                    const Eigen::Affine3f& sensor_pose, RangeImageLocal::CoordinateFrame coordinate_frame,
                                                    float noise_level, float min_range, int border_size)
{
  //MEASURE_FUNCTION_TIME;
  
  //std::cout << "Starting to create range image from "<<point_cloud.points.size ()<<" points.\n";
  
  // If the sensor pose is inside of the sphere we have to calculate the image the normal way
  if ((point_cloud_center-sensor_pose.translation()).norm() <= point_cloud_radius) {
    createFromPointCloud (point_cloud, angular_resolution_x, angular_resolution_y,
                          deg2radLocal (360.0f), deg2radLocal (180.0f),
                          sensor_pose, coordinate_frame, noise_level, min_range, border_size);
    return;
  }
  
  setAngularResolution (angular_resolution_x, angular_resolution_y);
  
  getCoordinateFrameTransformation (coordinate_frame, to_world_system_);
  to_world_system_ = sensor_pose * to_world_system_;
  to_range_image_system_ = to_world_system_.inverse (Eigen::Isometry);
  
  float max_angle_size = getMaxAngleSize (sensor_pose, point_cloud_center, point_cloud_radius);
  int pixel_radius_x = pcl_lrint (ceil (0.5f*max_angle_size*angular_resolution_x_reciprocal_)),
  pixel_radius_y = pcl_lrint (ceil (0.5f*max_angle_size*angular_resolution_y_reciprocal_));
  width  = 2*pixel_radius_x;
  height = 2*pixel_radius_y;
  is_dense = false;
  
  image_offset_x_ = image_offset_y_ = 0;  // temporary values for getImagePoint
  int center_pixel_x, center_pixel_y;
  getImagePoint (point_cloud_center, center_pixel_x, center_pixel_y);
  image_offset_x_ = (std::max) (0, center_pixel_x-pixel_radius_x);
  image_offset_y_ = (std::max) (0, center_pixel_y-pixel_radius_y);
  
  points.clear ();
  points.resize (width*height, unobserved_point);
  
  int top=height, right=-1, bottom=-1, left=width;
  doZBuffer (point_cloud, noise_level, min_range, top, right, bottom, left);
  
  cropImage (border_size, top, right, bottom, left);
  
  recalculate3DPointPositions ();
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudTypeWithViewpoints> void
RangeImageLocal::createFromPointCloudWithViewpoints (const PointCloudTypeWithViewpoints& point_cloud,
                                                     float angular_resolution,
                                                     float max_angle_width, float max_angle_height,
                                                     RangeImageLocal::CoordinateFrame coordinate_frame,
                                                     float noise_level, float min_range, int border_size)
{
  createFromPointCloudWithViewpoints (point_cloud, angular_resolution, angular_resolution,
                                      max_angle_width, max_angle_height, coordinate_frame,
                                      noise_level, min_range, border_size);
}

/////////////////////////////////////////////////////////////////////////
template <typename PointCloudTypeWithViewpoints> void
RangeImageLocal::createFromPointCloudWithViewpoints (const PointCloudTypeWithViewpoints& point_cloud,
                                                     float angular_resolution_x, float angular_resolution_y,
                                                     float max_angle_width, float max_angle_height,
                                                     RangeImageLocal::CoordinateFrame coordinate_frame,
                                                     float noise_level, float min_range, int border_size)
{
  Eigen::Vector3f average_viewpoint = getAverageViewPoint (point_cloud);
  Eigen::Affine3f sensor_pose = static_cast<Eigen::Affine3f> (Eigen::Translation3f (average_viewpoint));
  createFromPointCloud (point_cloud, angular_resolution_x, angular_resolution_y, max_angle_width, max_angle_height,
                        sensor_pose, coordinate_frame, noise_level, min_range, border_size);
}

/////////////////////////////////////////////////////////////////////////

template <typename PointCloudType> void
RangeImageLocal::doZBuffer (const PointCloudType& point_cloud, float noise_level, float min_range, int& top, int& right, int& bottom, int& left)
{
  typedef typename PointCloudType::PointType PointType2;
  const typename pcl::PointCloud<PointType2>::VectorType &points2 = point_cloud.points;
  
  unsigned int size = width*height;
  int* counters = new int[size];
  ERASE_ARRAY (counters, size);
  
  top=height; right=-1; bottom=-1; left=width;
  
  float x_real, y_real, range_of_current_point;
  int x, y;
  for (typename pcl::PointCloud<PointType2>::VectorType::const_iterator it=points2.begin (); it!=points2.end (); ++it)
  {
    //if (!isFinite (*it->get))  // Check for NAN etc
    //  continue;
    Eigen::Map<const Eigen::Vector3f > current_point = it->getVector3fMap ();
    
    this->getImagePoint (current_point, x_real, y_real, range_of_current_point);
    this->real2DToInt2D (x_real, y_real, x, y);
    
    if (range_of_current_point < min_range|| !isInImage (x, y))
      continue;
    //std::cout << " ("<<current_point[0]<<", "<<current_point[1]<<", "<<current_point[2]<<") falls into pixel "<<x<<","<<y<<".\n";
    
    // Do some minor interpolation by checking the three closest neighbors to the point, that are not filled yet.
    int floor_x = pcl_lrint (floor (x_real)), floor_y = pcl_lrint (floor (y_real)),
    ceil_x  = pcl_lrint (ceil (x_real)),  ceil_y  = pcl_lrint (ceil (y_real));
    
    int neighbor_x[4], neighbor_y[4];
    neighbor_x[0]=floor_x; neighbor_y[0]=floor_y;
    neighbor_x[1]=floor_x; neighbor_y[1]=ceil_y;
    neighbor_x[2]=ceil_x;  neighbor_y[2]=floor_y;
    neighbor_x[3]=ceil_x;  neighbor_y[3]=ceil_y;
    //std::cout << x_real<<","<<y_real<<": ";
    
    for (int i=0; i<4; ++i)
    {
      int n_x=neighbor_x[i], n_y=neighbor_y[i];
      //std::cout << n_x<<","<<n_y<<" ";
      if (n_x==x && n_y==y)
        continue;
      if (isInImage (n_x, n_y))
      {
        int neighbor_array_pos = n_y*width + n_x;
        if (counters[neighbor_array_pos]==0)
        {
          float& neighbor_range = points[neighbor_array_pos].range;
          neighbor_range = (pcl_isinf (neighbor_range) ? range_of_current_point : (std::min) (neighbor_range, range_of_current_point));
          top= (std::min) (top, n_y); right= (std::max) (right, n_x); bottom= (std::max) (bottom, n_y); left= (std::min) (left, n_x);
        }
      }
    }
    //std::cout <<std::endl;
    
    // The point itself
    int arrayPos = y*width + x;
    float& range_at_image_point = points[arrayPos].range;
    int& counter = counters[arrayPos];
    bool addCurrentPoint=false, replace_with_current_point=false;
    
    if (counter==0)
    {
      replace_with_current_point = true;
    }
    else
    {
      if (range_of_current_point < range_at_image_point-noise_level)
      {
        replace_with_current_point = true;
      }
      else if (fabs (range_of_current_point-range_at_image_point)<=noise_level)
      {
        addCurrentPoint = true;
      }
    }
    
    if (replace_with_current_point)
    {
      counter = 1;
      range_at_image_point = range_of_current_point;
      top= (std::min) (top, y); right= (std::max) (right, x); bottom= (std::max) (bottom, y); left= (std::min) (left, x);
      //std::cout << "Adding point "<<x<<","<<y<<"\n";
    }
    else if (addCurrentPoint)
    {
      ++counter;
      range_at_image_point += (range_of_current_point-range_at_image_point)/counter;
    }
  }
  
  delete[] counters;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (float x, float y, float z, float& image_x, float& image_y, float& range) const
{
  Eigen::Vector3f point (x, y, z);
  getImagePoint (point, image_x, image_y, range);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (float x, float y, float z, float& image_x, float& image_y) const
{
  float range;
  getImagePoint (x, y, z, image_x, image_y, range);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (float x, float y, float z, int& image_x, int& image_y) const
{
  float image_x_float, image_y_float;
  getImagePoint (x, y, z, image_x_float, image_y_float);
  real2DToInt2D (image_x_float, image_y_float, image_x, image_y);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (const Eigen::Vector3f& point, float& image_x, float& image_y, float& range) const
{
  Eigen::Vector3f transformedPoint = to_range_image_system_ * point;
  range = transformedPoint.norm ();
  float angle_x = atan2LookUp (transformedPoint[0], transformedPoint[2]),
  angle_y = asinLookUp (transformedPoint[1]/range);
  getImagePointFromAngles (angle_x, angle_y, image_x, image_y);
  //std::cout << " ("<<point[0]<<","<<point[1]<<","<<point[2]<<")"
  //<< " => ("<<transformedPoint[0]<<","<<transformedPoint[1]<<","<<transformedPoint[2]<<")"
  //<< " => "<<angle_x<<","<<angle_y<<" => "<<image_x<<","<<image_y<<"\n";
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (const Eigen::Vector3f& point, int& image_x, int& image_y, float& range) const {
  float image_x_float, image_y_float;
  getImagePoint (point, image_x_float, image_y_float, range);
  real2DToInt2D (image_x_float, image_y_float, image_x, image_y);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (const Eigen::Vector3f& point, float& image_x, float& image_y) const
{
  float range;
  getImagePoint (point, image_x, image_y, range);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePoint (const Eigen::Vector3f& point, int& image_x, int& image_y) const
{
  float image_x_float, image_y_float;
  getImagePoint (point, image_x_float, image_y_float);
  real2DToInt2D (image_x_float, image_y_float, image_x, image_y);
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::checkPoint (const Eigen::Vector3f& point, PointWithRange& point_in_image) const
{
  int image_x, image_y;
  float range;
  getImagePoint (point, image_x, image_y, range);
  if (!isInImage (image_x, image_y))
    point_in_image = unobserved_point;
  else
    point_in_image = getPoint (image_x, image_y);
  return range;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getRangeDifference (const Eigen::Vector3f& point) const
{
  int image_x, image_y;
  float range;
  getImagePoint (point, image_x, image_y, range);
  if (!isInImage (image_x, image_y))
    return -std::numeric_limits<float>::infinity ();
  float image_point_range = getPoint (image_x, image_y).range;
  if (pcl_isinf (image_point_range))
  {
    if (image_point_range > 0.0f)
      return std::numeric_limits<float>::infinity ();
    else
      return -std::numeric_limits<float>::infinity ();
  }
  return image_point_range - range;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getImagePointFromAngles (float angle_x, float angle_y, float& image_x, float& image_y) const
{
  image_x = (angle_x*cosLookUp (angle_y) + static_cast<float> (M_PI))*angular_resolution_x_reciprocal_ - static_cast<float> (image_offset_x_);
  image_y = (angle_y + 0.5f*static_cast<float> (M_PI))*angular_resolution_y_reciprocal_ - static_cast<float> (image_offset_y_);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::real2DToInt2D (float x, float y, int& xInt, int& yInt) const
{
  xInt = static_cast<int> (pcl_lrintf (x));
  yInt = static_cast<int> (pcl_lrintf (y));
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::isInImage (int x, int y) const
{
  return (x >= 0 && x < static_cast<int> (width) && y >= 0 && y < static_cast<int> (height));
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::isValid (int x, int y) const
{
  return isInImage (x,y) && pcl_isfinite (getPoint (x,y).range);
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::isValid (int index) const
{
  return pcl_isfinite (getPoint (index).range);
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::isObserved (int x, int y) const
{
  if (!isInImage (x,y) || (pcl_isinf (getPoint (x,y).range)&&getPoint (x,y).range<0.0f))
    return false;
  return true;
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::isMaxRange (int x, int y) const
{
  float range = getPoint (x,y).range;
  return pcl_isinf (range) && range>0.0f;
}

/////////////////////////////////////////////////////////////////////////
const PointWithRange&
RangeImageLocal::getPoint (int image_x, int image_y) const
{
  if (!isInImage (image_x, image_y))
    return unobserved_point;
  return points[image_y*width + image_x];
}

/////////////////////////////////////////////////////////////////////////
const PointWithRange&
RangeImageLocal::getPointNoCheck (int image_x, int image_y) const
{
  return points[image_y*width + image_x];
}

/////////////////////////////////////////////////////////////////////////
PointWithRange&
RangeImageLocal::getPointNoCheck (int image_x, int image_y)
{
  return points[image_y*width + image_x];
}

/////////////////////////////////////////////////////////////////////////
PointWithRange&
RangeImageLocal::getPoint (int image_x, int image_y)
{
  return points[image_y*width + image_x];
}


/////////////////////////////////////////////////////////////////////////
const PointWithRange&
RangeImageLocal::getPoint (int index) const
{
  return points[index];
}

/////////////////////////////////////////////////////////////////////////
const PointWithRange&
RangeImageLocal::getPoint (float image_x, float image_y) const
{
  int x, y;
  real2DToInt2D (image_x, image_y, x, y);
  return getPoint (x, y);
}

/////////////////////////////////////////////////////////////////////////
PointWithRange&
RangeImageLocal::getPoint (float image_x, float image_y)
{
  int x, y;
  real2DToInt2D (image_x, image_y, x, y);
  return getPoint (x, y);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getPoint (int image_x, int image_y, Eigen::Vector3f& point) const
{
  //std::cout << getPoint (image_x, image_y)<< " - " << getPoint (image_x, image_y).getVector3fMap ()<<"\n";
  point = getPoint (image_x, image_y).getVector3fMap ();
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getPoint (int index, Eigen::Vector3f& point) const
{
  point = getPoint (index).getVector3fMap ();
}

/////////////////////////////////////////////////////////////////////////
const Eigen::Map<const Eigen::Vector3f>
RangeImageLocal::getEigenVector3f (int x, int y) const
{
  return getPoint (x, y).getVector3fMap ();
}

/////////////////////////////////////////////////////////////////////////
const Eigen::Map<const Eigen::Vector3f>
RangeImageLocal::getEigenVector3f (int index) const
{
  return getPoint (index).getVector3fMap ();
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::calculate3DPoint (float image_x, float image_y, float range, Eigen::Vector3f& point) const
{
  float angle_x, angle_y;
  //std::cout << image_x<<","<<image_y<<","<<range;
  getAnglesFromImagePoint (image_x, image_y, angle_x, angle_y);
  
  float cosY = cosf (angle_y);
  point = Eigen::Vector3f (range * sinf (angle_x) * cosY, range * sinf (angle_y), range * cosf (angle_x)*cosY);
  point = to_world_system_ * point;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::calculate3DPoint (float image_x, float image_y, Eigen::Vector3f& point) const
{
  const PointWithRange& point_in_image = getPoint (image_x, image_y);
  calculate3DPoint (image_x, image_y, point_in_image.range, point);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::calculate3DPoint (float image_x, float image_y, float range, PointWithRange& point) const {
  point.range = range;
  Eigen::Vector3f tmp_point;
  calculate3DPoint (image_x, image_y, range, tmp_point);
  point.x=tmp_point[0];  point.y=tmp_point[1];  point.z=tmp_point[2];
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::calculate3DPoint (float image_x, float image_y, PointWithRange& point) const
{
  const PointWithRange& point_in_image = getPoint (image_x, image_y);
  calculate3DPoint (image_x, image_y, point_in_image.range, point);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getAnglesFromImagePoint (float image_x, float image_y, float& angle_x, float& angle_y) const
{
  angle_y = (image_y+static_cast<float> (image_offset_y_))*angular_resolution_y_ - 0.5f*static_cast<float> (M_PI);
  float cos_angle_y = cosf (angle_y);
  angle_x = (cos_angle_y==0.0f ? 0.0f : ( (image_x+ static_cast<float> (image_offset_x_))*angular_resolution_x_ - static_cast<float> (M_PI))/cos_angle_y);
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getImpactAngle (int x1, int y1, int x2, int y2) const
{
  if (!isInImage (x1, y1) || !isInImage (x2,y2))
    return -std::numeric_limits<float>::infinity ();
  return getImpactAngle (getPoint (x1,y1),getPoint (x2,y2));
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getImpactAngle (const PointWithRange& point1, const PointWithRange& point2) const {
  if ( (pcl_isinf (point1.range)&&point1.range<0) || (pcl_isinf (point2.range)&&point2.range<0))
    return -std::numeric_limits<float>::infinity ();
  
  float r1 = (std::min) (point1.range, point2.range),
  r2 = (std::max) (point1.range, point2.range);
  float impact_angle = static_cast<float> (0.5f * M_PI);
  
  if (pcl_isinf (r2))
  {
    if (r2 > 0.0f && !pcl_isinf (r1))
      impact_angle = 0.0f;
  }
  else if (!pcl_isinf (r1))
  {
    float r1Sqr = r1*r1,
    r2Sqr = r2*r2,
    dSqr  = squaredEuclideanDistance (point1, point2),
    d     = sqrtf (dSqr);
    float cos_impact_angle = (r2Sqr + dSqr - r1Sqr)/ (2.0f*r2*d);
    cos_impact_angle = (std::max) (0.0f, (std::min) (1.0f, cos_impact_angle));
    impact_angle = acosf (cos_impact_angle);  // Using the cosine rule
  }
  
  if (point1.range > point2.range)
    impact_angle = -impact_angle;
  
  return impact_angle;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getAcutenessValue (const PointWithRange& point1, const PointWithRange& point2) const
{
  float impact_angle = getImpactAngle (point1, point2);
  if (pcl_isinf (impact_angle))
    return -std::numeric_limits<float>::infinity ();
  float ret = 1.0f - float (fabs (impact_angle)/ (0.5f*M_PI));
  if (impact_angle < 0.0f)
    ret = -ret;
  //if (fabs (ret)>1)
  //std::cout << PVARAC (impact_angle)<<PVARN (ret);
  return ret;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getAcutenessValue (int x1, int y1, int x2, int y2) const
{
  if (!isInImage (x1, y1) || !isInImage (x2,y2))
    return -std::numeric_limits<float>::infinity ();
  return getAcutenessValue (getPoint (x1,y1), getPoint (x2,y2));
}

/////////////////////////////////////////////////////////////////////////
const Eigen::Vector3f
RangeImageLocal::getSensorPos () const
{
  return Eigen::Vector3f (to_world_system_ (0,3), to_world_system_ (1,3), to_world_system_ (2,3));
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getSurfaceAngleChange (int x, int y, int radius, float& angle_change_x, float& angle_change_y) const
{
  angle_change_x = angle_change_y = -std::numeric_limits<float>::infinity ();
  if (!isValid (x,y))
    return;
  Eigen::Vector3f point;
  getPoint (x, y, point);
  Eigen::Affine3f transformation = getTransformationToViewerCoordinateFrame (point);
  
  if (isObserved (x-radius, y) && isObserved (x+radius, y))
  {
    Eigen::Vector3f transformed_left;
    if (isMaxRange (x-radius, y))
      transformed_left = Eigen::Vector3f (0.0f, 0.0f, -1.0f);
    else
    {
      Eigen::Vector3f left;
      getPoint (x-radius, y, left);
      transformed_left = - (transformation * left);
      //std::cout << PVARN (transformed_left[1]);
      transformed_left[1] = 0.0f;
      transformed_left.normalize ();
    }
    
    Eigen::Vector3f transformed_right;
    if (isMaxRange (x+radius, y))
      transformed_right = Eigen::Vector3f (0.0f, 0.0f, 1.0f);
    else
    {
      Eigen::Vector3f right;
      getPoint (x+radius, y, right);
      transformed_right = transformation * right;
      //std::cout << PVARN (transformed_right[1]);
      transformed_right[1] = 0.0f;
      transformed_right.normalize ();
    }
    angle_change_x = transformed_left.dot (transformed_right);
    angle_change_x = (std::max) (0.0f, (std::min) (1.0f, angle_change_x));
    angle_change_x = acosf (angle_change_x);
  }
  
  if (isObserved (x, y-radius) && isObserved (x, y+radius))
  {
    Eigen::Vector3f transformed_top;
    if (isMaxRange (x, y-radius))
      transformed_top = Eigen::Vector3f (0.0f, 0.0f, -1.0f);
    else
    {
      Eigen::Vector3f top;
      getPoint (x, y-radius, top);
      transformed_top = - (transformation * top);
      //std::cout << PVARN (transformed_top[0]);
      transformed_top[0] = 0.0f;
      transformed_top.normalize ();
    }
    
    Eigen::Vector3f transformed_bottom;
    if (isMaxRange (x, y+radius))
      transformed_bottom = Eigen::Vector3f (0.0f, 0.0f, 1.0f);
    else
    {
      Eigen::Vector3f bottom;
      getPoint (x, y+radius, bottom);
      transformed_bottom = transformation * bottom;
      //std::cout << PVARN (transformed_bottom[0]);
      transformed_bottom[0] = 0.0f;
      transformed_bottom.normalize ();
    }
    angle_change_y = transformed_top.dot (transformed_bottom);
    angle_change_y = (std::max) (0.0f, (std::min) (1.0f, angle_change_y));
    angle_change_y = acosf (angle_change_y);
  }
}


//inline float RangeImageLocal::getSurfaceChange (const PointWithRange& point, const PointWithRange& neighbor1, const PointWithRange& neighbor2) const
//{
//if (!pcl_isfinite (point.range) || (!pcl_isfinite (neighbor1.range)&&neighbor1.range<0) || (!pcl_isfinite (neighbor2.range)&&neighbor2.range<0))
//return -std::numeric_limits<float>::infinity ();
//if (pcl_isinf (neighbor1.range))
//{
//if (pcl_isinf (neighbor2.range))
//return 0.0f;
//else
//return acosf ( (Eigen::Vector3f (point.x, point.y, point.z)-getSensorPos ()).normalized ().dot ( (Eigen::Vector3f (neighbor2.x, neighbor2.y, neighbor2.z)-Eigen::Vector3f (point.x, point.y, point.z)).normalized ()));
//}
//if (pcl_isinf (neighbor2.range))
//return acosf ( (Eigen::Vector3f (point.x, point.y, point.z)-getSensorPos ()).normalized ().dot ( (Eigen::Vector3f (neighbor1.x, neighbor1.y, neighbor1.z)-Eigen::Vector3f (point.x, point.y, point.z)).normalized ()));

//float d1_squared = squaredEuclideanDistance (point, neighbor1),
//d1 = sqrtf (d1_squared),
//d2_squared = squaredEuclideanDistance (point, neighbor2),
//d2 = sqrtf (d2_squared),
//d3_squared = squaredEuclideanDistance (neighbor1, neighbor2);
//float cos_surface_change = (d1_squared + d2_squared - d3_squared)/ (2.0f*d1*d2),
//surface_change = acosf (cos_surface_change);
//if (pcl_isnan (surface_change))
//surface_change = static_cast<float> (M_PI);
////std::cout << PVARN (point)<<PVARN (neighbor1)<<PVARN (neighbor2)<<PVARN (cos_surface_change)<<PVARN (surface_change)<<PVARN (d1)<<PVARN (d2)<<PVARN (d1_squared)<<PVARN (d2_squared)<<PVARN (d3_squared);

//return surface_change;
//}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getMaxAngleSize (const Eigen::Affine3f& viewer_pose, const Eigen::Vector3f& center, float radius)
{
  return 2.0f * asinf (radius/ (viewer_pose.translation ()-center).norm ());
}

/////////////////////////////////////////////////////////////////////////
Eigen::Vector3f
RangeImageLocal::getEigenVector3f (const PointWithRange& point)
{
  return Eigen::Vector3f (point.x, point.y, point.z);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::get1dPointAverage (int x, int y, int delta_x, int delta_y, int no_of_points, PointWithRange& average_point) const
{
  //std::cout << __PRETTY_FUNCTION__<<" called.\n";
  //MEASURE_FUNCTION_TIME;
  float weight_sum = 1.0f;
  average_point = getPoint (x,y);
  if (pcl_isinf (average_point.range))
  {
    if (average_point.range>0.0f)  // The first point is max range -> return a max range point
      return;
    weight_sum = 0.0f;
    average_point.x = average_point.y = average_point.z = average_point.range = 0.0f;
  }
  
  int x2=x, y2=y;
  Eigen::Map< Eigen::Vector4f, Eigen::Aligned > average_point_eigen = average_point.getVector4fMap ();
  //std::cout << PVARN (no_of_points);
  for (int step=1; step<no_of_points; ++step)
  {
    //std::cout << PVARC (step);
    x2+=delta_x;  y2+=delta_y;
    if (!isValid (x2, y2))
      continue;
    const PointWithRange& p = getPointNoCheck (x2, y2);
    average_point_eigen+=p.getVector4fMap (); average_point.range+=p.range;
    weight_sum += 1.0f;
  }
  if (weight_sum<= 0.0f)
  {
    average_point = unobserved_point;
    return;
  }
  float normalization_factor = 1.0f/weight_sum;
  average_point_eigen *= normalization_factor;
  average_point.range *= normalization_factor;
  //std::cout << PVARN (average_point);
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getEuclideanDistanceSquared (int x1, int y1, int x2, int y2) const
{
  if (!isObserved (x1,y1)||!isObserved (x2,y2))
    return -std::numeric_limits<float>::infinity ();
  const PointWithRange& point1 = getPoint (x1,y1),
  & point2 = getPoint (x2,y2);
  if (pcl_isinf (point1.range) && pcl_isinf (point2.range))
    return 0.0f;
  if (pcl_isinf (point1.range) || pcl_isinf (point2.range))
    return std::numeric_limits<float>::infinity ();
  return squaredEuclideanDistance (point1, point2);
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getAverageEuclideanDistance (int x, int y, int offset_x, int offset_y, int max_steps) const
{
  float average_pixel_distance = 0.0f;
  float weight=0.0f;
  for (int i=0; i<max_steps; ++i)
  {
    int x1=x+i*offset_x,     y1=y+i*offset_y;
    int x2=x+ (i+1)*offset_x, y2=y+ (i+1)*offset_y;
    float pixel_distance = getEuclideanDistanceSquared (x1,y1,x2,y2);
    if (!pcl_isfinite (pixel_distance))
    {
      //std::cout << x<<","<<y<<"->"<<x2<<","<<y2<<": "<<pixel_distance<<"\n";
      if (i==0)
        return pixel_distance;
      else
        break;
    }
    //std::cout << x<<","<<y<<"->"<<x2<<","<<y2<<": "<<sqrtf (pixel_distance)<<"m\n";
    weight += 1.0f;
    average_pixel_distance += sqrtf (pixel_distance);
  }
  average_pixel_distance /= weight;
  //std::cout << x<<","<<y<<","<<offset_x<<","<<offset_y<<" => "<<average_pixel_distance<<"\n";
  return average_pixel_distance;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getImpactAngleBasedOnLocalNormal (int x, int y, int radius) const
{
  if (!isValid (x,y))
    return -std::numeric_limits<float>::infinity ();
  const PointWithRange& point = getPoint (x, y);
  int no_of_nearest_neighbors = static_cast<int> (pow (static_cast<double> ( (radius + 1.0)), 2.0));
  Eigen::Vector3f normal;
  if (!getNormalForClosestNeighbors (x, y, radius, point, no_of_nearest_neighbors, normal, 1))
    return -std::numeric_limits<float>::infinity ();
  return deg2radLocal (90.0f) - acosf (normal.dot ( (getSensorPos ()-getEigenVector3f (point)).normalized ()));
}


/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getNormal (int x, int y, int radius, Eigen::Vector3f& normal, int step_size) const
{
  pcl::VectorAverage3f vector_average;
  for (int y2=y-radius; y2<=y+radius; y2+=step_size)
  {
    for (int x2=x-radius; x2<=x+radius; x2+=step_size)
    {
      if (!isInImage (x2, y2))
        continue;
      const PointWithRange& point = getPoint (x2, y2);
      if (!pcl_isfinite (point.range))
        continue;
      vector_average.add (Eigen::Vector3f (point.x, point.y, point.z));
    }
  }
  if (vector_average.getNoOfSamples () < 3)
    return false;
  Eigen::Vector3f eigen_values, eigen_vector2, eigen_vector3;
  vector_average.doPCA (eigen_values, normal, eigen_vector2, eigen_vector3);
  if (normal.dot ( (getSensorPos ()-vector_average.getMean ()).normalized ()) < 0.0f)
    normal *= -1.0f;
  return true;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getNormalBasedAcutenessValue (int x, int y, int radius) const
{
  float impact_angle = getImpactAngleBasedOnLocalNormal (x, y, radius);
  if (pcl_isinf (impact_angle))
    return -std::numeric_limits<float>::infinity ();
  float ret = 1.0f - static_cast<float> ( (impact_angle / (0.5f * M_PI)));
  //std::cout << PVARAC (impact_angle)<<PVARN (ret);
  return ret;
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getNormalForClosestNeighbors (int x, int y, int radius, const PointWithRange& point,
                                               int no_of_nearest_neighbors, Eigen::Vector3f& normal, int step_size) const
{
  return getNormalForClosestNeighbors (x, y, radius, Eigen::Vector3f (point.x, point.y, point.z), no_of_nearest_neighbors, normal, NULL, step_size);
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getNormalForClosestNeighbors (int x, int y, Eigen::Vector3f& normal, int radius) const
{
  if (!isValid (x,y))
    return false;
  int no_of_nearest_neighbors = static_cast<int> (pow (static_cast<double> (radius + 1.0), 2.0));
  return getNormalForClosestNeighbors (x, y, radius, getPoint (x,y).getVector3fMap (), no_of_nearest_neighbors, normal);
}

namespace
{  // Anonymous namespace, so that this is only accessible in this file
  struct NeighborWithDistance
  {  // local struct to help us with sorting
    float distance;
    const PointWithRange* neighbor;
    bool operator < (const NeighborWithDistance& other) const { return distance<other.distance;}
  };
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getSurfaceInformation (int x, int y, int radius, const Eigen::Vector3f& point, int no_of_closest_neighbors, int step_size,
                                        float& max_closest_neighbor_distance_squared,
                                        Eigen::Vector3f& normal, Eigen::Vector3f& mean, Eigen::Vector3f& eigen_values,
                                        Eigen::Vector3f* normal_all_neighbors, Eigen::Vector3f* mean_all_neighbors,
                                        Eigen::Vector3f* eigen_values_all_neighbors) const
{
  max_closest_neighbor_distance_squared=0.0f;
  normal.setZero (); mean.setZero (); eigen_values.setZero ();
  if (normal_all_neighbors!=NULL)
    normal_all_neighbors->setZero ();
  if (mean_all_neighbors!=NULL)
    mean_all_neighbors->setZero ();
  if (eigen_values_all_neighbors!=NULL)
    eigen_values_all_neighbors->setZero ();
  
  int blocksize = static_cast<int> (pow (static_cast<double> ( (2.0 * radius + 1.0)), 2.0));
  
  PointWithRange given_point;
  given_point.x=point[0];  given_point.y=point[1];  given_point.z=point[2];
  
  std::vector<NeighborWithDistance> ordered_neighbors (blocksize);
  int neighbor_counter = 0;
  for (int y2=y-radius; y2<=y+radius; y2+=step_size)
  {
    for (int x2=x-radius; x2<=x+radius; x2+=step_size)
    {
      if (!isValid (x2, y2))
        continue;
      NeighborWithDistance& neighbor_with_distance = ordered_neighbors[neighbor_counter];
      neighbor_with_distance.neighbor = &getPoint (x2, y2);
      neighbor_with_distance.distance = squaredEuclideanDistance (given_point, *neighbor_with_distance.neighbor);
      ++neighbor_counter;
    }
  }
  no_of_closest_neighbors = (std::min) (neighbor_counter, no_of_closest_neighbors);
  
  std::sort (ordered_neighbors.begin (), ordered_neighbors.begin () + neighbor_counter);  // Normal sort seems to be the fastest method (faster than partial_sort)
  //std::stable_sort (ordered_neighbors, ordered_neighbors+neighbor_counter);
  //std::partial_sort (ordered_neighbors, ordered_neighbors+no_of_closest_neighbors, ordered_neighbors+neighbor_counter);
  
  max_closest_neighbor_distance_squared = ordered_neighbors[no_of_closest_neighbors-1].distance;
  //float max_distance_squared = max_closest_neighbor_distance_squared;
  float max_distance_squared = max_closest_neighbor_distance_squared*4.0f;  // Double the allowed distance value
  //max_closest_neighbor_distance_squared = max_distance_squared;
  
  pcl::VectorAverage3f vector_average;
  //for (int neighbor_idx=0; neighbor_idx<no_of_closest_neighbors; ++neighbor_idx)
  int neighbor_idx;
  for (neighbor_idx=0; neighbor_idx<neighbor_counter; ++neighbor_idx)
  {
    if (ordered_neighbors[neighbor_idx].distance > max_distance_squared)
      break;
    //std::cout << ordered_neighbors[neighbor_idx].distance<<"\n";
    vector_average.add (ordered_neighbors[neighbor_idx].neighbor->getVector3fMap ());
  }
  
  if (vector_average.getNoOfSamples () < 3)
    return false;
  //std::cout << PVARN (vector_average.getNoOfSamples ());
  Eigen::Vector3f eigen_vector2, eigen_vector3;
  vector_average.doPCA (eigen_values, normal, eigen_vector2, eigen_vector3);
  Eigen::Vector3f viewing_direction = (getSensorPos ()-point).normalized ();
  if (normal.dot (viewing_direction) < 0.0f)
    normal *= -1.0f;
  mean = vector_average.getMean ();
  
  if (normal_all_neighbors==NULL)
    return true;
  
  // Add remaining neighbors
  for (int neighbor_idx2=neighbor_idx; neighbor_idx2<neighbor_counter; ++neighbor_idx2)
    vector_average.add (ordered_neighbors[neighbor_idx2].neighbor->getVector3fMap ());
  
  vector_average.doPCA (*eigen_values_all_neighbors, *normal_all_neighbors, eigen_vector2, eigen_vector3);
  //std::cout << PVARN (vector_average.getNoOfSamples ())<<".\n";
  if (normal_all_neighbors->dot (viewing_direction) < 0.0f)
    *normal_all_neighbors *= -1.0f;
  *mean_all_neighbors = vector_average.getMean ();
  
  //std::cout << viewing_direction[0]<<","<<viewing_direction[1]<<","<<viewing_direction[2]<<"\n";
  
  return true;
}

/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getSquaredDistanceOfNthNeighbor (int x, int y, int radius, int n, int step_size) const
{
  const PointWithRange& point = getPoint (x, y);
  if (!pcl_isfinite (point.range))
    return -std::numeric_limits<float>::infinity ();
  
  int blocksize = static_cast<int> (pow (static_cast<double> (2.0 * radius + 1.0), 2.0));
  std::vector<float> neighbor_distances (blocksize);
  int neighbor_counter = 0;
  for (int y2=y-radius; y2<=y+radius; y2+=step_size)
  {
    for (int x2=x-radius; x2<=x+radius; x2+=step_size)
    {
      if (!isValid (x2, y2) || (x2==x&&y2==y))
        continue;
      const PointWithRange& neighbor = getPointNoCheck (x2,y2);
      float& neighbor_distance = neighbor_distances[neighbor_counter++];
      neighbor_distance = squaredEuclideanDistance (point, neighbor);
    }
  }
  std::sort (neighbor_distances.begin (), neighbor_distances.begin () + neighbor_counter);  // Normal sort seems to be
  // the fastest method (faster than partial_sort)
  n = (std::min) (neighbor_counter, n);
  return neighbor_distances[n-1];
}


/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getNormalForClosestNeighbors (int x, int y, int radius, const Eigen::Vector3f& point, int no_of_nearest_neighbors,
                                               Eigen::Vector3f& normal, Eigen::Vector3f* point_on_plane, int step_size) const
{
  Eigen::Vector3f mean, eigen_values;
  float used_squared_max_distance;
  bool ret = getSurfaceInformation (x, y, radius, point, no_of_nearest_neighbors, step_size, used_squared_max_distance,
                                    normal, mean, eigen_values);
  
  if (ret)
  {
    if (point_on_plane != NULL)
      *point_on_plane = (normal.dot (mean) - normal.dot (point))*normal + point;
  }
  return ret;
}


/////////////////////////////////////////////////////////////////////////
float
RangeImageLocal::getCurvature (int x, int y, int radius, int step_size) const
{
  pcl::VectorAverage3f vector_average;
  for (int y2=y-radius; y2<=y+radius; y2+=step_size)
  {
    for (int x2=x-radius; x2<=x+radius; x2+=step_size)
    {
      if (!isInImage (x2, y2))
        continue;
      const PointWithRange& point = getPoint (x2, y2);
      if (!pcl_isfinite (point.range))
        continue;
      vector_average.add (Eigen::Vector3f (point.x, point.y, point.z));
    }
  }
  if (vector_average.getNoOfSamples () < 3)
    return false;
  Eigen::Vector3f eigen_values;
  vector_average.doPCA (eigen_values);
  return eigen_values[0]/eigen_values.sum ();
}


/////////////////////////////////////////////////////////////////////////
template <typename PointCloudTypeWithViewpoints> Eigen::Vector3f
RangeImageLocal::getAverageViewPoint (const PointCloudTypeWithViewpoints& point_cloud)
{
  Eigen::Vector3f average_viewpoint (0,0,0);
  int point_counter = 0;
  for (unsigned int point_idx=0; point_idx<point_cloud.points.size (); ++point_idx)
  {
    const typename PointCloudTypeWithViewpoints::PointType& point = point_cloud.points[point_idx];
    if (!pcl_isfinite (point.vp_x) || !pcl_isfinite (point.vp_y) || !pcl_isfinite (point.vp_z))
      continue;
    average_viewpoint[0] += point.vp_x;
    average_viewpoint[1] += point.vp_y;
    average_viewpoint[2] += point.vp_z;
    ++point_counter;
  }
  average_viewpoint /= point_counter;
  
  return average_viewpoint;
}

/////////////////////////////////////////////////////////////////////////
bool
RangeImageLocal::getViewingDirection (int x, int y, Eigen::Vector3f& viewing_direction) const
{
  if (!isValid (x, y))
    return false;
  viewing_direction = (getPoint (x,y).getVector3fMap ()-getSensorPos ()).normalized ();
  return true;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getViewingDirection (const Eigen::Vector3f& point, Eigen::Vector3f& viewing_direction) const
{
  viewing_direction = (point-getSensorPos ()).normalized ();
}

/////////////////////////////////////////////////////////////////////////
Eigen::Affine3f
RangeImageLocal::getTransformationToViewerCoordinateFrame (const Eigen::Vector3f& point) const
{
  Eigen::Affine3f transformation;
  getTransformationToViewerCoordinateFrame (point, transformation);
  return transformation;
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getTransformationToViewerCoordinateFrame (const Eigen::Vector3f& point, Eigen::Affine3f& transformation) const
{
  Eigen::Vector3f viewing_direction = (point-getSensorPos ()).normalized ();
  pcl::getTransformationFromTwoUnitVectorsAndOrigin (Eigen::Vector3f (0.0f, -1.0f, 0.0f), viewing_direction, point, transformation);
}

/////////////////////////////////////////////////////////////////////////
void
RangeImageLocal::getRotationToViewerCoordinateFrame (const Eigen::Vector3f& point, Eigen::Affine3f& transformation) const
{
  Eigen::Vector3f viewing_direction = (point-getSensorPos ()).normalized ();
  pcl::getTransformationFromTwoUnitVectors (Eigen::Vector3f (0.0f, -1.0f, 0.0f), viewing_direction, transformation);
}

/////////////////////////////////////////////////////////////////////////
inline void
RangeImageLocal::setAngularResolution (float angular_resolution)
{
  angular_resolution_x_ = angular_resolution_y_ = angular_resolution;
  angular_resolution_x_reciprocal_ = angular_resolution_y_reciprocal_ = 1.0f / angular_resolution;
}

/////////////////////////////////////////////////////////////////////////
inline void
RangeImageLocal::setAngularResolution (float angular_resolution_x, float angular_resolution_y)
{
  angular_resolution_x_ = angular_resolution_x;
  angular_resolution_x_reciprocal_ = 1.0f / angular_resolution_x_;
  angular_resolution_y_ = angular_resolution_y;
  angular_resolution_y_reciprocal_ = 1.0f / angular_resolution_y_;
}

inline void
RangeImageLocal::setTransformationToRangeImageSystem (const Eigen::Affine3f& to_range_image_system)
{
  to_range_image_system_ = to_range_image_system;
  to_world_system_ = to_range_image_system_.inverse ();
}

inline void
RangeImageLocal::getAngularResolution (float& angular_resolution_x, float& angular_resolution_y) const
{  
  angular_resolution_x = angular_resolution_x_;
  angular_resolution_y = angular_resolution_y_;
}

