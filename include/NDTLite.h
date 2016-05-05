/*
 * IRON Reference Implementation:
 * *************************************************************************************************************
 * Schmiedel, Th., Einhorn, E., Gross, H.-M.
 * IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching and Its Application to Robot Localization.
 * IEEE/RSJ Int. Conf. on Intelligent Robots and Systems (IROS), Hamburg, Germany, 2015
 * *************************************************************************************************************
 * Copyright (C) 2015, Thomas Schmiedel (thomas.schmiedel@tu-ilmenau.de)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/


#ifndef _NDTLITE_H_
#define _NDTLITE_H_


///////////////////////////////
// includes
///////////////////////////////
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <float.h>
#include <stdint.h>


///////////////////////////////
// namespaces
///////////////////////////////
namespace IRON
{

////////////////////////////////////////////////////
// To understand how to choose parameters, here's
// a short explanation of the NDTLite module:
// Basically what it does is to quickly transform
// a point cloud into an NDT map. With the current
// set of parameters - run on a modern i7 CPU -
// it takes only about 3ms to transform a Kinect
// point cloud with sufficient accuracy.
//
// This works as follows:
// - a fast xorshift random generator is used to
//   draw points from the input point cloud,
//   ("subsamplingFactor" tells what fraction of
//    the original amount of points is drawn)
// - if noise is enabled, it will now be added
//   to the point values in sensor direction
// - those noisy points are inserted into 3d
//   grid cells, as long as they have less than
//   "maxPointsPerCell"
// - after all points have been processed,
//   grid cells are selected that contain at least
//   "minPointsPerCell"; they are transformed into
//   NDT cells and later comprise the NDT map
// - the transformation is achieved by computing the
//   mean and covariance matrix of all points within
//   one cell
// - "cellSize" specifies the dimension of one NDT
//   cell, its default size is 0.1m x 0.1m x 0.1m
//
// PLEASE NOTE:
// This NDT module is currently assumed to process
// dense point clouds only, which means their points
// are gathered in a rather limited area. This is
// perfectly true for Kinect/AsusXtion point clouds
// and other data captured by depth sensors in indoor
// environments. However, for large scans from whole
// buildings etc., it will (depending on "cellSize")
// likely run out of memory.
// So far, there was no need to work with such large
// scans. If, however, there's demand to have this
// feature (or others) implemented,
// do not hesitate to contact me:
//
//   thomas.schmiedel@tu-ilmenau.de
//
////////////////////////////////////////////////////

struct NDTLiteConfig
{
    NDTLiteConfig() : maxPointsPerCell(100u),     // do not add more points within one cell
                      minPointsPerCell(5u),       // skip cells with fewer points
                      preAllocCells(0u),          // pre allocate grid during construction
                      normalDistTableSize(1000u), // reduce for faster initialization
                      normalDistLoops(12u),       // reduce for faster initialization
                      subsamplingFactor(0.2f),    // subsample point cloud (1.0=all points, 0.5=half of points, ...)
                      sensorNoise(0.001f),        // will be added in sensor direction, if noise is enabled
                      cellSize(0.1f),             // NDT cell size
                      clippingDistance(5.0f)      // throw away points further away from the sensor
    { }

    uint maxPointsPerCell,
         minPointsPerCell,
         preAllocCells,
         normalDistTableSize,
         normalDistLoops;
    float subsamplingFactor,
          sensorNoise,
          cellSize,
          clippingDistance;
};


///////////////////////////////
// NDT cell lite class
///////////////////////////////
template <typename NumericalType>
struct NDTCellLite
{
    Eigen::Matrix<NumericalType, 3, 3> cov;
    Eigen::Matrix<NumericalType, 3, 1> mu;
};


///////////////////////////////
// NDT map lite class
///////////////////////////////
template <typename NumericalType>
struct NDTMapLite
{
private:
    NDTMapLite(const NDTMapLite &other);
    NDTMapLite &operator=(const NDTMapLite &other);

public:
    NDTMapLite()
    { }
    std::vector<NDTCellLite<NumericalType> > cells;
};


typedef NDTMapLite<float>  NDTMapLitef;
typedef NDTMapLite<double> NDTMapLited;


///////////////////////////////
// simple random number gen.
///////////////////////////////
class XorShift
{
private:
    XorShift(const XorShift &other);
    XorShift &operator=(const XorShift &other);

public:
    XorShift(uint64_t seed)
    {
        setSeed(seed);
    }
    XorShift()
    {
        setSeed(time(NULL));
    }
    void setSeed(uint64_t seed)
    {
        mSeed = (seed == 0) ? 0xdeadbeefu : seed;
    }
    uint64_t rand()
    {
        mSeed ^= mSeed << 13;
        mSeed ^= mSeed >> 7;
        mSeed ^= mSeed << 17;
        return mSeed;
    }
private:
    uint64_t mSeed;
};


///////////////////////////////
// normal distribution gen
///////////////////////////////
template <typename NumericalType>
class NormalDistGenerator
{
private:
    NormalDistGenerator(const NormalDistGenerator &other);
    NormalDistGenerator &operator=(const NormalDistGenerator &other);

public:
    NormalDistGenerator(NumericalType mu    = (NumericalType)0.0,
                        NumericalType sigma = (NumericalType)1.0,
                        uint loops = 12u) : mMu(mu),
                                            mSigma(sigma),
                                            mLoops(loops),
                                            mScale(std::sqrt((NumericalType)12.0 / loops))
    { }

    void setMu(NumericalType val)
    {
        mMu = val;
    }

    void setSigma(NumericalType val)
    {
        mSigma = val;
    }

    void setLoops(uint val)
    {
        mLoops = val;
        mScale = std::sqrt((NumericalType)12.0 / val);
    }

    void setSeed(uint64_t seed)
    {
        mRnd.setSeed(seed);
    }

    NumericalType rand()
    {
        NumericalType sum = (NumericalType)0.0;
        for (uint i = 0; i < mLoops; ++i)
        {
            sum += uniform();
        }
        return (sum - (NumericalType)0.5 * mLoops) * mScale * mSigma + mMu;
    }

private:
    NumericalType uniform()
    {
        return (mRnd.rand() % 1000000u) * (NumericalType)0.000001;
    }

    XorShift mRnd;
    NumericalType mMu, mSigma;
    uint mLoops;
    NumericalType mScale;
};


///////////////////////////////
// simple look up table
///////////////////////////////
template <typename Type>
class CircularLUT
{
private:
    CircularLUT(const CircularLUT &other);
    CircularLUT &operator=(const CircularLUT &other);

public:
    CircularLUT(uint cnt = 0) : mPos(0)
    {
        reserve(cnt);
    }

    void reserve(uint cnt)
    {
        mBuffer.reserve(cnt);
    }

    void clear()
    {
        mBuffer.clear();
    }

    void insert(const Type &value)
    {
        mBuffer.push_back(value);
    }

    Type next()
    {
        if (mPos == mBuffer.size())
        {
            mPos = 0;
        }
        return mBuffer[mPos++];
    }

protected:
    std::vector<Type> mBuffer;
    uint mPos;
};


///////////////////////////////
// ndt map lite creator
///////////////////////////////
enum IRONNoise
{
    IRON_NO_NOISE,
    IRON_CONSTANT_NOISE,
    IRON_PRIMESENSE_NOISE
};
template <typename NumericalType,
          typename ContainerType,
          typename PointAccessor,
          IRONNoise NoiseType = IRON_NO_NOISE>
class NDTMapLiteCreator
{
private:
    NDTMapLiteCreator(const NDTMapLiteCreator &other);
    NDTMapLiteCreator &operator=(const NDTMapLiteCreator &other);

public:
    NDTMapLiteCreator(const NDTLiteConfig &config) : mConfig(config), mAlloc(0u)
    {
        // create LUT for normal distribution
        NormalDistGenerator<NumericalType> mNGen;
        // fix seeds to reproduce behavior
        mRnd.setSeed(0x12345678);
        mNGen.setSeed(0x87654321);
        // normal dist gen loops
        if (NoiseType == IRON_CONSTANT_NOISE || NoiseType == IRON_PRIMESENSE_NOISE)
        {
            mNGen.setLoops(config.normalDistLoops);
            mNlut.reserve(config.normalDistTableSize);
            for (uint i = 0; i < config.normalDistTableSize; ++i)
            {
                mNlut.insert(mNGen.rand());
            }
        }
        // do initial allocation
        realloc(config.preAllocCells);
    }

    void realloc(uint cells)
    {
        mGrid.resize(cells);
        for (uint i = 0; i < cells; ++i)
        {
            mGrid[i].reserve(mConfig.maxPointsPerCell);
        }
        mAlloc = cells;
    }

    void createMapFromPointValues(NDTMapLite<NumericalType> *map,
                                  const ContainerType &container,
                                  const Eigen::Transform<NumericalType, 3, Eigen::Affine> &sensorPose)
    {
        // clear
        map->cells.clear();

        // subsampling
        const uint limit = mConfig.subsamplingFactor * container.size();

        // alloc space for noisy points
        mNoiseVec.clear();
        mNoiseVec.reserve(limit);

        // find min and max values
        NumericalType minx = FLT_MAX,
                      miny = FLT_MAX,
                      minz = FLT_MAX,
                      maxx = -FLT_MAX,
                      maxy = -FLT_MAX,
                      maxz = -FLT_MAX;

        Eigen::Matrix<NumericalType, 3, 1> dir;
        const NumericalType clipDistSqr = mConfig.clippingDistance * mConfig.clippingDistance;
        for (uint i = 0; i < limit; ++i)
        {
            // get random point
            const uint idx = mRnd.rand() % container.size();
            const auto &pt = container[idx];

            // create vector from it
            const Eigen::Matrix<NumericalType, 3, 1> rawpt(PointAccessor::x(pt),
                                                           PointAccessor::y(pt),
                                                           PointAccessor::z(pt));

            // compute squared sensor distance and check clippingDist
            const NumericalType distSqr = norm_sqr(rawpt);
            if (distSqr > clipDistSqr) continue;

            // add noise; depending on template param, 2 of 3 will be optimized out
            if (NoiseType == IRON_CONSTANT_NOISE)
            {
                dir = rawpt;
                normalize_fast_from_sqr(dir, distSqr);
                mNoiseVec.push_back(sensorPose * (rawpt + mNlut.next() * mConfig.sensorNoise * dir));
            }
            if (NoiseType == IRON_PRIMESENSE_NOISE)
            {
                dir = rawpt;
                normalize_fast_from_sqr(dir, distSqr);
                const NumericalType sigma = mConfig.sensorNoise * distSqr;
                mNoiseVec.push_back(sensorPose * (rawpt + mNlut.next() * sigma * dir));
            }
            if (NoiseType == IRON_NO_NOISE)
            {
                // pretty boring here
                mNoiseVec.push_back(sensorPose * rawpt);
            }

            // update limits
            const Eigen::Matrix<NumericalType, 3, 1> &noisept = mNoiseVec.back();
            if (noisept(0) < minx) minx = noisept(0);
            if (noisept(1) < miny) miny = noisept(1);
            if (noisept(2) < minz) minz = noisept(2);
            if (noisept(0) > maxx) maxx = noisept(0);
            if (noisept(1) > maxy) maxy = noisept(1);
            if (noisept(2) > maxz) maxz = noisept(2);
        }

        // any point left?
        if (mNoiseVec.size() == 0)
        {
            std::cerr << "all points were filtered out, increase clipping distance and/or check sensor pose\n";
            return;
        }

        // grid size
        const NumericalType invcellsize = (NumericalType)1.0 / mConfig.cellSize;

        // align min values to global grid
        minx = std::floor(minx * invcellsize) * mConfig.cellSize;
        miny = std::floor(miny * invcellsize) * mConfig.cellSize;
        minz = std::floor(minz * invcellsize) * mConfig.cellSize;

        // compute cells
        const uint gridx = std::ceil((maxx - minx) * invcellsize),
                   gridy = std::ceil((maxy - miny) * invcellsize),
                   gridz = std::ceil((maxz - minz) * invcellsize),
                   cells = gridx * gridy * gridz;

        // alloc
        if (cells > mAlloc)
        {
            //std::cerr << "reallocating..." << gridx << "x" << gridy << "x" << gridz << std::endl;
            realloc(cells);
        }
        // clear ptr lists
        for (uint i = 0; i < cells; ++i)
        {
            mGrid[i].clear();
        }

        // gather points
        const uint layer = gridx * gridy;
        for (uint i = 0; i < mNoiseVec.size(); ++i)
        {
            // get noisy point
            const Eigen::Matrix<NumericalType, 3, 1> &noisept = mNoiseVec[i];

            // get cell coords
            const uint cx = static_cast<uint>((noisept(0) - minx) * invcellsize),
                       cy = static_cast<uint>((noisept(1) - miny) * invcellsize),
                       cz = static_cast<uint>((noisept(2) - minz) * invcellsize);

            // compute array index
            const uint gridindex = cz * layer + cy * gridx + cx;

            // add ptr
            if (mGrid[gridindex].size() == mConfig.maxPointsPerCell) continue;
            mGrid[gridindex].push_back(&noisept);
        }

        // reserve cells
        map->cells.reserve(cells);

        // compute mean and covariance
        for (uint i = 0; i < cells; ++i)
        {
            // skip with few points
            if (mGrid[i].size() < mConfig.minPointsPerCell) continue;
            map->cells.push_back(NDTCellLite<NumericalType>());
            map->cells.back().mu = mean(mGrid[i]);
            map->cells.back().cov = cov(mGrid[i], map->cells.back().mu);
        }
    }

private:
    // quick vector norm
    void normalize_fast(Eigen::Matrix<NumericalType, 3, 1> &vec) const
    {
        vec *= invsqrt_fast(vec(0) * vec(0) + vec(1) * vec(1) + vec(2) * vec(2));
    }

    // normalize from given dist squared
    void normalize_fast_from_sqr(Eigen::Matrix<NumericalType, 3, 1> &vec, const NumericalType normsqr) const
    {
        vec *= invsqrt_fast(normsqr);
    }

    // norm squared
    NumericalType norm_sqr(const Eigen::Matrix<NumericalType, 3, 1> &vec) const
    {
        return vec(0) * vec(0) + vec(1) * vec(1) + vec(2) * vec(2);
    }

    // must remain float (=32 bit)
    float invsqrt_fast(const float val) const
    {
        union { float f; uint32_t u; } tmp;
        tmp.f = val;
        tmp.u = 0x5f3759df - (tmp.u >> 1);
        return tmp.f;
    }

    // compute mean
    Eigen::Matrix<NumericalType, 3, 1> mean(const std::vector<const Eigen::Matrix<NumericalType, 3, 1> *> &samples) const
    {
        Eigen::Matrix<NumericalType, 3, 1> res((NumericalType)0.0, (NumericalType)0.0, (NumericalType)0.0);
        for (uint i = 0; i < samples.size(); ++i)
        {
            res += *(samples[i]);
        }
        return res / samples.size();
    }

    // compute cov
    Eigen::Matrix<NumericalType, 3, 3> cov(const std::vector<const Eigen::Matrix<NumericalType, 3, 1> *> &samples,
                                           const Eigen::Matrix<NumericalType, 3, 1> &mean) const
    {
        Eigen::Matrix<NumericalType, 3, 3> res;
        res.setZero();
        Eigen::Matrix<NumericalType, 3, 1> d;
        for (uint i = 0; i < samples.size(); ++i)
        {
            d = *(samples[i]) - mean;
            res += d * d.transpose();
        }
        return res / samples.size();
    }

protected:
    NDTLiteConfig mConfig;
    XorShift mRnd;
    CircularLUT<NumericalType> mNlut;
    std::vector<Eigen::Matrix<NumericalType, 3, 1> > mNoiseVec; // for Vector3 we can do this
    std::vector<std::vector<const Eigen::Matrix<NumericalType, 3, 1> *> > mGrid;
    uint mAlloc;
};


// end of namespace
}

#endif
