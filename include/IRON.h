/*
 * IRON Reference Implementation:
 * *************************************************************************************************************
 * Schmiedel, Th., Einhorn, E., Gross, H.-M.
 * IRON: A Fast Interest Point Descriptor for Robust NDT-Map Matching and Its Application to Robot Localization.
 * IEEE/RSJ Int. Conf. on Intelligent Robots and Systems (IROS), Hamburg, Germany, 2015
 * *************************************************************************************************************
 * Copyright (C) 2015, Thomas Schmiedel
 * 
 * THE BSD LICENSE
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


// TODO: uint -> size_t


#ifndef IRON_H_
#define IRON_H_


////////////////////////////////////////
// INCLUDES
////////////////////////////////////////
#include "NDTLite.h"
#include "KNNAdapter.h"
#include <algorithm>
#include <iomanip>


////////////////////////////////////////
// NAMESPACE
////////////////////////////////////////
namespace IRON
{
#define IRON_LOG(_x) (std::log((NumericalType)_x) * (NumericalType)3.32192809489) // base e


////////////////////////////////////////
// PARAMETERS
////////////////////////////////////////

// DIST_EXP: Influences mapping from
// distances to histogram bins
// 1: Linear mapping
// 2: Quadratic mapping [default]
// 4: Quartic mapping
// The linear model may lead to sparsely
// populated first distance rows
// > 2 is the fastest
#define DIST_EXP 2

// params
struct IRONConfig
{
    IRONConfig() : neighborSearchRadius(0.5f), // radius in m around the base NDT-cell
                                               // to be considered for neighbor search
                   entropyThreshold(0.55f),    // normalized entropy threshold used for
                                               // keypoint classification; a lower threshold
                                               // results in more keypoints but is slower
                   maxEV0EV1Ratio(0.9f),       // maximum ratio between eigenvalue 0 and 1 (EV0 is the smallest)
                   minEV1EV2Ratio(0.1f),       // minimum ratio between eigenvalue 1 and 2 (EV2 is the biggest)
                                               // > both params together will prevent spherically and needle shaped NDT cells
                   matchingTolerance(0.05f),   // maximum distance in m for a sample point to be considered an inlier
                                               // by RANSAC; 0.5 * NDT-cellSize - 1.0 * NDT-cellSize is a good choice
                   bfMatchesLimit(1.0f),       // take only N best matches for further processing (e.g. 0.5 = best half of matches)
                   distanceBins(3u),           // descriptor histogram bins for encoding distance values
                   angleBins(3u),              // descriptor histogram bins for encoding angular values
                   minNeighborCount(5u),       // skip NDT-cells with less neighbors than N inside the search radius
                   ransacLoops(2000u),         // amount of loops for outlier rejection
                   acosTableSize(5000u)        // acos look-up table has N entries
    { }

    float neighborSearchRadius,
          entropyThreshold,
          maxEV0EV1Ratio,
          minEV1EV2Ratio,
          matchingTolerance,
          bfMatchesLimit;
    uint  distanceBins,
          angleBins,
          minNeighborCount,
          ransacLoops,
          acosTableSize;
};


////////////////////////////////////////
// IRON DESCRIPTOR
////////////////////////////////////////
template <typename NumericalType>
struct IRONDescriptor
{
    const NDTCellLite<NumericalType> *cellPtr;
    NumericalType *h;
    Eigen::Matrix<NumericalType, 3, 1> normal;
    uint k;
    NumericalType entropy;
    IRONDescriptor(const NDTCellLite<NumericalType> *ptr) : cellPtr(ptr), h(NULL), k(0), entropy((NumericalType)0.0)
    { }
    const Eigen::Matrix<NumericalType, 3, 1> &mu()  const { return cellPtr->mu;  }
    const Eigen::Matrix<NumericalType, 3, 3> &cov() const { return cellPtr->cov; }
};

template <typename NumericalType>
struct IRONDescriptorVector
{
    typedef std::vector<IRONDescriptor<NumericalType> > type;
};

typedef IRONDescriptorVector<float>::type IRONDescriptorVectorf;
typedef IRONDescriptorVector<double>::type IRONDescriptorVectord;


////////////////////////////////////////
// ACCESSOR FOR DESCRIPTOR MU
////////////////////////////////////////
template <typename NumericalType>
struct IRONAccessor
{
private:
    IRONAccessor();
    IRONAccessor(const IRONAccessor &other);
    IRONAccessor &operator=(const IRONAccessor &other);

public:
    static NumericalType at(const IRONDescriptor<NumericalType> &item, uint i) { return item.mu()[i]; }
};


//////////////////////////////////////
// acos look up table
//////////////////////////////////////
template <typename NumericalType>
class ACosLUT
{
private:
    ACosLUT(const ACosLUT &other);
    ACosLUT &operator=(const ACosLUT &other);

public:
    ACosLUT(uint elements = 2000u)
    {
        setTableSize(elements);
    }

    void setTableSize(uint elements)
    {
        mElements = elements;
        mStep = (NumericalType)2.0 / elements;
        mStepInv = elements / (NumericalType)2.0;
    }

    void create()
    {
        mBuffer.resize(mElements + 1u);
        for (uint i = 0u; i <= mElements; ++i)
        {
            mBuffer[i] = std::acos(mStep * i - (NumericalType)1.0);
        }
    }

    NumericalType acos(NumericalType value) const
    {
        if (value < (NumericalType)-1.0) { value = (NumericalType)-1.0; }
        else if (value > (NumericalType)1.0) { value = (NumericalType)1.0; }
        return mBuffer[(uint)((value + (NumericalType)1.0) * mStepInv)];
    }

protected:
    uint mElements;
    NumericalType mStep,
                  mStepInv;
    std::vector<NumericalType> mBuffer;
};


typedef ACosLUT<float> ACosLUTf;
typedef ACosLUT<double> ACosLUTd;


//////////////////////////////////////
// TMP STRUCT FOR MATCHING
//////////////////////////////////////
template <typename NumericalType>
struct CostIdxPair
{
    NumericalType cost;
    uint from, to;
    CostIdxPair() : cost((NumericalType)0.0), from(0), to(0)
    { }
};


//////////////////////////////////////
// STRUCT FOR MATCHES
//////////////////////////////////////
template <typename NumericalType>
struct IRONMatch
{
    NumericalType cost;
    const IRONDescriptor<NumericalType> *from, *to;
    IRONMatch() : cost((NumericalType)0.0), from(NULL), to(NULL)
    { }
};


template <typename NumericalType>
struct IRONMatchVector
{
    typedef std::vector<IRONMatch<NumericalType> > type;
};


typedef IRONMatchVector<float>::type IRONMatchVectorf;
typedef IRONMatchVector<double>::type IRONMatchVectord;


//////////////////////////////////////
// TRANSFORM RESULT
//////////////////////////////////////
template <typename NumericalType>
struct IRONTransformResult
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    bool success;
    Eigen::Transform<NumericalType, 3, Eigen::Affine> tf;
    IRONTransformResult() : success(false)
    {
        tf.setIdentity();
    }
};


typedef IRONTransformResult<float> IRONTransformResultf;
typedef IRONTransformResult<double> IRONTransformResultd;


//////////////////////////////////////
// COMPUTE DESCRIPTORS AND KEYPOINTS
//////////////////////////////////////
template <typename NumericalType>
class IRONEngine
{
private:
    IRONEngine(const IRONEngine &other);
    IRONEngine &operator=(const IRONEngine &other);

public:
    // initialization
    IRONEngine(const IRONConfig &config) : mConfig(config)
    {
        // acos look up table
        mACLUT.setTableSize(config.acosTableSize);
        mACLUT.create();
        // fix random seed for deterministic behavior
        mRnd.setSeed(0xfedcba01u);
    }

    // compute descriptors and keypoints
    void computeDescriptors(typename IRONDescriptorVector<NumericalType>::type *vec,
                            const NDTMapLite<NumericalType> &map) const
    {
        // prepare
        releaseDescriptors(vec);
        const uint bins = 2u * mConfig.angleBins * mConfig.distanceBins;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<NumericalType, 3, 3> > esolver;
        const NumericalType minev = mConfig.minEV1EV2Ratio,
                            maxev = mConfig.maxEV0EV1Ratio;
        vec->reserve(map.cells.size());

        // discard degenerate cells
        for (uint i = 0; i < map.cells.size(); ++i)
        {
            // compute eigenvalues
            esolver.compute(map.cells[i].cov); //TODO: evaluate computeDirect()

            // check cell shape
            if (esolver.eigenvalues()(0) / esolver.eigenvalues()(1) > maxev ||
                esolver.eigenvalues()(1) / esolver.eigenvalues()(2) < minev)
            {
                continue;
            }

            // alloc
            vec->push_back(IRONDescriptor<NumericalType>(&(map.cells[i])));
            vec->back().h = new NumericalType[bins];
            memset(vec->back().h, 0, sizeof(NumericalType) * bins);
            vec->back().k = 0;
            vec->back().normal = esolver.eigenvectors().col(0);
        }

        // check enough valid cells
        //std::cerr << "got " << validCells.size() << " valid cells\n";
        if (vec->size() < 5u)
        {
            std::cerr << "not enough ndt cells left after shape filtering\n";
            return;
        }

        // prepare for nearest neighbor search
        GenericKNNAdapter<3, typename IRONDescriptorVector<NumericalType>::type,
                          IRONAccessor<NumericalType>, NumericalType> knn;
        typedef typename GenericKNNAdapter<3, typename IRONDescriptorVector<NumericalType>::type,
                                           IRONAccessor<NumericalType>, NumericalType>::IdxDistPairs IdxDistPairs;
        knn.attach(vec);

        // prepare some parameters
        const NumericalType radiusSqr = mConfig.neighborSearchRadius * mConfig.neighborSearchRadius;
#if   DIST_EXP == 1
        const NumericalType distmax = cfg.neighborSearchRadius;
#elif DIST_EXP == 2
        const NumericalType distmax = mConfig.neighborSearchRadius * mConfig.neighborSearchRadius;
#elif DIST_EXP == 4
        const NumericalType distmax = mConfig.neighborSearchRadius * mConfig.neighborSearchRadius *
                                      mConfig.neighborSearchRadius * mConfig.neighborSearchRadius;
#else
#error "Unsupported DIST_EXP"
#endif
        const NumericalType dfac = mConfig.distanceBins / (distmax + (NumericalType)0.000001),
                            afac = mConfig.angleBins / ((NumericalType)0.5 * M_PI + (NumericalType)0.000001);
        const uint minneigh = mConfig.minNeighborCount,
                   distBins = mConfig.distanceBins,
                   angBins  = mConfig.angleBins,
                   shift    = mConfig.angleBins * mConfig.distanceBins;

        // some temp vars
        Eigen::Matrix<NumericalType, 3, 1> dvec;
        uint *freq = new uint[distBins];

        // loop valid cells
        for (uint i = 0; i < vec->size(); ++i)
        {
            // get neighbors inside radius
            const IdxDistPairs &neighbors = knn.radiusSearch(vec->at(i), radiusSqr);
            if (neighbors.size() < minneigh + 1u) continue;
            vec->at(i).k = neighbors.size() - 1u;

            // reference to normal vector
            const Eigen::Matrix<NumericalType, 3, 1> &basenvec = vec->at(i).normal;

            // fill histograms
            memset(freq, 0, distBins * sizeof(uint));
            for (uint k = 0; k < neighbors.size(); ++k)
            {
                // check index
                const uint idx = neighbors[k].first;
                if (idx == i) continue; // this is the base node

                // now compute histogram values
#if   DIST_EXP == 1
                const NumericalType dist = std::sqrt(neighbors[k].second);
#elif DIST_EXP == 2
                const NumericalType dist = neighbors[k].second;
#elif DIST_EXP == 4
                const NumericalType dist = neighbors[k].second * neighbors[k].second;
#else
#error "Unsupported DIST_EXP"
#endif

                // direction vector
                dvec = vec->at(i).mu() - vec->at(idx).mu();
                dvec.normalize();

                // angles
                const NumericalType nangle = limitedAngle(basenvec, vec->at(idx).normal),
                                    dangle = limitedAngle(basenvec, dvec);

                // discretization
                // distance between base and neighbor mu
                const uint dindex = static_cast<uint>(dist * dfac);
                freq[dindex] += 1u;
                // angle between surface normals
                uint aindex = static_cast<uint>(nangle * afac);
                const uint tmpidx = dindex * angBins;
                vec->at(i).h[tmpidx + aindex] += (NumericalType)1.0;
                // angle between base normal and direction to neighbor
                aindex = static_cast<uint>(dangle * afac);
                vec->at(i).h[shift + tmpidx + aindex] += (NumericalType)1.0;
            }

            // normalize histograms
            uint fcnt = 0;
            for (uint x = 0; x < distBins; ++x)
            {
                // count filled rows
                if (freq[x] != 0) ++fcnt;
            }
            for (uint x = 0; x < distBins; ++x)
            {
                // skip empty row
                if (freq[x] == 0) continue;
                // norm value
                const NumericalType dnorm = (NumericalType)1.0 / (fcnt * freq[x]);
                for (uint y = 0; y < angBins; ++y)
                {
                    const uint tmpidx = x * angBins + y;
                    vec->at(i).h[tmpidx] *= dnorm;
                    vec->at(i).h[shift + tmpidx] *= dnorm;
                }
            }
        }
        delete[] freq;

        // nn search is done, now throw away nodes with k=0 or low a-histo-entropy
        NumericalType entropy, p;
        const NumericalType entthres = mConfig.entropyThreshold,
                            norment  = (NumericalType)-1.0 / IRON_LOG((NumericalType)shift);
        for (uint i = 0; i < vec->size(); )
        {
            // if not enough neighbors, skip
            if (vec->at(i).k == 0)
            {
                remove(vec, i);
                continue;
            }

            // compute entropy
            entropy = (NumericalType)0.0;
            for (uint e = 0; e < shift; ++e)
            {
                p = vec->at(i).h[e];
                if (p == (NumericalType)0.0) continue;
                entropy += p * IRON_LOG(p);
            }

            // normalize
            entropy *= norment;

            // check
            //std::cerr << "ent=" << entropy << std::endl;
            if (entropy < entthres)
            {
                remove(vec, i);
                continue;
            }

            // store entropy
            vec->at(i).entropy = entropy;
            ++i;
        }
    }

    // free memory
    void releaseDescriptors(typename IRONDescriptorVector<NumericalType>::type *vec) const
    {
        for (int i = (int)(vec->size()) - 1; i >= 0; --i)
        {
            remove(vec, i);
        }
    }

    // print descriptor information
    void printDescriptor(const IRONDescriptor<NumericalType> &vec) const
    {
        std::cerr << std::setw(4) << std::setprecision(2) << std::fixed;
        std::cerr << "MU    : " << vec.mu()(0) << ", " << vec.mu()(1) << ", " << vec.mu()(2) << std::endl
                  << "NORMAL: " << vec.normal(0) << ", " << vec.normal(1) << ", " << vec.normal(2) << std::endl
                  << "NEIGHB: " << vec.k << std::endl
                  << "ENTROP: " << vec.entropy << std::endl;
        if (vec.k == 0) return;
        std::cerr << "DESCRIPTOR=\n";
        const uint bins = mConfig.angleBins * mConfig.distanceBins;
        for (uint k = 0; k < mConfig.distanceBins; ++k)
        {
            for (uint i = 0; i < mConfig.angleBins; ++i)
            {
                std::cerr << vec.h[k * mConfig.angleBins + i] << " ";
            }
            for (uint i = 0; i < mConfig.angleBins; ++i)
            {
                std::cerr << vec.h[bins + k * mConfig.angleBins + i] << " ";
            }
            std::cerr << std::endl;
        }
    }

    // brute force matching
    void computeMatches(typename IRONMatchVector<NumericalType>::type *matches,
                        const typename IRONDescriptorVector<NumericalType>::type &from,
                        const typename IRONDescriptorVector<NumericalType>::type &to)
    {
        // prepare
        matches->clear();

        // alloc
        mCostVector.clear();
        mCostVector.reserve(from.size() * to.size());
        NumericalType tmp;
        const uint size = 2u * mConfig.angleBins * mConfig.distanceBins;

        // compute cross wise costs between all descriptors
        for (uint i = 0; i < from.size(); ++i)
        {
            for (uint k = 0; k < to.size(); ++k)
            {
                // compute cost
                CostIdxPair<NumericalType> cip;
                for (uint c = 0; c < size; ++c)
                {
                    tmp = from[i].h[c] - to[k].h[c];
                    cip.cost += tmp * tmp;
                }
                cip.from = i;
                cip.to = k;
                mCostVector.push_back(cip);
            }
        }

        // sort
        std::sort(mCostVector.begin(), mCostVector.end(), cmpfunc);

        // extract
        std::vector<bool> used1(from.size(), false),
                          used2(to.size(), false);
        const uint limit = (NumericalType)1.0 * mConfig.bfMatchesLimit * std::min(from.size(), to.size());
        matches->reserve(limit);
        for (uint i = 0; i < mCostVector.size(); ++i)
        {
            const uint i1 = mCostVector[i].from,
                       i2 = mCostVector[i].to;
            if (used1[i1] || used2[i2]) continue;
            IRONMatch<NumericalType> match;
            match.from = &(from[i1]);
            match.to   = &(to[i2]);
            match.cost = mCostVector[i].cost;
            matches->push_back(match);
            if (matches->size() == limit) break;
            used1[i1] = true;
            used2[i2] = true;
        }
    }


    // outlier detection
    IRONTransformResult<NumericalType> detectOutliers(typename IRONMatchVector<NumericalType>::type *inlierset,
                                                      const typename IRONMatchVector<NumericalType>::type &matches)
    {
        // best fitting transform
        IRONTransformResult<NumericalType> result;
        inlierset->clear();

        // check
        if (matches.size() < 3u)
        {
            std::cerr << "not enough matches to run RANSAC\n";
            return result;
        }

        // vars
        std::vector<uint> supporterVec, bestSupporterVec;
        supporterVec.reserve(matches.size());
        bestSupporterVec.reserve(matches.size());

        // main ransac loop
        const uint limit = mConfig.ransacLoops;
        const NumericalType tolerance = mConfig.matchingTolerance;
        NumericalType bestscore = (NumericalType)0.0;

        Eigen::Matrix<NumericalType, 3, 1> pt1, pt2, n1, n2;
        Eigen::Transform<NumericalType, 3, Eigen::Affine> tinv, treg, rot1, rot2, TF;
        tinv = treg = rot1 = rot2 = Eigen::Transform<NumericalType, 3, Eigen::Affine>::Identity();
        // break after "limit" loops with no improvement
        for (uint a = 0; a < limit; ++a)
        {
            // draw 3 random samples
            uint idx1 = mRnd.rand() % matches.size(),
                 idx2 = mRnd.rand() % matches.size(),
                 idx3 = mRnd.rand() % matches.size();
            if (idx1 == idx2 || idx2 == idx3 || idx1 == idx3) continue;

            // compute normals for each 3 points
            const Eigen::Matrix<NumericalType, 3, 1> &center1 = matches[idx1].from->mu(),
                                                     &center2 = matches[idx1].to->mu();
            pt1 = matches[idx2].from->mu() - center1;
            pt2 = matches[idx2].to->mu()   - center2;
            n1  = (pt1).cross(matches[idx3].from->mu() - center1);
            n2  = (pt2).cross(matches[idx3].to->mu()   - center2);

            // normalize for angle computation
            pt1.normalize();
            pt2.normalize();
            n1.normalize();
            n2.normalize();

            // construct transform
            tinv.translation() = -center1;
            treg.translation() = center2;
            rot1.linear() = rotFromVectors(n1, n2);
            rot2.linear() = rotFromVectors(rot1 * pt1, pt2);
            TF = treg * rot2 * rot1 * tinv;

            // check supporters
            supporterVec.clear();
            NumericalType tmpscore = (NumericalType)0.0;
            for (uint m = 0; m < matches.size(); ++m)
            {
                // get pair mu
                const Eigen::Matrix<NumericalType, 3, 1> &mu1 = matches[m].from->mu(),
                                                         &mu2 = matches[m].to->mu();

                // distance between descriptors
                const NumericalType dist = (mu2 - TF * mu1).norm();
                if (dist < tolerance)
                {
                    tmpscore += tolerance - dist;
                    supporterVec.push_back(m);
                }
            }

            // best score for now?
            if (tmpscore > bestscore)
            {
                a = 0; // reset loop counter
                std::swap(supporterVec, bestSupporterVec);
                bestscore = tmpscore;
                result.tf = TF;
                result.success = true;
            }
        }

        // now commit to inlierset
        inlierset->reserve(bestSupporterVec.size());
        for (uint i = 0; i < bestSupporterVec.size(); ++i)
        {
            inlierset->push_back(matches[bestSupporterVec[i]]);
        }

        // return best transform
        return result;
    }


    // compute transform from inlierset
    IRONTransformResult<NumericalType> computeTransform(const typename IRONMatchVector<NumericalType>::type &inlierset) const
    {
        IRONTransformResult<NumericalType> result;

        // check size
        if (inlierset.size() < 3u)
        {
            std::cerr << "need at least 3 matches to compute transform\n";
            return result;
        }

        // point means
        Eigen::Matrix<NumericalType, 3, 1> mean1((NumericalType)0.0, (NumericalType)0.0, (NumericalType)0.0),
                                           mean2((NumericalType)0.0, (NumericalType)0.0, (NumericalType)0.0);
        for (uint i = 0; i < inlierset.size(); ++i)
        {
            mean1 += inlierset[i].from->mu();
            mean2 += inlierset[i].to->mu();
        }
        NumericalType invSize = (NumericalType)1.0 / inlierset.size();
        mean1 *= invSize;
        mean2 *= invSize;

        // compute cov
        Eigen::Matrix<NumericalType, 3, 3> cov;
        Eigen::Matrix<NumericalType, 3, 1> da, db;
        cov.setZero();
        for (uint i = 0; i < inlierset.size(); ++i)
        {
            da = inlierset[i].from->mu() - mean1;
            db = inlierset[i].to->mu() - mean2;
            cov += da * db.transpose();
        }

        // SVD
        Eigen::JacobiSVD<Eigen::Matrix<NumericalType, 3, 3> > svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);

        // get rotation matrix
        const NumericalType det = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        Eigen::Matrix<NumericalType, 3, 3> S;
        S << (NumericalType)1.0, (NumericalType)0.0, (NumericalType)0.0,
             (NumericalType)0.0, (NumericalType)1.0, (NumericalType)0.0,
             (NumericalType)0.0, (NumericalType)0.0, det;
        const Eigen::Matrix<NumericalType, 3, 3> R = svd.matrixV() * S * svd.matrixU().transpose();

        // construct transform
        Eigen::Transform<NumericalType, 3, Eigen::Affine> tinv, treg, rot;
        tinv = treg = rot = Eigen::Transform<NumericalType, 3, Eigen::Affine>::Identity();
        tinv.translation() = -mean1;
        treg.translation() = mean2;
        rot.linear() = R;

        // done
        result.tf = treg * rot * tinv;
        result.success = true;
        return result;
    }


private:
    NumericalType limitedAngle(const Eigen::Matrix<NumericalType, 3, 1> &v1,
                               const Eigen::Matrix<NumericalType, 3, 1> &v2) const
    {
        // vectors MUST be normalized
        const NumericalType angle = mACLUT.acos(v1.dot(v2));
        if (angle > (NumericalType)0.5 * M_PI)
        {
            return M_PI - angle;
        }
        return angle;
    }

    // quickly remove descriptor
    void remove(typename IRONDescriptorVector<NumericalType>::type *vec,
                uint idx) const
    {
        // free memory
        delete[] vec->at(idx).h;
        std::swap(vec->at(idx), vec->back());
        vec->pop_back();
    }

    // function for comparison
    static bool cmpfunc(const CostIdxPair<NumericalType> &a,
                        const CostIdxPair<NumericalType> &b)
    {
        return a.cost < b.cost;
    }

    // get 3D rotation matrix from two normalized vectors
    Eigen::Matrix<NumericalType, 3, 3> rotFromVectors(const Eigen::Matrix<NumericalType, 3, 1> &v1,
                                                      const Eigen::Matrix<NumericalType, 3, 1> &v2) const
    {
        // vectors MUST be normalized
        const Eigen::Matrix<NumericalType, 3, 1> v = v1.cross(v2);
        const NumericalType s = v.norm(),
                            c = v1.dot(v2);
        Eigen::Matrix<NumericalType, 3, 3> vx;
        vx << (NumericalType)  0.0, (NumericalType)-v(2), (NumericalType) v(1),
              (NumericalType) v(2), (NumericalType)  0.0, (NumericalType)-v(0),
              (NumericalType)-v(1), (NumericalType) v(0), (NumericalType)  0.0;
        return Eigen::Matrix<NumericalType, 3, 3>::Identity() + vx + vx * vx * ((NumericalType)1.0 - c) / (s * s);
    }

protected:
    IRONConfig mConfig;
    ACosLUT<NumericalType> mACLUT;
    std::vector<CostIdxPair<NumericalType> > mCostVector;
    XorShift mRnd;
};


typedef IRONEngine<float> IRONEnginef;
typedef IRONEngine<double> IRONEngined;


// namespace
}


#endif
