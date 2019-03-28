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


#ifndef _KNNADAPTER_H_
#define _KNNADAPTER_H_


////////////////////////////////////////
// INCLUDES
////////////////////////////////////////
#include "_nanoflann.hpp"
#include <Eigen/Dense>


////////////////////////////////////////
// NAMESPACE
////////////////////////////////////////
namespace IRON
{
////////////////////////////////////////
// GenericKNNAdapter
////////////////////////////////////////
template <int VecDim,
          typename Container,
          typename Accessor,
          typename NumericalType>
class GenericKNNAdapter;


template <int VecDim,
          typename Container,
          typename Accessor,
          typename NumericalType>
struct KDTree
{
    typedef nanoflann::KDTreeSingleIndexAdaptor
            <nanoflann::L2_Adaptor<NumericalType, GenericKNNAdapter<VecDim, Container, Accessor, NumericalType> >,
            GenericKNNAdapter<VecDim, Container, Accessor, NumericalType>, VecDim, uint> type;
};


///////////////////////////////////////////////////
// VecDim...: number of dimensions of an item
//            -1 for dynamic allocation at run-time
//            (give attach() the correct dim)
// Container: e.g. std::vector<item>
//            must implement size(), operator[]
// Accessor.: must have a static function
//            at(item, dim) to get n-th element
// NumericalType: double / float / ...
///////////////////////////////////////////////////
template <int VecDim,
          typename Container,
          typename Accessor,
          typename NumericalType>
class GenericKNNAdapter
{
private:
    GenericKNNAdapter(const GenericKNNAdapter &other);
    GenericKNNAdapter &operator=(const GenericKNNAdapter &other);

public:
    // containers for search results (index, distance)
    typedef std::pair<uint, NumericalType> IdxDistPair;
    typedef std::vector<IdxDistPair> IdxDistPairs;

    // ctor, dtor
    GenericKNNAdapter() : mList(NULL), mTree(NULL), mDim(-1) { }
    ~GenericKNNAdapter()
    {
        if (mTree != NULL) delete mTree;
    }

    // attach a container, note that the
    // actual data is not stored here
    // if VecDim is -1, you MUST specify "allocDim" here
    void attach(const Container *list, int allocDim = VecDim)
    {
        mDim = allocDim;
        mList = list;
        // rebuild tree
        rebuild();
    }

    // refresh tree, when data has changed
    void rebuild()
    {
        if (mTree != NULL) delete mTree;
        mTree = new typename KDTree<VecDim,
                                    Container,
                                    Accessor,
                                    NumericalType>::type(mDim,
                                                         *this,
                                                         nanoflann::KDTreeSingleIndexAdaptorParams(10));
        mTree->buildIndex();
    }

    // has tree?
    bool hasTree() const
    {
        return mTree != NULL;
    }

    // get neighbors within squared radius
    const IdxDistPairs &radiusSearch(const typename Container::value_type &item,
                                     NumericalType radiusSqr)
    {
        mPairs.clear();
        nanoflann::SearchParams params;
        params.sorted = false;
        NumericalType *pt = new NumericalType[mDim];
        for (int i = 0; i < mDim; ++i)
        {
            pt[i] = Accessor::at(item, i);
        }
        mTree->radiusSearch(pt, radiusSqr, mPairs, params);
        delete[] pt;
        return mPairs;
    }

    // ...stripped down for publication...

    // get direct access to tree
    const typename KDTree<VecDim, Container, Accessor, NumericalType>::type *getTree()
    {
        return mTree;
    }

    // now some nanoflann stuff
    uint kdtree_get_point_count() const { return mList->size(); }
    NumericalType kdtree_get_pt(uint idx, int dim) const
    {
        return Accessor::at((*mList)[idx], dim);
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const { return false; }

private:
    // pointer to data
    const Container *mList;

    // actual k-d-tree
    typename KDTree<VecDim,
                    Container,
                    Accessor,
                    NumericalType>::type *mTree;
    int mDim;
    IdxDistPairs mPairs;
};


////////////////////////////////////////////
// The same with eigen vector
////////////////////////////////////////////
typedef std::vector<Eigen::Vector3f, Eigen::aligned_allocator<Eigen::Vector3f> > PointVector3f;
typedef std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > PointVector3d;

template <typename NumericalType>
struct PointAccessor
{
private:
    PointAccessor();
    PointAccessor(const PointAccessor &other);
    PointAccessor &operator=(const PointAccessor &other);

public:
    static NumericalType at(const Eigen::Matrix<NumericalType, 3, 1> &item, uint i) { return item[i]; }
};

typedef GenericKNNAdapter<3, PointVector3f, PointAccessor<float>, float> PointVectorKNNAdapter3f;
typedef GenericKNNAdapter<3, PointVector3d, PointAccessor<double>, double> PointVectorKNNAdapter3d;


// namespace
}

#endif
