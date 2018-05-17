/*
 * OctoMap - An Efficient Probabilistic 3D Mapping Framework Based on Octrees
 * http://octomap.github.com/
 *
 * Copyright (c) 2009-2013, K.M. Wurm and A. Hornung, University of Freiburg
 * All rights reserved.
 * License: New BSD
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef OCTOMAP_SALIENCY_OCTREE_H
#define OCTOMAP_SALIENCY_OCTREE_H


#include <iostream>
#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>

namespace octomap {

  // forward declaraton for "friend"
  class SaliencyOcTree;

  // Node of saliency voxel
  class SaliencyOcTreeNode : public OcTreeNode {
  public:
    friend class SaliencyOcTree; // needs access to node children (inherited)

    // Define paramaters stored in saliency voxel
    class Saliency {
    public:
      enum Type {
        VOXEL_NORMAL    = (uint8_t)0,
        VOXEL_SALIENCY  = (uint8_t)1,
        VOXEL_RETIRED   = (uint8_t)2
      };

      Saliency(): value(0),
                  value_buff(0),
                  type(VOXEL_NORMAL),
                  counter(0),
                  timestamp(0),
                  viewpoint(0),
                  density(0){}
      Saliency(uint8_t sal): value(0),
                  value_buff(0),
                  type(VOXEL_NORMAL),
                  counter(0),
                  timestamp(0),
                  viewpoint(0),
                  density(0){
        value = sal;
      }
      inline bool operator== (const Saliency &other) const {
        return (value==other.value && type==other.type);
      }

      inline bool operator!= (const Saliency &other) const {
        return (value!=other.value || type!=other.type);
      }

      uint8_t value;
      uint8_t value_buff;
      Type type;
      uint16_t counter;
      uint32_t timestamp;
      uint32_t viewpoint;
      uint32_t density;
    };

  public:
    SaliencyOcTreeNode() : OcTreeNode() {}

    SaliencyOcTreeNode(const SaliencyOcTreeNode& rhs) : OcTreeNode(rhs), saliency(rhs.saliency) {}

    bool operator==(const SaliencyOcTreeNode& rhs) const{
      return (rhs.value == value && rhs.saliency == saliency);
    }

    void copyData(const SaliencyOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->saliency =  from.getSaliency();
    }

    inline Saliency getSaliency() const { return saliency; }
    inline void  setSaliency(Saliency s) {this->saliency = s; }
    inline void  setSaliency(uint8_t sal) {this->saliency.value = sal; }

    Saliency& getSaliency() { return saliency; }

    inline bool isSaliencySet() const {
      return (saliency.timestamp);
    }

    void updateSaliencyChildren();

    SaliencyOcTreeNode::Saliency getAverageChildSaliency() const;

    // file I/O
    std::istream& readData(std::istream &s);
    std::ostream& writeData(std::ostream &s) const;

  protected:
    Saliency saliency;
  };

  // tree definition
  class SaliencyOcTree : public OccupancyOcTreeBase <SaliencyOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    SaliencyOcTree(double resolution);

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    SaliencyOcTree* create() const {return new SaliencyOcTree(resolution); }

    std::string getTreeType() const {return "SaliencyOcTree";}

     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(SaliencyOcTreeNode* node);

    virtual bool isNodeCollapsible(const SaliencyOcTreeNode* node) const;

    // set node color at given key or coordinate. Replaces previous color.
    SaliencyOcTreeNode* setNodeSaliency(const OcTreeKey& key, uint8_t sal);

    SaliencyOcTreeNode* setNodeSaliency(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeSaliency(key, sal);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    SaliencyOcTreeNode* averageNodeSaliency(const OcTreeKey& key, uint8_t sal);

    SaliencyOcTreeNode* averageNodeSaliency(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeSaliency(key, sal);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    SaliencyOcTreeNode* integrateNodeSaliency(const OcTreeKey& key, uint8_t sal);

    SaliencyOcTreeNode* integrateNodeSaliency(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeSaliency(key, sal);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeSaliencyHistogram(std::string filename);

  protected:
    void updateInnerOccupancyRecurs(SaliencyOcTreeNode* node, unsigned int depth);

    /**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           SaliencyOcTree* tree = new SaliencyOcTree(0.1);
           tree->clearKeyRays();
           AbstractOcTree::registerTreeType(tree);
         }

         /**
         * Dummy function to ensure that MSVC does not drop the
         * StaticMemberInitializer, causing this tree failing to register.
         * Needs to be called from the constructor of this octree.
         */
         void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer saliencyOcTreeMemberInit;

  };

  //! user friendly output in format (r g b)
  std::ostream& operator<<(std::ostream& out, SaliencyOcTreeNode::Saliency const& s);

} // end namespace

#endif
