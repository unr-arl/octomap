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

#ifndef OCTOMAP_ANOMALY_OCTREE_H
#define OCTOMAP_ANOMALY_OCTREE_H


#include <iostream>
#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>

#define R_BINS    3
#define YAW_BINS  16

namespace octomap {

  // forward declaraton for "friend"
  class AnomalyOcTree;

  // Node of anomaly voxel
  class AnomalyOcTreeNode : public OcTreeNode {
  public:
    friend class AnomalyOcTree; // needs access to node children (inherited)

    // Define paramaters stored in anomaly voxel
    class Anomaly {
    public:
      enum Type {
        VOXEL_NORMAL    = (uint8_t)0,
        VOXEL_ANOMALY  = (uint8_t)1,
        VOXEL_RETIRED   = (uint8_t)2
      };

      Anomaly(): value(0),
                  type(VOXEL_NORMAL),
                  logodds(0),
                  timestamp(0){
        for (int i=0; i < R_BINS; i++)
          for (int j=0; j < YAW_BINS; j++)
            hist[i][j] = 0;
      }
      Anomaly(uint8_t sal): value(0),
                  type(VOXEL_NORMAL),
                  logodds(0),
                  timestamp(0){
        value = sal;
        for (int i=0; i < R_BINS; i++)
          for (int j=0; j < YAW_BINS; j++)
            hist[i][j] = 0;
      }
      inline bool operator== (const Anomaly &other) const {
        return (value==other.value && type==other.type);
      }

      inline bool operator!= (const Anomaly &other) const {
        return (value!=other.value || type!=other.type);
      }

      uint8_t value; // value to scale the
      Type type;
      float logodds;
      uint32_t timestamp;
      uint8_t hist[R_BINS][YAW_BINS];
    };

  public:
    AnomalyOcTreeNode() : OcTreeNode() {}

    AnomalyOcTreeNode(const AnomalyOcTreeNode& rhs) : OcTreeNode(rhs), anomaly(rhs.anomaly) {}

    bool operator==(const AnomalyOcTreeNode& rhs) const{
      return (rhs.value == value && rhs.anomaly == anomaly);
    }

    void copyData(const AnomalyOcTreeNode& from){
      OcTreeNode::copyData(from);
      this->anomaly =  from.getAnomaly();
    }

    inline void getHistBins(int *bins) {
      bins[0] = R_BINS;
      bins[1] = YAW_BINS;
    }

    inline Anomaly getAnomaly() const { return anomaly; }
    inline void  setAnomaly(Anomaly s) {this->anomaly = s; }
    inline void  setAnomaly(uint8_t sal) {this->anomaly.value = sal; }

    Anomaly& getAnomaly() { return anomaly; }

    inline bool isAnomalySet() const {
      return (anomaly.timestamp);
    }

    void updateAnomalyChildren();

    AnomalyOcTreeNode::Anomaly getAverageChildAnomaly() const;

    // file I/O
    std::istream& readData(std::istream &s);
    std::ostream& writeData(std::ostream &s) const;

  protected:
    Anomaly anomaly;
  };

  // tree definition
  class AnomalyOcTree : public OccupancyOcTreeBase <AnomalyOcTreeNode> {

  public:
    /// Default constructor, sets resolution of leafs
    AnomalyOcTree(double resolution);

    /// virtual constructor: creates a new object of same type
    /// (Covariant return type requires an up-to-date compiler)
    AnomalyOcTree* create() const {return new AnomalyOcTree(resolution); }

    std::string getTreeType() const {return "AnomalyOcTree";}

     /**
     * Prunes a node when it is collapsible. This overloaded
     * version only considers the node occupancy for pruning,
     * different colors of child nodes are ignored.
     * @return true if pruning was successful
     */
    virtual bool pruneNode(AnomalyOcTreeNode* node);

    virtual bool isNodeCollapsible(const AnomalyOcTreeNode* node) const;

    // set node color at given key or coordinate. Replaces previous color.
    AnomalyOcTreeNode* setNodeAnomaly(const OcTreeKey& key, uint8_t sal);

    AnomalyOcTreeNode* setNodeAnomaly(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return setNodeAnomaly(key, sal);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    AnomalyOcTreeNode* averageNodeAnomaly(const OcTreeKey& key, uint8_t sal);

    AnomalyOcTreeNode* averageNodeAnomaly(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return averageNodeAnomaly(key, sal);
    }

    // integrate color measurement at given key or coordinate. Average with previous color
    AnomalyOcTreeNode* integrateNodeAnomaly(const OcTreeKey& key, uint8_t sal);

    AnomalyOcTreeNode* integrateNodeAnomaly(float x, float y, float z, uint8_t sal) {
      OcTreeKey key;
      if (!this->coordToKeyChecked(point3d(x,y,z), key)) return NULL;
      return integrateNodeAnomaly(key, sal);
    }

    // update inner nodes, sets color to average child color
    void updateInnerOccupancy();

    // uses gnuplot to plot a RGB histogram in EPS format
    void writeAnomalyHistogram(std::string filename);

  protected:
    void updateInnerOccupancyRecurs(AnomalyOcTreeNode* node, unsigned int depth);

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
           AnomalyOcTree* tree = new AnomalyOcTree(0.1);
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
    static StaticMemberInitializer anomalyOcTreeMemberInit;

  };

  //! user friendly output in format (r g b)
  std::ostream& operator<<(std::ostream& out, AnomalyOcTreeNode::Anomaly const& s);

} // end namespace

#endif
