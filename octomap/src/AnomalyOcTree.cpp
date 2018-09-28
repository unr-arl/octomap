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

#include <octomap/AnomalyOcTree.h>

namespace octomap {

  // node implementation  --------------------------------------
  std::ostream& AnomalyOcTreeNode::writeData(std::ostream &s) const {
    s.write((const char*) &(anomaly.value), sizeof(anomaly.value)); // anomaly value
    s.write((const char*) &(anomaly.type), sizeof(anomaly.type)); // anomaly type
    return s;
  }

  std::istream& AnomalyOcTreeNode::readData(std::istream &s) {
    s.read((char*) &(anomaly.value), sizeof(anomaly.value)); // anomaly value
    s.read((char*) &(anomaly.type), sizeof(anomaly.type)); // anomaly type
    return s;
  }

  AnomalyOcTreeNode::Anomaly AnomalyOcTreeNode::getAverageChildAnomaly() const {
    int m = 0;
    int c = 0;

    if (children != NULL){
      for (int i=0; i<8; i++) {
        AnomalyOcTreeNode* child = static_cast<AnomalyOcTreeNode*>(children[i]);

        if (child != NULL && child->isAnomalySet()) {
          m += child->getAnomaly().value;
          ++c;
        }
      }
    }

    if (c > 0) {
      m /= c;
      return Anomaly(m);
    }
    else { // no child had a anomaly other than white
      return Anomaly();
    }
  }

  void AnomalyOcTreeNode::updateAnomalyChildren() {
    anomaly = getAverageChildAnomaly();
  }

  // tree implementation  --------------------------------------
  AnomalyOcTree::AnomalyOcTree(double resolution)
  : OccupancyOcTreeBase<AnomalyOcTreeNode>(resolution) {
    anomalyOcTreeMemberInit.ensureLinking();
  };

  AnomalyOcTreeNode* AnomalyOcTree::setNodeAnomaly(const OcTreeKey& key,
                                             uint8_t sal) {
    AnomalyOcTreeNode* n = search (key);
    if (n != 0) {
      n->setAnomaly(sal);
    }
    return n;
  }

  bool AnomalyOcTree::pruneNode(AnomalyOcTreeNode* node) {
    if (!isNodeCollapsible(node))
      return false;

    // set value to children's values (all assumed equal)
    node->copyData(*(getNodeChild(node, 0)));

    if (node->isAnomalySet()) // TODO check
      node->setAnomaly(node->getAverageChildAnomaly());

    // delete children
    for (unsigned int i=0;i<8;i++) {
      deleteNodeChild(node, i);
    }
    delete[] node->children;
    node->children = NULL;

    return true;
  }

  bool AnomalyOcTree::isNodeCollapsible(const AnomalyOcTreeNode* node) const{

    // TODO tung: disable prune node now
    return false;

    // all children must exist, must not have children of
    // their own and have the same occupancy probability
    if (!nodeChildExists(node, 0))
      return false;

    const AnomalyOcTreeNode* firstChild = getNodeChild(node, 0);
    if (nodeHasChildren(firstChild))
      return false;

    for (unsigned int i = 1; i<8; i++) {
      // compare nodes only using their occupancy, ignoring anomaly for pruning
      if (!nodeChildExists(node, i) || nodeHasChildren(getNodeChild(node, i)) || !(getNodeChild(node, i)->getValue() == firstChild->getValue()))
        return false;
    }

    return true;
  }

  AnomalyOcTreeNode* AnomalyOcTree::averageNodeAnomaly(const OcTreeKey& key,
                                                 uint8_t sal) {
    // TODO tung: possible derive the average here
    AnomalyOcTreeNode* n = search(key);
    if (n != 0) {
      if (n->isAnomalySet()) {
        AnomalyOcTreeNode::Anomaly prev_anomaly = n->getAnomaly();
        n->setAnomaly((prev_anomaly.value + sal)/2);
      }
      else {
        n->setAnomaly(sal);
      }
    }
    return n;
  }

  AnomalyOcTreeNode* AnomalyOcTree::integrateNodeAnomaly(const OcTreeKey& key,
                                                   uint8_t sal) {
    AnomalyOcTreeNode* n = search (key);
    if (n != 0) {
      if (n->isAnomalySet()) {
        AnomalyOcTreeNode::Anomaly prev_anomaly = n->getAnomaly();
        double node_prob = n->getOccupancy();
        uint8_t new_sal = (uint8_t) ((double) prev_anomaly.value * node_prob
                                               +  (double) sal * (0.99-node_prob));
        n->setAnomaly(new_sal);
      }
      else {
        n->setAnomaly(sal);
      }
    }
    return n;
  }


  void AnomalyOcTree::updateInnerOccupancy() {
    this->updateInnerOccupancyRecurs(this->root, 0);
  }

  void AnomalyOcTree::updateInnerOccupancyRecurs(AnomalyOcTreeNode* node, unsigned int depth) {

    //TODO tung check this function updateAnomalyChildren
    // only recurse and update for inner nodes:
    if (nodeHasChildren(node)){
      // return early for last level:
      if (depth < this->tree_depth){
        for (unsigned int i=0; i<8; i++) {
          if (nodeChildExists(node, i)) {
            updateInnerOccupancyRecurs(getNodeChild(node, i), depth+1);
          }
        }
      }
      node->updateOccupancyChildren();
      node->updateAnomalyChildren();
    }
  }

  void AnomalyOcTree::writeAnomalyHistogram(std::string filename) {
  }

  std::ostream& operator<<(std::ostream& out, AnomalyOcTreeNode::Anomaly const& s) {
    return out << '(' << (uint8_t)s.value << ' ' << (uint8_t)s.type << ' ' << (int)s.timestamp << ')';
  }

  AnomalyOcTree::StaticMemberInitializer AnomalyOcTree::anomalyOcTreeMemberInit;

} // end namespace
