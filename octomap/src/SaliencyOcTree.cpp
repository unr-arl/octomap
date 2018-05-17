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

#include <octomap/SaliencyOcTree.h>

namespace octomap {

  // node implementation  --------------------------------------
  std::ostream& SaliencyOcTreeNode::writeData(std::ostream &s) const {
    s.write((const char*) &(saliency.value), sizeof(saliency.value)); // saliency value
    s.write((const char*) &(saliency.type), sizeof(saliency.type)); // saliency type
    s.write((const char*) &(saliency.viewpoint), sizeof(saliency.viewpoint)); //
    return s;
  }

  std::istream& SaliencyOcTreeNode::readData(std::istream &s) {
    s.read((char*) &(saliency.value), sizeof(saliency.value)); // saliency value
    s.read((char*) &(saliency.type), sizeof(saliency.type)); // saliency type
    s.read((char*) &(saliency.viewpoint), sizeof(saliency.viewpoint)); //
    return s;
  }

  SaliencyOcTreeNode::Saliency SaliencyOcTreeNode::getAverageChildSaliency() const {
    int m = 0;
    int c = 0;

    if (children != NULL){
      for (int i=0; i<8; i++) {
        SaliencyOcTreeNode* child = static_cast<SaliencyOcTreeNode*>(children[i]);

        if (child != NULL && child->isSaliencySet()) {
          m += child->getSaliency().value;
          ++c;
        }
      }
    }

    if (c > 0) {
      m /= c;
      return Saliency(m);
    }
    else { // no child had a saliency other than white
      return Saliency();
    }
  }

  void SaliencyOcTreeNode::updateSaliencyChildren() {
    saliency = getAverageChildSaliency();
  }

  // tree implementation  --------------------------------------
  SaliencyOcTree::SaliencyOcTree(double resolution)
  : OccupancyOcTreeBase<SaliencyOcTreeNode>(resolution) {
    saliencyOcTreeMemberInit.ensureLinking();
  };

  SaliencyOcTreeNode* SaliencyOcTree::setNodeSaliency(const OcTreeKey& key,
                                             uint8_t sal) {
    SaliencyOcTreeNode* n = search (key);
    if (n != 0) {
      n->setSaliency(sal);
    }
    return n;
  }

  bool SaliencyOcTree::pruneNode(SaliencyOcTreeNode* node) {
    if (!isNodeCollapsible(node))
      return false;

    // set value to children's values (all assumed equal)
    node->copyData(*(getNodeChild(node, 0)));

    if (node->isSaliencySet()) // TODO check
      node->setSaliency(node->getAverageChildSaliency());

    // delete children
    for (unsigned int i=0;i<8;i++) {
      deleteNodeChild(node, i);
    }
    delete[] node->children;
    node->children = NULL;

    return true;
  }

  bool SaliencyOcTree::isNodeCollapsible(const SaliencyOcTreeNode* node) const{

    // @tung: Disable prune node now
    return false;

    // all children must exist, must not have children of
    // their own and have the same occupancy probability
    if (!nodeChildExists(node, 0))
      return false;

    const SaliencyOcTreeNode* firstChild = getNodeChild(node, 0);
    if (nodeHasChildren(firstChild))
      return false;

    for (unsigned int i = 1; i<8; i++) {
      // compare nodes only using their occupancy, ignoring saliency for pruning
      if (!nodeChildExists(node, i) || nodeHasChildren(getNodeChild(node, i)) || !(getNodeChild(node, i)->getValue() == firstChild->getValue()))
        return false;
    }

    return true;
  }

  SaliencyOcTreeNode* SaliencyOcTree::averageNodeSaliency(const OcTreeKey& key,
                                                 uint8_t sal) {
    SaliencyOcTreeNode* n = search(key);
    if (n != 0) {
      if (n->isSaliencySet()) {
        SaliencyOcTreeNode::Saliency prev_saliency = n->getSaliency();
        n->setSaliency((prev_saliency.value + sal)/2);
      }
      else {
        n->setSaliency(sal);
      }
    }
    return n;
  }

  SaliencyOcTreeNode* SaliencyOcTree::integrateNodeSaliency(const OcTreeKey& key,
                                                   uint8_t sal) {
    SaliencyOcTreeNode* n = search (key);
    if (n != 0) {
      if (n->isSaliencySet()) {
        SaliencyOcTreeNode::Saliency prev_saliency = n->getSaliency();
        double node_prob = n->getOccupancy();
        uint8_t new_sal = (uint8_t) ((double) prev_saliency.value * node_prob
                                               +  (double) sal * (0.99-node_prob));
        n->setSaliency(new_sal);
      }
      else {
        n->setSaliency(sal);
      }
    }
    return n;
  }


  void SaliencyOcTree::updateInnerOccupancy() {
    this->updateInnerOccupancyRecurs(this->root, 0);
  }

  void SaliencyOcTree::updateInnerOccupancyRecurs(SaliencyOcTreeNode* node, unsigned int depth) {

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
      node->updateSaliencyChildren();
    }
  }

  void SaliencyOcTree::writeSaliencyHistogram(std::string filename) {

#ifdef _MSC_VER
    fprintf(stderr, "The saliency histogram uses gnuplot, this is not supported under windows.\n");
#else
    // build RGB histogram
    std::vector<int> histogram_saliency (256,0);
    for(SaliencyOcTree::tree_iterator it = this->begin_tree(),
          end=this->end_tree(); it!= end; ++it) {
      if (!it.isLeaf() || !this->isNodeOccupied(*it)) continue;
      SaliencyOcTreeNode::Saliency& s = it->getSaliency();
      ++histogram_saliency[s.value];
    }
    // plot data
    FILE *gui = popen("gnuplot ", "w");
    fprintf(gui, "set term postscript eps enhanced saliency\n");
    fprintf(gui, "set output \"%s\"\n", filename.c_str());
    fprintf(gui, "plot [-1:256] ");
    fprintf(gui,"'-' w filledcurve lt 1 lc 1 tit \"r\",");
    fprintf(gui, "'-' w filledcurve lt 1 lc 2 tit \"g\",");
    fprintf(gui, "'-' w filledcurve lt 1 lc 3 tit \"b\",");
    fprintf(gui, "'-' w l lt 1 lc 1 tit \"\",");
    fprintf(gui, "'-' w l lt 1 lc 2 tit \"\",");
    fprintf(gui, "'-' w l lt 1 lc 3 tit \"\"\n");

    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui,"0 0\n"); fprintf(gui, "e\n");
    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui, "e\n");
    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui, "e\n");
    for (int i=0; i<256; ++i) fprintf(gui,"%d %d\n", i, histogram_saliency[i]);
    fprintf(gui, "e\n");
    fflush(gui);
#endif
  }

  std::ostream& operator<<(std::ostream& out, SaliencyOcTreeNode::Saliency const& s) {
    return out << '(' << (uint8_t)s.value << ' ' << (uint8_t)s.type << ' ' << (int)s.timestamp << ')';
  }

  SaliencyOcTree::StaticMemberInitializer SaliencyOcTree::saliencyOcTreeMemberInit;

} // end namespace
