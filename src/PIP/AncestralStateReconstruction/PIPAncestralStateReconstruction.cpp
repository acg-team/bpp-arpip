/*******************************************************************************
 * Licensed Materials - Property of Gholamhossein Jowkar
 *
 *
 * Copyright (C) 2019-2023 by Gholamhossein Jowkar
 *
 * This work was supported by the Swiss National Science Foundation (SNF) grants $31003A\_176316$ to Dr. M. Anisimova.
 * The funding body did not play any role in the design of the study and collection, analysis, and interpretation of
 * data and in writing the code.
 *******************************************************************************
 *
 * This file is part of ARPIP project
 *
 * ARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process
 * ARPIP is a joint maximum likelihood approach for phylogenetic ancestral sequence reconstruction, capable of modeling
 * indels both biological and mathematically.
 *
 *
 * This software is based and extends the following libraries:
 *
 * - the Bio++ libraries
 *   developed by the Bio++ Development Team <http://biopp.univ-montp2.fr>
 *
 *
 * ARPIP is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 *
 * ARPIP is a free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with ARPIP. If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/

/**
 * @file PIPAncestralStateReconstruction.cpp
 * @author Gholamhossein Jowkar
 * @date 11 08 2020
 * @version 1.0.0
 * @maintainer Gholamhossein Jowkar
 * @email jowk@zhaw.ch
 * @status Development
 *
 * @brief
 * @details
 * @pre
 * @bug
 * @warning
 *
 * @see For more information visit:
 */
// From ARPIP:
#include "PIPAncestralStateReconstruction.h"

// Form bpp:
#include<Bpp/Numeric/VectorTools.h>
#include<Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

/************************************************* Constructor ********************************************************/
PIPAncestralStateReconstruction::PIPAncestralStateReconstruction(const PIPDRTreeLikelihood *lik,
                                                                 const PIPMLIndelPoints *mlIndelPoints,
                                                                 const std::string &probProfileType) :
        likelihood_(lik),
        mlIndelPoints_(mlIndelPoints),
        tree_(lik->getTree()),
        alphabet_(lik->getAlphabet()),
        nbInnerNodes_(lik->getTree().getInnerNodesId().size()),
        nbSites_(lik->getLikelihoodData()->getNumberOfSites()),
        nbDistinctSites_(lik->getLikelihoodData()->getNumberOfDistinctSites()),
        nbStates_(lik->getLikelihoodData()->getNumberOfStates()),
        rootPatternLinks_(lik->getLikelihoodData()->getRootArrayPositions()),
        probabilityProfile_(0),
        probProfileType_(probProfileType) {

    DLOG(INFO) << "[PIP ASR] PIPAncestralStateReconstruction object is constructed successfully.";
}

/********************************************* Copy constructor *******************************************************/
PIPAncestralStateReconstruction::PIPAncestralStateReconstruction(const PIPAncestralStateReconstruction &arpip):
        likelihood_(arpip.likelihood_),
        mlIndelPoints_(arpip.mlIndelPoints_),
        tree_(arpip.tree_),
        alphabet_(arpip.alphabet_),
        nbInnerNodes_(arpip.nbInnerNodes_),
        nbSites_(arpip.nbSites_),
        nbDistinctSites_(arpip.nbDistinctSites_),
        nbStates_(arpip.nbStates_),
        rootPatternLinks_(arpip.rootPatternLinks_),
        probabilityProfile_(arpip.probabilityProfile_),
        probProfileType_(arpip.probProfileType_){}

/********************************************* Operator overloading *******************************************************/

PIPAncestralStateReconstruction &PIPAncestralStateReconstruction::operator=(const PIPAncestralStateReconstruction &arpip){
    likelihood_         = arpip.likelihood_;
    mlIndelPoints_      = arpip.mlIndelPoints_;
    tree_               = arpip.tree_;
    alphabet_           = arpip.alphabet_;
    nbInnerNodes_       = arpip.nbInnerNodes_;
    nbSites_            = arpip.nbSites_;
    nbDistinctSites_    = arpip.nbDistinctSites_;
    nbStates_           = arpip.nbStates_;
    rootPatternLinks_   = arpip.rootPatternLinks_;
    probabilityProfile_= arpip.probabilityProfile_;
    probProfileType_    = arpip.probProfileType_;

    return *this;
}
/****************************************** PIP based functions **********************************************/
/**********************************************************************************************************************/

void PIPAncestralStateReconstruction::getAncestralStatesForSite(const Node *node,
                                                                std::map<int, std::vector<int> > &ancestors,
                                                                VVdouble &lArray, VVdouble &cArray,
                                                                const size_t siteNumber, const int newRootId) const {

    int nodeId = node->getId();
    if (nodeId == newRootId) {
//        std::cout << bpp::VectorTools::whichMax(lArray[nodeId]*10000);
//        std::cout << std::max_element(lArray[nodeId].begin(), lArray[nodeId].end());
        ancestors[nodeId][siteNumber] = bpp::VectorTools::whichMax(lArray[nodeId]);
    } else {
        const Node *father = node->getFather();
        int winnerIndex = ancestors[father->getId()][siteNumber];
        ancestors[nodeId][siteNumber] = cArray[nodeId][winnerIndex];
    }

}
/**********************************************************************************************************************/

std::tuple<std::map<int, std::vector<int>>, VVVdouble> PIPAncestralStateReconstruction::getAllAncestralStatesWGap() const{

    std::map<int, std::vector<int>> ancestors;

    // Pupko's Arrays:
    VVVdouble likelihoodArray;
    VVVdouble cLikelihoodArray;

    // Allocating the mentioned arrays:
    likelihoodArray.resize(nbDistinctSites_);
    cLikelihoodArray.resize(nbDistinctSites_);

    // Clone the data into a AlignedSequenceContainer for more efficiency:
     AlignedSequenceContainer *data = new AlignedSequenceContainer(*likelihood_->getLikelihoodData()->getShrunkData());
    DVLOG(2) << "[PIP ASR] Cloning the data with "<< data->getNumberOfSites()<<" sites.";

    // Indel points extracted using MLIndelPoints algorithm:
    std::map<size_t, std::vector<std::string>> insertionPoints = mlIndelPoints_->getInsertionPoints();
    std::map<size_t, std::vector<std::string>> deletionPoints = mlIndelPoints_->getDeletionPoints();
    DLOG(INFO) << "[PIP ASR] The insertion and deletion points copied from IndelPoints object containing " <<
             insertionPoints.size() << " insertion points and " << deletionPoints.size() << " deletion points.";

    // List of internal nodes:
    std::vector<int> ancestorNodeListNb = tree_.getInnerNodesId();

    for (size_t siteNb = 0; siteNb < nbDistinctSites_; siteNb++) {
        // Copying the existing tree:
        TreeTemplate<Node> *PIPSubTree = tree_.clone();
        if ((insertionPoints[siteNb].at(0).empty()))
            throw Exception(
                    "PIPAncestralStateReconstruction::getAllAncestralStatesWGap. Error, NO Insertion point for this site!");
        // Extracting the subtree based on insertion point:
        const Node *newRoot = PIPSubTree->getNode(std::stoi(insertionPoints[siteNb][0]), 0);
        DVLOG(2) << "[PIP ASR] The new root for site number " << siteNb << " is node number: " << newRoot->getId();

        // case newroot != treeRoot
        if(PIPSubTree->getRootId() != newRoot->getId()){
            PIPSubTree = PIPSubTree->cloneSubtree(newRoot->getId());
            DVLOG(2)<< "[PIP ASR] The new root is different from the tree root.";
        }

        // Initialize the ancestor array with gap for both scenario that there are deletions and an insertion different than the root.
        for(auto innerNodeNb:ancestorNodeListNb){
            ancestors[innerNodeNb].resize(nbDistinctSites_);
            ancestors[innerNodeNb][siteNb]= -1;//char for gap state
            DVLOG(2) << "[PIP ASR] The ancestor value for internal node nb " << innerNodeNb << " and siteNb " << siteNb
                     << " is initialized by:" << ancestors[innerNodeNb][siteNb];
        }
        std::vector<int> PIPDeletionPoints;
        if (!newRoot->isLeaf()) {
            for(auto &str:deletionPoints[siteNb]) {
                PIPDeletionPoints.push_back(std::stoi(str));
                if(PIPSubTree->isLeaf(std::stoi(str))){ //todo: check!
                    Node *fNode = PIPSubTree->getNode(std::stoi(str))->getFather();
                    Node *cNode = PIPSubTree->getNode(std::stoi(str));
                    fNode->removeSon(cNode);
                } else {

                    Node *fNode = PIPSubTree->getNode(std::stoi(str))->getFather();
                    Node *cNode = PIPSubTree->getNode(std::stoi(str));
                    fNode->removeSon(cNode);
                }
            }
            Node *newPIPRoot = PIPSubTree->getRootNode();
            recursiveJointAncestralStates(newPIPRoot, PIPSubTree, ancestors, likelihoodArray, cLikelihoodArray,
                                          *data, siteNb, newRoot->getId(), 0);
        }
        delete newRoot;
    }
    delete data;
    return std::make_tuple(ancestors, likelihoodArray);
}
/**********************************************************************************************************************/

void PIPAncestralStateReconstruction::recursiveJointAncestralStates(const Node *node, const TreeTemplate<Node> *tree,
                                                                    std::map<int, std::vector<int> > &ancestors,
                                                                    VVVdouble &likelihoodArray, VVVdouble &cLikelihoodArray,
                                                                    AlignedSequenceContainer &data,
                                                                    const size_t siteNumber,
                                                                    const int newRootId, bool isVisited) const {

    VVdouble *larray = &likelihoodArray[siteNumber];
    VVdouble *carray = &cLikelihoodArray[siteNumber];

    bool *visited = &isVisited;
    DVLOG(2) << "[PIP ASR] This node " << node->getId() << " has been visited: " << isVisited;

    if (!(*visited)) {
        *visited = true;
        recursiveComputeLikelihoodAtSite(node, tree, *larray,*carray, siteNumber, newRootId);
    }
    if (!node->isLeaf() || tree->getRootId()==node->getId()) {//check if the new subtree is not a single leaf node or is a rooted with one son.
        getAncestralStatesForSite(node, ancestors, *larray, *carray, siteNumber, newRootId);
        for (size_t i = 0; i < node->getNumberOfSons(); i++) {
            recursiveJointAncestralStates(node->getSon(i), tree, ancestors, likelihoodArray, cLikelihoodArray, data, siteNumber,
                                          newRootId, isVisited);
        }
    }
}
/**********************************************************************************************************************/

void PIPAncestralStateReconstruction::recursiveComputeLikelihoodAtSite(const Node *node,
                                                                       const TreeTemplate<Node> *tree,
                                                                       VVdouble &likelihoodArray,
                                                                       VVdouble &characterArray,
                                                                       const size_t siteNumber,
                                                                       const int newRootId) const {

    for (size_t i = 0; i < node->getNumberOfSons(); i++) {
        recursiveComputeLikelihoodAtSite(node->getSon(i), tree, likelihoodArray, characterArray, siteNumber, newRootId);
    }
    // Compute Pupko's likelihood using interface defined:
    likelihood_->computeLikelihoodAtSite(node, tree, likelihoodArray, characterArray,  siteNumber, newRootId);// From PIPDRHomogeneousTreeLikelihood
    DVLOG(2) << "[PIP ASR] The likelihood value for node number " << node->getId() << " with name " << node->getName()
             << " @site " << siteNumber << " is compute successfully.";


}


/**********************************************************************************************************************/
// Joint version //todo: Marginal should be implemented.
Sequence *PIPAncestralStateReconstruction::getAncestralSequenceForNode(int nodeId) const {

    // Check to see nodeId is a member of internal nodes:
    std::vector<int> InternalNodes = tree_.getInnerNodesId();
    if (!(std::count(InternalNodes.begin(), InternalNodes.end(), nodeId)))
    {
        throw Exception(
                "PIPAncestralStateReconstruction::getAncestralSequenceForNode(). nodeId is not member of internal nodes. (pick another node)");
    }

    // Do the joint reconstruction:
    std::map<int, std::vector<int> > allNodes{};
    VVVdouble probabilityProfile {};
    tie(allNodes, probabilityProfile) = getAllAncestralStatesWGap();
    //Allocating the new array for sequence:
    std::vector<int> allStates(nbSites_);

    std::string name = tree_.hasNodeName(nodeId) ? tree_.getNodeName(nodeId) : ("" + TextTools::toString(nodeId));


    DVLOG(2) << "[PIP ASR] Ancestral Sequence for node (" << nodeId << ") @node ";

    // Mapping the shrunk data to the normal one:
    for (int i = 0; i < nbSites_; ++i) {
        // Extracting the position of the site in the shrunkData
        allStates.at(i) = allNodes[nodeId][rootPatternLinks_[i]];
        DVLOG(2) << "[PIP ASR] " << "COL " << i << " is: " << (mlIndelPoints_->getShrunkData()->getSite(i)).toString()
                 << std::endl;
        DVLOG(2) << "[PIP ASR] Ancestral Sequence reconstruction for Site (" << i << ") in Shrunk data ("
                 << rootPatternLinks_[i] << ")" << " is: " << allStates.at(i);
    }
    return new BasicSequence(name, allStates, alphabet_);
}

/**********************************************************************************************************************/
//// Joint version //todo: Marginal should be implemented.
//std::vector<size_t> PIPAncestralStateReconstruction::getAncestralStatesForNode(int nodeId) const{
//
//    // Check to see nodeId is a member of internal nodes:
//    std::vector<int> InternalNodes = tree_.getInnerNodesId();
//    if (!(std::count(InternalNodes.begin(), InternalNodes.end(), nodeId)))
//    {
//        throw Exception(
//                "PIPAncestralStateReconstruction::getAncestralStatesForNode(). nodeId is not member of internal nodes. (pick another node)");
//    }
//    // Do the joint reconstruction:
//    std::map<int, std::vector<int> > allNodes = getAllAncestralStatesWGap);
//
//    return allNodes[nodeId];
//}




/**********************************************************************************************************************/

AlignedSequenceContainer *PIPAncestralStateReconstruction::getAncestralSequences() const {

    // Do the joint reconstruction:
    VVVdouble probPorfile {};
    std::map<int, std::vector<int>> allNodesAllStates{};
    tie(allNodesAllStates, probPorfile) = getAllAncestralStatesWGap();
    DLOG(INFO) << "[PIP ASR] Joint Ancestral Sequence Reconstruction is done successfully.";

    // Allocate the new array for Sequences:
    AlignedSequenceContainer *asc = new AlignedSequenceContainer(alphabet_);

    // List of internal nodes:
    std::vector<int> internalNodesIndex = tree_.getInnerNodesId();

    for (size_t nodeIId = 0; nodeIId < nbInnerNodes_; nodeIId++) {

        size_t currentNodeId = internalNodesIndex[nodeIId];

        DVLOG(2) << "[PIP ASR] Ancestral Sequence for node (" << nodeIId << ") @node ";
        std::string name = tree_.hasNodeName(currentNodeId) ? tree_.getNodeName(currentNodeId) : ("" +
                                                                                      TextTools::toString(currentNodeId));

        std::vector<int> allStates;
        allStates.empty();
        allStates.resize(nbSites_);
        // Mapping the shrunk data to the normal one:
        for (size_t i = 0; i < nbSites_; i++) {

            DVLOG(2) << "[PIP ASR] " << "COL " << i << " is: "
                     << (mlIndelPoints_->getShrunkData()->getSite(i)).toString() << std::endl;

            size_t rootIndex = rootPatternLinks_[i];

            // Extracting the position of the site in the shrunkData
            allStates[i] = allNodesAllStates[currentNodeId][rootIndex];

            DVLOG(2) << "[PIP ASR] Ancestral Sequence reconstruction for Site (" << i << ") in Shrunk data ("
                     << rootIndex << ")" << " is: " << allStates.at(i);
        }

        BasicSequence *seq = new BasicSequence(name, allStates, alphabet_);

        asc->addSequence(*seq);
        DVLOG(2) << "[PIP ASR] Ancestral states for sequence number (node)" << nodeIId << "is" << seq;
        delete seq;
    }
    DLOG(INFO)<< "[PIP ASR] Ancestral sequences were mapped to the original order with "<< asc->getNumberOfSites() << "sites.";

    return asc;
}

/**********************************************************************************************************************/

std::tuple<AlignedSequenceContainer *, VVVdouble> PIPAncestralStateReconstruction::getAncestralSequencesWithProbability() const {

    // Do the joint reconstruction:
    VVVdouble PupkoLiklihood {};
    std::map<int, std::vector<int>> allNodesAllStates{};
    tie(allNodesAllStates, PupkoLiklihood) = getAllAncestralStatesWGap();
    DLOG(INFO) << "[PIP ASR] Joint Ancestral Sequence Reconstruction is done successfully.";


    // compute the profile probability for ancestral nodes:
    VVVdouble probabilityProfile {};
    probabilityProfile.resize(nbDistinctSites_);
    computeProbabilityProfileForAllSites(PupkoLiklihood, probabilityProfile);

    // Allocate the new array for Sequences:
    AlignedSequenceContainer *asc = new AlignedSequenceContainer(alphabet_);

    // List of internal nodes:
    std::vector<int> internalNodesIndex = tree_.getInnerNodesId();

    for (size_t nodeIId = 0; nodeIId < nbInnerNodes_; nodeIId++) {

        size_t currentNodeId = internalNodesIndex[nodeIId];

        DVLOG(2) << "[PIP ASR] Ancestral Sequence for node (" << nodeIId << ") @node ";
        std::string name = tree_.hasNodeName(currentNodeId) ? tree_.getNodeName(currentNodeId) : ("" +
                                                                                                  TextTools::toString(currentNodeId));

        std::vector<int> allStates;
        allStates.empty();
        allStates.resize(nbSites_);
        // Mapping the shrunk data to the normal one:
        for (size_t i = 0; i < nbSites_; i++) {

            DVLOG(2) << "[PIP ASR] " << "COL " << i << " is: "
                     << (mlIndelPoints_->getShrunkData()->getSite(i)).toString() << std::endl;// check if it should be likelihood_

            size_t rootIndex = rootPatternLinks_[i];

            // Extracting the position of the site in the shrunkData
            allStates[i] = allNodesAllStates[currentNodeId][rootIndex];

            DVLOG(2) << "[PIP ASR] Ancestral Sequence reconstruction for Site (" << i << ") in Shrunk data ("
                     << rootIndex << ")" << " is: " << allStates.at(i);
        }

        BasicSequence *seq = new BasicSequence(name, allStates, alphabet_);

        asc->addSequence(*seq);
        DVLOG(2) << "[PIP ASR] Ancestral states for sequence number (node)" << nodeIId << "is" << seq;
        delete seq;
    }
    DLOG(INFO)<< "[PIP ASR] Ancestral sequences were mapped to the original order with "<< asc->getNumberOfSites() << "sites.";

    return make_tuple(asc, probabilityProfile);
}

/**********************************************************************************************************************/

VVVdouble PIPAncestralStateReconstruction::getProbabilityProfileForAllSites() const{
    return probabilityProfile_;
}

/**********************************************************************************************************************/

void PIPAncestralStateReconstruction::computeProbabilityProfileForAllSites(VVVdouble &likelihoodArray, VVVdouble &probability) const {

//    *prob &likelihoodArray;
//    // get the internal nodes:
    std::vector<int> nodeIds = tree_.getInnerNodesId();
    auto m = mlIndelPoints_->getDeletionPoints();

    if (probProfileType_ == "raw") {
        VVVdouble *oLik = &likelihoodArray;
        VVVdouble *oProb = &probability;
        for (size_t s = 0; s < nbDistinctSites_ ; s++) {
            VVdouble *likelihoodArray_s = &(*oLik)[s];
            VVdouble *probability_s = &(*oProb)[s];
            probability_s->resize(nbInnerNodes_);
            for (size_t n = 0; n < nbInnerNodes_ ; n++) {
                //check to see if the vector is empty or not:
                if (!(&(*likelihoodArray_s)[nodeIds[n]])->empty()) {
                    Vdouble *likelihoodArray_s_n = &(*likelihoodArray_s)[nodeIds[n]];
                    Vdouble *probability_s_n  = &(*probability_s)[n];
                    probability_s_n ->resize(nbStates_);
                    for (size_t x = 0; x < nbStates_; x++) {
                        (*probability_s_n)[x] = (*likelihoodArray_s_n)[x];
                    }
                }
            }
        }

    } else if(probProfileType_ == "normalized"){ // normalizing probabilities:
        VVVdouble *oLik = &likelihoodArray;
        VVVdouble *oProb = &probability;
        for (size_t s = 0; s < nbDistinctSites_ ; s++) {
            VVdouble *likelihoodArray_s = &(*oLik)[s];
            VVdouble *probability_s = &(*oProb)[s];
            probability_s->resize(nbInnerNodes_);
            for (size_t n = 0; n < nbInnerNodes_ ; n++) {
                //check to see if the vector is empty or not:
                if (!(&(*likelihoodArray_s)[nodeIds[n]])->empty()) {
                    Vdouble *likelihoodArray_s_n = &(*likelihoodArray_s)[nodeIds[n]];
                    Vdouble *probability_s_n  = &(*probability_s)[n];
                    probability_s_n ->resize(nbStates_);
                    // divide each element by the sum of likelihoodArray_s_n:
                    double sumLikelihoods = std::accumulate(likelihoodArray_s_n->begin(), likelihoodArray_s_n->end(), 0.0);
                    for (size_t x = 0; x < nbStates_; x++) {
                        (*probability_s_n)[x] = (*likelihoodArray_s_n)[x] / sumLikelihoods;
                    }
                }
            }
        }

    } else if(probProfileType_ == "naive_posterior"){ // compute naive posterior probabilities:

        // get background frequencies for the substitution model:
        const Vdouble *p = &mlIndelPoints_->getModel()->getFrequencies();
        VVVdouble *oLik = &likelihoodArray;
        VVVdouble *oProb = &probability;
        for (size_t s = 0; s < nbDistinctSites_ ; s++) {
            VVdouble *likelihoodArray_s = &(*oLik)[s];
            VVdouble *probability_s = &(*oProb)[s];
            probability_s->resize(nbInnerNodes_);
            for (size_t n = 0; n < nbInnerNodes_ ; n++) {
                //check to see if the vector is empty or not:
                if (!(&(*likelihoodArray_s)[nodeIds[n]])->empty()) {
                    Vdouble *likelihoodArray_s_n = &(*likelihoodArray_s)[nodeIds[n]];
                    Vdouble *probability_s_n = &(*probability_s)[n];
                    probability_s_n->resize(nbStates_);
                    // divide each element by the sum product of p and likelihoodArray_s_n:
                    double sum_product {0};
                    for (size_t x = 0; x < nbStates_; x++) {
                        (*probability_s_n)[x] = (*likelihoodArray_s_n)[x] * p->at(x);
                        sum_product += (*probability_s_n)[x];
                    }
                    for(size_t y=0; y<nbStates_; y++){
                        (*probability_s_n)[y] = (*probability_s_n)[y] / sum_product;
                    }
                }
            }
        }

    } else {

        throw Exception("PIPAncestralStateReconstruction::computeProbabilityProfileForAllSites. Error, unknown profile type!");
        LOG(ERROR)<< "[PIP ASR] Error, unknown profile type!";

    }
    DLOG(INFO) << "[PIP ASR] Probability profile for all sites is computed successfully.";
}

/**********************************************************************************************************************/
