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
 * @file PIPMLIndelPoints.cpp
 * @authors Massimo Maiolo and Gholamhossein Jowkar
 * @date 06 07 2020
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
#include <zconf.h>
#include "PIPMLIndelPoints.h"
using namespace bpp;

/************************************************* Constructor ********************************************************/
PIPMLIndelPoints::PIPMLIndelPoints(const PIPDRHomogeneousTreeLikelihood *lik) :
        likelihood_(lik),
        tree_(0),
        shrunkData_(0),
        model_(0),
        nbSites_(lik->getNumberOfSites()),
        nbStates_(lik->getNumberOfStates()),
        nbClasses_(lik->getNumberOfClasses()),

        nodeIota_(),nodeBeta_(),
        nodeSetG_(), nodeSetA_(),

        nodeText_(), nodeProb_(),
        nodeFelsensteinValueText_(), nodeFelsensteinValue_() , homoPath_(),
        insertionPoints_(), deletionPoints_()
        {
    // The same data used by likelihood function: sometimes the order of columns is changed.
    nbDistinctSites_ = lik->getLikelihoodData()->getNumberOfDistinctSites();
    shrunkData_ = lik->getLikelihoodData()->getShrunkData();
    tree_ = new TreeTemplate<Node> (lik->getTree());

    if(lik->isPip()) {
        model_ = lik->getModel();
        nodeIota_ = &(lik->getLikelihoodData()->getPipParam().getNodeIotaData());
        nodeBeta_ = &(lik->getLikelihoodData()->getPipParam().getNodeBetaData());
        nodeSetG_ = (lik->getLikelihoodData()->getNodeSetG());
        nodeSetA_ = (lik->getLikelihoodData()->getNodeSetA());
        DLOG(INFO) << "[Indel points] In constructor the PIP parameter is copied from the likelihood object.";
    }
    else {
        ApplicationTools::displayError(
                "[Indel points] This class only supports PIP based models. (meaning PIP13,...)");
        DLOG(ERROR) << "[Indel points] The likelihood object is not PIP-based.";
    }

    // Initializing the algorithm
    initMaximumLikelihoodIndelPoints_(tree_->getRootNode());
    DLOG(INFO) << "[Indel points] IndelPoints object is initialized successfully.";

}


/********************************************* Copy constructor *******************************************************/

PIPMLIndelPoints::PIPMLIndelPoints(const PIPMLIndelPoints &mlindel) :
        tree_(mlindel.tree_),
        shrunkData_(mlindel.shrunkData_),
        model_(mlindel.model_),
        nbSites_(mlindel.nbSites_),
        nbStates_(mlindel.nbStates_),
        nbClasses_(mlindel.nbClasses_),
        nbDistinctSites_(mlindel.nbDistinctSites_),

        nodeIota_(mlindel.nodeIota_),
        nodeBeta_(mlindel.nodeBeta_),
        nodeSetG_(mlindel.nodeSetG_),
        nodeSetA_(mlindel.nodeSetA_),

        nodeText_(mlindel.nodeText_),
        nodeProb_(mlindel.nodeProb_),
        nodeFelsensteinValueText_(mlindel.nodeFelsensteinValueText_),
        nodeFelsensteinValue_(mlindel.nodeFelsensteinValue_),
        homoPath_(mlindel.homoPath_),
        insertionPoints_(mlindel.insertionPoints_),
        deletionPoints_(mlindel.deletionPoints_) {
    if (mlindel.shrunkData_)
        shrunkData_ = mlindel.shrunkData_->clone();

}

/********************************************* Operator overloading *******************************************************/

PIPMLIndelPoints &PIPMLIndelPoints::operator=(const PIPMLIndelPoints &mlindel) {

    tree_                       = mlindel.tree_;
    shrunkData_                 = mlindel.shrunkData_;
    model_                      = mlindel.model_;
    nbSites_                    = mlindel.nbSites_;
    nbStates_                   = mlindel.nbStates_;
    nbClasses_                  = mlindel.nbClasses_;
    nbDistinctSites_            = mlindel.nbDistinctSites_;

    nodeIota_                   = mlindel.nodeIota_;
    nodeBeta_                   = mlindel.nodeBeta_;
    nodeSetA_                   = mlindel.nodeSetA_;
    nodeSetG_                   = mlindel.nodeSetG_;

    nodeText_                   = mlindel.nodeText_;
    nodeProb_                   = mlindel.nodeProb_;
    nodeFelsensteinValueText_   = mlindel.nodeFelsensteinValueText_;
    nodeFelsensteinValue_       = mlindel.nodeFelsensteinValue_;
    homoPath_                   = mlindel.homoPath_;

    insertionPoints_            = mlindel.insertionPoints_;
    deletionPoints_             = mlindel.deletionPoints_;

}

/******************************************* PIP specific functions ***************************************************/
void PIPMLIndelPoints::initMaximumLikelihoodIndelPoints_(const Node *node)  {


    /* Getting Total tree length: */
//    const double tau = TreeTemplateTools::getTotalLength(*tree_->getRootNode(), false);

    /* Getting the parameter mu: */
    ParameterList PIPpl = model_->getParameters();
    const float mu = PIPpl.getParameterValue("PIP13.mu");
    /* INFO: The prefix "PIP13." is the namespace of the parameters, to avoid confusion when more complicated models are used.
     * Parameter values can also be directly accessed using:
     * cout << model->getParameterValue("PIP13.mu") << endl;
     * This syntax is then independent of the namespace used.
     */


//    std::cout << std::endl << "-------------|" << "extractIndelPoints" << "|--------------" << std::endl;
    ApplicationTools::displayTask("Extracting the Indel Points");
    extractIndelPoint(node, mu);
    DLOG(INFO)<< "[Indel points] Indel points were inferred successfully.";
    ApplicationTools::displayTaskDone();
    //    MLHomoPath(nodeText_[tree_->getRootId()]);
    MLHomoPath(nodeText_[tree_->getRootId()]);
}

/******************************************** extractIndelPoint ********************************************/

void PIPMLIndelPoints::extractIndelPoint(const Node *node, const float mu) {
    for (size_t siteNumber = 0; siteNumber < nbDistinctSites_; siteNumber++) {
        extractIndelPointForSite_(node, mu, siteNumber);
    }
}

/****************************************** extractIndelPointForSite ******************************************/

void PIPMLIndelPoints::extractIndelPointForSite_(const Node *node, const float mu, size_t siteNumber) {

    // To have access to nodeAllGap_ element of current node:
    int nodeId = node->getId();
    std::map<size_t, std::vector<bool> > *node_all_gap_ = &nodeSetG_;
    std::vector<bool> *node_all_gap_n_ = &(*node_all_gap_)[nodeId];
    if (node->isLeaf()) {
        /******************************* Leaf Node ************************************/
        leafHomologyPath_(node,siteNumber);
    } else {
        /**************************** Internal Node **********************************/
        // We initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();
        if (nbSonNodes != 2) throw Exception("PIPDRTreeLikelihoodData::extractIndelPoint. The input tree should be binary tree.");

        const Node *leftChildNode = (node->getSon(0));
        const Node *rightChildNode = (node->getSon(1));

        extractIndelPointForSite_(leftChildNode, mu, siteNumber);
        extractIndelPointForSite_(rightChildNode, mu, siteNumber);



        if ((*node_all_gap_n_)[siteNumber]) {
            /**------------------ Both subtrees are gaps ---------------------------%
            %-----------------------------------------------------------------------%
            %    o
            %   / \
            %  -   -
            %-----------------------------------------------------------------------%
            %    o
            %   / \
            %  -   o
            %     / \
            %    -   -
            %-----------------------------------------------------------------------%
            %     o
            %    / \
            %   o   -
            %  / \
            % -   -
            %----------------------------------------------------------------------%
            %      o
            %     / \
            %    /   \
            %   o     o
            %  / \   / \
            % -   - -   -
            %----------------------------------------------------------------------%
            */
            GapOnlySubtreeHomologyPath_(node, leftChildNode, rightChildNode, mu, siteNumber);
        } else {
            /*************** At least one subtree contains a character ************/
            if (((*node_all_gap_)[leftChildNode->getId()][siteNumber]) == 0 &&
                ((*node_all_gap_)[rightChildNode->getId()][siteNumber] == 0)) {
                /*%-------------------------------------------------------------%
                %    o
                %   / \
                %  A   A
                %-------------------------------------------------------------%
                %      o
                %     / \
                %    A   o
                %       / \
                %      A   A
                %-------------------------------------------------------------%
                %      o
                %     / \
                %    o   A
                %   / \
                %  A   A
                %-------------------------------------------------------------%
                %      o
                %     / \
                %    /   \
                %   o     o
                %  / \   / \
                % A   A A   A
                %-------------------------------------------------------------%
                 */
                NeitherChildIsGapOnlySubtreeHomologyPath_(node, leftChildNode, rightChildNode, mu, siteNumber);
            } else {
                if (leftChildNode->isLeaf() && rightChildNode->isLeaf()) {
                    /******************* first internal node **********************/
                    /** %---------------------------------------------------------%
                    %    o
                    %   / \
                    %  -   A
                    %---------------------------------------------------------%
                    %    o
                    %   / \
                    %  A   -
                    %---------------------------------------------------------%
                     */
                    BothChildrenLeavesOrOneLeafSubtreeHomologyPath_(node, leftChildNode, rightChildNode, mu, siteNumber);
                } else if ((!leftChildNode->isLeaf()) != (!rightChildNode->isLeaf())) {//logical xor
                    /** %------------------------------------------------------%
                    %      o
                    %     / \
                    %    -   o
                    %       / \
                    %      A   A
                    %---------------------------------------------------------%
                    %      o
                    %     / \
                    %    A   o
                    %       / \
                    %      -   -
                    %---------------------------------------------------------%
                    %      o
                    %     / \
                    %    o   A
                    %   / \
                    %  -   -
                    %---------------------------------------------------------%
                    %     o
                    %    / \
                    %   o   -
                    %  / \
                    % A   A
                    %---------------------------------------------------------%
                     */
                    OneChildLeafSubtreeHomologyPath_(node, leftChildNode, rightChildNode, mu, siteNumber);
                } else {
                    /** %---------------------------------------------------------%
                    %      o
                    %     / \
                    %    /   \
                    %   o     o
                    %  / \   / \
                    % -   - A   A
                    %---------------------------------------------------------%
                    %      o
                    %     / \
                    %    /   \
                    %   o     o
                    %  / \   / \
                    % A   A -   -
                    %---------------------------------------------------------%
                    */
                    NoChildLeafHomologyPath_(node, leftChildNode, rightChildNode, mu, siteNumber);
                }
            }
        }
    }
}

/******************************************** leafHomologyPath_ ********************************************/

void PIPMLIndelPoints::leafHomologyPath_(const Node *node, size_t siteNumber) {

    size_t nodeId = node->getId();
    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);


//    Vdouble *node_iota_ = &nodeIota_;
//    Vdouble node_beta_ = nodeBeta_;

//    VVdouble* leavesLikelihoods_leaf = &leafData_[node->getId()].getLikelihoodArray();

    std::map<size_t, std::vector<bool>> *setA_node_ = &nodeSetA_;
    std::vector<bool> *setA_node_n_ = &(*setA_node_)[nodeId];

    if ((*setA_node_n_)[siteNumber]) {
        (*node_prob_n_)[siteNumber] = (*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId];
//        (*node_text_n_).push_back(node->getId() + ":I;");
        (*node_text_n_)[siteNumber] = std::to_string(node->getId()) + ":I;";
    } else {
        (*node_prob_n_)[siteNumber] = 0.;
        (*node_text_n_)[siteNumber] = "";
    }

    (*node_fv_n_)[siteNumber] = 1.;
    (*node_fv_text_n_)[siteNumber] = "";

}

/********************************************* GapOnlySubtreeHomologyPath_ *********************************/

void PIPMLIndelPoints:: GapOnlySubtreeHomologyPath_(
        const Node *node,
        const Node *childNode1,
        const Node *childNode2,
        const float mu,
        const size_t siteNumber)
        {
    /*%-----------------------------------------------------------------%
    %- both subtrees are gaps
    %-----------------------------------------------------------------%
    %    o
    %   / \
    %  -   -
    %-----------------------------------------------------------------%
    %    o
    %   / \
    %  -   o
    %     / \
    %    -   -
    %-----------------------------------------------------------------%
    %     o
    %    / \
    %   o   -
    %  / \
    % -   -
    %-----------------------------------------------------------------%
    %      o
    %     / \
    %    /   \
    %   o     o
    %  / \   / \
    % -   - -   -
    %-----------------------------------------------------------------%*/
    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

//    VVdouble* leavesLikelihoods_leaf = &leafData_[node->getId()].getLikelihoodArray();

    // Definition of local variables:
    double fv1{0}, fv2{0};
    fv1 = (childNode1->isLeaf() ? 0. : (*node_fv_)[child1NodeId][siteNumber]);
    fv2 = (childNode2->isLeaf() ? 0. : (*node_fv_)[child2NodeId][siteNumber]);

    Vdouble fvR, fvL;
    fvL.push_back(fv1 * exp((-mu) * (childNode1->getDistanceToFather())));
    fvL.push_back(1 - exp((-mu) * (childNode1->getDistanceToFather())));

    fvR.push_back(fv2 * exp((-mu) * (childNode2->getDistanceToFather())));
    fvR.push_back(1 - exp((-mu) * (childNode2->getDistanceToFather())));

    std::vector<std::string> fv_txt_1, fv_txt_2;
    fv_txt_1.push_back((*node_fv_text_)[child1NodeId][siteNumber]);
    fv_txt_2.push_back((*node_fv_text_)[child2NodeId][siteNumber]);

    fv_txt_1.push_back(std::to_string(childNode1->getId())+ ":X;");
    fv_txt_2.push_back(std::to_string(childNode2->getId())+ ":X;");

    // Inline definition with whichMax and max doesn't work!
    size_t idxl = VectorTools::whichMax(fvL);
    size_t idxr = VectorTools::whichMax(fvR);

    double fvL_max = VectorTools::max(fvL);
    double fvR_max = VectorTools::max(fvR);

    // Doesn't work and I don't know why!!!
//    std::cout << VectorTools::max(fvL);
//    std::cout << fv_txt_1[VectorTools::whichMax(fvL)];


    (*node_prob_n_)[siteNumber] = 0.;
    (*node_text_n_)[siteNumber] = "";
    (*node_fv_n_)[siteNumber] = fvL_max * fvR_max;
    (*node_fv_text_n_)[siteNumber] = fv_txt_1[idxl] + fv_txt_2[idxr];

}

/******************************************* NeitherChildIsGapOnlySubtreeHomologyPath_ ********************************/

void PIPMLIndelPoints::NeitherChildIsGapOnlySubtreeHomologyPath_(
        const Node *node,
        const Node *childNode1,
        const Node *childNode2,
        const float mu,
        const size_t siteNumber){
    /*%-------------------------------------------------------------------------%
    %    o
    %   / \
    %  A   A
    %-------------------------------------------------------------------------%
    %      o
    %     / \
    %    A   o
    %       / \
    %      A   A
    %-------------------------------------------------------------------------%
    %      o
    %     / \
    %    o   A
    %   / \
    %  A   A
    %-------------------------------------------------------------------------%
    %      o
    %     / \
    %    /   \
    %   o     o
    %  / \   / \
    % A   A A   A
    %-------------------------------------------------------------------------% */
    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

//    Vdouble *node_iota_ = &(nodeIota_);
//    Vdouble *node_beta_ = &(nodeBeta_);

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

    std::map<size_t, std::vector<bool>> *node_setA_ = &nodeSetA_;
    std::vector<bool> *node_setA_n_ = &(*node_setA_)[nodeId];

    // Definition of local variables:
    double fv1{0}, fv2{0};
    fv1 = (*node_fv_)[child1NodeId][siteNumber];
    fv2 = (*node_fv_)[child2NodeId][siteNumber];

    Vdouble fvR, fvL;
    fvL.push_back(fv1 * exp((-mu) * (childNode1->getDistanceToFather())));
    fvR.push_back(fv2 * exp((-mu) * (childNode2->getDistanceToFather())));

    std::vector<std::string> fv_txt_1, fv_txt_2;

    fv_txt_1.push_back((*node_fv_text_)[child1NodeId][siteNumber]);
    fv_txt_2.push_back((*node_fv_text_)[child2NodeId][siteNumber]);

    Vdouble fv;
    std::vector<std::string> txt;
    for (size_t i = 0; i < fv_txt_1.size(); i++) {
        for (size_t j = 0; j < fv_txt_2.size(); j++) {
            fv.push_back(fvL[i] * fvR[j]);
            txt.push_back(fv_txt_1[i] + fv_txt_2[j]);
        }
    }

    (*node_fv_n_)[siteNumber] = VectorTools::max(fv);
    (*node_fv_text_n_)[siteNumber] = txt[VectorTools::whichMax(fv)];

    if ((*node_setA_n_)[siteNumber] == 1) {
        (*node_prob_n_)[siteNumber] = (*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId] * (*node_fv_n_)[siteNumber];
        (*node_text_n_)[siteNumber] = (*node_fv_text_n_)[siteNumber] + std::to_string(node->getId()) + ":I;";
    } else {
        (*node_prob_n_)[siteNumber] = 0.;
        (*node_text_n_)[siteNumber] = "";
    }
}

/**************************************** BothChildrenLeavesOrOneLeafSubtreeHomologyPath_ ******************************/

void PIPMLIndelPoints::BothChildrenLeavesOrOneLeafSubtreeHomologyPath_(
        const Node *node,
        const Node *childNode1,
        const Node *childNode2,
        const float mu,
        const size_t siteNumber) {

    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

    // To have access to nodeAllGap_ element of current node:
    std::map<size_t, std::vector<bool> > *node_all_gap_ = &nodeSetG_;
    std::map<size_t, std::vector<bool> > *setA_node_ = &nodeSetA_;
    std::vector<bool> *node_all_gap_n_ = &(*node_all_gap_)[nodeId];
    std::vector<bool> *setA_node_n_ = &(*setA_node_)[nodeId];

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

//    Vdouble *node_iota_ = &(nodeIota_);
//    Vdouble *node_beta_ = &(nodeBeta_);

    // Local and tmp variables:
    Vdouble p;
    std::vector<std::string> txt;

    if (((*node_all_gap_)[child1NodeId])[siteNumber]) {
        OneChildGapHomologyPath_(node, childNode1, childNode2, mu, siteNumber);
        if ((*setA_node_n_)[siteNumber]) {

            p.push_back((*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId] * (*node_fv_n_)[siteNumber]);
            p.push_back((*node_prob_)[child2NodeId][siteNumber]);

            txt.push_back((*node_fv_text_n_)[siteNumber] + std::to_string(node->getId()) + ":I;");// siteNumber -> txt.push_back((*node_fv_text_n_)[siteNumber]
            txt.push_back((*node_text_)[child2NodeId][siteNumber]);

            (*node_prob_n_)[siteNumber] = VectorTools::max(p);
            (*node_text_n_)[siteNumber] = txt[VectorTools::whichMax(p)];

        } else {
            (*node_prob_n_)[siteNumber] = 0;
            (*node_text_n_)[siteNumber] = "";
        }
    } else {
        OneChildGapHomologyPath_(node, childNode2, childNode1, mu, siteNumber);//Children are swapped.
        if ((*setA_node_n_)[siteNumber]) {

            p.push_back((*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId] * (*node_fv_n_)[siteNumber]);
            p.push_back((*node_prob_)[child1NodeId][siteNumber]);

            txt.push_back((*node_fv_text_n_)[siteNumber] + std::to_string(node->getId()) + ":I;");// 0 or siteNumber -> txt.push_back((*node_fv_text_n_)[siteNumber]
            txt.push_back((*node_text_)[child1NodeId][siteNumber]);

            (*node_prob_n_)[siteNumber] = VectorTools::max(p);
            (*node_text_n_)[siteNumber] = txt[VectorTools::whichMax(p)];

        } else {
            (*node_prob_n_)[siteNumber] = 0;
            (*node_text_n_)[siteNumber] = "";
        }
    }
}

/******************************************* OneChildLeafSubtreeHomologyPath_ ***********************************/

void PIPMLIndelPoints::OneChildLeafSubtreeHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2, const float mu,
                                 const size_t siteNumber) {

    if (childNode1->isLeaf() || childNode2->isLeaf()) {
        /* %---------------------------------------------------------%
        %      o
        %     / \
        %    -   o
        %       / \
        %      A   A
        %---------------------------------------------------------%*/
        BothChildrenLeavesOrOneLeafSubtreeHomologyPath_(node, childNode1, childNode2, mu, siteNumber);
    } else {
        /*%---------------------------------------------------------%
        %      o
        %     / \
        %    o   A
        %   / \
        %  -   -
        %---------------------------------------------------------%*/
        NoChildLeafHomologyPath_(node, childNode1, childNode2, mu, siteNumber);
    }
}

/*********************************************** NoChildLeafHomologyPath_ *******************************/

void PIPMLIndelPoints::NoChildLeafHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2, const float mu,
                                 const size_t siteNumber) {

    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

    // To have access to nodeAllGap_ element of current node:
    std::map<size_t, std::vector<bool> > *node_all_gap_ = &nodeSetG_;
    std::map<size_t, std::vector<bool> > *node_setA_ = &nodeSetA_;
    std::vector<bool> *node_all_gap_n_ = &(*node_all_gap_)[nodeId];
    std::vector<bool> *setA_node_n_ = &(*node_setA_)[nodeId];

//    Vdouble *node_iota_ = &(nodeIota_);
//    Vdouble *node_beta_ = &(nodeBeta_);

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

    // Local and tmp variable:
    Vdouble p;
    std::vector<std::string> txt;

    if ((*node_all_gap_)[child1NodeId][siteNumber]) {
        /*%---------------------------------------------------------%
        %      o
        %     / \
        %    /   \
        %   o     o
        %  / \   / \
        % -   - A   A
        %---------------------------------------------------------%*/
        AllChildGapHomologyPath_(node, childNode2, childNode1, mu, siteNumber);//Children are swapped.

        if ((*setA_node_n_)[siteNumber]) {
            p.push_back((*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId] * (*node_fv_n_)[siteNumber]);
            p.push_back((*node_prob_)[child2NodeId][siteNumber]);

            txt.push_back((*node_fv_text_n_)[siteNumber] + std::to_string(node->getId()) + ":I;");
            txt.push_back((*node_text_)[child2NodeId][siteNumber]);

            (*node_prob_n_)[siteNumber] = VectorTools::max(p);
            (*node_text_n_)[siteNumber]= txt[VectorTools::whichMax(p)];

        } else {
            (*node_prob_n_)[siteNumber] = 0;
            (*node_text_n_)[siteNumber] = "";
        }
    } else {
        /*%---------------------------------------------------------%
        %      o
        %     / \
        %    /   \
        %   o     o
        %  / \   / \
        % A   A -   -
        %---------------------------------------------------------%*/
        AllChildGapHomologyPath_(node, childNode1, childNode2, mu, siteNumber);//Children are swapped.
        if ((*setA_node_n_)[siteNumber]) {

            p.push_back((*nodeIota_)[nodeId] * (*nodeBeta_)[nodeId] * (*node_fv_n_)[siteNumber]);
            p.push_back(((*node_prob_)[child1NodeId][siteNumber]));

            txt.push_back((*node_fv_text_n_)[siteNumber] + std::to_string(node->getId()) + ":I;");// siteNumber -> txt.push_back((*node_fv_text_n_)[siteNumber]
            txt.push_back((*node_text_)[child1NodeId][siteNumber]);

            (*node_prob_n_)[siteNumber] = VectorTools::max(p);
            (*node_text_n_)[siteNumber] = txt[VectorTools::whichMax(p)];
        } else {
            (*node_prob_n_)[siteNumber] = 0;
            (*node_text_n_)[siteNumber] = "";
        }
    }
}

/***************************************** OneChildGapHomologyPath_ *************************************/

void PIPMLIndelPoints::OneChildGapHomologyPath_(
        const Node *node,
        const Node *childNode1,
        const Node *childNode2,
        const float mu,
        const size_t siteNumber) {
    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

    std::map<size_t, std::vector<bool>> *node_setA_ = &nodeSetA_;
    std::vector<bool> *node_setA_n_ = &(*node_setA_)[nodeId];

    // Definition of local variables:
    double fv2{0};
//    fv1 = (*node_fv_)[child1NodeId][siteNumber];
    fv2 = (*node_fv_)[child2NodeId][siteNumber];

    Vdouble fvR, fvL;
    fvL.push_back(1- exp((-mu) * (childNode1->getDistanceToFather())));
    fvR.push_back(fv2 * exp((-mu) * (childNode2->getDistanceToFather())));

    std::vector<std::string> fv_txt_1, fv_txt_2;

    fv_txt_1.push_back(std::to_string(childNode1->getId()) + ":X;");
    fv_txt_2.push_back((*node_fv_text_)[child2NodeId][siteNumber]);

    Vdouble fv;
    std::vector<std::string> txt;
//    for (size_t i = 0; i < fv_txt_1.size(); i++) {
    for (size_t j = 0; j < fv_txt_2.size(); j++) {
        fv.push_back(fvL[0] * fvR[j]);
        txt.push_back(fv_txt_1[0] + fv_txt_2[j]);
    }
//    }

    (*node_fv_n_)[siteNumber] = VectorTools::max(fv);
    (*node_fv_text_n_)[siteNumber] = txt[VectorTools::whichMax(fv)];


}

/************************************* AllChildGapHomologyPath_ *****************************************/

void PIPMLIndelPoints::AllChildGapHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2,
                                 const float mu, const size_t siteNumber) {

    size_t nodeId = node->getId();
    size_t child1NodeId = childNode1->getId();
    size_t child2NodeId = childNode2->getId();

//    Vdouble *node_iota_ = &(nodeIota_);
//    Vdouble *node_beta_ = &(nodeBeta_);

    // Initialize Text and Prob vector:
    std::map<int, std::vector<std::string>> *node_text_ = &nodeText_;
    std::map<int, Vdouble> *node_prob_ = &nodeProb_;

    std::vector<std::string> *node_text_n_ = &(*node_text_)[nodeId];
    Vdouble *node_prob_n_ = &(*node_prob_)[nodeId];
    node_text_n_->resize(nbDistinctSites_);
    node_prob_n_->resize(nbDistinctSites_);

    std::map<int, std::vector<std::string>> *node_fv_text_ = &nodeFelsensteinValueText_;
    std::map<int, Vdouble> *node_fv_ = &nodeFelsensteinValue_;

    std::vector<std::string> *node_fv_text_n_ = &(*node_fv_text_)[nodeId];
    Vdouble *node_fv_n_ = &(*node_fv_)[nodeId];
    node_fv_n_->resize(nbDistinctSites_);
    node_fv_text_n_->resize(nbDistinctSites_);

    std::map<size_t, std::vector<bool>> *node_setA_ = &nodeSetA_;
    std::vector<bool> *node_setA_n_ = &(*node_setA_)[nodeId];

    // Definition of local variables:
    double fv1{0}, fv2{0};
    fv1 = (*node_fv_)[child1NodeId][siteNumber];
    fv2 = (*node_fv_)[child2NodeId][siteNumber];

    Vdouble fvR, fvL;
    fvR.push_back(fv2 * exp((-mu) * childNode2->getDistanceToFather()));
    fvR.push_back(1 - (exp((-mu) * childNode2->getDistanceToFather())));

    fvL.push_back(fv1 * exp((-mu) * childNode1->getDistanceToFather()));


    std::vector<std::string> fv_txt_1, fv_txt_2;

    fv_txt_1.push_back((*node_fv_text_)[child1NodeId][siteNumber]);
    fv_txt_2.push_back((*node_fv_text_)[child2NodeId][siteNumber]);

    std::vector<std::string> txt;//might be problematic
    txt = fv_txt_2;
    txt.push_back(std::to_string(childNode2->getId()) + ":X;");


    (*node_fv_n_)[siteNumber] = VectorTools::max(fvL) * VectorTools::max(fvR);
    (*node_fv_text_n_)[siteNumber] = fv_txt_1[VectorTools::whichMax(fvL)] + txt[VectorTools::whichMax(fvR)];
//
}

/****************************************** MLHomoPath ************************************/

void PIPMLIndelPoints::MLHomoPath(std::vector<std::string> &text) {
    homoPath_ = text;
    std::map<size_t, std::vector<std::string>> *insertion_points_ = &(insertionPoints_);
    std::map<size_t, std::vector<std::string >> *deletion_points_ = &(deletionPoints_);

    for (size_t nodeId = 0; nodeId < text.size(); nodeId++) {
        std::vector<std::string> *insertion_points_n_ = &(*insertion_points_)[nodeId];
        std::vector<std::string> *deletion_points_n_ = &(*deletion_points_)[nodeId];

        std::string del = "";
        std::string ins = "";
        std::string s = text[nodeId];
        size_t pos = 0;
        std::string delimiter = ";";
        std::string token;
        while ((pos = s.find(delimiter))!= std::string::npos) {
            token = s.substr(0, pos);
            if (!token.empty()) {
//                std::cout << token << std::endl;
                if (token.find(":X") != std::string::npos) {
                    deletion_points_n_->push_back(token.substr(0, token.find(":")));
                } else {
                    insertion_points_n_->push_back(token.substr(0, token.find(":")));
                }
            }
            s.erase(0, pos + delimiter.length());
        }
    }
}
/************************************* PIP getters and setters *************************************************/
const std::map<size_t, std::vector<bool>> &PIPMLIndelPoints::getNodeSetG() const {
    return nodeSetG_;
}

void PIPMLIndelPoints::setNodeSetG(const std::map<size_t, std::vector<bool>> &nodeSetG) {
    nodeSetG_ = nodeSetG;
}

const std::map<size_t, std::vector<bool>> &PIPMLIndelPoints::getNodeSetA() const {
    return nodeSetA_;
}

void PIPMLIndelPoints::setNodeSetA(const std::map<size_t, std::vector<bool>> &nodeSetA) {
    nodeSetA_ = nodeSetA;
}

const std::map<int, std::vector<std::string>> &PIPMLIndelPoints::getNodeText() const {
    return nodeText_;
}

void PIPMLIndelPoints::setNodeText(const std::map<int, std::vector<std::string>> &nodeText) {
    nodeText_ = nodeText;
}

const std::map<int, std::vector<std::string>> &PIPMLIndelPoints::getNodeFelsensteinValueText() const {
    return nodeFelsensteinValueText_;
}

void
PIPMLIndelPoints::setNodeFelsensteinValueText(const std::map<int, std::vector<std::string>> &nodeFelsensteinValueText) {
    nodeFelsensteinValueText_ = nodeFelsensteinValueText;
}

const std::map<int, Vdouble> &PIPMLIndelPoints::getNodeProb() const {
    return nodeProb_;
}

void PIPMLIndelPoints::setNodeProb(const std::map<int, Vdouble> &nodeProb) {
    nodeProb_ = nodeProb;
}

const std::map<int, Vdouble> &PIPMLIndelPoints::getNodeFelsensteinValue() const {
    return nodeFelsensteinValue_;
}

void PIPMLIndelPoints::setNodeFelsensteinValue(const std::map<int, Vdouble> &nodeFelsensteinValue) {
    nodeFelsensteinValue_ = nodeFelsensteinValue;
}

const std::vector<std::string> &PIPMLIndelPoints::getHomoPath() const {
    return homoPath_;
}

void PIPMLIndelPoints::setHomoPath(const std::vector<std::string> &homoPath) {
    homoPath_ = homoPath;
}

const std::map<size_t, std::vector<std::string>> &PIPMLIndelPoints::getInsertionPoints() const {
    return insertionPoints_;
}

void PIPMLIndelPoints::setInsertionPoints(const std::map<size_t, std::vector<std::string>> &insertionPoints) {
    insertionPoints_ = insertionPoints;
}

const std::map<size_t, std::vector<std::string>> &PIPMLIndelPoints::getDeletionPoints() const {
    return deletionPoints_;
}

void PIPMLIndelPoints::setDeletionPoints(const std::map<size_t, std::vector<std::string>> &deletionPoints) {
    deletionPoints_ = deletionPoints;
}





