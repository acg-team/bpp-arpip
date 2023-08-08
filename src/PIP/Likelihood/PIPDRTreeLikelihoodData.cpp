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
 * @file PIPDRTreeLikelihoodData.cpp
 * @author Gholamhossein Jowkar
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

// From STL
#include "vector"

#include "PIPDRTreeLikelihoodData.h"
//From PhylLib:
#include "Bpp/Phyl/PatternTools.h"

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <cstring>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

// From ARPIP
#include "../SubstitutionModel/PIP13.h"
#include "../Utils/ARPIPTools.h"
#include "PIPDRTreeLikelihoodData.h"
#include "PIPDRTreeLikelihood.h"

// From GLOG
#include <glog/logging.h>

using namespace bpp;



/******************************************************************************/
/************************ PIPDRTreeLikelihoodLeafData *************************/
/******************************************************************************/

/******************************** Constructor **********************************************/
PIPDRTreeLikelihoodLeafData::PIPDRTreeLikelihoodLeafData(): leafLikelihood_(), leaf_(0){}

/***************************** Copy Constructor ************************************************/
PIPDRTreeLikelihoodLeafData::PIPDRTreeLikelihoodLeafData(const PIPDRTreeLikelihoodLeafData &data):
        leafLikelihood_(data.leafLikelihood_), leaf_(data.leaf_){}

/******************************* Operator Overloading ***********************************************/
PIPDRTreeLikelihoodLeafData &PIPDRTreeLikelihoodLeafData::operator=(const PIPDRTreeLikelihoodLeafData &data) {
    leafLikelihood_ = data.leafLikelihood_;
    leaf_ = data.leaf_;

    return *this;
}

/******************************************************************************/
/************************ PIPDRTreeLikelihoodNodeData *************************/
/******************************************************************************/

/*********************************** Constructor **********************************************/
PIPDRTreeLikelihoodNodeData::PIPDRTreeLikelihoodNodeData()
        : nodeLikelihoods_(), node_(0), nodeDLikelihoods_(), nodeD2Likelihoods_()
        {}

/***************************** Copy Constructor ************************************************/
PIPDRTreeLikelihoodNodeData::PIPDRTreeLikelihoodNodeData(const PIPDRTreeLikelihoodNodeData &data) :
        nodeLikelihoods_(data.nodeLikelihoods_),
        node_(data.node_),

        // Derivatives
        nodeDLikelihoods_(data.nodeDLikelihoods_),
        nodeD2Likelihoods_(data.nodeD2Likelihoods_)
{}

/******************************* Operator Overloading *******************************************/
PIPDRTreeLikelihoodNodeData &PIPDRTreeLikelihoodNodeData::operator=(const PIPDRTreeLikelihoodNodeData &data)
{
    nodeLikelihoods_ = data.nodeLikelihoods_;
    node_ = data.node_;

    //Derivatives
    nodeDLikelihoods_ = data.nodeDLikelihoods_;
    nodeD2Likelihoods_ = data.nodeD2Likelihoods_;

    return *this;
}

/******************************************************************************/
/*********************** PIPDRTreeLikelihoodParameter *************************/
/******************************************************************************/

/*********************************** Constructor **********************************************/

PIPDRTreeLikelihoodParameters::PIPDRTreeLikelihoodParameters()
        : tree_(), tau_(0), lambda_(0), mu_(0), nu_(0), nodeIotaData_(),
          nodeBetaData_() {
}

/*********************************** Constructor **********************************************/

PIPDRTreeLikelihoodParameters::PIPDRTreeLikelihoodParameters(const TreeTemplate<Node> *tree,
                                                             const double lambda,
                                                             const double mu) :
        tree_(tree), tau_(0), lambda_(lambda), mu_(mu), nu_(0), nodeIotaData_(), nodeBetaData_() {

    initPIPParameter_();
}

/******************************** Copy Constructor *******************************************/

PIPDRTreeLikelihoodParameters::PIPDRTreeLikelihoodParameters(const PIPDRTreeLikelihoodParameters &data):
        tree_(data.tree_),
        tau_(data.tau_),
        lambda_(data.lambda_),
        mu_(data.mu_),
        nu_(data.nu_),
        nodeIotaData_(data.nodeIotaData_),
        nodeBetaData_(data.nodeBetaData_)
        {}


/*************************** Operator Overloading *******************************************/

PIPDRTreeLikelihoodParameters &PIPDRTreeLikelihoodParameters::operator=(const PIPDRTreeLikelihoodParameters &data) {
    tree_           = data.tree_;
    tau_            = data.tau_;
    lambda_         = data.lambda_;
    mu_             = data.mu_;
    nu_             = data.nu_;
    nodeIotaData_   = data.nodeIotaData_;
    nodeBetaData_   = data.nodeBetaData_;
    return *this;
}

/************************************* initPIPParameter_ *****************************************/

void PIPDRTreeLikelihoodParameters::initPIPParameter_() {
    /* Getting total tree length */
    tau_ = TreeTemplateTools::getTotalLength(*tree_->getRootNode(), false);

    if (lambda_ == std::string::npos || mu_ == std::string::npos) {
        //need to call estimator, but right now throw exception:
        throw Exception("PIPDRTreeLikelihoodParameters::initPIPParameter_(). Lambda or Mu is missing!");
    }

    ApplicationTools::displayTask("\nInitializing iota and beta parameters");

    // Initializing iota and beta parameters
    computePIPIota(tree_->getRootNode(), tau_, mu_, true);
    computePIPBeta(tree_->getRootNode(), mu_, true);
    DLOG(INFO) << "[PIP tree likelihood] Parameter iota and beta are initialized." << std::endl;
    ApplicationTools::displayTaskDone();

    // Initializing nu parameter
    computeNu_();
    DLOG(INFO) << "[PIP tree likelihood] The value of nu is: " << nu_;

}

/************************************** computeNu_ ****************************************/

void PIPDRTreeLikelihoodParameters::computeNu_(){


    if (fabs(mu_) < 1e-8)
        DLOG(WARNING) << "[PIP tree likelihood] Constraint match at parameter mu, badValue = " << mu_ << " [ 1e-08; 10000]";

    // Initializing nu parameter
    nu_ = lambda_ * (tau_ + 1 / mu_);

}

/************************************** computePIPIota ****************************************/

void PIPDRTreeLikelihoodParameters::computePIPIota(const Node *node, const double tau, const double mu, bool isRoot) {

    if (fabs(mu) < 1e-8)
        DLOG(WARNING) << "[PIP tree likelihood] Constraint match at parameter mu, badValue = " << mu << " [ 1e-08; 10000]";

    // Initialize iota vector:
    double iotaValue{0};
    double bl{0};
    Vdouble *iota_node_ = &nodeIotaData_;

    iota_node_->resize(tree_->getNumberOfNodes());

    if (isRoot) {
        isRoot = false;
        iotaValue = (1 / mu) / (tau + (1 / mu));
    } else {
        bl = node->getDistanceToFather();
        // This is one of the cases that user forget to make sure that the branch length of zero is not valid except for
        // the root. So we have to check it here, resolve it and warn the user.
        // To resolve it, we set the branch length to a very small value (1e-8).
        if (bl == 0) {

//            iotaValue = 0;
            int nId = node->getId();
            Node *theNode = const_cast<Node *>(tree_->getNode(nId));
            theNode->setDistanceToFather(1e-20);

            ApplicationTools::displayWarning("Branch length of node " + node->getName() +
                                             " is 0. It is not valid except for the root. now, we set branch length to very small value 1e-20.");

            DLOG(WARNING) << "[PIP tree likelihood] Branch length of node " << node->getName()
                          << " is 0. It is not valid except for the root. now, we set branch length to very small value 1e-20.";
        }
        iotaValue = bl / (tau + (1 / mu));
    }

    // setNodeIota
    (*iota_node_)[node->getId()] = iotaValue;

    if (!node->isLeaf()) {
        // we initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();
        for (size_t l = 0; l < nbSonNodes; l++) {
            // For each son node
            computePIPIota(node->getSon(l), tau, mu, isRoot);
        }
    }
}

/************************************** computePIPBeta ****************************************/

void PIPDRTreeLikelihoodParameters::computePIPBeta(const Node *node, const double mu, bool isRoot) {

    if (fabs(mu) < 1e-8)
        DLOG(WARNING) << "[PIP tree Llikelihood] Constraint match at parameter mu, badValue = " << mu << " [ 1e-08; 10000]";

    // Initialize Beta vector:
    double betaValue{0};
    double bl{0};
    Vdouble *beta_node_ = &nodeBetaData_;
    beta_node_->resize(tree_->getNumberOfNodes());

    if (isRoot) {
        isRoot = false;
        betaValue = 1.;
    } else {

        bl = node->getDistanceToFather();
        // see the comment in computePIPIota function
        if (bl== 0)
            DLOG(WARNING) << "[PIP tree likelihood] Branch length of node " << node->getName()
                          << " is 0. It is not valid except for the root. This should never happened as it is solved in iota.";

        betaValue = (1 - exp(-mu * bl)) / (mu * bl);
    }

    // setNodeBeta
    (*beta_node_)[node->getId()] = betaValue;

    if (!node->isLeaf()) {

        // we initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();
        for (size_t l = 0; l < nbSonNodes; l++)
        {
            // For each son node
            computePIPBeta(node->getSon(l), mu, isRoot);
        }
    }
}



/******************************************************************************/
/************************* PIPDRTreeLikelihoodData ****************************/
/******************************************************************************/
/*********************************** Constructor **********************************************/
PIPDRTreeLikelihoodData::PIPDRTreeLikelihoodData(const TreeTemplate<Node> *tree):
        AbstractTreeLikelihoodData(tree),
        pipParam_(),
        nodeData_(), leafData_(), rootLikelihoods_(),
        rootLikelihoodsS_(), rootLikelihoodsSR_(),
        shrunkData_(0), nbSites_(0), nbStates_(0), nbClasses_(0), nbDistinctSites_(0) {

}

/*********************************** Constructor **********************************************/
PIPDRTreeLikelihoodData::PIPDRTreeLikelihoodData(const TreeTemplate<Node> *tree, const double lambda, const double mu) :
        AbstractTreeLikelihoodData(tree),
        pipParam_(tree, lambda, mu),
        nodeData_(), leafData_(), rootLikelihoods_(),
        rootLikelihoodsS_(), rootLikelihoodsSR_(),
        shrunkData_(0), nbSites_(0), nbStates_(0), nbClasses_(0), nbDistinctSites_(0) {

}

/******************************** Copy Constructor *******************************************/

PIPDRTreeLikelihoodData::PIPDRTreeLikelihoodData(const PIPDRTreeLikelihoodData &data):
        AbstractTreeLikelihoodData(data),
        pipParam_(data.pipParam_),
        nodeData_(data.nodeData_), leafData_(data.leafData_),
        rootLikelihoods_(data.rootLikelihoods_),
        rootLikelihoodsS_(data.rootLikelihoodsS_),
        rootLikelihoodsSR_(data.rootLikelihoodsSR_),
        shrunkData_(0),
        nbSites_(data.nbSites_), nbStates_(data.nbStates_),
        nbClasses_(data.nbClasses_), nbDistinctSites_(data.nbDistinctSites_)
{
    if (data.shrunkData_)
        shrunkData_ = dynamic_cast<SiteContainer*>(data.shrunkData_->clone());
}

/*************************** Operator Overloading *******************************************/

PIPDRTreeLikelihoodData &PIPDRTreeLikelihoodData::operator=(const PIPDRTreeLikelihoodData &data)
{
    AbstractTreeLikelihoodData::operator=(data);
    pipParam_          = data.pipParam_;
    nodeData_          = data.nodeData_;
    leafData_          = data.leafData_;
    rootLikelihoods_   = data.rootLikelihoods_;
    rootLikelihoodsS_  = data.rootLikelihoodsS_;
    rootLikelihoodsSR_ = data.rootLikelihoodsSR_;
    nbSites_           = data.nbSites_;
    nbStates_          = data.nbStates_;
    nbClasses_         = data.nbClasses_;
    nbDistinctSites_   = data.nbDistinctSites_;

    if (shrunkData_) delete shrunkData_;
    if (data.shrunkData_)
        shrunkData_      = dynamic_cast<SiteContainer *>(data.shrunkData_->clone());
    else
        shrunkData_      = 0;
    return *this;
}

/*****************************************************************************************/

void PIPDRTreeLikelihoodData::setTree(const TreeTemplate<Node> *tree)
{
    tree_ = tree;
    for (std::map<int, PIPDRTreeLikelihoodNodeData>::iterator it = nodeData_.begin(); it != nodeData_.end(); it++)
    {
        int id = it->second.getNode()->getId();
        it->second.setNode(tree_->getNode(id));
    }
    for (std::map<int, PIPDRTreeLikelihoodLeafData>::iterator it = leafData_.begin(); it != leafData_.end(); it++)
    {
        int id = it->second.getNode()->getId();
        it->second.setNode(tree_->getNode(id));
    }
}

/***************************************************************************************/

void PIPDRTreeLikelihoodData::initLikelihoods(const SiteContainer &sites, const TransitionModel &model)
{
    if (sites.getNumberOfSequences() == 1)
        throw Exception("Error, only 1 sequence!");
    if (sites.getNumberOfSequences() == 0)
        throw Exception("Error, no sequence!");
    if (sites.getAlphabet()->getAlphabetType()
        != model.getAlphabet()->getAlphabetType())
        throw AlphabetMismatchException("PIPDRTreeLikelihoodData::initLikelihoods. Data and model must have the same alphabet type.",
                                        sites.getAlphabet(),
                                        model.getAlphabet());
    if(model.getName()!="PIP13")
        throw Exception("PIPDRTreeLikelihoodData::initLikelihoods. Error, You have to use PIP13 model!");

    alphabet_ = sites.getAlphabet();
    nbStates_ = model.getNumberOfStates();
    nbSites_  = sites.getNumberOfSites();

    SitePatterns pattern(&sites);
    if (shrunkData_)
        delete shrunkData_;
    shrunkData_       = pattern.getSites();
    rootWeights_      = pattern.getWeights();
    rootPatternLinks_ = pattern.getIndices();
    nbDistinctSites_  = shrunkData_->getNumberOfSites();

    const SiteContainer* sequences = new AlignedSequenceContainer(*shrunkData_);
    // Now initialize likelihood values and pointers:
    // Clone data for more efficiency on sequences access:
    initLikelihoods(tree_->getRootNode(), *sequences, model);


    delete sequences;

    // Now initialize root likelihoods with one: No derivatives in this version
    rootLikelihoods_.resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
        VVdouble* rootLikelihoods_i_ = &rootLikelihoods_[i];
        rootLikelihoods_i_->resize(nbClasses_);
        for (size_t c = 0; c < nbClasses_; c++)
        {
            Vdouble* rootLikelihoods_i_c_ = &(*rootLikelihoods_i_)[c];
            rootLikelihoods_i_c_->resize(nbStates_);
            for (size_t x = 0; x < nbStates_; x++)
            {
                (*rootLikelihoods_i_c_)[x] = 1.;
            }
        }
    }
}

/****************************************************************************************/

void PIPDRTreeLikelihoodData::initLikelihoods(const Node *node, const SiteContainer &sites, const TransitionModel &model)
{
    if (node->isLeaf())
    {
        // Init leaves likelihoods:
        const Sequence* seq;
        try
        {
            seq = &sites.getSequence(node->getName());
        }
        catch (SequenceNotFoundException& snfe)
        {
            throw SequenceNotFoundException("PIPDRTreeLikelihoodData::initlikelihoods. Leaf name in tree not found in site container: ", (node->getName()));
        }
        PIPDRTreeLikelihoodLeafData *leafData = &leafData_[node->getId()];
        VVdouble *leavesLikelihoods_leaf = &leafData->getLikelihoodArray();
        leafData->setNode(node);
        leavesLikelihoods_leaf->resize(nbDistinctSites_);
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
            Vdouble* leavesLikelihoods_leaf_i = &(*leavesLikelihoods_leaf)[i];
            leavesLikelihoods_leaf_i->resize(nbStates_);
            int state = seq->getValue(i);
            double test = 0.;

            for (size_t s = 0; s < nbStates_; s++)
            {
                // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
                // otherwise value set to 0:
                ( *leavesLikelihoods_leaf_i)[s] = model.getInitValue(s, state);
                test += ( *leavesLikelihoods_leaf_i)[s];
            }
            if (test < 0.000001)
                std::cerr << "WARNING!!! Likelihood will be 0 for site " << i << std::endl;
        }
    }

    // We initialize each son node first ('node' is an internal node):
    size_t nbSonNodes = node->getNumberOfSons();
    for (size_t l = 0; l < nbSonNodes; l++)
    {
        // For each son node,
        initLikelihoods(node->getSon(l), sites, model);
    }

    // Initialize likelihood vector:
    PIPDRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
    std::map<int, VVVdouble>* likelihoods_node_ = &nodeData->getLikelihoodArrays();
    nodeData->setNode(node);

    int nbSons = static_cast<int>(node->getNumberOfSons());

    for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
        const Node* neighbor = (*node)[n];
        VVVdouble* likelihoods_node_neighbor_ = &(*likelihoods_node_)[neighbor->getId()];

        likelihoods_node_neighbor_->resize(nbDistinctSites_);

        if (neighbor->isLeaf())
        {
            VVdouble* leavesLikelihoods_leaf_ = &leafData_[neighbor->getId()].getLikelihoodArray();
            for (size_t i = 0; i < nbDistinctSites_; i++)
            {
                Vdouble* leavesLikelihoods_leaf_i_ = &(*leavesLikelihoods_leaf_)[i];
                VVdouble* likelihoods_node_neighbor_i_ = &(*likelihoods_node_neighbor_)[i];
                likelihoods_node_neighbor_i_->resize(nbClasses_);
                for (size_t c = 0; c < nbClasses_; c++)
                {
                    Vdouble* likelihoods_node_neighbor_i_c_ = &(*likelihoods_node_neighbor_i_)[c];
                    likelihoods_node_neighbor_i_c_->resize(nbStates_);
                    for (size_t s = 0; s < nbStates_; s++)
                    {
                        (*likelihoods_node_neighbor_i_c_)[s] = (*leavesLikelihoods_leaf_i_)[s];
                    }
                }
            }
        }
        else
        {
            for (size_t i = 0; i < nbDistinctSites_; i++)
            {
                VVdouble* likelihoods_node_neighbor_i_ = &(*likelihoods_node_neighbor_)[i];
                likelihoods_node_neighbor_i_->resize(nbClasses_);
                for (size_t c = 0; c < nbClasses_; c++)
                {
                    Vdouble* likelihoods_node_neighbor_i_c_ = &(*likelihoods_node_neighbor_i_)[c];
                    likelihoods_node_neighbor_i_c_->resize(nbStates_);
                    for (size_t s = 0; s < nbStates_; s++)
                    {
                        (*likelihoods_node_neighbor_i_c_)[s] = 1.; // All likelihoods are initialized to 1.
                    }
                }
            }
        }
    }

    // Initialize d and d2 likelihoods:
    Vdouble* dLikelihoods_node_ = &nodeData->getDLikelihoodArray();
    Vdouble* d2Likelihoods_node_ = &nodeData->getD2LikelihoodArray();
    dLikelihoods_node_->resize(nbDistinctSites_);
    d2Likelihoods_node_->resize(nbDistinctSites_);
}

/***************************************************************************************/

void PIPDRTreeLikelihoodData::firParameterChanged(const double lambda, const double mu) {
    pipParam_.setLambda(lambda);
    pipParam_.setMu(mu);
    double tau = pipParam_.getTau();
    if(tau!=pipParam_.getTau())
        ApplicationTools::displayError("PIPDRTreeLikelihoodData::firParameterChanged. The tree is different.");

    // Recompute Nu_
    pipParam_.computeNu_();
    // Set all iotas and beta
    pipParam_.computePIPIota(tree_->getRootNode(),tau, mu, true);
    pipParam_.computePIPBeta(tree_->getRootNode(), mu, true);

}

/***************************************************************************************/

void PIPDRTreeLikelihoodData::reInit()
{
    reInit(tree_->getRootNode());
}

void PIPDRTreeLikelihoodData::reInit(const Node* node)
{
    if (node->isLeaf())
    {
        PIPDRTreeLikelihoodLeafData* leafData = &leafData_[node->getId()];
        leafData->setNode(node);
    }

    PIPDRTreeLikelihoodNodeData* nodeData = &nodeData_[node->getId()];
    nodeData->setNode(node);
    nodeData->eraseNeighborArrays();

    int nbSons = static_cast<int>(node->getNumberOfSons());

    for (int n = (node->hasFather() ? -1 : 0); n < nbSons; n++)
    {
        const Node* neighbor = (*node)[n];
        VVVdouble* array = &nodeData->getLikelihoodArrayForNeighbor(neighbor->getId());

        array->resize(nbDistinctSites_);
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
            VVdouble* array_i = &(*array)[i];
            array_i->resize(nbClasses_);
            for (size_t c = 0; c < nbClasses_; c++)
            {
                Vdouble* array_i_c = &(*array_i)[c];
                array_i_c->resize(nbStates_);
                for (size_t s = 0; s < nbStates_; s++)
                {
                    (*array_i_c)[s] = 1.; // All likelihoods are initialized to 1.
                }
            }
        }
    }

    // We re-initialize each son node:
    size_t nbSonNodes = node->getNumberOfSons();
    for (size_t l = 0; l < nbSonNodes; l++)
    {
        // For each son node,
        reInit(node->getSon(l));
    }

}
/*************************************setPIPTopologicalFlags**************************************************/

void PIPDRTreeLikelihoodData::setPIPTopologicalFlags(){
    // Topological flags:
    // Initializing Number Of Non-Gap Chars
    const std::vector<size_t> numberOfNonGapChar = initNumberOfNonGapChar();
    std::map<size_t, std::vector<int>> nodeSumCharacter;

    // Initializing set A members
    ApplicationTools::displayTask("Initializing set A for all sites");
    computePIPSetA(tree_->getRootNode(), numberOfNonGapChar, nodeSumCharacter);
    ApplicationTools::displayTaskDone();

    // Initializing set G members
    ApplicationTools::displayTask("Initializing set G for all sites");
    computePIPSetG(tree_->getRootNode());
    ApplicationTools::displayTaskDone();
}

/******************************************* computePIPSetAForSite ********************************************/

void PIPDRTreeLikelihoodData::computePIPSetAForSite(const Node *node, const Site &colSite,
                                                           const std::vector<size_t> numberOfNonGapChar,
                                                           std::map<size_t, std::vector<int>> &nodeSumCharacter, size_t nbSite){
    std::map<size_t, std::vector<bool>> *setA_node_ = &nodeSetA_;
    std::map<size_t, std::vector<int>> *sumChar_node_ = &nodeSumCharacter;
    int nodeId = node->getId();
    std::vector<bool> *setA_node_n_ = &(*setA_node_)[nodeId];
    std::vector<int> *sumChar_node_n_ = &(*sumChar_node_)[nodeId];
    setA_node_n_->resize(nbDistinctSites_);
    sumChar_node_n_->resize(nbDistinctSites_);

    if (node->isLeaf()) {
//        std::cout << std::endl << "Node " << node->getName() << " with number:" << node->getId();

        // Char state i is Gap:
        if (shrunkData_->getSequence(node->getName())[nbSite] == -1) {
//            if (colSite[nodeId] == -1) {
            (*sumChar_node_n_)[nbSite] = 0;
        } else {
            (*sumChar_node_n_)[nbSite] = 1.;
        }
        if ((*sumChar_node_n_)[nbSite] >= numberOfNonGapChar[nbSite]) {
            (*setA_node_n_)[nbSite] = 1.;
        } else {
            (*setA_node_n_)[nbSite] = 0;
        }
    } else {

        // we initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();
        for (size_t l = 0; l < nbSonNodes; l++) {
            // For each son node:
            computePIPSetAForSite(node->getSon(l), colSite, numberOfNonGapChar, nodeSumCharacter, nbSite);
        }
        //////////// By value !!!!!!!
        // To have access to node element:
        std::vector<int> sumChar_node_left_right_n_;
        sumChar_node_left_right_n_.resize(nbDistinctSites_);

        for (size_t l = 0; l < nbSonNodes; l++) {
            sumChar_node_left_right_n_[nbSite] += (nodeSumCharacter)[node->getSon(l)->getId()][nbSite];
        }
        // Assigning the value:
        (*sumChar_node_n_)[nbSite] = sumChar_node_left_right_n_[nbSite];
        if ((*sumChar_node_n_)[nbSite] >= numberOfNonGapChar[nbSite]) {
            (*setA_node_n_)[nbSite] = 1;
        } else {
            (*setA_node_n_)[nbSite] = 0;
        }
    }
}

/***************************************** initNumberOfNonGapChar *************************************/

std::vector<size_t> PIPDRTreeLikelihoodData::initNumberOfNonGapChar()
{
    std::vector<std::string> seqNames = shrunkData_->getSequencesNames();
    shrunkData_->getSequence(1).toString();
    std::vector<size_t> numChar;
    numChar.resize(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++) {
        Site columnSite = shrunkData_->getSite(i);
//        std::cout << "COL" << i << " is: " << columnSite.toString() ;

        for (size_t s = 0; s < columnSite.size(); s++) {
            // Char state is not Gap:
            if(columnSite[s]!=-1)
                numChar[i] += 1;
        }
//        std::cout << " and number of Non-gap char is: " << numChar[i] << std::endl;
    }
    if (std::find(numChar.begin(), numChar.end(), 0) !=numChar.end()){
        throw Exception("PIPDRTreeLikelihoodData::initNumberOfNonGapChar: Column full of gap is not allowed.");
    }
    return numChar;
}

/******************************************* computePIPSetGForSite ***********************************/

void PIPDRTreeLikelihoodData::computePIPSetGForSite(const Node *node, const Site &colSite, size_t siteNb) {


    std::map<size_t, std::vector<bool>> *node_all_gap_ = &nodeSetG_;
    int nodeId = node->getId();
    std::vector<bool> *node_all_gap_n_ = &(*node_all_gap_)[nodeId];
    node_all_gap_n_->resize(nbDistinctSites_);

    if (node->isLeaf()) {
//        std::cout << "Node " << node->getName() << " with number:" << node->getId();
//        std::cout << " and the " << siteNb << "th element is: " << data_->getSequence(node->getName()).toString()[siteNb]
//                  << std::endl;

        // Char state i is Gap:
        if (shrunkData_->getSequence(node->getName())[siteNb] == -1) {
            (*node_all_gap_n_)[siteNb] = 1.;
        } else {
            (*node_all_gap_n_)[siteNb] = 0.;
        }

    } else {
        // We initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();

        for (int l = 0; l < nbSonNodes; ++l) {
            computePIPSetGForSite(node->getSon(l), colSite, siteNb);
        }
        // To have access to node elements:
        if (((*node_all_gap_)[node->getSonsId()[0]])[siteNb] && ((*node_all_gap_)[node->getSonsId()[1]])[siteNb]) {
            (*node_all_gap_n_)[siteNb] = 1;
        }else{
            (*node_all_gap_n_)[siteNb] = 0;
        }
    }
}

/***************************************************************************************/



