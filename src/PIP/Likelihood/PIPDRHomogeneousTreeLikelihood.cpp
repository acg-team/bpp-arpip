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
 * @file PIPDRHomogeneousTreeLikelihood.cpp
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
#include "PIPDRHomogeneousTreeLikelihood.h"
#include "PIPDRTreeLikelihoodData.h"

// From PhyLib
#include <Bpp/Phyl/PatternTools.h>

// From SeqLib
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainer.h>

// From CoreLib
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From STL
#include <iostream>

using namespace std;

/********************************** Constructor ********************************************/

PIPDRHomogeneousTreeLikelihood::PIPDRHomogeneousTreeLikelihood(
        const Tree &tree,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, false, verbose),
        // checkRooted is override with "false" to prevent a twice deletion of last node.
        likelihoodData_(0), shrunkData_(0), isPIP_(1), p_c_0_(0), p_c_(0), phi_(0),
        nodeFTildeValues_(), nodeFValues_(), nodeSumFTilde_(),
        minusLogLik_(-1.)
{
    init_();
}


/******************************** Constructor **********************************************/

PIPDRHomogeneousTreeLikelihood::PIPDRHomogeneousTreeLikelihood(
        const Tree &tree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        bool checkRooted,
        bool verbose) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, false, verbose),
        // checkRooted is override with "false" to prevent a twice deletion of last node.
        likelihoodData_(0), shrunkData_(0), isPIP_(1), p_c_0_(0), p_c_(0), phi_(0),
        nodeFTildeValues_(), nodeFValues_(), nodeSumFTilde_(),
        minusLogLik_(-1.) {
    init_();
    setData(data);
    initialize();// from AbstractHomogeneousTreeLikelihood to fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes by computeAllTransitionProbabilities();

}

/******************************** Constructor Full PIP version **********************************************/

PIPDRHomogeneousTreeLikelihood::PIPDRHomogeneousTreeLikelihood(
        const Tree &tree,
        const SiteContainer &data,
        TransitionModel *model,
        DiscreteDistribution *rDist,
        const double lambda, const double mu,
        bool checkRooted,
        bool verbose) :
        AbstractHomogeneousTreeLikelihood(tree, model, rDist, false, verbose),
        // checkRooted is override with "false" to prevent a twice deletion of last node.
        likelihoodData_(0), shrunkData_(0), isPIP_(1), p_c_0_(0), p_c_(0), phi_(0),
        nodeFTildeValues_(), nodeFValues_(), nodeSumFTilde_(),
        minusLogLik_(-1.) {
    init_(lambda, mu);
    setData(data);
    initialize();// from AbstractHomogeneousTreeLikelihood to fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes by computeAllTransitionProbabilities();
//    computePIPTreeLikelihood(lambda, mu);
}

/***************************** Copy Constructor ************************************************/

PIPDRHomogeneousTreeLikelihood::PIPDRHomogeneousTreeLikelihood(const PIPDRHomogeneousTreeLikelihood &lik) :
        AbstractHomogeneousTreeLikelihood(lik),
        likelihoodData_(0), shrunkData_(0),
        nodeFTildeValues_(lik.nodeFTildeValues_),
        nodeFValues_(lik.nodeFValues_),
        nodeSumFTilde_(lik.nodeSumFTilde_),
        minusLogLik_(-1.)
{
    likelihoodData_ = dynamic_cast<PIPDRTreeLikelihoodData *>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);

    if (shrunkData_)
        delete shrunkData_;
    shrunkData_       = PatternTools::shrinkSiteSet(*data_);

    minusLogLik_ = lik.minusLogLik_;

}

/******************************* Operator Overloading ***********************************************/

PIPDRHomogeneousTreeLikelihood& PIPDRHomogeneousTreeLikelihood::operator=(const PIPDRHomogeneousTreeLikelihood& lik)
{
    AbstractHomogeneousTreeLikelihood::operator=(lik);
    if (likelihoodData_)
        delete likelihoodData_;
    likelihoodData_     = dynamic_cast<PIPDRTreeLikelihoodData*>(lik.likelihoodData_->clone());
    likelihoodData_->setTree(tree_);
    shrunkData_         = lik.shrunkData_;
    isPIP_              = lik.isPIP_;
    p_c_0_              = lik.p_c_0_;
    p_c_                = lik.p_c_;
    phi_                = lik.phi_;
    nodeFTildeValues_   = lik.nodeFTildeValues_;
    nodeFValues_        = lik.nodeFValues_;
    nodeSumFTilde_      = lik.nodeSumFTilde_;
    minusLogLik_        = lik.minusLogLik_;

    return *this;
}

/************************************* init_ *****************************************/

void PIPDRHomogeneousTreeLikelihood::init_() {

    //Check again if the tree is rooted (Tree should be rooted in PIP model):
    if (!tree_->isRooted()) throw Exception("PIPDRHomogeneousTreeLikelihood::init_(). Tree is unrooted!");
    if (TreeTemplateTools::isMultifurcating(*tree_->getRootNode())) throw Exception("PIPDRHomogeneousTreeLikelihood::init_(). Tree is multifurcation.");
    setMinimumBranchLength(0.);

    // init data structure ( Allocate the DS)
    likelihoodData_ = new PIPDRTreeLikelihoodData(tree_);

    computeFirstOrderDerivatives_ = false;
    computeSecondOrderDerivatives_ = false;
}
/************************************* init_ *****************************************/

void PIPDRHomogeneousTreeLikelihood::init_(const double lambda, const double mu) {

    //Check again if the tree is rooted (Tree should be rooted in PIP model):
    if (!tree_->isRooted()) throw Exception("PIPDRHomogeneousTreeLikelihood::init_(). Tree is unrooted!");
    if (TreeTemplateTools::isMultifurcating(*tree_->getRootNode())) throw Exception("PIPDRHomogeneousTreeLikelihood::init_(). Tree is multifurcation.");
    setMinimumBranchLength(0.);

    // init data structure ( Allocate the DS)
    likelihoodData_ = new PIPDRTreeLikelihoodData(
            tree_,
            lambda, mu);
    //Not implemented yet.  todo: implementing the derivatives.
    computeFirstOrderDerivatives_ = false;
    computeSecondOrderDerivatives_ = false;

}
/***************************************************************************************/
/******************************* PIP Specific Functions ********************************/
/***************************************************************************************/

/************************** computePIPLikelihoodEmptySite ******************************/
void PIPDRHomogeneousTreeLikelihood::computePIPLikelihoodForEmptySiteAtNode_(const Node *node){
    int nodeId = node->getId();
    /* In this scenario we only have one imaginary column full of gap
     * So, we just define the local variable with similar name.
     * The only thing which we store here is just the p_c_0!
     */
//    std::map<int ,Vdouble> fTildeValue_node;
    VVdouble *fTildeValue_node_n = &nodeFTildeValues_[nodeId];
    Vdouble *fTildeValue_node_n_i = &(*fTildeValue_node_n)[0];
    fTildeValue_node_n_i->resize(nbStates_);


//    std::map<int ,double> fValue_node;
    double *fValue_node_n_i = &nodeFValues_[nodeId][0];


//    Vdouble sumFTilde_node;
//    sumFTilde_node.resize(nbNodes_);
    double *sumFTilde_node_n_i = &nodeSumFTilde_[nodeId][0];


    const double *iota_node_n = &(likelihoodData_->getPipParam().getNodeIotaData(nodeId));
    const double *beta_node_n = &(likelihoodData_->getPipParam().getNodeBetaData(nodeId));


    if (node->isLeaf()) {
        // Init leaves likelihoods:
        for (int s = 0; s < nbStates_-1; ++s) {
            (*fTildeValue_node_n_i)[s] = 0;
        }
        (*fTildeValue_node_n_i)[nbStates_-1] = 1;

    } else {
        size_t nbSonNodes = node->getNumberOfSons();
        // Recursive call:
        for (int l = 0; l < nbSonNodes; ++l) {
            computePIPLikelihoodForEmptySiteAtNode_(node->getSon(l));
        }
        fTildeValue_node_n_i->assign(nbStates_, 1);
        //Computation:
        vector<const Vdouble *> iLik;
        vector<const VVdouble *> tProb;

        for (size_t n = 0; n < nbSonNodes; n++) {
            const Node *son = node->getSon(n);
            tProb.push_back(&pxy_[son->getId()][0]);// all zero AbstractHomogeneousTreeLikelihood-> solved in PIPDRHomogeneousTreeLikelihood::init
            iLik.push_back(&nodeFTildeValues_[son->getId()][0]);//
        }

        for (size_t n = 0; n < nbSonNodes; n++) {
            const VVdouble *pxy_n = tProb[n];
            const Vdouble *iLik_n = iLik[n];
//            const Vdouble *iLik_n_i = &(*iLik_n)[0];
            for (int x = 0; x < nbStates_; ++x) {
                // For each internal state
                double tmp_fv = 0;
                const Vdouble *pxy_n_x = &(*pxy_n)[x];
//                const double *iLik_n_x = &(*iLik_n)[x];
                for (size_t y = 0; y < nbStates_; ++y) {
//                    cout << "1." << (*pxy_n_x)[y] << endl;
//                    cout << "2." << (*iLik_n)[y] << endl;
                    tmp_fv += (*pxy_n_x)[y] * (*iLik_n)[y];
//                    cout  << y << x << "\t" << endl;
                }
                // We store this conditional likelihood into the corresponding array
                (*fTildeValue_node_n_i)[x] *= tmp_fv;
            }
        }
    }

    // Computing the Likelihood value: aka Felsenstein's value
    double tmpSumFV = 0;
    for (int x = 0; x < nbStates_; ++x) {
        double tmpValue = 0;
        tmpValue = rootFreqs_[x] * (*fTildeValue_node_n_i)[x];
        tmpSumFV += tmpValue;
    }
    *sumFTilde_node_n_i = tmpSumFV;


    if (tree_->getRootNode() == node) {
        (*fValue_node_n_i) = *sumFTilde_node_n_i;
    } else {
        (*fValue_node_n_i) = 1 + (*beta_node_n * (*sumFTilde_node_n_i - 1));
    }
    p_c_0_ += (*fValue_node_n_i) * (*iota_node_n);
}

/*********************************** computePIPPhi_ ***********************************/

void PIPDRHomogeneousTreeLikelihood::computePIPPhi_(){
    if (nbSites_==0){
        throw Exception("PIPDRHomogeneousTreeLikelihood::computePIPPhi_: The number of sites is zero!");
        LOG(FATAL)<<"[PIPDRHomogeneousTreeLikelihood] computePIPPhi_: The number of sites is zero.";
    }
    if(p_c_0_ ==0 || isnan(p_c_0_)){
        throw Exception(
                "PIPDRHomogeneousTreeLikelihood::computePIPPhi_: The likelihood value of empty columns is zero or nan.");
        LOG(FATAL)<<"[PIPDRHomogeneousTreeLikelihood] computePIPPhi_: The likelihood value of empty columns is zero or nan,"
                    "In the case of nan p_c_0_: plz check the branch lengths. Branches of length zero is not valid except for the root.";
    }
//            phi_ = (pow(nu_, nbSites_) * exp(((p_c_0_) - 1) * nu_)) /std::tgamma(nbSites_ + 1) ;
    double log_factorial_m;

    log_factorial_m = 0;
    for (int i = 1; i <= nbSites_; i++) {
        log_factorial_m += log(i);
    }

    double nu = likelihoodData_->getPipParam().getNu();

    phi_ = -log_factorial_m + nbSites_ * log(nu) + (nu * (p_c_0_ - 1));

}

/************************* computePIPLikelihoodNoneEmptySites **************************/

void PIPDRHomogeneousTreeLikelihood::computePIPLikelihoodForNonEmptySitesAtNode_(const Node *node) {
    // This section is independent from likelihood data structure defined in the bpp lib!!!

    int nodeId = node->getId();
    VVdouble *fTildeValue_node = &nodeFTildeValues_[nodeId];
    fTildeValue_node->resize(nbDistinctSites_);

    Vdouble *fValue_node = &nodeFValues_[nodeId];
    fValue_node->resize(nbDistinctSites_);

    Vdouble *sumFTilde_node = &nodeSumFTilde_[nodeId];
    sumFTilde_node->resize(nbDistinctSites_);



    const std::vector<bool> *setA_node = &(likelihoodData_->getNodeSetA(nodeId));
    const double *beta_node_n = &(likelihoodData_->getPipParam().getNodeBetaData(nodeId));
    const double *iota_node_n = &(likelihoodData_->getPipParam().getNodeIotaData(nodeId));

    //initialize the likelihood array
    if (node->isLeaf()) {

        // Init leaves likelihoods:
        const Sequence* seq;
        try
        {
            seq = &shrunkData_->getSequence(node->getName());
        }
        catch (SequenceNotFoundException& snfe)
        {
            LOG(FATAL)<<"[PIP homogeneous tree likelihood]"<< snfe.message();
            throw SequenceNotFoundException(
                    "PIPDRHomogeneousTreeLikelihood::computePIPLikelihoodForNonEmptySitesAtNode_. Leaf name in tree not found in site container: ",
                    (node->getName()));
        }

        for (size_t i = 0; i < nbDistinctSites_; i++) {
            Vdouble *fValueTilde_node_i = &(*fTildeValue_node)[i];
            fValueTilde_node_i->resize(nbStates_);
            int st_seq = seq->getValue(i);
//            cout << st_seq << endl;
            for (size_t s = 0; s < nbStates_; s++)
            {
                // Leaves likelihood are set to 1 if the char correspond to the site in the sequence,
                // otherwise value set to 0:
//                    ( *fValueTilde_node_i)[s] = model_->getInitValue(s, st_seq);
                (*fValueTilde_node_i)[s] = (s == st_seq)? (1):(0);
            }
            //gap as additional character:
            (*fValueTilde_node_i)[nbStates_-1] = (st_seq==-1)? (1):(0);
        }

    } else {

        // We initialize each son node:
        size_t nbSonNodes = node->getNumberOfSons();
        if (nbSonNodes == 0) {
            LOG(FATAL)<<"[PIP homogeneous tree likelihood] Post-order recursive call doesn't work.";
            throw Exception(
                    "PIPDRHomogeneousTreeLikelihood::computePIPLikelihoodForNonEmptySitesAtNode_: "
                    "Post-order recursive call doesn't work.");
        }
        // Recursive call:
        for (int l = 0; l < nbSonNodes; ++l) {
            computePIPLikelihoodForNonEmptySitesAtNode_(node->getSon(l));
        }

        // Set all likelihood array values to 1:
        for (size_t i = 0; i < nbDistinctSites_; i++) {
            Vdouble *fValueTilde_node_i = &(*fTildeValue_node)[i];
            fValueTilde_node_i->resize(nbStates_);
            fValueTilde_node_i->assign(nbStates_, 1);
        }

        //Computation:
        vector<const VVdouble *> iLik;
        vector<const VVdouble *> tProb;

        for (size_t n = 0; n < nbSonNodes; n++) {
            const Node *son = node->getSon(n);
            tProb.push_back(&pxy_[son->getId()][0]);// all zero AbstractHomogeneousTreeLikelihood-> solved in PIPDRHomogeneousTreeLikelihood::init
            iLik.push_back(&(nodeFTildeValues_)[son->getId()]);//
        }

        for (size_t n = 0; n < nbSonNodes; n++) {
            const VVdouble *pxy_n = tProb[n];
            const VVdouble *iLik_n = iLik[n];
            for (size_t i = 0; i < nbDistinctSites_; i++) {
                Vdouble *fValueTilde_node_i = &(*fTildeValue_node)[i];
                const Vdouble *iLik_n_i = &(*iLik_n)[i];
//                fValueTilde_node_i->resize(nbStates_);
//                fValueTilde_node_i->assign(nbStates_, 1);
                for (int x = 0; x < nbStates_; ++x) {
                    // For each internal state
                    double tmp_fv = 0;
                    const Vdouble *pxy_n_x = &(*pxy_n)[x];
                    for (size_t y = 0; y < nbStates_; ++y) {
//                        cout << "1." << (*pxy_n_x)[y] << endl;
//                        cout << "2." << (*iLik_n_i)[y] << endl;
                        tmp_fv += (*pxy_n_x)[y] * (*iLik_n_i)[y];
//                        cout << i << "\t" << y << x << "\t" << endl;
                    }
                    // We store this conditional likelihood into the corresponding array
                    (*fValueTilde_node_i)[x] *= tmp_fv;
                }
            }
        }
    }

    // Computing the intermediate variable for Likelihood value aka Felsenstein's value
    for (int i = 0; i < nbDistinctSites_; ++i) {
        double *sumFTilde_node_i = &(*sumFTilde_node)[i];
        Vdouble *fTildeValue_node_i = &(*fTildeValue_node)[i];
        double tmpSumFV = 0;
        for (int x = 0; x < nbStates_; ++x) {
            double tmpValue = 0;
            tmpValue = rootFreqs_[x] * (*fTildeValue_node_i)[x];
            tmpSumFV += tmpValue;
        }
        *sumFTilde_node_i = tmpSumFV;
    }

    // Computing the non-empty likelihood p(c):
    for (int i = 0; i < nbDistinctSites_; ++i) {

        double *sumFTilde_node_i = &(*sumFTilde_node)[i];
        double *fValue_node_i = &(*fValue_node)[i];
        double *p_c_i = &p_c_[i];

        if (tree_->getRootNode() == node) {
            (*fValue_node_i) = *sumFTilde_node_i;
        } else if ((*setA_node)[i]) {
            (*fValue_node_i) = (*sumFTilde_node_i) * (*beta_node_n);
        } else {
            (*fValue_node_i) = 0;
        }
        (*p_c_i) += (*fValue_node_i) * (*iota_node_n);
    }
}

/*********************************** resetFVArrays_ *************************************/

void PIPDRHomogeneousTreeLikelihood::resetFVArrays_(const Node *node)
{
    for (size_t n = 0; n < node->getNumberOfSons(); n++)
    {
        const Node *subNode = node->getSon(n);
        resetFVArrays_(subNode);
    }
    int nodeId = node->getId();
    VVdouble *fTildeValues_node_n = &nodeFTildeValues_[nodeId];
    Vdouble *sumFtilde_n = &nodeSumFTilde_[nodeId];
    sumFtilde_n->assign(nbDistinctSites_, 1);

    Vdouble *fValues_node_n = &nodeFValues_[nodeId];
    fValues_node_n->assign(nbStates_, 1);

    for (int i = 0; i < nbDistinctSites_; ++i) {
        Vdouble *fTildeValues_node_n_i = &(*fTildeValues_node_n)[i];
        fTildeValues_node_n_i->assign(nbStates_, 1);
    }

}

/*********************************** computePIPTreeLikelihood *************************************/

long double PIPDRHomogeneousTreeLikelihood::computePIPTreeLikelihood(const double lambda, const double mu) {

    // Initialise lk
    long double logLK;

    // Initialise vector of likelihood per each site
    std::vector<double> lk_sites(nbDistinctSites_);

    // Reset likelihood array:
    p_c_.assign(nbDistinctSites_, 0);
    p_c_0_ = 0;
//    resetFVArrays_(tree_->getRootNode());


    // Update the PIP parameters
    firePIPParameterChanged_(lambda, mu);
    DLOG(INFO)<< "[DR homogeneous tree likelihood] PIP's parameters (including transition probabilities) is now updated.";

    // For non-empty sites and all nodes in the MSA
    computePIPLikelihoodForNonEmptySitesAtNode_(tree_->getRootNode());

    // Reset Fv values since they are used again for empty sites
    resetFVArrays_(tree_->getRootNode());

    // For empty site and all nodes in the msa, compute the likelihood of empty column
    computePIPLikelihoodForEmptySiteAtNode_(tree_->getRootNode());

    // Compute PHi as PIP parameter
    computePIPPhi_();

    // We need to map distinct sites to normal site:
    // The number of times each unique site was found.
    const std::vector<unsigned int> *siteWeight = &likelihoodData_->getWeights();
    for (int i = 0; i < nbDistinctSites_; i++) {

        lk_sites[i] = log(p_c_[i]) * siteWeight->at(i);
        DVLOG(2) << "[DR homogeneous tree likelihood] site log_lk[" << i << "]=" << std::setprecision(18) << lk_sites[i] << std::endl;

    }

    // Sum all the values stored in the lk vector
    logLK = std::accumulate(lk_sites.begin(), lk_sites.end(), decltype(lk_sites)::value_type(0));
    if (isinf(logLK)) {
        throw Exception(
                "PIPDRHomogeneousTreeLikelihood::computePIPTreeLikelihood: Log likelihood value is not valid.");
        DLOG(INFO) << "[DR homogeneous tree likelihood] log likelihood value is " << logLK << " and not valid.";
    }
    DVLOG(2) << "LK Sites [BPP] " << std::setprecision(18) << logLK;

    // Add log phi to site likelihood
    logLK += phi_;

    // Calls the routine to compute the FV values
    minusLogLik_ = -logLK;
    DLOG(INFO) << "[DR homogeneous tree likelihood] PIP likelihood is computed now with value " << logLK;

    return logLK;

}

/*********************************** firePIPParameterChanged_ *************************************/

void PIPDRHomogeneousTreeLikelihood::firePIPParameterChanged_(const double lambda, const double mu) {

    if (lambda == std::string::npos || mu == std::string::npos) {
        //need to call estimator, but right now throw exception:
        LOG(FATAL)<<"[PIP homogeneous tree likelihood] Lambda or Mu is missing!";
        throw Exception("PIPDRHomogeneousTreeLikelihood::initPIPParameter_(). Lambda or Mu is missing!");
    }

    likelihoodData_->firParameterChanged(lambda, mu);
    
//    ParameterList parModel;
    if (!hasParameter("mu"))
        addParameter_(new Parameter("mu", mu));
    else {
        std::string parName = "mu";
        deleteParameter_(parName);
        addParameter_(new Parameter("mu", mu));
    }

    model_->setParametersValues(getParameters());

    // indel parameter changed, need to recompute all probs:
    computeAllTransitionProbabilities();
}

/*********************************** computeTreeLikelihood *************************************/

void PIPDRHomogeneousTreeLikelihood::computeTreeLikelihood()
{
    ApplicationTools::displayMessage("Not implemented: use computePIPTreeLikelihood instead.");
}

/***************************************************************************************/
/************************* Maximum Likelihood Indel Points *****************************/
/***************************************************************************************/




/********************************************************************************************************/
/******** The TreeLikelihood interface (other methods by AbstractTreeLikelihood class)*******************/
/************************************** setData *********************************************************/

void PIPDRHomogeneousTreeLikelihood::setData(const SiteContainer& sites) {
    if(data_)
        delete data_;
    data_ = PatternTools::getSequenceSubset(sites, *tree_->getRootNode());

    if(verbose_)
        ApplicationTools::displayTask("Initializing data structure");
    likelihoodData_->initLikelihoods(*data_, *model_);

    if(verbose_)
        ApplicationTools::displayTaskDone();

    if(isPIP_)// compute the SetA and SetG if true
        likelihoodData_->setPIPTopologicalFlags();

    // We need shrunk data:
    shrunkData_ = PatternTools::shrinkSiteSet(*data_);

//    // Fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes
//    computeAllTransitionProbabilities();

    nbSites_ = likelihoodData_->getNumberOfSites();
    nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
    nbStates_ = likelihoodData_->getNumberOfStates();

    if(verbose_)
        ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));
    initialized_ = false;
}

/******************************************************************************/

double PIPDRHomogeneousTreeLikelihood::getLikelihood() const {
    // Initialize lk
    double lik = 1.;
    ApplicationTools::displayMessage("The likelihood value is very small, try getLogLikelihood instead.");
    // The likelihood value is very small when we consider p_c_0

    // We need to map distinct sites to normal site:
    // The number of times each unique site was found.
    const std::vector<unsigned int> *w = &likelihoodData_->getWeights();
    for (int i = 0; i < nbDistinctSites_; i++) {
        lik *= std::pow(p_c_[i], w->at(i));
    }
    DVLOG(2) << "Likelihood value without empty columns is: "<< lik << std::endl;
    ApplicationTools::displayResult("The PIP Likelihood value without considering empty columns.", lik);

    return lik;
}

/******************************************************************************/

double PIPDRHomogeneousTreeLikelihood::getLogLikelihood() const {

    if(minusLogLik_ != -1)
        return -minusLogLik_;

    double logLik = 0;
    // We need to map distinct sites to normal site:
    // The number of times each unique site was found.
    const std::vector<unsigned int> *rootWeights = &likelihoodData_->getWeights();
    for (int i = 0; i < nbDistinctSites_; i++) {
        logLik += log(p_c_[i]) * rootWeights->at(i);
    }

    // Add logs phi to site likelihood
    logLik += phi_;


    return logLik;
}

/******************************************************************************/

double PIPDRHomogeneousTreeLikelihood::getLikelihoodForASite(size_t site) const
{
    ApplicationTools::displayMessage("The likelihood value is very small, try getLogLikelihoodForASite instead.");
    // The likelihood value is very small when we take account empty columns.

    // Extracting the position of the site in the shrunkData
    size_t rootPosition = likelihoodData_->getRootArrayPosition(site);

    p_c_[rootPosition];
    DVLOG(2) << "site lk[" << site << "]=" << std::setprecision(18) << p_c_[rootPosition] << std::endl;

    return p_c_[rootPosition];
}

/******************************************************************************/

double PIPDRHomogeneousTreeLikelihood::getLogLikelihoodForASite(size_t site) const
{
    if (!isInitialized())
        throw Exception("PIPDRHomogeneousTreeLikelihood::getLogLikelihoodForASite(). Instance is not initialized.");
    // Initialize lk
    double logLik= 1.;

    // Extracting the position of the site in the shrunkData
    size_t rootPosition = likelihoodData_->getRootArrayPosition(site);


    logLik = log(p_c_[rootPosition]);
    DVLOG(2) << "site log_lk[" << site << "]=" << std::setprecision(18) << logLik << std::endl;

    return logLik;
}




/*********************************** Not working yet!********************************************/
/********************** The DiscreteRatesAcrossSites interface implementation ********************/
/***********************************************************************************************/
double PIPDRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
    ApplicationTools::displayMessage(
            "getLikelihoodForASiteForARateClass for site " + to_string(site) + " and " + to_string(rateClass) +
            " is not implemented.");

    return 1;
}
/******************************************************************************/
double PIPDRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const
{
    ApplicationTools::displayMessage(
            "getLogLikelihoodForASiteForARateClass for site " + to_string(site) + " and " + to_string(rateClass) +
            " is not implemented.");
    return 1;
}
/******************************************************************************/
double PIPDRHomogeneousTreeLikelihood::getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
    ApplicationTools::displayMessage(
            "getLikelihoodForASiteForARateClassForAState for site " + to_string(site) + " and " + to_string(rateClass) +
            " is not implemented.");
    return 1;
}
/******************************************************************************/
double PIPDRHomogeneousTreeLikelihood::getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const
{
    ApplicationTools::displayMessage(
            "getLogLikelihoodForASiteForARateClassForAState for site " + to_string(site) + " and " + to_string(rateClass) +
            " is not implemented.");
    return 1;
}
/******************************************************************************/
void PIPDRHomogeneousTreeLikelihood::setParameters(const ParameterList &parameters)
{
    setParametersValues(parameters);
}
/******************************************************************************/
void PIPDRHomogeneousTreeLikelihood::fireParameterChanged(const ParameterList &params)
{
    applyParameters();

    if (rateDistribution_->getParameters().getCommonParametersWith(params).size() > 0
        || model_->getParameters().getCommonParametersWith(params).size() > 0)
    {
        // Rate parameter changed, need to recompute all probs:
        computeAllTransitionProbabilities();
    }
    else if (params.size() > 0)
    {
        // We may save some computations:
        for (size_t i = 0; i < params.size(); i++)
        {
            string s = params[i].getName();
            if (s.substr(0, 5) == "BrLen")
            {
                // Branch length parameter:
                computeTransitionProbabilitiesForNode(nodes_[TextTools::to < size_t > (s.substr(5))]);
            }
        }
    }

    // Check if we are using the pip model or not.
    if(isPIP_) {
        const double tmpLambda = likelihoodData_->getPipParam().getLambda();
        const double tmpMu = likelihoodData_->getPipParam().getMu();

        // To handle when one of the lambda or mu is zero!
        if (tmpLambda == 0 || tmpMu == 0) {
            computePIPTreeLikelihood();// call the default values.
        } else {
            computePIPTreeLikelihood(likelihoodData_->getPipParam().getLambda(),
                                     likelihoodData_->getPipParam().getMu());
        }
    } else{
            computeTreeLikelihood();
    }

    // todo: implementing these two methods.
    if (computeFirstOrderDerivatives_)
    {
        computeTreeDLikelihoods();
    }
    if (computeSecondOrderDerivatives_)
    {
        computeTreeD2Likelihoods();
    }

    minusLogLik_ = -getLogLikelihood();
}
/******************************************************************************/
double PIPDRHomogeneousTreeLikelihood::getValue() const {
    if (!isInitialized()) throw Exception("RHomogeneousTreeLikelihood::getValue(). Instance is not initialized.");
    return minusLogLik_;
}
/******************************************************************************/
void PIPDRHomogeneousTreeLikelihood::displayLikelihood(const Node* node)
{
    DVLOG(2) << "Likelihoods at node " << node->getName() << ": ";
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getFatherId(), node->getId()));
    DVLOG(2) << "                                         ***";
}

/******************************************************************************/

void PIPDRHomogeneousTreeLikelihood::resetLikelihoodArrays(const Node* node)
{
    for (size_t n = 0; n < node->getNumberOfSons(); n++)
    {
        const Node* subNode = node->getSon(n);
        resetLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), subNode->getId()));
    }
    if (node->hasFather())
    {
        const Node* father = node->getFather();
        resetLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), father->getId()));
    }
}



/***************************************************************************************/
/************************* ASR using Pupko's DP algorithm with indel **********************/
/***************************************************************************************/

/******************************************************************************/

void PIPDRHomogeneousTreeLikelihood::computeLikelihoodAtSite_(const Node *node,
                                                                   VVdouble &likelihoodArray,
                                                                   VVdouble &characterArray, std::size_t siteNumber,
                                                                   const int newRootId, const Node *sonNode) const {

    int nodeId = node->getId();
    likelihoodArray.resize(nbNodes_+1);
    characterArray.resize(nbNodes_+1);

    map<int, VVVdouble> *likelihoods_node = &likelihoodData_->getLikelihoodArrays(nodeId);

    Vdouble *pupko_candidate_ancestor_node = &characterArray[nodeId];
    Vdouble *pupko_likelihood_node = &likelihoodArray[nodeId];

    pupko_candidate_ancestor_node->resize(nbStates_);
    pupko_likelihood_node->resize(nbStates_);

    // Initialize likelihood array:
    if (node->isLeaf() && tree_->getRootId()!=nodeId) {
        vector<int> InnerNodes = tree_->getInnerNodesId();
//        VVdouble *leavesLikelihoods_node;
        if(std::find(InnerNodes.begin(), InnerNodes.end(), nodeId) != InnerNodes.end()) {
            //1.Set all likelihoods to 1 and Clikelihood -1 for a start:
//            for (size_t x = 0; x < nbStates_; x++) {
//                (*pupko_likelihood_node)[x] = 1.;
//                (*pupko_candidate_ancestor_node)[x] = -1;
//            }

            //2. Gap as additional char to backpropagate
            Vdouble tmpLikelihood;
            tmpLikelihood.resize(nbStates_);
            for (int tmpi = 0; tmpi < nbStates_; ++tmpi) {
                tmpLikelihood[tmpi] = 0;
            }
            tmpLikelihood[nbStates_-1] = 1;
            VVVdouble *pxy_n_i = &pxy_[nodeId];


            pupko_candidate_ancestor_node->resize(nbStates_);
            pupko_likelihood_node->resize(nbStates_);
            VVdouble *pxy_n_c = &(*pxy_n_i)[0];// in this simple case we only have one class
            size_t jj = VectorTools::whichMax(tmpLikelihood);
            for (size_t x = 0; x < nbStates_; x++) {
                (*pupko_likelihood_node)[x] = (*pxy_n_c)[x][jj];
                (*pupko_candidate_ancestor_node)[x] = VectorTools::whichMax(tmpLikelihood);
            }

        } else {
            VVdouble *leavesLikelihoods_node = &likelihoodData_->getLeafLikelihoods(nodeId);
            VVVdouble *pxy_n_i = &pxy_[nodeId];

            Vdouble *leavesLikelihoods_node_i = &(*leavesLikelihoods_node)[siteNumber];
            pupko_candidate_ancestor_node->resize(nbStates_);
            pupko_likelihood_node->resize(nbStates_);
            VVdouble *pxy_n_c = &(*pxy_n_i)[0];// in this simple case we only have one class
            size_t jj = VectorTools::whichMax(*leavesLikelihoods_node_i);
            for (size_t x = 0; x < nbStates_; x++) {
                (*pupko_likelihood_node)[x] = (*pxy_n_c)[x][jj];
                (*pupko_candidate_ancestor_node)[x] = VectorTools::whichMax(*leavesLikelihoods_node_i);
            }
        }

    } else {
        // Otherwise:
        // Set all likelihoods to 1 for a start:

        Vdouble *likelihoodArray_i = &likelihoodArray[nodeId];
        Vdouble *characterArray_i = &characterArray[nodeId];
        likelihoodArray_i->resize(nbStates_);
        characterArray_i->resize(nbStates_);
        for (size_t x = 0; x < nbStates_; x++) {
            (*likelihoodArray_i)[x] = 1.;
            (*characterArray_i)[x] = 1.;
        }

        size_t nbNodes = node->getNumberOfSons();

        if (sonNode) {
            throw ("PIPDRHomogeneousTreeLikelihood::computeLikelihoodAtSite_(...). 'sonNode' not found as a son of 'node'.");
        }
        vector<const Vdouble *> iLik;
        vector<const VVVdouble *> tProb;
        vector<const Vdouble *> cLik;
//        vector<vector<size_t>> *candidate_ancestor_node = &candidateAncestral[nodeId];
        for (size_t n = 0; n < nbNodes; n++) {
            const Node *son = node->getSon(n);
            cLik.push_back(&characterArray[son->getId()]);
            tProb.push_back(
                    &pxy_[son->getId()]);// all zero AbstractHomogeneousTreeLikelihood-> solved in PIPDRHomogeneousTreeLikelihood::init
            iLik.push_back(&likelihoodArray[son->getId()]);

        }

//        if (node->hasFather()) {
        if (node->getId() != newRootId) {
            const Node *father = node->getFather();
            computeLikelihoodFromArrays(iLik, tProb, cLik,
                                             &(*likelihoods_node)[father->getId()], &pxy_[nodeId],
                                             *likelihoodArray_i, *characterArray_i, siteNumber,
                                             nbNodes, nbDistinctSites_, nbClasses_, nbStates_);


        } else {
            const Vdouble *p = &model_->getFrequencies();
            computeLikelihoodFromArraysForRoot(iLik, tProb, cLik, p, //&(*likelihoods_node)[nodeId],
                                                    *likelihoodArray_i, *characterArray_i, siteNumber, nbNodes,
                                                    nbDistinctSites_, nbClasses_, nbStates_);

        }
    }
}

/*******************************************************************************/

void PIPDRHomogeneousTreeLikelihood::computeLikelihoodFromArraysForRoot(const vector<const Vdouble *> &iLik,
                                                                             const vector<const VVVdouble *> &tProb,
                                                                             const vector<const Vdouble *> &cLik,
                                                                             const Vdouble *p, Vdouble &oLik, Vdouble &oCLik,
                                                                             size_t siteNb, size_t nbNodes, size_t nbDistinctSites,
                                                                             size_t nbClasses,size_t nbStates) {

    Vdouble *oLik_i = &(oLik);
    Vdouble *oCLik_i = &(oCLik);
    // For each initial state
//    Vdouble likelihood {1};
    for (size_t n = 0; n < nbNodes; n++) {
        const Vdouble *iLik_n = iLik[n];
        const Vdouble *cLik_n = cLik[n];
        for (size_t x = 0; x < nbStates; x++) {
//            cout << n << x << (*iLik_n)[x] << endl;
            (*oLik_i)[x] *= (*iLik_n)[x];
//            cout << (*oLik_i)[x] << endl;
        }
    }
    for (size_t y = 0; y < nbStates; y++) {
//        cout << "2." << (*oLik_i)[y] << endl;
//        cout << "3." << (*p)[y] << endl;
        // We store this conditional likelihood into the corresponding array
        (*oLik_i)[y] *= (*p)[y];
//        cout << y << "\t" << (*oLik_i)[y] << endl;

//        (*oLik_i)[y] *= likelihood[y];
    }
    // We pick the argmax corresponding to the root
    (*oCLik_i)[0] = VectorTools::whichMax(*oLik_i);
}

/*******************************************************************************/
void PIPDRHomogeneousTreeLikelihood::computeLikelihoodFromArrays(
        const std::vector<const Vdouble*> &iLik,
        const std::vector<const VVVdouble*> &tProb,
        const std::vector<const Vdouble*> &cLik,
        const VVVdouble *iLikR,
        const VVVdouble *tProbR,
        Vdouble &oLik,
        Vdouble &oCLik,
        size_t siteNb,
        size_t nbNodes,
        size_t nbDistinctSites,
        size_t nbClasses,
        size_t nbStates) {

    Vdouble *oLik_i = &(oLik);
    Vdouble *oCLike_i = &(oCLik);
    const VVdouble *pxy_fa_c = &(*tProbR)[0];


    for (size_t x = 0; x < nbStates-1; x++) {
        Vdouble likelihood;
        likelihood.resize(nbStates);
        // Set all likelihoods to 1 for a start:
        for (size_t i = 0; i < nbStates; i++) {
            likelihood[i] = 1;
        }

        for (size_t n = 0; n < nbNodes; n++) {
            const Vdouble *iLik_n = iLik[n];
            const Vdouble *cLik_n = cLik[n];
//            const VVVdouble *pxy_n = tProb[n];

            for (size_t y = 0; y < nbStates; y++) {
//                cout << "1:" << (*iLik_n)[y] << endl;
//                cout << "2:" << (*cLik_n)[y] << endl;
                likelihood[y] *= (*iLik_n)[y];// left and right child multiplication
//                likelihood += (*pxy_n_c_x)[y] * (*iLik_n)[y];
//                cout << i << "\t" << c << "\t" << x << "\t" << y << "\t" <<  (* pxy__son_c_x)[y] << "\t" << (* likelihoods_root_son_i_c)[y] << endl;
            }
        }

        Vdouble lar;
        lar.resize(nbStates);
        for (size_t z = 0; z < nbStates; z++) {
            // We store this conditional likelihood into the corresponding array:
            lar[z] = (*pxy_fa_c)[x][z] * likelihood[z];
        }
        double winL = VectorTools::max(lar);
        (*oLik_i)[x] = winL;
        unsigned long winC = VectorTools::whichMax(lar);
        (*oCLike_i)[x] = winC;
    }
    //Gap state:
    (*oLik_i)[nbStates-1] = 0;
    (*oCLike_i)[nbStates-1] = nbStates-1;





}
/*******************************************************************************/
/**
void PIPDRHomogeneousTreeLikelihood::computeLikelihoodAtNode_(const Node *node, VVVdouble &likelihoodArray, const Node *sonNode) const {
    // const Node* node = tree_->getId();
    int nodeId = node->getId();
    likelihoodArray.resize(nbDistinctSites_);
    map<int, VVVdouble> *likelihoods_node = &likelihoodData_->getLikelihoodArrays(nodeId);

    // Initialize likelihood array:
    if (node->isLeaf()) {
        VVdouble *leavesLikelihoods_node = &likelihoodData_->getLeafLikelihoods(nodeId);
        for (size_t i = 0; i < nbDistinctSites_; i++) {
            VVdouble *likelihoodArray_i = &likelihoodArray[i];
            Vdouble *leavesLikelihoods_node_i = &(*leavesLikelihoods_node)[i];
            likelihoodArray_i->resize(nbClasses_);
            for (size_t c = 0; c < nbClasses_; c++) {
                Vdouble *likelihoodArray_i_c = &(*likelihoodArray_i)[c];
                likelihoodArray_i_c->resize(nbStates_);
                for (size_t x = 0; x < nbStates_; x++) {
                    (*likelihoodArray_i_c)[x] = (*leavesLikelihoods_node_i)[x];
                }
            }
        }

    } else {
        // Otherwise:
        // Set all likelihoods to 1 for a start:
        for (size_t i = 0; i < nbDistinctSites_; i++) {
            VVdouble *likelihoodArray_i = &likelihoodArray[i];
            likelihoodArray_i->resize(nbClasses_);
            for (size_t c = 0; c < nbClasses_; c++) {
                Vdouble *likelihoodArray_i_c = &(*likelihoodArray_i)[c];
                likelihoodArray_i_c->resize(nbStates_);
                for (size_t x = 0; x < nbStates_; x++) {
                    (*likelihoodArray_i_c)[x] = 1.;
                }
            }
        }
    }
    size_t nbNodes = node->getNumberOfSons();

    vector<const VVVdouble *> iLik;
    vector<const VVVdouble *> tProb;
    bool test = false;
    for (size_t n = 0; n < nbNodes; n++) {
        const Node *son = node->getSon(n);
        if (son != sonNode) {
            tProb.push_back(&pxy_[son->getId()]);// all zero AbstractHomogeneousTreeLikelihood-> solved in PIPDRHomogeneousTreeLikelihood::init
            iLik.push_back(&(*likelihoods_node)[son->getId()]);// all one
        } else {
            test = true;
        }
    }
    if (sonNode) {
        if(test)
            nbNodes--;
        else
            throw ("PIPDRHomogeneousTreeLikelihood::computeLikelihoodAtNode_(...). 'sonNode' not found as a son of 'node'.");
    }

    if (node->hasFather()) {
        const Node *father = node->getFather();
        computeLikelihoodFromArrays(iLik, tProb, &(*likelihoods_node)[father->getId()], &pxy_[nodeId], likelihoodArray,
                                    nbNodes, nbDistinctSites_, nbClasses_, nbStates_,
                                    false);

    } else {
        computeLikelihoodFromArrays(iLik, tProb, likelihoodArray, nbNodes, nbDistinctSites_, nbClasses_, nbStates_,
                                    false);
        // We have to account for the equilibrium frequencies:
        for (size_t i = 0; i < nbDistinctSites_; i++) {
            VVdouble *likelihoodArray_i = &likelihoodArray[i];
            for (size_t c = 0; c < nbClasses_; c++) {
                Vdouble *likelihoodArray_i_c = &(*likelihoodArray_i)[c];
                for (size_t x = 0; x < nbStates_; x++) {
                    (*likelihoodArray_i_c)[x] *= rootFreqs_[x];
                }
            }
        }
    }
}
**/
/******************************************************************************/

