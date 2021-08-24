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
 * @file PIPDRTreeLikelihoodData.h
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


#ifndef ARPIP_PIPDRTREELIKELIHOODDATA_H
#define ARPIP_PIPDRTREELIKELIHOODDATA_H

// From arpip:
#include "../SubstitutionModel/PIP13.h"

// From PhylLib:
#include "Bpp/Phyl/Likelihood/AbstractTreeLikelihoodData.h"
#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include "Bpp/Phyl/PatternTools.h"
#include "Bpp/Phyl/SitePatterns.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>

// From the STL:
#include <map>

namespace bpp {

    //------------------------------------------------------------------------------------------------------------------
    /**
    * @brief Likelihood data structure for a leaf.
    *
    * This class is for use with the PIPDRTreeLikelihoodData class.
    *
    * Store the likelihoods arrays associated to a leaf.
    *
    * @see PIPDRTreeLikelihoodLeafData
    */
    class PIPDRTreeLikelihoodLeafData:
            public virtual TreeLikelihoodNodeData {
    private:
        mutable VVdouble leafLikelihood_;
        const Node *leaf_;

    public:
        /**
        * @brief Constructor.
        */
        PIPDRTreeLikelihoodLeafData();

        /**
        * @brief Copy constructor.
        */
        PIPDRTreeLikelihoodLeafData(const PIPDRTreeLikelihoodLeafData &data);

        /**
        * @brief Operator overloading.
        */
        PIPDRTreeLikelihoodLeafData &operator=(const PIPDRTreeLikelihoodLeafData &data);

        /**
        * @brief Destructor.
        */
//        virtual ~PIPDRTreeLikelihoodLeafData();

        PIPDRTreeLikelihoodLeafData *clone() const { return new PIPDRTreeLikelihoodLeafData(*this); }

    public:
        const Node *getNode() const { return leaf_; }
        void setNode(const Node *node){ leaf_ = node; }

        VVdouble &getLikelihoodArray(){ return leafLikelihood_; }
    };

    //------------------------------------------------------------------------------------------------------------------
    /**
    * @brief Likelihood data structure for a node.
    *
    * This class is for use with the PIPDRTreeLikelihoodData class.
    *
    * Store for each neighbor node an array with conditional likelihoods.
    *
    * @see PIPDRTreeLikelihoodNodeData
    */
    class PIPDRTreeLikelihoodNodeData:
            public virtual TreeLikelihoodNodeData {
    private:
        /**
        * @brief This contains all likelihood values used for computation.
        *
        * <pre>
        * x[b][i][c][s]
        *   |------------> Neighbor node of n (id)
        *     |---------> Site i
        *         |------> Rate class c
        *            |---> state s
        * </pre>
        * We call this the <i>likelihood array</i> for each node.
        */
        mutable std::map<int, VVVdouble> nodeLikelihoods_;

        // All  first and second order likelihood derivatives values used for computation.
        mutable Vdouble nodeDLikelihoods_;
        mutable Vdouble nodeD2Likelihoods_;

        const Node *node_;

    public:
        /**
        * @brief Constructor.
        */
         PIPDRTreeLikelihoodNodeData() ;

        /**
        * @brief Copy constructor.
        */
        PIPDRTreeLikelihoodNodeData(const PIPDRTreeLikelihoodNodeData &data);

        /**
        * @brief Operator overloading.
        */
        PIPDRTreeLikelihoodNodeData &operator=(const PIPDRTreeLikelihoodNodeData &data);

        /**
        * @brief Destructor.
        */
        virtual ~PIPDRTreeLikelihoodNodeData(){}

        PIPDRTreeLikelihoodNodeData *clone() const{ return new PIPDRTreeLikelihoodNodeData(*this); }

    public:
        const Node *getNode() const { return node_; }

        void setNode(const Node *node) { node_ = node; }

        const std::map<int, VVVdouble> &getLikelihoodArrays() const { return nodeLikelihoods_; }

        std::map<int, VVVdouble> &getLikelihoodArrays() { return nodeLikelihoods_; }

        VVVdouble &getLikelihoodArrayForNeighbor(int neighborId)
        {
            return nodeLikelihoods_[neighborId];
        }

        const VVVdouble &getLikelihoodArrayForNeighbor(int neighborId) const
        {
            return nodeLikelihoods_[neighborId];
        }

        bool isNeighbor(int neighborId) const
        {
            return nodeLikelihoods_.find(neighborId) != nodeLikelihoods_.end();
        }

        void eraseNeighborArrays()
        {
            nodeLikelihoods_.erase(nodeLikelihoods_.begin(), nodeLikelihoods_.end());
        }

        // Derivatives

        Vdouble &getDLikelihoodArray() { return nodeDLikelihoods_;  }

        const Vdouble &getDLikelihoodArray() const  {  return nodeDLikelihoods_;  }

        Vdouble &getD2LikelihoodArray()  {  return nodeD2Likelihoods_; }

        const Vdouble &getD2LikelihoodArrayForNeighbor() const  { return nodeD2Likelihoods_; }

    };

    //------------------------------------------------------------------------------------------------------------------
    /**
    * @brief PIP data structure for likelihood computation.
    */
    class PIPDRTreeLikelihoodParameters{
    private:
        const TreeTemplate<Node>* tree_;

        /**
         *  @brief This contains total tree length \tau used for computation.
         *
         * We compute this <i>\tua</i> once for the tree.
         */
        mutable double tau_;

        /**
         *  @brief This contains insertion rate \lambda used for computation.
         *
         *
         * We compute this the <i>\lambda</i> once for the PIP.
         */
        mutable float lambda_;

        /**
        *  @brief This contains deletion rate \mu used for computation.
        *
        *
        * We compute this the <i>\mu</i> once for the PIP.
        */
        mutable float mu_;

        /**
        * @brief This contains all \nu values used for computation.
        *
        *
        * We call this the <i>\nu</i> for the PIP.
        *
        \f[
         $\nu = \lambda \left(||\tau|| + \nicefrac{1}{\mu}\right)$.
        \f]
        */
        mutable double nu_;

        /**
         * @brief This contains all Iota values used for computation.
         *
         * <pre>
         * x[n]
         *   |---------> node n
         * </pre>
         *
         * We call this the <i>iota array</i> for each node.
         * * For an internal node $v$, the probability $\iota(v)$ of inserting a single character on branch
         * $pa(v) \rightarrow v$, is proportional to branch length $b(v)$. For $v\neq \Omega$ it is given by
         * $\iota(v)=b(v)/(|\tau|+ \mu^{-1})$; at the root atomic mass point probability
         * $\iota(\Omega) = \mu^{-1}/(|\tau|+\mu^{-1})$ so that $\sum_{v\subset V} \iota (v) = 1$.
         */
        mutable Vdouble nodeIotaData_;

        /**
         *  @brief This contains all Beta values used for computation.
         *
         * <pre>
         * x[n]
         *   |---------> node n
         * </pre>
         *
         * We call this the <i>Beta array</i> for each node.
         \f[
             \beta(v) =
            \nicefrac{(1 - e^{{-\mu}b(v)})}{(\mu b(v))} & \text{if $v \neq \Omega$}\\
            1 & \text{if $v = \Omega$},
          \f]
         */
        mutable Vdouble nodeBetaData_;


    public:
        /**
        * @brief Constructor without parameter for scenario without PIP.
        */
        PIPDRTreeLikelihoodParameters();

        /**
        * @brief Constructor with PIP parameter for scenario with PIP.
        */
        PIPDRTreeLikelihoodParameters(const TreeTemplate<Node> *tree, const double lambda, const double mu);

        /**
        * @brief Copy constructor.
        */
        PIPDRTreeLikelihoodParameters(const PIPDRTreeLikelihoodParameters &data);

        /**
        * @brief Operator overloading.
        */
        PIPDRTreeLikelihoodParameters &operator=(const PIPDRTreeLikelihoodParameters &data);

        /**
        * @brief Destructor.
        */
        virtual ~PIPDRTreeLikelihoodParameters() { }

        PIPDRTreeLikelihoodParameters *clone() const { return new PIPDRTreeLikelihoodParameters(*this); }

    public:
        /**
        * @brief Method to initialize PIP model Parameter
        * called by constructor
        */
        void initPIPParameter_();

        /**
         * @brief Recursive formula associated to the node v that computes the character deletion probability on the
         * sub-tree rooted at node v.
         *
         */
         void computeNu_();

        /**
       * @brief PIP likelihood: Compute array of intensity of insertion at nodes.
       * @author Gholamhossein Jowkar
       *
       * This method is to be called when the PIP begins.
       *
       * For an internal node $v$, the probability $\iota(v)$ of inserting a single character on branch
       * $pa(v) \rightarrow v$, is proportional to branch length $b(v)$. For $v\neq \Omega$ it is given by
       * $\iota(v)=b(v)/(|\tau|+ \mu^{-1})$; at the root atomic mass point probability
       * $\iota(\Omega) = \mu^{-1}/(|\tau|+\mu^{-1})$ so that $\sum_{v\subset V} \iota (v) = 1$.
       *
       * NB: This method is recursive.
       *
       * @param node The node defining the subtree to analyse.
       * @param tau Total tree length.
       * @param mu The deletion rate.
       * @param isRoot The flag for recursive iteration.
       *
       * @see [Bouchard-Côté, Jordan, (2013), PNAS]
       */
        void computePIPIota(const Node *node, const double tau, const double mu, bool isRoot);

        /**
        * @brief PIP likelihood: Compute array of survival probability of inserted character at inner nodes.
        * @author Gholamhossein Jowkar
        *
        * This method is to be called when the PIP begins.
        *
        * For an internal node $v$,  the survival probability $\beta =(v)$ associated with an inserted character on
        * branch $pa(v) \rightarrow v$ is given by $ \beta(\Omega)= 1$ and $\beta(v) = (1-exp({-\mu}b(v)))/(\mu b(v))$.
        *
        * NB: This method is recursive.
        * @param node The node defining the subtree to analyse.
        * @param mu The deletion rate.
        * @param isRoot The flag for recursive iteration.
        *
        * @see [Bouchard-Côté, Jordan, (2013), PNAS]
        */
        void computePIPBeta(const Node *node, const double mu, bool isRoot);

        /********************************** Setter and Getter ********************************************/
        const TreeTemplate<Node> *getTree() const { return tree_;}

        void setTree(const TreeTemplate<Node> *tree) { tree_ = tree;}

        double getTau() const { return tau_; }

        void setTau(double tau) { tau_ = tau; }

        float getMu() const { return mu_; }

        void setMu(float mu) { mu_ = mu; }

        float getLambda() const { return lambda_; }

        void setLambda(float lambda) { lambda_ = lambda; }

        double getNu() const { return nu_; }

        void setNu(double nu) { nu_ = nu; }

        const Vdouble &getNodeIotaData() const { return nodeIotaData_; }
        const double &getNodeIotaData(int nodeId) const { return nodeIotaData_[nodeId]; }

        void setNodeIotaData(const Vdouble &nodeIotaData) { nodeIotaData_ = nodeIotaData; }

        const Vdouble &getNodeBetaData() const {return nodeBetaData_; }
        const double &getNodeBetaData(int nodeId) const {return nodeBetaData_[nodeId]; }

        void setNodeBetaData(const Vdouble &nodeBetaData) { nodeBetaData_ = nodeBetaData; }

    };


    //------------------------------------------------------------------------------------------------------------------
    /**
    * @brief Likelihood data structure for uniform rate sites models, using a double-recursive algorithm.
    */
    class PIPDRTreeLikelihoodData:
            public virtual AbstractTreeLikelihoodData {
    private:
        mutable std::map<int, PIPDRTreeLikelihoodNodeData> nodeData_;
        mutable std::map<int, PIPDRTreeLikelihoodLeafData> leafData_;
        mutable VVVdouble rootLikelihoods_;
        mutable VVdouble  rootLikelihoodsS_;
        mutable Vdouble   rootLikelihoodsSR_;

        mutable PIPDRTreeLikelihoodParameters pipParam_;

        /**
       * @brief This contains all gap flags used for computation.
       *
       * <pre>
       * x[n][i]
       *   |---------> node n
       *      |------> site i
       * </pre>
       * We call this the <i>Gaps flag array</i> for each node.
       */
        mutable std::map<size_t, std::vector<bool>> nodeSetG_;

        /**
         * @brief This contains all setA flags used for computation.
         *
         * <pre>
         * x[n][i]
         *   |---------> node n
         *      |------> site i
         * </pre>
         * We call this the <i> setA flag array</i> for each node.
         */
        mutable std::map<size_t, std::vector<bool>> nodeSetA_;

        SiteContainer* shrunkData_;
        size_t nbSites_;
        size_t nbStates_;
        size_t nbClasses_;
        size_t nbDistinctSites_;

    public:
        /**
        * @brief Constructor with PIP.
        */
        PIPDRTreeLikelihoodData(const TreeTemplate<Node> *tree, const double mu, const double lambda);

        /**
        * @brief Constructor without PIP.
        */
        PIPDRTreeLikelihoodData(const TreeTemplate<Node> *tree);

        /**
        * @brief Copy constructor.
        */
        PIPDRTreeLikelihoodData(const PIPDRTreeLikelihoodData &data);

        /**
        * @brief Operator overloading.
        */
        PIPDRTreeLikelihoodData &operator=(const PIPDRTreeLikelihoodData& data);

        /**
        * @brief Destructor.
        */
        virtual ~PIPDRTreeLikelihoodData() { delete shrunkData_; }

        PIPDRTreeLikelihoodData *clone() const { return new PIPDRTreeLikelihoodData(*this); }

    public:
        /******************************* Normal Likelihood  functions *************************************************/

        /**
        * @brief Set the tree associated to the data.
        *
        * All node data will be actualized accordingly by calling the setNode() method on the corresponding nodes.
        * @warning: the old tree and the new tree must be two clones! And particularly, they have to share the
        * same topology and nodes id.
        *
        * @param tree The tree to be associated to this data.
        */
        void setTree(const TreeTemplate<Node> *tree);

        /**
        * @brief Resize and initialize all likelihood arrays according to the given data set and substitution model.
        *
        * @param sites The sequences to use as data.
        * @param model The substitution model to use.
        * @throw Exception if an error occurs.
        */
        void initLikelihoods(const SiteContainer &sites, const TransitionModel &model);


        /**
        * @brief Rebuild likelihood arrays at inner nodes.
        *
        * This method is to be called when the topology of the tree has changed.
        * Node arrays relationship are rebuilt according to the new topology of the tree.
        * The leaves likelihood remain unchanged, so as for the first and second order derivatives.
        */
        void reInit();

        void reInit(const Node *node);

        /**
        * @brief PIP likelihood: Number of non-gap char states for each site
        * @author Gholamhossein Jowkar
        *
        * This method is used by PIPDRHomogeneousTreeLikelihood to initialize PIP topological flags/sets.
        *
        */
        void setPIPTopologicalFlags();


        /**
        * @brief PIP likelihood: Number of non-gap char states for each site
        * @author Gholamhossein Jowkar
        *
        * This method is used by AssignPIPSetA.
        *
        * @param node The node defining the subtree to analyse.
        * @param sites The sequence container to use.
        */
        std::vector<size_t> initNumberOfNonGapChar();

        /**
        * @brief PIP likelihood: This method called to put a flag on ancestral states candidate based on observed char
         * states in order to enhance computation at inner nodes. (The places that insertion might be happened.)
         *
        * @author Gholamhossein Jowkar
        *
        *
        * NB: This method is recursive.
        *
        * @param node The node defining the subtree to analyse.
        * @param sites The sequence container to use.
        * @param numberOfNonGapChar The number of character in sites that are not Gap.
        * @param nodeSumCharacter  The number of non-Gap characters in children of corresponding node.
        */
        void computePIPSetA(const Node *node, const std::vector<size_t> numberOfNonGapChar,
                            std::map<size_t, std::vector<int>> &nodeSumCharacter){
            for (size_t siteNumb = 0; siteNumb < nbDistinctSites_; siteNumb++) {
                Site columnSite = shrunkData_->getSite(siteNumb);
//        std::cout << std::endl <<"COL" << siteNumb << " is: " << columnSite.toString();
                computePIPSetAForSite(node, columnSite, numberOfNonGapChar, nodeSumCharacter, siteNumb);
            }
        }

        /**
        * @brief PIP likelihood: This method called to put a flag on ancestral states candidate based on observed char
         * states in order to enhance computation at inner nodes. (The places that insertion might be happened.)
         *
        * @author Gholamhossein Jowkar
        *
        *
        * NB: This method is recursive.
        *
        * @param node The node defining the subtree to analyse.
        * @param sites The sequence container to use.
        * @param numberOfNonGapChar The number of character in sites that are not Gap.
        * @param nodeSumCharacter  The number of non-Gap characters in children of corresponding node.
        * @param nbSite The site number.
        */
        void computePIPSetAForSite(const Node *node, const Site &colSite,
                                   const std::vector<size_t> numberOfNonGapChar,
                                   std::map<size_t, std::vector<int>> &nodeSumCharacter, size_t nbSite);

        /**
        * @brief PIP likelihood: This method called to put a flag on common ancestral states of gaps to enhance
        * computation at inner nodes.  (The places that deletion might be happened.)
        * @author Gholamhossein Jowkar
        *
        * This method is to be called when the PIP begins to put a flag on ancestral states of gap char for further
        * computation enhancement.
        *
        * NB: This method is recursive.
        *
        * @param node The node defining the subtree to analyse.
        * @param sites The sequence container to use.
        */
        void computePIPSetG(const Node *node){
            for (size_t siteNumb = 0; siteNumb < nbDistinctSites_; ++siteNumb) {
                Site columnSite = shrunkData_->getSite(siteNumb);
//                std::cout << std::endl << "COL" << siteNumb << " is: " << columnSite.toString();
                computePIPSetGForSite(node, columnSite, siteNumb);
            }
        }

        /**
        * @brief PIP likelihood: This method called to put a flag on common ancestral states of gaps to enhance
        * computation at inner nodes.  (The places that deletion might be happened.)
         *
        * @author Gholamhossein Jowkar
        *
        *
        * NB: This method is recursive.
        *
        * @param node The node defining the subtree to analyse.
        * @param sites The sequence container to use.
        * @param nbSite The site number.
        */
        void computePIPSetGForSite(const Node *node, const Site &colSite, size_t siteNb );

        /**
        * @brief Rebuild likelihood arrays by PIP.
        *
        * called by PIPDRHomogeneousTreeLikelihood::firePIPParameterChanged_
        */
        void firParameterChanged(const double lambda, const double mu);


        /****************************************** Normal getter  ****************************************************/

        const PIPDRTreeLikelihoodParameters &getPipParam() const { return pipParam_; };

        const  std::vector<bool> &getNodeSetG(int nodeId) const { return nodeSetG_[nodeId]; };

        const std::vector<bool> &getNodeSetA( int nodeId) const { return nodeSetA_[nodeId]; };

        PIPDRTreeLikelihoodNodeData &getNodeData(int nodeId) { return nodeData_[nodeId]; }

        const PIPDRTreeLikelihoodNodeData &getNodeData(int nodeId) const { return nodeData_[nodeId]; }

        PIPDRTreeLikelihoodLeafData &getLeafData(int nodeId)  {  return leafData_[nodeId];  }

        const PIPDRTreeLikelihoodLeafData &getLeafData(int nodeId) const    { return leafData_[nodeId]; }

        size_t getArrayPosition(int parentId, int sonId, size_t currentPosition) const { return currentPosition; }

        const std::map<int, VVVdouble> &getLikelihoodArrays(int nodeId) const
        {
            return nodeData_[nodeId].getLikelihoodArrays();
        }

        std::map<int, VVVdouble>& getLikelihoodArrays(int nodeId)
        {
            return nodeData_[nodeId].getLikelihoodArrays();
        }

        VVVdouble &getLikelihoodArray(int parentId, int neighborId)
        {
            return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
        }

        const VVVdouble &getLikelihoodArray(int parentId, int neighborId) const
        {
            return nodeData_[parentId].getLikelihoodArrayForNeighbor(neighborId);
        }

        // Derivatives
        Vdouble &getDLikelihoodArray(int nodeId) { return nodeData_[nodeId].getDLikelihoodArray(); }

        const Vdouble &getDLikelihoodArray(int nodeId) const { return nodeData_[nodeId].getDLikelihoodArray(); }

        Vdouble &getD2LikelihoodArray(int nodeId) { return nodeData_[nodeId].getD2LikelihoodArray(); }

        const Vdouble &getD2LikelihoodArray(int nodeId) const { return nodeData_[nodeId].getD2LikelihoodArray(); }

        VVdouble &getLeafLikelihoods(int nodeId) { return leafData_[nodeId].getLikelihoodArray(); }

        const VVdouble &getLeafLikelihoods(int nodeId) const { return leafData_[nodeId].getLikelihoodArray(); }

        VVVdouble &getRootLikelihoodArray() { return rootLikelihoods_; }
        const VVVdouble &getRootLikelihoodArray() const { return rootLikelihoods_; }

        VVdouble &getRootSiteLikelihoodArray() { return rootLikelihoodsS_; }
        const VVdouble &getRootSiteLikelihoodArray() const { return rootLikelihoodsS_; }

        Vdouble &getRootRateSiteLikelihoodArray() { return rootLikelihoodsSR_; }
        const Vdouble &getRootRateSiteLikelihoodArray() const { return rootLikelihoodsSR_; }

        size_t getNumberOfDistinctSites() const { return nbDistinctSites_; }

        size_t getNumberOfSites() const { return nbSites_; }

        size_t getNumberOfStates() const { return nbStates_; }

        size_t getNumberOfClasses() const { return nbClasses_; }

        const SiteContainer* getShrunkData() const { return shrunkData_; }

    protected:
        /**
         * @brief This method initializes the leaves according to a sequence container.
         *
         * Here the container shrunkData_ is used.
         * Likelihood is set to 1 for the state corresponding to the sequence site,
         * otherwise it is set to 0.
         *
         * All likelihood arrays at each nodes are initialized according to alphabet
         * size and sequences length, and filled with 1.
         *
         * NB: This method is recursive.
         *
         * @param node  The node defining the subtree to analyse.
         * @param sites The sequence container to use.
         * @param model The model, used for initializing leaves' likelihoods.
         */
        void initLikelihoods(const Node *node, const SiteContainer &sites, const TransitionModel &model);


    };



}// end of namespace
#endif //ARPIP_PIPDRTREELIKELIHOODDATA_H
