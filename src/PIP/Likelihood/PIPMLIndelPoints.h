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
 * @file PIPMLIndelPoints.h
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

#ifndef ARPIP_PIPMLINDELPOINTS_H
#define ARPIP_PIPMLINDELPOINTS_H

// From bpp-phy:
#include <Bpp/Phyl/AncestralStateReconstruction.h>

// From bpp-arpip
#include "PIPDRTreeLikelihood.h"
#include "PIPDRHomogeneousTreeLikelihood.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From STL:
#include<vector>

namespace bpp {
    /**
     * @brief MLIndelPoints data structure for a node.
     *
     * This class is for use with the PIPDRTreeLikelihoodData class.
     *
     * Performs the indel reconstruction on the tree given msa.
     *
     * @see PIPMLIndelPoints
     */
    class PIPMLIndelPoints {

    private:
        /**
         * @brief This contains the PIP marginal likelihood data structure and all the relative variables.
         */
        const PIPDRTreeLikelihood *likelihood_;

        const TreeTemplate<Node> *tree_;
//        const Tree *tree_;
        const SiteContainer *shrunkData_;
        const TransitionModel *model_;
        size_t nbSites_;
        size_t nbStates_;
        size_t nbClasses_;
        size_t nbDistinctSites_;

        /**
         * @brief This contains all Iota values used for computation.
         *
         * <pre>
         * x[n]
         *   |---------> node n
         * </pre>
         * We call this the <i>Iota array</i> for each node.
         */
        const Vdouble *nodeIota_;
        /**
         * @brief This contains all Beta values used for computation.
         *
         * <pre>
         * x[n]
         *   |---------> node n
         * </pre>
         * We call this the <i>Beta array</i> for each node.
         */
        const Vdouble *nodeBeta_;

        /**
         * @brief This contains all gap flags used for computation defined as the set of all internal nodes where all
         * descendant leaves have a gap in the respective MSA column.
         *
         * <pre>
         * x[n][i]
         *   |---------> node n
         *      |------> site i
         * </pre>
         *
         * We call this the <i>set of potential deletion points</i> for each node.
         */
        mutable std::map<size_t, std::vector<bool>> nodeSetG_;

        /**
         * @brief This contains all setA flags used for computation defined as follow:
         * We consider subset S of leaves that have a non-gap character,S={v∈L:mv6=ε}. Given set S, we define setA of
         * potential insertion points to include all nodes that are ancestral to all the leaves in S.
         *
         * <pre>
         * x[n][i]
         *   |---------> node n
         *      |------> site i
         * </pre>
         * We call this the <i> potential insertion points </i> for each node.
         */
        mutable std::map<size_t, std::vector<bool>> nodeSetA_;

        /**
        * @brief This contains all indel scenarios used for computation.
        *
        * <pre>
        * x[n][i][s]
        *   |------------> id of node (n)
        *      |---------> Site i
         *        |-----------> scenario s
        * </pre>
        * We call this the <i> indel text array </i> for each node.
        */
        mutable std::map<int, std::vector<std::string>> nodeText_;
        mutable std::map<int, std::vector<std::string>> nodeFelsensteinValueText_;

        mutable std::map<int, Vdouble> nodeProb_;
        mutable std::map<int, Vdouble> nodeFelsensteinValue_;

        mutable std::vector<std::string> homoPath_;
        mutable std::map<size_t, std::vector<std::string>> insertionPoints_;
        mutable std::map<size_t, std::vector<std::string>> deletionPoints_;

    public:
        /**
         * @brief Constructor
         */
        PIPMLIndelPoints(const PIPDRHomogeneousTreeLikelihood *lik);

        /**
         * @brief Copy Constructor
         */
        PIPMLIndelPoints(const PIPMLIndelPoints &mlindel);

        /**
         * @brief Operator overloading
         */
        PIPMLIndelPoints &operator=(const PIPMLIndelPoints &mlindel);

        /**
         * Destructor
         */
        virtual ~PIPMLIndelPoints() {delete shrunkData_;}

        PIPMLIndelPoints *clone() const { return new PIPMLIndelPoints(*this); }

        /******************************************* PIP specific functions ***************************************************/

        /**
         * @author Gholamhossein Jowkar
         *
         * @param node The root but generally the node defining the subtree to analyse.
         * @param sites The sequence container to use.
         * @param model The model, used for initializing leaves' likelihoods.
         */
        void initMaximumLikelihoodIndelPoints_(const Node *node);


        /**
        * @author Gholamhossein Jowkar/ Massimo Maiolo
        * @brief This method call the functions that extract Indel Points.
        *
        * NB: This method is recursive.
        *
        * @param node The root but generally the node defining the subtree to analyse.
        * @param mu deletion rate
        */
        void extractIndelPoint(const Node *node, const float mu);

        void extractIndelPointForSite_(const Node *node, const float mu, size_t siteNumber);

        /**
         * @brief This method, return insertion points and deletion points as text file.
         *
         * @param text
         */
        void MLHomoPath(std::vector<std::string> &text);


    protected:

        /**
         * @author Gholamhossein Jowkar/ Massimo Maiolo
         *
         * @brief PIP likelihood: This method compute likelihood function of a leaf node to inference the homology path
         * under PIP which progressively build the indel history.
         *
         * @param node The leaf node but generally the node defining the subtree to analyse.
         */
        void leafHomologyPath_(const Node *node, size_t siteNumber);

        /**
        * @brief PIP likelihood: This method compute likelihood function for subtree containing gaps only.
        * @author Gholamhossein Jowkar
        *
         * @note we have to infer whether the deletion happened at the node $v$ or its children. For this purpose, the
         * probability of homology path for each individual scenario would be computed by picking the maximum value of
         * pure survival/non-survival probability.
        *
        * @param node The internal node but generally the node defining the subtree to analyse.
        * @param childNode1 The first son of the parent node.
        * @param childNode2 The second son of the parent node.
        * @param mu The deletion rate.
         */
        void GapOnlySubtreeHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2,
                const float mu, const size_t siteNumber);

        /**
         * @brief In this case the sub-tree rooted at node $v$ contains no immediate child that is a member of the
         * potential deletion points set $\mathcal{G}$, so the character inserted at node $v$ must survive at least to
         * the grandchildren of node $v$.
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void NeitherChildIsGapOnlySubtreeHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2, const float mu,
                const size_t siteNumber);

        /**
         * @brief In this case, we are considering the sub-tree containing two leaves with one or more (called by
         * OneChildLeafSubtreeHomologyPath_) gap characters.
         *
         * @note function A
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void BothChildrenLeavesOrOneLeafSubtreeHomologyPath_(const Node *node, const Node *childNode1,
                const Node *childNode2, const float mu, const size_t siteNumber);


        /**
         * @brief In this case, the root of sub-tree contains one child leaf node. If the two children of root were also
         * leaves, we would apply the BothChildrenLeavesOrOneLeafSubtreeHomologyPath_ , otherwise, the sub-tree contains
         * no leaf children at all. Hence, the NoChildLeafHomologyPath_ should be called.
         *
         * @note function B
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void OneChildLeafSubtreeHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2,
                const float mu, const size_t siteNumber);


        /**
         * @brief In this subsection, it is assumed that the sub-tree contains no leaves so we have to check which one
         * of the sub-trees is a gap only sub-tree (full gap tree). And the probability of the homology path is computed.
         *
         * @note function C
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void NoChildLeafHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2,
                const float mu, const size_t siteNumber);

        /**
         * @brief function D
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void OneChildGapHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2, const float mu,
                const size_t siteNumber);


        /**
         * @brief function E
         *
         * @param node The node defining the subtree to analyse.
         * @param childNode1 The first son of the parent node.
         * @param childNode2 The second son of the parent node.
         * @param mu The deletion rate.
         * @param siteNumber The site number used for computation.
         */
        void AllChildGapHomologyPath_(const Node *node, const Node *childNode1, const Node *childNode2, const float mu,
                const size_t siteNumber);



        /************************************* PIP getters and setters *************************************************/
    public:

        const PIPDRTreeLikelihood *getLikelihood() const { return likelihood_; }

        void setLikelihood(const PIPDRTreeLikelihood *likelihood) { likelihood_ = likelihood; }

        const TreeTemplate<Node> *getTree() const { return tree_; }

        void setTree(const TreeTemplate<Node> *tree) { tree_ = tree; }

        const SiteContainer *getShrunkData() const { return shrunkData_; }

        void setShrunkData(const SiteContainer *shrunkData) { shrunkData_ = shrunkData; }

        const TransitionModel *getModel() const { return model_; }

        void setModel(const TransitionModel *model) { model_ = model; }

        size_t getNbSites() const { return nbSites_; }

        void setNbSites(size_t nbSites) { nbSites_ = nbSites; }

        const std::map<size_t, std::vector<bool>> &getNodeSetG() const;

        void setNodeSetG(const std::map<size_t, std::vector<bool>> &nodeSetG);

        const std::map<size_t, std::vector<bool>> &getNodeSetA() const;

        void setNodeSetA(const std::map<size_t, std::vector<bool>> &nodeSetA);

        const std::map<int, std::vector<std::string>> &getNodeText() const;

        void setNodeText(const std::map<int, std::vector<std::string>> &nodeText);

        const std::map<int, std::vector<std::string>> &getNodeFelsensteinValueText() const;

        void setNodeFelsensteinValueText(const std::map<int, std::vector<std::string>> &nodeFelsensteinValueText);

        const std::map<int, Vdouble> &getNodeProb() const;

        void setNodeProb(const std::map<int, Vdouble> &nodeProb);

        const std::map<int, Vdouble> &getNodeFelsensteinValue() const;

        void setNodeFelsensteinValue(const std::map<int, Vdouble> &nodeFelsensteinValue);

        const std::vector<std::string> &getHomoPath() const;

        void setHomoPath(const std::vector<std::string> &homoPath);

        const std::map<size_t, std::vector<std::string>> &getInsertionPoints() const;

        void setInsertionPoints(const std::map<size_t, std::vector<std::string>> &insertionPoints);

        const std::map<size_t, std::vector<std::string>> &getDeletionPoints() const;

        void setDeletionPoints(const std::map<size_t, std::vector<std::string>> &deletionPoints);




    };

}// end namespace
#endif //ARPIP_PIPMLINDELPOINTS_H
