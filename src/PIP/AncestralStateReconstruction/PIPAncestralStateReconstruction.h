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
 * @file PIPAncestralStateReconstruction.h
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
#ifndef ARPIP_PIPANCESTRALSTATERECONSTRUCTION_H
#define ARPIP_PIPANCESTRALSTATERECONSTRUCTION_H

// From Phyl:
#include <Bpp/Phyl/AncestralStateReconstruction.h>

// From ARPIP:
#include "../Likelihood/PIPDRTreeLikelihood.h"
#include "../Likelihood/PIPMLIndelPoints.h"

// From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From STL:
#include<vector>

namespace bpp {
    class PIPAncestralStateReconstruction :
            public virtual AncestralStateReconstruction {
    private:
        const PIPDRTreeLikelihood *likelihood_;
        const PIPMLIndelPoints *mlIndelPoints_;
        TreeTemplate<Node> tree_;
        const Alphabet *alphabet_;
        size_t nbInnerNodes_;
        size_t nbSites_;
        size_t nbDistinctSites_;
        size_t nbStates_;
        std::vector<size_t> rootPatternLinks_;
        VVVdouble probabilityProfile_;
        std::string probProfileType_;
//        std::vector<double> r_;
//        std::vector<double> l_;

    public:
        /**
         * @brief Constructor
         */
        PIPAncestralStateReconstruction(const PIPDRTreeLikelihood *lik, const PIPMLIndelPoints *mlIndelPoints,
                                        const std::string &probProfileType = "none");

        /**
         * @brief Copy Constructor
         */
        PIPAncestralStateReconstruction(const PIPAncestralStateReconstruction &arpip);

        /**
         * @brief Operator Overloading
         */
        PIPAncestralStateReconstruction &operator=(const PIPAncestralStateReconstruction &arpip);

        /**
         * @brief Destructor
         */
         virtual ~PIPAncestralStateReconstruction()= default;

         PIPAncestralStateReconstruction *clone()const { return new PIPAncestralStateReconstruction(*this); }

    public:

        /****************************************** PIP base functions **********************************************/


        /**
         * @brief Get ancestral states for a given node at specific site using Pupko's method.
         *
         * The size of the vector is the number of distinct sites in the container
         * associated to the likelihood object.
         * This method is mainly for efficient internal use in other classes.
         * Consider using the recursiveGetAncestralStatesForAllNodes() method for a more general output.
         *
         * @param node The node at which the states must be reconstructed.
         * @param ancestors The array containing the ancestor matrix.
         * @param lArray The likelihood array containing stored likelihood values.
         * @param cArray The character state corresponding to likelihood values, read FastML paper for more info.
         * @param siteNumber  Site number that all the computation is done on.
         * @param newRootId The id of new root for the subtree we work on.
         * @see [Pupko, et al., 2000, Molecular biology and evolution]
         */
        void getAncestralStatesForSite(const Node *node, std::map<int, std::vector<int>> &ancestors, VVdouble &lArray,
                                       VVdouble &cArray, const size_t siteNumber, const int newRootId) const;


        /**
        * @brief Get all ancestral states for all nodes including gap (aka joint ancestral sequence reconstruction using DP approach with indel).
        *
         * @note This method considers gap as an additional character.
        *
        * Call the recursiveJointAncestralStates() method on each subtree per site.
        *
        * @return A map with nodes id as key, and a vector of states indices as value.
        * @see JointAncestralSequenceReconstruction and [Pupko, et al., 2000, Molecular biology and evolution]
        */
        std::tuple<std::map<int, std::vector<int>>, VVVdouble>  getAllAncestralStatesWGap() const;


        /**
        * @brief Get Joint ancestral sequences for all nodes including gap (aka joint ancestral sequence reconstruction using DP approach with indel).
        *
        * @note This method considers gap as an additional character. Moreover, you  can choose to have it with/without probability profile.
        *
        * Call the getAncestralSequences() method.
        * @return
        */
        AlignedSequenceContainer *JointAncestralSequencesReconstruction(){
            AlignedSequenceContainer *asc = new AlignedSequenceContainer(alphabet_);
            if (probProfileType_ == "none"){
                DLOG(INFO) << "[PIP ASR] Joint ancestral sequence reconstruction is called without probability profile.";
                asc = getAncestralSequences();
            }
            else if (probProfileType_ != "none") {
                // in the of probability profile we call the method with probability profile.
                DLOG(INFO) << "[PIP ASR] Joint ancestral sequence reconstruction is called with probability profile.";
                std::tie(asc,probabilityProfile_) = getAncestralSequencesWithProbability();
            }
            return asc;
        }



        void recursiveComputeLikelihoodAtSite(const Node *node, const TreeTemplate<Node> *tree,
                                              VVdouble &likelihoodArray, VVdouble &characterArray,
                                              const size_t siteNumber, const int newRootId) const;

        /****************************************** Base class functions **********************************************/

        std::map<int, std::vector<size_t>> getAllAncestralStates() const override { std::map<int, std::vector<size_t>>  ret_val; return ret_val;};
        // We can not use this definition because of "std::vector<size_t>": for gap state we have -1 which need int type.

        /**
        * @brief Get the ancestral sequence for a given node.
        *
        * This method is mainly for efficient internal use in other classes.
        * Consider using the recursiveGetAncestralStatesForAllNodes() method for a more general output.
        *
        * @param nodeId The id of the node at which the sequence must be reconstructed.
        * @return A sequence object.
        */
        Sequence *getAncestralSequenceForNode(int nodeId) const;

        /**
        * @brief Get ancestral states for a given node as a vector of int.
        *
        * The size of the vector depends on the implementation.
        * This method is mainly for efficient internal use in other classes.
        * Consider using the getAncestralSequenceForNode() method for a more
        * general output.
        *
        * @param nodeId the id of the node at which the states must be reconstructed.
        * @return A vector of states indices.
        * @see getAncestralSequenceForNode
        */
        std::vector<size_t> getAncestralStatesForNode(int nodeId) const final {std::vector<size_t> ret_val; return ret_val;};// todo: Marginal should be implemented and this bug should be fixed.
        // We can not use this definition cause "std::vector<size_t>": for gap state we have -1 which need int type.



        /**
        * @brief Get Joint ancestral sequences for all nodes including gap (aka joint ancestral sequence reconstruction using DP approach with indel).
        *
        * @note This method considers gap as an additional character.
        *
        * Call the getAllAncestralStates() method and map the sites to the original one.
        * @return A new SiteContainer object containing joint ancestral sequences.
        */
        AlignedSequenceContainer *getAncestralSequences() const ;


        /**
        * @brief Get Joint ancestral sequences for all nodes including gap (aka joint ancestral sequence reconstruction using DP approach with indel).
        *
        * @note This method considers gap as an additional character.
        *
        * Call the getAllAncestralStates() method and map the sites to the original one.
        * @return A new SiteContainer object containing joint ancestral sequences.
        */
        std::tuple<AlignedSequenceContainer *, VVVdouble> getAncestralSequencesWithProbability() const ;


        /**
         * @brief Get the probability profile for a given all sites.
         */
         VVVdouble getProbabilityProfileForAllSites() const;

        /**
         * @brief Set the probability profile for a given site.
         */
        void computeProbabilityProfileForAllSites(VVVdouble &likelihoodArray, VVVdouble &ProbabilityProfile) const;



        /**
         * @brief Get the tree associated to this object.
         * @return A pointer toward the tree associated to this object.
         */
        TreeTemplate<Node> getTree() const{ return tree_; }



    private:

        /**
         *  This method is called when the joint ancestral state reconstruction is called by getAllAncestralStates().
         *
         *  This method is designed for internal usage @see getAllAncestralStates().
         *
         * @param node  The node at which ancestral is reconstructed.
         * @param tree  The subtree at which the computation is done.
         * @param ancestors A vector of states indices.
         * @param likelihoodArray  The likelihood matrix.
         * @param cLikelihoodArray  The character states corresponding to likelihood array.
         * @param data The sequences that the computation is done on.
         * @param siteNumber Site number that all the computation is done on.
         * @param newRootId The id of new root for the subtree we work on.
         * @param isVisited A flag to indicate whether the node is visited or not.
         */
        void recursiveJointAncestralStates( const Node *node, const TreeTemplate<Node> *tree,
                                            std::map<int, std::vector<int> > &ancestors,
                                            VVVdouble &likelihoodArray, VVVdouble &cLikelihoodArray,
                                            AlignedSequenceContainer &data, const size_t siteNumber,
                                            const int newRootId, bool isVisited) const;


    };
} // end of namespace bpp.
#endif //ARPIP_PIPANCESTRALSTATERECONSTRUCTION_H
