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
 * @file PIPDRTreeLikelihood.h
 * @author Gholamhossein Jowkar
 * @date 23 07 2020
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


#ifndef ARPIP_PIPDRTREELIKELIHOOD_H
#define ARPIP_PIPDRTREELIKELIHOOD_H

// From bpp-phyl
#include "Bpp/Phyl/Likelihood/AbstractNonHomogeneousTreeLikelihood.h"

// From bpp-arpip
#include "PIPDRTreeLikelihoodData.h"


// From bpp-core
#include "Bpp/Numeric/VectorTools.h"
#include "Bpp/Numeric/Prob/DiscreteDistribution.h"

namespace bpp{

    /**
        * @brief Reimplementation of the interface for double-recursive (DR) implementation of the likelihood computation "under PIP".
        *
        * In the DR implementation, three conditional likelihoods are stored at each node, corresponding to the three connected subtrees
        * The DR implementation hence uses 3x more memory than the simple recursive (R) implementation.
        * However, the likelihood of the tree can be computed independently at each node of the tree,
        * which is very convenient for topology estimation (the likelihood change of each NNI movement can be computed directly),
        * ancestral state reconstruction (all marginal ancestral states can be computed in one pass over the tree), or for
        * substitution mapping.
        *
        * This interface provides
        * - a method to access the DR likelihood data structure,
        * - a method to compute the likelihood array at each node.
        *
        * For now, this interface inherits from DiscreteRatesAcrossSitesTreeLikelihood and not TreeLikelihood,
        * since the data structure available accounts for rate across site variation.
        * This may change in the future.
        *
        * @see DRTreeLikelihoodTools
        * @see DRTreeLikelihood
        */
    class PIPDRTreeLikelihood:
            public virtual DiscreteRatesAcrossSitesTreeLikelihood{
    public:
        PIPDRTreeLikelihood() {}
        virtual ~PIPDRTreeLikelihood(){}

        PIPDRTreeLikelihood *clone() const = 0;

        /**
         * @name Get the likelihood data structure associated to this class.
         *
         * @{
         */
        virtual PIPDRTreeLikelihoodData *getLikelihoodData() = 0;
        virtual const PIPDRTreeLikelihoodData *getLikelihoodData() const = 0;
        /** @} */


        /**
         * @brief Compute the likelihood array at a given site.
         *
         * @param node The node of the tree to be considered as root of subtree.
         * @param tree The tree.
         * @param likelihoodArray The array where to store the results.
         * @param characterArray The array where to store the candidate ancestral values.
         * @param siteNumber The number of site which we operate on.
         * @param newRootId The id of new root which the MLIndelPoints algorithm returns.
         */

        virtual void
        computeLikelihoodAtSite(const Node *node, const TreeTemplate<Node> *tree, VVdouble &likelihoodArray,
                                     VVdouble &characterArray, size_t siteNumber, const int newRootId) const = 0;






    };


}// end namespace bpp
#endif //ARPIP_PIPDRTREELIKELIHOOD_H
