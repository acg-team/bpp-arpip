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
 * @file PIPDRHomogeneousTreeLikelihood.h
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

#ifndef ARPIP_PIPDRHOMOGENEOUSTREELIKELIHOOD_H
#define ARPIP_PIPDRHOMOGENEOUSTREELIKELIHOOD_H

// From bpp-phyl
#include "Bpp/Phyl/Likelihood/AbstractHomogeneousTreeLikelihood.h"
#include "Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h"

// From bpp-arpip
#include "PIPDRTreeLikelihood.h"
#include "PIPDRTreeLikelihoodData.h"
#include "PIPParameterData.h"


// From bpp-core
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>

// From Std
#include <math.h>
#include <glog/logging.h>
#include <numeric>

namespace bpp {

    /**
     * @brief This class implements the likelihood computation for a tree using PIP model
     *
     * <p> The substitution model is the same over the tree (homogeneous model).
     *
     *  A uniform distribution of rates among the sites is allowed.
     *
     * This class uses an instance of the PIPDRTreeLikelihoodData for conditional likelihood storage.
     *
     * All nodes share the same site patterns.</p>
     */

    class PIPDRHomogeneousTreeLikelihood:
            public AbstractHomogeneousTreeLikelihood,
            public PIPDRTreeLikelihood{

    private:
        /**
         * @brief This contains the likelihood values and PIP parameters
         */
        mutable PIPDRTreeLikelihoodData *likelihoodData_;

        /**
         * @brief This contains the sites
         */
        mutable SiteContainer *shrunkData_;

        /**
         * @brief This is the flag to indicate we are using PIP model.
         */
        bool isPIP_;


        /**
        * @brief This contains all $\tilde{f}_v(n)$ Felsenstein's temporary vector of values used for likelihood computation.
        *
        * <pre>
        * x[n][i][s]
        *   |------------> Id of node n
        *     |---------> Site i
        *            |---> State s
        * </pre>
        * We call this the <i> $\tilde{f}_v$ Felsenstein's tmp array</i> for each node.
        */
        mutable std::map<int, VVdouble> nodeFTildeValues_;

        /**
        * @brief This contains all $\tilde{f}_v$ Felsenstein's temporary values used for likelihood computation.
        *
        * <pre>
        * x[n][i][s]
        *   |------------> Id of node n
        *     |---------> Site i
        * </pre>
        * We call this the <i> $\tilde{f}_v$ Felsenstein's tmp array</i> for each node.
        *  $\tilde{f}_v = \sum_{n \in \Sigma_\e} \pi^\e_n \tilde{f}_v(n)$
        */
        mutable std::map<int, Vdouble> nodeSumFTilde_;

        /**
        * @brief This contains all $f_v$ Felsenstein's values used for likelihood computation.
        *
        * <pre>
        * x[n][i][s]
        *   |------------> Id of node n
        *     |---------> Site i
        * </pre>
        * We call this the <i> $f_v$ Felsenstein's array</i> for each node.
        * For an internal node $v$, the probability of homology path $f_v$ is equal to $\tilde{f}_v$ for root
        * and for $v \in \mathcal{A}$ is equal to $\beta_v * \tiled{f}_v$ and for the rest is equal to zero.
        */
        mutable std::map<int, Vdouble> nodeFValues_;

        /**
         * @brief Likelihood values of non-empty columns
         */
        mutable Vdouble p_c_;

        /**
         * @brief Likelihood value of single empty column
         */
        mutable double p_c_0_;

        /**
        * @brief This contains all Iota values used for computation.
        *
        *
        * We call this the <i>Nu array</i> for each node.
        *
	    $\varphi = \frac{ \nu^{|\sequence|} e^{(p(c_\emptyset)-1)\nu}}{|\sequence|!}$\;
        */
        mutable long double phi_;

    protected:
        double minusLogLik_;

    public:
        /**
        * @brief Build a new PIPDRHomogeneousTreeLikelihood object without data.
        *
        * This constructor only initialize the parameters.
        * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
        *
        * @param tree The tree to use.
        * @param data Sequences to use.
        * @param model The substitution model to use.
        * @param rDist The rate across sites distribution to use.
        * @param checkRooted Tell if we have to see if the tree is rooted or not.
        * @param verbose Should I display some info?
        * @throw Exception in an error occurred.
        */
        PIPDRHomogeneousTreeLikelihood(
                const Tree& tree,
                TransitionModel* model,
                DiscreteDistribution* rDist,
                bool checkRooted = true,
                bool verbose = true);


        /**
        * @brief Build a new PIPDRHomogeneousTreeLikelihood object and compute the corresponding likelihood.
        *
        * This constructor initializes all parameters, data, and likelihood arrays.
        *
        * @param tree The tree to use.
        * @param data Sequences to use.
        * @param model The substitution Model to use.
        * @param rDist The rate across sites distribution to use.
        * @param checkRooted Tell if we have to see if the tree is rooted or not.
        * @param verbose Should I display some info?
        * @throw Exception in an error occured.
        */
        PIPDRHomogeneousTreeLikelihood(
                const Tree &tree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                bool checkRooted = true,
                bool verbose = true
        );

        /**
        * @brief Build a new PIPDRHomogeneousTreeLikelihood object and compute the corresponding likelihood.
        *
        * This constructor initializes all parameters, data, and likelihood arrays.
        *
        * @param tree The tree to use.
        * @param data Sequences to use.
        * @param model The substitution Model to use.
        * @param rDist The rate across sites distribution to use.
        * @param checkRooted Tell if we have to see if the tree is rooted or not.
        * @param verbose Should I display some info?
        * @throw Exception in an error occured.
        */
        PIPDRHomogeneousTreeLikelihood(
                const Tree &tree,
                const SiteContainer &data,
                TransitionModel *model,
                DiscreteDistribution *rDist,
                const double mu,
                const double lambda,
                bool checkRooted = true,
                bool verbose = true
        );
        /**
         * @brief Copy constructor.
         */
        PIPDRHomogeneousTreeLikelihood(const PIPDRHomogeneousTreeLikelihood &lik);

        /**
          * @brief  Operator overloading.
          */
        PIPDRHomogeneousTreeLikelihood &operator=(const PIPDRHomogeneousTreeLikelihood &lik);

        /**
          * @brief destructor.
          */
        virtual ~PIPDRHomogeneousTreeLikelihood(){ delete likelihoodData_; };

        PIPDRHomogeneousTreeLikelihood *clone() const { return new PIPDRHomogeneousTreeLikelihood(*this); }

    private:
        /**
        * @brief Method called by constructor.
        */
        void init_();

        /**
        * @brief Method called by constructor for PIP model.
        */
        void init_(const double mu, const double lambda);

    public:

        /**
         * @name The TreeLikelihood interface.
         *
         * Other methods are implemented in the AbstractTreeLikelihood class.
         *
         * @{
         */
        void setData(const SiteContainer &sites);
        double getLikelihood () const;
        double getLogLikelihood() const;
        double getLikelihoodForASite (size_t site) const;
        double getLogLikelihoodForASite(size_t site) const;
        size_t getSiteIndex(size_t site) const { return likelihoodData_->getRootArrayPosition(site); }
        /** @} */

        /**
         * @name The DiscreteRatesAcrossSites interface implementation:
         *
         * @{
         */
        double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
        double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
        double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
        double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
        /** @} */


        /**
         * @brief Implements the Function interface.
         *
         * Update the parameter list and call the fireParameterChanged() method.
         *
         * If a subset of the whole parameter list is passed to the function,
         * only these parameters are updated and the other remain constant (i.e.
         * equal to their last value).
         *
         * @param parameters The parameter list to pass to the function.
         */
        void setParameters(const ParameterList &parameters);

        /**
         * @brief Function and NNISearchable interface.
         */
        double getValue() const;

        /**
        * @name DerivableFirstOrder interface.
        *
        * @{
        */
        double getFirstOrderDerivative(const std::string& variable) const{ return 0; } // Not implemented for now.
        /** @{ */

        /**
         * @name UPKODerivableSecondOrder interface.
         *
         * @{
         */
        double getSecondOrderDerivative(const std::string& variable) const{ return 0; } // Not implemented for now.
        double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const { return 0; } // Not implemented for now.
        /** @} */

        virtual void computeTreeDLikelihoodAtNode(const Node* node) { std::cout << "not implemented yet!"; }
        virtual void computeTreeDLikelihoods() { std::cout << "not implemented yet!"; }

        virtual void computeTreeD2LikelihoodAtNode(const Node* node) { std::cout << "not implemented yet!"; }
        virtual void computeTreeD2Likelihoods() { std::cout << "not implemented yet!"; }



        virtual void fireParameterChanged(const ParameterList& params);

        virtual void resetLikelihoodArrays(const Node* node);

        /**
         * @brief This method is mainly for debugging purpose.
         *
         * @param node The node at which likelihood values must be displayed.
         */
        virtual void displayLikelihood(const Node* node);


        /**
         * @brief This method is compute the tree likelihood.
         */
        void computeTreeLikelihood();



        /**
        * @name PIP specific likelihood functions.
        * @author Gholamhossein Jowkar
        * Other methods are implemented in the private part of this section.
        *
        * @{
        */
        virtual long double computePIPTreeLikelihood_(const double lambda, const double mu);
        void firePIPParameterChanged_(const double lambda, const double mu);

        /** @} */


        /**
         * @name PIPDRTreeLikelihood interface.
         *
         * @{
         */
        PIPDRTreeLikelihoodData *getLikelihoodData() { return likelihoodData_; }
        const PIPDRTreeLikelihoodData *getLikelihoodData() const  { return likelihoodData_; }

        virtual void computeLikelihoodAtSite(
                const Node *node,
                const TreeTemplate<Node> *tree,
                VVdouble &likelihoodArray,
                VVdouble &characterArray,
                size_t siteNumber,
                const int newRootId) const {
            computeLikelihoodAtSite_(node, likelihoodArray, characterArray, siteNumber, newRootId);
        }
        /** @} */

    protected:

        /**
         * @brief Marginal likelihood of single empty column (non-observable): p(c(∅)).
         *
         * @param node The node defining the subtree to analyse (Usually it is the root).
         *
         * @return p_c_0 Likelihood of a single column full of gaps.
         */
        void computePIPLikelihoodForEmptySiteAtNode_(const Node *node);

        /**
         * @brief Marginal likelihood of non-observable empty columns.
         *
         * Marginal likelihood of non-observable empty columns for an alignment of length |m| where p(c(∅),|m| ) is the
         * likelihood of a single MSA column full of gaps.
         *
         * @note the p_c_0 and nu have to be computed in advance.
         */
        void computePIPPhi_();

        /**
         * @brief Marginal likelihood of observable columns.
         *
         * @param node The node defining the subtree to analyse (Usually it is the root).
         */
        void computePIPLikelihoodForNonEmptySitesAtNode_(const Node *node);

        /**
         * @brief Reset all the PIP arrays.
         *
         * Usually used for PIP param optimization.
         * @param node The node defining the subtree to analyse (Usually it is the root).
         */
        void resetFVArrays_(const Node *node);


        virtual void
        computeLikelihoodAtSite_(
                const Node *node,
                VVdouble &likelihoodArray,
                VVdouble &characterArray,
                size_t siteNumber,
                const int newRootId,
                const Node *sonNode = 0) const;

        virtual void computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode = 0) const;


        /**
         * @brief Compute conditional likelihoods using Pupko et al..
         *
         * This method is the "core" likelihood computation function, performing all the product upon all nodes, the summation for each ancestral state and each rate class.
         * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
         * Use with care!
         *
         * @param iLik A vector of likelihood arrays, one for each conditional node.
         * @param tProb A vector of transition probabilities, one for each node.
         * @param oLik The likelihood array to store the computed likelihoods.
         * @param cLik A vector of candidate char, one for each conditional node.
         * @param nbNodes The number of nodes = the size of the input vectors.
         * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
         * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
         * @param nbStates The number of states (the third dimension of the likelihood array).
         * If true, the resetLikelihoodArray method will be called.
         */
        static void computeLikelihoodFromArraysForRoot(
                const std::vector<const Vdouble*>& iLik,
                const std::vector<const VVVdouble*>& tProb,
                const std::vector<const Vdouble*>& cLik,
                const Vdouble *p,
                Vdouble &oLik,
                Vdouble &oClik,
                size_t siteNb,
                size_t nbNodes,
                size_t nbDistinctSites,
                size_t nbClasses,
                size_t nbStates);


        /**
         * @brief Compute conditional likelihoods.
         *
         * This method is the "core" likelihood computation function, performing all the product upon all nodes, the summation for each ancestral state and each rate class.
         * This function is specific to non-reversible models: the subtree containing the root is specified separately.
         * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
         * Use with care!
         *
         * @param iLik A vector of likelihood arrays, one for each conditional node.
         * @param tProb A vector of transition probabilities, one for each node.
         * @param iLikR The likelihood array for the subtree containing the root of the tree.
         * @param tProbR The transition probabilities for thr subtree containing the root of the tree.
         * @param oLik The likelihood array to store the computed likelihoods.
         * @param oClik The candidate ancestral char array to store the computed candidates.
         * @param nbNodes The number of nodes = the size of the input vectors.
         * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
         * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
         * @param nbStates The number of states (the third dimension of the likelihood array).
         * If true, the resetLikelihoodArray method will be called.
         */
        static void computeLikelihoodFromArrays(
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
                size_t nbStates);







    };
}// end of namespace bpp.

#endif //ARPIP_PIPDRHOMOGENEOUSTREELIKELIHOOD_H
