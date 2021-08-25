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
 * @file PIPIndelRateInference.h
 * @author Gholamhossein Jowkar
 * @date 11 06 2021
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

#ifndef ARPIP_PIPINDELRATEINFERENCE_H
#define ARPIP_PIPINDELRATEINFERENCE_H


// From bpp-arpip
#include "../Likelihood/PIPDRHomogeneousTreeLikelihood.h"


// From bpp-core
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Phyl/Model/Protein/WAG01.h>

namespace bpp {

    class OptimizationFunction:
            public virtual Function,
            public AbstractParametrizable{
    private:

        /**
         * @brief This contains the objective function we are trying to optimize.
         */
        long double fValue_;

        /**
        * @brief This contains the likelihood values and PIP parameters
        */
        mutable PIPDRHomogeneousTreeLikelihood *likeFunData_;

    public:
        /**
         * @brief constructor.
         */
         OptimizationFunction(PIPDRHomogeneousTreeLikelihood *likeFunObj);

        /**
         * @brief Copy constructor.
         */
        OptimizationFunction(const OptimizationFunction &obj);

        /**
         * @brief Operator overloading.
         */
        OptimizationFunction &operator= (const OptimizationFunction &obj);

        /**
         * @brief destructor.
         */
        virtual ~OptimizationFunction(){};

        OptimizationFunction *clone() const { return new OptimizationFunction(*this); }

    public:
        void setParameters(const ParameterList &pl);

        double getValue()const { return fValue_; }

        void fireParameterChanged(const ParameterList &pl);


    };

//-----------------------------------------------------------------------------------------------------------------------

    class PIPIndelRateInference {
    private:
        mutable double lambda_;
        mutable double mu_;
        mutable PIPDRHomogeneousTreeLikelihood *likelihood_;

    public:
        /**
         * @brief Constructor when the likelihood object is already there.
         *
         * @param likeFunObj object from the class implements the likelihood computation for a tree using PIP model.
         */
        PIPIndelRateInference(PIPDRHomogeneousTreeLikelihood &likeFunObj,
                              int MaxIteration = 100000,
                              double tolerance = 0.0001);

        /**
         * @brief Constructor which build the likelihood array inside itself
         *
         * @param likeFunObj object from the class implements the likelihood computation for a tree using PIP model.
         */
        PIPIndelRateInference(const SiteContainer &sites,
                              const TreeTemplate<Node> &ttree,
                              std::map<std::string,
                              std::string> modelMap,
                              double *lambda,
                              double *mu,
                              int MaxIteration = 100000,
                              double tolerance = 0.0001);

        /**
         * @brief Copy constructor.
         */
        PIPIndelRateInference(const PIPIndelRateInference &inf);

        /**
         * @brief Operator overloading.
         */
        PIPIndelRateInference &operator=(const PIPIndelRateInference &inference);

        /**
         * @brief Destructor.
         */
//        virtual ~PIPIndelRateInference();

        /**
        * @brief Method called by constructor
        */
        void
        init_(const SiteContainer &sites, const TreeTemplate<Node> &ttree, std::map<std::string, std::string> modelMap);

        /**
        * @brief Build a new PIPInferenceIndelRates object and infer the corresponding parameters.
        *
        * Inference of the @param lambda and @param mu using method described in paper.
        *
        * @param likelihood The instance of class that store and compute the tree likelihood using PIP.
        * @param lambda The evolutionary parameter we are inferring.
        * @param mu The evolutionary parameter we are inferring.
        * @throw Exception in an error occured.
        */
        void inferIndelRateFromSequences(PIPDRHomogeneousTreeLikelihood &likelihood,
                                         int maxIteration,
                                         double tolerance);

        /***************************** Setters and Getters ***************************************************/

        double getLambda() const { return lambda_; }

        void setLambda(double lambda) { lambda_ = lambda; }

        double getMu() const { return mu_; }

        void setMu(double mu) { mu_ = mu; }

        PIPDRHomogeneousTreeLikelihood *getLikelihood() const { return likelihood_; }

        void setLikelihood(PIPDRHomogeneousTreeLikelihood *likelihood) { likelihood_ = likelihood; }


    };

}// end of namespace
#endif //ARPIP_PIPINDELRATEINFERENCE_H
