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
 * @file PIP13.h
 * @author Gholamhossein Jowkar
 * @date 01 07 2020
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

#ifndef _PIP13_H_
#define _PIP13_H_

#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include "Bpp/Phyl/Model/AbstractSubstitutionModel.h"
#include "Bpp/Phyl/Model/Nucleotide/NucleotideSubstitutionModel.h"
#include "Bpp/Phyl/Model/Protein/ProteinSubstitutionModel.h"
#include "Bpp/Phyl/Model/Codon/CodonSubstitutionModel.h"

namespace bpp{

    /**
 * @brief The PIP substitution model with gap characters.
 *
 * This model expends any reversible substitution model with gaps as an additional state.
 * Although the conditionnal subtitution process is reversible, the insertion/deletion process
 * needs not be. The model hence adds one parameter just for deletions, @f$\mu@f$.
 * If we note @f$Q@f$ the (simple) transition matrix (= Markov generator) and @f$Q^\epsilon@f$ the extended one, we have:
 * @f[
 * Q^\epsilon =
 * \left(
 * \begin{array}{ccc|c}
 * & & & \mu \\
 * & \rule{0cm}{1cm}\rule{1cm}{0cm}Q-\mu\delta_{ij} & & \vdots \\
 * & & & \mu \\
 * \hline
 * \lambda\pi_1 & \ldots & \lambda\pi_n & -\lambda \\
 * \end{array}
 * \right)
 * @f]
 * where @f$n@f$ is the number of states of the simple model (in most case equal to the size of the alphabet) and @f$(\pi_1,\ldots,\pi_n)@f$ is the vector of equilibrium frequencies of the conditional model.
 * @f$\delta_{ij}@f$ is 1 if i=j, 0 otherwise.
 * Note that in the original paper @f$Q@f$ is noted as @f$R@f$, and @f$Q_t@f$ is used for the probability matrix, which is referred here as @f$P^\epsilon(t)@f$ for consistency with the documentation of other models.
 *
 * The extended Markov model is reversible, and the equilibrium frequencies are
 * @f[
 * \pi^\epsilon = \left( \pi \cdot \frac{\lambda}{\lambda + \mu}, \frac{\mu}{\lambda + \mu}\right).
 * @f]
 * The corresponding exchangeability matrix is:
 * @f[
 * S^\epsilon =
 * \left(
 * \begin{array}{ccc|c}
 * & & & \lambda + \mu \\
 * & \rule{0cm}{1cm}\rule{1cm}{0cm}(S - \frac{\mu\delta_{ij}}{\pi_i})\frac{\lambda+\mu}{\lambda} & & \vdots \\
 * & & & \lambda + \mu \\
 * \hline
 * \lambda + \mu & \ldots & \lambda + \mu & - (\lambda + \mu) \\
 * \end{array}
 * \right)
 * @f]
 * The eigen values and vectors are computed numerically, but the transition probabilities are computed analytically from the simple substitution model.
 *
 * Reference:
 * - Bouchard-Côté A, Jordan MI. Evolutionary inference via the Poisson Indel Process. Proc NatlAcad Sci USA. 2013;110(4):1160-1166. doi:10.1073/pnas.1220450110 -.
 */



    class PIP13: public AbstractReversibleSubstitutionModel {
    private:
        std::unique_ptr<ReversibleSubstitutionModel> simpleModel_;
        RowMatrix<double> simpleGenerator_;
        RowMatrix<double> simpleExchangeabilities_;
        mutable double exp_;
        mutable RowMatrix<double> p_;
        double mu_;
        std::string nestedPrefix_;

    public:
        /**
        * @brief Build a new PIP model from a standard substitution model.
        *
        * The alphabet and number of states for the extended model will be derived from the simple one.
        *
        * @param simpleModel The simple model to use to build the extended one.
        * THE PIP13 class will own the simple one, meaning that it will be destroyed together with the PIP13 instance, and cloned when cloning the PIP13 instance.
        * To prevent the original simple model to be destroyed, you should make a copy of it before creating the PIP13 instance.
        * @param mu     Deletion rate.
        */
        PIP13(ReversibleSubstitutionModel* simpleModel, double mu = 0.1);

        PIP13(const PIP13& model):
        AbstractParameterAliasable(model),
        AbstractReversibleSubstitutionModel(model),
        simpleModel_(dynamic_cast<ReversibleSubstitutionModel*>(model.simpleModel_->clone())),
        simpleGenerator_(model.simpleGenerator_),
        simpleExchangeabilities_(model.simpleExchangeabilities_),
        exp_(model.exp_),
        p_(model.p_),
        mu_(model.mu_),
        nestedPrefix_(model.nestedPrefix_)
        {}

        PIP13& operator=(const PIP13& model)
        {
            AbstractParameterAliasable::operator=(model);
            AbstractSubstitutionModel::operator=(model);
            AbstractReversibleSubstitutionModel::operator=(model);
            simpleModel_.reset(dynamic_cast<ReversibleSubstitutionModel*>(model.simpleModel_->clone()));
            simpleGenerator_         = model.simpleGenerator_;
            simpleExchangeabilities_ = model.simpleExchangeabilities_;
            exp_                     = model.exp_;
            p_                       = model.p_;
            mu_                      = model.mu_;
            nestedPrefix_            = model.nestedPrefix_;
            return *this;
        }

        virtual ~PIP13(){}

        PIP13* clone() const { return new PIP13(*this); }

        /**
         *@brief All probabilities of change from state i to state j during time t.
         * @param i starting state
         * @param j ending state
         * @param d distance or phylogeny time
         * @return probabilities of change from state i to state j
         */
        double Pij_t(size_t i, size_t j, double d) const;

        const Matrix<double> &getPij_t(double d) const;

        std::string getName() const{ return "PIP13"; }


        /**
        * @brief This method is forwarded to the simple model.
        *
        * @param data The data to be passed to the simple model (gaps will be ignored).
        * @param pseudoCount A (typically small) value to add to each count to avoid 0 estimates.
        */
        void setFreqFromData(const SequenceContainer &data, double pseudoCount = 0) {
            simpleModel_->setFreqFromData(data, pseudoCount);
        }

        void fireParameterChanged(const ParameterList &parameters) {
            AbstractParameterAliasable::fireParameterChanged(parameters);
            simpleModel_->matchParametersValues(parameters);
            mu_ = getParameter_(0).getValue();
            updateMatrices();
        }

        /**
         * @brief Get the number of states. For most models, this equals the size of the alphabet.
         * @return     The number of different states in the model.
         */
        size_t getNumberOfStates() const{ return size_; }

        /**
         * @brief This method is used to initialize likelihoods in recursions.
         * It typically sends 1 if i = state, 0 otherwise, where i is one of the possible states
         * of the alphabet allowed in the model and state is the observed state in the considered
         * sequence/site.
         *
         * @param the index of the state in the model.
         * @param state An observed state in the sequence/site.
         * @return 1 or 0 depending if the two states are compatible.
         */
        double getInitValue(size_t i, int state) const;

        /**
         * @brief  Overloading the default function: This method is used to initialize likelihoods in recursions.
         * For each state i: P_ij (t_y), where t_y is the branch length between y and its father.
         *
         * @param the index of the state in the model.
         * @param state An observed state in the sequence/site.
         * @param branchLength
         * @return Transition probability from observed state to all existing states.
         *
         * @see [Pupko, et al., 2000, Molecular biology and evolution]
         */
        using TransitionModel::getInitValue;
        double getInitValue(size_t i, int state, double branchLength) const;

        void setNamespace(const std::string &prefix);

        const ReversibleSubstitutionModel* getNestedModel() const { return simpleModel_.get(); }

    protected:

        void updateMatrices();


        void setParametersValues(const ParameterList &parameters);
    };

} //end of namespace bpp.
#endif // _PIP13_H_