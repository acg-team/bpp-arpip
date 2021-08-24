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
 * @file PIP13.cpp
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


#include "PIP13.h"

using namespace bpp;

#include <cmath>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

PIP13::PIP13(ReversibleSubstitutionModel *simpleModel, const double mu) :
        AbstractParameterAliasable("PIP13"),
        AbstractReversibleSubstitutionModel(simpleModel->getAlphabet(),
                new CanonicalStateMap(simpleModel->getStateMap(), true), "PIP13"),
        simpleModel_(simpleModel),
        simpleGenerator_(),
        simpleExchangeabilities_(),
        exp_(), p_(), mu_(mu),
        nestedPrefix_("model_" + simpleModel->getNamespace()) {

    addParameter_(new Parameter("PIP13.mu", mu, &Parameter::R_PLUS));
    simpleModel_->setNamespace("PIP13." + nestedPrefix_);
    addParameters_(simpleModel->getParameters());


    //We need to override this from the AbstractSubstitutionModel constructor,
    //since the number of states in the model is no longer equal to the size of the alphabet.
    /*****/
    size_ = simpleModel->getNumberOfStates() + 1;
    generator_.resize(size_, size_);
    exchangeability_.resize(size_, size_);
    freq_.resize(size_);
    eigenValues_.resize(size_);
    leftEigenVectors_.resize(size_, size_);
    rightEigenVectors_.resize(size_, size_);
    p_.resize(size_, size_);
    updateMatrices();
    /*****/
}

/******************************************************************************/
void PIP13::updateMatrices() {
//    double f = (mu_ == 0) ? 1 : mu_;
    for(size_t i = 0; i < size_ - 1; i++)
        freq_[i] = simpleModel_->freq(i);
    //The extended Markov model is reversible, and the equilibrium frequencies: freq = [freq 0];
    freq_[size_ - 1] = 0;

    //The generator matrix $Q$ of the model
    simpleGenerator_ = simpleModel_->getGenerator();
    simpleExchangeabilities_ = simpleModel_->getExchangeabilityMatrix();

    // Generator and exchangeabilities:
    for (size_t i = 0; i < size_ - 1; ++i) {
        for (size_t j = 0; j < size_ - 1; ++j) {
            generator_(i, j) = simpleGenerator_(i, j);
            exchangeability_(i, j) = simpleExchangeabilities_(i, j);
            if (i == j) {
                generator_(i, j) -= mu_;
                exchangeability_(i, j) = simpleExchangeabilities_(i, j) - mu_;
            }
        }
        generator_(i, size_ - 1) = mu_;
        generator_(size_ - 1, i) = 0;
        exchangeability_(i, size_ - 1) = mu_;
//        exchangeability_(size_ - 1, i) = mu_;
    }
    generator_(size_-1, size_-1) = 0;
//    exchangeability_(size_ - 1, size_ - 1) = -0;

    // Normalization:
//    setDiagonal();

    // Compute eigen values and vectors:
    // In this version of bpp AbstractSubstitutionModel::updateMatrices() is not working.
    if (enableEigenDecomposition()) {
        EigenValue<double> ev(generator_);
        rightEigenVectors_ = ev.getV();
        eigenValues_ = ev.getRealEigenValues();
        iEigenValues_ = ev.getImagEigenValues();
        try {
            MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
            isNonSingular_ = true;
            isDiagonalizable_ = true;
            for (size_t i = 0; i < size_ && isDiagonalizable_; i++) {
                if (abs(iEigenValues_[i]) > NumConstants::TINY())
                    isDiagonalizable_ = false;
            }
        }
        catch (ZeroDivisionException &e) {
            ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

            isNonSingular_ = false;
            isDiagonalizable_ = false;
            MatrixTools::Taylor(generator_, 30, vPowGen_);
        }
    }


    //It is very likely that we are able to compute the eigen values and vector from the one of the simple model.
    //For now however, we will use a numerical diagonalization:
//    AbstractSubstitutionModel::updateMatrices();
    //We do not use the one from  AbstractReversibleSubstitutionModel, since we already computed the generator.
    // With this version of bpp lib, the AbstractSubstitutionModel::updateMatrices() is not working correctly,
    // since it is unable to find eigenvector for eigenvalue 0. Taylor series used instead.

}

/******************************************************************************/
////This function is implemented according to ARPIPv1 and is quiet different from RE08.
double PIP13::Pij_t(size_t i, std::size_t j, double d) const {
    return  getPij_t(d) (i, j);
}

/******************************************************************************/
const Matrix<double>& PIP13::getPij_t(double d) const {

    if (d == 0) {
        MatrixTools::getId(size_, pijt_);
    } else if (isNonSingular_) {
        if (isDiagonalizable_) {
            MatrixTools::mult<double>(rightEigenVectors_, VectorTools::exp(eigenValues_ * (rate_ * d)),
                                      leftEigenVectors_, pijt_);
        } else {
            std::vector<double> vdia(size_);
            std::vector<double> vup(size_ - 1);
            std::vector<double> vlo(size_ - 1);
            double c = 0, s = 0;
            double l = rate_ * d;
            for (size_t i = 0; i < size_; i++) {
                vdia[i] = std::exp(eigenValues_[i] * l);
                if (iEigenValues_[i] != 0) {
                    s = std::sin(iEigenValues_[i] * l);
                    c = std::cos(iEigenValues_[i] * l);
                    vup[i] = vdia[i] * s;
                    vlo[i] = -vup[i];
                    vdia[i] *= c;
                    vdia[i + 1] = vdia[i]; // trick to avoid computation
                    i++;
                } else {
                    if (i < size_ - 1) {
                        vup[i] = 0;
                        vlo[i] = 0;
                    }
                }
            }
            MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, pijt_);
        }
    } else {
        MatrixTools::getId(size_, pijt_);
        double s = 1.0;
        double v = rate_ * d;
        size_t m = 0;
        while (v > 0.5)    // exp(r*t*A)=(exp(r*d/(2^m) A))^(2^m)
        {
            m += 1;
            v /= 2;
        }
        for (size_t i = 1; i < vPowGen_.size(); i++) {
            s *= v / static_cast<double>(i);
            MatrixTools::add(pijt_, s, vPowGen_[i]);
        }
        while (m > 0)  // recover the 2^m
        {
            MatrixTools::mult(pijt_, pijt_, tmpMat_);
            MatrixTools::copy(tmpMat_, pijt_);
            m--;
        }
    }
//  MatrixTools::print(pijt_);
    return pijt_;
}
/******************************************************************************/

double PIP13::getInitValue(size_t i, int state) const
{
    if (i >= size_) throw IndexOutOfBoundsException("PIP13::getInitValue", i, 0, size_ - 1);
    if (state < -1 || !getAlphabet()->isIntInAlphabet(state))
        throw BadIntException(state, "PIP13::getInitValue. Character " + getAlphabet()->intToChar(state) + " is not allowed in model.");
    if (i == size_ - 1 && state == -1) return 1.;
    vector<int> states = getAlphabet()->getAlias(state);
    for (size_t j = 0; j < states.size(); j++)
        if ((int)i == states[j])
            return 1.;
    return 0.;
}

/******************************************************************************/

double PIP13::getInitValue(size_t i, int state, double branchLength) const {
    if (i >= size_) throw IndexOutOfBoundsException("PIP13::getInitValue", i, 0, size_ - 1);
    if (state < -1 || !getAlphabet()->isIntInAlphabet(state))
        throw BadIntException(state, "PIP13::getInitValue. Character " + getAlphabet()->intToChar(state) + " is not allowed in model.");
    if (i == size_ - 1 && state == -1) return 1.;

    vector<int> states = getAlphabet()->getAlias(state);

    for (size_t j = 0; j < states.size(); j++)
    if ((int)i == states[j])
        return 1.;
    return 0.;

}


/******************************************************************************/

void PIP13::setNamespace(const string &prefix)
{
    AbstractSubstitutionModel::setNamespace(prefix);
    //We also need to update the namespace of the nested model:
    simpleModel_->setNamespace(prefix + nestedPrefix_);
}

/******************************************************************************/

void PIP13::setParametersValues(const ParameterList &parameters){
    mu_ = parameters.getParameterValue("mu");
    updateMatrices();
}






