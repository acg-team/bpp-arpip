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
 * @file PIPIndelRateInference.cpp
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

// From ARPIP
#include "PIPIndelRateInference.h"
#include "../Likelihood/PIPDRHomogeneousTreeLikelihood.h"
#include "../SubstitutionModel/PIP13.h"


// From bpp-core
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Numeric/Function/DirectionFunction.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

// From bpp-phy
#include <Bpp/Phyl/Model/Protein/LG08.h>
#include <Bpp/Phyl/Model/Protein/LG08.h>
#include <Bpp/Phyl/Model/Protein/JCprot.h>
#include <Bpp/Phyl/Model/Protein/JTT92.h>
#include <Bpp/Phyl/Model/Protein/LG08.h>
#include <Bpp/Phyl/Model/Protein/WAG01.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>


using namespace bpp;

/******************************************************************************/
/**************************** OptimizationFunction ****************************/
/******************************************************************************/
/******************************** Constructor *********************************/
OptimizationFunction::OptimizationFunction(PIPDRHomogeneousTreeLikelihood *likeFunObj) :
        AbstractParametrizable(""), fValue_(0), likeFunData_(likeFunObj) {

            //We declare parameters here:
            addParameter_(new Parameter("lambda", bpp::RandomTools::giveRandomNumberBetweenZeroAndEntry(10),
                                        &Parameter::R_PLUS_STAR));
//            addParameter_(new Parameter("lambda", bpp::RandomTools::giveRandomNumberBetweenZeroAndEntry(10)));
            addParameter_(new Parameter("mu", bpp::RandomTools::giveRandomNumberBetweenZeroAndEntry(10),
                                new IntervalConstraint(0.0001, 100, false, true)));

            fireParameterChanged(getParameters());
        }

/******************************** Copy constructor *********************************/

OptimizationFunction::OptimizationFunction(const OptimizationFunction &obj) :
        AbstractParametrizable(""),
        fValue_(obj.fValue_),
        likeFunData_(obj.likeFunData_) {}

/****************************** Operator overloading ******************************/

OptimizationFunction &OptimizationFunction::operator=(const OptimizationFunction &obj) {
    AbstractParametrizable::operator=(obj);
    fValue_         = obj.fValue_;
    likeFunData_    = obj.likeFunData_;
    return *this;
}

/****************************** setParameters ******************************/

void OptimizationFunction::setParameters(const ParameterList &pl)
{
    matchParametersValues(pl);
}

/****************************** Operator overloading ******************************/

void OptimizationFunction::fireParameterChanged(const ParameterList &pl) {

    double lambda = getParameterValue("lambda");
    double mu = getParameterValue("mu");

    /*
     * We used -log likelihood to change maximization to minimization problem.
     */
    fValue_ = -(likeFunData_->computePIPTreeLikelihood(lambda, mu));
}


/******************************************************************************/
/**************************** PIPIndelRateInference ***************************/
/******************************************************************************/
/******************************** Constructor *********************************/

PIPIndelRateInference::PIPIndelRateInference(PIPDRHomogeneousTreeLikelihood &likeFunObj, int MaxIteration,
                                             double tolerance) :
        likelihood_(&likeFunObj), lambda_(0), mu_(0) {
    inferIndelRateFromSequences(*likelihood_, MaxIteration, tolerance);
}

/******************************** Constructor *********************************/

PIPIndelRateInference::PIPIndelRateInference(const SiteContainer &sites,
                                             const TreeTemplate<Node> &ttree,
                                             std::map<std::string, std::string> modelMap,
                                             double *lambda,
                                             double *mu,
                                             int MaxIteration,
                                             double tolerance) : lambda_(*lambda), mu_(*mu) {

    init_(sites, ttree, modelMap);
    inferIndelRateFromSequences(*likelihood_, MaxIteration, tolerance);

}

/******************************** Copy constructor *********************************/

PIPIndelRateInference::PIPIndelRateInference(const PIPIndelRateInference &inf) :
        likelihood_(inf.likelihood_), lambda_(inf.lambda_), mu_(inf.mu_) {}

/******************************* Operator overloading *******************************/

PIPIndelRateInference &PIPIndelRateInference::operator=(const PIPIndelRateInference &inference) {
    likelihood_     = inference.likelihood_;
    lambda_         = inference.lambda_;
    mu_             = inference.mu_;
    return *this;
}

/************************************** init_ ******************************************/

void PIPIndelRateInference::init_(const SiteContainer &sites,
                                  const TreeTemplate<Node> &ttree,
                                  std::map<std::string, std::string> modelMap) {

    // Genetic code
    unique_ptr<GeneticCode> gCode;

    // Extract alphabet to decide about the substitution model
    const Alphabet *inputAlphabet = sites.getAlphabet();

    string alphabetType = inputAlphabet->getAlphabetType();

    ReversibleSubstitutionModel *baseModel = nullptr;

    std::string baseModelName;

    std::map<std::string, std::string> baseModelMap;

    KeyvalTools::parseProcedure(modelMap["model"], baseModelName, baseModelMap);

    bool computeFrequenciesFromData = false;

    std::vector<std::string> keys;
    for (auto it = baseModelMap.begin(); it != baseModelMap.end(); ++it) keys.push_back(it->first);

    if (!keys.empty()) {
        baseModelName += "(";
        for (auto &key:keys) {
            if (key != "initFreqs") {
                baseModelName += key + "=" + baseModelMap[key];
            } else {
                if (baseModelMap[key] == "observed") {
                    computeFrequenciesFromData = true;
                }
            }
            baseModelName += ",";
        }
        baseModelName.pop_back();
        baseModelName += ")";
        modelMap["model"] = baseModelName;
    }

//    const ProteicAlphabet *alphabet = new ProteicAlphabet;
//    ReversibleSubstitutionModel *wagModel = new WAG01(alphabet);

//    if (alphabetType == "Proteic") {
        // The substitution model should be ReversibleSubstitutionModel.

        baseModel = dynamic_cast<ReversibleSubstitutionModel *>(PhylogeneticsApplicationTools::getSubstitutionModel(
                inputAlphabet, gCode.get(), &sites,
                modelMap, "", true, false, 0));


//        if (baseModelName == "WAG01") {
//            const ProteicAlphabet *alphabet = new ProteicAlphabet;
//            baseModel = new WAG01(alphabet);
//        } else if (modelMap["model"] == "LG08") {
//            const ProteicAlphabet *alphabet = new ProteicAlphabet;
//            baseModel = new LG08(alphabet);
//        } else if (modelMap["model"] == "JTT92") {
//            const ProteicAlphabet *alphabet = new ProteicAlphabet;
//            baseModel = new JTT92(alphabet);
//        } else if (modelMap["model"] == "JCprot") {
//            const ProteicAlphabet *alphabet = new ProteicAlphabet;
//            baseModel = new JCprot(alphabet);
//        } else if (modelMap["model"] == "DSO78") {
//            const ProteicAlphabet *alphabet = new ProteicAlphabet;
//            baseModel = new JCprot(alphabet);
//        } else {
//            ApplicationTools::displayError(
//                    "PIPIndelRateInference::init_. The substitution model should be ReversibleSubstitutionModel.");
//        }
//    } else if (alphabetType == "DNA") {
//
//    }

    TransitionModel *TrModel = new PIP13(baseModel, mu_);

    DiscreteDistribution *rDist = new ConstantRateDistribution();


    // Create joint likelihood function instance
    likelihood_ = new PIPDRHomogeneousTreeLikelihood(ttree, sites, TrModel, rDist, lambda_, mu_,
                                                     false, false);

}


/*************************** inferIndelRateFromSequences ****************************/
void PIPIndelRateInference::inferIndelRateFromSequences(PIPDRHomogeneousTreeLikelihood &likelihood,
                                                        int maxIteration,
                                                        double tolerance) {
    OptimizationFunction func(&likelihood);

    ApplicationTools::displayResult("Starting -log Likelihood value", func.getValue());

    // Brent method:
    SimpleMultiDimensions optimizer(&func);
    optimizer.setVerbose(0);
    optimizer.setProfiler(0); // showing the likelihood log 0: no profiling. in the case of debugging should be commented.
    optimizer.getOptimizationProgressCharacter();
    optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO); // to set the constraints of parameters to avoid out of range values.
    optimizer.setMaximumNumberOfEvaluations(maxIteration);// maxIteration=10000
    optimizer.getStopCondition()->setTolerance(tolerance);// tolerance=0.001
    optimizer.init(func.getParameters());
    optimizer.optimize();
    optimizer.getProfiler(); // activate the profiler again to display the results
    ApplicationTools::displayMessage("\nThe Brent Multidimensional optimization used.");

    double minf = func.getValue();

    mu_ = func.getParameterValue("mu");
    lambda_ = func.getParameterValue("lambda");

    DLOG(INFO) << "The Brent Multi-Dimensional optimization method inferred the indel parameters.";

    ApplicationTools::displayResult("The new Mu", mu_);
    ApplicationTools::displayResult("The new lambda", lambda_);
    ApplicationTools::displayResult("-logLikelihood value", minf);

}















