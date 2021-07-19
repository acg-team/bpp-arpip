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
 * @file PIPInferenceIndelRates.cpp
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
/*
 * From the STL:
 */
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <limits>
#include <map>
#include <stdio.h>

/*
* From the BPP library:
*/
// From bpp-core
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabetState.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>//Sequence containers
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOSequenceReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Io/BppOSequenceWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

//From bpp-phyl
#include <Bpp/Phyl/Model/Protein/WAG01.h>
//#include <Bpp/Phyl/Model/FrequencySet/ProteinFrequencySet.h>
//#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include "Bpp/Phyl/SitePatterns.h"

//Phylogenetic trees
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
//#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
//#include <Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.h>
//#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Phyl/Io/Newick.h>
//#include <Bpp/Phyl/Io/BppOFrequencySetFormat.h>
//#include <Bpp/Phyl/Mapping/StochasticMapping.h>
//#include <Bpp/Phyl/Likelihood/JointLikelihoodFunction.h>
#include <Bpp/Phyl/Io/Newick.h>
//Sequence evolution models
#include <Bpp/Phyl/Model/RE08.h>
#include <Bpp/Phyl//Model/Protein/JTT92.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/NeighborJoining.h>
#include <Bpp/Phyl/OptimizationTools.h>
//Likelihood methods, including ancestral states reconstruction
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h> //Rate distributions
#include <Bpp/Phyl/Likelihood/RHomogeneousClockTreeLikelihood.h>
#include <Bpp/Phyl/NodeTemplate.h>

// From bpp-arpip
#include "PIP/Utils/ARPIPApplication.hpp"
#include "PIP/Utils/Version.hpp"

/*
* From GSL:
*/
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

/*
* From GLog:
*/
#include <glog/logging.h>
/*
* From Boost
*/
#include <boost/log/trivial.hpp>

using namespace bpp;

/******************************************************************************/
/*********************************** Main *************************************/
/******************************************************************************/


int main(int argc, char *argv[]) {

    std::cout << "Hello, World!" << std::endl;
        try {
            std::cout << "hello" << endl;

        ARPIPApplication arpipapp(argc,
                              argv,
                              std::string(software::name + " " + software::version),
                              std::string(software::releasegitbranch + " " + software::releasegitref),
                              std::string(software::releasedate + ", " + software::releasetime));

        if (argc < 2) {
//            arpipapp.help();
            exit(0);
        } else {
//            arpipapp.banner();
//            arpipapp.startTimer();
        }

        ApplicationTools::displayResult("Random seed set to", arpipapp.getSeed());
        ApplicationTools::displayResult("Log files location", std::string("current execution path"));
    /************************************************** INPUT *****************************************************/
    /* ***************************************************
     * Standard workflow
     * [INPUT]
     * 1. tree + alignment => (1.1) just parse everything
     * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using bioNJ
     * (For future release!!!)
     * 3. sequences + tree => (3.1) Perform alignment using ProPIP => (3.2) generate tree using bioNJ
     */
    double lambda = 0;
    double mu = 0;
    /********************************************* Process the MSA ************************************************/
    // ALPHABET
    // The alphabet object contains the not-extended alphabet as requested by the user,
    // while in computation we use the extended version of the same alphabet.
        std::string defaultAlphabet = ApplicationTools::getStringParameter("alphabet", arpipapp.getParams(),
                                                                         "DNA", "", true, true);

    // Alphabet without gaps
        bpp::Alphabet *alphabetNoGaps = SequenceApplicationTools::getAlphabet(arpipapp.getParams(), "", false,
                                                                                   false);

    /* Process the tree*/
    cout << "Testing parsing a directory for multiple trees..." << endl;
//    Newick tReader;
//    string treeFileName ="../test/data/input/test_01/sim-0_new.newick";
//    Tree *tree = tReader.read(treeFileName);
//    bpp::Fasta seqReader;

    /************************************ Deleting the pointers **********************************************/


        std::cout << "Hello, ARPIP!" << std::endl;
        return 0;
    } catch (exception &exception) {
        cout << "ARPIP exception:" << endl;
        std::cout << exception.what() << std::endl;
        exit(1);
    }
    return 0;
}
