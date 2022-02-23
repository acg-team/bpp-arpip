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
 * @file mlIndelPoints.cpp
 * @author Gholamhossein Jowkar
 * @date 30 09 2021
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

#include "gtest/gtest.h"
#include "../src/PIP/Likelihood/PIPMLIndelPoints.h"

#include <memory>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <gsl/gsl_blas_types.h>


class TestIndelPoints: public testing::Test{
protected:
    void SetUp() override {
        lambda = 10;
        mu = 0.01;
        bpp::Fasta seqReader;
        inputSequencePath = "../../test/input/msa.fasta";
        alphabet = new bpp::DNA;
        sequences = seqReader.readSequences(inputSequencePath, alphabet);
        sites = new bpp::VectorSiteContainer(*sequences);
        bpp::Site objectSite = sites->getSite(0);
        inputTreePath = "../../test/input/tree.nwk";
        bpp::Newick tReader;
        tree = tReader.read(inputTreePath);
        ttree = new bpp::TreeTemplate<bpp::Node>(*tree);
        rDist = new bpp::ConstantRateDistribution();

        sModel = new bpp::JCnuc(&bpp::AlphabetTools::DNA_ALPHABET);

        extendedSModel = new bpp::PIP13(sModel, mu);

        likFunctionPIP20 = new bpp::PIPDRHomogeneousTreeLikelihood(*ttree, *sites,
                                                                   extendedSModel, rDist,
                                                                   lambda, mu,
                                                                   false);


    }
    double lambda;
    double mu;
    bpp::Alphabet *alphabet;
    std::string inputSequencePath;
    std::string inputTreePath;
    bpp::SequenceContainer *sequences;
    bpp::Tree *tree;
    bpp::SiteContainer *sites;
    bpp::TreeTemplate<bpp::Node> *ttree;
    bpp::ReversibleSubstitutionModel *sModel;
    bpp::ReversibleSubstitutionModel *extendedSModel;
    bpp::DiscreteDistribution *rDist ;
    bpp::PIPDRHomogeneousTreeLikelihood *likFunctionPIP20;

    void TearDown() override {
        delete sequences;
        delete sites;
    }

};

TEST_F(TestIndelPoints, IndelPointsPre){
    ASSERT_TRUE(ttree != nullptr);
    ASSERT_TRUE(sequences != nullptr);
    EXPECT_EQ(likFunctionPIP20->getLikelihoodData()->getNumberOfDistinctSites(), 17);
    EXPECT_NEAR(likFunctionPIP20->getLikelihoodData()->getPipParam().getNodeBetaData(1),0.999375, 10e-3);
    EXPECT_NEAR(likFunctionPIP20->getLikelihoodData()->getPipParam().getNodeIotaData(1),0.00122399, 10e-3);
    bpp::ApplicationTools::displayTask("bpp::Likelihood module");
    bpp::ApplicationTools::displayTaskDone();

//    EXPECT_EQ(likFunctionPIP20,5);
}






