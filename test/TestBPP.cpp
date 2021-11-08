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
 * @file test_bpp.cpp
 * @author Gholamhossein Jowkar
 * @date 09 08 2021
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

// From STL:
#include <iostream>
#include <filesystem>
#include <memory>

// From Gtest:
#include <gtest/gtest.h>

// From BPP
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>



class TestSequenceReader: public testing::Test {
protected:
    void SetUp()  override {
        std::string currentPath = get_current_dir_name();
        std::string inputSequencePath = "../data/input/test_00/";

    }

    bpp::Alphabet *alphabet;
    std::string inputSequencePath;
    std::string outputSequncePath;

    void TearDown() override {
        delete alphabet;
    }

};

TEST_F(TestSequenceReader, ReadmsaFile){
    bpp::Fasta seqReader;
    bpp::SequenceContainer *sequences = seqReader.readSequences(inputSequencePath, alphabet);
    bpp::SiteContainer *sites = new bpp::VectorSiteContainer(*sequences);
    ASSERT_TRUE(sequences != nullptr);
    EXPECT_EQ(sites->getNumberOfSequences(),8);
    EXPECT_EQ(sites->getNumberOfSites(), 14);
    EXPECT_EQ(sites->getSequencesNames()[0],"leaf_0");
    EXPECT_EQ(sites->getSequence(0).toString(),"ATGTC----ACTGT");

    delete sequences;
    delete sites;
}

class TestTreeReader: public testing::Test{
protected:
    void SetUp() override {

    }
//    bpp::Alphabet *alphabet;
    std::string inputTreePath;
    std::string outputTreePath;
    void TearDown() override {}
};

TEST_F(TestTreeReader, ReadTreeFile){
    bpp::Newick *newickReader;
    bpp::TreeTemplate<bpp::Node> *ttree = newickReader->read(inputTreePath);
    ASSERT_TRUE(ttree != nullptr);
    EXPECT_EQ(ttree->getNumberOfNodes(),15);
    EXPECT_EQ(ttree->getNumberOfLeaves(),8);
    EXPECT_EQ(ttree->getTotalLength(),7.21381);
}