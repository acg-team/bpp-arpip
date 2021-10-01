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
 * @file ARPIPTools.cpp
 * @author Gholamhossein Jowkar
 * @date 16 08 2021
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


#include <Bpp/Phyl/TreeTemplate.h>
#include "ARPIPTools.h"
#include <iostream>
#include <iterator>


using namespace bpp;

/******************************************************************************/
/****************************** ARPIPTreeTools ********************************/
/******************************************************************************/
 void ARPIPTreeTools::treeAncestorRelation(Node *node, std::vector<std::string> &relation){

    int nodeId = node->getId();
    if (!node->isLeaf()) {
        int nbSons = node->getNumberOfSons();
        for (int i = 0; i < nbSons; ++i) {
            treeAncestorRelation(node->getSon(i), relation);
        }
    }
    std::string *relation_n = &relation[nodeId];

    if(node->hasFather())
        *relation_n = node->getName() + "     " + node->getFather()->getName();
    else
        *relation_n = node->getName() + "        " + "Root!";
    int nbSons = node->getNumberOfSons();
    for (int i = 0; i < nbSons; ++i) {
        *relation_n = *relation_n + "   " + node->getSon(i)->getName();
    }
}

/******************************************************************************/

void ARPIPTreeTools::renameInternalNodes(bpp::TreeTemplate<Node> *ttree, std::string prefix) {

    // Rename internal nodes with standard Vxx * where xx is a progressive number
    for (auto &nodeId:ttree->getNodesId()) {

        if (!ttree->hasNodeName(nodeId)) {

        std::string stringId;
        std::string stringName;

        stringId = std::to_string(nodeId);
        stringName = prefix + stringId;

        ttree->setNodeName(nodeId, stringName);

        }
    }

}

/******************************************************************************/



/******************************************************************************/
/******************************* ARPIPIOTools *********************************/
/******************************************************************************/

void ARPIPIOTools::writeNodeRelationToFile(std::vector<std::string> &params, std::string path) {
    std::ofstream outTreeAncestorFile;
    outTreeAncestorFile.open(path);
    outTreeAncestorFile << "Name \t Parent \t Child \n";
    std::ostream_iterator<std::string> output_iterator(outTreeAncestorFile, "\n");
    copy(params.begin(), params.end(), output_iterator);
    outTreeAncestorFile.close();
    bpp::ApplicationTools::displayResult("Output node relation name file", path);

}
/******************************************************************************/

void ARPIPIOTools::writeMLIndelPointsToFile(bpp::PIPMLIndelPoints *mlindelpoint,  const std::string &path,  bool overwrite) {
    try {
        // Open file in specified mode
        std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));

        size_t nbSite = mlindelpoint->getNbSites();
        std::vector<std::string> homoPathPtr;
        homoPathPtr.resize(nbSite);
        bpp::SitePatterns shrunkData = mlindelpoint->getShrunkData();

        // Print the Shrunk sequences:
        SiteContainer *columnSite = shrunkData.getSites();
        for (size_t i = 0; i < columnSite->getNumberOfSequences(); i++) {
            std::cout << "Seq " << i << "= " << columnSite->getName(i) << " is:" << columnSite->toString(i)
                      << " ****** " << std::endl;
        }


        std::vector<size_t> index = shrunkData.getIndices();
        for (int i = 0; i < nbSite; ++i) {
            // Extracting the position of the site in the shrunkData
            size_t rootPosition = mlindelpoint->getLikelihood()->getLikelihoodData()->getRootArrayPosition(i);
//            std::cout << rootPosition << std::endl;
            std::vector<unsigned int> arrayWeights = mlindelpoint->getLikelihood()->getLikelihoodData()->getWeights();
            std::vector<unsigned int> weight_pos = mlindelpoint->getLikelihood()->getLikelihoodData()->getWeights();
            homoPathPtr.at(i) = mlindelpoint->getHomoPath()[rootPosition];
        }

        size_t homSize = homoPathPtr.size();
        for (size_t i = 0; i < homSize; i++) {
            std::cout << i << ":" << homoPathPtr[i] << std::endl;
            output << homoPathPtr.at(i) << std::endl;
        }
        output << std::endl;
        output.close();



    }catch (IOException &e) {
        std::stringstream ss;
        ss << e.what() << "\nProblem writing mlIndels to file " << path
           << "\n Is the file path correct and do you have the proper authorizations? ";
        throw (IOException(ss.str()));
    }

}

