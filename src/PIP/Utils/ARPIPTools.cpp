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
 void ARPIPTreeTools::TreeAncestorRelation(Node *node, std::vector<std::string> &relation){

    int nodeId = node->getId();
    if (!node->isLeaf()) {
        int nbSons = node->getNumberOfSons();
        for (int i = 0; i < nbSons; ++i) {
            TreeAncestorRelation(node->getSon(i), relation);
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

void ARPIPIOTools::WriteNodeRelationToFile(std::vector<std::string> &params, std::string path) {
    std::ofstream outTreeAncestorFile;
    outTreeAncestorFile.open(path);
    outTreeAncestorFile << "Name \t Parent \t Child \n";
    std::ostream_iterator<std::string> output_iterator(outTreeAncestorFile, "\n");
    copy(params.begin(), params.end(), output_iterator);
    outTreeAncestorFile.close();
    bpp::ApplicationTools::displayResult("Output node relation name file", path);

}


