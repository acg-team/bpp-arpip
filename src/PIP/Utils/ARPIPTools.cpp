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
#include <dirent.h>


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
        *relation_n = node->getName() + "\t" + node->getFather()->getName();
    else
        *relation_n = node->getName() + "\t" + "--";
    int nbSons = node->getNumberOfSons();
    for (int i = 0; i < nbSons; ++i) {
        *relation_n = *relation_n + "\t" + node->getSon(i)->getName();
    }
}

/******************************************************************************/

bool ARPIPTreeTools::detectZeroBranchLength(bpp::Tree *tree){
    std::vector<double> branchLens = tree->getBranchLengths();
    int ans = std::count(branchLens.begin(), branchLens.end(), 0);
    if (ans > 1)
        return 1;
    return 0;
}

/******************************************************************************/

void ARPIPTreeTools::renameInternalNodes(bpp::TreeTemplate<Node> *ttree, std::string prefix) {

    // Rename internal nodes with standard Vxx * where xx is a progressive number
//    size_t i = 1; // We can have a ordinal number: 1, 2, 3 instead of progressive one.

    for (auto &nodeId:ttree->getNodesId()) {

        if (!ttree->hasNodeName(nodeId)) {

            std::string stringId;
            std::string stringName;

            stringId = std::to_string(nodeId);
            stringName = prefix + stringId;
//            stringName = prefix + std::to_string(i);

            ttree->setNodeName(nodeId, stringName);

//            i+=1;

        }
    }

}

/******************************************************************************/

size_t ARPIPTreeTools::getLongestBranchesNodeId(TreeTemplate<Node> *ttree){

    double branchLength = 0;
    size_t nodeWLongestBranch = ttree->getRootId();

    for (auto &nodeId:ttree->getNodesId()) {

        if(ttree->hasFather(nodeId))
        {
            double newBranchLength = ttree->getNode(nodeId)->getDistanceToFather();

            if (newBranchLength > branchLength) {
                nodeWLongestBranch = nodeId;
            }

        }

    }
    return nodeWLongestBranch;

}

/******************************************************************************/
void ARPIPTreeTools::scaleBranches(bpp::TreeTemplate<bpp::Node> *ttree, std::string strScale) {

    int dScale = std::stoi(strScale);
    std::vector<int> nodeIds = ttree->getNodesId();
    DLOG(INFO) << "[ARPIPTreeTools] The total tree length was " << ttree->getTotalLength() << ".";

    for (size_t i = 0; i < nodeIds.size(); i++) {
    if( nodeIds[i] != ttree->getRootId())
        ttree->setDistanceToFather(nodeIds[i], ttree->getDistanceToFather(nodeIds[i]) * dScale);
    }

    DLOG(INFO) << "[ARPIPTreeTools] Branches of the tree is now scaled." << std::endl;
    DLOG(INFO) << "[ARPIPTreeTools] The total tree length now is " << ttree->getTotalLength() << ".";

}


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

void ARPIPIOTools::writeInferredPIPParams(double lambda, double mu, long double logLikelihood, const std::string path,
                                          bool overwrite) {
    try {
        // Open file and writing the content
        std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
        output << "Lambda:" << lambda << std::endl;
        output << "Mu:" << mu << std::endl;
        output << "Loglikelihood:" << logLikelihood << std::endl;
        output << std::endl;
        output.close();
        bpp::ApplicationTools::displayResult("Output estimated PIPParameter file", path);


    } catch (IOException &e) {
        std::stringstream ss;
        ss << e.what() << "\nProblem writing PIPParameter to file " << path
           << "\n Is the file path correct and do you have the proper authorizations? ";
        throw (IOException(ss.str()));
    }
}

/******************************************************************************/

void
ARPIPIOTools::writeMLIndelPointsToFile(bpp::PIPMLIndelPoints *mlindelpoint, const std::string &path, bool overwrite,
                                       bool verbose) {

    try {
        std::cout << "Display the progress...................: " << std::boolalpha << verbose << "." << std::endl;
        size_t nbSites = mlindelpoint->getNbSites();
        std::vector<std::string> homoPathPtr;
        homoPathPtr.resize(nbSites);
        bpp::SitePatterns shrunkData = mlindelpoint->getShrunkData();

        // Print the Shrunk sequences:
        if (verbose) {
            SiteContainer *columnSite = shrunkData.getSites();
            bpp::ApplicationTools::displayMessage("The input MSA:");
            for (size_t i = 0; i < columnSite->getNumberOfSequences(); i++) {
                std::cout << "Seq " << i << "= " << columnSite->getName(i) << " is:" << columnSite->toString(i)
                          << " ****** " << std::endl;
            }
            std::cout << std::endl;
        }

        // Mapping the shrunk data to the normal one:
        std::vector<size_t> index = shrunkData.getIndices();
        for (int i = 0; i < nbSites; ++i) {
            // Extracting the position of the site in the shrunkData
            size_t rootPosition = mlindelpoint->getLikelihood()->getLikelihoodData()->getRootArrayPosition(i);
//            std::cout << rootPosition << std::endl;
            std::vector<unsigned int> arrayWeights = mlindelpoint->getLikelihood()->getLikelihoodData()->getWeights();
            std::vector<unsigned int> weight_pos = mlindelpoint->getLikelihood()->getLikelihoodData()->getWeights();
//            homoPathPtr.at(i) = mlindelpoint->getHomoPath()[rootPosition];

            // Processing the indel string:
            std::string str_homo = mlindelpoint->getHomoPath()[rootPosition];
            // New storage for new homo_path
            std::string new_str_homo{};
            std::string delimiter_semicolon = ";";
            size_t pos_h = 0;
            // removing the ";" at the end of each string and finding position of ";" using find()
            while ((pos_h = str_homo.find(delimiter_semicolon)) != std::string::npos) {
                // copy the substring before ";"
                std::string sub = str_homo.substr(0, pos_h);
                str_homo = str_homo.substr(pos_h + delimiter_semicolon.length());
                if (!sub.empty()) {
                    // resetting the pos value for the next usage
                    size_t pos = 0;
                    std::string delimiter_colon = ":";
                    std::string token;
                    //find the position of ":"
                    while ((pos = sub.find(delimiter_colon)) != std::string::npos) {
                        token = sub.substr(0, pos);
                        if (!token.empty()) {
                            new_str_homo = new_str_homo + mlindelpoint->getTree()->getNodeName(std::stoi(token)) + ":" +
                                           sub.substr(pos + 1, sub.length()) + ";";
                        }
                        sub.erase(0, pos + delimiter_semicolon.length());
                    }
                }
            }
            homoPathPtr.at(i) = new_str_homo;
        }

        // Open file and writing the content
        std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));

        size_t homSize = homoPathPtr.size();
        for (size_t i = 0; i < homSize; i++) {
            if(verbose)
                std::cout << "Site " << i << ":" << homoPathPtr[i] << std::endl;
            output << homoPathPtr.at(i) << std::endl;
        }
        output << std::endl;
        output.close();
        bpp::ApplicationTools::displayResult("Output IndelPoint file", path);


    } catch (IOException &e) {
        std::stringstream ss;
        ss << e.what() << "\nProblem writing Indels to file " << path
           << "\n Is the file path correct and do you have the proper authorizations? ";
        throw (IOException(ss.str()));
    }

}

/******************************************************************************/

std::vector<char*> ARPIPIOTools::ReadDirectory(char const *path){

    DIR *dir;
//    dir = path;
    struct dirent *diread;
    std::vector<char *> files;
    if ((dir = opendir(path)) != nullptr) {
        while ((diread = readdir(dir)) != nullptr) {
            files.push_back(diread->d_name);
        }
        closedir (dir);
    } else {
        perror ("opendir");
        std::cout << "FAILED TO READ!!!" << std::endl;
    }
    for (auto file : files)
        std::cout << file << "| ";
    std::cout << std::endl;

    return files;
}