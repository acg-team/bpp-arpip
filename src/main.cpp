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

/*
* From the BPP library:
*/
// From bpp-core

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabetState.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>//Sequence containers
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/Fasta.h>



//From bpp-phyl
//Sequence evolution models
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/RE08.h>
#include <Bpp/Phyl/Model/Protein/WAG01.h>
#include <Bpp/Phyl/Model/FrequenciesSet/ProteinFrequenciesSet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
//Phylogenetic trees
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplate.h>



// From bpp-arpip
#include "PIP/Utils/ARPIPApplication.h"
#include "PIP/Utils/ARPIPVersion.h"
#include "PIP/Utils/ARPIPTools.h"
#include "PIP/Likelihood/PIPDRHomogeneousTreeLikelihood.h"
#include "PIP/Utils/PIPIndelRateInference.h"
#include "PIP/Likelihood/PIPMLIndelPoints.h"
#include "PIP/AncestralStateReconstruction/PIPAncestralStateReconstruction.h"

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
//#include <boost/log/trivial.hpp>

/*
* From Gtest
*/
#include <gtest/gtest.h>


/**********************************************************************************************************************/
/****************************************************** Main **********************************************************/
/**********************************************************************************************************************/

int main(int argc, char *argv[]) {
    FLAGS_log_dir = "../logs/";
    ::google::InitGoogleLogging(software::name.c_str());
    ::google::InstallFailureSignalHandler();


    try {
        bpp::ARPIPApplication arpipapp
                (argc,
                 argv,
                 std::string(software::name + " " + software::version),
                 std::string(software::releasegitbranch + " " + software::releasegitref),
                 std::string(software::releasedate + ", " + software::releasetime));

        if (argc < 2) {
            arpipapp.help();
            exit(0);
        } else {
            arpipapp.banner();
            arpipapp.startTimer();
        }

//        bool PAR = bpp::ApplicationTools::getBooleanParameter("", arpipapp.getParams(), false);
        /*
         * Read all the files in directory:
         */
//        std::string inPath = "../data/input/";
//        std::vector<char *> files = ReadDirectory(in_path.c_str());

        bpp::ApplicationTools::displayResult("Random seed set to", arpipapp.getSeed());
        bpp::ApplicationTools::displayResult("Log files location", std::string("current execution path"));
        /************************************************** INPUT *****************************************************/
        /* ***************************************************
         * Standard workflow
         * [INPUT]
         * 1. tree + alignment => (1.1) just parse everything
         * 2. alignment  => (2.1) parse alignment => (2.2) generate tree using one of supported methods
         * 3. sequences + tree => (3.1) Perform alignment using ProPIP => (3.2) generate tree using one of supported methods (For future release!!!)
         */

        double lambda = 0.01;
        double mu = 0.01;
        bool estimatedPIPParameters = false;

        std::string App_model_substitution = bpp::ApplicationTools::getStringParameter("model", arpipapp.getParams(),
                                                                                       "JC69", "",
                                                                                       true, true);

        std::string modelStringName;
        std::map<std::string, std::string> modelMap;
        bpp::KeyvalTools::parseProcedure(App_model_substitution, modelStringName, modelMap);
        bool App_model_indels = modelStringName == "PIP";

        // AMONG-SITE-RATE-VARIATION
        // Markov-modulated Markov model: constant!
        bpp::DiscreteDistribution *rateDist = new bpp::ConstantRateDistribution();
        bpp::ApplicationTools::displayMessage("This model works with constant rate distribution for DNA/AA.");

        /********************************************* Process the MSA ************************************************/
        // ALPHABET
        // The alphabet object contains the not-extended alphabet as requested by the user,
        // while in computation we use the extended version of the same alphabet.
        std::string App_alphabet = bpp::ApplicationTools::getStringParameter("alphabet", arpipapp.getParams(),
                                                                         "DNA", "", true, true);

        // Alphabet without gaps
        bpp::Alphabet *alphabetNoGaps = bpp::SequenceApplicationTools::getAlphabet(arpipapp.getParams(), "", false,
                                                                                       false);

        // Genetic code
        unique_ptr<bpp::GeneticCode> gCode;

        // Alphabet used for all the computational purpose (it can allows for gap extension)
        bpp::Alphabet *alphabet;
        if (App_alphabet.find("DNA") != std::string::npos ) {
            alphabet = new bpp::DNA();
        } else if (App_alphabet.find("Protein") != std::string::npos) {
            alphabet = new bpp::ProteicAlphabet();
        }

        // Check what is the alphabet:
        bpp::ApplicationTools::displayResult("The alphabet is", alphabet->getAlphabetType());

        bpp::ApplicationTools::displayMessage("\n[Preparing input Files]");
        std::string input_sequences = bpp::ApplicationTools::getAFilePath("input.sequence.file", arpipapp.getParams(),
                                                                         true, true, "",
                                                                         false, "", 1);

        bpp::SequenceContainer *sequences = nullptr;
        bpp::SiteContainer *sites = nullptr;

        try {

            // Read aligned sequences
            bpp::Fasta seqReader;
            sequences = seqReader.readSequences(input_sequences, alphabet);


            bpp::VectorSiteContainer *allSites = bpp::SequenceApplicationTools::getSiteContainer(alphabet, arpipapp.getParams());

            sites = bpp::SequenceApplicationTools::getSitesToAnalyse(*allSites, arpipapp.getParams(), "", true,
                                                                0, true, 1);

            // Removing gap within sites:
            //sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, arpipapp.getParams(), "", true,
            //                                                    1, true, 1);

            delete allSites;

            bpp::ApplicationTools::displayResult("Number of sequences", sites->getNumberOfSequences());
            bpp::ApplicationTools::displayResult("Number of sites", sites->getNumberOfSites());


        } catch (bpp::Exception &e) {
            LOG(FATAL) << "Error when reading sequence file due to: " << e.message();
        }



        /**************************************** Process the tree ****************************************************/
        // TREE
        bpp::ApplicationTools::displayMessage("\n[Preparing initial tree]");
        bpp::Tree *tree = nullptr;
        string initTreeOpt = bpp::ApplicationTools::getStringParameter("init.tree", arpipapp.getParams(),
                                                                       "user", "", false, 1);
        bpp::ApplicationTools::displayResult("Initial tree", initTreeOpt);
//        DLOG(INFO) << "Tree Description:\n";
        if (initTreeOpt == "user") {

            DLOG(INFO) << "[Input tree parser] Tree Provided by User.\n";
            tree = bpp::PhylogeneticsApplicationTools::getTree(arpipapp.getParams());
            DLOG(INFO) << "[Input tree parser] Tree name is:" << tree->getName() << std::endl;
        } else {
            // If there is no tree as input
            std::string App_distance_method = bpp::ApplicationTools::getStringParameter("init.tree.method",
                                                                                        arpipapp.getParams(), "nj");
            bpp::ApplicationTools::displayResult("Initial tree reconstruction method", App_distance_method);
            DLOG(INFO) << "[Tree reconstruction] Reconstruction Method: " << App_distance_method;

            bpp::AgglomerativeDistanceMethod *distMethod = nullptr;

            std::string token = App_distance_method.substr(0, App_distance_method.find("-"));

            if (token == "wpgma") {
                auto *wpgma = new bpp::PGMA(true);
                distMethod = wpgma;
            } else if (token == "upgma") {
                auto *upgma = new bpp::PGMA(false);
                distMethod = upgma;
            } else if (token == "nj") {
                auto *nj = new bpp::NeighborJoining();
                nj->outputPositiveLengths(true);
                distMethod = nj;
            } else if (token == "bionj") {
                auto *bionj = new bpp::BioNJ();
                bionj->outputPositiveLengths(true);
                distMethod = bionj;
            } else throw bpp::Exception("Tree reconstruction method is not supported.");

            // Building tree using one of the mentioned methods
            /* ***************************************************
             * Standard workflow
             * 1. (1.1) Build temporary substitution model (1.2) build a distance matrix using the model (1.3) tree reconstruction.
            */

            // PIP parameter is computed using known tree and for computing the tree we have to assume that the PIP param is given

            std::string tmpBaseModel;

            std::map<std::string, std::string> tmpBaseModelMap;
            bpp::KeyvalTools::parseProcedure(modelMap["model"], tmpBaseModel, tmpBaseModelMap);

            bpp::SubstitutionModel *tmpSModel = bpp::PhylogeneticsApplicationTools::getSubstitutionModel(alphabetNoGaps,
                                                                                                         gCode.get(),
                                                                                                         sites,
                                                                                                         modelMap,
                                                                                                         "",
                                                                                                         true,
                                                                                                         false,
                                                                                                         0);
            bpp::ApplicationTools::displayResult(
                    "Default model without gap is used", tmpSModel->getName());
            DLOG(INFO)
                    << "[Tree reconstruction] Normal substitution model is used, meaning no gap were considered for tree reconstruction. ";
            // Cloning the site for building the tree
            bpp::SiteContainer *tmpSites = sites->clone();

            // Removing gappy regions
            bpp::ApplicationTools::displayTask("Removing gap only sites");
            bpp::SiteContainerTools::removeGapOnlySites(*tmpSites);
            DLOG(INFO) << "[Tree reconstruction] Sites with full gap column were removed.";

            bpp::ApplicationTools::displayTask("Changing gaps to unknown characters");
            bpp::SiteContainerTools::changeGapsToUnknownCharacters(*tmpSites);
            DLOG(INFO) << "[Tree reconstruction] Gap states were changed to unknown characters.";


            // Computing distance matrix
            bpp::DistanceEstimation distanceMethod(tmpSModel, rateDist, tmpSites);
            DLOG(INFO) << "[Tree reconstruction] Distance matrix was computed.";


            // Retrieve the computed distances
            bpp::DistanceMatrix *distances = distanceMethod.getMatrix();

            distMethod->setDistanceMatrix(*distances);
            distMethod->computeTree();
            tree = distMethod->getTree();
            DLOG(INFO) << "[Tree reconstruction] Reconstructed by: " << distMethod->getName() << std::endl;

            delete tmpSites;
            delete distances;
            delete distMethod;
            std::cout << endl;
        }


        bpp::TreeTemplate<bpp::Node> *ttree_ = new bpp::TreeTemplate<bpp::Node>(*tree);

        // In the case of scaling the tree:
        string initTreeScale = bpp::ApplicationTools::getStringParameter("opt.tree.scale", arpipapp.getParams(),
                                                                       "1", "", false, 1);
        if (initTreeScale != "1") {
            bpp::ApplicationTools::displayResult("The tree scale is activated with value ", initTreeScale);
            DLOG(INFO) << "[Input tree parser] Tree is now scaling " << initTreeScale << " times.";
            ARPIPTreeTools::scaleBranches(ttree_, initTreeScale);
        }

        // If tree is mulifurcation, then resolve it with midpoint rooting
        if (ttree_->isMultifurcating()) {
            try{
                bpp::ApplicationTools::displayMessage("Tree is multifuracted.");
                bpp::TreeTemplateTools::midRoot(*(ttree_), bpp::TreeTemplateTools::MIDROOT_VARIANCE, false);
                DLOG(INFO) << "[Input tree parser] Tree is now binary by Midpoint rooting method.";

            } catch (bpp::Exception e){
                bpp::ApplicationTools::displayError("Error when multifuracting the tree");
                LOG(FATAL) << "[Input tree parser] Error when multifuracting the tree" << e.message();
            }
        }

        // If tree is not rooted, resolve it by picking a random node
        if (ttree_->isRooted()) {
            int root = tree->getRootId();
            bpp::ApplicationTools::displayResult("The tree root is node number", root);
            DLOG(INFO) << "[Input tree parser] Tree is rooted and now the root is " << root << endl;
        } else {
            bpp::ApplicationTools::displayMessage("Tree is not rooted: the tree must have a root in PIP model!!!!");
            DLOG(INFO) << "[Input tree parser] The input tree is not rooted, the tree must have a root in PIP model!!!!" << endl;
//            int root = ttree_->getRootId();
            std::vector<double> nodeBr = ttree_->getBranchLengths();
            size_t newRoot = ARPIPTreeTools::getLongestBranchesNodeId(ttree_);
            bpp::Node *n = ttree_->getNode(newRoot);
            ttree_->newOutGroup(newRoot);// it should be a random node!
            if(ttree_->isRooted()){
                bpp::ApplicationTools::displayResult("New random root is", newRoot);
                DLOG(INFO) << "[Input tree parser] Now, the tree is rooted....." << "(node " << newRoot << " is the new root)" << endl;
            } else throw bpp::Exception("Tree is not rooted yet.");
        }


        // Tree description
        bpp::ApplicationTools::displayResult("Number of nodes", ttree_->getNumberOfNodes());
        bpp::ApplicationTools::displayResult("Number of leaves", ttree_->getNumberOfLeaves());
        bpp::ApplicationTools::displayResult("Total tree length, aka tau", ttree_->getTotalLength());
        DLOG(INFO) << "[Input tree parser] Number of Nodes: " << tree->getNumberOfNodes() << std::endl;
        DLOG(INFO) << "[Input tree parser] Number of leaves: " << tree->getNumberOfLeaves() << std::endl;
        DLOG(INFO) << "[Input tree parser] Total tree length: " << tree->getTotalLength() << std::endl;


        DLOG(INFO) << "[Input tree parser] Initial tree topology is: " << bpp::TreeTools::treeToParenthesis(*tree, true);
        // Rename internal nodes with standard Vxx * where xx is a progressive number
        ttree_->setNodeName(tree->getRootId(), "root");
        ARPIPTreeTools::renameInternalNodes(ttree_);

        // Write down the reconstructed tree
        bpp::PhylogeneticsApplicationTools::writeTree(*ttree_, arpipapp.getParams());
        DLOG(INFO) << "[Input tree parser] Initial tree topology after topological modifications: "
                   << bpp::TreeTools::treeToParenthesis(*ttree_, true);

        // Write down the relation between nodes
        vector<string> strTreeFileAncestor;
        strTreeFileAncestor.resize(ttree_->getNumberOfNodes());
        // Store relation between node in a string vector
        ARPIPTreeTools::treeAncestorRelation(ttree_->getRootNode(), strTreeFileAncestor);
        ARPIPIOTools::writeNodeRelationToFile(strTreeFileAncestor, arpipapp.getParam("output.node_rel.file"));



        /*********************************** Process the substitution model *******************************************/

        // SUBSTITUTION MODEL
        bpp::ApplicationTools::displayMessage("\n[Setting up substitution model]");

        bpp::ReversibleSubstitutionModel *sModel = nullptr;
//        bpp::TransitionModel *tModel = nullptr;


        // Instantiate a substitution model and extend it with PIP
        bool computeFreqFromData = false;
        std::string baseModel{};

        std::map<std::string, std::string> baseModelMap;
        bpp::KeyvalTools::parseProcedure(modelMap["model"], baseModel, baseModelMap);

        std::vector<std::string> modelKeys{};
        for (auto k = baseModelMap.begin(); k != baseModelMap.end(); ++k)modelKeys.push_back(k->first);

        if (!modelKeys.empty()) {
            baseModel += "(";
            for (auto &key:modelKeys) {
                if (key != "initFreqs") {
                    baseModel += key + "=" + baseModelMap[key];
                } else {
                    if (baseModelMap[key] == "observed") {
                        computeFreqFromData = true;
                    }
                }
                baseModel += ",";
            }
            baseModel.pop_back();
            baseModel += ")";
            modelMap["model"] = baseModel;
        }

        if (App_alphabet.find("Protein") != std::string::npos || App_alphabet.find("DNA") != std::string::npos) {
            sModel = dynamic_cast<bpp::ReversibleSubstitutionModel *> (bpp::PhylogeneticsApplicationTools::getSubstitutionModel(
                    alphabetNoGaps, gCode.get(), sites, modelMap, "", true, false, 0));
        }else{
            sModel = dynamic_cast<bpp::ReversibleSubstitutionModel *>(bpp::PhylogeneticsApplicationTools::getSubstitutionModel(
                    alphabet, gCode.get(), sites, modelMap,
                    "", true, false, 0));
        }
        DLOG(INFO) << "[Substitution model] Base substitution model is " << sModel->getName() << " with "
                   << (int) sModel->getNumberOfStates() << " number of states.";
        bpp::ApplicationTools::displayResult("Base substitution model", sModel->getName());


        // If PIP, then check if lambda/mu initial values are estimated from the data
        if (modelMap.find("estimated") != modelMap.end()) {
            bpp::ApplicationTools::displayError(
                    "The use of the tag [observed] is obsolete. Use the tag [initFreqs] instead.");
            exit(1);
        } else if (modelMap.find("initFreqs") != modelMap.end()) {
            if (modelMap["initFreqs"] == "observed") {
                estimatedPIPParameters = true;
            }
        } else if (modelMap.find("lambda") == modelMap.end() ||
                   modelMap.find("mu") == modelMap.end()) {
            // One of the PIP parameters is missing
            estimatedPIPParameters = true;
        }

        if (computeFreqFromData) {
            sModel->setFreqFromData(*sites, 0.01);
            DLOG(INFO) << "[Substitution model] Base substitution model used frequency from data" ;
            bpp::ApplicationTools::displayMessage("Set the frequency from data");
        }


        /////////////////////////////
        // Estimate PIP Parameters: inference of Indel rates
        if (estimatedPIPParameters) {
            bpp::PIPIndelRateInference *PIPIndelParam = new bpp::PIPIndelRateInference(*sites, *ttree_, modelMap, &lambda,
                                                                                       &mu);
            lambda = PIPIndelParam->getLambda();
            mu = PIPIndelParam->getMu();

            DLOG(INFO) << "[PIP model] Estimated PIP parameters from data using input sequences (lambda=" <<
                       lambda << ",mu=" << mu << "," "I=" << lambda * mu << ")";
            delete PIPIndelParam;
        } else {
            lambda = (modelMap.find("lambda") == modelMap.end()) ? 0.1 : std::stod(modelMap["lambda"]);
            mu = (modelMap.find("mu") == modelMap.end()) ? 0.2 : std::stod(modelMap["mu"]);

            //sim-0 True:
            // lambda = 10 ;
            // mu = 0.1

            // prank
            //mu = 0.210789;
            //lambda = 163.872;

            // ProPIP
            //mu = 0.242833;
            //lambda = 210.235;
            DLOG(INFO) << "[Substitution model][PIP model] Fixed PIP parameters to (lambda=" << lambda << ",mu=" << mu << "," "I="
                       << lambda * mu << ")";
        }




        // New PIP substitution model:
        sModel = new bpp::PIP13(sModel, mu);

        bpp::ApplicationTools::displayResult("Substitution model", sModel->getName());

        if (App_model_indels)
            bpp::ApplicationTools::displayResult("Indel parameter initial value",
                                            (estimatedPIPParameters) ? "estimated" : "fixed");

        bpp::ParameterList parameters = sModel->getParameters();
        for (size_t i = 0; i < parameters.size(); i++) {
            bpp::ApplicationTools::displayResult(parameters[i].getName(), bpp::TextTools::toString(parameters[i].getValue()));
        }

        double pip_intensity = lambda * parameters.getParameter("PIP13.mu").getValue();
        bpp::ApplicationTools::displayResult("PIP13.intensity", bpp::TextTools::toString(pip_intensity));

        DLOG(INFO) << "[Substitution model] Substitution model is " << sModel->getName() << " model with mu="
                   << parameters.getParameter("PIP13.mu").getValue();

        bpp::ApplicationTools::displayMessage("Equilibrium frequency of substitution model:");
        for (size_t i = 0; i < sModel->getFrequencies().size(); i++) {

            bpp::ApplicationTools::displayResult("eq.freq(" + sModel->getAlphabet()->getName(i) + ")",
                                            bpp::TextTools::toString(sModel->getFrequencies()[i], 4));
        }

        bpp::StdStr s1;
        bpp::PhylogeneticsApplicationTools::printParameters(sModel, s1, 1, true);
        DLOG(INFO) <<"[Substitution model] Value from the model is " << s1.str();

        /////////////////////////
        // Among site rate variation (ASVR)
        bpp::DiscreteDistribution *rDist = new bpp::ConstantRateDistribution();
        bpp::ApplicationTools::displayResult("Among site rate variation (ASVR)", rDist->getName());
        DLOG(INFO) << "[Substitution model] Constant distribution rate is used in this model.";

        /************************************* Tree Likelihood Computation ********************************************/

        bpp::ApplicationTools::displayMessage("\n[Setting up Tree Likelihood]");

        bpp::PIPDRHomogeneousTreeLikelihood *likFunctionPIP20 = new bpp::PIPDRHomogeneousTreeLikelihood(*ttree_, *sites,
                                                                                                        sModel, rDist,
                                                                                                        lambda, mu,
                                                                                                        false);
        DLOG(INFO) << "[PIP tree likelihood] Likelihood object under PIP is constructed successfully.";

        //////////////////////////////
        // Estimate the PIPLogLikelihood value:
        string optPIPLikelihoodCmp = bpp::ApplicationTools::getStringParameter("opt.likelihood", arpipapp.getParams(),
                                                                               "0", "", false, 1);
        if (stoi(optPIPLikelihoodCmp) == 1) {
            long double logLikelihood = likFunctionPIP20->computePIPTreeLikelihood(lambda, mu);
            bpp::ApplicationTools::displayResult("The insertion rate", lambda);
            bpp::ApplicationTools::displayResult("The deletion rate", mu);
            bpp::ApplicationTools::displayResult("The log likelihood under PIP", logLikelihood);
            bpp::ApplicationTools::displayWarning(
                    "By activating \"opt.likelihood\" in config file the operation would be terminated!");
            bpp::ApplicationTools::displayMessage(
                    "By running the ASR module, the PIP parameters including logLikelihood would be printed.");
            DLOG(INFO) << "[PIP tree likelihood] Just Log likelihood asked to be printed.";
            return 1;
        }
        /////////////////////////
        bpp::ApplicationTools::displayMessage("\n[Extracting Indel Points]");

        // Maximum Likelihood Indel Points
        bpp::PIPMLIndelPoints *mlIndePoints = new bpp::PIPMLIndelPoints(likFunctionPIP20);
        DLOG(INFO) << "[Indel Points] Indel points using Maximum Likelihood are inferred successfully.";

        // Write down the output of the algorithm to a single file
//        ARPIPTreeTools::treeAncestorRelation(ttree_->getRootNode(), strTreeFileAncestor);
        ARPIPIOTools::writeMLIndelPointsToFile(mlIndePoints, arpipapp.getParam("output.mlindelpoints.file"));
        DLOG(INFO) << "[Indel Points] File is written successfully.";


        /************************************* Ancestral Sequence Reconstruction **************************************/

        bpp::ApplicationTools::displayMessage("\n[Computing Ancestral Sequences Reconstruction]");


        bpp::PIPAncestralStateReconstruction jarPIP(likFunctionPIP20, mlIndePoints);

        bpp::ApplicationTools::displayTask("Ancestral sequence reconstruction");
        bpp::AlignedSequenceContainer *asr = jarPIP.JointAncestralSequencesReconstruction();
        bpp::ApplicationTools::displayTaskDone();
        DLOG(INFO) << "[PIP ASR] Ancestral sequence were successfully reconstructed.";

        DVLOG(2) << "[PIP ASR] This sequence is coded with a" << asr->getAlphabet()->getAlphabetType();
        for (size_t nbseq = 0; nbseq < asr->getNumberOfSequences(); nbseq++) {
            // To change the Alphabet of a sequence, we need to decode and recode it:
            try {
                bpp::Sequence *sequence = new bpp::BasicSequence(asr->getSequence(nbseq).getName(), asr->getSequence(nbseq).toString(), alphabet);
                DVLOG(2) << "[PIP ASR] This sequence is now coded with a" << sequence->getAlphabet()->getAlphabetType();
                DVLOG(3) << "[PIP ASR] The sequence is " << sequence->toString();
                delete sequence;
            } catch (bpp::Exception &ex) {
                cerr << ex.what() << endl;
                LOG(FATAL) << "Error when " << ex.message();
            }
        }

        /************************************ Writing the result **********************************************/

        bpp::Fasta fastaWtiter;
        fastaWtiter.writeSequences(arpipapp.getParam("output.ancestral.file"), *asr);
        LOG(INFO) << "[PIP ASR] File is written successfully.";
        bpp::ApplicationTools::displayResult("Output ASR file", arpipapp.getParam("output.ancestral.file"));

        /**************************************** Deleting the pointers ***********************************************/

        delete mlIndePoints;
        delete sequences;
        delete asr;
        delete rDist;
        delete rateDist;

        arpipapp.done();
        google::ShutdownGoogleLogging();
        exit(0);

    } catch (exception &ex) {
        cout << "ARPIP exception:" << endl;
        std::cout << ex.what() << std::endl;
        LOG(FATAL) << "Error when " << ex.what();
        google::ShutdownGoogleLogging();
        exit(1);
    }
    return 0;
}
