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
 * @file ARPIPApplication.hpp
 * @author Gholamhossein Jowkar
 * @date 12 07 2021
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
#ifndef ARPIP_ARPIPAPPLICATION_HPP
#define ARPIP_ARPIPAPPLICATION_HPP

// From the STL:
#include <iostream>
#include <string>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>
#include <iostream>

// From BPP
#include <Bpp/Exceptions.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Utils/AttributesTools.h>

// From Glog
#include <glog/logging.h>

#include <boost/asio/ip/host_name.hpp>



namespace bpp{
    class ARPIPApplication {
    private:
        std::string appName_;
        std::string appBuild_;
        std::string appVersion_;
        mutable std::map<std::string, std::string> params_;
        bool timerStarted_;
        long seed_;

    public:
        ARPIPApplication(int argc, char *argv[], const std::string &name, const std::string &strVersion, const std::string &build_date):
                appName_(name), appBuild_(appBuild_), appVersion_(strVersion), params_(), timerStarted_(false) {
            params_ = bpp::AttributesTools::parseOptions(argc, argv);
            bool showversion = bpp::ApplicationTools::getBooleanParameter("version", params_, false, "", true, 3);
            bpp::ApplicationTools::warningLevel = bpp::ApplicationTools::getIntParameter("warning", params_, 0, "",
                                                                                         true, 3);
            bool noint = bpp::ApplicationTools::getBooleanParameter("noninteractive", params_, false, "", true, 3);
            bpp::ApplicationTools::interactive = !noint;
            seed_ = bpp::ApplicationTools::getParameter<long>("seed", params_, -1, "", true, 3);
            if (seed_ >= 0) {
                bpp::RandomTools::setSeed(seed_);

            } else {
                unsigned int time_ui = (unsigned int) std::time(NULL);
                seed_ = time_ui;
                bpp::RandomTools::setSeed(seed_);
            }

            // Print version header

            if (showversion) {
                this->version();
                exit(0);
            }
        }

    public:

        void done() {
            DLOG(INFO) << appName_ << "'s done.";
            if (timerStarted_)
                bpp::ApplicationTools::displayTime("Total execution time:");
        }

        std::map<std::string, std::string> &getParams() { return params_; }

        void startTimer(){
            ApplicationTools::startTimer();
            timerStarted_ = true;
        }

        const std::string &getParam(const std::string &name) const {
            if (params_.find(name) == params_.end()) throw bpp::Exception("BppApplication::getParam(). Parameter '" + name + "' not found.");
            return params_[name];
        }

        std::string &getParam(const std::string &name) { return params_[name]; }

        long getSeed() {return seed_;}

        void help() {
            std::cout << appName_ << std::endl << std::endl;
            std::cout << "Usage: Arpip [arguments] or [params=file.txt]" << std::endl;
            std::cout << "Documentation can be found at ..." << std::endl;
        }


        void banner() {

//            char *hostname[1024];
//            gethostname(*hostname, 1024);
            auto host_name = boost::asio::ip::host_name();


            bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
            bpp::ApplicationTools::displayMessage(appName_);
            bpp::ApplicationTools::displayMessage("PARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process");
            bpp::ApplicationTools::displayMessage("Authors: Gholamhossein Jowkar");
            bpp::ApplicationTools::displayMessage("Build on commit: " + appVersion_);
            bpp::ApplicationTools::displayMessage("On date: "+ appBuild_);
            bpp::ApplicationTools::displayMessage("------------------------------------------------------------------------------");
            bpp::ApplicationTools::displayResult("Execution started on:", host_name);


        }

        void version() {
            std::cout << appName_ << std::endl;
            std::cout << appVersion_ << std::endl;
            std::cout << appBuild_ << std::endl;
        }


    };


}// end of namespace bpp.

#endif //ARPIP_ARPIPAPPLICATION_HPP
