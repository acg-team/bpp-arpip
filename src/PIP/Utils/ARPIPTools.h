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
 * @file ARPIPTools.hpp
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

#ifndef ARPIP_ARPIPTOOLS_H
#define ARPIP_ARPIPTOOLS_H


#include <Bpp/Phyl/Node.h>


namespace PIPVectorTools {
    typedef std::vector<std::string> Vstring;
    typedef std::vector<Vstring> VVstring;
    typedef std::vector<VVstring> VVVstring;
    typedef std::vector<VVVstring> VVVVVstring;
}



namespace ARPIPIOTools {
    void WriteNodeRelationToFile(std::vector< std::string>& params, std::string path);

}



namespace ARPIPTreeTools {

    /**
    * @brief Father-child relation is extracted.
    * NB: This method is recursive.
    */
    void TreeAncestorRelation(bpp::Node *node, std::vector <std::string> &relation);

    /**
    * @brief Rename internal nodes with standard Vxx * where xx is a progressive number during post-order traverse of the tree.
    */
    void renameInternalNodes(bpp::TreeTemplate<bpp::Node> *ttree, std::string prefix = "V");



}

#endif //ARPIP_ARPIPTOOLS_H