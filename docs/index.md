<img src="ARPIP_logo_fine_3.jpg" width="100" >
 

# ARPIP: Ancestral sequence Reconstruction under the Poisson Indel Process

ARPIP is an Ancestral Sequence Reconstruction software - a method for ancestral sequence reconstruction that uses 
Poisson Indel Process (Bouchard-Côté & Jordan, PNAS, 2013) to model insertions and deletions on a phylogenetic tree 
assuming independence among sites. ARPIP consists of two main algorithms. The first algorithm (IndelPoints) infers for
each site the insertion and deletion points on the tree topology by highest probability. The second algorithm reconstructs 
ancestral characters on the pruned subtrees in a manner similar to FastML (Pupko et al. 2000). 

---

#### Documentation

1. **[Compiling ARPIP](arpip_ancestral_sequence_reconstruction_under_poisson_indel_proccess_compile.md)**
2. **[Features and project structure](arpip_ancestral_sequence_reconstruction_under_poisson_indel_proccess_features.md)**
3. **[Executing an analysis](arpip_ancestral_sequence_reconstruction_under_poisson_indel_proccess_inference_examples.md)**
4. **[Interpret the result](arpip_ancestral_sequence_reconstruction_under_poisson_indel_proccess_result.md)**

---

#### Citation

Gh. Jowkar, et al. 2022.
**ARPIP: Ancestral sequence Reconstruction with insertions and deletions under the Poisson Indel Process.** 
*Systematic Biology, 2022;* syac050, https://doi.org/10.1093/sysbio/syac050
---

#### Licence

     * ARPIP is a computer program whose purpose is to infer Ancestral Sequence 
     * under Poisson Indel Process for nucleotide and protein datasets.
     *
     * This software is based and extends the following libraries:
     *
     * - Bio++ libraries Released version 2.4.1.
     *   developed by the Bio++ development team <http://biopp.univ-montp2.fr>
     *
     * - Google's C++ test framework version 1.11.0.
     *   developed by the google development team <https://github.com/google/googletest>
     *
     * - Google Logging Library version 0.5.0
     *   developed by the google development team <https://github.com/google/glog>
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

---

#### Reference

[1] Alexandre Bouchard-Côté and Michael I. Jordan. Evolutionary inference via the Poisson Indel Process.Proc. Natl. Acad. Sci. U.S.A., 110(4):1160, 2013.

[2] Tal Pupko and et al.  A fast algorithm for joint reconstruction of ancestral amino acid sequences Molecular biology and evolution. 17(6). 890–896, 2000.

---

#### Support
In case of bugs or improvement suggestions feel free to:
    
- Write to [Gholam-Hossein Jowkar](mailto:jowk@zhaw.ch)
    
---   

#### Funding
This work was supported by the Swiss National Science Foundation (SNF) grant [31003A_176316](https://p3.snf.ch/Project-176316) to M.A. The funding body did not play any role in the design of 
the study and collection, analysis, and interpretation of data and in writing.

<img src="sib_logo_trans_background.png" width="100" >
<img src="SNF_logo_standard_web_color_pos_e.png" width="200">
