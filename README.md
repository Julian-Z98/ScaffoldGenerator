# Scaffold Generator
##### A CDK-based library for generating Scaffold Trees and Scaffold Networks

## Contents of this document
* [Description](#Description)
* [Functionalities and options](#Functionalities-and-options)
  * [Available functionalities](#Available-functionalities)
  * [Deviations from the Scaffold Tree prioritization rules](#Deviations-from-the-Scaffold-Tree-prioritization-rules)
  * [Aromaticity handling](#Aromaticity-handling)
  * [Available settings and options](#Available-settings-and-options)
  * [Notes about the implementation](#Notes-about-the-implementation)
* [Contents of this repository](#Contents-of-this-repository)
  * [Sources](#Sources)
  * [Resources](#Resources)
* [Installation](#Installation)
* [Dependencies](#Dependencies)
* [References and useful links](#References-and-useful-links)

## Description
The Scaffold Generator library is designed to make molecular scaffold-related functionalities available in applications 
and workflows based on the [Chemistry Development Kit (CDK)](https://cdk.github.io/). Building upon the works by 
[Bemis and Murcko](https://doi.org/10.1021/jm9602928), [Schuffenhauer et al.](https://doi.org/10.1021/ci600338x), 
and [Varin et al.](https://doi.org/10.1021/ci2000924), it offers scaffold 
perception and dissection based on single molecules and molecule collections. 
From the latter, Scaffold Trees and Scaffold Networks can also be constructed, represented in data structures and visualised 
using the [GraphStream library](https://graphstream-project.org). Multiple options to fine-tune and adapt the routines are available.

## Functionalities and options
### Available functionalities
From a molecule represented by a CDK IAtomContainer object, the molecular scaffold can be extracted. According to
[Bemis and Murcko](https://doi.org/10.1021/jm9602928), this is constituted by its rings and the non-cyclic structures 
connecting them (linkers). Terminal side-chains are excluded. Different scaffold types based on this first definition
can be selected in Scaffold Generator (see below). Additionally, the separate building blocks of the scaffold, rings and 
linkers, and the side-chains removed for scaffold generation can be extracted. Ring perception is based on the CDK "relevant"
cycle finder algorithm that extracts the smallest set of uniquely defined short cycles. Fused ring systems will therefore 
be dissected into their constituting separate smallest rings.<p>
Extracted scaffolds can be further dissected into their smaller parent scaffolds using two different methods. All possible 
parent scaffolds can be enumerated that would result 
from the step-wise removal of terminal rings, exploring all possible combinations of such removal steps. This dissection
is the basis for Scaffold Networks, as described by [Varin et al.](https://doi.org/10.1021/ci2000924).
<br>[Schuffenhauer et al.](https://doi.org/10.1021/ci600338x) built upon a similar scaffold dissection but introduced 
13 chemical rules to prioritize one specific parent scaffold at every step. Applying these rules, only one specifically 
determined terminal ring is chosen to be removed at every stage to generate one specifically chosen parent scaffold. This
procedure is the basis for [Schuffenhauer et al.'s](https://doi.org/10.1021/ci600338x) Scaffold Trees.
<br>Both scaffold dissection methods can be applied to a molecule using Scaffold Generator. But it is also possible 
to directly apply them to a collection of molecules and thus generate Scaffold Networks based on [Varin et al.](https://doi.org/10.1021/ci2000924)
and Scaffold Trees based on [Schuffenhauer et al.](https://doi.org/10.1021/ci600338x). For both approaches, data structures 
are available in this library to manage the collection of resulting parent/child scaffolds and their connections in a 
graph-based way that can be visualised using the [GraphStream library](https://graphstream-project.org).
<br>Other functionalities of the data structures include the retrieval of chemical information from the scaffolds and 
their origin molecules and the determination of virtual scaffolds, i.e. scaffolds that were only produced by dissection of 
bigger child scaffolds and not present in the original data set.

### Deviations from the Scaffold Tree prioritization rules
non-single bonds
different scaffolds/frameworks
ring perception

### Aromaticity handling


### Available settings and options

### Notes about the implementation


## Contents of this repository
### Sources


### Resources


## Installation
This is a Maven project. In order to use the source code for your own software, download or clone the repository and
open it in a Maven-supporting IDE (e.g. IntelliJ) as a Maven project and execute the pom.xml file. Maven will then take
care of installing all dependencies. A Java Development Kit (JDK) of version 11 or higher must also be pre-installed.
<br>To run the COCONUT-analysing tests, an SD file of the database needs to be placed in the test resources folder. 
The respective file can be downloaded at
[https://coconut.naturalproducts.net/download](https://coconut.naturalproducts.net/download). Add the "COCONUT_DB.sdf" 
file to <i>src\test\resources</i>. 

## Dependencies
**Needs to be pre-installed:**
* Java Development Kit (JDK) version 11
  * [AdoptOpenJDK](https://adoptopenjdk.net) (as one possible source of the JDK)
* Apache Maven version 4
  * [Apache Maven](http://maven.apache.org)

**Managed by Maven**
* Chemistry Development Kit (CDK) version 2.5
  * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
  * License: GNU Lesser General Public License 2.1 
* GraphStream version 2.0
  * [GraphStream project](https://graphstream-project.org)
  * License: CeCILL-C and GNU Lesser General Public License 3
* JUnit version 4.13.2
  * [JUnit 4](https://junit.org/junit4/)
  * License: Eclipse Public License 1.0
* Google Guava version 30.1.1-jre
  * [Guava](https://guava.dev)
  * License: Apache License 2.0

## References and useful links
**Conceptual Scaffold, Scaffold Tree, and Scaffold Network articles**
* [G. W. Bemis and M. A. Murcko, “The Properties of Known Drugs. 1. Molecular Frameworks,” J. Med. Chem., vol. 39, no. 15, pp. 2887–2893, Jan. 1996, doi: 10.1021/jm9602928.](https://doi.org/10.1021/jm9602928)
* [S. J. Wilkens, J. Janes, and A. I. Su, “HierS: Hierarchical Scaffold Clustering Using Topological Chemical Graphs,” J. Med. Chem., vol. 48, no. 9, pp. 3182–3193, May 2005, doi: 10.1021/jm049032d.](https://doi.org/10.1021/jm049032d)
* [M. A. Koch et al., “Charting biologically relevant chemical space: A structural classification of natural products (SCONP),” Proceedings of the National Academy of Sciences, vol. 102, no. 48, pp. 17272–17277, Nov. 2005, doi: 10.1073/pnas.0503647102.](https://doi.org/10.1073/pnas.0503647102)
* [A. Schuffenhauer, P. Ertl, S. Roggo, S. Wetzel, M. A. Koch, and H. Waldmann, “The Scaffold Tree − Visualization of the Scaffold Universe by Hierarchical Scaffold Classification,” J. Chem. Inf. Model., vol. 47, no. 1, pp. 47–58, Jan. 2007, doi: 10.1021/ci600338x.](https://doi.org/10.1021/ci600338x)
* [T. Varin et al., “Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data,” J. Chem. Inf. Model., vol. 50, no. 12, pp. 2067–2078, Dec. 2010, doi: 10.1021/ci100203e.](https://doi.org/10.1021/ci100203e)
* [T. Varin, A. Schuffenhauer, P. Ertl, and S. Renner, “Mining for Bioactive Scaffolds with Scaffold Networks: Improved Compound Set Enrichment from Primary Screening Data,” J. Chem. Inf. Model., vol. 51, no. 7, pp. 1528–1538, Jul. 2011, doi: 10.1021/ci2000924.](https://doi.org/10.1021/ci2000924)

**Chemistry Development Kit (CDK)**
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
* [Groovy Cheminformatics with the Chemistry Development Kit](https://github.com/egonw/cdkbook)

**GraphStream**
* [GraphStream project](https://graphstream-project.org)
* [Antoine Dutot, Frédéric Guinand, Damien Olivier, Yoann Pigné. GraphStream: A Tool for bridging the gap between Complex Systems and Dynamic Graphs. Emergent Properties in Natural and Artificial Complex Systems. Satellite Conference within the 4th European Conference on Complex Systems (ECCS’2007), Oct 2007, Dresden, Germany. ffhal-00264043](https://hal.archives-ouvertes.fr/hal-00264043/)

**COlleCtion of Open NatUral producTs (COCONUT)**
* [COCONUT Online home page](https://coconut.naturalproducts.net)
* [Sorokina, M., Merseburger, P., Rajan, K. et al. COCONUT online: Collection of Open Natural Products database. J Cheminform 13, 2 (2021). https://doi.org/10.1186/s13321-020-00478-9](https://doi.org/10.1186/s13321-020-00478-9)
* [Sorokina, M., Steinbeck, C. Review on natural products databases: where to find data in 2020. J Cheminform 12, 20 (2020).](https://doi.org/10.1186/s13321-020-00424-9)











## ScaffoldGeneratorTest
__Class ScaffoldGeneratorTest shows the funktions of the ScaffoldGenerator__

  A set of mol files is loaded from the resources folder, processed with the ScaffoldGenerator and saved as image in the generated scaffoldTestOutput. 
  The ScaffoldGenerator method getSchuffenhauerScaffold() transforms the molecules in Schuffenhauer scaffolds.

https://coconut.naturalproducts.net/download

## MurckoFragmenterTest

__Class MurckoFragmenterTest shows the funktions of the CDK MurckoFragmenter__
  
  A set of mol files is loaded from the resources folder and processed with the MurckoFragmenter. The following molecules are saved as image in the generated scaffoldTestOutput:
  * Original: The unchanged molecule
  * Fragments: The fragments generated by the MurckoFragmenter
  * Rings: The rings generated by the MurckoFragmenter
  * Frameworks: The frameworks generated by the MurckoFragmenter
  
  Examples:
  * *Test1* shows that linkers are not further decomposed
  * *Test3* shows that ring systems are separated from one another
  * *Test6* shows that side chains of the ring systems are removed and ring systems are not further broken down into individual rings

## Required ressources
COCONUT DB for speed test: https://coconut.naturalproducts.net/download
Add the COCONUT_DB.sdf file to src\test\resources
