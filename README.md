# Scaffold Generator
##### A CDK-based library for generating Scaffold Trees and Scaffold Networks

## Description
The Scaffold Generator library is designed to make molecular scaffold-related functionalities available in applications 
and workflows based on the [Chemistry Development Kit (CDK)](https://cdk.github.io/). Building mainly upon the works by 
[Bemis and Murcko](https://doi.org/10.1021/jm9602928), [Schuffenhauer et al.](https://doi.org/10.1021/ci600338x), 
and [Varin et al.](https://doi.org/10.1021/ci2000924), it offers scaffold perception and dissection based on single 
molecules and molecule collections. From the latter, Scaffold Trees and Scaffold Networks can be constructed, 
represented in data structures, and visualised using the [GraphStream library](https://graphstream-project.org). 
Multiple options to fine-tune and adapt the routines are available.

## Contents of this repository
### Sources
* ScaffoldGenerator: Core class of the library making its main functionalities available through convenient, high-level methods.
* ScaffoldNodeCollectionBase: Abstract base class of Scaffold Trees and Scaffold Networks. Manages a collection of scaffold nodes and their connections in a tree-based representation that can be exported as an adjacency matrix.
* ScaffoldNetwork: Class used to represent Scaffold Network structures.
* ScaffoldTree: Class used to represent Scaffold Tree structures.
* ScaffoldNodeBase: Abstract base class for nodes in a Scaffold Tree or Scaffold Network that maintains structural information about its scaffold and the molecules it originated from.
* NetworkNode: A class for nodes of a Scaffold Network that can have multiple parents.
* TreeNode: A class for nodes of a Scaffold Tree that can have only one parent.
* GraphStreamUtility: Utility class for displaying Scaffold Trees and Scaffold Networks using the [GraphStream library](https://graphstream-project.org). It is primarily employed in ScaffoldGeneratorTest (see below) for visual inspection.


### Tests
* MurckoFragmenterTest: Illustrates the functions and options of CDK’s MurckoFragmenter that is used in Scaffold Generator for initial scaffold creation. Example molecules from resources (see below) are imported and image files of the generated framework structures saved in an output folder.
* ScaffoldGeneratorTest: Illustrates and tests the functionalities of Scaffold Generator; the correct output of its basic methods like scaffold generation, the more advanced functions to build Scaffold Trees and Scaffold Networks, the correct application of [Schuffenhauer et al.'s](https://doi.org/10.1021/ci600338x) prioritization rules (based on the schemata given in their publication), and the correct workings of the available settings and options. Some examples of Scaffold Trees and Scaffold Networks are displayed for visual inspection using the [GraphStream library](https://graphstream-project.org) and examples for the basic functionalities are visualised using example molecules imported from the resource folder (see below) and saved as image files in an output folder. Additionally, performance tests are included that apply specific routines of Scaffold Generator to the whole [COCONUT database](https://doi.org/10.1186/s13321-020-00478-9).


### Resources
The test resources folder at path <i>src\test\resources\ </i> contains MDL MOL files of 23 test molecules used to illustrate the basic functionalities of Scaffold Generator. They are imported in multiple test methods and the results saved as image files in respective molecule-specific output folders.
<br>An SD file of the [COCONUT database](https://doi.org/10.1186/s13321-020-00478-9) to run the performance tests, is not included in the repository (see below).


## Installation
This is a Maven project. In order to use the source code for your own software, download or clone the repository and
open it in a Maven-supporting IDE (e.g. IntelliJ) as a Maven project and execute the pom.xml file. Maven will then take
care of installing all dependencies. A Java Development Kit (JDK) of version 11 or higher must also be pre-installed.
<br>To run the COCONUT-analysing tests, an SD file of the database needs to be placed in the test "resources" folder
at path <i>src\test\resources\COCONUT_DB.sdf</i>. 
The respective file can be downloaded at [https://coconut.naturalproducts.net/download](https://coconut.naturalproducts.net/download).

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
* [C. Manelfi et al., “‘Molecular Anatomy’: a new multi-dimensional hierarchical scaffold analysis tool,” J Cheminform, vol. 13, no. 1, p. 54, Dec. 2021, doi: 10.1186/s13321-021-00526-y.](https://doi.org/10.1186/s13321-021-00526-y)

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