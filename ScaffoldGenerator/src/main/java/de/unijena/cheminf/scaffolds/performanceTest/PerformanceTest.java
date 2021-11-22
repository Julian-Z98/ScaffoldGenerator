 /**
 * Performance test for
 * ScaffoldGenerator for CDK
 * Copyright (C) 2021 Julian Zander
 *
 * Source code is available at <https://github.com/Julian-Z98/ScaffoldGenerator>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.unijena.cheminf.scaffolds.performanceTest;

import de.unijena.cheminf.scaffolds.ScaffoldGenerator;
import de.unijena.cheminf.scaffolds.ScaffoldNetwork;
import de.unijena.cheminf.scaffolds.ScaffoldNodeBase;
import de.unijena.cheminf.scaffolds.ScaffoldTree;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

/**
 * An application for testing the performance of the ScaffoldGenerator methods.
 *
 * @author Julian Zander(zanderjulian@gmx.de)
 */
public class PerformanceTest {

    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    /**
     * Name of file for logging occurred exceptions
     */
    private static final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log.txt";

    /**
     * Name of file for writing results
     */
    private static final String RESULTS_FILE_NAME = "Results";

    /**
     * Name of CSV file for writing time processing time
     */
    private static final String CSV_TIME_FILE_NAME = "CSV_Processing_Time";

    /**
     * Name of CSV file for the network SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_NETWORK_FILE_NAME = "CSV_Origin_Network.";

    /**
     * Name of CSV file for the forest SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_FOREST_FILE_NAME = "CSV_Origin_Forest";
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Instantiates and starts the application.
     * It first loads all molecules from a given SD file into memory And then applies the methods of the ScaffoldGenerator to them.
     *
     * First, all molecules with more than 10 rings are removed.
     * The Schuffenhauer fragments and the enumerative fragments of all molecules are then generated and the time of generation is measured.
     * The molecules are then divided into 100 subsets in ascending order of size and the generation of networks and forests is performed with these subsets.
     * The size of the subsets is 1/100 to 100/100 of the dataset and the respective runtime for each subset is stored in a CSV file.
     *
     * Results and logs are generated and saved in different txt and csv files in the .jar folder:
     * -Results.txt contains the summary of all tests run and their runtimes
     * -Exceptions_Log.txt contains the exceptions that may have been thrown together with their description
     * -CSV_Origin_Forest.csv contains the SMILES of each node of the forest and the number of their origins.
     * -CSV_Origin_Network.csv contains the SMILES of each node of the network and the number of their origins.
     * -CSV_Processing_Time.csv contains the size of all subsets that were used to generate a forest/network and their specific runtime in ms.
     *
     * It exits the system if an unexpected exception occurs that prevents the application from working, e.g. an IllegalArgumentException
     * (will be logged to a file, not printed on the console).
     *
     * @param anArgs the command line arguments, anArgs[0] must be the name of the SD file to load (must be located in
     * the same directory as the application's JAR file)
     * @throws IOException if the constructor is unable to open a text file for logging occurred exceptions
     */
    public PerformanceTest(String anArgs) throws IOException {
        /*Set up exception log file*/
        this.workingPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpProcessingTime = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        File tmpExceptionsLogFile = new File(this.workingPath
                + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        PrintWriter tmpExceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        tmpExceptionsPrintWriter.println("#########################################################################");
        tmpExceptionsPrintWriter.println("Processing Time: " + tmpProcessingTime);
        tmpExceptionsPrintWriter.println();
        tmpExceptionsPrintWriter.flush();
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
        try {
            /*Load SD file*/
            File tmpDBFile = new File(this.workingPath + anArgs);
            FileInputStream tmpDBFileInputStream;
            try {
                tmpDBFileInputStream = new FileInputStream(tmpDBFile);
            } catch (FileNotFoundException | SecurityException anException) {
                throw new IllegalArgumentException("The database file (name) is not valid: " + anException.getMessage());
            }
            /*Results file*/
            File tmpResultsLogFile = new File(this.workingPath + PerformanceTest.RESULTS_FILE_NAME + tmpProcessingTime + ".txt");
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            PrintWriter tmpResultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            tmpResultsPrintWriter.println("#########################################################################");
            tmpResultsPrintWriter.println("Processing Time: " + tmpProcessingTime);
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.println("Application initialized. Loading database file named " + anArgs + ".");
            tmpResultsPrintWriter.flush();
            System.out.println("\nApplication initialized. Loading database file named " + anArgs + ".");
            /*ProcessingTime file*/
            File tmpCSVTimeFile = new File(this.workingPath + PerformanceTest.CSV_TIME_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpCSVTimeFileWriter = new FileWriter(tmpCSVTimeFile, false);
            PrintWriter tmpCSVTimePrintWriter = new PrintWriter(tmpCSVTimeFileWriter);
            tmpCSVTimePrintWriter.println("Number of the molecules,Calculation time for the network in ms,Calculation time for the forest in ms");
            tmpCSVTimePrintWriter.flush();
            /*Network origin file*/
            File tmpOriginNetworkOriginFile = new File(this.workingPath + PerformanceTest.CSV_ORIGIN_NETWORK_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpOriginNetworkOriginFileWriter = new FileWriter(tmpOriginNetworkOriginFile, false);
            PrintWriter tmpOriginNetworkOriginPrintWriter = new PrintWriter(tmpOriginNetworkOriginFileWriter);
            tmpOriginNetworkOriginPrintWriter.println("SMILES of the node,Number of origins");
            tmpOriginNetworkOriginPrintWriter.flush();
            /*Forest origin file*/
            File tmpOriginForestOriginFile = new File(this.workingPath + PerformanceTest.CSV_ORIGIN_FOREST_FILE_NAME + tmpProcessingTime + ".csv");
            FileWriter tmpOriginForestOriginFileWriter = new FileWriter(tmpOriginForestOriginFile, false);
            PrintWriter tmpOriginForestOriginPrintWriter = new PrintWriter(tmpOriginForestOriginFileWriter);
            tmpOriginForestOriginPrintWriter.println("SMILES of the node,Number of origins");
            tmpOriginForestOriginPrintWriter.flush();
            IteratingSDFReader tmpDBReader = new IteratingSDFReader(tmpDBFileInputStream, DefaultChemObjectBuilder.getInstance(), true);
            /*Load molecules*/
            List<IAtomContainer> tmpMoleculesList = new LinkedList<>();
            while (tmpDBReader.hasNext()) {
                try {
                    IAtomContainer tmpMolecule = (IAtomContainer) tmpDBReader.next();
                    tmpMoleculesList.add(tmpMolecule);
                } catch (Exception anException) {
                    tmpExceptionsPrintWriter.println("Load molecules ERROR");
                    tmpExceptionsPrintWriter.flush();
                    this.appendToLogfile(anException);
                }
            }
            try {
                tmpDBReader.close();
            } catch (IOException anException) {
                tmpExceptionsPrintWriter.println("DBReader close IO ERROR");
                tmpExceptionsPrintWriter.flush();
                this.appendToLogfile(anException);
            }
            tmpResultsPrintWriter.println("\nDone Loading database. Found and processed " + tmpMoleculesList.size() + " valid molecules.");
            System.out.println("Done Loading database. Found and processed " + tmpMoleculesList.size() + " valid molecules.");
            /*Remove all molecules with more than 10 rings from list*/
            ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
            tmpScaffoldGenerator.setSmilesGeneratorSetting(tmpSmilesGenerator);
            for(int tmpIndex = 0 ; tmpIndex < tmpMoleculesList.size(); tmpIndex++) {
                try {
                    if(tmpScaffoldGenerator.getRings(tmpMoleculesList.get(tmpIndex), false, false).size() > 10) {
                        tmpMoleculesList.remove(tmpIndex);
                        tmpIndex--;
                    }
                } catch(Exception anException) {
                    tmpExceptionsPrintWriter.println("Remove molecules with more than 10 rings ERROR");
                    tmpExceptionsPrintWriter.flush();
                    this.appendToLogfile(anException);
                }
            }
            int tmpListSize = tmpMoleculesList.size();
            tmpResultsPrintWriter.println("\nNumber of molecules with less than 11 rings: " + tmpListSize);
            System.out.println("Number of molecules with less than 11 rings: " + tmpListSize);
            /*Generate all Schuffenhauer fragments of the molecule*/
            tmpResultsPrintWriter.println("\nGenerate all Schuffenhauer fragments.");
            System.out.println("Generate all Schuffenhauer fragments.");
            long tmpStartTime = System.currentTimeMillis();
            /*Real measured process*/
            for(IAtomContainer tmpMolecule : tmpMoleculesList) {
                try {
                    tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
                } catch (Exception anException) {
                    try {
                        tmpExceptionsPrintWriter.println("applySchuffenhauerRules ERROR with molecule: " + tmpSmilesGenerator.create(tmpMolecule));
                        tmpExceptionsPrintWriter.flush();
                        this.appendToLogfile(anException);
                    } catch(Exception anExceptionException) {
                        tmpExceptionsPrintWriter.println("applySchuffenhauerRules ERROR. Cant generate SMILES of the molecule");
                        tmpExceptionsPrintWriter.flush();
                        this.appendToLogfile(anException);
                    }
                }
            }
            long tmpEndTime = System.currentTimeMillis();
            tmpResultsPrintWriter.println("Schuffenhauer fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            System.out.println("Schuffenhauer fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            /*Generate all enumerative fragments of the molecule*/
            tmpResultsPrintWriter.println("\nGenerate all enumerative removal fragments.");
            System.out.println("Generate all enumerative removal fragments.");
            tmpStartTime = System.currentTimeMillis();
            /*Real measured process*/
            for(IAtomContainer tmpMolecule : tmpMoleculesList) {
                try {
                    tmpScaffoldGenerator.applyEnumerativeRemoval(tmpMolecule);
                } catch (Exception anException) {
                    try {
                        tmpExceptionsPrintWriter.println("applyEnumerativeRemoval ERROR with molecule: " + tmpSmilesGenerator.create(tmpMolecule));
                        tmpExceptionsPrintWriter.flush();
                        this.appendToLogfile(anException);
                    } catch(Exception anExceptionException) {
                        tmpExceptionsPrintWriter.println("applyEnumerativeRemoval ERROR. Cant generate SMILES of the molecule");
                        tmpExceptionsPrintWriter.flush();
                        this.appendToLogfile(anException);
                    }
                }
            }
            tmpEndTime = System.currentTimeMillis();
            tmpResultsPrintWriter.println("Enumerative removal fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            System.out.println("Enumerative removal fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            /*Log the skipped molecules*/
            Logger tmpLogger = Logger.getLogger("");
            FileHandler tmpFileHandler = null;
            try {
                tmpFileHandler = new FileHandler("Exceptions_Log.txt", true);
                tmpLogger.addHandler(tmpFileHandler);
            } catch(Exception anException) {
                tmpExceptionsPrintWriter.println("FileHandler ERROR");
                tmpExceptionsPrintWriter.flush();
                this.appendToLogfile(anException);
            }
            /*Generate 100 random subsets from 1/100 to 100/100 and generate a network and a forest out of them*/









            //Change to 100
            int tmpNumberOfRounds = 1;
            for (int tmpRound = 1; tmpRound < (tmpNumberOfRounds + 1); tmpRound++) {
                try {
                    int tmpRate = (int) (tmpListSize / (float)tmpNumberOfRounds * tmpRound);
                    tmpResultsPrintWriter.println("\nProcess " + tmpRate + " valid molecules.");
                    System.out.println("Process " + tmpRate + " valid molecules.");
                    Collections.shuffle(tmpMoleculesList, new Random(42));
                    List<IAtomContainer> tmpMoleculeSubList = tmpMoleculesList.subList(0, tmpRate);
                    tmpResultsPrintWriter.flush();
                    tmpStartTime = System.currentTimeMillis();
                    /*Real measured process*/
                    ScaffoldNetwork tmpScaffoldNetwork = tmpScaffoldGenerator.generateScaffoldNetwork(tmpMoleculeSubList);
                    tmpEndTime = System.currentTimeMillis();
                    long tmpNetworkOriginTime = tmpEndTime - tmpStartTime;
                    tmpResultsPrintWriter.println("Network generation took " + tmpNetworkOriginTime + " ms.");
                    System.out.println("Network generation took " + tmpNetworkOriginTime + " ms.");
                    tmpStartTime = System.currentTimeMillis();
                    /*Real measured process*/
                    List<ScaffoldTree> tmpTreeList = tmpScaffoldGenerator.generateSchuffenhauerForest(tmpMoleculeSubList);
                    tmpEndTime = System.currentTimeMillis();
                    long tmpForestOriginTime = tmpEndTime - tmpStartTime;
                    tmpResultsPrintWriter.println("Forest generation took " + tmpForestOriginTime + " ms.");
                    System.out.println("Forest generation took " + tmpForestOriginTime + " ms.");
                    tmpCSVTimePrintWriter.println(tmpRate + "," + tmpNetworkOriginTime + "," + tmpForestOriginTime);
                    tmpCSVTimePrintWriter.flush();
                    if(tmpRound == tmpNumberOfRounds) {
                        try {
                            List<ScaffoldNodeBase> tmpNetworkNodeList = tmpScaffoldNetwork.getAllNodes();
                            tmpResultsPrintWriter.println("Number of NetworkNodes: " + tmpNetworkNodeList.size());
                            System.out.println("Number of NetworkNodes: " + tmpNetworkNodeList.size());
                            for(ScaffoldNodeBase tmpNode : tmpNetworkNodeList) {
                                try {
                                    String tmpSMILES = tmpSmilesGenerator.create((IAtomContainer) tmpNode.getMolecule());
                                    tmpOriginNetworkOriginPrintWriter.println(tmpSMILES + "," +  tmpNode.getOriginCount());
                                    tmpOriginNetworkOriginPrintWriter.flush();
                                } catch (Exception anException) {
                                    tmpExceptionsPrintWriter.println("SmilesGenerator ERROR. Cant generate SMILES of the network node");
                                    tmpExceptionsPrintWriter.flush();
                                    this.appendToLogfile(anException);
                                }
                            }
                        } catch (Exception anException) {
                            tmpExceptionsPrintWriter.println("Network node list ERROR");
                            tmpExceptionsPrintWriter.flush();
                            this.appendToLogfile(anException);
                        }
                        tmpResultsPrintWriter.println("Number of trees: " + tmpTreeList.size());
                        System.out.println("Number of trees: " + tmpTreeList.size());
                        try {
                            long tmpNumberOfTreeNodes = 0;
                            for(ScaffoldTree tmpTree : tmpTreeList) {
                                ScaffoldNodeBase tmpRootNode = tmpTree.getRoot();
                                tmpNumberOfTreeNodes = tmpNumberOfTreeNodes + tmpTree.getAllNodes().size();
                                for(ScaffoldNodeBase tmpNode : tmpTree.getAllNodes()) {
                                    try {
                                        String tmpSMILES = tmpSmilesGenerator.create((IAtomContainer) tmpNode.getMolecule());
                                        tmpOriginForestOriginPrintWriter.print(tmpSMILES + "," +  tmpNode.getOriginCount());
                                        if(tmpNode.equals(tmpRootNode)) {
                                            if(tmpTree.getAllNodes().size() < 50) {
                                                for(ScaffoldNodeBase tmpOriginNode : tmpTree.getAllNodes()) {
                                                    for(Object tmpOrigin : tmpOriginNode.getNonVirtualOriginSmilesList()) {
                                                        tmpOriginForestOriginPrintWriter.print("," + tmpOrigin);
                                                    }
                                                }
                                            } else {
                                                tmpOriginForestOriginPrintWriter.print("," + "-");
                                            }
                                        }
                                        tmpOriginForestOriginPrintWriter.println();
                                        tmpOriginForestOriginPrintWriter.flush();
                                    } catch (Exception anException) {
                                        tmpExceptionsPrintWriter.println("SmilesGenerator ERROR. Cant generate SMILES of the forest node");
                                        tmpExceptionsPrintWriter.flush();
                                        this.appendToLogfile(anException);

                                    }
                                }
                            }
                            tmpResultsPrintWriter.println("Number of tree nodes: " + tmpNumberOfTreeNodes);
                            System.out.println("Number of tree nodes: " + tmpNumberOfTreeNodes);
                        } catch (Exception anException) {
                            tmpExceptionsPrintWriter.println("Forest node list ERROR");
                            tmpExceptionsPrintWriter.flush();
                            this.appendToLogfile(anException);
                        }
                    }
                } catch(Exception anException) {
                    tmpExceptionsPrintWriter.println("Subset ERROR");
                    tmpExceptionsPrintWriter.flush();
                    this.appendToLogfile(anException);
                }
            }
            String tmpFinTime = LocalDateTime.now().format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
            tmpResultsPrintWriter.println("---FIN at " + tmpFinTime  + "---");
            tmpResultsPrintWriter.flush();
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.flush();
            tmpResultsPrintWriter.close();
            tmpCSVTimePrintWriter.flush();
            tmpCSVTimePrintWriter.close();
            tmpOriginNetworkOriginPrintWriter.flush();
            tmpOriginNetworkOriginPrintWriter.close();
            tmpOriginForestOriginPrintWriter.flush();
            tmpOriginForestOriginPrintWriter.close();
            tmpFileHandler.close();
            tmpExceptionsPrintWriter.close();
        } catch (Exception anException) {
            this.appendToLogfile(anException);
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }
    //</editor-fold>

    /**
     * Appends the given exception's stack trace to a log file.
     *
     * @param anException the exception to log
     */
    private void appendToLogfile(Exception anException) {
        if (anException == null) {
            return;
        }
        PrintWriter tmpPrintWriter = null;
        try {
            FileWriter tmpFileWriter = new FileWriter(this.workingPath
                    + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME, true);
            tmpPrintWriter = new PrintWriter(tmpFileWriter);
            StringWriter tmpStringWriter = new StringWriter();
            anException.printStackTrace(new PrintWriter(tmpStringWriter));
            String tmpStackTrace = tmpStringWriter.toString();
            tmpPrintWriter.println(tmpStackTrace);
        } catch (IOException anIOException) {
            anIOException.printStackTrace(System.err);
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.println();
                tmpPrintWriter.flush();
                tmpPrintWriter.close();
            }
        }
    }
}
