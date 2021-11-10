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
    private static final String RESULTS_FILE_NAME = "Results.txt";

    /**
     * Name of CSV file for writing time stamps
     */
    private static final String CSV_TIME_FILE_NAME = "CSV_Time_Stamp.csv";

    /**
     * Name of CSV file for the network SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_NETWORK_FILE_NAME = "CSV_Origin_Network.csv";

    /**
     * Name of CSV file for the forest SMILES and the number of origins
     */
    private static final String CSV_ORIGIN_FOREST_FILE_NAME = "CSV_Origin_Forest.csv";
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
        String tmpTimeStamp = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        File tmpExceptionsLogFile = new File(this.workingPath
                + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        PrintWriter tmpExceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        tmpExceptionsPrintWriter.println("#########################################################################");
        tmpExceptionsPrintWriter.println("Time-stamp: " + tmpTimeStamp);
        tmpExceptionsPrintWriter.println();
        tmpExceptionsPrintWriter.flush();
        tmpExceptionsPrintWriter.close();
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
            File tmpResultsLogFile = new File(this.workingPath + PerformanceTest.RESULTS_FILE_NAME);
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            PrintWriter tmpResultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            tmpResultsPrintWriter.println("#########################################################################");
            tmpResultsPrintWriter.println("Time-stamp: " + tmpTimeStamp);
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.println("Application initialized. Loading database file named " + anArgs + ".");
            tmpResultsPrintWriter.flush();
            System.out.println("\nApplication initialized. Loading database file named " + anArgs + ".");
            /*Time stamp file*/
            File tmpCSVTimeFile = new File(this.workingPath + PerformanceTest.CSV_TIME_FILE_NAME);
            FileWriter tmpCSVTimeFileWriter = new FileWriter(tmpCSVTimeFile, false);
            PrintWriter tmpCSVTimePrintWriter = new PrintWriter(tmpCSVTimeFileWriter);
            tmpCSVTimePrintWriter.println("Molecule,Network,Forest");
            tmpCSVTimePrintWriter.flush();
            /*Network origin file*/
            File tmpOriginNetworkOriginFile = new File(this.workingPath + PerformanceTest.CSV_ORIGIN_NETWORK_FILE_NAME);
            FileWriter tmpOriginNetworkOriginFileWriter = new FileWriter(tmpOriginNetworkOriginFile, false);
            PrintWriter tmpOriginNetworkOriginPrintWriter = new PrintWriter(tmpOriginNetworkOriginFileWriter);
            tmpOriginNetworkOriginPrintWriter.println("SMILES,Number of Origins");
            tmpOriginNetworkOriginPrintWriter.flush();
            /*Forest origin file*/
            File tmpOriginForestOriginFile = new File(this.workingPath + PerformanceTest.CSV_ORIGIN_FOREST_FILE_NAME);
            FileWriter tmpOriginForestOriginFileWriter = new FileWriter(tmpOriginForestOriginFile, false);
            PrintWriter tmpOriginForestOriginPrintWriter = new PrintWriter(tmpOriginForestOriginFileWriter);
            tmpOriginForestOriginPrintWriter.println("SMILES,Number of Origins");
            tmpOriginForestOriginPrintWriter.flush();
            IteratingSDFReader tmpDBReader = new IteratingSDFReader(tmpDBFileInputStream, DefaultChemObjectBuilder.getInstance(), true);
            /*Load molecules*/
            List<IAtomContainer> tmpMoleculesList = new LinkedList<>();
            while (tmpDBReader.hasNext()) {
                try {
                    IAtomContainer tmpMolecule = (IAtomContainer) tmpDBReader.next();
                    tmpMoleculesList.add(tmpMolecule);
                } catch (Exception anException) {
                    this.appendToLogfile(anException);
                }
            }
            try {
                tmpDBReader.close();
            } catch (IOException ex) {
                this.appendToLogfile(ex);
            }
            tmpResultsPrintWriter.println("\nDone Loading database. Found and processed " + tmpMoleculesList.size() + " valid molecules.");
            System.out.println("Done Loading database. Found and processed " + tmpMoleculesList.size() + " valid molecules.");
            /*Remove all molecules with more than 10 Rings from list*/
            ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
            for(int tmpIndex = 0 ; tmpIndex < tmpMoleculesList.size(); tmpIndex++) {
                System.out.println(tmpIndex);
                IAtomContainer tmpScaffold = tmpScaffoldGenerator.getScaffold(tmpMoleculesList.get(tmpIndex), false);
                if(tmpScaffoldGenerator.getRings(tmpScaffold, false, false).size() > 10) {
                    tmpMoleculesList.remove(tmpIndex);
                    tmpIndex--;
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
                tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
            }
            long tmpEndTime = System.currentTimeMillis();
            tmpResultsPrintWriter.println("Schuffenhauer fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            System.out.println("Schuffenhauer fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            /*Generate all Schuffenhauer fragments of the molecule*/
            tmpResultsPrintWriter.println("\nGenerate all enumerative removal fragments.");
            System.out.println("Generate all enumerative removal fragments.");
            tmpStartTime = System.currentTimeMillis();
            /*Real measured process*/
            for(IAtomContainer tmpMolecule : tmpMoleculesList) {
                tmpScaffoldGenerator.applyEnumerativeRemoval(tmpMolecule);
            }
            tmpEndTime = System.currentTimeMillis();
            /*Log the skipped molecules*/
            Logger tmpLogger = Logger.getLogger("");
            FileHandler tmpFileHandler = null;
            try {
                tmpFileHandler = new FileHandler("Exceptions_Log.txt", true);
            } catch (IOException e) {
                e.printStackTrace();
            }
            tmpLogger.addHandler(tmpFileHandler);
            tmpResultsPrintWriter.println("Enumerative removal fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            System.out.println("Enumerative removal fragment generation took " + (tmpEndTime - tmpStartTime) + " ms.");
            /*Generate 100 random subsets from 1/100 to 10/100 and generate a network and a forest out of them*/
            int tmpNumberOfRounds = 10;
            for (int tmpRound = 1; tmpRound < (tmpNumberOfRounds + 1); tmpRound++) {
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
                long tmpNetworkOrigin = tmpEndTime - tmpStartTime;
                tmpResultsPrintWriter.println("Network generation took " + tmpNetworkOrigin + " ms.");
                System.out.println("Network generation took " + tmpNetworkOrigin + " ms.");
                tmpStartTime = System.currentTimeMillis();
                /*Real measured process*/
                List<ScaffoldTree> tmpTreeList = tmpScaffoldGenerator.generateSchuffenhauerForest(tmpMoleculeSubList);
                tmpEndTime = System.currentTimeMillis();
                long tmpForestOrigin = tmpEndTime - tmpStartTime;
                tmpResultsPrintWriter.println("Forest generation took " + tmpForestOrigin + " ms.");
                System.out.println("Forest generation took " + tmpForestOrigin + " ms.");
                tmpCSVTimePrintWriter.println(tmpRate + "," + tmpNetworkOrigin + "," + tmpForestOrigin);
                tmpCSVTimePrintWriter.flush();
                if(tmpRound == tmpNumberOfRounds) {
                    SmilesGenerator tmpSmilesGenerator = new SmilesGenerator(SmiFlavor.Stereo);
                    tmpResultsPrintWriter.println("Number of NetworkNodes: " + tmpScaffoldNetwork.getAllNodes().size());
                    System.out.println("Number of NetworkNodes: " +tmpScaffoldNetwork.getAllNodes().size());
                    for(ScaffoldNodeBase tmpNode : tmpScaffoldNetwork.getAllNodes()) {
                        String tmpSMILES = tmpSmilesGenerator.create((IAtomContainer) tmpNode.getMolecule());
                        tmpOriginNetworkOriginPrintWriter.println(tmpSMILES + "," +  tmpNode.getOriginCount());
                        tmpOriginNetworkOriginPrintWriter.flush();
                    }
                    tmpResultsPrintWriter.println("Number of trees: " + tmpTreeList.size());
                    System.out.println("Number of trees: " + tmpTreeList.size());
                    long tmpNumberOfTreeNodes = 0;
                    for(ScaffoldTree tmpTree : tmpTreeList) {
                        tmpNumberOfTreeNodes = tmpNumberOfTreeNodes + tmpTree.getAllNodes().size();
                        for(ScaffoldNodeBase tmpNode : tmpTree.getAllNodes()) {
                            String tmpSMILES = tmpSmilesGenerator.create((IAtomContainer) tmpNode.getMolecule());
                            tmpOriginForestOriginPrintWriter.println(tmpSMILES + "," +  tmpNode.getOriginCount());
                            tmpOriginForestOriginPrintWriter.flush();
                        }
                    }
                    tmpResultsPrintWriter.println("Number of tree nodes: " + tmpNumberOfTreeNodes);
                    System.out.println("Number of tree nodes: " + tmpNumberOfTreeNodes);
                }
            }
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
                    + PerformanceTest.EXCEPTIONS_LOG_FILE_NAME,
                    true);
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
