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
import de.unijena.cheminf.scaffolds.ScaffoldTree;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

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
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    //</editor-fold>

    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Instantiates and starts the application. It first loads all molecules from a given SD file into memory and then
     * distributes them equally on the given number of different threads to use. It measures the time it takes for all
     * threads to complete the extraction of functional groups using the ScaffoldGenerator. It exits the system
     * if an unexpected exception occurs that prevents the application from working, e.g. an IllegalArgumentException
     * (will be logged to a file, not printed on the console).
     *
     * @param anArgs the command line arguments, anArgs[0] must be the name of the SD file to load (must be located in
     * the same directory as the application's JAR file) and anArgs[1] must be the number of different threads to use
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

            File tmpResultsLogFile = new File(this.workingPath
                    + PerformanceTest.RESULTS_FILE_NAME);
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            PrintWriter tmpResultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            tmpResultsPrintWriter.println("#########################################################################");
            tmpResultsPrintWriter.println("Time-stamp: " + tmpTimeStamp);
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.println("Application initialized. Loading database file named " + anArgs + ".");
            tmpResultsPrintWriter.flush();
            System.out.println("\nApplication initialized. Loading database file named " + anArgs + ".");
            IteratingSDFReader tmpDBReader = new IteratingSDFReader(tmpDBFileInputStream, DefaultChemObjectBuilder.getInstance(), true);
            /*Load molecules*/
            List<IAtomContainer> tmpMoleculesList = new LinkedList<>();
            while (tmpDBReader.hasNext()) {
                try {
                    IAtomContainer tmpMolecule = (IAtomContainer) tmpDBReader.next();
                    tmpMoleculesList.add(tmpMolecule);
                } catch (Exception anException) {
                    /*If an IllegalArgumentException is thrown in applyFiltersAndPreprocessing (meaning that the molecule
                    should be filtered) the molecule is skipped by catching this exception*/
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
            for (int tmpDeci = 1; tmpDeci < 101; tmpDeci++) {
                int tmpRate = (int) (tmpListSize / 100.0 * tmpDeci);
                tmpResultsPrintWriter.println("\nProcess " + tmpRate + " valid molecules.");
                System.out.println("Process " + tmpRate + " valid molecules.");
                long tmpSeed = System.nanoTime();
                Collections.shuffle(tmpMoleculesList, new Random(tmpSeed));
                List<IAtomContainer> tmpMoleculeSubList = tmpMoleculesList.subList(0, tmpRate);
                tmpResultsPrintWriter.flush();
                tmpStartTime = System.currentTimeMillis();
                /*Real measured process*/
                tmpScaffoldGenerator.generateScaffoldNetwork(tmpMoleculeSubList);
                tmpEndTime = System.currentTimeMillis();
                tmpResultsPrintWriter.println("Network generation took " + (tmpEndTime - tmpStartTime) + " ms.");
                System.out.println("Network generation took " + (tmpEndTime - tmpStartTime) + " ms.");
                tmpStartTime = System.currentTimeMillis();
                /*Real measured process*/
                List<ScaffoldTree> tmpTreeList = tmpScaffoldGenerator.generateSchuffenhauerForest(tmpMoleculeSubList);
                tmpEndTime = System.currentTimeMillis();
                tmpResultsPrintWriter.println("Forest generation took " + (tmpEndTime - tmpStartTime) + " ms.");
                System.out.println("Forest generation took " + (tmpEndTime - tmpStartTime) + " ms.");
                System.out.println("Number of trees: " + tmpTreeList.size());
            }
            tmpResultsPrintWriter.flush();
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.flush();
            tmpResultsPrintWriter.close();
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
