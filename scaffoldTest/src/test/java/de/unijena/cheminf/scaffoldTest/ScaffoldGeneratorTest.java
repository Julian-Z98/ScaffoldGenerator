/*
 * Copyright (c) 2021 Julian Zander, Jonas Schaub,  Achim Zielesny
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, version 2.1.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */

package de.unijena.cheminf.scaffoldTest;

import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.TimeUnit;

import static junit.framework.TestCase.assertEquals;

/**
 * JUnit test class for the ScaffoldGenerator
 */
public class ScaffoldGeneratorTest {
    private ScaffoldGenerator scaffoldGenerator;
    @Before
    public void initFragmenter() {
        scaffoldGenerator = new ScaffoldGenerator();
    }

    //<editor-fold desc="Fundamental Method Tests">
    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold() with V2000 and V3000 mol files.
     * Loads the 22 Test(Test1.mol-Test23.mol) molfiles from the Resources folder and creates the SchuffenhauerScaffolds with getSchuffenhauerScaffold().
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerScaffoldTest() throws CloneNotSupportedException, CDKException, IOException {
        for (int tmpCount = 1; tmpCount < 23; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            /*Generate picture of the original*/
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
            BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
            /*Save the original picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Original.png").mkdirs();
            File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Original.png");
            ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
            //Generate picture of the SchuffenhauerScaffold
            BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png").mkdirs();
            File tmpOutputSchuff = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png");
            ImageIO.write(tmpImgSchuff, "png" ,tmpOutputSchuff);
        }
    }
    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold with SMILES.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerNonCTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("[S+]#S1CCCCC1"); //Triple bond
        //IAtomContainer tmpMolecule = tmpParser.parseSmiles("P=[P+]1CCCCC1"); // P=P Test
        tmpMolecule = tmpParser.parseSmiles("S=P12N=P3(OC=4C=CC(C=NN(C)P(=S)(C=5C=CC=CC5)N(N=CC=6C=CC(OP(=NP(=S)(OC=7C=CC(C=NN(C)P(=S)(C=8C=CC=CC8)N(N=CC=9C=CC(O1)=CC9)C)=CC7)OC=%10C=CC(C=NN(C)P(=S)(C=%11C=CC=CC%11)N(N=CC=%12C=CC(O2)=CC%12)C)=CC%10)(OC=%13C=CC(C=NN(C)P(=S)(C=%14C=CC=CC%14)N(N=CC=%15C=CC(O3)=CC%15)C)=CC%13)C=%16C=CC=CC%16)=CC6)C)=CC4)C=%17C=CC=CC%17"); // Test
        //IAtomContainer tmpMolecule = tmpParser.parseSmiles("S=S1CCCCC1"); //S=S Test
        /*Generate picture of the Original*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold*/
        IAtomContainer tmpSchuff = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpSchuff).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Schuffenhauer.png").mkdirs();
        File tmpOutputSchuff = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Schuffenhauer.png");
        ImageIO.write(tmpImgSchuff, "png" ,tmpOutputSchuff);
        //Generate rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk()),true);
        /*Generate pictures of the rings*/
        int tmpCounter = 1;
        for (IAtomContainer tmpRing : tmpRings) {
            BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Ring" + tmpCounter + ".png").mkdirs();
            File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCDoubleBond/Ring" + tmpCounter + ".png");
            ImageIO.write(tmpImgRing, "png", tmpOutputRing);
            tmpCounter++;
        }
    }

    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold with SMILES.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void removeRingNonCTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("c2cc[p+]1CCCCc1c2"); //Triple bond
        tmpMolecule = tmpParser.parseSmiles("c1cc2CCCc3cccc(c1)c23");//Test
        //tmpMolecule = tmpParser.parseSmiles("c1ccc3c(c1)c2ccccc2c4ccccc34");//Test
        //tmpMolecule = tmpParser.parseSmiles("c1ccc2c(c1)c7cccc6CC3CCc4cccc5cc2c(c3c45)c67");//Test
        tmpMolecule = tmpParser.parseSmiles("CCN(C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C#N)C(=O)C");//Test
        //tmpMolecule = tmpParser.parseSmiles("c7ccc(c5c(c1ccccc1)c(c2ccccc2)c(c3ccccc3)c(c4ccccc4)c5c6ccccc6)cc7");//Test
        //tmpMolecule = tmpParser.parseSmiles("c1c3CCC4CCC5CCC6CCc7cc2CCc1c2c8c3C4=C5C6c78");//Test
        /*Generate picture of the Original*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold*/
        IAtomContainer tmpSchuff = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpSchuff).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Schuffenhauer.png").mkdirs();
        File tmpOutputSchuff = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Schuffenhauer.png");
        ImageIO.write(tmpImgSchuff, "png" ,tmpOutputSchuff);
        //Generate rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk()),true);
        /*Generate pictures of the rings*/
        int tmpCounter = 1;
        for (IAtomContainer tmpRing : tmpRings) {
            IAtomContainer tmpRingRemoved = scaffoldGenerator.removeRing(tmpSchuff, tmpRing);
            BufferedImage tmpImgRingRemoved = tmpGenerator.depict(tmpRingRemoved).toImg();
            BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Ring" + tmpCounter + "Removed.png").mkdirs();
            File tmpOutputRingRemoved = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Ring" + tmpCounter + "Removed.png");
            ImageIO.write(tmpImgRingRemoved, "png", tmpOutputRingRemoved);
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Ring" + tmpCounter + ".png").mkdirs();
            File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/NonCRemove/Ring" + tmpCounter + ".png");
            ImageIO.write(tmpImgRing, "png", tmpOutputRing);
            tmpCounter++;
        }
    }
    /**
     * Test of Cycles.mcb() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates the rings of the SchuffenhauerScaffold with getRings().
     * All generated Rings are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRingsTest() throws IOException, CDKException, CloneNotSupportedException {
        for (int tmpCount = 1; tmpCount < 23; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate the SchuffenhauerScaffold
            tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
            //Generate rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
            /*Generate pictures of the rings*/
            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                /*Save the picture*/
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/GeneratedRing" + tmpCounter + ".png").mkdirs();
                File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/GeneratedRing" + tmpCounter + ".png");
                ImageIO.write(tmpImgRing, "png", tmpOutputRing);
                tmpCounter++;
            }
        }
    }
    /**
     * Test of removeRing() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates for each generated ring, the corresponding total molecule with removed ring.
     * All generated molecules are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void removeRingTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 2; tmpCount < 23; tmpCount++) {
            String tmpFileName = "Test" + tmpCount ;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
            //Generate Rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                /*Generate SchuffenhauerScaffold with removed ring*/
                //tmpSchuffenhauerScaffold = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                /*Generate picture of the SchuffenhauerScaffold with removed ring*/
                DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
                /*Save the picture*/
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/RingRemove" + tmpCounter + ".png").mkdirs();
                File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/RingRemove" + tmpCounter + ".png");
                ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                tmpCounter++;
            }
        }
    }

    /**
     * Example of the errors caused by AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms().
     * All molecules containing B throw a NullPointerException in line 347 of the CDKAtomMatcher.
     * All molecules containing P with 4 bonds throw a NullPointerException in line 1373 of the CDKAtomMatcher.
     * In both cases Atom.getFormalCharge() = null there and triggers this error.
     * In this method, both an example molecule with B and with P are given. Both trigger an error. With tmpBypassError the FormalCharges can be set to 0 and the errors can be bypassed.
     */
    @Test
    public void percieveAtomTypesErrorTest() throws CDKException {
        boolean tmpBypassError = false; //Pass the error or not
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = null;
        //P Examples:
        tmpMolecule = tmpParser.parseSmiles("C1C2CC3CC1CC(C2)C3(CC(=O)C(=P(C4=CC=CC=C4)(C5=CC=CC=C5)C6=CC=CC=C6)C#N)C7=CC=C(C=C7)F"); //PubChem CID: 89041793
        //tmpMolecule = tmpParser.parseSmiles("C1=CC=C(C=C1)C(=O)OCCCC[P+](C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4"); //PubChem CID: 2755300
        //tmpMolecule = tmpParser.parseSmiles("C1CCC(=C(C(=P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C(=O)C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)C1"); //PubChem CID: 2752540
        //B Examples:
        //tmpMolecule = tmpParser.parseSmiles("[B-](C1=CC=CC=C1)(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4"); //PubChem CID: 8934
        //tmpMolecule = tmpParser.parseSmiles("[B-]123OC4C(=O)OC(CCC=CC=CCCC(CC(=O)C(C5CCC(C(O1)(O5)C(O2)C(=O)OC(CCC=CC=CCCC(CC(=O)C(C6CCC(C4(O3)O6)(C)O)C)O)C)(C)O)C)O)C"); //PubChem CID: 637168
        //tmpMolecule = tmpParser.parseSmiles("[B-]123OC4C(=O)OC5CC(C=CCC(C(C6CCC(C(O1)(O6)C(O2)C(=O)OC7CC(C=CCC(C(C8CCC(C4(O3)O8)C)(C)C)O)OC7C)C)(C)C)O)OC5C"); //PubChem CID: 43587
        //AtomContainerManipulator.clearAtomConfigurations(tmpMolecule);
        for(IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpAtom.setHybridization((IAtomType.Hybridization) CDKConstants.UNSET);
        }
        if(tmpBypassError) {
            /*Avoid the error by setting the FormalCharge to 0*/
            for(IAtom tmpAtom : tmpMolecule.atoms()) {
                if(tmpAtom.getSymbol() == "B" || (tmpAtom.getSymbol() == "P" && tmpMolecule.getConnectedBondsList(tmpAtom).size() == 4)) {
                    tmpAtom.setFormalCharge(0);
                }
            }
        }
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule); //Line that causes the error
    }

    /**
     * Test of isRingTerminal() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates for each generated terminal ring, the corresponding total molecule with removed ring.
     * All generated molecules are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void isRingTerminalTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 2; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
            //Generate Rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                //Check that rings are terminal
                if(scaffoldGenerator.isRingTerminal(tmpSchuffenhauerScaffold, tmpRing)) {
                    //Generate SchuffenhauerScaffold with removed ring
                    IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                    /*Generate picture of the SchuffenhauerScaffold with removed ring*/
                    DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                    BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
                    /*Save the picture*/
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TerminalRingRemove" + tmpCounter + ".png").mkdirs();
                    File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TerminalRingRemove" + tmpCounter + ".png");
                    ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                }
                tmpCounter++;
            }
        }
    }
    //</editor-fold>

    //<editor-fold desc="Output Method Test">
    /**
     * Test of getIterativeRemoval() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and iteratively removes the terminal rings.
     * All generated molecules are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getIterativeRemovalTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 1; tmpCount < 24; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate a list of molecules with iteratively removed terminal rings
            List<IAtomContainer> tmpMolecules = scaffoldGenerator.getIterativeRemoval(tmpMolecule, ElectronDonation.cdk());
            int tmpCounter = 1;
            for (IAtomContainer tmpIterative : tmpMolecules) {
                /*Generate picture of the molecule*/
                DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
                /*Save the picture*/
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Iterative" + tmpCounter + ".png").mkdirs();
                File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Iterative" + tmpCounter + ".png");
                ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                tmpCounter++;
            }
        }
    }

    /**
     * Test of getRemovalTree() with V2000 and V3000 mol files.
     * Loads the 13 Test(Test1.mol-Test13.mol) molfiles from the Resources folder, iteratively removes the terminal rings and saves the molecules in a tree.
     * All generated molecules are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRemovalTreeTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 1; tmpCount < 23; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate a tree of molecules with iteratively removed terminal rings
            TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule, ElectronDonation.cdk()));
            int tmpCounter = 0;
            while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
                TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
                /*Save the picture*/
                DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpMoleculeNode.getMolecule()).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TreeTest" + tmpCounter  + "Level" + tmpMoleculeNode.getLevel() + ".png").mkdirs();
                File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TreeTest" + tmpCounter +  "Level" + tmpMoleculeNode.getLevel() + ".png");
                ImageIO.write(tmpSecImgRemove, "png", tmpSecOutputRemove);
                tmpCounter++;
            }
        }
    }

    /**
     * Test of getRemovalTree() with V2000 and V3000 mol files.
     * Loads one molecule(insert in tmpFileName)  from the Resources folder, iteratively removes the terminal rings and saves the molecules in a tree.
     * Saves the parent and the children of one Node and saves them as images.
     * Set file with: tmpFileName
     * Set Node with: tmpTestNumber
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRemovalTreeStructureTest() throws CDKException, CloneNotSupportedException, IOException {
        String tmpFileName = "Test11"; //File to be tested
        int tmpTestNumber = 2; //Node whose children and parent are to be displayed
        //Load molecule from molfile
        IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
        //Generate a tree of molecules with iteratively removed terminal rings
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule, ElectronDonation.cdk()));
        int tmpCounter = 0;
        while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            if(tmpCounter == tmpTestNumber -1){
                /*Save the picture of the test Node*/
                DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                BufferedImage tmpNodeImg = tmpGenerator.depict(tmpMoleculeNode.getMolecule()).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TestNode.png").mkdirs();
                File tmpNodeFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TestNode.png");
                ImageIO.write(tmpNodeImg, "png", tmpNodeFile);
                /*Save the picture of the parent*/
                BufferedImage tmpParentImg = tmpGenerator.depict(tmpMoleculeNode.getParent().getMolecule()).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ParentNode.png").mkdirs();
                File tmpParentFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ParentNode.png");
                ImageIO.write(tmpParentImg, "png", tmpParentFile);
                /*Save pictures of the children*/
                int tmpChildCounter = 0;
                for(TreeNode<IAtomContainer> tmpNode : tmpMoleculeNode.getChildren()) {
                    tmpChildCounter++;
                    BufferedImage tmpChildImg = tmpGenerator.depict(tmpNode.getMolecule()).toImg();
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ChildNode" + tmpChildCounter + ".png").mkdirs();
                    File tmpChildFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ChildNode" + tmpChildCounter + ".png");
                    ImageIO.write(tmpChildImg, "png", tmpChildFile);
                }
                break;
            }
            tmpCounter++;
        }
    }

    /**
     * Tests the methods of the ScaffoldTree class with a V2000 or V3000 mol file as test molecule.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void scaffoldTreeTest() throws CDKException, CloneNotSupportedException, IOException {
        String tmpFileName = "Test11";
        //Load molecule from molfile
        IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
        //Generate a tree of molecules with iteratively removed terminal rings
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule, ElectronDonation.cdk()));
        /*Build ScaffoldTree*/
        int tmpCounter = 0;
        ScaffoldTree tmpScaffoldTree = new ScaffoldTree();
        while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
            tmpCounter++;
            //Add nodes to tree
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            tmpScaffoldTree.addNode(tmpMoleculeNode);
        }
        /*isMoleculeInTree and getMatrixNode checker*/
        System.out.println("---isMoleculeInTree and getMatrixNode check---");
        /*Should be false*/
        System.out.println("Is the original molecule in the tree: " + tmpScaffoldTree.isMoleculeInTree(tmpMolecule));
        assertEquals(false,tmpScaffoldTree.isMoleculeInTree(tmpMolecule));
        /*Should be true*/
        System.out.println("Is the molecule of node 14 in the tree: " + tmpScaffoldTree.isMoleculeInTree((IAtomContainer) tmpScaffoldTree.getMatrixNode(14).getMolecule()));
        assertEquals(true, tmpScaffoldTree.isMoleculeInTree((IAtomContainer) tmpScaffoldTree.getMatrixNode(14).getMolecule()));
        /*getTreeNode checker*/
        System.out.println("---getTreeNode check---");
        IAtomContainer tmpTestMolecule = (IAtomContainer) tmpScaffoldTree.getMatrixNode(10).getMolecule();
        System.out.println("Size of the test molecule: " + tmpTestMolecule.getAtomCount());
        IAtomContainer tmpResultMolecule = (IAtomContainer) tmpScaffoldTree.getTreeNode(tmpTestMolecule).getMolecule();
        System.out.println("Size of the result molecule: " + tmpResultMolecule.getAtomCount()); //Should be the same size
        assertEquals(tmpTestMolecule.getAtomCount(), tmpResultMolecule.getAtomCount());
        /*getAllNodes, getLevel and getMaxLevel checker*/
        System.out.println("---getAllNodes, getLevel and getMaxLevel check---");
        System.out.println("Total number of nodes: " + tmpScaffoldTree.getAllNodes().size());
        assertEquals(15, tmpScaffoldTree.getAllNodes().size());
        System.out.println("getMaxLevel:" + tmpScaffoldTree.getMaxLevel());
        assertEquals(3, tmpScaffoldTree.getMaxLevel());
        System.out.println("get number of nodes on level 2: " + tmpScaffoldTree.getAllNodesOnLevel(2).size());
        assertEquals(4, tmpScaffoldTree.getAllNodesOnLevel(2).size());
        /*Matrix check*/
        System.out.println("---getTreeAsMatrix check---");
        int tmpTreeCounter = 0;
        System.out.println("Is tree connected: " + tmpScaffoldTree.isTreeConnected());
        Integer[][] tmpMatrix = tmpScaffoldTree.getTreeAsMatrix();
        for(int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            for(int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) {
                System.out.print(tmpMatrix[tmpRow][tmpCol] + " ");
            }
            System.out.println(" " + tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow));
        }
        for(int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            if(tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow) < 10) {
                System.out.print(tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow) + " ");
            } else {
                System.out.print(tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow));
            }
        }
        /*Store all molecules in the order in which they appear in the matrix as images*/
        for(TreeNode tmpNode : tmpScaffoldTree.getMatrixNodes().values()) {
            IAtomContainer tmpNodeMolecule = (IAtomContainer) tmpNode.getMolecule();
            /*Save the picture*/
            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
            BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpNodeMolecule).toImg();
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/MatrixTest" + "/MatrixTest" + tmpTreeCounter  + ".png").mkdirs();
            File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/MatrixTest" + "/MatrixTest" + tmpTreeCounter  + ".png");
            ImageIO.write(tmpSecImgRemove, "png", tmpSecOutputRemove);
            tmpTreeCounter++;
        }
        /*RemoveNode, getRoot, hasOneSingleRootNode and isTreeConnected check*/
        System.out.println("---RemoveNode, getRoot, hasOneSingleRootNode and isTreeConnected check---");
        IAtomContainer tmpRemovedMolecule = (IAtomContainer) tmpScaffoldTree.getMatrixNode(14).getMolecule();
        tmpScaffoldTree.removeNode(tmpScaffoldTree.getMatrixNode(14));
        //Because the molecule occurs twice in the tree, it must still be present after the removal of a single node
        assertEquals(true, tmpScaffoldTree.isMoleculeInTree(tmpRemovedMolecule));
        System.out.println("Size of the tree after on Node is removed: "+ tmpScaffoldTree.getAllNodes().size());
        assertEquals(14, tmpScaffoldTree.getAllNodes().size());
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
        System.out.println("Root of the tree: " + tmpSmilesGenerator.create((IAtomContainer) tmpScaffoldTree.getRoot().getMolecule()));
        assertEquals(22,  ((IAtomContainer) tmpScaffoldTree.getRoot().getMolecule()).getAtomCount());
        System.out.println("Is tree connected: " + tmpScaffoldTree.isTreeConnected());
        assertEquals(true, tmpScaffoldTree.isTreeConnected());
        System.out.println("Has the tree one single root Node: " + tmpScaffoldTree.hasOneSingleRootNode());
        assertEquals(true, tmpScaffoldTree.hasOneSingleRootNode());
        tmpScaffoldTree.removeNode(tmpScaffoldTree.getRoot()); //Remove one node
        System.out.println("Is the tree connected after the root has been removed: " + tmpScaffoldTree.isTreeConnected());
        assertEquals(false, tmpScaffoldTree.isTreeConnected());
        System.out.println("Has the tree one single root node after the root has been removed: " + tmpScaffoldTree.hasOneSingleRootNode());
        assertEquals(false, tmpScaffoldTree.hasOneSingleRootNode());
    }

    /**
     * Creates a ScaffoldTree from a V2000 or V3000 mol file and displays it as a tree with GraphStream.
     * @throws InterruptedException if a problem with sleep occurs
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Ignore
    @Test
    public void graphStreamTest() throws InterruptedException, CDKException, IOException, CloneNotSupportedException {
        String tmpFileName = "Test11" ;
        //Load molecule from molfile
        IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
        //Generate a tree of molecules with iteratively removed terminal rings
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule, ElectronDonation.cdk()));
        int tmpCounter = 0;
        /*Create ScaffoldTree*/
        ScaffoldTree tmpScaffoldTree = new ScaffoldTree();
        while(tmpNodeIter.hasNext()) { //As long as there are still other molecules in the tree
            tmpCounter++;
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            tmpScaffoldTree.addNode(tmpMoleculeNode);
        }
        /*Remove some nodes*/
        //tmpScaffoldTree.removeNode(tmpScaffoldTree.getMatrixNode(24));
        //tmpScaffoldTree.removeNode(tmpScaffoldTree.getMatrixNode(22));
        //tmpScaffoldTree.removeNode(tmpScaffoldTree.getMatrixNode(23));
        /*Create a graph from the ScaffoldTree*/
        Graph tmpGraph = new SingleGraph("TestGraph");
        tmpGraph.setAttribute("ui.stylesheet", "node { size: 100px, 100px; }");
        System.setProperty("org.graphstream.ui", "swing");
        /*Add edges and nodes*/
        int tmpEdgeCount = 0;
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        Integer[][] tmpMatrix = tmpScaffoldTree.getTreeAsMatrix(); //Create the adjacency matrix
        for(int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) { //Create a node for each row
            /*Add the ScaffoldTree nodes to the graph*/
            tmpGraph.addNode(String.valueOf(tmpRow));
            Node tmpNode = tmpGraph.getNode(String.valueOf(tmpRow));
            tmpNode.setAttribute("Node", tmpScaffoldTree.getMatrixNode(tmpRow));
            /*Add a label to each node that corresponds to the position in the matrix*/
            tmpNode.setAttribute("ui.label", tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow));
            /*Add the images*/
            TreeNode tmpTreeNode =  tmpScaffoldTree.getMatrixNode(tmpScaffoldTree.getMatrixNodesNumbers().get(tmpRow));
            IAtomContainer tmpTreeNodeMolecule = (IAtomContainer) tmpTreeNode.getMolecule();
            BufferedImage tmpNodeImg = tmpGenerator.withSize(100,100).depict(tmpTreeNodeMolecule).toImg();
            //The images are stored temporarily, as I have not found a way to use them directly
            new File(System.getProperty("user.dir") + "//target/test-classes/GraphStream" + tmpRow + ".png").mkdirs();
            File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "//target/test-classes/GraphStream" + tmpRow + ".png");
            ImageIO.write(tmpNodeImg, "png", tmpSecOutputRemove);
            //set the images
            tmpNode.setAttribute("ui.style", "fill-mode: image-scaled-ratio-max;" + "fill-image: url('GraphStream" + tmpRow + ".png');");
            /*Add edges*/
            for(int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) { //Go through each column of the row
                if(tmpRow < tmpCol) { //Skip a diagonal half to get edges in one direction only.
                    continue;
                }
                if(tmpMatrix[tmpRow][tmpCol] == 1) { //Insert an edge if there is a 1 in it
                    tmpGraph.addEdge("Edge" + tmpEdgeCount, tmpRow, tmpCol);
                    tmpEdgeCount++;
                }
            }
        }
        /*Display graph*/
        System.setProperty("org.graphstream.ui", "swing");
        tmpGraph.display();
        TimeUnit.SECONDS.sleep(300);
    }
    //</editor-fold>

    //<editor-fold desc="Schuffenhauer Rules Tests">
    /**
     * Test of ScaffoldGenerator.applySchuffenhauerRules() with V2000 and V3000 mol files.
     * Loads the Test(Test1.mol-Test21.mol) molfiles from the Resources folder and creates the SchuffenhauerScaffolds with getSchuffenhauerScaffold().
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void applySchuffenhauerRulesTest() throws CloneNotSupportedException, CDKException, IOException {
        for (int tmpCount = 1; tmpCount < 23; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            /*Generate picture of molecule*/
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
            List<IAtomContainer> tmpSchuffenhauerFragments = scaffoldGenerator.applySchuffenhauerRules(tmpMolecule, true, ElectronDonation.cdk(), true);
            int tmpCounter = 0;
            for(IAtomContainer tmpFragment : tmpSchuffenhauerFragments) {
                tmpCounter++;
                /*Generate picture*/
                BufferedImage tmpImgFragment = tmpGenerator.depict(tmpFragment).toImg();
                /*Save the picture*/
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerRules/Fragment"+  tmpCounter + ".png").mkdirs();
                File tmpOutputFragment = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerRules/Fragment"+  tmpCounter + ".png");
                ImageIO.write(tmpImgFragment, "png", tmpOutputFragment);
            }
        }
    }

    /**
     * Test of ScaffoldGenerator.applySchuffenhauerRules() with SMILES.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void applySchuffenhauerRulesSMILESTest() throws CloneNotSupportedException, CDKException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C5CCCCCCCCCCCCCC1NC1CCCCCCCCCCCCCC(C3CC2CC2C4NC34)CC5");//Original
        //tmpMolecule = tmpParser.parseSmiles("C6CCC(CCC4CC(CC1CCCC1)C(C2CCCC2)C5C(C3CCCC3)CCC45)C6");
        //tmpMolecule = tmpParser.parseSmiles("OCC1OC2OC3C(O)C(O)C(OC3CO)OC4C(O)C(O)C(OC4CO)OC5C(O)C(O)C(OC5CO)OC6C(O)C(O)C(OC6CO)OC7C(O)C(O)C(OC7CO)OC8C(O)C(O)C(OC8CO)OC9C(O)C(O)C(OC9CO)OC%10C(O)C(O)C(OC%10CO)OC%11C(O)C(O)C(OC%11CO)OC%12C(O)C(O)C(OC%12CO)OC1C(O)C2O");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[F+]C1CC2");
        //tmpMolecule = tmpParser.parseSmiles("C2CCC1SC1C2");
        //tmpMolecule = tmpParser.parseSmiles("C2CCC(C1C[Br+]1)C2");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[I+]C1CCC2");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[I+]C1=CCC2");
        //tmpMolecule = tmpParser.parseSmiles("C2CC(C1C[Br+]1)CC2C3C[F+]3");
        //tmpMolecule = tmpParser.parseSmiles("[Cl+]2C3C1[I+]C1C4[Cl+]C234");
        //tmpMolecule = tmpParser.parseSmiles("C=1C=CC=2C(C1)=C3C=CC4=C5C=CC=CC5=C6C=CC2C3=C46");
        tmpMolecule = tmpParser.parseSmiles("c2ccc1[nH]ccc1c2");
        tmpMolecule = tmpParser.parseSmiles("O=C4C(=O)C3CC2C1CNNC1NC2C3C4=O");
        tmpMolecule= tmpParser.parseSmiles("C1NOCC2CSNCC12");
        tmpMolecule = tmpParser.parseSmiles("C2CC1NOCC1C3NSCC23");
        tmpMolecule = tmpParser.parseSmiles("C2NNNC3CC1SNNCC1CC23");
        tmpMolecule = tmpParser.parseSmiles("O=C2CC(=[Br+])C1CCCCCC1CC2=P");
        //tmpMolecule = tmpParser.parseSmiles("CCN(C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C#N)C(=O)C");
        tmpMolecule = tmpParser.parseSmiles("N1=C2C=CC=CC2=NC3=C1C(=NN3C=4C=CC=CC4)N");
        tmpMolecule = tmpParser.parseSmiles("N1=CC=2C(=NN(C=3C=CC=CC3)S2N1C)C=4C=CC=CC4");
        tmpMolecule = tmpParser.parseSmiles("c1ccc5c(c1)c2ccccc2C6c3ccccc3c4ccccc4C56");
        tmpMolecule = tmpParser.parseSmiles("N=1C=CN(C1)C2=NC=CO2"); //CNP0001107 example rule 7 [N] problem
        tmpMolecule = tmpParser.parseSmiles("O=CN1C=CC2=CC=NN21"); //CNP0000568 rule 7 example
        tmpMolecule = tmpParser.parseSmiles("N=1C=2C=CC=CC2N=C3C1C=4C=CC=C(C34)C"); //CNP0058597 rule 7 example
        //tmpMolecule = tmpParser.parseSmiles("c1c[nH]cn1");
        //tmpMolecule = tmpParser.parseSmiles("C5CCN(CCC3CN(CCC1CCNC1)CN(CCN2CCCC2)C3CCC4CCCN4)C5"); //Rule 12
        //tmpMolecule = tmpParser.parseSmiles("C5CCN(C3CN(C1CCCN1)C(C2CCNC2)N3N4CCCC4)C5"); //Rule 12
        //tmpMolecule = tmpParser.parseSmiles("[I+]=C(NC1CCCN1)C4C(CCN2CCCC2)CN(CCC3CCNC3)CN4CCN5CCCC5"); //Rule 12
        tmpMolecule = tmpParser.parseSmiles("O=C(O)C1CCC(=CO)C(CCCNC2=CC=C(C=[NH+]2)C=3C=CC=C(C3)CC4=C5C=CC=CC5=CC6C4=CC78C=CCC9(C)C(O)CCC6(C%10=C(C7)C%11%12CC(O)C%131C=CC%14C%15%16CCC%14%12C(C=CC%16CC(CC=%17C=CC=CC%17CCCCC)C%15)CC(C%10)C%13%11C)C89)C%18CCCCC%18");
        tmpMolecule = tmpParser.parseSmiles("[H]OC(=O)C1([H])C23C([H])=C([H])C4([H])C56C7(C8=C(C9%10C%11([H])C(C(=C%12C(C([H])=C([H])C([H])=C%12[H])=C%11[H])C([H])([H])C%13=C([H])C(C=%14C([H])=C([H])C(N([H])C([H])([H])C([H])([H])C([H])([H])C([H])(C(=C([H])O[H])C([H])([H])C1([H])[H])C%15([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C%15([H])[H])=[N+]([H])C%14[H])=C([H])C([H])=C%13[H])=C([H])C%16(C([H])=C([H])C([H])([H])C(C([H])(O[H])C([H])([H])C9([H])[H])(C%16%10[H])C([H])([H])[H])C8([H])[H])C([H])([H])C([H])(C72C([H])([H])[H])C([H])([H])C6([H])C([H])=C([H])C%17([H])C4(C([H])([H])C([H])(C([H])([H])C=%18C(=C([H])C([H])=C([H])C%18[H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H])C%17([H])[H])C([H])([H])C5([H])[H])C([H])([H])C3([H])O[H]");
        tmpMolecule = tmpParser.parseSmiles("S=P12N=P3(OC=4C=CC(C=NN(C)P(=S)(C=5C=CC=CC5)N(N=CC=6C=CC(OP(=NP(=S)(OC=7C=CC(C=NN(C)P(=S)(C=8C=CC=CC8)N(N=CC=9C=CC(O1)=CC9)C)=CC7)OC=%10C=CC(C=NN(C)P(=S)(C=%11C=CC=CC%11)N(N=CC=%12C=CC(O2)=CC%12)C)=CC%10)(OC=%13C=CC(C=NN(C)P(=S)(C=%14C=CC=CC%14)N(N=CC=%15C=CC(O3)=CC%15)C)=CC%13)C=%16C=CC=CC%16)=CC6)C)=CC4)C=%17C=CC=CC%17");
        //tmpMolecule = tmpParser.parseSmiles("OCC1OC2OC3C(O)C(O)C(OC3CO)OC4C(O)C(O)C(OC4CO)OC5C(O)C(O)C(OC5CO)OC6C(O)C(O)C(OC6CO)OC7C(O)C(O)C(OC7CO)OC8C(O)C(O)C(OC8CO)OC9C(O)C(O)C(OC9CO)OC%10C(O)C(O)C(OC%10CO)OC%11C(O)C(O)C(OC%11CO)OC1C(O)C2O");
        //tmpMolecule = tmpParser.parseSmiles("O=C(OCC1OC2OC3C(O)C(O)C(OC3CO)OC4C(O)C(O)C(OC4COC(=O)C)OC5C(O)C(O)C(OC5COC(=O)C)OC6C(O)C(O)C(OC6COC(=O)C)OC7C(O)C(O)C(OC7CO)OC1C(O)C2O)C");
        tmpMolecule = tmpParser.parseSmiles("C1C2CC3CC1CC(C2)C3");
        //tmpMolecule = tmpParser.parseSmiles("OCC1OC2OC3C(O)C(O)C(OC3CO)OC4C(O)C(O)C(OC4CO)OC5C(O)C(O)C(OC5CO)OC6C(O)C(O)C(OC6CO)OC7C(O)C(O)C(OC7CO)OC8C(O)C(O)C(OC8CO)OC9C(O)C(O)C(OC9CO)OC1C(O)C2O");
        /*Generate picture of molecule*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        List<IAtomContainer> tmpSchuffenhauerFragments = scaffoldGenerator.applySchuffenhauerRules(tmpMolecule, true, ElectronDonation.cdk(), false);
        int tmpCounter = 0;
        for(IAtomContainer tmpFragment : tmpSchuffenhauerFragments) {
            tmpCounter++;
            /*Generate picture*/
            BufferedImage tmpImgFragment = tmpGenerator.depict(tmpFragment).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/RulesSMILESTest/Fragment" + tmpCounter + ".png").mkdirs();
            File tmpOutputFragment = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/RulesSMILESTest/Fragment" + tmpCounter + ".png");
            ImageIO.write(tmpImgFragment, "png", tmpOutputFragment);
        }
    }
    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold() with SMILES.
     * Loads Scheme 1 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES and generates SchuffenhauerScaffold.
     * Flucloxacillin is generated from the SMILES and all terminal side chains are removed. Rings, linkers and double bonds on these structures are obtained.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme1Test() throws IOException, CDKException, CloneNotSupportedException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");
        //tmpMolecule = tmpParser.parseSmiles("C=C1CCCCNC2=CC=C(C=[NH+]2)C=3C=CC=C(C3)CC4=CC=CC5C4=CC6CC7=C(CC8CC9C=CCCC%10C=CC(CCC1)C8C7C%109)C%115CCCCC6%11");
        /*Generate picture of the Original molecule*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture of the original*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme1/Original.png").mkdirs();
        File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme1/Original.png");
        ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerSMILES = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme1/Schuffenhauer.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme1/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpSchuffenhauerSMILES), "O=C(NC1C(=O)N2CCSC21)C3=CON=C3C=4C=CC=CC4");
    }

    /**
     * Loads Scheme 2b from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Adamantane is generated from the SMILES and it is checked whether rings can be removed. This should not be the case.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme2bTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1C2CC3CC1CC(C2)C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme2b/Schuffenhauer.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme2b/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold with removed ring*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        //Get rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
        int tmpCounter = 0;
        for(IAtomContainer tmpRing : tmpRings) {
            boolean tmpIsRingRemovable = scaffoldGenerator.isRingRemovable(tmpSchuffenhauerScaffold, tmpRings, tmpRing);
            /*Remove rings*/
            IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
            /*Generate picture*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme2b/RingRemovable" + tmpIsRingRemovable + tmpCounter + ".png").mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme2b/RingRemovable" + tmpIsRingRemovable + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
            /*Check boolean*/
            assertEquals(tmpIsRingRemovable, false);
        }
    }

    /**
     * Loads Scheme 3a from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * A molecule consisting of two rings is generated from a SMILES.
     * One of these rings is aromatic and has to be removed.
     * At the point where this aromatic ring was bound to the other ring, a double bond should now be formed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme3aTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("c2ccc1CNCCc1c2");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1CCNCC1=CCC2");

        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3a/Schuffenhauer.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3a/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold with removed ring*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        //Get rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
        //Remove Ring
        IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRings.get(1));
        /*Generate picture*/
        BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3a/RemovedRing.png").mkdirs();
        File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3a/RemovedRing.png");
        ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRemovedSchuff), "C1=CCCNC1");
    }

    /**
     * Loads Scheme 3b from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * A molecule consisting of three rings is generated from a SMILES. One of these rings is aromatic.
     * It is tested whether this aromatic ring can be removed. This should not be the case.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme3bTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("c1cc2CCCc3c[nH]c(c1)c23");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3b/Schuffenhauer.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3b/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold with removed ring*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        //Get rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
        int tmpCounter = 0;
        for(IAtomContainer tmpRing : tmpRings) {
            boolean tmpIsRemovable = scaffoldGenerator.isRingRemovable(tmpRing, tmpRings, tmpSchuffenhauerScaffold);
            /*Remove rings*/
            IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
            /*Generate picture*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3b/IsRingRemovable_" + tmpIsRemovable + tmpCounter + ".png").mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme3b/IsRingRemovable_" + tmpIsRemovable + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Check booleans*/
        assertEquals(scaffoldGenerator.isRingRemovable(tmpRings.get(0), tmpRings, tmpSchuffenhauerScaffold), false);
        assertEquals(scaffoldGenerator.isRingRemovable(tmpRings.get(1), tmpRings, tmpSchuffenhauerScaffold), false);
        assertEquals(scaffoldGenerator.isRingRemovable(tmpRings.get(2), tmpRings, tmpSchuffenhauerScaffold), true);
    }

    /**
     * Loads Scheme 4 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Epothilone A is generated from a SMILES and the ring consisting of 3 atoms is removed.
     * The removal of this hetero ring should result in a double bond at the removed position.
     * In addition, some similar example molecules are also commented out.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme4Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1CCCC2C(O2)CC(OC(=O)CC(C(C(=O)C(C1O)C)(C)C)O)C(=CC3=CSC(=N3)C)C");//Original
        //tmpMolecule = tmpParser.parseSmiles("C3CCC1[F+]C1C2CC2C3");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[F+]C1CC2");
        //tmpMolecule = tmpParser.parseSmiles("C2=c1[nH]c1=CCC2");
        //tmpMolecule = tmpParser.parseSmiles("C2CCC1SC1C2");
        //tmpMolecule = tmpParser.parseSmiles("C2CCC(C1C[Br+]1)C2");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[I+]C1CCC2");
        //tmpMolecule = tmpParser.parseSmiles("C2=C1[I+]C1=CCC2");
        //tmpMolecule = tmpParser.parseSmiles("C2CC(C1C[Br+]1)CC2C3C[F+]3");
        //tmpMolecule = tmpParser.parseSmiles("[Cl+]2C3C1[I+]C1C4[Cl+]C234");
        //tmpMolecule = tmpParser.parseSmiles("O=C1C2CCCCC12");
        //tmpMolecule = tmpParser.parseSmiles("C2CCC(C1NO1)C2");
        //tmpMolecule = tmpParser.parseSmiles("C2CC(C1NO1)C3NC23");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        //Get rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
        int tmpCounter = 0;
        for(IAtomContainer tmpRing : tmpRings) {
            /*Remove rings*/
            IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
            /*Generate picture*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/RemovedRing" + tmpCounter + ".png").mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/RemovedRing" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Generate picture of the SchuffenhauerRuleOne*/
        List<IAtomContainer> tmpRuleOne = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRuleOne = tmpGenerator.depict(tmpRuleOne.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/RuleOne.png").mkdirs();
        File tmpOutputRuleOne = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme4/RuleOne.png");
        ImageIO.write(tmpImgRuleOne, "png" ,tmpOutputRuleOne);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRuleOne.get(1)), "O=C1OC(C=CC=2N=CSC2)CC=CCCCCCCC(=O)CCC1");
    }

    /**
     * Loads Scheme5 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Seglitide is generated from a SMILES and the two single rings connected via linker are removed.
     * Then, according to the second rule, the aromatic 6 ring is removed to obtain the macroring.
     * In addition, some similar example molecules are also commented out.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme5Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N1C)CC2=CC=CC=C2)C(C)C)CCCCN)CC3=CNC4=CC=CC=C43)CC5=CC=C(C=C5)O");//Original
        //tmpMolecule = tmpParser.parseSmiles("Sc2cccc3[nH]c(CC1CCCCCCCCCCCCCCCCC1)cc23");
        //tmpMolecule = tmpParser.parseSmiles("Sc3cc1CCCc1c4[nH]c(CC2CCCCCCCCCCCCCCCCC2)cc34");
        //tmpMolecule = tmpParser.parseSmiles("Sc3cccc4[nH]c(CC1CCCCCCCCCCCC2(CCCCC1)CC2)cc34");
        //tmpMolecule = tmpParser.parseSmiles("O=C3CCCC2NC(CC1CCCCCCCCCCCCCCCCC1)=CC2C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerScaffold with removed ring*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        //Get rings
        List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true).subList(1,3);
        tmpSchuffenhauerScaffold = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRings.get(0));
        tmpSchuffenhauerScaffold = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRings.get(1));
        /*Generate picture*/
        BufferedImage tmpImgRemove = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/ModifiedMolecule" + ".png").mkdirs();
        File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/ModifiedMolecule" + ".png");
        ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/RuleTwo.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme5/RuleTwo.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "O=C1NCC(=O)NCC(=O)NC(C(=O)NCC(=O)NCC(=O)NC1)CC=2C=CNC2");
    }

    /**
     * Loads Scheme 6  from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Flucloxacillin is generated from a SMILES and the ring consisting of 6 atoms is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme6Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");//Original
        //tmpMolecule = tmpParser.parseSmiles("O=C(CC(=S)C(=O)C2CCCC(CCCCCCC1CCCCC1)C2)C(=[I+])C3CCCCC3");//Linker with double bond
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme6/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme6/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme6/RuleThree.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme6/RuleThree.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "O=C(NC1C(=O)N2CCSC21)C=3C=NOC3");
    }

    /**
     * Loads Scheme 7 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Pentazocine is generated from a SMILES and the aromatic ring consisting of 6 atoms is removed.
     * A double bond is inserted at the point where it was removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme7Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1C2CC3=C(C1(CCN2CC=C(C)C)C)C=C(C=C3)O");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme7/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme7/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme7/RuleFour.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme7/RuleFour.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "C1=CC2CCNC(C1)C2");
    }

    /**
     * Loads Scheme 8 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Sophocarpin is generated from a SMILES and the ring that only has an overlapping bond is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme8Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1CC2CN3C(CC=CC3=O)C4C2N(C1)CCC4");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme8/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme8/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme8/RuleFour.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme8/RuleFour.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "N1CC2CCCN3CCCC(C1)C32");
    }

    /**
     * Loads Scheme 9 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Rhynchophylline is generated from a SMILES and the aromatic ring is removed.
     * The 6 ring is now removed from this fragment.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme9Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCC1CN2CCC3(C2CC1C(=COC)C(=O)OC)C4=CC=CC=C4NC3=O");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/RuleFour.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme9/RuleFour.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "O=C1NC=CC12CNCC2");
    }

    /**
     * Loads Scheme 10 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Cafestol is generated from a SMILES and the aromatic ring is removed.
     * The 6 ring is now removed from this fragment.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme10Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC12CCC3=C(C1CCC45C2CCC(C4)C(C5)(CO)O)C=CO3");
        //tmpMolecule = tmpParser.parseSmiles("C1CCC23CCC(CCC2C1)C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(3)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/RuleFive.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme10/RuleFive.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(3)), "C1CC2CCC(C1)C2");
    }
    /**
     * Loads Scheme 11a from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Flucloxacillin is generated from a SMILES and the rings connected via linkers are removed.
     * Then, according to the sixth rule, the ring of size 5 is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme11aTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");
        //tmpMolecule = tmpParser.parseSmiles("C1CCC23CCC(CCC2C1)C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(3)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/RuleSix.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11a/RuleSix.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(3)), "O=C1NCC1");
    }

    /**
     * Loads Scheme 11b from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Epinastine is generated from a SMILES and the aromatic rings are removed.
     * Then, according to the sixth rule, the ring of size 5 is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme11bTest() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1C2C3=CC=CC=C3CC4=CC=CC=C4N2C(=N1)N");
        //tmpMolecule = tmpParser.parseSmiles("C1CCC23CCC(CCC2C1)C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(3)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/RuleSix.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme11b/RuleSix.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(3)), "C1=CCC=CCN1");
    }

    /**
     * Loads Scheme 12 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Zaleplon is generated from a SMILES and the aromatic 6 ring is removed.
     * Then, according to the seventh rule, the ring of size 6 is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme12Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCN(C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C#N)C(=O)C"); //Original
        //tmpMolecule = tmpParser.parseSmiles("c2cnc1ccnn1c2");
        //tmpMolecule = tmpParser.parseSmiles("c1(ccn2)n2cccn1");
        //tmpMolecule = tmpParser.parseSmiles("O=C(OC)N1C=NC=2C(=NC=NC21)N"); //better example for CDK
        //tmpMolecule = tmpParser.parseSmiles("[S-]C1=NC=NC2=C1C=NN2C=3C=CC=CC3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpMolecule, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/RuleSeven.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme12/RuleSeven.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        //assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "C1=CCC=CCN1");
    }

    /**
     * Test that illustrates aromaticity detection for pyrimidine, pyrazole, and their combination in
     * pyrazolo[1,5-a]pyrimidine, which is used in Scheme 12 for illustrating rule 7. Different cycle finder algorithms
     * and electron donation models are combined and the results printed as SMILES strings with aromaticity encoded.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void pyrazoloPyrimidineElectronDonationAndCycleFinderTest() throws Exception {
        String tmpPyrimidineSMILES = "C1=CC=NC=N1";
        String tmpPyrazoleSMILES = "C1=CC=NN1";
        String tmpPyrazoloPyrimidineSMILES = "N1(C=CC=N2)C2=CC=N1";
        HashMap<String, String> tmpMoleculesMap = new HashMap<>(5, 1);
        tmpMoleculesMap.put(tmpPyrimidineSMILES, "pyrimidine");
        tmpMoleculesMap.put(tmpPyrazoleSMILES, "pyrazole");
        tmpMoleculesMap.put(tmpPyrazoloPyrimidineSMILES, "Pyrazolo[1,5-a]pyrimidine");
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols);
        HashMap<ElectronDonation, String> tmpElectronDonationModelsMap = new HashMap(6, 1);
        tmpElectronDonationModelsMap.put(ElectronDonation.cdk(), "CDK");
        tmpElectronDonationModelsMap.put(ElectronDonation.cdkAllowingExocyclic(), "CDK exocyclic");
        tmpElectronDonationModelsMap.put(ElectronDonation.daylight(), "Daylight");
        tmpElectronDonationModelsMap.put(ElectronDonation.piBonds(), "pi bonds");
        CycleFinder[] tmpCycleFinders = {Cycles.all(), Cycles.cdkAromaticSet(), Cycles.mcb(), Cycles.relevant()};
        for (String tmpSMILES : tmpMoleculesMap.keySet()) {
            System.out.println("\n" + tmpMoleculesMap.get(tmpSMILES));
            for (CycleFinder tmpCF : tmpCycleFinders) {
                System.out.println("\n\t" + tmpCF);
                for (ElectronDonation tmpEDmodel : tmpElectronDonationModelsMap.keySet()) {
                    Aromaticity tmpAromaticityModel = new Aromaticity(tmpEDmodel, tmpCF);
                    IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles(tmpSMILES);
                    AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
                    tmpAromaticityModel.apply(tmpMolecule);
                    int tmpNumberOfRings = tmpCF.find(tmpMolecule).numberOfCycles();
                    System.out.println("\t\t" + tmpElectronDonationModelsMap.get(tmpEDmodel) + " " + tmpNumberOfRings + " " + tmpSmiGen.create(tmpMolecule));
                }
            }
        }
    }


    /**
     * Loads Scheme 13 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * A double ring system is generated from a SMILES.
     * According to the eighth rule the ring with the least heterocycle is removed
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme13Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("c2ccc1[nH]ccc1c2");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme13/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme13/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme13/RuleEight.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme13/RuleEight.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "C=1C=CNC1");
    }

    /**
     * Loads Scheme 14 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Ticlopidine is generated from a SMILES.
     * According to the ninth rule the ring with the S is removed first
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme14Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1CN(CC2=C1SC=C2)CC3=CC=CC=C3Cl.Cl");
        //tmpMolecule = tmpParser.parseSmiles("C2CC1NCOCC1NS2");
        //tmpMolecule = tmpParser.parseSmiles("N=C2CC1NCOCC1NS2");
        //tmpMolecule = tmpParser.parseSmiles("P=S2C1NPPPC1S(=P)S(=P)S2=P");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/RuleNine.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme14/RuleNine.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "C1=CCCNC1");
    }
    /**
     * Loads a molecule to check the rule ten from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * A molecule with two 7 rings and one 8 ring is generated from a SMILES.
     * According to the tenth rule the 7 rings are removed first
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRule10Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("O=C3CC2C(CCCC1CCCCCCC12)C(=O)C(=O)C3=O");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/RuleTen.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule10/RuleTen.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "C1CCCCCCC1");
    }

    /**
     * Loads Scheme 15 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Sertraline is generated from a SMILES and the linker bonded 6 ring is removed.
     * According to the eleventh rule the aromatic ring is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme15Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/RuleEleven.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme15/RuleEleven.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "C1=CCCCC1");
    }

    /**
     * Loads Scheme 16 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Deferasirox is generated from a SMILES.
     * According to the twelfth rule the aromatic ring bond to the N is removed.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme16Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1=CC=C(C(=C1)C2=NN(C(=N2)C3=CC=CC=C3O)C4=CC=C(C=C4)C(=O)O)O");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme16/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme16/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme16/RuleTwelve.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme16/RuleTwelve.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(1)), "N=1[N]C(=NC1C=2C=CC=CC2)C=3C=CC=CC3");
    }

    /**
     * Loads Scheme 17 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Ormeloxifene is generated from a SMILES and the linker bonded 6 rings are removed.
     * The generated scaffold "Thirteen" does not correspond to the illustration in the paper.
     * This is due to the fact that unique SMILES are generated for rule 13, although canical SMILES are used in the paper.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme17Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1(C(C(C2=C(O1)C=C(C=C2)OC)C3=CC=C(C=C3)OCCN4CCCC4)C5=CC=CC=C5)C");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpMod = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMod.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Modified.png").mkdirs();
        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Modified.png");
        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Thirteen.png").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/Thirteen.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate picture of the last step*/
        List<IAtomContainer> tmpLast = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgLast = tmpGenerator.depict(tmpLast.get(3)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/LastStep.png").mkdirs();
        File tmpOutputLast = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme17/LastStep.png");
        ImageIO.write(tmpImgLast, "png" ,tmpOutputLast);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpRule.get(2)), "O1C=2C=CC=CC2C(C=3C=CC=CC3)CC1");
    }

    /**
     * Loads Scheme 18 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Diazepam, Bromazepam, Zolazepam and Clotiazepam are generated from SMILES.
     * The Schuffenhauer rules are then applied to these molecules.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme18Test() throws CDKException, CloneNotSupportedException, IOException {
        /*-----Diazepam-----*/
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMoleculeDiazepam = tmpParser.parseSmiles("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffoldDiazepam = scaffoldGenerator.getSchuffenhauerScaffold(tmpMoleculeDiazepam, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILESDiazepam = tmpGenerator.depict(tmpSchuffenhauerScaffoldDiazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamOriginal.png").mkdirs();
        File tmpOutputSMILESDiazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamOriginal.png");
        ImageIO.write(tmpImgSMILESDiazepam, "png" ,tmpOutputSMILESDiazepam);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpStep1MolDiazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldDiazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep1Diazepam = tmpGenerator.depict(tmpStep1MolDiazepam.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamStep1.png").mkdirs();
        File tmpOutputStep1Diazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamStep1.png");
        ImageIO.write(tmpImgStep1Diazepam, "png" ,tmpOutputStep1Diazepam);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpStep2MolDiazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldDiazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep2Diazepam = tmpGenerator.depict(tmpStep2MolDiazepam.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamStep2.png").mkdirs();
        File tmpOutputStep2Diazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/DiazepamStep2.png");
        ImageIO.write(tmpImgStep2Diazepam, "png" ,tmpOutputStep2Diazepam);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolDiazepam.get(1)), "O=C1NC=2C=CC=CC2C=NC1");
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolDiazepam.get(2)), "O=C1NC=CC=NC1");

        /*-----Bromazepam-----*/
        //SMILES to IAtomContainer
        IAtomContainer tmpMoleculeBromazepam = tmpParser.parseSmiles("C1C(=O)NC2=C(C=C(C=C2)Br)C(=N1)C3=CC=CC=N3");
        /*Generate picture of the SchuffenhauerScaffold*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffoldBromazepam = scaffoldGenerator.getSchuffenhauerScaffold(tmpMoleculeBromazepam, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILESBromazepam = tmpGenerator.depict(tmpSchuffenhauerScaffoldBromazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamOriginal.png").mkdirs();
        File tmpOutputSMILESBromazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamOriginal.png");
        ImageIO.write(tmpImgSMILESBromazepam, "png" ,tmpOutputSMILESBromazepam);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpStep1MolBromazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldBromazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep1Bromazepam = tmpGenerator.depict(tmpStep1MolBromazepam.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamStep1.png").mkdirs();
        File tmpOutputStep1Bromazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamStep1.png");
        ImageIO.write(tmpImgStep1Bromazepam, "png" ,tmpOutputStep1Bromazepam);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpStep2MolBromazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldBromazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep2Bromazepam = tmpGenerator.depict(tmpStep2MolBromazepam.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamStep2.png").mkdirs();
        File tmpOutputStep2Bromazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/BromazepamStep2.png");
        ImageIO.write(tmpImgStep2Bromazepam, "png" ,tmpOutputStep2Bromazepam);
        /*Generate and check SMILES*/
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolBromazepam.get(1)), "O=C1NC=2C=CC=CC2C=NC1");
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolBromazepam.get(2)), "O=C1NC=CC=NC1");

        /*-----Zolazepam-----*/
        //SMILES to IAtomContainer
        IAtomContainer tmpMoleculeZolazepam = tmpParser.parseSmiles("CC1=NN(C2=C1C(=NCC(=O)N2C)C3=CC=CC=C3F)C");
        /*Generate picture of the SchuffenhauerScaffold*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffoldZolazepam = scaffoldGenerator.getSchuffenhauerScaffold(tmpMoleculeZolazepam, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILESZolazepam = tmpGenerator.depict(tmpSchuffenhauerScaffoldZolazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamOriginal.png").mkdirs();
        File tmpOutputSMILESZolazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamOriginal.png");
        ImageIO.write(tmpImgSMILESZolazepam, "png" ,tmpOutputSMILESZolazepam);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpStep1MolZolazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldZolazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep1Zolazepam = tmpGenerator.depict(tmpStep1MolZolazepam.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamStep1.png").mkdirs();
        File tmpOutputStep1Zolazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamStep1.png");
        ImageIO.write(tmpImgStep1Zolazepam, "png" ,tmpOutputStep1Zolazepam);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpStep2MolZolazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldZolazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep2Zolazepam = tmpGenerator.depict(tmpStep2MolZolazepam.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamStep2.png").mkdirs();
        File tmpOutputStep2Zolazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ZolazepamStep2.png");
        ImageIO.write(tmpImgStep2Zolazepam, "png" ,tmpOutputStep2Zolazepam);
        /*Generate and check SMILES*/
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolZolazepam.get(1)), "O=C1NC=2NN=CC2C=NC1");
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolZolazepam.get(2)), "O=C1NC=CC=NC1");

        /*-----Clotiazepam-----*/
        //SMILES to IAtomContainer
        IAtomContainer tmpMoleculeClotiazepam = tmpParser.parseSmiles("CCC1=CC2=C(S1)N(C(=O)CN=C2C3=CC=CC=C3Cl)C");
        /*Generate picture of the SchuffenhauerScaffold*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffoldClotiazepam = scaffoldGenerator.getSchuffenhauerScaffold(tmpMoleculeClotiazepam, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILESClotiazepam = tmpGenerator.depict(tmpSchuffenhauerScaffoldClotiazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamOriginal.png").mkdirs();
        File tmpOutputSMILESClotiazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamOriginal.png");
        ImageIO.write(tmpImgSMILESClotiazepam, "png" ,tmpOutputSMILESClotiazepam);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpStep1MolClotiazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldClotiazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep1Clotiazepam = tmpGenerator.depict(tmpStep1MolClotiazepam.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamStep1.png").mkdirs();
        File tmpOutputStep1Clotiazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamStep1.png");
        ImageIO.write(tmpImgStep1Clotiazepam, "png" ,tmpOutputStep1Clotiazepam);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpStep2MolClotiazepam = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldClotiazepam, true, ElectronDonation.cdk(), true);
        BufferedImage tmpImgStep2Clotiazepam = tmpGenerator.depict(tmpStep2MolClotiazepam.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamStep2.png").mkdirs();
        File tmpOutputStep2Clotiazepam = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme18/ClotiazepamStep2.png");
        ImageIO.write(tmpImgStep2Clotiazepam, "png" ,tmpOutputStep2Clotiazepam);
        /*Generate and check SMILES*/
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolClotiazepam.get(1)), "O=C1NC=2SC=CC2C=NC1");
        assertEquals(tmpSmilesGenerator.create(tmpStep2MolClotiazepam.get(2)), "O=C1NC=CC=NC1");
    }

    /**
     * Loads Scheme 19 from the "The Scaffold Tree" Paper by Schuffenhauer et al as SMILES.
     * Baccatin III is generated from a SMILES and decomposed according to the Schuffenahauer rules.
     * -Step 1: The aromatic 6 ring is removed according to rule 3
     * -Step 2: The 4 ring is removed according to rule 4
     * -Step 3: The 6 ring without DB is removed according to rule 4
     * -Step 4: The 6 ring with DB is removed according to rule 6
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getScheme19Test() throws CDKException, CloneNotSupportedException, IOException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1O)O)OC(=O)C5=CC=CC=C5)(CO4)OC(=O)C)O)C)OC(=O)C");
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Original.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Original.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate the scaffolds*/
        List<IAtomContainer> tmpScaffolds = scaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffold, true, ElectronDonation.cdk(), true);
        /*Generate picture of the modified molecule*/
        BufferedImage tmpImgStep1 = tmpGenerator.depict(tmpScaffolds.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step1.png").mkdirs();
        File tmpOutputStep1 = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step1.png");
        ImageIO.write(tmpImgStep1, "png" ,tmpOutputStep1);
        /*Generate picture of the modified molecule*/
        BufferedImage tmpImgStep2 = tmpGenerator.depict(tmpScaffolds.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step2.png").mkdirs();
        File tmpOutputStep2 = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step2.png");
        ImageIO.write(tmpImgStep2, "png" ,tmpOutputStep2);
        /*Generate picture of the modified molecule*/
        BufferedImage tmpImgStep3 = tmpGenerator.depict(tmpScaffolds.get(3)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step3.png").mkdirs();
        File tmpOutputStep3 = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step3.png");
        ImageIO.write(tmpImgStep3, "png" ,tmpOutputStep3);
        /*Generate picture of the modified molecule*/
        BufferedImage tmpImgStep4 = tmpGenerator.depict(tmpScaffolds.get(4)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step4.png").mkdirs();
        File tmpOutputStep4 = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Scheme19/Step4.png");
        ImageIO.write(tmpImgStep4, "png" ,tmpOutputStep4);
        /*Generate and check SMILES*/
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        assertEquals(tmpSmilesGenerator.create(tmpScaffolds.get(1)), "O=C1CC2=CCCC(C2)CC3C1CCC4OCC43");
        assertEquals(tmpSmilesGenerator.create(tmpScaffolds.get(2)), "O=C1CC2=CCCC(C2)CC3CCCCC13");
        assertEquals(tmpSmilesGenerator.create(tmpScaffolds.get(3)), "O=C1CC2=CCCC(C2)CCC1");
        assertEquals(tmpSmilesGenerator.create(tmpScaffolds.get(4)), "O=C1CCCCCCC1");
    }
    //</editor-fold>

    //<editor-fold desc="Speed Tests">
    /**
     * Speed test for the getSchuffenhauerScaffold() Method with over 400000 molecules from the COCONUT DB.
     * @throws CDKException problem with the AtomContainerManipulator
     * @throws IOException if COCONUT DB not found
     */
    @Ignore
    @Test
    public void calculateSchuffenhauerSpeedTest() throws IOException, CDKException {
        this.ScaffoldGeneratorSpeedTest(true, false, false, false, true, 4242);
    }

    /**
     * Speed test for the getRings() Method with over 400000 molecules from the COCONUT DB.
     * getSchuffenhauerScaffold() must also be executed for all molecules.
     * @throws CDKException problem with the AtomContainerManipulator
     * @throws IOException if COCONUT DB not found
     */
    @Ignore
    @Test
    public void calculateRingsSpeedTest() throws IOException, CDKException {
        this.ScaffoldGeneratorSpeedTest(false, true, false, false, true, 4242);
    }

    /**
     * Speed test for the removeRing() Method with over 400000 molecules from the COCONUT DB.
     * getSchuffenhauerScaffold() and getRings() must also be executed for all molecules.
     * Skips all molecules with more than 1000 rings.
     * In this case, these are the molecules with the COCONUT IDs: CNP0022608, CNP0029543, CNP0065312 and CNP0103752.
     * Runtime
     * @throws CDKException problem with the AtomContainerManipulator
     * @throws IOException if COCONUT DB not found
     */
    @Ignore
    @Test
    public void calculateRemoveRingsSpeedTest() throws IOException, CDKException {
        this.ScaffoldGeneratorSpeedTest(false, false, true, false, true, 4242);
    }

    /**
     * Speed test for the removeRing() Method with over 400000 molecules from the COCONUT DB.
     * getSchuffenhauerScaffold() and getRings() must also be executed for all molecules.
     * Skips all molecules with more than 1000 rings.
     * In this case, these are the molecules with the COCONUT IDs: CNP0022608, CNP0029543, CNP0065312 and CNP0103752.
     * Runtime
     * @throws CDKException problem with the AtomContainerManipulator
     * @throws IOException if COCONUT DB not found
     */
    @Ignore
    @Test
    public void calculateApplySchuffenhauerRulesSpeedTest() throws CDKException, IOException {
        this.ScaffoldGeneratorSpeedTest(false, false, false, true, true, 4242);
    }

    /**
     * Speed test for the getSchuffenhauerScaffold(), getRing() and removeRing() Method with over 400000 molecules from the COCONUT DB.
     * Which methods are tested can be set via the booleans.
     * To perform the test download the COCONUT DB(https://coconut.naturalproducts.net/download) and add the COCONUT_DB.sdf file to src\test\resources
     * @param anIsSchuffenhauerScaffoldCalculated Generate SchuffenhauerScaffolds
     * @param anIsRingCalculated Calculate Rings and Schuffenhauer scaffolds.
     * @param anIsRemoveRingCalculated The molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated. The Schuffenhauer scaffolds and the Rings are also calculated for this.
     * @param anIsApplySchuffenhauerCalculated Creates all molecule fragments generated by the schuffenhauer rules
     * @param anIsPictureCreated Show control pictures from one molecule.
     * @param aPictureNumber Number of the molecule from which control images are to be taken(from 0 to 406000)
     * @throws FileNotFoundException if file not found
     */
    private void ScaffoldGeneratorSpeedTest(boolean anIsSchuffenhauerScaffoldCalculated, boolean anIsRingCalculated, boolean anIsRemoveRingCalculated,
                                            boolean anIsApplySchuffenhauerCalculated, boolean anIsPictureCreated, int aPictureNumber) throws IOException, CDKException {
        /*Counter*/
        int tmpExceptionCounter = 0;
        int tmpNumberCounter = 0;
        int tmpSkipCounter = 0;
        int tmpCounter = 0;
        //Number of molecules that have more than 10 rings and were skipped in the getIterativeRemoval speed test
        //Number of molecules where more than 1000 rings were detected and skipped in the calculaterRemoveRingsSpeedTest().
        int tmpFusedRingCounter = 0;
        /*Loading and reading the library*/
        File tmpResourcesDirectory = new File("src/test/resources/COCONUT_DB.sdf");
        IteratingSDFReader tmpReader = new IteratingSDFReader( new FileInputStream(tmpResourcesDirectory), DefaultChemObjectBuilder.getInstance());
        //Start timer
        long tmpStartTime = System.nanoTime();
        /*Start report*/
        System.out.println("-----START REPORT-----");
        if(anIsSchuffenhauerScaffoldCalculated == true && anIsRingCalculated == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds are calculated for all molecules.");
        }
        if(anIsRingCalculated == true && anIsRemoveRingCalculated == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
        }
        if(anIsRemoveRingCalculated == true && anIsApplySchuffenhauerCalculated == false){
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
            System.out.println("In addition, the molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated.");
        }
        if(anIsApplySchuffenhauerCalculated == true){
            System.out.println("In this test, the molecules are decomposed according to the Schuffenhauer rules (except rule 7)");
        }
        /*Going through the library*/
        while (tmpReader.hasNext()) {
            String tmpCoconutID = null;
            IAtomContainer tmpMolecule = (IAtomContainer) tmpReader.next();
            tmpCoconutID = tmpMolecule.getProperty("coconut_id");
            try {
                /*Calculate SchuffenhauerScaffolds*/
                if(anIsSchuffenhauerScaffoldCalculated == true && anIsRingCalculated == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
                    /*Generate control picture*/
                    if(anIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                        /*Generate and save molecule picture*/
                        BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                        File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                    }
                }
                /*Calculate SchuffenhauerScaffolds and Rings*/
                if(anIsRingCalculated == true && anIsRemoveRingCalculated == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
                    /*Generate control pictures*/
                    if(anIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        IAtomContainer tmpRing = tmpRings.get(tmpRings.size()-1);
                        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                        /*Generate and save molecule picture*/
                        BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                        File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                        /*Generate and save ring picture*/
                        BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png").mkdirs();
                        File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png");
                        ImageIO.write(tmpImgRing, "png", tmpOutputRing);
                    }
                }
                /*Calculate SchuffenhauerScaffolds, Rings and the molecules for which the rings have been removed from the Schuffenhauer scaffolds*/
                if(anIsRemoveRingCalculated == true && anIsApplySchuffenhauerCalculated == false) {
                    IAtomContainer tmpSchuff = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, true, ElectronDonation.cdk());
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuff,true);
                    /*Skip all molecules with more than 1000 rings*/
                    if(tmpRings.size() > 100) {
                        System.out.println(tmpCoconutID);
                        System.out.println("Skipped Ring size:  " + tmpRings.size());
                        tmpNumberCounter++;
                        tmpSkipCounter++;
                        continue;
                    }
                    boolean tmpIsFused = false;
                    boolean tmpIsRuleSeven = false;
                    List<IAtomContainer> tmpRemovableRings = new ArrayList<>(tmpRings.size());
                    for(IAtomContainer tmpRing : tmpRings) {
                        if(scaffoldGenerator.isRingRemovable(tmpRing, tmpRings, tmpSchuff) && scaffoldGenerator.isRingTerminal(tmpSchuff, tmpRing)) {
                            tmpRemovableRings.add(tmpRing);
                        }
                        /*Detect molecules with aromatic fused ring systems*/
                        if(scaffoldGenerator.hasFusedAromaticRings(tmpRing,tmpRings,tmpSchuff) == true) {
                            tmpIsFused = true;
                        }
                        IAtomContainer tmpRemoveMol = scaffoldGenerator.removeRing(tmpSchuff, tmpRing);
                        /*Generate control pictures*/
                        if(anIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                            /*Generate and save molecule picture*/
                            BufferedImage tmpImgMol = tmpGenerator.depict(tmpSchuff).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                            File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                            ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                            /*Generate and save ring picture*/
                            BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png").mkdirs();
                            File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png");
                            ImageIO.write(tmpImgRing, "png", tmpOutputRing);
                            /*Generate and save removed ring picture*/
                            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemoveMol).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png").mkdirs();
                            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png");
                            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                        }
                    }

                    if(scaffoldGenerator.ruleSevenTest(tmpSchuff, tmpRemovableRings, true, ElectronDonation.cdk())) {
                        tmpIsRuleSeven = true;
                    }


                    /*Detected molecule with aromatic fused ring system*/
                    if(tmpIsFused == true) {
                        tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
                        //CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
                        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMolecule).toImg();
                        /*Save the picture*/
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/AromaticFused/" + tmpCoconutID + ".png").mkdirs();
                        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/AromaticFused/" + tmpCoconutID + ".png");
                        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);

                        System.out.println("Fused ring detected " + tmpCoconutID);
                        tmpFusedRingCounter++;
                    }
                    if(tmpIsRuleSeven == true) {
                        tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
                        //CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
                        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                        BufferedImage tmpImgMod = tmpGenerator.depict(tmpMolecule).toImg();
                        /*Save the picture*/
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/RuleSeven/" + tmpCoconutID + ".png").mkdirs();
                        File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/RuleSeven/" + tmpCoconutID + ".png");
                        ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
                        System.out.println("Rule Seven Positive" + tmpCoconutID);
                        tmpFusedRingCounter++;
                    }
                }

                /*Calculate a list of molecules with iteratively removed terminal rings*/
                if(anIsApplySchuffenhauerCalculated == true) {
                    //if(scaffoldGenerator.getRings(tmpMolecule,false).size() > 100000) {
                     //   tmpSkipCounter++;
                     //   continue;
                    //}
                    boolean tmpIsConspicuous = false;
                    IAtomContainer tmpSchuff = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule, false, null);
                    int tmpRingNumber =  scaffoldGenerator.getRings(tmpSchuff, false).size();
                    float tmpRingAtomRatio = (float) tmpSchuff.getAtomCount() / tmpRingNumber;
                    if(tmpRingAtomRatio < 1.0 ) {
                        tmpIsConspicuous = true;
                    }

                    List<IAtomContainer> tmpIterations = scaffoldGenerator.applySchuffenhauerRules(tmpMolecule, true, ElectronDonation.cdk(), false);
                    if(anIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        int tmpIterationCounter = 0;
                        for(IAtomContainer tmpIteration : tmpIterations) {
                            /*Generate control pictures*/
                            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                            /*Generate and save molecule picture*/
                            BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                            File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                            ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                            /*Generate picture of the Iteration*/
                            BufferedImage tmpImgIter = tmpGenerator.depict(tmpIteration).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration" + tmpIterationCounter + ".png").mkdirs();
                            File tmpOutputIter = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration" + tmpIterationCounter + ".png");
                            ImageIO.write(tmpImgIter, "png", tmpOutputIter);
                            tmpIterationCounter++;
                        }
                    }
                    if(tmpIsConspicuous == true) {
                        if(tmpNumberCounter != 141586){
                            tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
                            //CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
                            DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                            BufferedImage tmpImgMod = tmpGenerator.depict(tmpMolecule).toImg();
                            /*Save the picture*/
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Conspicuous/" + tmpCoconutID + tmpIsMoleculeRemovable + ".png").mkdirs();
                            File tmpOutputMod = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Conspicuous/" + tmpCoconutID + tmpIsMoleculeRemovable + ".png");
                            ImageIO.write(tmpImgMod, "png" ,tmpOutputMod);
                            System.out.println("Conspicuous: " + tmpCoconutID);
                            tmpFusedRingCounter++;
                        }
                    }
                }

                if(tmpCounter == 10000){
                    tmpCounter = 0;
                    System.out.println("At molecule Number: " + tmpNumberCounter);
                }


                /*First status report*/
                if(tmpNumberCounter == 100000) {
                    System.out.println("-----STATUS REPORT(1/4)-----");
                    System.out.println("A quarter of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
                /*Second status report*/
                if(tmpNumberCounter == 200000) {
                    System.out.println("-----STATUS REPORT(2/4)-----");
                    System.out.println("A half of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
                /*Third status report*/
                if(tmpNumberCounter == 300000) {
                    System.out.println("-----STATUS REPORT(3/4)-----");
                    System.out.println("Two thirds of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
            }
            /*Count exceptions*/
            catch(Exception e) {
                System.out.println("Exception at number: " + tmpNumberCounter);
                System.out.println("COCONUT ID: " + tmpCoconutID);
                SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
                System.out.println("SMILES: " + tmpSmilesGenerator.create(tmpMolecule));
                /*Print out the exception stack trace*/
                StringWriter sw = new StringWriter();
                PrintWriter pw = new PrintWriter(sw);
                e.printStackTrace(pw);
                String sStackTrace = sw.toString();
                System.out.println(sStackTrace);

                tmpExceptionCounter++;


                /*Generate control pictures*/
                DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                /*Generate and save molecule picture*/
                tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
                BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Exception/" + tmpCoconutID + ".png").mkdirs();
                File tmpOutputMol = new File(System.getProperty("user.dir") +  "/scaffoldTestOutput/Exception/" + tmpCoconutID + ".png");
                ImageIO.write(tmpImgMol, "png", tmpOutputMol);
            }
            tmpCounter++;
            tmpNumberCounter++;
        }
        /*End report*/
        System.out.println("-----END REPORT-----");
        System.out.println("All molecules completed");
        System.out.println("Total number of exceptions: " + tmpExceptionCounter);
        System.out.println("total Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
        System.out.println("Number of skipped molecules" + tmpSkipCounter);
        System.out.println("Number of fused molecules" + tmpFusedRingCounter);
    }
    //</editor-fold>

    /**
     * Loads a mol file of a specific path and returns it as IAtomContainer.
     * Supports V2000 and V3000 mol files.
     * @param aFilePath Path of the molecule to be loaded
     * @return IAtomContainer of the charged molecule
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     */
    private IAtomContainer loadMolFile(String aFilePath) throws IOException, CDKException {
        /*Get molecule path*/
        File tmpResourcesDirectory = new File(aFilePath);
        BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory));
        /*Get mol file version*/
        FormatFactory tmpFactory = new FormatFactory();
        IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
        IAtomContainer tmpMolecule = new AtomContainer();
        /*Load V2000 mol file*/
        if(tmpFormat.getReaderClassName().contains("V2000")) {
            MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
            IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
            tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            /*Load V3000 mol file*/
        } else if(tmpFormat.getReaderClassName().contains("V3000")) {
            MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
            IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
            tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
        }
        return tmpMolecule;
    }
}
