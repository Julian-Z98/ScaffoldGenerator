/*
 *
 * MIT License
 *
 * Copyright (c) 2021 Julian Zander, Jonas Schaub,  Achim Zielesny
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 *
 */

package de.unijena.cheminf.scaffoldTest;

import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
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
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
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

    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates the SchuffenhauerScaffolds with getSchuffenhauerScaffold().
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerScaffoldTest() throws CloneNotSupportedException, CDKException, IOException {
        for (int tmpCount = 1; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            /*Generate picture of the original*/
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
            CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder()).addImplicitHydrogens(tmpMolecule);
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
            /*Save the original picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Original.png").mkdirs();
            File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Original.png");
            ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate picture of the SchuffenhauerScaffold
            BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png").mkdirs();
            File tmpOutputSchuff = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png");
            ImageIO.write(tmpImgSchuff, "png" ,tmpOutputSchuff);
        }
    }
    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold() with SMILES.
     * Loads Flucloxacillin(Test11) as SMILES and generates SchuffenhauerScaffold.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format can not be detected
     * @throws CDKException if file can not be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerSMILESTest() throws IOException, CDKException, CloneNotSupportedException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerSMILES = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
        /*Generate picture of the SchuffenhauerScaffold*/
        DepictionGenerator tmpGenerator = new DepictionGenerator();
        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerSMILES).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Test11/SchuffenhauerSMILES.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Test11/SchuffenhauerSMILES.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);

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
        for (int tmpCount = 1; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate the SchuffenhauerScaffold
            tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
            /*Generate pictures of the rings*/
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
        for (int tmpCount = 2; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount ;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate Rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpSchuffenhauerScaffold, true);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                /*Generate SchuffenhauerScaffold with removed ring*/
                //tmpSchuffenhauerScaffold = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                /*Generate picture of the SchuffenhauerScaffold with removed ring*/
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate Rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule, true);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                //Check that rings are terminal
                if(scaffoldGenerator.isRingTerminal(tmpSchuffenhauerScaffold, tmpRing)) {
                    //Generate SchuffenhauerScaffold with removed ring
                    IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                    /*Generate picture of the SchuffenhauerScaffold with removed ring*/
                    DepictionGenerator tmpGenerator = new DepictionGenerator();
                    tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
        for (int tmpCount = 1; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate a list of molecules with iteratively removed terminal rings
            List<IAtomContainer> tmpMolecules = scaffoldGenerator.getIterativeRemoval(tmpMolecule);
            int tmpCounter = 1;
            for (IAtomContainer tmpIterative : tmpMolecules) {
                /*Generate picture of the molecule*/
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
        for (int tmpCount = 1; tmpCount < 21; tmpCount++) {
            String tmpFileName = "Test" + tmpCount ;
            //Load molecule from molfile
            IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
            //Generate a tree of molecules with iteratively removed terminal rings
            TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
            int tmpCounter = 0;
            while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
                TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
                /*Save the picture*/
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpMoleculeNode.getData()).toImg();
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
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
        int tmpCounter = 0;
        while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            if(tmpCounter == tmpTestNumber -1){
                /*Save the picture of the test Node*/
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                BufferedImage tmpNodeImg = tmpGenerator.depict(tmpMoleculeNode.getData()).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TestNode.png").mkdirs();
                File tmpNodeFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/TestNode.png");
                ImageIO.write(tmpNodeImg, "png", tmpNodeFile);
                /*Save the picture of the parent*/
                BufferedImage tmpParentImg = tmpGenerator.depict(tmpMoleculeNode.getParent().getData()).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ParentNode.png").mkdirs();
                File tmpParentFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Tree" + "/ParentNode.png");
                ImageIO.write(tmpParentImg, "png", tmpParentFile);
                /*Save pictures of the children*/
                int tmpChildCounter = 0;
                for(TreeNode<IAtomContainer> tmpNode : tmpMoleculeNode.getChildren()) {
                    tmpChildCounter++;
                    BufferedImage tmpChildImg = tmpGenerator.depict(tmpNode.getData()).toImg();
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
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
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
        System.out.println("Is the molecule of node 14 in the tree: " + tmpScaffoldTree.isMoleculeInTree((IAtomContainer) tmpScaffoldTree.getMatrixNode(14).getData()));
        assertEquals(true, tmpScaffoldTree.isMoleculeInTree((IAtomContainer) tmpScaffoldTree.getMatrixNode(14).getData()));
        /*getTreeNode checker*/
        System.out.println("---getTreeNode check---");
        IAtomContainer tmpTestMolecule = (IAtomContainer) tmpScaffoldTree.getMatrixNode(10).getData();
        System.out.println("Size of the test molecule: " + tmpTestMolecule.getAtomCount());
        IAtomContainer tmpResultMolecule = (IAtomContainer) tmpScaffoldTree.getTreeNode(tmpTestMolecule).getData();
        System.out.println("Size of the result molecule: " + tmpResultMolecule.getAtomCount());; //Should be the same size
        assertEquals(tmpTestMolecule.getAtomCount(), tmpResultMolecule.getAtomCount());
        /*getUniqueTreeNodes checker*/
        System.out.println("---getUniqueTreeNodes check---");
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        for(TreeNode tmpNode : tmpScaffoldTree.getUniqueTreeNodes()) {
            IAtomContainer tmpNodeMolecule = (IAtomContainer) tmpNode.getData();
            System.out.println(tmpSmilesGenerator.create(tmpNodeMolecule));
        }
        System.out.println("Number of unique tree nodes: " +  tmpScaffoldTree.getUniqueTreeNodes().size());
        assertEquals(10,  tmpScaffoldTree.getUniqueTreeNodes().size());
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
            IAtomContainer tmpNodeMolecule = (IAtomContainer) tmpNode.getData();
            /*Save the picture*/
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpNodeMolecule).toImg();
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/MatrixTest" + "/MatrixTest" + tmpTreeCounter  + ".png").mkdirs();
            File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/MatrixTest" + "/MatrixTest" + tmpTreeCounter  + ".png");
            ImageIO.write(tmpSecImgRemove, "png", tmpSecOutputRemove);
            tmpTreeCounter++;
        }
        /*RemoveNode, getRoot, hasOneSingleRootNode and isTreeConnected check*/
        System.out.println("---RemoveNode, getRoot, hasOneSingleRootNode and isTreeConnected check---");
        tmpScaffoldTree.removeNode(tmpScaffoldTree.getMatrixNode(3));
        System.out.println("Size of the tree after on Node is removed: "+ tmpScaffoldTree.getAllNodes().size());
        assertEquals(14, tmpScaffoldTree.getAllNodes().size());
        System.out.println("Root of the tree: " + tmpSmilesGenerator.create((IAtomContainer) tmpScaffoldTree.getRoot().getData()));
        assertEquals(22,  ((IAtomContainer) tmpScaffoldTree.getRoot().getData()).getAtomCount());
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
    @Test
    public void graphStreamTest() throws InterruptedException, CDKException, IOException, CloneNotSupportedException {
        String tmpFileName = "Test11" ;
        //Load molecule from molfile
        IAtomContainer tmpMolecule = this.loadMolFile("src/test/resources/" + tmpFileName + ".mol");
        //Generate a tree of molecules with iteratively removed terminal rings
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
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
        DepictionGenerator tmpGenerator = new DepictionGenerator();
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
            IAtomContainer tmpTreeNodeMolecule = (IAtomContainer) tmpTreeNode.getData();
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

    /**
     * Speed test for the getSchuffenhauerScaffold() Method with over 400000 molecules from the COCONUT DB.
     * @throws FileNotFoundException if COCONUT DB not found
     */
    @Test
    public void calculateSchuffenhauerSpeedTest() throws FileNotFoundException {
        this.ScaffoldGeneratorSpeedTest(true, false, false, false, true, 4242);
    }

    /**
     * Speed test for the getRings() Method with over 400000 molecules from the COCONUT DB.
     * getSchuffenhauerScaffold() must also be executed for all molecules.
     * @throws FileNotFoundException if COCONUT DB not found
     */
    @Test
    public void calculateRingsSpeedTest() throws FileNotFoundException {
        this.ScaffoldGeneratorSpeedTest(false, true, false, false, true, 4242);
    }

    /**
     * Speed test for the removeRing() Method with over 400000 molecules from the COCONUT DB.
     * getSchuffenhauerScaffold() and getRings() must also be executed for all molecules.
     * @throws FileNotFoundException if COCONUT DB not found
     */
    @Test
    public void calculateRemoveRingsSpeedTest() throws FileNotFoundException {
        this.ScaffoldGeneratorSpeedTest(false, false, true, false, true, 4242);
    }

    /**
     * Speed test for the getSchuffenhauerScaffold(), getRing() and removeRing() Method with over 400000 molecules from the COCONUT DB.
     * Which methods are tested can be set via the booleans.
     * To perform the test download the COCONUT DB(https://coconut.naturalproducts.net/download) and add the COCONUT_DB.sdf file to src\test\resources
     * @param aIsSchuffenhauerScaffoldCalculated Generate SchuffenhauerScaffolds
     * @param aIsRingCalculated Calculate Rings and Schuffenhauer scaffolds.
     * @param aIsRemoveRingCalculated The molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated. The Schuffenhauer scaffolds and the Rings are also calculated for this.
     * @param aIsIterativeRemovalCalculated Creates all molecules that are generated by iterative terminal removal of rings.
     * @param aIsPictureCreated Show control pictures from one molecule.
     * @param aPictureNumber Number of the molecule from which control images are to be taken(from 0 to 406000)
     * @throws FileNotFoundException if file not found
     */
    private void ScaffoldGeneratorSpeedTest(boolean aIsSchuffenhauerScaffoldCalculated, boolean aIsRingCalculated, boolean aIsRemoveRingCalculated,
                                            boolean aIsIterativeRemovalCalculated, boolean aIsPictureCreated, int aPictureNumber) throws FileNotFoundException {
        /*Counter*/
        int tmpExceptionCounter = 0;
        int tmpNumberCounter = 0;
        int tmpToManyRingsCounter = 0;
        /*Loading and reading the library*/
        File tmpResourcesDirectory = new File("src/test/resources/COCONUT_DB.sdf");
        IteratingSDFReader tmpReader = new IteratingSDFReader( new FileInputStream(tmpResourcesDirectory), DefaultChemObjectBuilder.getInstance());
        //Start timer
        long tmpStartTime = System.nanoTime();
        /*Start report*/
        System.out.println("-----START REPORT-----");
        if(aIsSchuffenhauerScaffoldCalculated == true && aIsRingCalculated == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds are calculated for all molecules.");
        }
        if(aIsRingCalculated == true && aIsRemoveRingCalculated == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
        }
        if(aIsRemoveRingCalculated == true && aIsIterativeRemovalCalculated == false){
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
            System.out.println("In addition, the molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated.");
        }
        if(aIsIterativeRemovalCalculated == true){
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
            System.out.println("The molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated.");
            System.out.println("In addition, all molecules with with iteratively removed terminal rings are calculated");
        }
        /*Going through the library*/
        while (tmpReader.hasNext()) {
            String tmpCoconutID = null;
            try {
                IAtomContainer tmpMolecule = (IAtomContainer) tmpReader.next();
                tmpCoconutID = tmpMolecule.getProperty("coconut_id");
                /*Calculate SchuffenhauerScaffolds*/
                if(aIsSchuffenhauerScaffoldCalculated == true && aIsRingCalculated == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    /*Generate control picture*/
                    if(aIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        DepictionGenerator tmpGenerator = new DepictionGenerator();
                        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                        /*Generate and save molecule picture*/
                        BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                        File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                    }
                }
                /*Calculate SchuffenhauerScaffolds and Rings*/
                if(aIsRingCalculated == true && aIsRemoveRingCalculated == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
                    /*Generate control pictures*/
                    if(aIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        IAtomContainer tmpRing = tmpRings.get(tmpRings.size()-1);
                        DepictionGenerator tmpGenerator = new DepictionGenerator();
                        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
                if(aIsRemoveRingCalculated == true && aIsIterativeRemovalCalculated == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
                    for(IAtomContainer tmpRing : tmpRings) {
                        IAtomContainer tmpRemoveMol = scaffoldGenerator.removeRing(tmpMolecule, tmpRing);
                        /*Generate control pictures*/
                        if(aIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                            DepictionGenerator tmpGenerator = new DepictionGenerator();
                            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
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
                            /*Generate and save removed ring picture*/
                            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemoveMol).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png").mkdirs();
                            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png");
                            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                        }
                    }
                }
                /*Calculate a list of molecules with iteratively removed terminal rings*/
                if(aIsIterativeRemovalCalculated == true) {
                    if(scaffoldGenerator.getRings(tmpMolecule,false).size() > 10) {
                        tmpToManyRingsCounter++;
                        continue;
                    }
                    List<IAtomContainer> tmpIterations = scaffoldGenerator.getIterativeRemoval(tmpMolecule);
                    if(aIsPictureCreated && (aPictureNumber) == tmpNumberCounter) {
                        for(IAtomContainer tmpIteration : tmpIterations) {
                        /*Generate control pictures*/
                            DepictionGenerator tmpGenerator = new DepictionGenerator();
                            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                            /*Generate and save molecule picture*/
                            BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                            File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                            ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                            /*Generate picture of the Iteration*/
                            BufferedImage tmpImgIter = tmpGenerator.depict(tmpIteration).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration.png").mkdirs();
                            File tmpOutputIter = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration.png");
                            ImageIO.write(tmpImgIter, "png", tmpOutputIter);
                        }
                    }
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
                /*Print out the exception stack trace*/
                StringWriter sw = new StringWriter();
                PrintWriter pw = new PrintWriter(sw);
                e.printStackTrace(pw);
                String sStackTrace = sw.toString();
                System.out.println(sStackTrace);
                tmpExceptionCounter++;
            }
            tmpNumberCounter++;
        }
        /*End report*/
        System.out.println("-----END REPORT-----");
        System.out.println("All molecules completed");
        System.out.println("Total number of exceptions: " + tmpExceptionCounter);
        System.out.println("total Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
    }

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
