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

import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmilesParser;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.List;
import java.util.concurrent.TimeUnit;

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
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates the SchuffenhauserScaffolds with getSchuffenhauerScaffold().
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerScaffoldTest() throws IOException, CDKException, CloneNotSupportedException {
        for (int tmpCount = 1; tmpCount < 13; tmpCount++) {
            String tmpFileName = "Test"+ tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory));
            //Get mol file version
            FormatFactory tmpFactory = new FormatFactory();
            IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
            IAtomContainer tmpMolecule = new AtomContainer();
            //Load V2000 mol file
            if(tmpFormat.getReaderClassName().contains("V2000")) {
                MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            //Load V3000 mol file
            } else if(tmpFormat.getReaderClassName().contains("V3000")) {
                MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            }
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate picture of the SchuffenhauerScaffold
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
            BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpSchuffenhauerScaffold).toImg();
            //Save the picture
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png").mkdirs();
            File tmpOutputSchuff = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png");
            ImageIO.write(tmpImgSchuff, "png" ,tmpOutputSchuff);
        }
    }
    /**
     * Test of ScaffoldGenerator.getSchuffenhauerScaffold() with SMILES.
     * Loads Flucloxacillin(Test11) as SMILES and generates SchuffenhauerScaffold.
     * All generated scaffolds are saved as images in a subfolder of the scaffoldTestOutput folder.
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getSchuffenhauerSMILESTest() throws IOException, CDKException, CloneNotSupportedException {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerSMILES = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
        //Generate picture of the SchuffenhauerScaffold
        DepictionGenerator tmpGenerator = new DepictionGenerator();
        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerSMILES).toImg();
        //Save the picture
        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Test11/SchuffenhauerSMILES.png").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Test11/SchuffenhauerSMILES.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);

    }
    /**
     * Test of Cycles.mcb() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates the rings of the SchuffenhauerScaffold with getRings().
     * All generated Rings are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRingsTest() throws IOException, CDKException, CloneNotSupportedException {
        for (int tmpCount = 1; tmpCount < 13; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory));
            //Get mol file version
            FormatFactory tmpFactory = new FormatFactory();
            IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
            IAtomContainer tmpMolecule = new AtomContainer();
            //Load V2000 mol file
            if (tmpFormat.getReaderClassName().contains("V2000")) {
                MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            //Load V3000 mol file
            } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            }
            //Generate the SchuffenhauerScaffold
            tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
            //Generate pictures of the rings
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                //Save the picture
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
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void removeRingTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 2; tmpCount < 13; tmpCount++) {
            String tmpFileName = "Test"+tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            IAtomContainer tmpMolecule;
            try (BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory))) {
                //Get mol file version
                FormatFactory tmpFactory = new FormatFactory();
                IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
                tmpMolecule = new AtomContainer();
                //Load V2000 mol file
                if (tmpFormat.getReaderClassName().contains("V2000")) {
                    MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                //Load V3000 mol file
                } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                    MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                }
            }
            //Generate SchuffenhauerScaffold
            IAtomContainer tmpSchuffenhauerScaffold = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
            //Generate Rings
            List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule, true);
            int tmpCounter = 1;
            for (IAtomContainer tmpRing : tmpRings) {
                //Generate SchuffenhauerScaffold with removed ring
                IAtomContainer tmpRemovedSchuff = scaffoldGenerator.removeRing(tmpSchuffenhauerScaffold, tmpRing);
                //Generate picture of the SchuffenhauerScaffold with removed ring
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
                //Save the picture
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/RingRemove" + tmpCounter + ".png").mkdirs();
                File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/RingRemove" + tmpCounter + ".png");
                ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                tmpCounter++;
            }
        }
    }

    /**
     * Test of isRingTerminal() with V2000 and V3000 mol files.
     * Loads the 12 Test(Test1.mol-Test12.mol) molfiles from the Resources folder and creates for each generated terminal ring, the corresponding total molecule with removed ring.
     * All generated molecules are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the input file.
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void isRingTerminalTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 2; tmpCount < 13; tmpCount++) {
            String tmpFileName = "Test"+tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            IAtomContainer tmpMolecule;
            try (BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory))) {
                //Get mol file version
                FormatFactory tmpFactory = new FormatFactory();
                IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
                tmpMolecule = new AtomContainer();
                //Load V2000 mol file
                if (tmpFormat.getReaderClassName().contains("V2000")) {
                    MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                //Load V3000 mol file
                } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                    MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                }
            }
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
                    //Generate picture of the SchuffenhauerScaffold with removed ring
                    DepictionGenerator tmpGenerator = new DepictionGenerator();
                    tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                    BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemovedSchuff).toImg();
                    //Save the picture
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
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getIterativeRemovalTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 1; tmpCount < 14; tmpCount++) {
            String tmpFileName = "Test"+ tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            IAtomContainer tmpMolecule;
            try (BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory))) {
                //Get mol file version
                FormatFactory tmpFactory = new FormatFactory();
                IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
                tmpMolecule = new AtomContainer();
                //Load V2000 mol file
                if (tmpFormat.getReaderClassName().contains("V2000")) {
                    MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                //Load V3000 mol file
                } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                    MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                }
            }
            //Generate a list of molecules with iteratively removed terminal rings
            List<IAtomContainer> tmpMolecules = scaffoldGenerator.getIterativeRemoval(tmpMolecule);
            int tmpCounter = 1;
            for (IAtomContainer tmpIterative : tmpMolecules) {
                //Generate picture of the molecule
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
                //Save the picture
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
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRemovalTreeTest() throws CDKException, CloneNotSupportedException, IOException {
        for (int tmpCount = 1; tmpCount < 14; tmpCount++) {
            String tmpFileName = "Test" + tmpCount;
            //Get molecule path
            //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
            File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
            IAtomContainer tmpMolecule;
            try (BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory))) {
                //Get mol file version
                FormatFactory tmpFactory = new FormatFactory();
                IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
                tmpMolecule = new AtomContainer();
                //Load V2000 mol file
                if (tmpFormat.getReaderClassName().contains("V2000")) {
                    MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                    //Load V3000 mol file
                } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                    MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                    IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                    tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                }
            }
            //Generate a tree of molecules with iteratively removed terminal rings
            TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
            int tmpCounter = 0;
            while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
                tmpCounter++;
                TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
                //Save the picture
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpMoleculeNode.data).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter  + "Level" + tmpMoleculeNode.getLevel() + ".png").mkdirs();
                File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter +  "Level" + tmpMoleculeNode.getLevel() + ".png");
                ImageIO.write(tmpSecImgRemove, "png", tmpSecOutputRemove);
            }
            /**
            //Code to check the tree structure
            TreeNode<IAtomContainer> tmpMoleculeTree = scaffoldGenerator.getRemovalTree(tmpMolecule);
            //First Level
            for(int tmpCounter = 0;tmpCounter < tmpMoleculeTree.children.size(); tmpCounter++) {
                IAtomContainer tmpFirstMolecule = tmpMoleculeTree.children.get(tmpCounter).data;
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                BufferedImage tmpImgRemove = tmpGenerator.depict(tmpFirstMolecule).toImg();
                //Save the picture
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter + ".png").mkdirs();
                File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter + ".png");
                ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                //Second Level
                for(int tmpSecondCounter = 0;tmpSecondCounter < tmpMoleculeTree.children.get(tmpCounter).children.size(); tmpSecondCounter++ ){
                    IAtomContainer tmpSecondMolecule = tmpMoleculeTree.children.get(tmpCounter).children.get(tmpSecondCounter).data;
                    BufferedImage tmpSecImgRemove = tmpGenerator.depict(tmpSecondMolecule).toImg();
                    //Save the picture
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter + tmpSecondCounter + ".png").mkdirs();
                    File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TreeTest" + tmpCounter + tmpSecondCounter +  ".png");
                    ImageIO.write(tmpSecImgRemove, "png", tmpSecOutputRemove);
                }
            }*/
        }
    }

    /**
     * Test of getRemovalTree() with V2000 and V3000 mol files.
     * Loads one molecule(insert in tmpFileName)  from the Resources folder, iteratively removes the terminal rings and saves the molecules in a tree.
     * Saves the parent and the children of one Node and saves them as images.
     * Set file with: tmpFileName
     * Set Node with: tmpTestNumber
     * @throws IOException if file format cant be detected
     * @throws CDKException if file cant be read
     * @throws CloneNotSupportedException if cloning is not possible
     */
    @Test
    public void getRemovalTreeStructureTest() throws CDKException, CloneNotSupportedException, IOException {
        String tmpFileName = "Test13"; //File to be tested
        int tmpTestNumber = 2; //Node whose children and parent are to be displayed
        //Get molecule path
        //InputStream tmpInputStream = ScaffoldGenerator.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
        File tmpResourcesDirectory = new File("src/test/resources/" + tmpFileName + ".mol");
        IAtomContainer tmpMolecule;
        try (BufferedInputStream tmpInputStream = new BufferedInputStream(new FileInputStream(tmpResourcesDirectory))) {
            //Get mol file version
            FormatFactory tmpFactory = new FormatFactory();
            IChemFormat tmpFormat = tmpFactory.guessFormat(tmpInputStream);
            tmpMolecule = new AtomContainer();
            //Load V2000 mol file
            if (tmpFormat.getReaderClassName().contains("V2000")) {
                MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                //Load V3000 mol file
            } else if (tmpFormat.getReaderClassName().contains("V3000")) {
                MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = DefaultChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            }
        }
        //Generate a tree of molecules with iteratively removed terminal rings
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(scaffoldGenerator.getRemovalTree(tmpMolecule));
        int tmpCounter = 0;
        while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            if(tmpCounter == tmpTestNumber -1){
                //Save the picture of the test Node
                DepictionGenerator tmpGenerator = new DepictionGenerator();
                BufferedImage tmpNodeImg = tmpGenerator.depict(tmpMoleculeNode.data).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TestNode.png").mkdirs();
                File tmpNodeFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/TestNode.png");
                ImageIO.write(tmpNodeImg, "png", tmpNodeFile);
                //Save the picture of the parent
                BufferedImage tmpParentImg = tmpGenerator.depict(tmpMoleculeNode.parent.data).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/ParentNode.png").mkdirs();
                File tmpParentFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/ParentNode.png");
                ImageIO.write(tmpParentImg, "png", tmpParentFile);
                //Save pictures of the children
                int tmpChildCounter = 0;
                for(TreeNode<IAtomContainer> tmpNode : tmpMoleculeNode.children) {
                    tmpChildCounter++;
                    BufferedImage tmpChildImg = tmpGenerator.depict(tmpNode.data).toImg();
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/ChildNode" + tmpChildCounter + ".png").mkdirs();
                    File tmpChildFile = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/ChildNode" + tmpChildCounter + ".png");
                    ImageIO.write(tmpChildImg, "png", tmpChildFile);
                }
                break;
            }
            tmpCounter++;
        }
    }
    /**
     * Speed test for the getSchuffenhauerScaffold(), getRing() and removeRing() Method with over 400000 molecules from the COCONUT DB.
     * Which methods are tested can be set via the booleans.
     * To perform the test download the COCONUT DB(https://coconut.naturalproducts.net/download) and add the COCONUT_DB.sdf file to src\test\resources
     * @throws FileNotFoundException if file not found
     */
    @Test
    public void ScaffoldGeneratorSpeedTest() throws FileNotFoundException {
        //Options
        boolean tmpCalculateSchuffenhauerScaffold = true; //Generate SchuffenhauerScaffolds
        boolean tmpCalculateRings = false; //Calculate Rings. The Schuffenhauer scaffolds are also calculated for this.
        boolean tmpCalculateRemoveRings = true; //The molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated. The Schuffenhauer scaffolds and the Rings are also calculated for this.
        boolean tmpIterativeRemoval = false; //Creates all molecules that are generated by iterative terminal removal of rings.
        boolean tmpGetPicture = false; //Show control pictures from one molecule.
        int tmpPictureNumber = 2461; //Number of the molecule from which control images are to be taken(from 0 to 406000)
        //Counter
        int tmpExceptionCounter = 0;
        int tmpNumberCounter = 0;
        int tmpToManyRingsCounter = 0;
        //Loading and reading the library
        File tmpResourcesDirectory = new File("src/test/resources/COCONUT_DB.sdf");
        IteratingSDFReader tmpReader = new IteratingSDFReader( new FileInputStream(tmpResourcesDirectory), DefaultChemObjectBuilder.getInstance());
        //Start timer
        long tmpStartTime = System.nanoTime();
        //Start report
        System.out.println("-----START REPORT-----");
        if(tmpCalculateSchuffenhauerScaffold == true && tmpCalculateRings == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds are calculated for all molecules.");
        }
        if(tmpCalculateRings == true && tmpCalculateRemoveRings == false) {
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
        }
        if(tmpCalculateRemoveRings == true && tmpIterativeRemoval == false){
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
            System.out.println("In addition, the molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated.");
        }
        if(tmpIterativeRemoval == true){
            System.out.println("In this test, the Schuffenhauer scaffolds and their rings are calculated for all molecules.");
            System.out.println("The molecules for which the rings have been removed from the Schuffenhauer scaffolds are calculated.");
            System.out.println("In addition, all molecules with with iteratively removed terminal rings are calculated");
        }
        // Going through the library
        while (tmpReader.hasNext()) {
            try {
                IAtomContainer tmpMolecule = (IAtomContainer) tmpReader.next();
                //Calculate SchuffenhauerScaffolds
                if(tmpCalculateSchuffenhauerScaffold == true && tmpCalculateRings == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    //Generate control picture
                    if(tmpGetPicture && (tmpPictureNumber) == tmpNumberCounter) {
                        DepictionGenerator tmpGenerator = new DepictionGenerator();
                        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                        //Generate and save molecule picture
                        BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                        File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                    }
                }
                //Calculate SchuffenhauerScaffolds and Rings
                if(tmpCalculateRings == true && tmpCalculateRemoveRings == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
                    //Generate control pictures
                    if(tmpGetPicture && (tmpPictureNumber) == tmpNumberCounter) {
                        IAtomContainer tmpRing = tmpRings.get(tmpRings.size()-1);
                        DepictionGenerator tmpGenerator = new DepictionGenerator();
                        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                        //Generate and save molecule picture
                        BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                        File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                        //Generate and save ring picture
                        BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                        new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png").mkdirs();
                        File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png");
                        ImageIO.write(tmpImgRing, "png", tmpOutputRing);
                    }
                }
                //Calculate SchuffenhauerScaffolds, Rings and the molecules for which the rings have been removed from the Schuffenhauer scaffolds
                if(tmpCalculateRemoveRings == true && tmpIterativeRemoval == false) {
                    tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                    List<IAtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
                    for(IAtomContainer tmpRing : tmpRings) {
                        IAtomContainer tmpRemoveMol = scaffoldGenerator.removeRing(tmpMolecule, tmpRing);
                        //Generate control pictures
                        if(tmpGetPicture && (tmpPictureNumber) == tmpNumberCounter) {
                            DepictionGenerator tmpGenerator = new DepictionGenerator();
                            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                            //Generate and save molecule picture
                            BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                            File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                            ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                            //Generate and save ring picture
                            BufferedImage tmpImgRing = tmpGenerator.depict(tmpRing).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png").mkdirs();
                            File tmpOutputRing = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRing.png");
                            ImageIO.write(tmpImgRing, "png", tmpOutputRing);
                            //Generate and save removed ring picture
                            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpRemoveMol).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png").mkdirs();
                            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestRingRemoved.png");
                            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
                        }
                    }
                }
                //Calculate a list of molecules with iteratively removed terminal rings
                if(tmpIterativeRemoval == true) {
                    if(scaffoldGenerator.getRings(tmpMolecule,false).size() > 10) {
                        tmpToManyRingsCounter++;
                        continue;
                    }
                    List<IAtomContainer> tmpIterations = scaffoldGenerator.getIterativeRemoval(tmpMolecule);
                    if(tmpGetPicture && (tmpPictureNumber) == tmpNumberCounter) {
                        for(IAtomContainer tmpIteration : tmpIterations) {
                        //Generate control pictures
                            DepictionGenerator tmpGenerator = new DepictionGenerator();
                            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                            //Generate and save molecule picture
                            BufferedImage tmpImgMol = tmpGenerator.depict(tmpMolecule).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png").mkdirs();
                            File tmpOutputMol = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestMol.png");
                            ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                            //Generate picture of the Iteration
                            BufferedImage tmpImgIter = tmpGenerator.depict(tmpIteration).toImg();
                            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration.png").mkdirs();
                            File tmpOutputIter = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest/SpeedTestIteration.png");
                            ImageIO.write(tmpImgIter, "png", tmpOutputIter);
                        }
                    }
                }

                //First status report
                if(tmpNumberCounter == 100000) {
                    System.out.println("-----STATUS REPORT(1/4)-----");
                    System.out.println("A quarter of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
                //Second status report
                if(tmpNumberCounter == 200000) {
                    System.out.println("-----STATUS REPORT(2/4)-----");
                    System.out.println("A half of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
                //Third status report
                if(tmpNumberCounter == 300000) {
                    System.out.println("-----STATUS REPORT(3/4)-----");
                    System.out.println("Two thirds of all molecules completed");
                    System.out.println("Number of exceptions: " + tmpExceptionCounter);
                    System.out.println("Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
                }
                tmpNumberCounter++;
            }
            //Count exceptions
            catch(Exception e) {
                tmpExceptionCounter++;
            }
        }
        //End report
        System.out.println("-----END REPORT-----");
        System.out.println("All molecules completed");
        System.out.println("Total number of exceptions: " + tmpExceptionCounter);
        System.out.println("total Runtime: " + TimeUnit.NANOSECONDS.toSeconds((System.nanoTime() - tmpStartTime)) + " seconds");
    }
}
