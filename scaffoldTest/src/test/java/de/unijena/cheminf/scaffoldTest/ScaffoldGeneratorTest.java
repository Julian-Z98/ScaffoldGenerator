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
            List<AtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule,true);
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
     * Test of Cycles.mcb() with V2000 and V3000 mol files.
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
            List<AtomContainer> tmpRings = scaffoldGenerator.getRings(tmpMolecule, true);
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
     * Speed test for the getSchuffenhauerScaffold() Method with over 400000 molecules from the COCONUT DB.
     * @throws FileNotFoundException if file not found
     */
    @Test
    public void getSchuffenhauerScaffoldSpeedTest() throws FileNotFoundException {
        int tmpExceptionCounter = 0;
        int tmpNumberCounter = 0;
        File tmpResourcesDirectory = new File("src/test/resources/COCONUT_DB.sdf");
        IteratingSDFReader reader = new IteratingSDFReader( new FileInputStream(tmpResourcesDirectory), DefaultChemObjectBuilder.getInstance());
        while (reader.hasNext()) {
            try {
                IAtomContainer tmpMolecule = (IAtomContainer) reader.next();
                tmpMolecule = scaffoldGenerator.getSchuffenhauerScaffold(tmpMolecule);
                ///**
                 tmpNumberCounter++;
                 if(tmpNumberCounter == 4242) {
                 //Generate picture of the SchuffenhauerScaffold with removed ring
                 DepictionGenerator tmpGenerator = new DepictionGenerator();
                 tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
                 BufferedImage tmpImgRemove = tmpGenerator.depict(tmpMolecule).toImg();
                 //Save the picture
                 new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest1.png").mkdirs();
                 File tmpOutputRemove = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/SpeedTest1.png");
                 ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);

                 }
                 //*/
            }
            catch(Exception e) {
                tmpExceptionCounter++;
            }
        }
        System.out.println("Number of Exceptions:"+ tmpExceptionCounter);
    }
}
