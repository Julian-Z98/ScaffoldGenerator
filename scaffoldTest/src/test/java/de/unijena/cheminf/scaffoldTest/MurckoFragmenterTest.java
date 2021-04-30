/*
 *
 *  * MIT License
 *  *
 *  * Copyright (c) 2021 Julian Zander, Jonas Schaub,  Achim Zielesny
 *  *
 *  * Permission is hereby granted, free of charge, to any person obtaining a copy
 *  * of this software and associated documentation files (the "Software"), to deal
 *  * in the Software without restriction, including without limitation the rights
 *  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  * copies of the Software, and to permit persons to whom the Software is
 *  * furnished to do so, subject to the following conditions:
 *  *
 *  * The above copyright notice and this permission notice shall be included in all
 *  * copies or substantial portions of the Software.
 *  *
 *  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  * SOFTWARE.
 *
 */

package de.unijena.cheminf.scaffoldTest;

import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

/**
  * JUnit test class for the org.openscience.cdk.fragment.MurckoFragmenter
  */
public class MurckoFragmenterTest {
    private MurckoFragmenter fragmenter;
    @Before
    public void initFragmenter() {
        fragmenter = new MurckoFragmenter(false,1);
    }

     /**
      * Tests the MurckoFragmenter: Loads the 11 Test(Test1.mol-Test11.mol) molfiles(V2000 and V3000) from the Resources folder and creates all fragments that can be generated with the MurckoFragmenter.
      * The unmodified molecule (original) and all generated fragments are saved as images in a subfolder of the generated scaffoldTestOutput folder.
      * The subfolder has the name of the input file.
      * @throws IOException if file format cant be detected
      * @throws CDKException if file cant be read
      */
    @Test
    public void testFragmenter() throws IOException, CDKException {
        for (int tmpCount = 1; tmpCount < 12; tmpCount++) {
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
                IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
                //Load V3000 mol file
            } else if(tmpFormat.getReaderClassName().contains("V3000")){
                MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
                IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
                tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
            }
            //Generate picture of the original molecule
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
            BufferedImage tmpImgOri = tmpGenerator.depict(tmpMolecule).toImg();
            //Save picture of the original molecule
            new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Original.png").mkdirs();
            File tmpOutputOri = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Original.png");
            ImageIO.write(tmpImgOri, "png" ,tmpOutputOri);
            //Generate fragments, rings and frameworks
            fragmenter.setComputeRingFragments(true);
            fragmenter.generateFragments(tmpMolecule);
            IAtomContainer[] tmpFragments = fragmenter.getFragmentsAsContainers();
            IAtomContainer[] tmpRings = fragmenter.getRingSystemsAsContainers();
            IAtomContainer[] tmpFrameworks = fragmenter.getFrameworksAsContainers();
            //Generate and save pictures of the fragments
            int tmpCountFra = 0;
            for(IAtomContainer tmpFragment : tmpFragments) {
                tmpCountFra++;
                BufferedImage tmpImgFra = tmpGenerator.depict(tmpFragment).toImg();
                new File(System.getProperty("user.dir")+"/src/test/scaffoldTestOutput/" + tmpFileName + "/Fragment"+tmpCountFra + ".png").mkdirs();
                File tmpOutputFra = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Fragment" + tmpCountFra+".png");
                ImageIO.write(tmpImgFra, "png" ,tmpOutputFra);
            }
            //Generate and save pictures of the rings
            int tmpCountRgs = 0;
            for(IAtomContainer tmpRing : tmpRings) {
                tmpCountRgs++;
                BufferedImage tmpImgRgs = tmpGenerator.depict(tmpRing).toImg();
                new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Ring" + tmpCountRgs + ".png").mkdirs();
                File tmpOutputRgs = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName+"/Ring" + tmpCountRgs + ".png");
                ImageIO.write(tmpImgRgs, "png" ,tmpOutputRgs);
            }
            //Generate and save pictures of the frameworks
            int tmpCountFrw = 0;
            for(IAtomContainer tmpFramework : tmpFrameworks) {
                tmpCountFrw++;
                BufferedImage tmpImgFrw = tmpGenerator.depict(tmpFramework).toImg();
                new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName+"/Framework" + tmpCountFrw+".png").mkdirs();
                File tmpOutputFrw = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName+"/Framework" + tmpCountFrw+".png");
                ImageIO.write(tmpImgFrw, "png" ,tmpOutputFrw);
            }
        }
    }
}