/*
 * Copyright (c) 2022 Julian Zander, Jonas Schaub, Achim Zielesny, Christoph Steinbeck
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

package de.unijena.cheminf.scaffolds;

import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.io.formats.IChemFormat;

import javax.imageio.ImageIO;
import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
  * JUnit test class for the org.openscience.cdk.fragment.MurckoFragmenter
  *
 * @author Julian Zander, Jonas Schaub (zanderjulian@gmx.de, jonas.schaub@uni-jena.de)
  * @version 1.0.0.1
  */
public class MurckoFragmenterTest {
    private MurckoFragmenter fragmenter;
    @Before
    public void initFragmenter() {
        fragmenter = new MurckoFragmenter(false,1);
    }

     /**
      * Tests the MurckoFragmenter: Loads the 23 Test(Test1.mol-Test23.mol) molfiles(V2000 and V3000) from the Resources folder and creates all fragments that can be generated with the MurckoFragmenter.
      * The unmodified molecule (original) and all generated fragments are saved as images in a subfolder of the generated scaffoldTestOutput folder.
      * The subfolder has the name of the input file.
      * @throws IOException if file format cant be detected
      * @throws CDKException if file cant be read
      */
    @Test
    public void testFragmenter() throws IOException, CDKException {
        for (int tmpCount = 1; tmpCount < 8; tmpCount++) {
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
            //Generate picture of the original molecule
            DepictionGenerator tmpGenerator = new DepictionGenerator();
            tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
            BufferedImage tmpImgOri = tmpGenerator.depict(tmpMolecule).toImg();
            //Save picture of the original molecule
            new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName).mkdirs();
            File tmpOutputOri = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Original.png");
            ImageIO.write(tmpImgOri, "png" ,tmpOutputOri);
            //Generate fragments, rings and frameworks
            fragmenter.setComputeRingFragments(true);
            fragmenter.generateFragments(tmpMolecule);
            IAtomContainer[] tmpFragments = fragmenter.getFragmentsAsContainers();
            IAtomContainer[] tmpRings = fragmenter.getRingSystemsAsContainers();
            IAtomContainer[] tmpFrameworks = fragmenter.getFrameworksAsContainers();
            //Generate and save pictures of the fragments
            int tmpCountFra = 1;
            for(IAtomContainer tmpFragment : tmpFragments) {
                BufferedImage tmpImgFra = tmpGenerator.depict(tmpFragment).toImg();
                new File(System.getProperty("user.dir")+"/scaffoldTestOutput/" + tmpFileName).mkdirs();
                File tmpOutputFra = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName + "/Fragment" + tmpCountFra+".png");
                ImageIO.write(tmpImgFra, "png" ,tmpOutputFra);
                tmpCountFra++;
            }
            //Generate and save pictures of the rings
            int tmpCountRgs = 1;
            for(IAtomContainer tmpRing : tmpRings) {
                BufferedImage tmpImgRgs = tmpGenerator.depict(tmpRing).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName).mkdirs();
                File tmpOutputRgs = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName+"/Ring" + tmpCountRgs + ".png");
                ImageIO.write(tmpImgRgs, "png" ,tmpOutputRgs);
                tmpCountRgs++;
            }
            //Generate and save pictures of the frameworks
            int tmpCountFrw = 1;
            for(IAtomContainer tmpFramework : tmpFrameworks) {
                BufferedImage tmpImgFrw = tmpGenerator.depict(tmpFramework).toImg();
                new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName).mkdirs();
                File tmpOutputFrw = new File(System.getProperty("user.dir") + "/scaffoldTestOutput/" + tmpFileName+"/Framework" + tmpCountFrw+".png");
                ImageIO.write(tmpImgFrw, "png" ,tmpOutputFrw);
                tmpCountFrw++;
            }
        }
    }
}