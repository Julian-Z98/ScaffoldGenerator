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

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.MDLV3000Reader;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class murckoFragmenter {
    public static void main(String[] args) throws CDKException {
        //Load molecule
        InputStream inputStream= murckoFragmenter.class.getClassLoader().getSystemResourceAsStream("Test8.mol");// Try 6,8,9 for Schuffenhauer Test
        MDLV3000Reader reader = new MDLV3000Reader(inputStream);
        IAtomContainer testMol = reader.read(new AtomContainer());
        //Mark the carbons to which an oxygen is double bonded and save the oxygens with double bonds.
        List<Integer> addAtomList = new ArrayList<>();//Stores positions of double bounded O in the testMol
        for (IAtom tmpAtom : testMol.atoms()){
            if(tmpAtom.getSymbol().equals("O")){
                for (IBond tmpBond : tmpAtom.bonds()){
                    if(tmpBond.getElectronCount() == 4){
                        addAtomList.add(testMol.getAtomNumber(tmpAtom));
                        tmpBond.getAtom(0).setFlag(CDKConstants.DUMMY_POINTER,true);//Flag all Carbons with double bounded O
                    }
                }
            }
        }
        //Generate picture of the original molecule
        DepictionGenerator generator = new DepictionGenerator();
        generator.withSize(300, 350).withMolTitle().withTitleColor(Color.BLACK);
        BufferedImage imgOri = generator.depict(testMol).toImg();
        ImageIcon iconOri = new ImageIcon(imgOri);
        JFrame frame= new JFrame();
        frame.setLayout(new FlowLayout());
        frame.setSize(1000,800);
        JLabel lblOri = new JLabel();
        lblOri.setIcon(iconOri);
        lblOri.setText("Original");
        frame.add(lblOri);
        //Generate fragments, rings and frameworks
        MurckoFragmenter murckoFragmenter = new MurckoFragmenter(false,1);
        murckoFragmenter.setComputeRingFragments(true);
        murckoFragmenter.generateFragments(testMol);
        IAtomContainer[] fragments = murckoFragmenter.getFragmentsAsContainers();
        IAtomContainer[] rings = murckoFragmenter.getRingSystemsAsContainers();
        IAtomContainer[] frameworks = murckoFragmenter.getFrameworksAsContainers();
        //Generate pictures of the fragments
        int tmpCountFra = 0;
        for(IAtomContainer tmpFragment : fragments) {
            tmpCountFra++;
            BufferedImage tmpImgFra = generator.withBackgroundColor(Color.LIGHT_GRAY).depict(tmpFragment).toImg();
            ImageIcon tmpIconFra = new ImageIcon(tmpImgFra);
            JFrame tmpFrameFra = new JFrame();
            tmpFrameFra.setLayout(new FlowLayout());
            tmpFrameFra.setSize(320,370);
            JLabel tmpLblFra = new JLabel();
            tmpLblFra.setIcon(tmpIconFra);
            tmpLblFra.setText("Fragment "+tmpCountFra);
            frame.add(tmpLblFra);
        }
        //Generate pictures of the rings
        int tmpCountRgs = 0;
        for(IAtomContainer tmpRing : rings) {
            tmpCountRgs++;
            BufferedImage tmpImgRgs = generator.withBackgroundColor(Color.GRAY).depict(tmpRing).toImg();
            ImageIcon tmpIconRgs = new ImageIcon(tmpImgRgs);
            JFrame tmpFrameRgs = new JFrame();
            tmpFrameRgs.setLayout(new FlowLayout());
            tmpFrameRgs.setSize(320,370);
            JLabel tmpLblRgs = new JLabel();
            tmpLblRgs.setIcon(tmpIconRgs);
            tmpLblRgs.setText("Ring "+tmpCountRgs);
            frame.add(tmpLblRgs);
        }
        //Generate pictures of the frameworks
        int tmpCountFrw = 0;
        for(IAtomContainer tmpFramework : frameworks) {
            tmpCountFrw++;
            BufferedImage tmpImgFrw = generator.withBackgroundColor(Color.PINK).depict(tmpFramework).toImg();
            ImageIcon tmpIconFrw = new ImageIcon(tmpImgFrw);
            JFrame tmpFrameFrw = new JFrame();
            tmpFrameFrw.setLayout(new FlowLayout());
            tmpFrameFrw.setSize(320,370);
            JLabel tmpLblFrw = new JLabel();
            tmpLblFrw.setIcon(tmpIconFrw);
            tmpLblFrw.setText("Framework "+tmpCountFrw);
            frame.add(tmpLblFrw);
        }
        //Adds previously removed double-bonded oxygens to the framework to obtain Schuffenhauer frameworks.
        IAtomContainer tmpSFramework = frameworks[0];
        int tmpListNumber = 0;
            for (IAtom tmpAtom : tmpSFramework.atoms()){
                if(tmpAtom.getFlag(CDKConstants.DUMMY_POINTER)){
                    double tmpMinDist = 1000;//Normal distances lower by a factor of 100
                    //Oxygen with the smallest distance to a flagged carbon is selected
                    for(int tmpONumber: addAtomList){
                        double tmpDist = testMol.getAtom(tmpONumber).getPoint2d().distance(tmpAtom.getPoint2d());
                        if(tmpDist<tmpMinDist){
                            tmpMinDist = tmpDist;
                            tmpListNumber = tmpONumber;
                        }
                    }
                    tmpSFramework.addAtom(testMol.getAtom(tmpListNumber));//Oxygen with the smallest distance to a flagged carbon is added
                    tmpSFramework.addBond(tmpAtom.getIndex(), tmpSFramework.getAtomNumber(testMol.getAtom(tmpListNumber)), IBond.Order.DOUBLE);//Double bond inserted between carbon and oxygen

                    //Alternative solution: Oxygen is newly generated and added to the DUMMY_POINTER carbons.
                    //Problem: Positions of the oxygens do not correspond to the original ones.
                    //tmpSFramework.addAtom(new Atom("O"));
                    //tmpSFramework.addBond(tmpSFramework.getLastAtom().getIndex(),tmpAtom.getIndex(), IBond.Order.DOUBLE);
                }
            }

        //Generate picture of the SchuffenhauerFramework
        BufferedImage tmpImgSFrw = generator.withBackgroundColor(Color.PINK).depict(tmpSFramework).toImg();
        ImageIcon tmpIconSFrw = new ImageIcon(tmpImgSFrw);
        JFrame tmpFrameSFrw = new JFrame();
        tmpFrameSFrw.setLayout(new FlowLayout());
        tmpFrameSFrw.setSize(320,370);
        JLabel tmpLblSFrw = new JLabel();
        tmpLblSFrw.setIcon(tmpIconSFrw);
        tmpLblSFrw.setText("Schuffenhauer");
        frame.add(tmpLblSFrw);

        //Show window
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    /**
     * Loads a molfile from the Resources folder and creates all fragments that can be generated with the MurckoFragmenter.
     * The unmodified molecule (original) and all generated fragments are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the file.
     * @param tmpFileName Name of the file to be loaded. File is loaded from the resources folder.
     * @throws CDKException if file cant be read as IAtomContainer
     */
    public void getMurckoFragments(String tmpFileName) throws IOException, CDKException {
        //Load molecule
        InputStream tmpInputStream= murckoFragmenter.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
        MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
        IAtomContainer tmpMolecule = tmpReader.read(new AtomContainer());
        //Generate picture of the original molecule
        DepictionGenerator tmpGenerator = new DepictionGenerator();
        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
        BufferedImage tmpImgOri = tmpGenerator.depict(tmpMolecule).toImg();
        new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Original.png").mkdirs();
        File tmpOutputOri = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Original.png");
        ImageIO.write(tmpImgOri, "png" ,tmpOutputOri);
        //Generate fragments, rings and frameworks
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(false,1);
        tmpMurckoFragmenter.setComputeRingFragments(true);
        tmpMurckoFragmenter.generateFragments(tmpMolecule);
        IAtomContainer[] tmpFragments = tmpMurckoFragmenter.getFragmentsAsContainers();
        IAtomContainer[] tmpRings = tmpMurckoFragmenter.getRingSystemsAsContainers();
        IAtomContainer[] tmpFrameworks = tmpMurckoFragmenter.getFrameworksAsContainers();
        //Generate pictures of the fragments
        int tmpCountFra = 0;
        for(IAtomContainer tmpFragment : tmpFragments) {
            tmpCountFra++;
            BufferedImage tmpImgFra = tmpGenerator.depict(tmpFragment).toImg();
            new File(System.getProperty("user.dir")+"/src/test/scaffoldTestOutput/" + tmpFileName + "/Fragment"+tmpCountFra + ".png").mkdirs();
            File tmpOutputFra = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Fragment" + tmpCountFra+".png");
            ImageIO.write(tmpImgFra, "png" ,tmpOutputFra);
        }
        //Generate pictures of the rings
        int tmpCountRgs = 0;
        for(IAtomContainer tmpRing : tmpRings) {
            tmpCountRgs++;
            BufferedImage tmpImgRgs = tmpGenerator.depict(tmpRing).toImg();
            new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/Ring" + tmpCountRgs + ".png").mkdirs();
            File tmpOutputRgs = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName+"/Ring" + tmpCountRgs + ".png");
            ImageIO.write(tmpImgRgs, "png" ,tmpOutputRgs);
        }
        //Generate pictures of the frameworks
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
