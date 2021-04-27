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

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class murckoFragmenter {
    public static void main(String[] args) {

    }
    /**
     * @param tmpFileName Name of the file to be loaded. File is loaded from the resources folder.
     * @throws CDKException if file cant be read as IAtomContainer
     * @throws IOException if Image cant saved in directory
     */
    public void getSchuffenhauerScaffold(String tmpFileName) throws CDKException, IOException {
        //Load molecule
        InputStream tmpInputStream= murckoFragmenter.class.getClassLoader().getSystemResourceAsStream(tmpFileName+".mol");
        MDLV2000Reader tmpReader = new MDLV2000Reader(tmpInputStream);
        IChemObjectBuilder tmpBuilder = SilentChemObjectBuilder.getInstance();
        IAtomContainer tmpMolecule = tmpReader.read(tmpBuilder.newAtomContainer());
        //IAtomContainer tmpMolecule = tmpReader.read(new AtomContainer()); Dont use new AtomContainer. Makes trouble: https://github.com/cdk/cdk/issues/687

        //Mark each atom with ascending number
        tmpMolecule.getBondCount();
        Integer tmpCounter = 0;
        List<Integer> tmpAddAtomList = new ArrayList<>();//Stores positions of double bounded O in the testMol
        for(IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpCounter++;
            tmpAtom.setProperty("AtomCounter", tmpCounter);
        }
        for(IBond tmpBond: tmpMolecule.bonds()){
            if(tmpBond.getOrder() == IBond.Order.DOUBLE){
                if(tmpBond.getAtom(0).getSymbol() == "C" && tmpBond.getAtom(1).getSymbol() == "O"){
                    //addAtomList.add(Integer.parseInt(tmpBond.getAtom(0).getProperty("AtomCounter")));
                    tmpAddAtomList.add(tmpBond.getAtom(0).getProperty("AtomCounter"));
                }
                if(tmpBond.getAtom(0).getSymbol() == "O" && tmpBond.getAtom(1).getSymbol() == "C"){
                    //addAtomList.add(Integer.parseInt(tmpBond.getAtom(1).getProperty("AtomCounter")));
                    tmpAddAtomList.add(tmpBond.getAtom(1).getProperty("AtomCounter"));
                }
            }
        }
        System.out.println(tmpAddAtomList);
        //Store the number of each carbon with an oxygen double bond
        /*for(IAtom tmpAtom : tmpMolecule.atoms()){
            if(tmpAtom.getSymbol().equals("O")) {
                for (IBond tmpBond : tmpAtom.bonds()) {
                    if (tmpBond.getElectronCount() == 4) {
                        if(tmpBond.getAtom(0).getSymbol() == "C"){
                            //addAtomList.add(Integer.parseInt(tmpBond.getAtom(0).getProperty("AtomCounter")));
                            tmpAddAtomList.add(tmpBond.getAtom(0).getProperty("AtomCounter"));
                        }
                        if(tmpBond.getAtom(1).getSymbol() == "C"){
                            //addAtomList.add(Integer.parseInt(tmpBond.getAtom(1).getProperty("AtomCounter")));
                            tmpAddAtomList.add(tmpBond.getAtom(1).getProperty("AtomCounter"));
                        }
                    }
                }
            }
        }*/
        //Generate the murckoFragment
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        tmpMurckoFragmenter.generateFragments(tmpMolecule);
        IAtomContainer[] tmpFrameworks = tmpMurckoFragmenter.getFrameworksAsContainers();
        //Get the longest murckoFragment
        int tmpAtomCount = 0;
        int tmpFragmentCounter = 0;
        int tmpFragmentNumber = 0;
        for(IAtomContainer tmpFragment : tmpFrameworks){
            if(tmpFragment.getAtomCount()>tmpAtomCount){
                tmpAtomCount = tmpFragment.getAtomCount();
                tmpFragmentNumber = tmpFragmentCounter;
            }
            tmpFragmentCounter++;
        }
        IAtomContainer tmpLongFragment = tmpFrameworks[tmpFragmentNumber];
        //Generate SchuffenhauerFragment
        for(IAtom tmpAtom : tmpLongFragment.atoms()){
            if(tmpAddAtomList.contains(tmpAtom.getProperty("AtomCounter"))){
                tmpLongFragment.addAtom(new Atom("O"));
                for(IAtom tmpNewAtom : tmpLongFragment.atoms()){
                    if(tmpNewAtom.getProperty("AtomCounter") == null){
                        tmpCounter++;
                        tmpNewAtom.setProperty("AtomCounter", tmpCounter);
                        tmpLongFragment.addBond(tmpNewAtom.getIndex(),tmpAtom.getIndex(), IBond.Order.DOUBLE);
                    }
                }
            }
        }
        //Add back hydrogens removed by the MurckoFragmenter
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpLongFragment);
        CDKHydrogenAdder.getInstance(tmpLongFragment.getBuilder()).addImplicitHydrogens(tmpLongFragment);
        //AtomContainerManipulator.convertImplicitToExplicitHydrogens(tmpFrameworks[0]);
        //System.out.println(tmpTestGenerator.create(tmpFrameworks[0]));

        //Generate picture of the SchuffenhauerScaffold
        DepictionGenerator tmpGenerator = new DepictionGenerator();
        tmpGenerator.withSize(600, 600).withTitleColor(Color.BLACK);
        BufferedImage tmpImgSchuff = tmpGenerator.depict(tmpLongFragment).toImg();
        new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png").mkdirs();
        File tmpOutputOri = new File(System.getProperty("user.dir") + "/src/test/scaffoldTestOutput/" + tmpFileName + "/SchuffenhauerScaffold.png");
        ImageIO.write(tmpImgSchuff, "png" ,tmpOutputOri);
    }

    /**
     * Loads a molfile from the Resources folder and creates all fragments that can be generated with the MurckoFragmenter.
     * The unmodified molecule (original) and all generated fragments are saved as images in a subfolder of the scaffoldTestOutput folder.
     * The subfolder has the name of the file.
     * @param tmpFileName Name of the file to be loaded. File is loaded from the resources folder.
     * @throws CDKException if file cant be read as IAtomContainer
     * @throws IOException if Image cant saved in directory
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
