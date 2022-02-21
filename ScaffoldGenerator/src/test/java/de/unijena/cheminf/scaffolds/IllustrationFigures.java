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

import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import javax.imageio.ImageIO;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Class that creates sample images.
 *
 * @author Julian Zander, Jonas Schaub (zanderjulian@gmx.de, jonas.schaub@uni-jena.de)
 * @version 1.1.0.0
 */
public class IllustrationFigures extends ScaffoldGenerator {

    /**
     * Loads Flucloxacillin (CID 21319) out of a SMILES.
     * Generates all types of Scaffolds for this molecule and saves them as pictures.
     * This molecule is also used in Scheme 1 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void scaffoldTypesTest() throws Exception {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1=C(C(=NO1)C2=C(C=CC=C2Cl)F)C(=O)NC3C4N(C3=O)C(C(S4)(C)C)C(=O)O");
        /*Generate picture of the Original molecule*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture of the original*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/Original.png");
        ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
        //Generate SchuffenhauerScaffold
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.SCHUFFENHAUER_SCAFFOLD);
        IAtomContainer tmpSchuffenhauerSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate Murcko Scaffold*/
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.MURCKO_FRAMEWORK);
        IAtomContainer tmpMurckoSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgMurcko = tmpGenerator.depict(tmpMurckoSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputMurcko = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/Murcko.png");
        ImageIO.write(tmpImgMurcko, "png" ,tmpOutputMurcko);
        /*Generate Basic Wire Frame*/
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.BASIC_WIRE_FRAME);
        IAtomContainer tmpBWFSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgBWF = tmpGenerator.depict(tmpBWFSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputBWF = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/BasicWireFrame.png");
        ImageIO.write(tmpImgBWF, "png" ,tmpOutputBWF);
        /*Generate Element Wire Frame*/
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.ELEMENTAL_WIRE_FRAME);
        IAtomContainer tmpEWFSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgEWF = tmpGenerator.depict(tmpEWFSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputEWF = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/ElementWireFrame.png");
        ImageIO.write(tmpImgEWF, "png" ,tmpOutputEWF);
        /*Generate Basic Framework*/
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.BASIC_FRAMEWORK);
        IAtomContainer tmpBFSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgBF = tmpGenerator.depict(tmpBFSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest").mkdirs();
        File tmpOutputBF = new File(System.getProperty("user.dir") + "/Figure/scaffoldTypesTest/BasicFramework.png");
        ImageIO.write(tmpImgBF, "png" ,tmpOutputBF);
    }

    /**
     * Loads Clotiazepam (CID 2811) out of a SMILES.
     * Generates all levels of a schuffenhauer fragmentation out of the molecule and saves them as picture.
     * This molecule is also used in Scheme 18 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void schuffenhauerFragmentationTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        //SMILES to IAtomContainer
        IAtomContainer tmpMoleculeClotiazepam = tmpParser.parseSmiles("CCC1=CC2=C(S1)N(C(=O)CN=C2C3=CC=CC=C3Cl)C");
        /*Generate picture of the Original*/
        BufferedImage tmpImgOriginalClotiazepam = tmpGenerator.depict(tmpMoleculeClotiazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest").mkdirs();
        File tmpOutputOriginalClotiazepam = new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest/ClotiazepamOriginal.png");
        ImageIO.write(tmpImgOriginalClotiazepam, "png" ,tmpOutputOriginalClotiazepam);
        /*Generate picture of the SchuffenhauerScaffold*/
        //Generate SchuffenhauerScaffold
        IAtomContainer tmpSchuffenhauerScaffoldClotiazepam =tmpScaffoldGenerator.getScaffold(tmpMoleculeClotiazepam, true);
        BufferedImage tmpImgScaffoldClotiazepam = tmpGenerator.depict(tmpSchuffenhauerScaffoldClotiazepam).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest").mkdirs();
        File tmpOutputScaffoldClotiazepam = new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest/ClotiazepamStep0.png");
        ImageIO.write(tmpImgScaffoldClotiazepam, "png" ,tmpOutputScaffoldClotiazepam);
        /*Generate picture of the modified molecule*/
        List<IAtomContainer> tmpStep1MolClotiazepam = tmpScaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldClotiazepam);
        BufferedImage tmpImgStep1Clotiazepam = tmpGenerator.depict(tmpStep1MolClotiazepam.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest").mkdirs();
        File tmpOutputStep1Clotiazepam = new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest/ClotiazepamStep1.png");
        ImageIO.write(tmpImgStep1Clotiazepam, "png" ,tmpOutputStep1Clotiazepam);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpStep2MolClotiazepam = tmpScaffoldGenerator.applySchuffenhauerRules(tmpSchuffenhauerScaffoldClotiazepam);
        BufferedImage tmpImgStep2Clotiazepam = tmpGenerator.depict(tmpStep2MolClotiazepam.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest").mkdirs();
        File tmpOutputStep2Clotiazepam = new File(System.getProperty("user.dir") + "/Figure/schuffenhauerFragmentationTest/ClotiazepamStep2.png");
        ImageIO.write(tmpImgStep2Clotiazepam, "png" ,tmpOutputStep2Clotiazepam);
    }

    /**
     * Loads Clotiazepam (CID 2811) out of a SMILES.
     * Generates all outcomes of an enumerative removal out of the molecule and saves them as picture.
     * This molecule is also used in Scheme 18 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void enumerativeRemovalTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCC1=CC2=C(S1)N(C(=O)CN=C2C3=CC=CC=C3Cl)C");
        BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
        new File(System.getProperty("user.dir") + "/Figure/enumerativeRemovalTest").mkdirs();
        File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/Figure/enumerativeRemovalTest/ClotiazepamOriginal.png");
        ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
        //Generate a list of molecules with iteratively removed terminal rings
        List<IAtomContainer> tmpMolecules = tmpScaffoldGenerator.applyEnumerativeRemoval(tmpMolecule);
        int tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpMolecules) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/enumerativeRemovalTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/enumerativeRemovalTest/ClotiazepamStep" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
    }

    /**
     * Loads Zaleplon (CID 5719) out of a SMILES.
     * Generates the Schuffenhauer fragments of this molecule with ruleSevenAppliedSetting turned on/off and saves them as picture.
     * In this case, rule 10 (Smaller Rings are Removed First) leads to the decision when rule 7 is deactivated.
     * This molecule is also used in Scheme 12 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void ruleSevenTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CCN(C1=CC=CC(=C1)C2=CC=NC3=C(C=NN23)C#N)C(=O)C");
        //Generate Schuffenhauer fragments
        List<IAtomContainer> tmpStep1 = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpStep1.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest/Step1.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate picture of the SchuffenhauerRule*/
        List<IAtomContainer> tmpRule = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        BufferedImage tmpImgRule = tmpGenerator.depict(tmpRule.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest").mkdirs();
        File tmpOutputRule = new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest/Step2RuleSevenTrue.png");
        ImageIO.write(tmpImgRule, "png" ,tmpOutputRule);
        /*Generate picture of the SchuffenhauerRule without Rule 7*/
        tmpScaffoldGenerator.setRuleSevenAppliedSetting(false);
        List<IAtomContainer> tmpRuleFalse = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        BufferedImage tmpImgRuleFalse = tmpGenerator.depict(tmpRuleFalse.get(2)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest").mkdirs();
        File tmpOutputRuleFalse = new File(System.getProperty("user.dir") + "/Figure/ruleSevenTest/Step2RuleSevenFalse.png");
        ImageIO.write(tmpImgRuleFalse, "png" ,tmpOutputRuleFalse);
    }

    /**
     * Loads 1,2,3,4,6,7-Hexahydroisoquinoline (CID	89002720) out of a SMILES.
     * Generates the Schuffenhauer fragments of this molecule with RetainOnlyHybridisationsAtAromaticBondsSetting on/off and saves them as picture.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void retainOnlyHybridisationsAtAromaticBondsSettingTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C2=C1CCNCC1=CCC2");
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(true);
        List<IAtomContainer> tmpSchuffenhauerFragments = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        /*Generate picture of the Original*/
        BufferedImage tmpImgFragment = tmpGenerator.depict(tmpSchuffenhauerFragments.get(0)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest").mkdirs();
        File tmpOutputFragment = new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest/Original.png");
        ImageIO.write(tmpImgFragment, "png", tmpOutputFragment);
        /*Generate picture with NonAromaticDBObtainedSetting turned off*/
        BufferedImage tmpImgFragmentFalse = tmpGenerator.depict(tmpSchuffenhauerFragments.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest").mkdirs();
        File tmpOutputFragmentFalse = new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest/DoNotKeepNonAromaticDB.png");
        ImageIO.write(tmpImgFragmentFalse, "png", tmpOutputFragmentFalse);
        /*Generate picture with NonAromaticDBObtainedSetting turned on*/
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(false);
        tmpSchuffenhauerFragments = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        BufferedImage tmpImgFragmentTrue = tmpGenerator.depict(tmpSchuffenhauerFragments.get(1)).toImg();
        /*Save the picture*/
        new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest").mkdirs();
        File tmpOutputFragmentTrue = new File(System.getProperty("user.dir") + "/Figure/retainOnlyHybridisationsAtAromaticBondsSettingTest/KeepNonAromaticDB.png");
        ImageIO.write(tmpImgFragmentTrue, "png", tmpOutputFragmentTrue);
    }

    /**
     * Loads 2-[Cyclohexyl(cyclopentylmethyl)amino]ethanol (CID 62016022) out of a SMILES.
     * Generates the SchuffenhauerScaffold sidechains, linkers and rings of this molecule and saves them as pictures.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void scaffoldFragmentationTest() throws Exception {
        //SMILES to IAtomContainer
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        //IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1CCC(CC1)N(CCO)CC2CCCC2");
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CC1(C(=O)C(=C(O1)C2=CC=C(C=C2)S(=O)(=O)N)C3=CC(=CC=C3)F)C");
        /*Generate picture of the Original molecule*/
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        BufferedImage tmpImgOriginal = tmpGenerator.depict(tmpMolecule).toImg();
        /*Save the picture of the original*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest").mkdirs();
        File tmpOutputOriginal = new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest/Original.png");
        ImageIO.write(tmpImgOriginal, "png" ,tmpOutputOriginal);
        //Generate SchuffenhauerScaffold
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        tmpScaffoldGenerator.setScaffoldModeSetting(ScaffoldGenerator.ScaffoldModeOption.SCHUFFENHAUER_SCAFFOLD);
        IAtomContainer tmpSchuffenhauerSMILES = tmpScaffoldGenerator.getScaffold(tmpMolecule, true);
        /*Generate picture of the SchuffenhauerScaffold*/
        BufferedImage tmpImgSMILES = tmpGenerator.depict(tmpSchuffenhauerSMILES).toImg();
        /*Save the picture of the schuffenhauer scaffold*/
        new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest").mkdirs();
        File tmpOutputSMILES = new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest/Schuffenhauer.png");
        ImageIO.write(tmpImgSMILES, "png" ,tmpOutputSMILES);
        /*Generate the sidechains*/
        List<IAtomContainer> tmpSideChainList = tmpScaffoldGenerator.getSideChains(tmpMolecule, true);
        int tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpSideChainList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest/Sidechain" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Generate the linkers*/
        List<IAtomContainer> tmpLinkerList = tmpScaffoldGenerator.getLinkers(tmpMolecule, true);
        tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpLinkerList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest/Linker" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Generate the rings*/
        List<IAtomContainer> tmpRingList = tmpScaffoldGenerator.getRings(tmpMolecule, true, true);
        tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpRingList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/scaffoldFragmentationTest/Ring" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
    }

    /**
     * Loads Pyrene (CID 31423) out of a SMILES.
     * Generates the Schuffenhauer fragments of this molecule with determineAromaticitySetting turned on/off and saves them as picture.
     * @throws Exception if anything goes wrong
     */
    @Test
    public void determineAromaticitySettingTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit().withAromaticDisplay(); //Vllt überall AromaticityDisplay aktivieren?
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1=CC2=C3C(=C1)C=CC4=CC=CC(=C43)C=C2");
        /*Apply Schuffenhauer with aromaticity detection*/
        /*Apply Schuffenhauer with aromaticity detection and without non-aromatic DBs*/
        tmpScaffoldGenerator.setDetermineAromaticitySetting(true);
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(true);
        List<IAtomContainer> tmpSchuffenhauerWithoutAromaticityList = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        /*Generate the rings*/
        int tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpSchuffenhauerWithoutAromaticityList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest/FragmentAromaticity1OnlyAromaticity1Number" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Apply Schuffenhauer without aromaticity detection and with non-aromatic DBs*/
        tmpScaffoldGenerator.setDetermineAromaticitySetting(true);
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(false);
        tmpSchuffenhauerWithoutAromaticityList = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        /*Generate the rings*/
        tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpSchuffenhauerWithoutAromaticityList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest/FragmentAromaticity1OnlyAromaticity0Number" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Apply Schuffenhauer with aromaticity detection and without non-aromatic DBs*/
        tmpScaffoldGenerator.setDetermineAromaticitySetting(false);
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(true);
        tmpSchuffenhauerWithoutAromaticityList = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        tmpCounter = 0;
        /*Generate the rings*/
        for (IAtomContainer tmpIterative : tmpSchuffenhauerWithoutAromaticityList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest/FragmentAromaticity0OnlyAromaticity1Number" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
        /*Apply Schuffenhauer without aromaticity detection and without non-aromatic DBs*/
        tmpScaffoldGenerator.setDetermineAromaticitySetting(false);
        tmpScaffoldGenerator.setRetainOnlyHybridisationsAtAromaticBondsSetting(false);
        tmpSchuffenhauerWithoutAromaticityList = tmpScaffoldGenerator.applySchuffenhauerRules(tmpMolecule);
        /*Generate the rings*/
        tmpCounter = 0;
        for (IAtomContainer tmpIterative : tmpSchuffenhauerWithoutAromaticityList) {
            /*Generate picture of the molecule*/
            BufferedImage tmpImgRemove = tmpGenerator.depict(tmpIterative).toImg();
            /*Save the picture*/
            new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest" ).mkdirs();
            File tmpOutputRemove = new File(System.getProperty("user.dir") + "/Figure/determineAromaticitySettingTest/FragmentAromaticity0OnlyAromaticity0Number" + tmpCounter + ".png");
            ImageIO.write(tmpImgRemove, "png", tmpOutputRemove);
            tmpCounter++;
        }
    }

    /**
     * Depicts the 20 most frequent scaffolds in DrugBank and COCONUT (for forest and network), imported from CSV files
     * that were created using the performance test command-line application.
     *
     * @throws Exception if anything goes wrong
     * @author Jonas Schaub
     */
    @Test
    public void gridFiguresTest() throws Exception {
        String tmpOutputFolderPath = new File(System.getProperty("user.dir")).getAbsolutePath() + File.separator
                + "Figure" + File.separator + "GridFigures" + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        } else {
            for (File tmpFile : tmpOutputFolderFile.listFiles()) {
                if (!tmpFile.isDirectory()) {
                    tmpFile.delete();
                }
            }
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        HashMap<String, String> tmpFilesMap = new HashMap<>(10);
        tmpFilesMap.put("COCONUT_Forest", "CSV_Origin_Forest_COCONUT.csv");
        tmpFilesMap.put("COCONUT_Network", "CSV_Origin_Network_COCONUT.csv");
        tmpFilesMap.put("DrugBank_Forest", "CSV_Origin_Forest_DrugBank.csv");
        tmpFilesMap.put("DrugBank_Network", "CSV_Origin_Network_DrugBank.csv");
        for (String tmpKey : tmpFilesMap.keySet()) {
            String tmpSpecificOutputFolderPath = tmpOutputFolderPath + tmpKey + File.separator;
            File tmpSpecificOutputFolderFile = new File(tmpSpecificOutputFolderPath);
            if (!tmpSpecificOutputFolderFile.exists()) {
                tmpSpecificOutputFolderFile.mkdirs();
            } else {
                for (File tmpFile : tmpSpecificOutputFolderFile.listFiles()) {
                    if (!tmpFile.isDirectory()) {
                        tmpFile.delete();
                    }
                }
            }
            System.out.println("Output directory: " + tmpSpecificOutputFolderFile);
            File tmpScaffoldsFile = new File(this.getClass().getResource(tmpFilesMap.get(tmpKey)).getFile());
            BufferedReader tmpReader = new BufferedReader(new FileReader(tmpScaffoldsFile));
            //assuming that there are no doublets among the scaffold SMILES
            HashMap<String, Integer> tmpSMILESToFrequencyMap = new HashMap<>(420000);
            String tmpLine = tmpReader.readLine(); //skip header;
            while ((tmpLine = tmpReader.readLine()) != null) {
                String[] tmpLineSplit = tmpLine.split(",");
                tmpSMILESToFrequencyMap.put(tmpLineSplit[0],Integer.parseInt(tmpLineSplit[1]));
            }
            List <Map.Entry<String, Integer>> tmpEntriesList = new LinkedList<>(tmpSMILESToFrequencyMap.entrySet());
            //map sorting based on: https://www.programiz.com/java-programming/examples/sort-map-values
            // call the sort() method of Collections
            Collections.sort(tmpEntriesList, (l1, l2) -> (-1)*l1.getValue().compareTo(l2.getValue()));
            // create a new map
            LinkedHashMap<String, Integer> tmpSortedMap = new LinkedHashMap();
            // get entry from list to the map
            for (Map.Entry<String, Integer> entry : tmpEntriesList) {
                tmpSortedMap.put(entry.getKey(), entry.getValue());
            }
            int i = 0;
            for (Map.Entry tmpEntry : tmpSortedMap.entrySet()) {
                int tmpFrequency = (Integer) tmpEntry.getValue();
                String tmpSMILESCode = (String) tmpEntry.getKey();
                IAtomContainer tmpScaffold = tmpSmiPar.parseSmiles(tmpSMILESCode);
                /*Generate picture of the molecule*/
                BufferedImage tmpImg = tmpDepictionGenerator.withSize(300,350).withZoom(2).depict(tmpScaffold).toImg();
                Graphics2D tmpGraphics2d = tmpImg.createGraphics();
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_ALPHA_INTERPOLATION, RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_COLOR_RENDERING, RenderingHints.VALUE_COLOR_RENDER_QUALITY);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_DITHERING, RenderingHints.VALUE_DITHER_ENABLE);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
                tmpGraphics2d.setRenderingHint(RenderingHints.KEY_STROKE_CONTROL, RenderingHints.VALUE_STROKE_PURE);
                tmpGraphics2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON));
                tmpGraphics2d.drawImage(tmpImg, 0, 0,null);
                tmpGraphics2d.setColor(Color.BLACK);
                tmpGraphics2d.setFont(new Font("Calibri", Font.PLAIN, 40));
                FontMetrics tmpFontMetric = tmpGraphics2d.getFontMetrics();
                int tmpTextWidth = tmpFontMetric.stringWidth(String.valueOf(tmpFrequency));
                tmpGraphics2d.drawString(String.valueOf(tmpFrequency), (tmpImg.getWidth() / 2) - tmpTextWidth / 2, tmpImg.getHeight()-20);
                tmpGraphics2d.dispose();
                /*Save the picture*/
                File tmpOutput = new File(tmpSpecificOutputFolderPath + tmpKey + i + ".png");
                ImageIO.write(tmpImg, "png", tmpOutput);
                i++;
                if (i == 20) {
                    break;
                }
            }
        }

    }

    /**
     * Loads Sertraline (CID 68617) out of a SMILES.
     * Generates the Schuffenhauer tree of this molecule and displays it with GraphStream.
     * This molecule is also used in Scheme 15 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void graphStreamTreeTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl");
        ScaffoldTree tmpScaffoldTree = tmpScaffoldGenerator.generateSchuffenhauerTree(tmpMolecule);
        /*Display the Tree*/
        GraphStreamUtility.displayWithGraphStream(tmpScaffoldTree, true);
    }

    /**
     * Loads Sertraline (CID 68617) out of a SMILES.
     * Generates the enumerative network of this molecule and displays it with GraphStream.
     * This molecule is also used in Scheme 15 in the <a href="https://doi.org/10.1021/ci600338x">"The scaffold Tree"</a> paper.
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void graphStreamNetworkTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("CNC1CCC(C2=CC=CC=C12)C3=CC(=C(C=C3)Cl)Cl");
        ScaffoldNetwork tmpScaffoldNetwork = tmpScaffoldGenerator.generateScaffoldNetwork(tmpMolecule);
        /*Display the Tree*/
        GraphStreamUtility.displayWithGraphStream(tmpScaffoldNetwork, true);
    }

    /**
     * Loads four molecules out of a SMILES.
     * Generates schuffenhauer trees of these molecules, merges them to one tree and displays the tree with GraphStream.
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void graphStreamTreeMergeTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("C1CCC2C(C1)C3=CN=CN=C3S2");
        IAtomContainer tmpMolecule1 = tmpParser.parseSmiles("C1CC2=C(C1)SC3=C2C(=NC=N3)Cl");
        IAtomContainer tmpMolecule2 = tmpParser.parseSmiles("CC1=C(SC=[N+]1CC2=CN=C(N=C2N)C)CCO.[Cl-]"); //6042
        IAtomContainer tmpMolecule3 = tmpParser.parseSmiles("COCCOC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C)OCCOC"); //176870
        List<IAtomContainer> tmpMoleculeList = new ArrayList<>();
        tmpMoleculeList.add(tmpMolecule);
        tmpMoleculeList.add(tmpMolecule1);
        tmpMoleculeList.add(tmpMolecule2);
        tmpMoleculeList.add(tmpMolecule3);
        List<ScaffoldTree> tmpScaffoldTreeList = tmpScaffoldGenerator.generateSchuffenhauerForest(tmpMoleculeList);
        /*Display the Tree*/
        GraphStreamUtility.displayWithGraphStream(tmpScaffoldTreeList.get(0), true);
    }

    /**
     * Loads two molecules out of a SMILES.
     * Generates networks of these molecules, merges them to one network and displays the network with GraphStream.
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void graphStreamNetworkMergeTest() throws Exception {
        SmilesParser tmpParser  = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        ScaffoldGenerator tmpScaffoldGenerator = new ScaffoldGenerator();
        //SMILES to IAtomContainer
        IAtomContainer tmpMolecule = tmpParser.parseSmiles("ClC2NC1SCNN1N2");
        IAtomContainer tmpMolecule1 = tmpParser.parseSmiles("c2ccc(C1NCNN1)cc2");
        List<IAtomContainer> tmpMoleculeList = new ArrayList<>();
        tmpMoleculeList.add(tmpMolecule);
        tmpMoleculeList.add(tmpMolecule1);
        ScaffoldNetwork tmpScaffoldTree = tmpScaffoldGenerator.generateScaffoldNetwork(tmpMoleculeList);
        /*Display the Tree*/
        GraphStreamUtility.displayWithGraphStream(tmpScaffoldTree, true);
    }
}
