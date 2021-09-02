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

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.*;


public class ScaffoldGenerator {

    //<editor-fold desc="Enum">
    /**
     * Enum with which the type of scaffolds to be generated can be set. These scaffolds are then used for the rest of the processing.
     */
    public enum ScaffoldModeOption {

        /**
         * Schuffenhauer scaffolds are generated. Based on the "The Scaffold Tree" Paper by Schuffenhauer et al. 2006.
         * The terminal side chains of the molecule are removed, but any atoms non-single bonded to linkers or rings are retained.
         */
        SCHUFFENHAUER_FRAMEWORK(),

        /**
         * Murcko scaffolds are generated. Based on the "The Properties of Known Drugs. 1. Molecular Frameworks" Paper by Bemis and Murcko 1996.
         * All terminal side chains are removed and only linkers and rings are retained.
         */
        MURCKO_FRAGMENT(),
        
        /**
         * Basic wire frames are generated based on the "Molecular Anatomy: a new multi‑dimensional hierarchical scaffold analysis tool"
         * Paper by Beccari et al. 2021.
         * It is a very abstract form of representation.
         * All side chains are removed, all bonds are converted into single bonds and all atoms are converted into carbons.
         */
        BECCARI_BASIC_WIRE_FRAME(),

        /**
         * All side chains are removed and multiple bonds are converted to single bonds, but the atomic elements remain.
         */
        ELEMENTAL_WIRE_FRAME(),

        /**
         * Basic Frameworks are generated based on the "Molecular Anatomy: a new multi‑dimensional hierarchical scaffold analysis tool"
         * Paper by Beccari et al. 2021.
         * All side chains are removed and all atoms are converted into carbons. The order of the remaining bonds is not changed.
         */
        BECCARI_BASIC_FRAMEWORK();
    }
    //</editor-fold>

    //<editor-fold desc="Public static final constants">
    /**
     * Property of the atoms according to which they are counted and identified.
     */
    public static final String SCAFFOLD_ATOM_COUNTER_PROPERTY = "SCAFFOLD_ATOM_COUNTER_PROPERTY";

    /**
     * Cycle finder used to detect rings.
     */
    public static final CycleFinder CYCLE_FINDER = Cycles.relevant();

    /**
     * Property is true if the backup cycle finder is to be used instead of the normal cycle finder.
     */
    public static final String CYCLE_FINDER_BACKUP_PROPERTY = "CYCLE_FINDER_BACKUP_PROPERTY";

    /**
     * Backup cycle finder used to detect rings.
     * The relevant cycle finder has problems with a few molecules and also finds too many rings in some molecules.
     * Therefore, the mcb is used in these cases.
     */
    public static final CycleFinder CYCLE_FINDER_BACKUP = Cycles.mcb();

    /**
     * Default setting for whether the aromaticity should be determined.
     * By default, the aromaticity is determined.
     */
    public static final boolean DETERMINE_AROMATICITY_SETTING_DEFAULT = true;

    /**
     * Default setting for which aromaticity model should be used.
     * By default, Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet()) is used.
     */
    public static final Aromaticity AROMATICITY_MODEL_SETTING_DEFAULT = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());

    /**
     * Default setting for whether rule 7 should be applied. By default, rule 7 is applied.
     */
    public static final boolean RULE_SEVEN_APPLIED_SETTING_DEFAULT = true;

    /**
     * Default setting for whether hybridisation should only be retained for aromatic bonds.
     * By default, the hybridisation of all bonds is retained.
     */
    public static final boolean RETAIN_ONLY_HYBRIDISATIONS_AT_AROMATIC_BONDS_SETTING_DEFAULT = true;

    /**
     * Default setting for which scaffold mode should be used.
     * By default, ScaffoldModeOption.SCHUFFENHAUER_FRAMEWORK is used.
     */
    public static final ScaffoldModeOption SCAFFOLD_MODE_OPTION_DEFAULT = ScaffoldModeOption.SCHUFFENHAUER_FRAMEWORK;
    //</editor-fold>

    //<editor-fold desc="Private variables">
    /**
     * MurckoFragmenter with the default settings for this class.
     */
    private MurckoFragmenter murckoFragmenter = new MurckoFragmenter(true, 1);

    /**
     * Specifies whether the aromaticity is to be taken into account.
     */
    private boolean determineAromaticitySetting;

    /**
     * Aromaticity model used to determine the aromaticity of the molecules.
     */
    private Aromaticity aromaticityModelSetting;

    /**
     * Indicates whether rule 7 is executed.
     * It can be useful to turn off rule 7 explicitly,
     * as it is only relevant for a relatively small number of molecules, but it increases the computing time.
     */
    private boolean ruleSevenAppliedSetting;

    /**
     * Currently used ScaffoldMode.
     */
    private ScaffoldModeOption scaffoldModeSetting;

    /**
     * With this setting, only the hybridisation of aromatic atoms can be obtained.
     */
    private boolean retainOnlyHybridisationsAtAromaticBondsSetting;
    //</editor-fold>

    //<editor-fold desc="Constructors">
    /**
     * The only constructor of this class. Sets all settings to their default values.
     */
    public ScaffoldGenerator() {
        this.restoreDefaultSettings();
        this.murckoFragmenter.setComputeRingFragments(false);
    }
    //</editor-fold>

    //<editor-fold desc="Public properties">
    //<editor-fold desc="get/are/is">
    /**
     * Specifies whether the aromaticity is to be taken into account.
     * @return true if the aromaticity is determined
     */
    public boolean isAromaticityDetermined() {
        return this.determineAromaticitySetting;
    }

    /**
     * Returns the currently applied Aromaticity model.
     * This consists of the CycleFinder and the ElectronDonation Model.
     * @return the Aromaticity model
     */
    public Aromaticity getAromaticityModel() {
        return this.aromaticityModelSetting;
    }

    /**
     * Indicates whether rule 7 is executed.
     * It can be useful to turn off rule 7 explicitly,
     * as it is only relevant for a relatively small number of molecules, but it increases the computing time.
     * @return true if Rule 7 is applied
     */
    public boolean isRuleSevenApplied() {
        return this.ruleSevenAppliedSetting;
    }

    /**
     * Returns the currently applied scaffold mode.
     * @return the now used scaffold mode
     */
    public ScaffoldModeOption getScaffoldModeSetting() {
        return this.scaffoldModeSetting;
    }

    /**
     * With this setting, only the hybridisation of aromatic atoms can be obtained.
     * @return true if only the hybridisation of aromatic atoms is obtained
     */
    public boolean areOnlyHybridisationsAtAromaticBondsRetained() {
        return this.retainOnlyHybridisationsAtAromaticBondsSetting;
    }

    //</editor-fold>

    //<editor-fold desc="set/clear">
    /**
     * Sets the option to not determine the aromaticity.
     * @param anIsAromaticitySet if true the aromaticity is determined
     */
    public void setDetermineAromaticitySetting(boolean anIsAromaticitySet) {
        this.determineAromaticitySetting = anIsAromaticitySet;
    }

    /**
     * Sets the applied aromaticity model. This consists of the CycleFinder and the ElectronDonation Model.
     * Must not be null. However, the aromaticity model is also not used if {@link ScaffoldGenerator#determineAromaticitySetting} == false.
     * @param anAromaticity the new Aromaticity model
     */
    public void setAromaticityModelSetting(Aromaticity anAromaticity) {
        Objects.requireNonNull(anAromaticity, "Given aromaticity model must not be null. " +
                "The aromaticity detection can instead be deactivated via setDetermineAromaticitySetting(false).");
        this.aromaticityModelSetting = anAromaticity;
    }

    /**
     * Sets the option to skip rule 7.
     * It can be useful to turn off rule 7 explicitly,
     * as it is only relevant for a relatively small number of molecules, but it increases the computing time.
     * @param anIsRuleSevenApplied if true rule 7 is applied
     */
    public void setRuleSevenAppliedSetting(boolean anIsRuleSevenApplied) {
        this.ruleSevenAppliedSetting = anIsRuleSevenApplied;
    }

    /**
     * Sets the now used scaffold mode.
     * @param anScaffoldMode the scaffold mode to use
     */
    public void setScaffoldModeSetting(ScaffoldModeOption anScaffoldMode) {
        Objects.requireNonNull(anScaffoldMode, "Given scaffold mode is null");
        this.scaffoldModeSetting = anScaffoldMode;

    }

    /**
     * Sets the setting that, only the hybridisation of aromatic atoms is obtained.
     * @param anIsOnlyHybridisationsAtAromaticBondsRetained true if only the hybridisation of aromatic atoms is obtained.
     */
    public void setRetainOnlyHybridisationsAtAromaticBondsSetting(boolean anIsOnlyHybridisationsAtAromaticBondsRetained) {
        this.retainOnlyHybridisationsAtAromaticBondsSetting = anIsOnlyHybridisationsAtAromaticBondsRetained;
    }
    /**
     * All settings are set to their default values. Automatically executed by the constructor.
     */
    public void restoreDefaultSettings() {
        this.setDetermineAromaticitySetting(this.DETERMINE_AROMATICITY_SETTING_DEFAULT);
        this.setAromaticityModelSetting(this.AROMATICITY_MODEL_SETTING_DEFAULT);
        this.setRuleSevenAppliedSetting(this.RULE_SEVEN_APPLIED_SETTING_DEFAULT);
        this.setRetainOnlyHybridisationsAtAromaticBondsSetting(this.RETAIN_ONLY_HYBRIDISATIONS_AT_AROMATIC_BONDS_SETTING_DEFAULT);
        this.setScaffoldModeSetting(this.SCAFFOLD_MODE_OPTION_DEFAULT);
    }
    //</editor-fold>
    //</editor-fold>

    //<editor-fold desc="Public methods">
    //<editor-fold desc="Fundamental methods">
    /**
     * Generates the selected fragment type for the entered molecule and returns it. You can choose from the types available in ScaffoldModeOption.
     * Depending on the internal settings via {@link ScaffoldGenerator#aromaticityModelSetting},
     * a specific aromaticity model is applied to determine the aromaticity of the individual atoms of the fragment.
     * {@link ScaffoldGenerator#determineAromaticitySetting} allows you to determine whether the aromaticity is to be determined.
     * All stereochemistry information is deleted.
     * The {@link ScaffoldGenerator#getScaffoldInternal(IAtomContainer, boolean, Aromaticity)} method is called internally.
     * @param aMolecule molecule whose scaffold is produced.
     * @return scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a scaffold of the used type.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present or problem with aromaticity.apply()
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public IAtomContainer getScaffold(IAtomContainer aMolecule) throws CDKException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        IAtomContainer tmpMolecule =  this.getScaffoldInternal(aMolecule, this.determineAromaticitySetting, this.aromaticityModelSetting);
        return tmpMolecule;
    }

    /**
     * Generates a set of rings depending on the CycleFinder selected by {@link ScaffoldGenerator#getCycleFinder(IAtomContainer)}.
     * Can optional add non-single bounded atoms to the rings and returns them.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for aMolecule.
     * @param aMolecule molecule whose rings are produced.
     * @param anIsKeepingNonSingleBonds if true, non-single bonded atoms are retained on the ring.
     * @return rings of the inserted molecule.
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present or problem with aromaticity.apply()
     */
    public List<IAtomContainer> getRings(IAtomContainer aMolecule, boolean anIsKeepingNonSingleBonds) throws CloneNotSupportedException, CDKException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        for(IAtom tmpAtom : aMolecule.atoms()) {
            Objects.requireNonNull(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY),
                    "SCAFFOLD_ATOM_COUNTER_PROPERTY must be set. Use getScaffold first");
        }
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        /*Generate cycles*/
        Cycles tmpNewCycles = this.getCycleFinder(tmpClonedMolecule).find(tmpClonedMolecule);
        IRingSet tmpRingSet = tmpNewCycles.toRingSet();
        List<IAtomContainer> tmpCycles = new ArrayList<>(tmpNewCycles.numberOfCycles());
        int tmpCycleNumber = tmpNewCycles.numberOfCycles();
        //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        HashSet<IBond> tmpAddBondSet = new HashSet<>((tmpClonedMolecule.getAtomCount() / 2), 1);
        /*Store non sinlge bonded atoms*/
        if(anIsKeepingNonSingleBonds == true) { //Only needed if non single bonded atoms are retained
            /*Generate the murckoFragment*/
            IAtomContainer tmpMurckoFragment = this.murckoFragmenter.scaffold(tmpClonedMolecule);
            /*Store the number of each Atom of the murckoFragment*/
            HashSet<Integer> tmpMurckoAtomNumbers = new HashSet<>(tmpClonedMolecule.getAtomCount(), 1);
            for (IAtom tmpMurckoAtom : tmpMurckoFragment.atoms()) {
                tmpMurckoAtomNumbers.add(tmpMurckoAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
            /*Store the number of each Atom that is not single bonded and the respective bond*/
            for (IBond tmpBond : tmpClonedMolecule.bonds()) {
                if (tmpBond.getOrder() != IBond.Order.SINGLE) {//Consider non-single bonds
                    //If both atoms of the bond are in the Murcko fragment, they are taken over anyway
                    if (tmpMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
                            && tmpMurckoAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        continue;
                    }
                    //The binding has not yet been added to the list
                    if(!tmpAddBondSet.contains(tmpBond)) {
                        /*Add the bond*/
                        tmpAddBondSet.add(tmpBond);
                    }
                }
            }
        }
        /*Add Cycles*/
        for(int tmpCount = 0; tmpCount < tmpCycleNumber; tmpCount++) { //Go through all generated rings
            IAtomContainer tmpCycle = tmpRingSet.getAtomContainer(tmpCount); //Store rings as AtomContainer
            if(anIsKeepingNonSingleBonds == true) {
                /*Add the missing atom and the respective bond*/
                HashMap<Integer, IAtom> tmpMurckoAtomMap = new HashMap<>(tmpCycle.getAtomCount(), 1);
                for(IAtom tmpAtom : tmpCycle.atoms()) {
                    /*Save the properties of the murcko fragment*/
                    int tmpAtomPropertyNumber = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    tmpMurckoAtomMap.put(tmpAtomPropertyNumber, tmpAtom);
                }
                for(IBond tmpBond : tmpAddBondSet) { //Go thought all saved bonds
                    /*If both atoms of the bond are contained in the murcko fragment, this bond does not need to be added anymore*/
                    if(tmpMurckoAtomMap.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) &&
                            tmpMurckoAtomMap.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        continue; //Skip this bond
                    }
                    /*Atom 1 of the bond is in the Murcko fragment*/
                    if(tmpMurckoAtomMap.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(0).clone();
                        tmpCycle.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoAtomMap.get(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 1);
                        tmpNewBond.setAtom(tmpClonedAtom, 0); //Set the second atom
                        tmpCycle.addBond(tmpNewBond); //Add the whole bond
                        continue; //Next bond
                    }
                    /*Atom 0 of the bond is in the Murcko fragment*/
                    if(tmpMurckoAtomMap.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(1).clone();
                        tmpCycle.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoAtomMap.get(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 0);
                        tmpNewBond.setAtom(tmpClonedAtom, 1); //Set the second atom
                        tmpCycle.addBond(tmpNewBond); //Add the whole bond
                    }
                }
            }
            tmpCycles.add(tmpCycle); //Add rings to list
        }
        return tmpCycles;
    }
    //</editor-fold>

    //<editor-fold desc="Advanced methods">
    /**
     * Iteratively removes the terminal rings. All resulting Scaffold are returned. Duplicates are not permitted.
     * The Scaffold of the entire entered molecule is stored first in the list.
     * @param aMolecule Molecule to be disassembled.
     * @return List with all resulting Scaffold.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> getIterativeRemoval(IAtomContainer aMolecule) throws CDKException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        IAtomContainer tmpScaffoldOriginal = this.getScaffoldInternal(aMolecule, true, this.aromaticityModelSetting);
        int tmpRingCount = this.getRings(tmpScaffoldOriginal, true).size();
        List<String> tmpAddedSMILESList = new ArrayList<>(tmpRingCount * 45);
        //List of all fragments already created and size estimated on the basis of an empirical value
        List<IAtomContainer> tmpIterativeRemovalList = new ArrayList<>(tmpRingCount * 45);
        tmpIterativeRemovalList.add(tmpScaffoldOriginal); //Add origin Scaffold
        for(int tmpCounter = 0 ; tmpCounter < tmpIterativeRemovalList.size(); tmpCounter++) {//Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterativeRemovalList.get(tmpCounter); //Take the next molecule from the list
            List<IAtomContainer> tmpAllRingsList = this.getRings(tmpIterMol,true);
            int tmpRingSize = tmpAllRingsList.size();
            for(IAtomContainer tmpRing : tmpAllRingsList) { //Go through all rings
                //Skip molecule if it has less than 2 rings or the ring is not removable
                if(tmpRingSize < 2) {
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing) && this.isRingRemovable(tmpRing, tmpAllRingsList, tmpIterMol)) { //Consider all terminal rings
                    boolean tmpIsInList = false;
                    IAtomContainer tmpRingRemoved = this.getScaffoldInternal(this.removeRing(tmpIterMol, tmpRing), false, null); //Remove next ring
                    String tmpRingRemovedSMILES = tmpGenerator.create(tmpRingRemoved); //Generate unique SMILES
                    if(tmpAddedSMILESList.contains(tmpRingRemovedSMILES)) { //Check if the molecule has already been added to the list
                        tmpIsInList = true;
                    }
                    if(tmpIsInList == false) { //Add the molecule only if it is not already in the list
                        tmpIterativeRemovalList.add(tmpRingRemoved);
                        tmpAddedSMILESList.add(tmpRingRemovedSMILES);
                    }
                }
            }
        }
        return tmpIterativeRemovalList;
    }

    /**
     * Iteratively removes the terminal rings. All resulting Scaffold are saved in a ScaffoldTree.
     * A new level is created with each removal step. Duplicates are permitted.
     * The Scaffold of the entire entered molecule is the root of the tree.
     * @param aMolecule Molecule to be disassembled.
     * @return ScaffoldTree with all resulting Scaffold.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public ScaffoldTree getRemovalTree(IAtomContainer aMolecule) throws CDKException, CloneNotSupportedException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        ScaffoldTree tmpScaffoldTree = new ScaffoldTree();
        IAtomContainer tmpScaffoldOriginal = this.getScaffoldInternal(aMolecule, true, this.aromaticityModelSetting);
        int tmpRingCount = this.getRings(tmpScaffoldOriginal, true).size();
        //List of all fragments already created and size estimated on the basis of an empirical value
        List<IAtomContainer> tmpIterativeRemovalList = new ArrayList<>(tmpRingCount * 45);
        List<TreeNode> tmpAllNodesList = new ArrayList<>(); //List of all TreeNodes
        tmpIterativeRemovalList.add(tmpScaffoldOriginal); //Add origin Scaffold
        TreeNode<IAtomContainer> tmpParentNode = new TreeNode<IAtomContainer>(tmpScaffoldOriginal); //Set origin Scaffold as root
        tmpAllNodesList.add(tmpParentNode);
        int tmpLevelCounter = 0; //Shows which level of the tree we are currently on.
        for(int tmpCounter = 0 ; tmpCounter < tmpIterativeRemovalList.size(); tmpCounter++) { //Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterativeRemovalList.get(tmpCounter); //Take the next molecule from the list
            List<IAtomContainer> tmpRings = this.getRings(tmpIterMol,true);
            int tmpRingSize = tmpRings.size();
            for(IAtomContainer tmpRing : tmpRings) { //Go through all rings
                if(tmpRingSize < 2) { //Skip molecule if it has less than 2 rings
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing) && this.isRingRemovable(tmpRing, tmpRings, tmpIterMol)) { //Consider all terminal rings
                    IAtomContainer tmpRingRemoved = this.getScaffoldInternal(this.removeRing(tmpIterMol, tmpRing), false, null); //Remove next ring
                    tmpAllNodesList.add(tmpAllNodesList.get(tmpLevelCounter).addChild(tmpRingRemoved)); //Add next node to current Level
                    tmpIterativeRemovalList.add(tmpRingRemoved); // The molecule added to the tree is added to the list
                }
            }
            tmpLevelCounter++; //Increases when a level is completed
        }
        TreeNodeIter<IAtomContainer> tmpNodeIter = new TreeNodeIter<>(tmpParentNode);
        /*Add generated nodes to the ScaffoldTree */
        while(tmpNodeIter.hasNext()) { // As long as there are still other molecules in the tree
            TreeNode<IAtomContainer> tmpMoleculeNode = tmpNodeIter.next(); // Next molecule in tree
            tmpScaffoldTree.addNode(tmpMoleculeNode);
        }
        return tmpScaffoldTree;
    }

    /**
     * Iteratively removes the rings of the molecule according to specific rules that are queried hierarchically.
     * Based on the rules from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * Rule 7 {@link ScaffoldGenerator#applySchuffenhauerRuleSeven(IAtomContainer, List)} is only applied
     * if {@link ScaffoldGenerator#ruleSevenAppliedSetting} is true
     * and the aromaticity is also redetermined by {@link ScaffoldGenerator#determineAromaticitySetting}.
     * @param aMolecule Molecule that is to be broken down into its fragments
     * @return Fragments of the molecule according to the Schuffenhauer rules
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> applySchuffenhauerRules(IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpScaffold = this.getScaffoldInternal(tmpClonedMolecule, this.determineAromaticitySetting ,this.aromaticityModelSetting);
        /*All molecules with an atom-to-ring ratio of less than 1.0 are assigned the CYCLE_FINDER_BACKUP_PROPERTY = true property,
         since too many rings were probably detected. The fact that a molecule has more rings than atoms seems concerning. That is why this value was chosen.*/
        int tmpRingNumber = this.getRings(tmpScaffold, false).size();
        float tmpRingAtomRatio = (float) tmpScaffold.getAtomCount() / tmpRingNumber;
        if(tmpRingAtomRatio < 1.0 ) {
            /*Change the property of all atoms of the molecule*/
            for(IAtom tmpAtom : tmpClonedMolecule.atoms()) {
                tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
            }
            /*Apply the new Cyclefinder to the molecules*/
            tmpRingNumber = this.getRings(tmpScaffold, false).size();
            tmpScaffold = this.getScaffoldInternal(tmpClonedMolecule, false ,null);
        }
        //List of all generated fragments
        List<IAtomContainer> tmpScaffoldFragments = new ArrayList<>(tmpRingNumber);
        tmpScaffoldFragments.add(tmpScaffold);
        /*Go through all the fragments generated and try to break them down further*/
        for(int tmpCounter = 0 ; tmpCounter < tmpScaffoldFragments.size(); tmpCounter++) {
            List<IAtomContainer> tmpRings = this.getRings(tmpScaffoldFragments.get(tmpCounter), true);
            /*If the fragment has only one ring or no ring, it does not need to be disassembled further*/
            if(tmpRings.size() == 1 || tmpRings.size() == 0) {
                break;
            }
            /*Only the removable terminal rings are further investigated*/
            List<IAtomContainer> tmpRemovableRings = new ArrayList<>(tmpRings.size());
            for (IAtomContainer tmpRing : tmpRings) {
                if (this.isRingTerminal(tmpScaffoldFragments.get(tmpCounter), tmpRing)
                        && this.isRingRemovable(tmpRing, tmpRings, tmpScaffoldFragments.get(tmpCounter))) {
                    tmpRemovableRings.add(tmpRing); //Add the candidate rings
                }
            }
            /*If the fragment has no  candidate ring, it does not need to be disassembled further*/
            if(tmpRemovableRings.size() == 0) {
                break;
            }
            /*Apply rule number one*/
            tmpRemovableRings = this.applySchuffenhauerRuleOne(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number two*/
            tmpRemovableRings = this.applySchuffenhauerRuleTwo(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number three*/
            tmpRemovableRings = this.applySchuffenhauerRuleThree(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number four and five*/
            tmpRemovableRings = this.applySchuffenhauerRuleFourAndFive(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number six*/
            tmpRemovableRings = this.applySchuffenhauerRuleSix(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            //Rule seven is only useful when aromaticity is redetermined
            if(this.ruleSevenAppliedSetting && this.determineAromaticitySetting) {
                int tmpRuleSevenRingNumber = tmpRemovableRings.size();
                /*Apply rule number seven*/
                tmpRemovableRings = this.applySchuffenhauerRuleSeven(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);

                /*Store molecules in which the number of rings to be examined is reduced*/
                //if(tmpRemovableRings.size() < tmpRuleSevenRingNumber) {
                if(false) {
                    System.out.println("COCONUT ID: " + aMolecule.getProperty("coconut_id"));
                    /*Generate control pictures*/
                    DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                    /*Generate and save molecule picture*/
                    aMolecule = AtomContainerManipulator.removeHydrogens(aMolecule);
                    BufferedImage tmpImgMol = tmpGenerator.depict(aMolecule).toImg();
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule7COCONUT/" + aMolecule.getProperty("coconut_id") + ".png").mkdirs();
                    File tmpOutputMol = new File(System.getProperty("user.dir") +  "/scaffoldTestOutput/Rule7COCONUT/" + aMolecule.getProperty("coconut_id") + ".png");
                    try {
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

                if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                    this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                    //After a new fragment has been added, the next one is investigated
                    continue;
                }
            }
            /*Apply rule number eight*/
            tmpRemovableRings = this.applySchuffenhauerRuleEight(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number nine*/
            tmpRemovableRings = this.applySchuffenhauerRuleNine(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number ten*/
            tmpRemovableRings = this.applySchuffenhauerRuleTen(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number eleven*/
            tmpRemovableRings = this.applySchuffenhauerRuleEleven(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number twelve*/
            tmpRemovableRings = this.applySchuffenhauerRuleTwelve(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRule(tmpRemovableRings.get(0), tmpScaffoldFragments);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number thirteen, the tiebreaking rule */
            tmpScaffoldFragments.add(this.applySchuffenhauerRuleThirteen(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings));
        }
        return tmpScaffoldFragments;
    }


    /**
     * Iteratively removes the rings of the molecule according to specific rules that are queried hierarchically.
     * A tree is built from the resulting fragments.
     * Based on the rules from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * Rule 7 {@link ScaffoldGenerator#applySchuffenhauerRuleSeven(IAtomContainer, List)} is only applied
     * if {@link ScaffoldGenerator#ruleSevenAppliedSetting} is true
     * and the aromaticity is also redetermined by {@link ScaffoldGenerator#determineAromaticitySetting}.
     * @param aMolecule Molecule that is to be broken down into its fragments
     * @return A tree consisting of fragments of the molecule according to the Schuffenhauer rules
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public ScaffoldTree applySchuffenhauerRulesTree(IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
        Objects.requireNonNull(aMolecule, "Input molecule must be non null");
        ScaffoldTree tmpScaffoldTree = new ScaffoldTree();List<TreeNode> tmpAllNodesList = new ArrayList<>(); //List of all TreeNodes
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpScaffold = this.getScaffoldInternal(tmpClonedMolecule, this.determineAromaticitySetting ,this.aromaticityModelSetting);
        /*All molecules with an atom-to-ring ratio of less than 1.0 are assigned the CYCLE_FINDER_BACKUP_PROPERTY = true property,
         since too many rings were probably detected. The fact that a molecule has more rings than atoms seems concerning. That is why this value was chosen.*/
        int tmpRingNumber = this.getRings(tmpScaffold, false).size();
        float tmpRingAtomRatio = (float) tmpScaffold.getAtomCount() / tmpRingNumber;
        if(tmpRingAtomRatio < 1.0 ) {
            /*Change the property of all atoms of the molecule*/
            for(IAtom tmpAtom : tmpClonedMolecule.atoms()) {
                tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
            }
            /*Apply the new Cyclefinder to the molecules*/
            tmpRingNumber = this.getRings(tmpScaffold, false).size();
            tmpScaffold = this.getScaffoldInternal(tmpClonedMolecule, false ,null);
        }
        //List of all generated fragments
        List<IAtomContainer> tmpScaffoldFragments = new ArrayList<>(tmpRingNumber);
        tmpScaffoldFragments.add(tmpScaffold);
        TreeNode<IAtomContainer> tmpParentNode = new TreeNode<IAtomContainer>(tmpScaffold); //Set origin Scaffold as root
        tmpAllNodesList.add(tmpParentNode);
        /*Go through all the fragments generated and try to break them down further*/
        for(int tmpCounter = 0 ; tmpCounter < tmpScaffoldFragments.size(); tmpCounter++) {
            List<IAtomContainer> tmpRings = this.getRings(tmpScaffoldFragments.get(tmpCounter), true);
            /*If the fragment has only one ring or no ring, it does not need to be disassembled further*/
            if(tmpRings.size() == 1 || tmpRings.size() == 0) {
                break;
            }
            /*Only the removable terminal rings are further investigated*/
            List<IAtomContainer> tmpRemovableRings = new ArrayList<>(tmpRings.size());
            for (IAtomContainer tmpRing : tmpRings) {
                if (this.isRingTerminal(tmpScaffoldFragments.get(tmpCounter), tmpRing)
                        && this.isRingRemovable(tmpRing, tmpRings, tmpScaffoldFragments.get(tmpCounter))) {
                    tmpRemovableRings.add(tmpRing); //Add the candidate rings
                }
            }
            /*If the fragment has no  candidate ring, it does not need to be disassembled further*/
            if(tmpRemovableRings.size() == 0) {
                break;
            }
            /*Apply rule number one*/
            tmpRemovableRings = this.applySchuffenhauerRuleOne(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number two*/
            tmpRemovableRings = this.applySchuffenhauerRuleTwo(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number three*/
            tmpRemovableRings = this.applySchuffenhauerRuleThree(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number four and five*/
            tmpRemovableRings = this.applySchuffenhauerRuleFourAndFive(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number six*/
            tmpRemovableRings = this.applySchuffenhauerRuleSix(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            //Rule seven is only useful when aromaticity is redetermined
            if(this.ruleSevenAppliedSetting && this.determineAromaticitySetting) {
                int tmpRuleSevenRingNumber = tmpRemovableRings.size();
                /*Apply rule number seven*/
                tmpRemovableRings = this.applySchuffenhauerRuleSeven(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);

                /*Store molecules in which the number of rings to be examined is reduced*/
                if(false) {
                    System.out.println("COCONUT ID: " + aMolecule.getProperty("coconut_id"));
                    /*Generate control pictures*/
                    DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(512,512).withFillToFit();
                    /*Generate and save molecule picture*/
                    aMolecule = AtomContainerManipulator.removeHydrogens(aMolecule);
                    BufferedImage tmpImgMol = tmpGenerator.depict(aMolecule).toImg();
                    new File(System.getProperty("user.dir") + "/scaffoldTestOutput/Rule7COCONUT/" + aMolecule.getProperty("coconut_id") + ".png").mkdirs();
                    File tmpOutputMol = new File(System.getProperty("user.dir") +  "/scaffoldTestOutput/Rule7COCONUT/" + aMolecule.getProperty("coconut_id") + ".png");
                    try {
                        ImageIO.write(tmpImgMol, "png", tmpOutputMol);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }

                if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                    this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                    //After a new fragment has been added, the next one is investigated
                    continue;
                }
            }
            /*Apply rule number eight*/
            tmpRemovableRings = this.applySchuffenhauerRuleEight(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number nine*/
            tmpRemovableRings = this.applySchuffenhauerRuleNine(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number ten*/
            tmpRemovableRings = this.applySchuffenhauerRuleTen(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number eleven*/
            tmpRemovableRings = this.applySchuffenhauerRuleEleven(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number twelve*/
            tmpRemovableRings = this.applySchuffenhauerRuleTwelve(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                this.removeRingForSchuffenhauerRuleTree(tmpRemovableRings.get(0), tmpScaffoldFragments, tmpAllNodesList);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number thirteen, the tiebreaking rule */
            IAtomContainer tmpLastRuleMolecule =this.applySchuffenhauerRuleThirteen(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1), tmpRemovableRings);
            tmpScaffoldFragments.add(tmpLastRuleMolecule);
            //Add the Node to the list of Nodes
            tmpAllNodesList.add(new TreeNode<IAtomContainer>(tmpScaffoldFragments.get(tmpScaffoldFragments.size() - 1))); //Add next node
        }
        IAtomContainer tmpReverseStartMolecule = (IAtomContainer) tmpAllNodesList.get((tmpAllNodesList.size() - 1)).getMolecule();
        /*Set the root for the ScaffoldTree*/
        TreeNode tmpReverseParentNode =  new TreeNode<IAtomContainer>(tmpReverseStartMolecule);
        tmpScaffoldTree.addNode(tmpReverseParentNode);
        /*Build the ScaffoldTree with the smallest fragment as root*/
        for(int i = 1; i < tmpAllNodesList.size(); i++) {
            TreeNode tmpNewNode = tmpAllNodesList.get((tmpAllNodesList.size() - 1) - i);
            IAtomContainer tmpTestMol = (IAtomContainer) tmpNewNode.getMolecule();
            tmpScaffoldTree.getAllNodesOnLevel(i - 1).get(0).addChild(tmpTestMol);
            TreeNode tmpTestNode = (TreeNode) tmpScaffoldTree.getAllNodesOnLevel(i - 1).get(0).getChildren().get(0);
            tmpScaffoldTree.addNode(tmpTestNode);
        }
        return tmpScaffoldTree;
    }

    /**
     * Decomposes the entered molecules into SchuffenhauerScaffolds, creates ScaffoldTrees from them and then assembles these trees if possible.
     * If trees have the same root (smallest fragment), they are joined together so that the same fragments are no longer duplicated.
     * In this way, no fragment created is lost when it is joined together.
     * @param aMoleculeList Molecules to be transferred into list of trees
     * @return List of ScaffoldTrees consisting of the fragments of the entered molecules.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible
     */
    public List<ScaffoldTree> mergeAllTrees (List<IAtomContainer> aMoleculeList) throws CDKException, CloneNotSupportedException {
        Objects.requireNonNull(aMoleculeList, "Input molecule list must be non null");
        /*Prepare the output list*/
        List<ScaffoldTree> tmpOutputForest = new ArrayList<>();
        ScaffoldTree tmpFirstTree = new ScaffoldTree();
        tmpOutputForest.add(tmpFirstTree);
        /*Go through all molecules*/
        for(IAtomContainer tmpMolecule : aMoleculeList) {
            boolean isMoleculeMerged = false;
            ScaffoldTree tmpOldTree = this.applySchuffenhauerRulesTree(tmpMolecule);
            /*Go through each newly created tree*/
            for(ScaffoldTree tmpNewTree : tmpOutputForest) {
                /*When one of the new trees has been joined with one of the old trees, move on to the next molecule*/
                if(tmpNewTree.mergeTree(tmpOldTree)) {
                    isMoleculeMerged = true;
                    break;
                }
            }
            /*If the molecule could not be included in a tree add the tree of the molecule*/
            if(isMoleculeMerged == false) {
                tmpOutputForest.add(tmpOldTree);
            }
        }
        return tmpOutputForest;
    }
    //</editor-fold>
    //</editor-fold>

    //<editor-fold desc="Protected methods">
    //<editor-fold desc="General processing">
    /**
     * Generates the selected fragment type for the entered molecule and returns it. You can choose from the types available in ScaffoldModeOption.
     * Depending on the internal settings via {@link ScaffoldGenerator#aromaticityModelSetting},
     * a specific aromaticity model is applied to determine the aromaticity of the individual atoms of the fragment.
     * {@link ScaffoldGenerator#determineAromaticitySetting} allows you to determine whether the aromaticity is to be determined.
     * All stereochemistry information is deleted.
     * @param aMolecule molecule whose Schuffenhauer scaffold is produced.
     * @param anIsAromaticitySet Indicates whether the aromaticity is to be set.
     * @param anAromaticity anAromaticity Model to be used to determine aromaticity. Can be null if anIsAromaticitySet == false.
     * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present or problem with aromaticity.apply()
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected IAtomContainer getScaffoldInternal(IAtomContainer aMolecule, boolean anIsAromaticitySet, Aromaticity anAromaticity) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        /*Clear the stereo chemistry of the molecule*/
        List<IStereoElement> tmpStereo = new ArrayList<>();
        tmpClonedMolecule.setStereoElements(tmpStereo);
        /*Basic wire frames and element wire frames will be numbered later, as their number will be deleted immediately by anonymization and skeleton*/
        if(ScaffoldModeOption.BECCARI_BASIC_WIRE_FRAME != this.getScaffoldModeSetting() && ScaffoldModeOption.ELEMENTAL_WIRE_FRAME != this.getScaffoldModeSetting()) {
            /*Mark each atom with ascending number*/
            Integer tmpCounter = 0;
            for(IAtom tmpAtom : tmpClonedMolecule.atoms()) {
                tmpAtom.setProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY, tmpCounter);
                tmpCounter++;
            }
        }
        /*Generate the murckoFragment*/
        IAtomContainer tmpMurckoFragment = this.murckoFragmenter.scaffold(tmpClonedMolecule);
        switch (this.scaffoldModeSetting) {
            /*Generate the Murcko scaffold*/
            case MURCKO_FRAGMENT:
                break;
            /*Generate the basic wire frame*/
            case BECCARI_BASIC_WIRE_FRAME:
                tmpMurckoFragment = AtomContainerManipulator.anonymise(tmpMurckoFragment);
                /*Mark each atom with ascending number after anonymization because all properties are removed*/
                Integer tmpCounterBWF = 0;
                for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY, tmpCounterBWF);
                    tmpCounterBWF++;
                }
                break;
            /*Generate the Schuffenhauer scaffold*/
            case SCHUFFENHAUER_FRAMEWORK:
                /*Store the number of each Atom of the murckoFragment*/
                HashSet<Integer> tmpMurckoAtomNumbers = new HashSet<>(tmpClonedMolecule.getAtomCount(), 1);
                for (IAtom tmpMurckoAtom : tmpMurckoFragment.atoms()) {
                    tmpMurckoAtomNumbers.add(tmpMurckoAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                }
                /*Store the number of each Atom that is not single bonded and the respective bond*/
                //HashMap cannot be larger than the total number of atoms. Key = Atom and Val = Bond
                HashSet<IBond> tmpAddBondSet = new HashSet<>((tmpClonedMolecule.getAtomCount() / 2), 1);
                for (IBond tmpBond : tmpClonedMolecule.bonds()) {
                    if (tmpBond.getOrder() != IBond.Order.SINGLE && tmpBond.getOrder() != IBond.Order.UNSET) {//Consider non-single bonds
                        //If both atoms of the bond are in the Murcko fragment, they are taken over anyway
                        if (tmpMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
                                && tmpMurckoAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                            continue;
                        }
                        //The binding has not yet been added to the list
                        if (!tmpAddBondSet.contains(tmpBond)) {
                            /*Add the bond*/
                            tmpAddBondSet.add(tmpBond);
                        }
                    }
                }
                /*Add the missing atom and the respective bond*/
                HashMap<Integer, IAtom> tmpMurckoAtomMap = new HashMap<>(tmpMurckoFragment.getAtomCount(), 1);
                for (IAtom tmpAtom : tmpMurckoFragment.atoms()) {
                    /*Save the properties of the murcko fragment*/
                    int tmpAtomProperty = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    tmpMurckoAtomMap.put(tmpAtomProperty, tmpAtom);
                }
                for (IBond tmpBond : tmpAddBondSet) { //Go thought all saved bonds
                    /*If both atoms of the bond are contained in the murcko fragment, this bond does not need to be added anymore*/
                    if (tmpMurckoAtomMap.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) &&
                            tmpMurckoAtomMap.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        continue; //Skip this bond
                    }
                    /*Atom 1 of the bond is in the Murcko fragment*/
                    if (tmpMurckoAtomMap.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(0).clone();
                        tmpMurckoFragment.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoAtomMap.get(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 1);
                        tmpNewBond.setAtom(tmpClonedAtom, 0); //Set the second atom
                        tmpMurckoFragment.addBond(tmpNewBond); //Add the whole bond
                        continue; //Next bond
                    }
                    /*Atom 0 of the bond is in the Murcko fragment*/
                    if (tmpMurckoAtomMap.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(1).clone();
                        tmpMurckoFragment.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoAtomMap.get(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 0);
                        tmpNewBond.setAtom(tmpClonedAtom, 1); //Set the second atom
                        tmpMurckoFragment.addBond(tmpNewBond); //Add the whole bond
                    }
                }
                break;
            /*Generate the element wire frame*/
            case ELEMENTAL_WIRE_FRAME:
                tmpMurckoFragment = AtomContainerManipulator.skeleton(tmpMurckoFragment);
                /*Mark each atom with ascending number after anonymization because all properties are removed*/
                Integer tmpCounterEWF = 0;
                for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY, tmpCounterEWF);
                    tmpCounterEWF++;
                }
                break;
            /*Generate the basic framework*/
            case BECCARI_BASIC_FRAMEWORK:
                for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
                    if(tmpAtom.getSymbol() != "C") {
                        tmpAtom.setSymbol("C");
                    }
                }
                break;
            }
        /*The Murcko fragmenter class does not adjust the hybridisation when the atoms are removed.
        Therefore, this is deleted and determined again.*/
        /*Add back hydrogens removed by the MurckoFragmenter class*/
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMurckoFragment);
        CDKHydrogenAdder.getInstance(tmpMurckoFragment.getBuilder()).addImplicitHydrogens(tmpMurckoFragment);
        /*Set aromaticity if necessary*/
        if (anIsAromaticitySet) {
            Objects.requireNonNull(anAromaticity, "If anIsAromaticitySet == true, anAromaticity must be non null");
            //Set aromaticity
            anAromaticity.apply(tmpMurckoFragment);
        }
        return tmpMurckoFragment;
    }

    /**
     * Removes the given ring from the total molecule and returns it.
     * Preserves the sp2 hybridisation of a border atom when an aromatic ring is removed.
     * Preserves the hybridisation of all molecules if {@link ScaffoldGenerator#retainOnlyHybridisationsAtAromaticBondsSetting} == true
     * With the removal of heterocycles of size 3 a double bond is inserted if it is directly adjacent to another ring.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for aMolecule/aRing and match.
     * @param aMolecule Molecule whose ring is to be removed.
     * @param aRing Ring to be removed.
     * @return Molecule whose ring has been removed.
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected IAtomContainer removeRing(IAtomContainer aMolecule, IAtomContainer aRing) throws CloneNotSupportedException, CDKException {
        /*Clone original molecules*/
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        IAtomContainer tmpRingClone = aRing.clone();
        boolean tmpIsRingAromatic = true;
        HashSet<Integer> tmpIsNotRing = new HashSet<>(aMolecule.getAtomCount(), 1);
        HashSet<Integer> tmpDoNotRemove = new HashSet<>(aMolecule.getAtomCount(), 1);
        int tmpBoundNumber = 0;
        /*Preparation for insertion of double bonds with removal of aromatic rings*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet<>(tmpMoleculeClone.getAtomCount(), 1);
        /*Store the number of each atom in the molecule*/
        for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
            tmpIsNotRing.add(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        /*Remove all numbers of the ring that is to be removed*/
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            tmpIsNotRing.remove(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        /*Get the number of bonds of the ring to other atoms*/
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                //All atoms of the ring in the original molecule
                if(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) {
                    for(IBond tmpBond : tmpMolAtom.bonds()){
                        //Bond between ring and non ring atom
                        if(tmpIsNotRing.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) ||
                                tmpIsNotRing.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                            tmpBoundNumber++;
                        }
                    }
                }
            }
        }
        /*Add all atoms of rings that are not to be removed to tmpDoNotRemove*/
        //Get all cycles of the molecule
        Cycles tmpCycles = this.getCycleFinder(tmpMoleculeClone).find(tmpMoleculeClone);
        HashSet<Integer> tmpRingPropertySet = new HashSet<>(tmpRingClone.getAtomCount(), 1);
        //Save the properties of the ring
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            tmpRingPropertySet.add(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        for(IAtomContainer tmpCycle : tmpCycles.toRingSet().atomContainers()) {
            boolean tmpIsRingToRemove = true;
            /*Check if it is the ring to be removed*/
            for(IAtom tmpCycleAtom : tmpCycle.atoms()) {
                //If one of the atoms of the ring to be removed is not included, it is not this ring
                if(!tmpRingPropertySet.contains(tmpCycleAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                    tmpIsRingToRemove = false;
                }
            }
            /*If it is not the ring you want to remove, add its atoms to the tmpDoNotRemove list*/
            if(tmpIsRingToRemove == false) {
                for(IAtom tmpCycleAtom : tmpCycle.atoms()) {
                    Integer tmpPropertyNumber = tmpCycleAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    if(!tmpDoNotRemove.contains(tmpPropertyNumber)) {
                        tmpDoNotRemove.add(tmpPropertyNumber);
                    }
                }
            }
        }
        if(tmpBoundNumber < 2) { //Remove all ring atoms, as there are less than two bonds to other atoms
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All atoms of the ring in the original molecule
                    if (tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) {
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms. tmpMoleculeCone.remove() not possible
                        /*Saturate the molecule with hydrogens after removal. Important for Scheme 16*/
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
                        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
                    }
                }
            }
        } else { //Remove only the ring atoms that are not bound to the rest of the molecule
            if(tmpRingClone.getAtomCount() == 3) {
                /* Rings consisting of 3 atoms are specially treated*/
                int tmpNonCCounter = 0;
                IAtom tmpNonCAtom = null;
                /*Count the heteroatoms*/
                for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                    if(tmpRingAtom.getSymbol() != "C") {
                        tmpNonCCounter++;
                        tmpNonCAtom = tmpRingAtom;
                    }
                }
                /*If the ring contains one heteroatom, it is treated specially.*/
                if (tmpNonCCounter == 1) {
                    tmpRingClone.removeAtom(tmpNonCAtom); //remove the heteroatom from the ring
                    IAtom tmpBondAtom1 = null;
                    IAtom tmpRemoveAtom = null;
                    for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) { //go through the whole molecule
                        /* Find the second atom to which the heteroatom was bonded if it was sp3 hybridised*/
                        if((tmpMolAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingClone.getAtom(0).getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) ||
                                tmpMolAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingClone.getAtom(1).getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY)) && tmpBondAtom1 != null
                                && tmpMolAtom.getHybridization() == IAtomType.Hybridization.SP3) {
                            //insert a double bond between the two atoms
                            tmpMoleculeClone.getBond(tmpBondAtom1 , tmpMolAtom).setOrder(IBond.Order.DOUBLE);
                        }
                        /*Find the first atom to which the heteroatom was bonded if it was sp3 hybridised*/
                        if((tmpMolAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingClone.getAtom(0).getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) ||
                                tmpMolAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingClone.getAtom(1).getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY)) && tmpBondAtom1 == null
                                && tmpMolAtom.getHybridization() == IAtomType.Hybridization.SP3) {
                            //Save this atom
                            tmpBondAtom1 = tmpMolAtom;
                        }
                        /*The heteroatom is to be removed*/
                        if(tmpNonCAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpMolAtom.getProperty(SCAFFOLD_ATOM_COUNTER_PROPERTY)) {
                            tmpRemoveAtom = tmpMolAtom;
                        }
                    }
                    //remove the heteroatom
                    tmpMoleculeClone.removeAtom(tmpRemoveAtom);
                }
            }
            /*To test whether the ring is aromatic, exocyclic atoms should not be included*/
            IAtomContainer tmpExocyclicRemovedRing = this.getRings(aRing.clone(), false).get(0);
            tmpIsRingAromatic = this.isAtomContainerAromatic(tmpExocyclicRemovedRing);
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    /*All atoms of the ring in the original molecule that are not bound to the rest of the molecule*/
                    if ((tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)
                            == tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) && !tmpDoNotRemove.contains(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms
                        /*Saturate the molecule with hydrogens after removal*/
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
                        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
                    }
                }
            }
            /*Store the number of all atoms from which an aromatic ring has been removed.
             * In these atoms, a double bond was removed without changing the hybridisation from sp2 to sp3.*/
            //Perform calculation only if the ring to be removed is aromatic or if non-aromatic atom hybridisation should also be preserved
            if(tmpIsRingAromatic || !this.retainOnlyHybridisationsAtAromaticBondsSetting) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All Atoms that are sp2 hybridised and in the ring to be removed
                    if (tmpMolAtom.getHybridization() == IAtomType.Hybridization.SP2
                            && tmpRingPropertySet.contains(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        boolean tmpIsSp3 = true;
                        for (IBond tmpBond : tmpMolAtom.bonds()) { //All bonds of the Atom
                            if (tmpBond.getOrder() != IBond.Order.SINGLE) { //If it contains a non-single bond it cannot be sp3
                                tmpIsSp3 = false;
                            }
                        }
                        if (tmpIsSp3) { //If the Atom contains only single bonds, it must be a wanted atom
                            tmpEdgeAtomNumbers.add(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                        }
                    }
                }
            }
            for(IBond tmpBond : tmpMoleculeClone.bonds()) {
                /*If both atoms of a bond were previously part of an aromatic ring, insert a double bond*/
                if(tmpEdgeAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
                        && tmpEdgeAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                    tmpBond.setOrder(IBond.Order.DOUBLE);
                    //Remove the atoms that have already been treated from the set
                    tmpEdgeAtomNumbers.remove(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                    tmpEdgeAtomNumbers.remove(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                }
            }
            /*Increase the number of hydrogens by 1 for all previously untreated edge atoms to compensate for the removed atom.*/
            for(IAtom tmpAtom : tmpMoleculeClone.atoms()) {
                if(tmpEdgeAtomNumbers.contains(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                    tmpAtom.setImplicitHydrogenCount(tmpAtom.getImplicitHydrogenCount() + 1);
                }
            }
        }
        /*Clear hybridisation. The hybridisation must be reset later by percieveAtomTypesAndConfigureAtoms, as the hybridisation is not changed on its own when the atoms are removed.
        sp2 atoms whose double bonds have been removed must be declared as sp3.*/
        for(IAtom tmpAtom : tmpMoleculeClone.atoms()) {
            tmpAtom.setHybridization((IAtomType.Hybridization) CDKConstants.UNSET);
        }
        /*Add back hydrogens removed by the MurckoFragmenter class*/
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Checks whether the tmpRing in the tmpMolecule is terminal. This means whether it can be removed without creating several unconnected parts.
     * Rings that lead to spiro ring systems when removed are also considered non-terminal.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for aMolecule/aRing and match.
     * @param aMolecule Molecule whose ring is to be checked
     * @param aRing Ring to check
     * @return true if the tmpRing is terminal
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected boolean isRingTerminal(IAtomContainer aMolecule, IAtomContainer aRing) throws CloneNotSupportedException {
        /*Clone molecule and ring*/
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpClonedRing = aRing.clone();
        /*Remove ring atoms from original molecule*/
        HashMap<Integer, IAtom> tmpMoleculeCounterMap = new HashMap<>((aMolecule.getAtomCount()), 1);
        for(IAtom tmpMolAtom : tmpClonedMolecule.atoms()) { //Save all atoms of the molecule
            tmpMoleculeCounterMap.put(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY), tmpMolAtom);
        }
        for(IAtom tmpRingAtom : tmpClonedRing.atoms()) { // Go through the ring
            if(tmpMoleculeCounterMap.containsKey(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) { //Is ring atom in molecule
                tmpClonedMolecule.removeAtom(tmpMoleculeCounterMap.get(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))); //Remove them
            }
        }
        /*Check if there is more than one molecule in the IAtomContainer*/
        ConnectivityChecker tmpChecker = new ConnectivityChecker();
        boolean tmpRingIsTerminal = tmpChecker.isConnected(tmpClonedMolecule);
        return tmpRingIsTerminal;
    }

    /**
     * Checks whether rings may be removed.
     * If the ring does not contain atoms that are not present in any other rings, it is not removable.
     * Furthermore, removal is impossible when it is an aromatic ring, that borders two consecutive rings.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for aMolecule/aRings/aRing and match.
     * @param aRing Ring being tested for its removability
     * @param aRings All Rings of the molecule
     * @param aMolecule Whole molecule
     * @return Whether the ring is removable
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected boolean isRingRemovable(IAtomContainer aRing, List<IAtomContainer> aRings, IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpClonedRing = aRing.clone();

        /*---Recognition of rings in which no atom belongs to another ring---*/
        List<IAtomContainer> tmpClonedRings = new ArrayList<>(aRings.size());
        HashSet<Integer> tmpRingsNumbers = new HashSet<>(aRings.size() * tmpClonedRing.getAtomCount(), 1);
        boolean isAnIndependentRing = false;
        /*Store all ring atoms of the whole molecule without the tested ring*/
        for(IAtomContainer tmpRing : aRings) {
            if(tmpRing == aRing) { //Skip the tested ring
                continue;
            }
            tmpClonedRings.add(tmpRing.clone()); //Store the rings
            for(IAtom tmpAtom : tmpRing.atoms()) { //Store the atoms of the rings
                tmpRingsNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
        }
        /*Investigate whether the ring contains atoms that do not occur in any other ring*/
        for (IAtom tmpSingleRingAtom : aRing.atoms()) {
            if(!tmpRingsNumbers.contains(tmpSingleRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))){
                isAnIndependentRing = true;
            }
        }
        /*If the ring does not contain atoms that are not present in any other rings, it is not removable*/
        if(isAnIndependentRing == false) {
            return false;
        }
        /*---If it is an aromatic ring that borders two consecutive rings, its removal is not possible.---*/
        /*Is it an aromatic ring at all*/
        //Remove exocyclic atoms
        IAtomContainer tmpRemovedRing = this.getRings(tmpClonedRing, false).get(0);
        if (!this.isAtomContainerAromatic(tmpRemovedRing)) {
            return true;
        }
        /*Store the number of all ring atoms*/
        HashSet<Integer> tmpMoleculeNumbers = new HashSet<>(tmpClonedMolecule.getAtomCount(), 1);
        for(IAtomContainer tmpRing : tmpClonedRings) {
            for(IAtom tmpRingAtom : tmpRing.atoms()) {
                int tmpAtomNumber = tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                if(!tmpMoleculeNumbers.contains(tmpAtomNumber)) {
                    tmpMoleculeNumbers.add(tmpAtomNumber);
                }
            }
        }
        /*Store all the atoms of the other rings bordering the aromatic ring*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet<>(tmpClonedMolecule.getAtomCount(), 1);
        for(IAtom tmpRingAtom : tmpClonedRing.atoms()) {
            //Skip the atom if it is already in the HashSet
            if(tmpMoleculeNumbers.contains(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                tmpEdgeAtomNumbers.add(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
        }
        /*At least 3 edge atoms are needed to cause a problem*/
        if(tmpEdgeAtomNumbers.size() < 3) {
            return true;
        }
        /*If one of the edge atoms occurs in more than one other ring, it is not possible to remove the ring*/
        for(Integer tmpEdgeAtomNumber : tmpEdgeAtomNumbers) {
            int tmpRingCounter = 0;
            for(IAtomContainer tmpRing : aRings) {
                for(IAtom tmpRingAtom : tmpRing.atoms()) {
                    //If one of the atoms of the ring to be tested matches one of the edge atoms
                    if(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpEdgeAtomNumber) {
                        tmpRingCounter++;
                        if(tmpRingCounter > 1) { //More than one bordering ring
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    /**
     * Removes the selected ring from the last fragment in the list and adds the resulting fragment to this list.
     * Specially designed for {@link ScaffoldGenerator#applySchuffenhauerRules(IAtomContainer)}
     * @param aRing Ring to be removed
     * @param aFragmentList List of all fragments created so far
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected void removeRingForSchuffenhauerRule(IAtomContainer aRing, List<IAtomContainer> aFragmentList) throws CDKException, CloneNotSupportedException {
        //Remove the ring from the fragment currently being treated
        IAtomContainer tmpRingRemoved = this.removeRing(aFragmentList.get(aFragmentList.size() - 1), aRing);
        //Remove the linkers
        IAtomContainer tmpSchuffRingRemoved = this.getScaffoldInternal(tmpRingRemoved, false, null);
        //Add the fragment to the list of fragments
        aFragmentList.add(tmpSchuffRingRemoved);
    }
    /**
     * Removes the selected ring from the last fragment in the list and adds the resulting fragment to this list.
     * In addition, the fragments are added to the node list as nodes.
     * Specially designed for {@link ScaffoldGenerator#applySchuffenhauerRulesTree(IAtomContainer)}
     * @param aRing Ring to be removed
     * @param aFragmentList List of all fragments created so far
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected void removeRingForSchuffenhauerRuleTree(IAtomContainer aRing, List<IAtomContainer> aFragmentList, List<TreeNode> aNodeList) throws CDKException, CloneNotSupportedException {
        //Remove the ring from the fragment currently being treated
        IAtomContainer tmpRingRemoved = this.removeRing(aFragmentList.get(aFragmentList.size() - 1), aRing);
        //Remove the linkers
        IAtomContainer tmpSchuffRingRemoved = this.getScaffoldInternal(tmpRingRemoved, false, null);
        //Add the fragment to the list of fragments
        aFragmentList.add(tmpSchuffRingRemoved);
        //Add the node to the list of nodes
        aNodeList.add(new TreeNode<IAtomContainer>(tmpSchuffRingRemoved));
    }
    /**
     * Selects the correct CycleFinder based on ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY.
     * @param aMolecule Molecule for which a CycleFinder is to be generated
     * @return CycleFinder that matches the properties of the molecule
     */
    protected CycleFinder getCycleFinder(IAtomContainer aMolecule) {
        boolean tmpIsBackupFinderUsed = false;
        /*Check whether the backup Cycle finder should be used*/
        for(IAtom tmpAtom : aMolecule.atoms()) {
            /*If no CycleFinder has been assigned to the Atom yet*/
            if(tmpAtom.getProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY) == null) {
                //The ScaffoldGenerator.CYCLE_FINDER is used by default
                tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, false);
                continue;
            }
            /*If one of the atoms has been assigned to the ScaffoldGenerator.CYCLE_FINDER_BACKUP cycle finder, it is used.*/
            if(tmpAtom.getProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY).equals(true)) {
                tmpIsBackupFinderUsed = true;
                break;
            }
        }
        /*Return the right CycleFinder*/
        if(tmpIsBackupFinderUsed == true) {
            return ScaffoldGenerator.CYCLE_FINDER_BACKUP;
        } else {
            return ScaffoldGenerator.CYCLE_FINDER;
        }
    }

    /**
     * Checks whether the ring of a molecule is in an aromatic fused ring system.
     * These systems cannot easily be further disassembled.
     * @param aRing Ring tested to see if it is in an aromatic fused ring system
     * @param aRings All Rings of the molecule
     * @param aMolecule Whole molecule
     * @return Whether the ring is part of an aromatic fused ring system
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected boolean hasFusedAromaticRings(IAtomContainer aRing, List<IAtomContainer> aRings, IAtomContainer aMolecule) throws CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpClonedRing = aRing.clone();
        List<IAtomContainer> tmpClonedRings = new ArrayList<>(aRings.size());
        List<Integer> tmpRingNumbers = new ArrayList<>(aMolecule.getAtomCount());
        List<Integer> tmpRingsNumbers = new ArrayList<>(aMolecule.getAtomCount());
        /*If the examined ring itself is not aromatic, it is not such a case*/
        for(IAtom tmpAtom : tmpClonedRing.atoms()){
            if(!tmpAtom.isAromatic()) {
                return false;
            }
            //Save the property numbers of the ring
            tmpRingNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        /*Store all ring atoms of the whole molecule without the tested ring*/
        for(IAtomContainer tmpRing : aRings) {
            if(tmpRing == aRing) { //Skip the tested ring
                continue;
            }
            if(this.isAtomContainerAromatic(tmpRing)) { //Skip non aromatic rings
                continue;
            }
            tmpClonedRings.add(tmpRing.clone()); //Store the aromatic rings
            for(IAtom tmpAtom : tmpRing.atoms()) { //Store the atoms of the aromatic rings
                tmpRingsNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
        }
        /*Store all the atoms of the other rings bordering the aromatic ring*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet<>(tmpClonedMolecule.getAtomCount(), 1);
        for(IAtom tmpRingAtom : tmpClonedRing.atoms()) {
            //Skip the atom if it is already in the HashSet
            if(tmpRingsNumbers.contains(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                tmpEdgeAtomNumbers.add(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
        }
        /*At least 3 edge atoms are needed to cause a problem*/
        if(tmpEdgeAtomNumbers.size() < 3) {
            return false;
        }
        /*If one of the edge atoms occurs in more than one other aromatic ring, it is not possible to remove the ring*/
        for(Integer tmpEdgeAtomNumber : tmpEdgeAtomNumbers) {
            int tmpRingCounter = 0;
            for(IAtomContainer tmpRing : tmpClonedRings) {
                for(IAtom tmpRingAtom : tmpRing.atoms()) {
                    //If one of the atoms of the ring to be tested matches one of the edge atoms
                    if(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpEdgeAtomNumber) {
                        tmpRingCounter++;
                        if(tmpRingCounter > 1) { //More than one bordering ring
                            return true; //Ring is in an aromatic fused ring system
                        }
                    }
                }
            }
        }
        return false;
    }

    /**
     * Checks the aromaticity of each atom of the input molecule.
     * If one of the atoms is not aromatic, the whole molecule is not aromatic.
     * @param aMolecule Molecule whose aromaticity is to be determined
     * @return true if the molecule is completely aromatic
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected boolean isAtomContainerAromatic(IAtomContainer aMolecule) throws CloneNotSupportedException {
        /*Check the aromaticity of each atom*/
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        for(IAtom tmpAtom : tmpClonedMolecule.atoms()) {
            //The cycle is not aromatic if one atom is not aromatic
            if(!tmpAtom.isAromatic()) {
                return false;
            }
        }
        return true;
    }
    //</editor-fold>

    //<editor-fold desc="Schuffenhauer rules">
    /**
     * Sort out the rings according to the first Schuffenhauer rule.
     * Based on the first rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Remove Heterocycles of Size 3 First.
     * Therefore, size 3 hetero rings are preferred when available.
     * Only these rings will be returned if present. If none are present, all rings entered will be returned.
     * @param aRings Rings to which the first rule is to be applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleOne(List<IAtomContainer> aRings) {
        int tmpHeteroCyclesCounter = 0; //Number of size 3 heterocycles
        List<IAtomContainer> tmpHeteroRingList = new ArrayList<>(aRings.size()); //Saved size 3 heterocycles
        /*Investigate how many size 3 heterocycles there are*/
        for(IAtomContainer tmpRing : aRings) {
            //All rings of size 3
            if(tmpRing.getAtomCount() == 3 ) {
                int tmpHeteroAtomCounter = 0;
                for(IAtom tmpAtom : tmpRing.atoms()) { //Atoms of the ring
                    if(tmpAtom.getSymbol() != "C"){
                        tmpHeteroAtomCounter++; //Increase if the ring contains a heteroatom
                    }
                }
                if(tmpHeteroAtomCounter == 1) { //If it is a heterocycle
                    tmpHeteroCyclesCounter++;
                    tmpHeteroRingList.add(tmpRing); //Save this ring
                }
            }
        }
        if(tmpHeteroCyclesCounter == 0) { //If there is no heterocycles of size 3
            return aRings; //Unchanged ring list
        } else { //If there are heterocycles of size 3
            return (tmpHeteroRingList); //Only the heterocycles of size 3
        }
    }

    /**
     * Sort out the rings according to the second Schuffenhauer rule.
     * Based on the second rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Do not remove rings with >= 12 Atoms if there are still smaller rings to remove.
     * Therefore, this method prefers smaller rings when macro rings are present.
     * If no macro rings are present, all rings entered will be returned.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected List<IAtomContainer> applySchuffenhauerRuleTwo(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpSmallRings = new ArrayList<>(aRings.size()); //Rings smaller 12
        /*Identify macrocycles and smaller rings*/
        boolean tmpHasRemovableMacroCycle = false;
        for(IAtomContainer tmpRing : aRings) {
            /*To determine the ring size, the exocyclic atoms must be removed*/
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            /*Check whether there are any removable macrocycles at all*/
            if(tmpRemovedExocyclic.toRingSet().getAtomContainer(0).getAtomCount() > 11 ) {
                tmpHasRemovableMacroCycle = true;
            } else { //All removable non macro rings
                tmpSmallRings.add(tmpRing);
            }
        }
        /*Return the unchanged ring list if there are no macrocycles or small rings*/
        if(tmpHasRemovableMacroCycle == false || tmpSmallRings.size() == 0) {
            return aRings;
        }
        /*Return the small rings if there are any*/
        return tmpSmallRings;
    }

    /**
     * Sort out the rings according to the third Schuffenhauer rule.
     * Based on the third rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Choose the Parent Scaffold Having the Smallest Number of Acyclic Linker Bonds.
     * Therefore, linked rings are given priority over fused rings.
     * The rings that are connected to the rest of the molecule via the longest linkers have priority in the removal process.
     * The number of atoms of the linkers is calculated here. The number of linkers is directly dependent on the number of atoms:
     * LinkerBonds = LinkerAtoms - 1
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleThree(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        List<IAtomContainer> tmpRemoveRings = new ArrayList<>(aRings.size()); //Rings with the longest linker
        List<Integer> tmpLinkerSize = new ArrayList<>(aRings.size()); //Linker length of each ring
        /*Generate the murcko fragment, as this removes the multiple bonded atoms at the linkers*/
        Integer tmpMoleculeAtomCount = this.murckoFragmenter.scaffold(tmpClonedMolecule).getAtomCount();
        /*Calculate the linker length of each ring. Negative integers are fused rings*/
        for(IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRemovedRing = this.removeRing(tmpClonedMolecule, tmpRing);
            //Generate the murcko fragment, as this removes the multiple bonded atoms at the linkers
            IAtomContainer tmpRemovedRingMurckoFragment = this.murckoFragmenter.scaffold(tmpRemovedRing);
            //The number of atoms of the removed ring and the molecule from which the ring and the linker were removed are subtracted from the atomic number of the whole molecule
            //This leaves only the atomic number of the linker
            tmpLinkerSize.add(tmpMoleculeAtomCount - (tmpRing.getAtomCount() + tmpRemovedRingMurckoFragment.getAtomCount()));
        }
        //Get the maximum linker size
        Integer tmpMaxList = tmpLinkerSize.stream().mapToInt(v->v).max().orElseThrow(NoSuchElementException::new);
        /*Save the linked rings if available*/
        if(tmpMaxList > -1) { //Is there is a linked ring
            for(int tmpCounter = 0 ; tmpCounter < tmpLinkerSize.size(); tmpCounter++) {
                if(tmpLinkerSize.get(tmpCounter) == tmpMaxList) { //Get the rings with the longest linkers
                    tmpRemoveRings.add(aRings.get(tmpCounter));
                }
            }
            return tmpRemoveRings;
        }
        /*Return the unchanged ring list if there are only fused rings*/
        return aRings;
    }

    /**
     * Sort out the rings according to the fourth and fifth Schuffenhauer rule.
     * Based on the fourth and fifth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The fourth rule says: Retain Bridged Rings, Spiro Rings, and Nonlinear Ring Fusion Patterns with Preference.
     * Therefore, delta is calculated as follows: |nrrb - (nR - 1)|
     * nrrb: number of bonds being a member in more than one ring
     * nR: number of rings
     * The rings with the highest absolute delta are returned
     * The artificial creation of spiro ring systems (see scheme 10 and rule 5) is not possible in our implementation,
     * because such a ring would not be detected as terminal (and only terminal rings are considered for removal)
     * The fifth rule says: Bridged Ring Systems Are Retained with Preference over Spiro Ring Systems.
     * Therefore, the rings with the positive maximum delta are preferred over the rings with the negative one.
     * Through the isRingTerminal() method, a removal that leads to spiro ring systems is not available for selection anyway.
     * For performance reasons, rules four and five are combined. This way, delta only has to be calculated once.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleFourAndFive(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        List<IAtomContainer> tmpRingsReturn = new ArrayList<>(aRings.size()); //Rings that are returned
        List<Integer> tmpDeltaList = new ArrayList<>(aRings.size()); //Delta values of all rings
        List<Integer> tmpDeltaListAbs = new ArrayList<>(aRings.size()); //Absolute Delta values of all rings
        /*Calculate the delta values for all rings*/
        for(IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRingRemoved = this.removeRing(tmpClonedMolecule, tmpRing); //Remove the ring

            //-----Eliminate Cycle Error-----
            Cycles tmpCycles = null;
            Iterable<IAtomContainer> tmpCycleIterable = null;
            /*With a few molecules, an error occurs with the relevant CycleFinder. Then the mcb CycleFinder is automatically used.*/
            try {
                tmpCycles = this.getCycleFinder(tmpRingRemoved).find(tmpRingRemoved); //get cycle number(nR)
                tmpCycleIterable = tmpCycles.toRingSet().atomContainers();
            } catch (NegativeArraySizeException e) {
                /*Save as a property of the atoms that from now on the CYCLE_FINDER_BACKUP is to be used*/
                for(IAtomContainer tmpOriginalRing : aRings) {
                    for(IAtom tmpAtom : tmpOriginalRing.atoms()) {
                        tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                    }
                }
                for(IAtom tmpAtom : aMolecule.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                }
                for(IAtom tmpAtom : tmpRingRemoved.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                }
                tmpCycles = this.getCycleFinder(tmpRingRemoved).find(tmpRingRemoved); //get cycle number(nR)
                tmpCycleIterable = tmpCycles.toRingSet().atomContainers();
            }
            HashSet<IBond> tmpCycleBonds = new HashSet<>(aRings.size());
            int tmpFusedRingBondCounter = 0; // Number of bonds being a member in more than one ring(nrrb)
            /*Count nrrb*/
            for(IAtomContainer tmpCycle : tmpCycleIterable) { //Go through all cycle
                for(IBond tmpBond : tmpCycle.bonds()) { //Go through all bonds of each cycle
                    //If the bond is already included in the list, it occurs in several rings
                    if(tmpCycleBonds.contains(tmpBond)) {
                        //It is assumed that a bond can be in a maximum of two rings
                        tmpFusedRingBondCounter++;
                    }
                    tmpCycleBonds.add(tmpBond);
                }
            }
            //Calculate the delta
            int tmpDelta = (tmpFusedRingBondCounter - (tmpCycles.numberOfCycles() - 1));
            tmpDeltaListAbs.add(Math.abs(tmpDelta));
            tmpDeltaList.add(tmpDelta);
        }
        //Get the maximum delta
        Integer tmpDeltaAbsMax = tmpDeltaListAbs.stream().mapToInt(v->v).max().getAsInt();
        Integer tmpDeltaMax = tmpDeltaList.stream().mapToInt(v->v).max().getAsInt();
        if(tmpDeltaAbsMax > 0) {
            /* In case the delta and the absolute delta are equal we jump to rule five.
            Rule five: if there is a positive maximum delta, only get the maximum positive deltas*/
            if(tmpDeltaAbsMax == tmpDeltaMax) {
                /* Add all rings that have the highest delta to the list*/
                for(int tmpCounter = 0 ; tmpCounter < tmpDeltaList.size(); tmpCounter++) {
                    if(tmpDeltaList.get(tmpCounter) == tmpDeltaAbsMax) {
                        tmpRingsReturn.add(aRings.get(tmpCounter));
                    }
                }
                return tmpRingsReturn; //All rings that have the highest delta
            }
            /*Rule four: Add all rings that have the highest absolute delta to the list*/
            for(int tmpCounter = 0 ; tmpCounter < tmpDeltaListAbs.size(); tmpCounter++) {
                if(tmpDeltaListAbs.get(tmpCounter).equals(tmpDeltaAbsMax)) {
                    tmpRingsReturn.add(aRings.get(tmpCounter));
                }
            }
            return tmpRingsReturn; //All rings that have the highest absolute delta
        }
        /*Return the unchanged ring list*/
        return aRings;
    }

    /**
     * Sort out the rings according to the third Schuffenhauer rule.
     * Based on the sixth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Remove Rings of Sizes 3, 5, and 6 First.
     * Therefore, the exocyclic atoms are removed and the size of the ring is determined.
     * Rings of size 3, 5 and 6 are preferred.
     * If no ring of these sizes is present, the original list is returned.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected List<IAtomContainer> applySchuffenhauerRuleSix(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        /*Size 3, 5 and 6 rings will be added to the list if present*/
        for(IAtomContainer tmpRing : aRings) {
            //To determine the ring size, the exocyclic atoms must be removed
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            IAtomContainer tmpRemovedExoRing = tmpRemovedExocyclic.toRingSet().getAtomContainer(0);
            if(tmpRemovedExoRing.getAtomCount() == 3 || tmpRemovedExoRing.getAtomCount() == 5 || tmpRemovedExoRing.getAtomCount() == 6) {
                tmpReturnRingList.add(tmpRing);
            }
        }
        /*If there are rings of the sizes searched for, they are returned*/
        if(tmpReturnRingList.size() != 0) {
            return tmpReturnRingList;
        }
        /*If there are no rings of the searched sizes, the original rings are returned*/
        return aRings;
    }

    /**
     * Sort out the rings according to the seventh Schuffenhauer rule.
     * Based on the seventh rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: A Fully Aromatic Ring System Must Not Be Dissected in a Way That the Resulting System Is Not Aromatic anymore.
     * It was changed to: The number of aromatic rings should be reduced by a maximum of one, when a ring is removed. Therefore, no additional aromatic rings should be deleted by removing a ring.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleSeven(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        SmilesGenerator tmpSmilesGenerator = new SmilesGenerator((SmiFlavor.Unique));
        List<IAtomContainer> tmpReturnRings = new ArrayList<>(aRings.size());
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        /*Check the number of aromatic rings in the original molecule*/
        int tmpOriginalAromaticRingCounter = 0;
        //Get all cycles without exocyclic atoms
        Cycles tmpOriginalCycles = this.getCycleFinder(tmpClonedMolecule).find(tmpClonedMolecule);
        for(IAtomContainer tmpCycle : tmpOriginalCycles.toRingSet().atomContainers()) {
            /*Count the aromatic rings*/
            if(this.isAtomContainerAromatic(tmpCycle)) {
                tmpOriginalAromaticRingCounter++;
            }
        }
        /*Remove each ring and count the number of remaining aromatic rings*/
        for(IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRemovedRing = this.removeRing(aMolecule, tmpRing);
            tmpRemovedRing = this.getScaffoldInternal(tmpRemovedRing, false, null);
            this.aromaticityModelSetting.apply(tmpRemovedRing);
            /*Check the number of aromatic rings*/
            int tmpRemovedAromaticRingCounter = 0;
            //-----Eliminate Cycle Error-----
            Cycles tmpRemovedCycles = null;
            Iterable<IAtomContainer> tmpCycleIterable = null;
            /*With a few molecules, an error occurs with the relevant CycleFinder. Then the mcb CycleFinder is automatically used.*/
            try {
                tmpRemovedCycles = this.getCycleFinder(tmpRemovedRing).find(tmpRemovedRing);
                tmpCycleIterable = tmpRemovedCycles.toRingSet().atomContainers();
            } catch (Exception e) {
                /*Save as a property of the atoms that from now on the CYCLE_FINDER_BACKUP is to be used*/
                for(IAtomContainer tmpOriginalRing : aRings) {
                    for(IAtom tmpAtom : tmpOriginalRing.atoms()) {
                        tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                    }
                }
                for(IAtom tmpAtom : aMolecule.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                }
                for(IAtom tmpAtom : tmpRemovedRing.atoms()) {
                    tmpAtom.setProperty(ScaffoldGenerator.CYCLE_FINDER_BACKUP_PROPERTY, true);
                }
                tmpRemovedCycles = this.getCycleFinder(tmpRemovedRing).find(tmpRemovedRing);
                tmpCycleIterable = tmpRemovedCycles.toRingSet().atomContainers();
            }
            /*Check the aromaticity of each Cycle*/
            for(IAtomContainer tmpCycle : tmpCycleIterable) {
                /*Count the aromatic rings*/
                if(this.isAtomContainerAromatic(tmpCycle)) {
                    tmpRemovedAromaticRingCounter++;
                }
            }
            /*Only fragments whose aromatic rings have decreased by a maximum of 1 are of interest.*/
            if((tmpOriginalAromaticRingCounter - tmpRemovedAromaticRingCounter) < 2 ) {
                tmpReturnRings.add(tmpRing);
            }
        }
        /*If the number of rings has changed due to this rule, return the changed number*/
        if(tmpReturnRings.size() < aRings.size() && tmpReturnRings.size() != 0) {
            return tmpReturnRings;
        }
        return aRings;
    }

    /**
     * Sort out the rings according to the eighth Schuffenhauer rule.
     * Based on the eighth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Remove Rings with the Least Number of Heteroatoms First
     * Therefore, the exocyclic atoms are removed and the number of cyclic heteroatoms is counted
     * Rings with the smallest number of heteroatoms are preferred
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if all rings have the same size of heteroatoms.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected List<IAtomContainer> applySchuffenhauerRuleEight(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        Integer tmpMinNumberOfHeteroAtoms = null;
        /*Store the rings with the lowest number of cyclic heteroatoms*/
        for(IAtomContainer tmpRing : aRings) {
            //get the cyclic atoms
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            //Number of heteroatoms in the ring
            int tmpNumberOfHeteroAtoms = 0;
            /*Count the heteroatoms*/
            for(IAtom tmpAtom : tmpRemovedExocyclic.toRingSet().getAtomContainer(0).atoms()) {
                if(tmpAtom.getSymbol() != "C") {
                    tmpNumberOfHeteroAtoms++;
                }
            }
            //Set the value of the first ring as starting value
            if(tmpMinNumberOfHeteroAtoms == null) {
                tmpMinNumberOfHeteroAtoms = tmpNumberOfHeteroAtoms;
            }
            /*If the number of heteroatoms matches the number of the least heteroatoms so far, add the ring to the list*/
            if(tmpNumberOfHeteroAtoms == tmpMinNumberOfHeteroAtoms) {
                tmpReturnRingList.add(tmpRing);
                continue;
            }
            /*If the ring has fewer heteroatoms, clear the list and add this ring to it*/
            if(tmpNumberOfHeteroAtoms < tmpMinNumberOfHeteroAtoms) {
                tmpMinNumberOfHeteroAtoms = tmpNumberOfHeteroAtoms;
                tmpReturnRingList.clear();
                tmpReturnRingList.add(tmpRing);
            }
        }
        return tmpReturnRingList;
    }

    /**
     * Sort out the rings according to the ninth Schuffenhauer rule.
     * Based on the ninth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: If the Number of Heteroatoms Is Equal, the Priority of Heteroatoms to Retain is N > O > S.
     * Therefore, the number of cyclic N, O and S of each ring is counted
     * The rings that have the lowest value of heteroatoms according to this rule are selected.
     * If two rings have the same number of N, their amount of O is considered.
     * Heteroatoms that are not N, O or S are ignored.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if all rings have the same size of heteroatoms.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected List<IAtomContainer> applySchuffenhauerRuleNine(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        /*Calculate the maximum number of heteroatoms that can occur*/
        Integer tmpMinNCount = null;
        Integer tmpMinOCount = null;
        Integer tmpMinSCount = null;
        /*Get the rings with the smallest value of heteroatoms*/
        for(IAtomContainer tmpRing : aRings) {
            //Only cyclic heteroatoms count
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            Integer tmpNCounter = 0;
            Integer tmpOCounter = 0;
            Integer tmpSCounter = 0;
            /*Record the composition of the heteroatoms for each ring */
            for(IAtom tmpAtom : tmpRemovedExocyclic.toRingSet().getAtomContainer(0).atoms()) {
                if(tmpAtom.getSymbol() == "N") {
                    tmpNCounter++;
                }
                if(tmpAtom.getSymbol() == "O") {
                    tmpOCounter++;
                }
                if(tmpAtom.getSymbol() == "S") {
                    tmpSCounter++;
                }
            }
            /*Search for the ring with the lowest value of heteroatoms*/
            //Set the values of the first ring as starting values
            if(tmpMinNCount == null) {
                tmpMinNCount = tmpNCounter;
                tmpMinOCount = tmpOCounter;
                tmpMinSCount = tmpSCounter;
            }
            //If the ring contains more N than the previous minimum, it is not eligible for removal
            if(tmpNCounter > tmpMinNCount) {
                continue;
            }
            //If this ring contains less N than the previous minimum, it is considered for removal
            if(tmpNCounter < tmpMinNCount) {
                tmpReturnRingList.clear();
                tmpReturnRingList.add(tmpRing);
                tmpMinNCount = tmpNCounter;
                tmpMinOCount = tmpOCounter;
                tmpMinSCount = tmpSCounter;
                continue;
            }
            /*If this ring has exactly as many N as the previous minimum*/
            //If the ring contains more O than the previous minimum, it is not eligible for removal
            if(tmpOCounter > tmpMinOCount) {
                continue;
            }
            //If this ring contains less O than the previous minimum, it is considered for removal
            if(tmpOCounter < tmpMinOCount) {
                tmpReturnRingList.clear();
                tmpReturnRingList.add(tmpRing);
                tmpMinNCount = tmpNCounter;
                tmpMinOCount = tmpOCounter;
                tmpMinSCount = tmpSCounter;
                continue;
            }
            /*If this ring has exactly as many N and O as the previous minimum*/
            //If the ring contains more S than the previous minimum, it is not eligible for removal
            if(tmpSCounter > tmpMinSCount) {
                continue;
            }
            //If this ring contains less S than the previous minimum, it is considered for removal
            if(tmpSCounter < tmpMinSCount) {
                tmpReturnRingList.clear();
                tmpReturnRingList.add(tmpRing);
                tmpMinNCount = tmpNCounter;
                tmpMinOCount = tmpOCounter;
                tmpMinSCount = tmpSCounter;
                continue;
            }
            /*If this ring has exactly as many N, O and S as the previous minimum*/
            tmpReturnRingList.add(tmpRing);
        }
        return tmpReturnRingList;
    }

    /**
     * Sort out the rings according to the tenth Schuffenhauer rule.
     * Based on the tenth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Smaller Rings are Removed First
     * Exocyclic atoms are not observed
     * Therefore, the exocyclic atoms are removed and the number of cyclic atoms is counted
     * Rings with the smallest number of atoms are preferred
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if all rings have the same size of heteroatoms.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    protected List<IAtomContainer> applySchuffenhauerRuleTen(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        Integer tmpMinimumAtomNumber = null;
        /*Store the rings with the lowest number of atoms*/
        for(IAtomContainer tmpRing : aRings) {
            //Remove the exocyclic atoms
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            IAtomContainer tmpCycle = tmpRemovedExocyclic.toRingSet().getAtomContainer(0);
            int tmpAtomNumber = tmpCycle.getAtomCount();
            /*Set the values of the first ring as starting values*/
            if(tmpMinimumAtomNumber == null) {
                tmpMinimumAtomNumber = tmpAtomNumber;
            }
            /*If the number of atoms matches the number of the least atoms so far, add the ring to the list*/
            if(tmpAtomNumber == tmpMinimumAtomNumber) {
                tmpReturnRingList.add(tmpRing);
                continue;
            }
            /*If the ring has fewer atoms, clear the list and add this ring to it*/
            if(tmpAtomNumber < tmpMinimumAtomNumber) {
                tmpMinimumAtomNumber = tmpAtomNumber;
                tmpReturnRingList.clear();
                tmpReturnRingList.add(tmpRing);
            }
        }
        return tmpReturnRingList;
    }

    /**
     * Sort out the rings according to the eleventh Schuffenhauer rule.
     * Based on the eleventh rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: For Mixed Aromatic/Nonaromatic Ring Systems, Retain Nonaromatic Rings with Priority.
     * Therefore, all rings are tested for aromaticity and the nonaromatic ones are preferably removed.
     * If it is not a mixed system, all rings will be returned.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list,
     * if the molecule is not a mixed ring system.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleEleven(List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        /*Add all fully aromatic rings to the list*/
        for(IAtomContainer tmpRing  : aRings) {
            //Remove the exocyclic atoms
            Cycles tmpRemovedExocyclic = this.getCycleFinder(tmpRing).find(tmpRing);
            IAtomContainer tmpCycle = tmpRemovedExocyclic.toRingSet().getAtomContainer(0);
            /*The ring is only fully aromatic, if all cyclic atoms are aromatic*/
            /*Add aromatic rings to the list*/
            if(this.isAtomContainerAromatic(tmpCycle)) {
                tmpReturnRingList.add(tmpRing);
            }
        }
        //Return aromatic rings if any are present
        if(tmpReturnRingList.size() > 0){
            return tmpReturnRingList;
        }
        //Return all rings if non-aromatic rings are present
        return aRings;
    }

    /**
     * Sort out the rings according to the twelfth Schuffenhauer rule.
     * Based on the twelfth rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Remove Rings First Where the Linker Is Attached
     * to a Ring Heteroatom at Either End of the Linker.
     * Therefore, rings that attached to a linker that have a heteroatom at al. least one end are prioritised.
     *
     * Two cases are treated differently.
     * In the first case, linkers consisting of only one bond are selected.
     * In this case, the ring to be examined is directly linked to the Murcko fragment from which this ring was removed.
     * This bond is found and it is checked whether it contains at least one heteroatom.
     * In the second case, all other linkers are treated. These consist of at least one atom.
     * Here, the linker atoms are filtered out by subtracting the atoms of the ring to be examined and
     * the atoms of the murcko fragment in which this ring was removed from the total molecule.
     * The remaining atoms are the linker atoms. Now it is checked whether their atoms are bound to heteroatoms of the rest of the molecule.
     *
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected List<IAtomContainer> applySchuffenhauerRuleTwelve(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        List<IAtomContainer> tmpRemoveRings = new ArrayList<>(aRings.size()); //Rings with the longest linker
        /*Check for each ring whether it is attached to a linker with a heteroatom at the end*/
        for(IAtomContainer tmpRing : aRings) {
            if(this.isRingAttachedToHeteroatomLinker(tmpClonedMolecule, tmpRing, this.murckoFragmenter)) {
                //If the ring is bound to such a linker add it to the list
                tmpRemoveRings.add(tmpRing);
            }
        }
        /*Return the rings attached to a linker with a heteroatom at the end if available*/
        if(tmpRemoveRings.size() > 0) {
            return tmpRemoveRings;
        }
        /*Return the unchanged ring list if the rule cannot be applied to the rings*/
        return aRings;
    }

    /**
     *
     * The ring to be examined is checked to determine whether it is attached to a linker that has a heteroatom at least one end.
     *
     * Two cases are treated differently.
     * In the first case, linkers consisting of only one bond are selected.
     * In this case, the ring to be examined is directly linked to the Murcko fragment from which this ring was removed.
     * This bond is found and it is checked whether it contains at least one heteroatom.
     * In the second case, all other linkers are treated. These consist of at least one atom.
     * Here, the linker atoms are filtered out by subtracting the atoms of the ring to be examined and
     * the atoms of the Murcko fragment in which this ring was removed from the total molecule.
     * The remaining atoms are the linker atoms. Now it is checked whether their atoms are bound to heteroatoms of the rest of the molecule.
     *
     * Designed for the {@link ScaffoldGenerator#applySchuffenhauerRuleTwelve(IAtomContainer, List)} method.
     * @param aMolecule Molecule from which a ring is to be removed
     * @param aRing rings of the molecule to which the rule is applied
     * @param aMurckoFragmenter MurckoFragmenter with which the molecules are treated
     * @return Whether it is one of the rings sought for
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected boolean isRingAttachedToHeteroatomLinker(IAtomContainer aMolecule, IAtomContainer aRing, MurckoFragmenter aMurckoFragmenter) throws CDKException, CloneNotSupportedException {
        HashSet<Integer> tmpRingPropertyNumbers = new HashSet<>(aRing.getAtomCount(), 1);
        HashSet<Integer> tmpRemovedMurckoAtomNumbers = new HashSet<>(aMolecule.getAtomCount(), 1);
        /*Save all numbers of the ring atoms*/
        for(IAtom tmpAtom : aRing.atoms()) {
            tmpRingPropertyNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        //Remove the examined ring
        IAtomContainer tmpRemovedRing = this.removeRing(aMolecule, aRing);
        //Generate the murcko fragment, as this removes the multiple bonded atoms at the linkers and exocyclic bonds
        IAtomContainer tmpRemovedRingMurckoFragment = aMurckoFragmenter.scaffold(tmpRemovedRing);
        /*Save all numbers of the murcko fragment atoms*/
        for(IAtom tmpAtom : tmpRemovedRingMurckoFragment.atoms()) {
            tmpRemovedMurckoAtomNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        /*Treatment of linkers consisting of only one bond*/
        for(IAtom tmpAtom : aMolecule.atoms()) {
            /*Get all ring atoms*/
            if(tmpRingPropertyNumbers.contains(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                /*Go thought all bonds of the ring atoms*/
                for(IBond tmpBond : tmpAtom.bonds()) {
                    /*Bond that connects ring atom and murcko fragment, so a linker*/
                    if(tmpRemovedMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        /*If the atom of the murcko fragment is a heteroatom, it is one of the rings we are looking for*/
                        if(tmpBond.getAtom(0).getSymbol() != "C") {
                            return true;
                        }
                        /*If the atom of the ring is a heteroatom, it is one of the rings we are looking for*/
                        if(tmpAtom.getSymbol() != "C") {
                            return true;
                        }
                    }
                    /*Bond that connects ring atom and murcko fragment, so a linker.*/
                    if(tmpRemovedMurckoAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        /*If the atom of the murcko fragment is a heteroatom, it is one of the rings we are looking for*/
                        if(tmpBond.getAtom(1).getSymbol() != "C") {
                            return true;
                        }
                        /*If the atom of the ring is a heteroatom, it is one of the rings we are looking for*/
                        if(tmpAtom.getSymbol() != "C") {
                            return true;
                        }
                    }
                }
            }
        }
        /*Treatment for linkers that consist of more than one bond, i.e. at least one atom*/
        //Generate the murcko fragment, as this removes the multiple bonded atoms at the linkers and exocyclic bonds
        IAtomContainer tmpMurcko = aMurckoFragmenter.scaffold(aMolecule);
        for(IAtom tmpAtom : tmpMurcko.atoms()) {
            /*Atom is not part of the murcko fragment from which the ring was removed, nor is it part of the ring under investigation.
            It is therefore a linker atom.*/
            if(!tmpRemovedMurckoAtomNumbers.contains(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) &&
                    !tmpRingPropertyNumbers.contains(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                /*Investigate all bonds of the atom*/
                for(IBond tmpBond : tmpAtom.bonds()) {
                    /*Check if atom 0 of the bond is a heteroatom*/
                    if(tmpBond.getAtom(0).getSymbol() != "C") {
                        /*If the heteroatom is in the ring or in the Murcko fragment with the ring removed, it must be a terminal linker atom*/
                        if(tmpRemovedMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) ||
                                tmpRingPropertyNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                            return true;
                        }
                    }
                    /*Check if atom 1 of the bond is a heteroatom*/
                    if(tmpBond.getAtom(1).getSymbol() != "C") {
                        /*If the heteroatom is in the ring or in the Murcko fragment with the ring removed, it must be a terminal linker atom*/
                        if(tmpRemovedMurckoAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) ||
                                tmpRingPropertyNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                            return true;
                        }
                    }
                }
            }
        }
        //If none of the cases apply, it is not one of the rings we are looking for
        return false;
    }

    /**
     * Remove a ring according to the thirteenth Schuffenhauer rule.
     * Based on rule number 13 from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * In contrast to the paper, unique SMILES are used here instead of canonical SMILES.
     * The entered rings are sorted alphabetically by their unique SMILES. The last ring of this sort is returned.
     * If two structures are the same, one is selected arbitrary.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return Molecule from which the ring selected by the rule has been removed
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    protected IAtomContainer applySchuffenhauerRuleThirteen(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        //Strings are stored in a sorted map. The natural order is alphabetical
        TreeMap<String, IAtomContainer> tmpRingRemovedMap = new TreeMap();//Sorted map
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRingRemoved = this.removeRing(tmpClonedMolecule, tmpRing);
            //Remove linker
            IAtomContainer tmpSchuff = this.getScaffoldInternal(tmpRingRemoved, false, null);
            //A few structures do not produce a truly unique SMILES. These are overwritten and are therefore not considered for further selection.
            tmpRingRemovedMap.put(tmpGenerator.create(tmpSchuff), tmpSchuff);
        }
        //The first key in the map is automatically the SMILES key, which has the lower rank in alphabetical order
        IAtomContainer tmpReturnedStructure = tmpRingRemovedMap.get(tmpRingRemovedMap.firstKey());
        return tmpReturnedStructure;
    }
    //</editor-fold>
    //</editor-fold>
}