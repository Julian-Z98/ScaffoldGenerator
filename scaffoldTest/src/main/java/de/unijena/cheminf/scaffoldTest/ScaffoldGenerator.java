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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class ScaffoldGenerator {
    /**
     * Property of the atoms according to which they are counted and identified
     */
    public static String scaffoldGeneratorAtomCounter = "ScaffoldGeneratorAtomCounter";
    /**
     * Generates the Schuffenhauer scaffold for the entered molecule and returns it.
     * @param tmpMolecule molecule whose Schuffenhauer scaffold is produced.
     * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
     * @throws CDKException if the CDKHydrogenAdder has a problem.
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public IAtomContainer getSchuffenhauerScaffold(IAtomContainer tmpMolecule) throws CDKException, CloneNotSupportedException {
        //Mark each atom with ascending number
        Integer tmpCounter = 0;
        for(IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpAtom.setProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter, tmpCounter);
            tmpCounter++;
        }
        //Store the number of each C that is double-bonded to O and the respective bond
        HashMap<Integer, IBond> tmpAddAtomMap = new HashMap((tmpMolecule.getAtomCount()/2), 1); //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        for(IBond tmpBond: tmpMolecule.bonds()) {
            if(tmpBond.getOrder() == IBond.Order.DOUBLE) {
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")) { //C in first position in the bond and O in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter),tmpBond);
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")) { //O in first position in the bond and C in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter),tmpBond);
                }
            }
        }
        //Generate the murckoFragment
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        IAtomContainer tmpMurckoFragment = tmpMurckoFragmenter.scaffold(tmpMolecule);
        //Add the missing O and the respective bond
        for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
            if(tmpAddAtomMap.containsKey(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) { //Every C that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                IBond tmpNewBond = null;
                if(tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).getSymbol() == "C") { //C in first position in the bond and O in second position
                    IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(1).clone();
                    tmpMurckoFragment.addAtom(tmpClonedO); //Add cloned O to the molecule
                    tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).clone(); //Clone the bond from the original molecule
                    tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                    tmpNewBond.setAtom(tmpClonedO,1); //Add cloned O to the bond
                }
                if(tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).getSymbol() == "O") { //O in first position in the bond and C in second position
                    IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).clone();
                    tmpMurckoFragment.addAtom(tmpClonedO); //Add cloned O to the molecule
                    tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).clone(); //Clone the bond from the original molecule
                    tmpNewBond.setAtom(tmpAtom,1); //Add tmpMurckoFragment C to the bond
                    tmpNewBond.setAtom(tmpClonedO,0); //Add cloned O to the bond
                }
                tmpMurckoFragment.addBond(tmpNewBond); //Add the new bond
            }
        }
        //Add back hydrogens removed by the MurckoFragmenter
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMurckoFragment);
        CDKHydrogenAdder.getInstance(tmpMurckoFragment.getBuilder()).addImplicitHydrogens(tmpMurckoFragment);
        return tmpMurckoFragment;
    }

    /**
     * Generates the smallest set of smallest rings(SSSR) with Cycles.mcb() for the entered molecule, adds double bounded oxygens and returns it.
     * @param tmpMolecule molecule whose rings are produced.
     * @param tmpWithDoubleO if true, double bonded oxygens are retained on the ring.
     * @return rings of the inserted molecule.
     */
    public List<IAtomContainer> getRings(IAtomContainer tmpMolecule, boolean tmpWithDoubleO) throws CloneNotSupportedException {
        //Store the number of each C that is double-bonded to O and the respective bond
        HashMap<Integer, IBond> tmpAddAtomMap = new HashMap((tmpMolecule.getAtomCount()/2), 1); //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        for(IBond tmpBond: tmpMolecule.bonds()) {
            if(tmpBond.getOrder() == IBond.Order.DOUBLE) {
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")) { //C in first position in the bond and O in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter),tmpBond);
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")) { //O in first position in the bond and C in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter),tmpBond);
                }
            }
        }
        //Generate cycles
        Cycles tmpNewCycles = Cycles.mcb(tmpMolecule);
        IRingSet tmpRingSet = tmpNewCycles.toRingSet();
        List<IAtomContainer> tmpCycles = new ArrayList<>(tmpNewCycles.numberOfCycles());
        int tmpCycleNumber = tmpNewCycles.numberOfCycles();
        for(int tmpCount = 0; tmpCount < tmpCycleNumber; tmpCount++) { //Go through all generated rings
            IAtomContainer tmpCycle = tmpRingSet.getAtomContainer(tmpCount); //Store rings as AtomContainer
            //Add the missing O and the respective bond
            if(tmpWithDoubleO == true) {
                for(IAtom tmpAtom : tmpCycle.atoms()) {
                    if(tmpAddAtomMap.containsKey(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) { //Every C that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                        IBond tmpNewBond = null;
                        if(tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).getSymbol() == "C") { //C in first position in the bond and O in second position
                            IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(1).clone();
                            tmpCycle.addAtom(tmpClonedO); //Add cloned O to the molecule
                            tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).clone(); //Clone the bond from the original molecule
                            tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                            tmpNewBond.setAtom(tmpClonedO,1); //Add cloned O to the bond
                        }
                        if(tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).getSymbol() == "O") { //O in first position in the bond and C in second position
                            IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).getAtom(0).clone();
                            tmpCycle.addAtom(tmpClonedO); //Add cloned O to the molecule
                            tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)).clone(); //Clone the bond from the original molecule
                            tmpNewBond.setAtom(tmpAtom,1); //Add tmpMurckoFragment C to the bond
                            tmpNewBond.setAtom(tmpClonedO,0); //Add cloned O to the bond
                        }
                        tmpCycle.addBond(tmpNewBond); //Add the new bond
                    }
                }
            }
            tmpCycles.add(tmpCycle); //Add rings to list
        }
        return tmpCycles;
    }


    /**
     * Removes the given ring from the total molecule and returns it. Important: Property (ScaffoldGenerator.scaffoldGeneratorAtomCounter) must be set for tmpMolecule/tmpRing and match.
     * @param tmpMolecule Molecule whose ring is to be removed.
     * @param tmpRing Ring to be removed.
     * @return Molecule whose ring has been removed.
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException if the CDKHydrogenAdder has a problem.
     */
    public IAtomContainer removeRing(IAtomContainer tmpMolecule, IAtomContainer tmpRing) throws CloneNotSupportedException, CDKException {
        //Clone original molecules
        IAtomContainer tmpMoleculeClone = tmpMolecule.clone();
        IAtomContainer tmpRingClone = tmpRing.clone();
        HashSet<Integer> tmpIsNotRing = new HashSet(tmpMolecule.getAtomCount(), 1);
        HashSet<Integer> tmpDontRemove = new HashSet(tmpMolecule.getAtomCount(), 1);
        int tmpBoundNumber = 0;
        //Store the number of each atom in the molecule
        for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
            tmpIsNotRing.add(tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter));
        }
        //Remove all numbers of the ring that is to be removed
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            tmpIsNotRing.remove(tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter));
        }
        //Get the number of bonds of the ring to other atoms and store the bond atoms of the ring
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                if(tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter) == tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)) { //All atoms of the ring in the original molecule
                    for(IBond tmpBond : tmpMolAtom.bonds()){
                        if(tmpIsNotRing.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) { //Bond between ring an non ring atom
                            tmpBoundNumber++;
                            tmpDontRemove.add(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)); //Store ring atom
                        }
                        if(tmpIsNotRing.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) { //Bond between ring an non ring atom
                            tmpBoundNumber++;
                            tmpDontRemove.add(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)); //Store ring atom
                        }
                    }
                }
            }
        }
        if(tmpBoundNumber < 2) { //Remove all ring atoms, as there are less than two bonds to other atoms
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    if (tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter) == tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)) {//All atoms of the ring in the original molecule
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms. tmpMoleculeCone.remove() not possible
                    }
                }
            }
        }
        else { //Remove only the ring atoms that are not bound to the rest of the molecule
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All atoms of the ring in the original molecule that are not bound to the rest of the molecule
                    if ((tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter) == tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter)) && !tmpDontRemove.contains(tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) {
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms
                    }
                }
            }
            //Set the hybridisation to sp3 for all C that have only single bonds
            for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                if(tmpMolAtom.getSymbol() == "C" && tmpMolAtom.getHybridization() != IAtomType.Hybridization.SP3) { //All C that are not sp3 hybridised
                    boolean tmpIsSp3 = true;
                    for(IBond tmpBond : tmpMolAtom.bonds()) { //All bonds of the C
                        if(tmpBond.getOrder() != IBond.Order.SINGLE) { //If it contains a non single bond it cannot be sp3
                            tmpIsSp3 = false;
                        }
                    }
                    if(tmpIsSp3) { //If the C contains only single bonds, it must be sp3
                        tmpMolAtom.setHybridization(IAtomType.Hybridization.SP3); //Set sp3
                    }
                }
            }
        }
        //Add back hydrogens removed by the MurckoFragmenter
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Checks whether the tmpRing in the tmpMolecule is terminal. This means whether it can be removed without creating several unconnected parts.
     * Important: Property (ScaffoldGenerator.scaffoldGeneratorAtomCounter) must be set for tmpMolecule/tmpRing and match.
     * @param tmpMolecule Molecule whose ring is to be checked
     * @param tmpRing Ring to check
     * @return true if the tmpRing is terminal
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public boolean isRingTerminal(IAtomContainer tmpMolecule, IAtomContainer tmpRing) throws CloneNotSupportedException {
        //Clone moleclule and ring
        IAtomContainer tmpClonedMolecule = tmpMolecule.clone();
        IAtomContainer tmpClonedRing = tmpRing.clone();
        //Remove ring atoms from orininal molecule
        HashMap<Integer, IAtom> tmpMoleculeCounterMap = new HashMap((tmpMolecule.getAtomCount()), 1);
        for(IAtom tmpMolAtom : tmpClonedMolecule.atoms()) { //Save all atoms of the molecule
            tmpMoleculeCounterMap.put(tmpMolAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter), tmpMolAtom);
        }
        for(IAtom tmpRingAtom : tmpClonedRing.atoms()) { // Go through the ring
            if(tmpMoleculeCounterMap.containsKey(tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))) { //Is ring atom in molecule
                tmpClonedMolecule.removeAtom(tmpMoleculeCounterMap.get(tmpRingAtom.getProperty(ScaffoldGenerator.scaffoldGeneratorAtomCounter))); //Remove them
            }
        }
        //Check if there is more than one molecule in the IAtomContainer
        ConnectivityChecker tmpChecker = new ConnectivityChecker();
        boolean tmpRingisTerminal = tmpChecker.isConnected(tmpClonedMolecule);
        return tmpRingisTerminal;
    }

    /**
     * Iteratively removes the terminal rings. All resulting Schuffenhauer scaffolds are returned. Duplicates are not permitted.
     * The Schuffenhauer scaffold of the entire entered molecule is stored first in the list.
     * @param tmpMolecule Molecule to be disassembled.
     * @return List with all resulting Schuffenhauer scaffolds.
     * @throws CDKException if the CDKHydrogenAdder has a problem.
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> getIterativeRemoval(IAtomContainer tmpMolecule) throws CDKException, CloneNotSupportedException {
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        int tmpRingCount = this.getRings(this.getSchuffenhauerScaffold(tmpMolecule), true).size();
        List<String> tmpStringList = new ArrayList<>(tmpRingCount * 45);
        List<IAtomContainer> tmpIterationList = new ArrayList<>(tmpRingCount * 45);//List of all fragments already created and size estimated on the basis of an empirical value
        tmpIterationList.add(this.getSchuffenhauerScaffold(tmpMolecule)); //Add origin SchuffenhauerScaffold
        for(int tmpCounter = 0 ; tmpCounter < tmpIterationList.size(); tmpCounter++) {//Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterationList.get(tmpCounter); //Take the next molecule from the list
            int tmpRingSize = this.getRings(tmpIterMol,true).size();
            for(IAtomContainer tmpRing : this.getRings(tmpIterMol,true)) { //Go through all rings
                if(tmpRingSize < 2) { //Skip molecule if it has less than 2 rings
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing)) { //Consider all terminal rings
                    boolean tmpIsInList = false;
                    IAtomContainer tmpRingRemoved = this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing)); //Remove next ring
                    String tmpRingRemovedSMILES = tmpGenerator.create(tmpRingRemoved); //Generate unique SMILES
                    if(tmpStringList.contains(tmpRingRemovedSMILES)) { //Check if the molecule has already been added to the list
                        tmpIsInList = true;
                    }
                    if(tmpIsInList == false) { //Add the molecule only if it is not already in the list
                        tmpIterationList.add(tmpRingRemoved);
                        tmpStringList.add(tmpRingRemovedSMILES);
                    }
                }
            }
        }
        return tmpIterationList;
    }

    /**
     * Iteratively removes the terminal rings. All resulting Schuffenhauer scaffolds are saved in a tree.
     * A new level is created with each removal step. Duplicates are permitted.
     * The Schuffenhauer scaffold of the entire entered molecule is the root of the tree.
     * @param tmpMolecule Molecule to be disassembled.
     * @return List with all resulting Schuffenhauer scaffolds.
     * @throws CDKException if the CDKHydrogenAdder has a problem.
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public TreeNode<IAtomContainer> getRemovalTree(IAtomContainer tmpMolecule) throws CDKException, CloneNotSupportedException {
        int tmpRingCount = this.getRings(this.getSchuffenhauerScaffold(tmpMolecule), true).size();
        List<IAtomContainer> tmpIterationList = new ArrayList<>(tmpRingCount * 45);// List of all fragments already created and size estimated on the basis of an empirical value
        List<TreeNode> tmpNodeList = new ArrayList<>(); //List of all TreeNodes
        tmpIterationList.add(this.getSchuffenhauerScaffold(tmpMolecule)); //Add origin SchuffenhauerScaffold
        TreeNode<IAtomContainer> tmpParentTree = new TreeNode<IAtomContainer>(this.getSchuffenhauerScaffold(tmpMolecule)); //Set origin Schuffenhauer as root
        tmpNodeList.add(tmpParentTree);
        int tmpLevelCounter = 0; //Shows which level of the tree we are currently on.
        for(int tmpCounter = 0 ; tmpCounter < tmpIterationList.size(); tmpCounter++) { //Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterationList.get(tmpCounter); //Take the next molecule from the list
            List<IAtomContainer> tmpRings = this.getRings(tmpIterMol,true);
            int tmpRingSize = tmpRings.size();
            for(IAtomContainer tmpRing : tmpRings) { //Go through all rings
                if(tmpRingSize < 2) { //Skip molecule if it has less than 2 rings
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing)) { //Consider all terminal rings
                    IAtomContainer tmpRingRemoved = this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing)); //Remove next ring
                    tmpNodeList.add(tmpNodeList.get(tmpLevelCounter).addChild(tmpRingRemoved)); //Add next node to current Level
                    tmpIterationList.add(tmpRingRemoved); // The molecule added to the tree is added to the list
                }
            }
            tmpLevelCounter++; //Increases when a level is completed
        }
        return tmpParentTree;
    }
}