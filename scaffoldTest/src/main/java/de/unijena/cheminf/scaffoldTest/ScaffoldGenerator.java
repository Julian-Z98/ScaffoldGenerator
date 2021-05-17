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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerComparator;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class ScaffoldGenerator {
    /**
     * Property of the atoms according to which they are counted and identified
     */
    private static String atomCounter = "ScaffoldGeneratorAtomCounter";
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
            tmpAtom.setProperty(this.atomCounter, tmpCounter);
            tmpCounter++;
        }
        //Store the number of each C that is double-bonded to O and the respective bond
        HashMap<Object, IBond> tmpAddAtomMap = new HashMap((tmpMolecule.getAtomCount()/2), 1); //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        for(IBond tmpBond: tmpMolecule.bonds()) {
            if(tmpBond.getOrder() == IBond.Order.DOUBLE) {
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")) { //C in first position in the bond and O in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(this.atomCounter),tmpBond);
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")) { //O in first position in the bond and C in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(this.atomCounter),tmpBond);
                }
            }
        }
        //Generate the murckoFragment
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        IAtomContainer tmpMurckoFragment = tmpMurckoFragmenter.scaffold(tmpMolecule);
        //Add the missing O and the respective bond
        for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
            if(tmpAddAtomMap.containsKey(tmpAtom.getProperty(this.atomCounter))) { //Every C that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                IBond tmpNewBond = null;
                if(tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).getSymbol() == "C") { //C in first position in the bond and O in second position
                    IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(1).clone();
                    tmpMurckoFragment.addAtom(tmpClonedO); //Add cloned O to the molecule
                    tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).clone(); //Clone the bond from the original molecule
                    tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                    tmpNewBond.setAtom(tmpClonedO,1); //Add cloned O to the bond
                }
                if(tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).getSymbol() == "O") { //O in first position in the bond and C in second position
                    IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).clone();
                    tmpMurckoFragment.addAtom(tmpClonedO); //Add cloned O to the molecule
                    tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).clone(); //Clone the bond from the original molecule
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
        HashMap<Object, IBond> tmpAddAtomMap = new HashMap((tmpMolecule.getAtomCount()/2), 1); //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        for(IBond tmpBond: tmpMolecule.bonds()) {
            if(tmpBond.getOrder() == IBond.Order.DOUBLE) {
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")) { //C in first position in the bond and O in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(this.atomCounter),tmpBond);
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")) { //O in first position in the bond and C in second position
                    tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(this.atomCounter),tmpBond);
                }
            }
        }
        //Generate cycles
        Cycles tmpNewCycles = Cycles.mcb(tmpMolecule);
        IRingSet tmpRingSet = tmpNewCycles.toRingSet();
        List<IAtomContainer> tmpCycles = new ArrayList<>(tmpNewCycles.numberOfCycles());
        for(int tmpCount = 0; tmpCount < tmpNewCycles.numberOfCycles(); tmpCount++) { //Go through all generated rings
            IAtomContainer tmpCycle = tmpRingSet.getAtomContainer(tmpCount); //Store rings as AtomContainer
            //Add the missing O and the respective bond
            if(tmpWithDoubleO == true) {
                for(IAtom tmpAtom : tmpCycle.atoms()) {
                    if(tmpAddAtomMap.containsKey(tmpAtom.getProperty(this.atomCounter))) { //Every C that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                        IBond tmpNewBond = null;
                        if(tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).getSymbol() == "C") { //C in first position in the bond and O in second position
                            IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(1).clone();
                            tmpCycle.addAtom(tmpClonedO); //Add cloned O to the molecule
                            tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).clone(); //Clone the bond from the original molecule
                            tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                            tmpNewBond.setAtom(tmpClonedO,1); //Add cloned O to the bond
                        }
                        if(tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).getSymbol() == "O") { //O in first position in the bond and C in second position
                            IAtom tmpClonedO = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).getAtom(0).clone();
                            tmpCycle.addAtom(tmpClonedO); //Add cloned O to the molecule
                            tmpNewBond = tmpAddAtomMap.get(tmpAtom.getProperty(this.atomCounter)).clone(); //Clone the bond from the original molecule
                            tmpNewBond.setAtom(tmpAtom,1); //Add tmpMurckoFragment C to the bond
                            tmpNewBond.setAtom(tmpClonedO,0); //Add cloned O to the bond
                        }
                        tmpCycle.addBond(tmpNewBond); //Add the new bond
                    }
                }
            }
            tmpCycles.add((IAtomContainer) tmpCycle); //Add rings to list
        }
        return tmpCycles;
    }


    /**
     * Removes the given ring from the total molecule and returns it. Important: Property (this.atomCounter) must be set for tmpMolecule/tmpRing and match.
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
        HashSet<Integer> tmpIsNotRing = new HashSet(tmpMolecule.getAtomCount());
        HashSet<Integer> tmpDontRemove = new HashSet(tmpMolecule.getAtomCount());
        int tmpBoundNumber = 0;
        //Store the number of each atom in the molecule
        for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
            tmpIsNotRing.add(tmpMolAtom.getProperty(this.atomCounter));
        }
        //Remove all numbers of the ring that is to be removed
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            tmpIsNotRing.remove(tmpRingAtom.getProperty(this.atomCounter));
        }
        //Get the number of bonds of the ring to other atoms and store the bond atoms of the ring
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                if(tmpMolAtom.getProperty(this.atomCounter) == tmpRingAtom.getProperty(this.atomCounter)) { //All atoms of the ring in the original molecule
                    for(IBond tmpBond : tmpMolAtom.bonds()){
                        if(tmpIsNotRing.contains(tmpBond.getAtom(0).getProperty(this.atomCounter))) { //Bond between ring an non ring atom
                            tmpBoundNumber++;
                            tmpDontRemove.add(tmpBond.getAtom(1).getProperty(this.atomCounter)); //Store ring atom
                        }
                        if(tmpIsNotRing.contains(tmpBond.getAtom(1).getProperty(this.atomCounter))) { //Bond between ring an non ring atom
                            tmpBoundNumber++;
                            tmpDontRemove.add(tmpBond.getAtom(0).getProperty(this.atomCounter)); //Store ring atom
                        }
                    }
                }
            }
        }
        if(tmpBoundNumber < 2) { //Remove all ring atoms, as there are less than two bonds to other atoms
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    if (tmpMolAtom.getProperty(this.atomCounter) == tmpRingAtom.getProperty(this.atomCounter)) {//All atoms of the ring in the original molecule
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms. tmpMoleculeCone.remove() not possible
                    }
                }
            }
        }
        else { //Remove only the ring atoms that are not bound to the rest of the molecule
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All atoms of the ring in the original molecule that are not bound to the rest of the molecule
                    if ((tmpMolAtom.getProperty(this.atomCounter) == tmpRingAtom.getProperty(this.atomCounter)) && !tmpDontRemove.contains(tmpMolAtom.getProperty(this.atomCounter))) {
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms
                    }
                }
            }
        }
        //Add back hydrogens
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Checks whether the tmpRing in the tmpMolecule is terminal. This means whether it can be removed without creating several unconnected parts.
     * Important: Property (this.atomCounter) must be set for tmpMolecule/tmpRing and match.
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
        for(IAtom tmpRingAtom : tmpClonedRing.atoms()) {
            for(IAtom tmpMolAtom : tmpClonedMolecule.atoms()){
                if(tmpMolAtom.getProperty(this.atomCounter) == tmpRingAtom.getProperty(this.atomCounter)) { //All atoms of the ring in the original molecule
                    tmpClonedMolecule.removeAtom(tmpMolAtom); //Remove them
                }
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
        List<IAtomContainer> tmpIterationList = new ArrayList<>();//Größe??
        tmpIterationList.add(this.getSchuffenhauerScaffold(tmpMolecule)); //Add origin SchuffenhauerScaffold
        AtomContainerComparator tmpComparator = new AtomContainerComparator();
        for(int tmpCounter = 0 ; tmpCounter < tmpIterationList.size(); tmpCounter++) { //Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterationList.get(tmpCounter); //Take the next molecule from the list
            for(IAtomContainer tmpRing : this.getRings(tmpIterMol,true)) { //Go through all rings
                if(this.getRings(tmpIterMol,true).size() < 2) { //Skip molcule if it has less than 2 rings
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing)) { //Consider all terminal rings
                    boolean tmpIsInList = false;
                    for(IAtomContainer tmpSavedMol : tmpIterationList) { //Go through the list of all molecules already created
                        if(tmpComparator.compare(tmpSavedMol,this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing)))== 0) {// Check if the molecule has already been added to the list
                            tmpIsInList = true;
                        }
                    }
                    if(tmpIsInList == false) { //Add the molecule only if it is not already in the list
                        tmpIterationList.add(this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing)));
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
        List<IAtomContainer> tmpIterationList = new ArrayList<>();//List of all fragments already created
        List<TreeNode> tmpNodeList = new ArrayList<>();
        tmpIterationList.add(this.getSchuffenhauerScaffold(tmpMolecule)); //Add origin SchuffenhauerScaffold
        TreeNode<IAtomContainer> tmpParentTree = new TreeNode<IAtomContainer>(this.getSchuffenhauerScaffold(tmpMolecule)); //Set origin Schuffenhauer as root
        tmpNodeList.add(tmpParentTree);
        int tmpLevelCounter = 0; //Shows which level of the tree we are currently on.
        for(int tmpCounter = 0 ; tmpCounter < tmpIterationList.size(); tmpCounter++) { //Go through all the molecules created
            IAtomContainer tmpIterMol = tmpIterationList.get(tmpCounter); //Take the next molecule from the list
            for(IAtomContainer tmpRing : this.getRings(tmpIterMol,true)) { //Go through all rings
                if(this.getRings(tmpIterMol,true).size() < 2) { //Skip molcule if it has less than 2 rings
                    continue;
                }
                if(this.isRingTerminal(tmpIterMol, tmpRing)) { //Consider all terminal rings
                    tmpNodeList.add(tmpNodeList.get(tmpLevelCounter).addChild(this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing)))); //Add next node to current Level
                    tmpIterationList.add(this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing))); // The molecule added to the tree is added to the list
                }
            }
            tmpLevelCounter++; //Increases when a level is completed
        }
        return tmpParentTree;
    }

}