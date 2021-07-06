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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

import static junit.framework.TestCase.assertEquals;

public class ScaffoldGenerator {
    //<editor-fold desc="Public Static Final Constants">
    /**
     * Property of the atoms according to which they are counted and identified
     */
    public static final String SCAFFOLD_ATOM_COUNTER_PROPERTY = "SCAFFOLD_ATOM_COUNTER_PROPERTY";
    /**
     * Cycle finder used to detect rings and to determine aromaticity
     */
    public static final CycleFinder CYCLE_FINDER = Cycles.relevant();
    //</editor-fold>

    //<editor-fold desc="Public Methods">
    //<editor-fold desc="Fundamental Methods">
    /**
     * Generates the Schuffenhauer scaffold for the entered molecule and returns it. All stereochemistry information is deleted.
     * The aromaticity can be set with a selected electron-donation model. Setting the aromaticity is essential for many of the following methods.
     * @param aMolecule molecule whose Schuffenhauer scaffold is produced.
     * @param anIsAromaticitySet Indicates whether the aromaticity is to be set.
     * @param anElectronDonation ElectronDonation Model to be used to determine aromaticity. Can be null if anIsAromaticitySet == false.
     * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present or problem with aromaticity.apply()
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public IAtomContainer getSchuffenhauerScaffold(IAtomContainer aMolecule, boolean anIsAromaticitySet, ElectronDonation anElectronDonation) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        /*Determine aromaticity*/
        if(anIsAromaticitySet) {
            /*Set aromaticity*/
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpClonedMolecule);
            CDKHydrogenAdder.getInstance(tmpClonedMolecule.getBuilder()).addImplicitHydrogens(tmpClonedMolecule);
            Objects.requireNonNull(anElectronDonation, "If anIsAromaticitySet == true, anElectronDonation must be non null");
            Aromaticity tmpAromaticity = new Aromaticity(anElectronDonation, ScaffoldGenerator.CYCLE_FINDER);
            tmpAromaticity.apply(tmpClonedMolecule);
        }
        /*Clear the stereo chemistry of the molecule*/
        List<IStereoElement> tmpStereo = new ArrayList<>();
        tmpClonedMolecule.setStereoElements(tmpStereo);
        /*Mark each atom with ascending number*/
        Integer tmpCounter = 0;
        for(IAtom tmpAtom : tmpClonedMolecule.atoms()) {
            tmpAtom.setProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY, tmpCounter);
            tmpCounter++;
        }
        /*Generate the murckoFragment*/
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        IAtomContainer tmpMurckoFragment = tmpMurckoFragmenter.scaffold(tmpClonedMolecule);
        /*Store the number of each Atom of the murckoFragment*/
        HashSet<Integer> tmpMurckoAtomNumbers = new HashSet(tmpClonedMolecule.getAtomCount(), 1);
        for(IAtom tmpMurckoAtom : tmpMurckoFragment.atoms()) {
            tmpMurckoAtomNumbers.add(tmpMurckoAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        /*Store the number of each Atom that is not single bonded and the respective bond*/
        //HashMap cannot be larger than the total number of atoms. Key = Atom and Val = Bond
        HashSet<IBond> tmpAddBondSet = new HashSet((tmpClonedMolecule.getAtomCount()/2), 1);
        for(IBond tmpBond: tmpClonedMolecule.bonds()) {
            if(tmpBond.getOrder() != IBond.Order.SINGLE && tmpBond.getOrder() != IBond.Order.UNSET) {//Consider non single bonds
                //If both atoms of the bond are in the Murcko fragment, they are taken over anyway
                if(tmpMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
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
        /*Add the missing atom and the respective bond*/
        HashMap<Integer, IAtom> tmpMurckoProperties = new HashMap(tmpMurckoFragment.getAtomCount(), 1);
        for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
            /*Save the properties of the murcko fragment*/
            int tmpAtomProperty = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
            tmpMurckoProperties.put(tmpAtomProperty, tmpAtom);
        }
        for(IBond tmpBond : tmpAddBondSet) { //Go though all saved bonds
            /*If both atoms of the bond are contained in the murcko fragment, this bond does not need to be added any more*/
            if(tmpMurckoProperties.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) &&
                    tmpMurckoProperties.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                continue; //Skip this bond
            }
            /*Atom 1 of the bond is in the Murcko fragment*/
            if(tmpMurckoProperties.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                IAtom tmpClonedAtom = tmpBond.getAtom(0).clone();
                tmpMurckoFragment.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                IBond tmpNewBond = tmpBond.clone();
                //Set the first atom
                tmpNewBond.setAtom(tmpMurckoProperties.get(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 1);
                tmpNewBond.setAtom(tmpClonedAtom, 0); //Set the second atom
                tmpMurckoFragment.addBond(tmpNewBond); //Add the whole bond
                continue; //Next bond
            }
            /*Atom 0 of the bond is in the Murcko fragment*/
            if(tmpMurckoProperties.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                IAtom tmpClonedAtom = tmpBond.getAtom(1).clone();
                tmpMurckoFragment.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                IBond tmpNewBond = tmpBond.clone();
                //Set the first atom
                tmpNewBond.setAtom(tmpMurckoProperties.get(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 0);
                tmpNewBond.setAtom(tmpClonedAtom, 1); //Set the second atom
                tmpMurckoFragment.addBond(tmpNewBond); //Add the whole bond
            }
        }
        /*Add back hydrogens removed by the MurckoFragmenter*/
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMurckoFragment);
        CDKHydrogenAdder.getInstance(tmpMurckoFragment.getBuilder()).addImplicitHydrogens(tmpMurckoFragment);
        return tmpMurckoFragment;
    }

    /**
     * Generates the smallest set of smallest rings(SSSR) with Cycles.relevant() for the entered molecule, adds double bounded oxygens and returns it.
     * @param aMolecule molecule whose rings are produced.
     * @param anIsKeepingNonSingleBonds if true, double bonded atoms are retained on the ring.
     * @return rings of the inserted molecule.
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws Intractable thrown if problem could not be solved within some predefined bounds.
     */
    public List<IAtomContainer> getRings(IAtomContainer aMolecule, boolean anIsKeepingNonSingleBonds) throws CloneNotSupportedException, Intractable {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        /*Generate cycles*/
        CycleFinder tmpCycleFinder = ScaffoldGenerator.CYCLE_FINDER;
        Cycles tmpNewCycles = tmpCycleFinder.find(tmpClonedMolecule);
        //Cycles tmpNewCycles = Cycles.relevant(tmpClonedMolecule);
        IRingSet tmpRingSet = tmpNewCycles.toRingSet();
        List<IAtomContainer> tmpCycles = new ArrayList<>(tmpNewCycles.numberOfCycles());
        int tmpCycleNumber = tmpNewCycles.numberOfCycles();
        //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
        HashSet<IBond> tmpAddBondSet = new HashSet((tmpClonedMolecule.getAtomCount() / 2), 1);
        /*Store double bonded atoms*/
        if(anIsKeepingNonSingleBonds == true) { //Only needed if double bonded atoms are retained
            /*Generate the murckoFragment*/
            MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true, 1);
            tmpMurckoFragmenter.setComputeRingFragments(false);
            IAtomContainer tmpMurckoFragment = tmpMurckoFragmenter.scaffold(tmpClonedMolecule);
            /*Store the number of each Atom of the murckoFragment*/
            HashSet<Integer> tmpMurckoAtomNumbers = new HashSet(tmpClonedMolecule.getAtomCount(), 1);
            for (IAtom tmpMurckoAtom : tmpMurckoFragment.atoms()) {
                tmpMurckoAtomNumbers.add(tmpMurckoAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
            /*Store the number of each Atom that is not single bonded and the respective bond*/
            for (IBond tmpBond : tmpClonedMolecule.bonds()) {
                if (tmpBond.getOrder() != IBond.Order.SINGLE) {//Consider non single bonds
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
                HashMap<Integer, IAtom> tmpMurckoProperties = new HashMap(tmpCycle.getAtomCount(), 1);
                for(IAtom tmpAtom : tmpCycle.atoms()) {
                    /*Save the properties of the murcko fragment*/
                    int tmpAtomProperty = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    tmpMurckoProperties.put(tmpAtomProperty, tmpAtom);
                }
                for(IBond tmpBond : tmpAddBondSet) { //Go though all saved bonds
                    /*If both atoms of the bond are contained in the murcko fragment, this bond does not need to be added any more*/
                    if(tmpMurckoProperties.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) &&
                            tmpMurckoProperties.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        continue; //Skip this bond
                    }
                    /*Atom 1 of the bond is in the Murcko fragment*/
                    if(tmpMurckoProperties.containsKey(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(0).clone();
                        tmpCycle.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoProperties.get(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 1);
                        tmpNewBond.setAtom(tmpClonedAtom, 0); //Set the second atom
                        tmpCycle.addBond(tmpNewBond); //Add the whole bond
                        continue; //Next bond
                    }
                    /*Atom 0 of the bond is in the Murcko fragment*/
                    if(tmpMurckoProperties.containsKey(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        IAtom tmpClonedAtom = tmpBond.getAtom(1).clone();
                        tmpCycle.addAtom(tmpClonedAtom); //Add the atom that is not yet in the murcko fragment
                        IBond tmpNewBond = tmpBond.clone();
                        //Set the first atom
                        tmpNewBond.setAtom(tmpMurckoProperties.get(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)), 0);
                        tmpNewBond.setAtom(tmpClonedAtom, 1); //Set the second atom
                        tmpCycle.addBond(tmpNewBond); //Add the whole bond
                    }
                }
            }
            tmpCycles.add(tmpCycle); //Add rings to list
        }
        return tmpCycles;
    }


    /**
     * Removes the given ring from the total molecule and returns it.
     * Preserves the sp2 hybridisation of a border atom when an aromatic ring is removed.
     * With the removal of heterocycles of size 3 a double bond is inserted if it is directly adjacent to another ring.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for tmpMolecule/tmpRing and match.
     * @param aMolecule Molecule whose ring is to be removed.
     * @param aRing Ring to be removed.
     * @return Molecule whose ring has been removed.
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
    public IAtomContainer removeRing(IAtomContainer aMolecule, IAtomContainer aRing) throws CloneNotSupportedException, CDKException {
        /*Clone original molecules*/
        IAtomContainer tmpMoleculeClone = aMolecule.clone();
        IAtomContainer tmpRingClone = aRing.clone();
        boolean tmpIsRingAromatic = true;
        HashSet<Integer> tmpIsNotRing = new HashSet(aMolecule.getAtomCount(), 1);
        HashSet<Integer> tmpDoNotRemove = new HashSet(aMolecule.getAtomCount(), 1);
        int tmpBoundNumber = 0;
        /*Preparation for insertion of double bonds with removal of aromatic rings*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet(tmpMoleculeClone.getAtomCount(), 1);
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
                        //Bond between ring an non ring atom
                        if(tmpIsNotRing.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) ||
                                tmpIsNotRing.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                            tmpBoundNumber++;
                        }
                    }
                }
            }
        }
        /*Add all atoms of rings that are not to be removed to tmpDoNotRemove*/
        CycleFinder tmpCycleFinder = ScaffoldGenerator.CYCLE_FINDER;
        //Get all cycles of the molecule
        Cycles tmpCycles = tmpCycleFinder.find(tmpMoleculeClone);
        HashSet<Integer> tmpRingProperties = new HashSet(tmpRingClone.getAtomCount(), 1);
        //Save the properties of the ring
        for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
            tmpRingProperties.add(tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
        }
        for(IAtomContainer tmpCycle : tmpCycles.toRingSet().atomContainers()) {
            boolean tmpIsRingToRemove = true;
            /*Check if it is the ring to be removed*/
            for(IAtom tmpCycleAtom : tmpCycle.atoms()) {
                //If one of the atoms of the ring to be removed is not included, it is not this ring
                if(!tmpRingProperties.contains(tmpCycleAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                    tmpIsRingToRemove = false;
                }
            }
            /*If it is not the ring you want to remove, add its atoms to the tmpDoNotRemove list*/
            if(tmpIsRingToRemove == false) {
                for(IAtom tmpCycleAtom : tmpCycle.atoms()) {
                    Integer tmpProperty = tmpCycleAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    if(!tmpDoNotRemove.contains(tmpProperty)) {
                        tmpDoNotRemove.add(tmpProperty);
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
                    for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) { //go though the whole molecule
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
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                if(!tmpRingAtom.isAromatic()) {
                    tmpIsRingAromatic = false; //If one atom is non aromatic the whole ring is not aromatic
                }
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    /*All atoms of the ring in the original molecule that are not bound to the rest of the molecule*/
                    if ((tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)
                            == tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) && !tmpDoNotRemove.contains(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        tmpMoleculeClone.removeAtom(tmpMolAtom); //Remove atoms
                    }
                }
            }
            /*Store the number of all atoms from which an aromatic ring has been removed.
            * In these atoms, a double bond was removed without changing the hybridisation from sp2 to sp3.*/
            if(tmpIsRingAromatic) { //Perform calculation only if the ring to be removed is aromatic
                for (IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All Atoms that are sp2 hybridised and in the ring to be removed
                    if (tmpMolAtom.getHybridization() == IAtomType.Hybridization.SP2
                            && tmpRingProperties.contains(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        boolean tmpIsSp3 = true;
                        for (IBond tmpBond : tmpMolAtom.bonds()) { //All bonds of the Atom
                            if (tmpBond.getOrder() != IBond.Order.SINGLE) { //If it contains a non single bond it cannot be sp3
                                tmpIsSp3 = false;
                            }
                        }
                        if (tmpIsSp3) { //If the Atom contains only single bonds, it must be a wanted atom
                            tmpEdgeAtomNumbers.add(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                        }
                    }
                }
            }
        }
        if(tmpIsRingAromatic) { //Perform calculation only if the ring to be removed is aromatic
            for(IBond tmpBond : tmpMoleculeClone.bonds()) {
                /*If both atoms of a bond were previously part of an aromatic ring, insert a double bond*/
                if(tmpEdgeAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
                    && tmpEdgeAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                    tmpBond.setOrder(IBond.Order.DOUBLE);
                }
            }
        }
        /*Clear hybridisation. The hybridisation must be reset later by percieveAtomTypesAndConfigureAtoms, as the hybridisation is not changed on its own when the atoms are removed.
        sp2 atoms whose double bonds have been removed must be declared as sp3.*/
        for(IAtom tmpAtom : tmpMoleculeClone.atoms()) {
            tmpAtom.setHybridization((IAtomType.Hybridization) CDKConstants.UNSET);
        }
        /*Add back hydrogens removed by the MurckoFragmenter*/
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoleculeClone);
        CDKHydrogenAdder.getInstance(tmpMoleculeClone.getBuilder()).addImplicitHydrogens(tmpMoleculeClone);
        return tmpMoleculeClone;
    }

    /**
     * Checks whether the tmpRing in the tmpMolecule is terminal. This means whether it can be removed without creating several unconnected parts.
     * Important: Property (ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) must be set for tmpMolecule/tmpRing and match.
     * @param aMolecule Molecule whose ring is to be checked
     * @param aRing Ring to check
     * @return true if the tmpRing is terminal
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public boolean isRingTerminal(IAtomContainer aMolecule, IAtomContainer aRing) throws CloneNotSupportedException {
        /*Clone molecule and ring*/
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpClonedRing = aRing.clone();
        /*Remove ring atoms from original molecule*/
        HashMap<Integer, IAtom> tmpMoleculeCounterMap = new HashMap((aMolecule.getAtomCount()), 1);
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
     * Checks whether rings may be removed. These are rings in no atom belongs to another ring.
     * Furthermore, removal is impossible when it is an aromatic ring, that borders two consecutive rings.
     * getSchuffenhauerScaffold() must have been previously applied to the molecule so that an aromaticity status has been assigned to the atoms.
     * @param aRing Ring being tested for its removability
     * @param aRings All Rings of the molecule
     * @param aMolecule Whole molecule
     * @return Whether the ring is removable
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public boolean isRingRemovable(IAtomContainer aRing, List<IAtomContainer> aRings, IAtomContainer aMolecule) throws CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        IAtomContainer tmpClonedRing = aRing.clone();

        /*---Recognition of rings in which no atom belongs to another ring---*/
        List<IAtomContainer> tmpClonedRings = new ArrayList<>(aRings.size());
        HashSet<Integer> tmpRingsNumbers = new HashSet(aRings.size() * tmpClonedRing.getAtomCount(), 1);
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
        for(IAtom tmpAtom : tmpClonedRing.atoms()) {
            if(!tmpAtom.isAromatic()) {
                return true;
            }
        }
        /*Store the number of all ring atoms*/
        HashSet<Integer> tmpMoleculeNumbers = new HashSet(tmpClonedMolecule.getAtomCount(), 1);
        for(IAtomContainer tmpRing : tmpClonedRings) {
            for(IAtom tmpRingAtom : tmpRing.atoms()) {
                int tmpAtomNumber = tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                if(!tmpMoleculeNumbers.contains(tmpAtomNumber)) {
                    tmpMoleculeNumbers.add(tmpAtomNumber);
                }
            }
        }
        /*Store all the atoms of the other rings bordering the aromatic ring*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet(tmpClonedMolecule.getAtomCount(), 1);
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
    //</editor-fold>

    /**
     * Checks whether the ring of a molecule is in an aromatic fused ring system. These systems cannot easily be further disassembled.
     * @param aRing Ring tested to see if it is in an aromatic fused ring system
     * @param aRings All Rings of the molecule
     * @param aMolecule Whole molecule
     * @return Whether the ring is part of an aromatic fused ring system
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public boolean hasFusedRings(IAtomContainer aRing, List<IAtomContainer> aRings, IAtomContainer aMolecule) throws CloneNotSupportedException {
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
            boolean tmpIsRingAromatic = true;
            /*Store the atoms of the aromatic rings*/
            for(IAtom tmpAtom : tmpRing.atoms()) {
                if(!tmpAtom.isAromatic()) {
                    tmpIsRingAromatic = false;
                }
            }
            if(tmpIsRingAromatic == false) { //Skip non aromatic rings
                continue;
            }
            tmpClonedRings.add(tmpRing.clone()); //Store the aromatic rings
            for(IAtom tmpAtom : tmpRing.atoms()) { //Store the atoms of the aromatic rings
                tmpRingsNumbers.add(tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
            }
        }

        /*Store all the atoms of the other rings bordering the aromatic ring*/
        HashSet<Integer> tmpEdgeAtomNumbers = new HashSet(tmpClonedMolecule.getAtomCount(), 1);
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

    //<editor-fold desc="Output Methods">
    /**
     * Iteratively removes the terminal rings. All resulting Schuffenhauer scaffolds are returned. Duplicates are not permitted.
     * The Schuffenhauer scaffold of the entire entered molecule is stored first in the list.
     * @param aMolecule Molecule to be disassembled.
     * @param anElectronDonation used electron donation model
     * @return List with all resulting Schuffenhauer scaffolds.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> getIterativeRemoval(IAtomContainer aMolecule, ElectronDonation anElectronDonation) throws CDKException, CloneNotSupportedException {
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        IAtomContainer tmpSchuffenhauerOriginal = this.getSchuffenhauerScaffold(aMolecule, true, anElectronDonation);
        int tmpRingCount = this.getRings(tmpSchuffenhauerOriginal, true).size();
        List<String> tmpAddedSMILESList = new ArrayList<>(tmpRingCount * 45);
        //List of all fragments already created and size estimated on the basis of an empirical value
        List<IAtomContainer> tmpIterativeRemovalList = new ArrayList<>(tmpRingCount * 45);
        tmpIterativeRemovalList.add(tmpSchuffenhauerOriginal); //Add origin SchuffenhauerScaffold
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
                    IAtomContainer tmpRingRemoved = this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing), false, null); //Remove next ring
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
     * Iteratively removes the terminal rings. All resulting Schuffenhauer scaffolds are saved in a tree.
     * A new level is created with each removal step. Duplicates are permitted.
     * The Schuffenhauer scaffold of the entire entered molecule is the root of the tree.
     * @param aMolecule Molecule to be disassembled.
     * @param anElectronDonation used electron donation model
     * @return List with all resulting Schuffenhauer scaffolds.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public TreeNode<IAtomContainer> getRemovalTree(IAtomContainer aMolecule, ElectronDonation anElectronDonation) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpSchuffenhauerOriginal = this.getSchuffenhauerScaffold(aMolecule, true, anElectronDonation);
        int tmpRingCount = this.getRings(tmpSchuffenhauerOriginal, true).size();
        //List of all fragments already created and size estimated on the basis of an empirical value
        List<IAtomContainer> tmpIterativeRemovalList = new ArrayList<>(tmpRingCount * 45);
        List<TreeNode> tmpAllNodesList = new ArrayList<>(); //List of all TreeNodes
        tmpIterativeRemovalList.add(tmpSchuffenhauerOriginal); //Add origin SchuffenhauerScaffold
        TreeNode<IAtomContainer> tmpParentNode = new TreeNode<IAtomContainer>(tmpSchuffenhauerOriginal); //Set origin Schuffenhauer as root
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
                    IAtomContainer tmpRingRemoved = this.getSchuffenhauerScaffold(this.removeRing(tmpIterMol, tmpRing), false, null); //Remove next ring
                    tmpAllNodesList.add(tmpAllNodesList.get(tmpLevelCounter).addChild(tmpRingRemoved)); //Add next node to current Level
                    tmpIterativeRemovalList.add(tmpRingRemoved); // The molecule added to the tree is added to the list
                }
            }
            tmpLevelCounter++; //Increases when a level is completed
        }
        return tmpParentNode;
    }
    //</editor-fold>
    //<editor-fold desc="Schuffenhauer Rules">
    /**
     * Iteratively removes the rings of the molecule according to specific rules that are queried hierarchically.
     * Based on the rules from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * @param aMolecule Molecule that is to be broken down into its fragments
     * @param anElectronDonation Electron donation model that should be used
     * @return Fragments of the molecule according to the Schuffenhauer rules
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> applySchuffenhauerRules(IAtomContainer aMolecule, ElectronDonation anElectronDonation) throws CloneNotSupportedException, CDKException {
        IAtomContainer tmpSchuffenhauer = this.getSchuffenhauerScaffold(aMolecule,true ,anElectronDonation );
        //List of all generated fragments
        List<IAtomContainer> tmpSchuffenhauerFragments = new ArrayList<>(this.getRings(tmpSchuffenhauer, false).size());
        tmpSchuffenhauerFragments.add(tmpSchuffenhauer);
        /*Go through all the fragments generated and try to break them down further*/
        for(int tmpCounter = 0 ; tmpCounter < tmpSchuffenhauerFragments.size(); tmpCounter++) {
            List<IAtomContainer> tmpRings = this.getRings(tmpSchuffenhauerFragments.get(tmpCounter), true);
            /*If the fragment has only one ring or no ring, it does not need to be disassembled further*/
            if(tmpRings.size() == 1 || tmpRings.size() == 0) {
                break;
            }
            /*Only the removable terminal rings are further investigated*/
            List<IAtomContainer> tmpRemovableRings = new ArrayList<>(tmpRings.size());
            for (IAtomContainer tmpRing : tmpRings) {
                if (this.isRingTerminal(tmpSchuffenhauerFragments.get(tmpCounter), tmpRing)
                        && this.isRingRemovable(tmpRing, tmpRings, tmpSchuffenhauerFragments.get(tmpCounter))) {
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
                //Remove the ring from from the fragment currently being treated
                IAtomContainer tmpRingRemoved = this.removeRing(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings.get(0));
                //Remove the linkers
                IAtomContainer tmpSchuffRingRemoved = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
                //Add the fragment to the list of fragments
                tmpSchuffenhauerFragments.add(tmpSchuffRingRemoved);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number two*/
            tmpRemovableRings = this.applySchuffenhauerRuleTwo(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                //Remove the ring from from the fragment currently being treated
                IAtomContainer tmpRingRemoved = this.removeRing(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings.get(0));
                //Remove the linkers
                IAtomContainer tmpSchuffRingRemoved = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
                //Add the fragment to the list of fragments
                tmpSchuffenhauerFragments.add(tmpSchuffRingRemoved);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number three*/
            tmpRemovableRings = this.applySchuffenhauerRuleThree(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                //Remove the ring from from the fragment currently being treated
                IAtomContainer tmpRingRemoved = this.removeRing(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings.get(0));
                //Remove the linkers
                IAtomContainer tmpSchuffRingRemoved = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
                //Add the fragment to the list of fragments
                tmpSchuffenhauerFragments.add(tmpSchuffRingRemoved);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number four and five*/
            tmpRemovableRings = this.applySchuffenhauerRuleFourAndFive(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                //Remove the ring from from the fragment currently being treated
                IAtomContainer tmpRingRemoved = this.removeRing(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings.get(0));
                //Remove the linkers
                IAtomContainer tmpSchuffRingRemoved = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
                //Add the fragment to the list of fragments
                tmpSchuffenhauerFragments.add(tmpSchuffRingRemoved);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number six*/
            tmpRemovableRings = this.applySchuffenhauerRuleSix(tmpRemovableRings);
            if (tmpRemovableRings.size() == 1) { //If only one eligible ring remains, it can be removed
                //Remove the ring from from the fragment currently being treated
                IAtomContainer tmpRingRemoved = this.removeRing(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings.get(0));
                //Remove the linkers
                IAtomContainer tmpSchuffRingRemoved = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
                //Add the fragment to the list of fragments
                tmpSchuffenhauerFragments.add(tmpSchuffRingRemoved);
                //After a new fragment has been added, the next one is investigated
                continue;
            }
            /*Apply rule number thirteen, the tiebreaking rule */
            tmpSchuffenhauerFragments.add(this.applySchuffenhauerRuleThirteen(tmpSchuffenhauerFragments.get(tmpSchuffenhauerFragments.size() - 1), tmpRemovableRings));
        }
        return tmpSchuffenhauerFragments;
    }

    /**
     * Sort out the rings according to the first Schuffenhauer rule.
     * Based on the first rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * The rule says: Remove Heterocycles of Size 3 First.
     * Therefore, size 3 hetero rings are preferred when available.
     * Only these rings will be returned if present. If none are present, all rings entered will be returned.
     * @param aRings Rings to which the first rule is to be applied
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     */
    public List<IAtomContainer> applySchuffenhauerRuleOne(List<IAtomContainer> aRings) {
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
                if(tmpHeteroAtomCounter == 1) { //If it is an heterocycle
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
    public List<IAtomContainer> applySchuffenhauerRuleTwo(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpSmallRings = new ArrayList<>(aRings.size()); //Rings smaller 12
        /*Identify macrocycles and smaller rings*/
        boolean tmpHasRemovableMacroCycle = false;
        for(IAtomContainer tmpRing : aRings) {
            /*To determine the ring size, the exocyclic atoms must be removed*/
            CycleFinder tmpCycleFinder = ScaffoldGenerator.CYCLE_FINDER;
            Cycles tmpRemovedExocyclic = tmpCycleFinder.find(tmpRing);
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
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> applySchuffenhauerRuleThree(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        List<IAtomContainer> tmpRemoveRings = new ArrayList<>(aRings.size()); //Rings with the longest linker
        List<Integer> tmpLinkerSize = new ArrayList<>(aRings.size()); //Linker length of each ring
        Integer tmpMoleculeAtomCount = tmpClonedMolecule.getAtomCount();
        /*Calculate the linker length of each ring. Negative integers are fused rings*/
        for(IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRemovedRing = this.removeRing(tmpClonedMolecule, tmpRing);
            IAtomContainer tmpSchuff = this.getSchuffenhauerScaffold(tmpRemovedRing, false, null);
            //The number of atoms of the removed ring and the molecule from which the ring and the linker were removed are subtracted from the atomic number of the whole molecule
            //This leaves only the atomic number of the linker
            tmpLinkerSize.add(tmpMoleculeAtomCount - (tmpRing.getAtomCount() + tmpSchuff.getAtomCount()));
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
     * The fifth rule says: Bridged Ring Systems Are Retained with Preference over Spiro Ring Systems.
     * Therefore, the rings with the positive maximum delta are preferred over the rings with the negative one.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return List of rings to be removed first according to the rule. Returns the unchanged list if the rule cannot be applied to the rings.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public List<IAtomContainer> applySchuffenhauerRuleFourAndFive(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        List<IAtomContainer> tmpRingsReturn = new ArrayList<>(aRings.size()); //Rings that are returned
        List<Integer> tmpDeltaList = new ArrayList<>(aRings.size()); //Delta values of all rings
        List<Integer> tmpDeltaListAbs = new ArrayList<>(aRings.size()); //Absolute Delta values of all rings
        /*Calculate the delta values for all rings*/
        for(IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRingRemoved = this.removeRing(tmpClonedMolecule, tmpRing); //Remove the ring
            CycleFinder tmpCycleFinder = ScaffoldGenerator.CYCLE_FINDER;
            Cycles tmpCycles = tmpCycleFinder.find(tmpRingRemoved); //get cycle number(nR)
            List<IBond> tmpCycleBonds = new ArrayList<>(aRings.size());
            int tmpFusedRingBondCounter = 0; // Number of bonds being a member in more than one ring(nrrb)
            /*Count nrrb*/
            for(IAtomContainer tmpCycle : tmpCycles.toRingSet().atomContainers()) { //Go through all cycle
                for(IBond tmpBond : tmpCycle.bonds()) { //Go through all bonds of each cycle
                    //If the bond is already included in the list, it occurs in several rings
                    if(tmpCycleBonds.contains(tmpBond)) {
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
        Integer tmpMaxAbsList = tmpDeltaListAbs.stream().mapToInt(v->v).max().orElseThrow(NoSuchElementException::new);
        Integer tmpMaxList = tmpDeltaList.stream().mapToInt(v->v).max().orElseThrow(NoSuchElementException::new);
        if(tmpMaxAbsList > 0) {
            /*Rule five: if there is a positive maximum delta, only get the maximum positive deltas*/
            if(tmpMaxAbsList == tmpMaxList) {
                /* Add all rings that have the highest delta to the list*/
                for(int tmpCounter = 0 ; tmpCounter < tmpDeltaList.size(); tmpCounter++) {
                    if(tmpDeltaList.get(tmpCounter) == tmpMaxAbsList) {
                        tmpRingsReturn.add(aRings.get(tmpCounter));
                    }
                }
                return tmpRingsReturn; //All rings that have the highest delta
            }
            /*Rule four: Add all rings that have the highest absolute delta to the list*/
            for(int tmpCounter = 0 ; tmpCounter < tmpDeltaListAbs.size(); tmpCounter++) {
                if(tmpDeltaListAbs.get(tmpCounter) == tmpMaxAbsList) {
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
    public List<IAtomContainer> applySchuffenhauerRuleSix(List<IAtomContainer> aRings) throws CDKException {
        List<IAtomContainer> tmpReturnRingList = new ArrayList<>(aRings.size());
        /*Size 3, 5 and 6 rings will be added to the list if present*/
        CycleFinder tmpCycleFinder = ScaffoldGenerator.CYCLE_FINDER;
        for(IAtomContainer tmpRing : aRings) {
            //To determine the ring size, the exocyclic atoms must be removed
            Cycles tmpRemovedExocyclic = tmpCycleFinder.find(tmpRing);
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
     * Remove a ring according to the thirteenth Schuffenhauer rule.
     * Based on rule number 13 from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * In contrast to the paper, unique SMILES are used here instead of canonical SMILES.
     * The entered rings are sorted alphabetically by their unique SMILES. The last ring of this sort is returned.
     * A few structures do not produce a truly unique SMILES. These are overwritten and are therefore not considered for further selection.
     * @param aRings Removable rings of the molecule to which the rule is applied
     * @param aMolecule Molecule from which a ring is to be removed
     * @return Molecule from which the ring selected by the rule has been removed
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     * @throws CloneNotSupportedException if cloning is not possible.
     */
    public IAtomContainer applySchuffenhauerRuleThirteen(IAtomContainer aMolecule, List<IAtomContainer> aRings) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpClonedMolecule = aMolecule.clone();
        TreeMap<String, IAtomContainer> tmpRingRemovedMap = new TreeMap();//Sorted map
        SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
        for (IAtomContainer tmpRing : aRings) {
            IAtomContainer tmpRingRemoved = this.removeRing(tmpClonedMolecule, tmpRing);
            //Remove linker
            IAtomContainer tmpSchuff = this.getSchuffenhauerScaffold(tmpRingRemoved, false, null);
            //A few structures do not produce a truly unique SMILES. These are overwritten and are therefore not considered for further selection.
            tmpRingRemovedMap.put(tmpGenerator.create(tmpSchuff), tmpSchuff);
        }
        //The last key in the map is automatically the SMILES key, which is at the end of the alphabetical order
        return tmpRingRemovedMap.get(tmpRingRemovedMap.lastKey());
    }
    //</editor-fold>
    //</editor-fold>
}