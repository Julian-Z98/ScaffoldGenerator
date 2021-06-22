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

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.*;

public class ScaffoldGenerator {
    /**
     * Property of the atoms according to which they are counted and identified
     */
    public static final String SCAFFOLD_ATOM_COUNTER_PROPERTY = "SCAFFOLD_ATOM_COUNTER_PROPERTY";

    //<editor-fold desc="Fundamental Methods">
        /**
         * Generates the Schuffenhauer scaffold for the entered molecule and returns it. All stereochemistry information is deleted.
         * The aromaticity can be set with a selected electron-donation model. Setting the aromaticity is essential for many of the following methods.
         * @param aMolecule molecule whose Schuffenhauer scaffold is produced.
         * @param aIsAromaticitySet Indicates whether the aromaticity is to be set.
         * @param aElectronDonation ElectronDonation Model to be used to determine aromaticity. Can be null if aISAromaticitySet == false.
         * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
         * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
         * @throws CloneNotSupportedException if cloning is not possible.
         */
        public IAtomContainer getSchuffenhauerScaffold(IAtomContainer aMolecule, boolean aIsAromaticitySet, ElectronDonation aElectronDonation) throws CDKException, CloneNotSupportedException {
            IAtomContainer tmpClonedMolecule = aMolecule.clone();
            /*Determine aromaticity*/
            if(aIsAromaticitySet) {
                /*Set aromaticity*/
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpClonedMolecule);
                CDKHydrogenAdder.getInstance(tmpClonedMolecule.getBuilder()).addImplicitHydrogens(tmpClonedMolecule);
                Objects.requireNonNull(aElectronDonation, "If aIsAromaticitySet == true, aElectronDonation must be non null");
                Aromaticity tmpAromaticity = new Aromaticity(aElectronDonation, Cycles.relevant());
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
            HashMap<Integer, IBond> tmpAddAtomMap = new HashMap((tmpClonedMolecule.getAtomCount()/2), 1);
            for(IBond tmpBond: tmpClonedMolecule.bonds()) {
                if(tmpBond.getOrder() != IBond.Order.SINGLE) {//Consider non single bonds
                    //If both atoms of the bond are in the Murcko fragment, they are taken over anyway
                    if(tmpMurckoAtomNumbers.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))
                            && tmpMurckoAtomNumbers.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                        continue;
                    }
                    //The binding has not yet been added to the list
                    if(!tmpAddAtomMap.containsValue(tmpBond)) {
                        /*Add the bond with both atoms*/
                        tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY),tmpBond);
                        tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY),tmpBond);
                    }
                }
            }
            /*Add the missing atom and the respective bond*/
            HashSet<IBond> tmpBondSaved = new HashSet(aMolecule.getBondCount(), 1);
            for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
                int tmpAtomProperty = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                //Every atom that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                if(tmpAddAtomMap.containsKey(tmpAtomProperty)) {
                    //Skip the bond if it has been added before
                    if(tmpBondSaved.contains(tmpAddAtomMap.get(tmpAtomProperty))) {
                        continue;
                    }
                    /*Select the atom that is no longer in the MurckoFragment*/
                    IAtom tmpClonedAtom = null;
                    int tmpAtom1Property = tmpAddAtomMap.get(tmpAtomProperty).getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                    if(tmpAtom1Property == tmpAtomProperty) {
                        tmpClonedAtom = tmpAddAtomMap.get(tmpAtomProperty).getAtom(0).clone();
                    } else {
                        tmpClonedAtom = tmpAddAtomMap.get(tmpAtomProperty).getAtom(1).clone();
                    }
                    /*Add all components to the Murcko fragment*/
                    tmpMurckoFragment.addAtom(tmpClonedAtom); //Add cloned Atom to the molecule
                    //Save the bonds that have been added
                    tmpBondSaved.add(tmpAddAtomMap.get(tmpAtomProperty));
                    //Clone the bond from the original molecule
                    IBond tmpNewBond = tmpAddAtomMap.get(tmpAtomProperty).clone();
                    tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                    tmpNewBond.setAtom(tmpClonedAtom,1); //Add cloned Atom to the bond
                    tmpMurckoFragment.addBond(tmpNewBond); //Add the new bond
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
         * @param aIsKeepingNonSingleBonds if true, double bonded atoms are retained on the ring.
         * @return rings of the inserted molecule.
         */
        public List<IAtomContainer> getRings(IAtomContainer aMolecule, boolean aIsKeepingNonSingleBonds) throws CloneNotSupportedException {
            IAtomContainer tmpClonedMolecule = aMolecule.clone();
            /*Generate cycles*/
            Cycles tmpNewCycles = Cycles.relevant(tmpClonedMolecule);
            IRingSet tmpRingSet = tmpNewCycles.toRingSet();
            List<IAtomContainer> tmpCycles = new ArrayList<>(tmpNewCycles.numberOfCycles());
            int tmpCycleNumber = tmpNewCycles.numberOfCycles();
            //HashMap cannot be larger than the total number of atoms. Key = C and Val = Bond
            HashMap<Integer, IBond> tmpAddAtomMap = new HashMap((tmpClonedMolecule.getAtomCount() / 2), 1);
            /*Store double bonded atoms*/
            if(aIsKeepingNonSingleBonds == true) { //Only needed if double bonded atoms are retained
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
                        if (!tmpAddAtomMap.containsValue(tmpBond)) {
                            /*Add the bond with both atoms*/
                            tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY), tmpBond);
                            tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY), tmpBond);
                        }
                    }
                }
            }
            /*Add Cycles*/
            for(int tmpCount = 0; tmpCount < tmpCycleNumber; tmpCount++) { //Go through all generated rings
                IAtomContainer tmpCycle = tmpRingSet.getAtomContainer(tmpCount); //Store rings as AtomContainer
                if(aIsKeepingNonSingleBonds == true) {
                    /*Add the missing atom and the respective bond*/
                    HashSet<IBond> tmpBondSaved = new HashSet(aMolecule.getBondCount(), 1);
                    for(IAtom tmpAtom : tmpCycle.atoms()) {
                        int tmpAtomProperty = tmpAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                        //Every atom that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                        if(tmpAddAtomMap.containsKey(tmpAtomProperty)) {
                            //Skip the bond if it has been added before
                            if(tmpBondSaved.contains(tmpAddAtomMap.get(tmpAtomProperty))) {
                                continue;
                            }
                            /*Select the atom that is no longer in the MurckoFragment*/
                            IAtom tmpClonedAtom = null;
                            int tmpAtom1Property = tmpAddAtomMap.get(tmpAtomProperty).getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY);
                            if(tmpAtom1Property == tmpAtomProperty) {
                                tmpClonedAtom = tmpAddAtomMap.get(tmpAtomProperty).getAtom(0).clone();
                            } else {
                                tmpClonedAtom = tmpAddAtomMap.get(tmpAtomProperty).getAtom(1).clone();
                            }
                            /*Add all components to the Cycle*/
                            tmpCycle.addAtom(tmpClonedAtom); //Add cloned Atom to the molecule
                            //Save the bonds that have been added
                            tmpBondSaved.add(tmpAddAtomMap.get(tmpAtomProperty));
                            //Clone the bond from the original molecule
                            IBond tmpNewBond = tmpAddAtomMap.get(tmpAtomProperty).clone();
                            tmpNewBond.setAtom(tmpAtom,0); //Add tmpMurckoFragment C to the bond
                            tmpNewBond.setAtom(tmpClonedAtom,1); //Add cloned Atom to the bond
                            tmpCycle.addBond(tmpNewBond); //Add the new bond
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
            /*Get the number of bonds of the ring to other atoms and store the bond atoms of the ring*/
            for(IAtom tmpRingAtom : tmpRingClone.atoms()) {
                for(IAtom tmpMolAtom : tmpMoleculeClone.atoms()) {
                    //All atoms of the ring in the original molecule
                    if(tmpMolAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY) == tmpRingAtom.getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY)) {
                        for(IBond tmpBond : tmpMolAtom.bonds()){
                            //Bond between ring an non ring atom
                            if(tmpIsNotRing.contains(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                                tmpBoundNumber++;
                                //Store ring atom
                                tmpDoNotRemove.add(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                            }
                            //Bond between ring an non ring atom
                            if(tmpIsNotRing.contains(tmpBond.getAtom(1).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY))) {
                                tmpBoundNumber++;
                                //Store ring atom
                                tmpDoNotRemove.add(tmpBond.getAtom(0).getProperty(ScaffoldGenerator.SCAFFOLD_ATOM_COUNTER_PROPERTY));
                            }
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
                        if (tmpMolAtom.getHybridization() == IAtomType.Hybridization.SP2) { //All Atoms that are sp2 hybridised
                            boolean tmpIsSp3 = true;
                            for (IBond tmpBond : tmpMolAtom.bonds()) { //All bonds of the Atom
                                if (tmpBond.getOrder() != IBond.Order.SINGLE) { //If it contains a non single bond it cannot be sp3
                                    tmpIsSp3 = false;
                                }
                            }
                            if (tmpIsSp3) { //If the Atom contains only single bonds, it must be a wanted Atom atom
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
         * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
         */
        public boolean isRingRemovable(IAtomContainer aRing, List<IAtomContainer> aRings, IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
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

    //<editor-fold desc="Output Methods">
        /**
         * Iteratively removes the terminal rings. All resulting Schuffenhauer scaffolds are returned. Duplicates are not permitted.
         * The Schuffenhauer scaffold of the entire entered molecule is stored first in the list.
         * @param aMolecule Molecule to be disassembled.
         * @param aElectronDonation used electron donation model
         * @return List with all resulting Schuffenhauer scaffolds.
         * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
         * @throws CloneNotSupportedException if cloning is not possible.
         */
        public List<IAtomContainer> getIterativeRemoval(IAtomContainer aMolecule, ElectronDonation aElectronDonation) throws CDKException, CloneNotSupportedException {
            SmilesGenerator tmpGenerator = new SmilesGenerator(SmiFlavor.Unique);
            IAtomContainer tmpSchuffenhauerOriginal = this.getSchuffenhauerScaffold(aMolecule, true, aElectronDonation);
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
         * @param aElectronDonation used electron donation model
         * @return List with all resulting Schuffenhauer scaffolds.
         * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
         * @throws CloneNotSupportedException if cloning is not possible.
         */
        public TreeNode<IAtomContainer> getRemovalTree(IAtomContainer aMolecule, ElectronDonation aElectronDonation) throws CDKException, CloneNotSupportedException {
            IAtomContainer tmpSchuffenhauerOriginal = this.getSchuffenhauerScaffold(aMolecule, true, aElectronDonation);
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
         * If there is only one heterocycle of size 3, it is removed.
         * Any ring containing an atom that is not a C is considered a heterocycle.
         * Based on the first rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
         * @param aMolecule Molecule to which the first rule is to be applied
         * @return Molecule to which the first rule was applied
         * @throws CloneNotSupportedException if cloning is not possible.
         * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
         */
        public IAtomContainer applySchuffenhauerRuleOne(IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
            IAtomContainer tmpClonedMolecule = aMolecule;
            List<IAtomContainer> tmpRings = this.getRings(tmpClonedMolecule, true);
            IAtomContainer tmpHeteroRing = null; //Saved size 3 heterocycle
            int tmpHeteroCyclesCounter = 0; //Number of size 3 heterocycles
            /*Investigate how many size 3 heterocycles there are*/
            for(IAtomContainer tmpRing : tmpRings) {
                //All removable and terminal rings of size 3
                if(tmpRing.getAtomCount() == 3 && this.isRingTerminal(tmpClonedMolecule, tmpRing)
                        && this.isRingRemovable(tmpRing, tmpRings, tmpClonedMolecule)) {
                    int tmpHeteroAtomCounter = 0;
                    for(IAtom tmpAtom : tmpRing.atoms()) { //Atoms of the ring
                        if(tmpAtom.getSymbol() != "C"){
                            tmpHeteroAtomCounter++; //Increase if the ring contains a heteroatom
                        }
                    }
                    if(tmpHeteroAtomCounter > 0) { //If it is an heterocycle
                        tmpHeteroCyclesCounter++;
                        tmpHeteroRing = tmpRing; //Save this ring
                    }
                }
            }
            if(tmpHeteroCyclesCounter == 1) { //If there is exactly one heterocycle of size 3
                return this.removeRing(tmpClonedMolecule, tmpHeteroRing); //Molecule with one ring removed
            }
            return tmpClonedMolecule; //Unchanged molecule
        }

    /**
     * If there is only one removable ring with < 12 atoms it will be removed.
     * If there is no ring that has >= 12 atoms or there is not only one small ring, the original molecule is returned.
     * Otherwise the smaller ring is removed.
     * The rule says: Do not remove rings with >= 12 Atoms if there are still smaller rings to remove.
     * Based on the second rule from the "The Scaffold Tree" Paper by Schuffenhauer et al.
     * @param aMolecule Molecule to be treated with the rule.
     * @return Molecule that has been modified or not modified according to the rule
     * @throws CloneNotSupportedException if cloning is not possible.
     * @throws CDKException problem with CDKHydrogenAdder: Throws if insufficient information is present
     */
        public IAtomContainer applySchuffenhauerRuleTwo(IAtomContainer aMolecule) throws CloneNotSupportedException, CDKException {
            IAtomContainer tmpClonedMolecule = this.getSchuffenhauerScaffold(aMolecule, true, ElectronDonation.cdk());
            List<IAtomContainer> tmpRings = this.getRings(tmpClonedMolecule, true);
            List<IAtomContainer> tmpRemovableRings = new ArrayList<>();
            boolean tmpHasMacroCycle = false;
            /*Check whether the molecule has a macrocycle at all*/
            for(IAtomContainer tmpRing : tmpRings) {
                //To determine the ring size, the double-bonded atoms must be removed
                if(this.getRings(tmpRing,false).get(0).getAtomCount() > 11 ) {
                    tmpHasMacroCycle = true; //True if the molecule has a macrocycle
                }
            }
            /*Return the unchanged molecule if there is no macrocycle*/
            if(tmpHasMacroCycle == false) {
                return tmpClonedMolecule; //Unchanged molecule
            }
            /*Identify removable macrocycles and smaller removable rings*/
            boolean tmpHasRemovableMacroCycle = false;
            for(IAtomContainer tmpRing : tmpRings) {
                //Only terminal and removable rings
                if(this.isRingTerminal(tmpClonedMolecule, tmpRing) && this.isRingRemovable(tmpRing, tmpRings, tmpClonedMolecule)) {
                    //To determine the ring size, the double-bonded atoms must be removed
                    if(this.getRings(tmpRing,false).get(0).getAtomCount() > 11 ) {
                        tmpHasRemovableMacroCycle = true;
                    } else { //All removable non macro rings
                        tmpRemovableRings.add(tmpRing);
                    }
                }
            }
            /*Return the unchanged molecule if there is no removable macrocycle*/
            if(tmpHasRemovableMacroCycle == false) {
                return tmpClonedMolecule;
            }
            /*Remove the small ring if there is only one*/
            if (tmpRemovableRings.size() == 1) {
                return this.removeRing(tmpClonedMolecule, tmpRemovableRings.get(0));
            }
            //Return the unchanged molecule if there is more than one small ring that could be removed
            return tmpClonedMolecule;
        }
    //</editor-fold>
}