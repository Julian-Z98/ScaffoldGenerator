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

import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ScaffoldGenerator {
    /**
     * Generates the Schuffenhauer scaffold for the entered molecule and returns it.
     * @param tmpMolecule molecule whose Schuffenhauer scaffold is produced.
     * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
     * @throws CDKException if the CDKHydrogenAdder has a problem.
     * @throws CloneNotSupportedException if cloning is not possible
     */
    public IAtomContainer getSchuffenhauerScaffold(IAtomContainer tmpMolecule) throws CDKException, CloneNotSupportedException {
        //Mark each atom with ascending number
        Integer tmpCounter = 0;
        for(IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpAtom.setProperty("AtomCounter", tmpCounter);
            tmpCounter++;
        }
        //Store all pairs of double bounded C and O in HashMap
        HashMap<Integer,Integer> tmpAddAtomMap = new HashMap(tmpMolecule.getAtomCount());//HashMap cannot be larger than the total number of atoms. Key = C and Val = O
        for(IBond tmpBond: tmpMolecule.bonds()) {
            if(tmpBond.getOrder() == IBond.Order.DOUBLE) {
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")) {
                    tmpAddAtomMap.put(tmpBond.getAtom(0).getProperty("AtomCounter"),tmpBond.getAtom(1).getProperty("AtomCounter"));
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")) {
                    tmpAddAtomMap.put(tmpBond.getAtom(1).getProperty("AtomCounter"),tmpBond.getAtom(0).getProperty("AtomCounter"));
                }
            }
        }
        //Generate the murckoFragment
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        tmpMurckoFragmenter.generateFragments(tmpMolecule);
        IAtomContainer tmpMurckoFragment = tmpMurckoFragmenter.scaffold(tmpMolecule);
        //Generate SchuffenhauerScaffold
        for(IAtom tmpAtom : tmpMurckoFragment.atoms()) {
            if(tmpAddAtomMap.containsKey(tmpAtom.getProperty("AtomCounter"))) { //Every C that occurs in the tmpAddAtomMap and in the SchuffenhauerScaffold
                for(IAtom tmpOriginalAtom : tmpMolecule.atoms()) {
                    if(tmpOriginalAtom.getProperty("AtomCounter") == tmpAddAtomMap.get(tmpAtom.getProperty("AtomCounter"))) { //Get the O belonging to the C
                        tmpMurckoFragment.addAtom(tmpOriginalAtom.clone()); //Add the clone of the associated O
                    }
                }
                for(IAtom tmpNewAtom : tmpMurckoFragment.atoms()) {
                    if(tmpNewAtom.getProperty("AtomCounter") == tmpAddAtomMap.get(tmpAtom.getProperty("AtomCounter"))) { //Get the newly added O
                        IAtom tmpOriginalO = null;
                        IAtom tmpOriginalC = null;
                        for(IAtom tmpOriginalAtom : tmpMolecule.atoms()){
                            if(tmpOriginalAtom.getProperty("AtomCounter") == tmpAddAtomMap.get(tmpAtom.getProperty("AtomCounter"))) { //Get the O from the original molecule
                                tmpOriginalO = tmpOriginalAtom;
                            }
                            if(tmpOriginalAtom.getProperty("AtomCounter") == tmpAtom.getProperty("AtomCounter")) { //Get the C from the original molecule
                                tmpOriginalC = tmpOriginalAtom;
                            }
                        }
                        IBond tmpNewBond = tmpMolecule.getBond(tmpOriginalO,tmpOriginalC).clone(); //Clone the bond from the original molecule
                        if(tmpMolecule.getBond(tmpOriginalO,tmpOriginalC).getAtom(0).getSymbol() == "C") { //The first atom of the bond is a C
                            //Add the new atoms to the cloned bond
                            tmpNewBond.setAtom(tmpNewAtom,1);
                            tmpNewBond.setAtom(tmpAtom,0);
                        }
                        if(tmpMolecule.getBond(tmpOriginalO,tmpOriginalC).getAtom(0).getSymbol() == "O") { //The first atom of the bond is a O
                            //Add the new atoms to the cloned bond
                            tmpNewBond.setAtom(tmpNewAtom,0);
                            tmpNewBond.setAtom(tmpAtom,1);
                        }
                        tmpMurckoFragment.addBond(tmpNewBond); //Add the new bond
                    }
                }
            }
        }
        //Add back hydrogens removed by the MurckoFragmenter
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMurckoFragment);
        CDKHydrogenAdder.getInstance(tmpMurckoFragment.getBuilder()).addImplicitHydrogens(tmpMurckoFragment);
        return tmpMurckoFragment;
    }

    /**
     * Generates the smallest set of smallest rings(SSSR) with Cycles.mcb() for the entered molecule and returns it.
     * @param tmpMolecule molecule whose rings are produced.
     * @return rings of the inserted molecule.
     */
    public List<AtomContainer> getRings(IAtomContainer tmpMolecule) {
        Cycles tmpNewCycles = Cycles.mcb(tmpMolecule); //Generate cycles
        List<AtomContainer> tmpCycles = new ArrayList<>();
        for(int tmpCount = 0; tmpCount < tmpNewCycles.numberOfCycles(); tmpCount++) { //Go through all generated rings
            IAtomContainer tmpCycle = tmpNewCycles.toRingSet().getAtomContainer(tmpCount); //Store rings as AtomContainer
            tmpCycles.add((AtomContainer) tmpCycle); //Add rings to list
        }
        return tmpCycles;
    }
}
