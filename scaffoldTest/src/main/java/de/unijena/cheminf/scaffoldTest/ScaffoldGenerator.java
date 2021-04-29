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
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.ArrayList;
import java.util.List;

public class ScaffoldGenerator {
    public static void main(String[] args) {

    }
    /**
     * Generates the Schuffenhauer scaffold for the entered molecule and returns it.
     * @param tmpMolecule molecule whose Schuffenhauer scaffold is produced.
     * @return Schuffenhauer scaffold of the inserted molecule. It can be an empty molecule if the original molecule does not contain a Schuffenhauer scaffold.
     * @throws CDKException if an error occurs when using the MurckoFragmenter.
     */
    public IAtomContainer getSchuffenhauerScaffold(IAtomContainer tmpMolecule) throws CDKException {
        //Mark each atom with ascending number
        Integer tmpCounter = 0;
        for(IAtom tmpAtom : tmpMolecule.atoms()) {
            tmpCounter++;
            tmpAtom.setProperty("AtomCounter", tmpCounter);
        }
        //Store every C to which an O is double-bonded
        List<Integer> tmpAddAtomList = new ArrayList<>();
        for(IBond tmpBond: tmpMolecule.bonds()){
            if(tmpBond.getOrder() == IBond.Order.DOUBLE){
                if(tmpBond.getAtom(0).getSymbol().equals("C") && tmpBond.getAtom(1).getSymbol().equals("O")){
                    tmpAddAtomList.add(tmpBond.getAtom(0).getProperty("AtomCounter"));
                }
                if(tmpBond.getAtom(0).getSymbol().equals("O") && tmpBond.getAtom(1).getSymbol().equals("C")){
                    tmpAddAtomList.add(tmpBond.getAtom(1).getProperty("AtomCounter"));
                }
            }
        }
        //Generate the murckoFragment
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(true,1);
        tmpMurckoFragmenter.setComputeRingFragments(false);
        tmpMurckoFragmenter.generateFragments(tmpMolecule);
        IAtomContainer[] tmpFrameworks = tmpMurckoFragmenter.getFrameworksAsContainers();
        if (tmpFrameworks.length == 0){
            return new AtomContainer(); //Returns an empty molecule because the original molecule does not contain a Schuffenhauer scaffold
        }
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
        //Generate SchuffenhauerScaffold
        for(IAtom tmpAtom : tmpLongFragment.atoms()){
            if(tmpAddAtomList.contains(tmpAtom.getProperty("AtomCounter"))){
                tmpLongFragment.addAtom(new Atom("O")); //Add an O for each remaining C with double bounded O
                for(IAtom tmpNewAtom : tmpLongFragment.atoms()){
                    if(tmpNewAtom.getProperty("AtomCounter") == null){ //The Atom without an AtomCounter is the new added O
                        tmpCounter++;
                        tmpNewAtom.setProperty("AtomCounter", tmpCounter);
                        tmpLongFragment.addBond(tmpNewAtom.getIndex(),tmpAtom.getIndex(), IBond.Order.DOUBLE); //Bond the new O and the C from list
                    }
                }
            }
        }
        //Add back hydrogens removed by the MurckoFragmenter
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpLongFragment);
        CDKHydrogenAdder.getInstance(tmpLongFragment.getBuilder()).addImplicitHydrogens(tmpLongFragment);
        return tmpLongFragment;
    }
}
