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

package de.unijena.cheminf.scaffolds;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;

/**
 * Base class of node objects.
 * Contains the basic functionality of nodes needed for a tree or a network.
 *
 * @param <MoleculeType> As MoleculeType, any data type can be defined.
 *                      In our scenario, the nodes contain molecules.
 *
 * Inspired by: https://github.com/gt4dev/yet-another-tree-structure
 */
public abstract class ScaffoldNodeBase<MoleculeType> {

    /**
     * Molecule that can be stored in each node
     */
    protected MoleculeType molecule;

    /**
     * List of SMILES of the molecules from which this fragment originates.
     *
     * If additional information of the origin is needed,
     * it can be stored in a matrix with the IAtomContainer. The SMILES stored here can then be used as a key.
     */
    protected ArrayList<String> OriginSmilesList;

    /**
     * Child of the Node
     */
    protected List<ScaffoldNodeBase<MoleculeType>> children;

    /**
     * List of indices of all elements
     */
    protected List<ScaffoldNodeBase<MoleculeType>> elementsIndex;

    /**
     * Indicates whether it has children
     * @return true if it has no children
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    /**
     * Shows if the node has parents
     * @return Whether the node has parents
     */
    public abstract boolean isOrphan();

    /**
     * Constructor
     * @param aMolecule molecule of the ScaffoldNodeBase
     */
    public ScaffoldNodeBase(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
        this.children = new LinkedList<ScaffoldNodeBase<MoleculeType>>();
        this.OriginSmilesList = new ArrayList<String>();
        this.elementsIndex = new LinkedList<ScaffoldNodeBase<MoleculeType>>();
        this.elementsIndex.add(this);
    }

    /**
     * Adds a child to the ScaffoldNodeBase, i.e. links it to a ScaffoldNodeBase on the level below
     * @param aMolecule Molecule of the child leave
     * @return Node of the child leave
     */
    public abstract ScaffoldNodeBase<MoleculeType> addChild(MoleculeType aMolecule);

    /**
     * Adds another string to the OriginSmilesList if it is not already present.
     * @param aString String to be added
     */
    public void addOriginSmiles(String aString) {
        Objects.requireNonNull(aString, "Given SMILES of the molecule is 'null'");
        if(!this.OriginSmilesList.contains(aString)) {
            OriginSmilesList.add(aString);
        }
    }

    //<editor-fold desc="get/set">
    /**
     * Outputs the level on which the node is located in the entire tree
     * @return level of the node in the entire tree
     */
    public abstract int getLevel();

    /**
     * Get the node molecule.
     * @return node molecule
     */
    public MoleculeType getMolecule() {
        return this.molecule;
    }

    /**
     * Set the node molecule.
     * @param aMolecule molecule that are set
     */
    public void setMolecule(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
    }

    /**
     * Get the children nodes.
     * @return children nodes
     */
    public List<ScaffoldNodeBase<MoleculeType>> getChildren() {
        return this.children;
    }

    /**
     * Set the children node.
     * @param aChildren children that are set
     */
    public void setChildren(List<ScaffoldNodeBase<MoleculeType>> aChildren) {
        Objects.requireNonNull(aChildren, "Given ScaffoldNodeBase List is 'null'");
        this.children = aChildren;
    }

    /**
     * Get the OriginSmilesList
     * @return List of SMILES of the molecules from which this fragment originates
     */
    public ArrayList<String> getOriginSmilesList() {
        return this.OriginSmilesList;
    }

    /**
     * Get the size of the OriginSmilesList
     * @return size of the OriginSmilesList
     */
    public Integer getOriginCount() {
        return this.OriginSmilesList.size();
    }
    /**
     * Set the entire OriginSmilesList
     * @param aOriginSmilesList SMILES of molecules that are set
     */
    public void setOriginSmilesList(ArrayList<String> aOriginSmilesList) {
        Objects.requireNonNull(aOriginSmilesList, "Given SMILES of the molecule List is 'null'");
        this.OriginSmilesList = aOriginSmilesList;
    }
    //</editor-fold>
}
