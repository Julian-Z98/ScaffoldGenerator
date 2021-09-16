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
 * The NetworkNodes are nodes from which a Network can be built.
 * It is used to organise the IAtomContainers and enables a relationship between the different objects.
 * As MoleculeType, any data type can be defined. In our scenario, the nodes contain molecules.
 *
 * Inspired by: https://github.com/gt4dev/yet-another-tree-structure
 */
public class NetworkNode <MoleculeType> {

    /**
     * Molecule that can be stored in each node
     */
    private MoleculeType molecule;

    /**
     * List of SMILES of the molecules from which this fragment originates.
     *
     * If additional information of the origin is needed,
     * it can be stored in a matrix with the IAtomContainer. The SMILES stored here can then be used as a key.
     */
    private ArrayList<String> OriginSmilesList;

    /**
     * parents of the node
     */
    private List<NetworkNode<MoleculeType>> parents;

    /**
     * Children of the Node
     */
    private List<NetworkNode<MoleculeType>> children;

    /**
     * List of indices of all elements
     */
    private List<NetworkNode<MoleculeType>> elementsIndex;

    /**
     * Shows if the node has parents
     * @return Whether the node has parents
     */
    public boolean isOrphan() {
        return parents.isEmpty();
    }

    /**
     * Indicates whether it is a leaf, i.e. whether there are deeper nodes
     * True if there are no deeper nodes
     * @return Whether it is a leaf
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    /**
     * Creates a NetworkNode
     * @param aMolecule molecule of the NetworkNode
     */
    public NetworkNode(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
        this.parents =  new ArrayList<NetworkNode<MoleculeType>>();
        this.children = new LinkedList<NetworkNode<MoleculeType>>();
        this.OriginSmilesList = new ArrayList<String>();
        this.elementsIndex = new LinkedList<NetworkNode<MoleculeType>>();
        this.elementsIndex.add(this);
    }

    /**
     * Adds a child to the NetworkNode, i.e. links it to a NetworkNode on the level below
     * @param aMolecule Molecule of the child leave
     * @return Node of the child leave
     */
    public NetworkNode<MoleculeType> addChild(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        NetworkNode<MoleculeType> tmpChildNode = new NetworkNode<MoleculeType>(aMolecule);
        this.children.add(tmpChildNode);
        return tmpChildNode;
    }

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
    /**
     * Outputs the level on which the node is located in the entire network
     * @return level of the node in the entire network
     */
    public int getLevel() {
        if (this.isOrphan())
            return 0;
        else
            return parents.get(0).getLevel() + 1;
    }

    //<editor-fold desc="get/set">
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
     * Get the parents node.
     * @return parents node
     */
    public List<NetworkNode<MoleculeType>> getParents() {
        return this.parents;
    }

    /**
     * Set the parents node.
     * @param aParents parents that are set
     */
    public void setParents(List<NetworkNode<MoleculeType>> aParents) {
        Objects.requireNonNull(aParents, "Given NetworkNode is 'null'");
        this.parents = aParents;
    }

    /**
     * Add the parents node and add this node as child to the parent node if not already done.
     * @param aParent parent that are added
     */
    public void addParent(NetworkNode<MoleculeType> aParent) {
        Objects.requireNonNull(aParent, "Given NetworkNode is 'null'");
        /*Add child if not already added*/
        boolean tmpIsAlreadyChild = false;
        for(NetworkNode tmpNode : aParent.getChildren()) {
            if(tmpNode.getMolecule() == this.getMolecule()){
                tmpIsAlreadyChild = true;
            }
        }
        if(tmpIsAlreadyChild == false) {
            aParent.addChild(this.getMolecule());
        }
        //Add parent
        this.parents.add(aParent);
    }

    /**
     * Get the children nodes.
     * @return children nodes
     */
    public List<NetworkNode<MoleculeType>> getChildren() {
        return this.children;
    }

    /**
     * Set the children node.
     * @param aChildren children that are set
     */
    public void setChildren(List<NetworkNode<MoleculeType>> aChildren) {
        Objects.requireNonNull(aChildren, "Given NetworkNode List is 'null'");
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
     * Set the entire OriginSmilesList
     * @param aOriginSmilesList SMILES of molecules that are set
     */
    public void setOriginSmilesList(ArrayList<String> aOriginSmilesList) {
        Objects.requireNonNull(aOriginSmilesList, "Given SMILES of the molecule List is 'null'");
        this.OriginSmilesList = aOriginSmilesList;
    }
}
