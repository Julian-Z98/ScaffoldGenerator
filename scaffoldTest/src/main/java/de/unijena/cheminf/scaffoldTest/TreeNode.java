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

import java.util.*;

/**
 * The TreeNodes are nodes from which a tree can be built.
 * It is used to organise the IAtomContainers and enables a relationship between the different objects.
 * As MoleculeType, any data type can be defined. In our scenario, the nodes contain molecules.
 *
 * Inspired by: https://github.com/gt4dev/yet-another-tree-structure
 */
public class TreeNode<MoleculeType> implements Iterable<TreeNode<MoleculeType>> {

    /**
     * Molecule that can be stored in each node
     */
    private MoleculeType molecule;

    /**
     * List of SMILES of the molecules from which this fragment originates
     */
    private ArrayList<String> OriginSmilesList;

    /**
     * Parent of the node
     */
    private TreeNode<MoleculeType> parent;

    /**
     * Child of the Node
     */
    private List<TreeNode<MoleculeType>> children;

    /**
     * List of indices of all elements
     */
    private List<TreeNode<MoleculeType>> elementsIndex;

    /**
     * Shows if the node has parents
     * @return Whether the node has parents
     */
    public boolean isOrphan() {
        return parent == null;
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
     * Creates a TreeNode
     * @param aMolecule molecule of the TreeNode
     */
    public TreeNode(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
        this.children = new LinkedList<TreeNode<MoleculeType>>();
        this.OriginSmilesList = new ArrayList<String>();
        this.elementsIndex = new LinkedList<TreeNode<MoleculeType>>();
        this.elementsIndex.add(this);
    }

    /**
     * Adds a child to the TreeNode, i.e. links it to a TreeNode on the level below
     * @param aMolecule Molecule of the child leave
     * @return Node of the child leave
     */
    public TreeNode<MoleculeType> addChild(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        TreeNode<MoleculeType> tmpChildNode = new TreeNode<MoleculeType>(aMolecule);
        tmpChildNode.parent = this;
        this.children.add(tmpChildNode);
        this.registerChildForSearch(tmpChildNode);
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
     * Outputs the level on which the node is located in the entire tree
     * @return level of the node in the entire tree
     */
    public int getLevel() {
        if (this.isOrphan())
            return 0;
        else
            return parent.getLevel() + 1;
    }

    /**
     * Child is registered for the index search.
     * The node can only be found with findTreeNode if it was previously registered in this way.
     * @param aNode Node to be registered
     */
    private void registerChildForSearch(TreeNode<MoleculeType> aNode) {
        Objects.requireNonNull(aNode, "Given Tree Node is 'null'");
        elementsIndex.add(aNode);
        if (parent != null)
            parent.registerChildForSearch(aNode);
    }

    /**
     * Find a specific TreeNode in the tree
     * Node must have been previously registered via registerChildForSearch
     * @param aCompMolecule Molecule on the basis of which the node is to be searched for
     * @return Element of the searched TreeNode
     */
    public TreeNode<MoleculeType> findTreeNode(Comparable<MoleculeType> aCompMolecule) {
        Objects.requireNonNull(aCompMolecule, "Given comparable Molecule is 'null'");
        for (TreeNode<MoleculeType> tmpElement : this.elementsIndex) {
            MoleculeType tmpMolecule = tmpElement.molecule;
            if (aCompMolecule.compareTo(tmpMolecule) == 0)
                return tmpElement;
        }
        return null;
    }

    /**
     * Iterator to iterate through the entire tree
     * @return the iterator
     */
    @Override
    public Iterator<TreeNode<MoleculeType>> iterator() {
        TreeNodeIter<MoleculeType> tmpIter = new TreeNodeIter<MoleculeType>(this);
        return tmpIter;
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
     * Get the parent node.
     * @return parent node
     */
    public TreeNode<MoleculeType> getParent() {
        return this.parent;
    }

    /**
     * Set the parent node.
     * @param aParent parent that are set
     */
    public void setParent(TreeNode<MoleculeType> aParent) {
        Objects.requireNonNull(aParent, "Given TreeNode is 'null'");
        this.parent = aParent;
    }

    /**
     * Get the children nodes.
     * @return children nodes
     */
    public List<TreeNode<MoleculeType>> getChildren() {
        return this.children;
    }

    /**
     * Set the children node.
     * @param aChildren children that are set
     */
    public void setChildren(List<TreeNode<MoleculeType>> aChildren) {
        Objects.requireNonNull(aChildren, "Given TreeNode List is 'null'");
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
    //</editor-fold>
}