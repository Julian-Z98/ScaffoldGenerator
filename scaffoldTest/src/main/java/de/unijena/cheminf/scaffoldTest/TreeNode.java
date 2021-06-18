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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;

/**
 * The TreeNodes are nodes from which a tree can be built.
 * It is used to organise the IAtomContainers and enables a relationship between the different objects.
 *
 * Inspired by: https://github.com/gt4dev/yet-another-tree-structure
 */
public class TreeNode<IAtomContainer> implements Iterable<TreeNode<IAtomContainer>> {

    /**
     * Molecule that can be stored in each node
     */
    private IAtomContainer molecule;

    /**
     * Parent of the node
     */
    private TreeNode<IAtomContainer> parent;

    /**
     * Child of the Node
     */
    private List<TreeNode<IAtomContainer>> children;

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
     * List of indices of all elements
     */
    private List<TreeNode<IAtomContainer>> elementsIndex;

    /**
     * Creates a TreeNode
     * @param aMolecule molecule of the TreeNode
     */
    public TreeNode(IAtomContainer aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
        this.children = new LinkedList<TreeNode<IAtomContainer>>();
        this.elementsIndex = new LinkedList<TreeNode<IAtomContainer>>();
        this.elementsIndex.add(this);
    }

    /**
     * Adds a child to the TreeNode, i.e. links it to a TreeNode on the level below
     * @param aMolecule Molecule of the child leave
     * @return Node of the child leave
     */
    public TreeNode<IAtomContainer> addChild(IAtomContainer aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        TreeNode<IAtomContainer> tmpChildNode = new TreeNode<IAtomContainer>(aMolecule);
        tmpChildNode.parent = this;
        this.children.add(tmpChildNode);
        this.registerChildForSearch(tmpChildNode);
        return tmpChildNode;
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
    private void registerChildForSearch(TreeNode<IAtomContainer> aNode) {
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
    public TreeNode<IAtomContainer> findTreeNode(Comparable<IAtomContainer> aCompMolecule) {
        Objects.requireNonNull(aCompMolecule, "Given comparable Molecule is 'null'");
        for (TreeNode<IAtomContainer> tmpElement : this.elementsIndex) {
            IAtomContainer tmpMolecule = tmpElement.molecule;
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
    public Iterator<TreeNode<IAtomContainer>> iterator() {
        TreeNodeIter<IAtomContainer> tmpIter = new TreeNodeIter<IAtomContainer>(this);
        return tmpIter;
    }

    //<editor-fold desc="get/set">
    /**
     * Get the node molecule.
     * @return node molecule
     */
    public IAtomContainer getMolecule() {
        return this.molecule;
    }

    /**
     * Set the node molecule.
     * @param aMolecule molecule that are set
     */
    public void setMolecule(IAtomContainer aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        this.molecule = aMolecule;
    }

    /**
     * Get the parent node.
     * @return parent node
     */
    public TreeNode<IAtomContainer> getParent() {
        return this.parent;
    }

    /**
     * Set the parent node.
     * @param aParent parent that are set
     */
    public void setParent(TreeNode<IAtomContainer> aParent) {
        Objects.requireNonNull(aParent, "Given TreeNode is 'null'");
        this.parent = aParent;
    }

    /**
     * Get the children nodes.
     * @return children nodes
     */
    public List<TreeNode<IAtomContainer>> getChildren() {
        return this.children;
    }

    /**
     * Set the children node.
     * @param aChildren children that are set
     */
    public void setChildren(List<TreeNode<IAtomContainer>> aChildren) {
        Objects.requireNonNull(aChildren, "Given TreeNode List is 'null'");
        this.children = aChildren;
    }
    //</editor-fold>
}