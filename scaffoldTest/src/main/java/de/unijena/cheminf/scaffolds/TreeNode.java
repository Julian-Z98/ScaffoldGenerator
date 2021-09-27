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

import java.util.Objects;

/**
 * The TreeNodes are nodes from which a tree can be built.
 * It is used to organise the IAtomContainers and enables a relationship between the different objects.
 * A TreeNode can have different children but only one parent.
 *
 * @param <MoleculeType> As MoleculeType, any data type can be defined.
 *                      In our scenario, the nodes contain molecules.
 */
public class TreeNode<MoleculeType> extends ScaffoldNodeBase<MoleculeType> {

    /**
     * Parent of the node
     */
    private TreeNode<MoleculeType> parent;

    /**
     * Shows if the node has parents
     * @return Whether the node has parents
     */
    public boolean isOrphan() {
        return parent == null;
    }

    /**
     * Constructor
     * @param aMolecule molecule of the TreeNode
     */
    public TreeNode(MoleculeType aMolecule) {
        super(aMolecule);
    }

    /**
     * Adds a child to the TreeNode, i.e. links it to a TreeNode on the level below
     * @param aMolecule Molecule of the child leave
     * @return Node of the child leave
     */
    public TreeNode<MoleculeType> addChild(MoleculeType aMolecule) {
        Objects.requireNonNull(aMolecule, "Given molecule is 'null'");
        TreeNode<MoleculeType> tmpChildNode = new TreeNode<MoleculeType>(aMolecule);
        this.children.add(tmpChildNode);
        tmpChildNode.parent = this;
        return tmpChildNode;
    }

    //<editor-fold desc="get/set">
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
    //</editor-fold>
}