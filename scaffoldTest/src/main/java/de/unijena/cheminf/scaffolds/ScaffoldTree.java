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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.Objects;


/**
 * Top-level class to organise the TreeNodes
 */
public class ScaffoldTree extends ScaffoldNodeCollectionBase {

    /**
     * Constructor
     * @param aSmilesGenerator Used SMILES Generator
     */
    public ScaffoldTree(SmilesGenerator aSmilesGenerator) {
        super(aSmilesGenerator);
    }

    /**
     * Default Constructor
     */
    public ScaffoldTree() {
        this(new SmilesGenerator(SmiFlavor.Unique));
    }

    /**
     * Add TreeNode to the ScaffoldTree
     * @param aNode TreeNode to be added
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public void addNode(TreeNode aNode) throws CDKException {
        Objects.requireNonNull(aNode, "Given TreeNode is 'null'");
        Objects.requireNonNull(aNode.getMolecule(), "Given Data is 'null'");
        //Add to nodeMap
        this.nodeMap.put(this.nodeCounter, aNode);
        //Add to reverseNodeMap
        this.reverseNodeMap.put(aNode, this.nodeCounter);
        /*Add to smilesMap*/
        IAtomContainer tmpMolecule = (IAtomContainer) aNode.getMolecule();
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
        this.smilesMap.put(tmpSmiles, aNode);
        //Add to levelMap
        this.levelMap.put(aNode.getLevel(), aNode);
        /*Add origins from the new node to parent nodes*/
        TreeNode tmpNode = aNode;
        for(int tmpCount = 0; tmpCount < aNode.getLevel(); tmpCount++) {
            TreeNode tmpNextNode = tmpNode.getParent();
            for(Object tmpString : tmpNode.getOriginSmilesList()) {
                tmpNextNode.addOriginSmiles((String) tmpString);
            }
            tmpNode = tmpNextNode;
        }
        //Increase nodeCounter
        this.nodeCounter++;
    }

    /**
     * Removes a Node. This does not change the order. The numbering does not move up.
     * @param aNode Node to remove
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the Scaffold
     */
    public void removeNode(TreeNode aNode) throws CDKException, IllegalArgumentException {
        Objects.requireNonNull(aNode, "Given ScaffoldNode is 'null'");
        if(!this.reverseNodeMap.containsKey(aNode)) { //Check if the node exists in the Scaffold
            throw new IllegalArgumentException("Node is not in Scaffold");
        }
        /*Remove from nodeMap and reverseNodeMap*/
        int tmpNumberInNodeMap = this.reverseNodeMap.get(aNode); //get number in nodeMap
        this.nodeMap.remove(tmpNumberInNodeMap);
        this.reverseNodeMap.remove(aNode);
        /*Remove from smilesMap*/
        String tmpSmiles = this.smilesGenerator.create((IAtomContainer) aNode.getMolecule()); //Convert molecule to SMILES
        this.smilesMap.remove(tmpSmiles, aNode);
        /*Remove from levelMap*/
        levelMap.remove(aNode.getLevel(), aNode);
    }

    /**
     * Checks if the parent of each node in the tree is present (root excluded).
     * If this is not the case, the tree is disjointed. Then it is actually more than one tree.
     * @return whether the tree is connected
     */
    public boolean isTreeConnected() {
        int tmpCounter = 0;
        boolean tmpIsTreeConnected = true;
        for(ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            TreeNode tmpNode = (TreeNode) tmpNodeBase;
            /*The root has no parent and is therefore skipped*/
            if(tmpCounter == 0) {
                tmpCounter++;
                continue;
            }
            /*Is the parent of the node in the tree*/
            if(!this.reverseNodeMap.containsKey(tmpNode.getParent())) {
                tmpIsTreeConnected = false;
            }
            tmpCounter++;
        }
        return tmpIsTreeConnected;
    }

    /**
     * Checks whether the tree has a single root. Returns false in any other case.
     * @return Whether the tree has only one root
     */
    public boolean hasOneSingleRootNode() {
        int tmpRootCounter = 0;
        for(ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            TreeNode tmpNode = (TreeNode) tmpNodeBase;
            /*Is the parent of the node in the tree*/
            if(!this.reverseNodeMap.containsKey(tmpNode.getParent())) {
                tmpRootCounter++;
            }
            if(tmpRootCounter < 1) { //If the tree has more than one root
                return false;
            }
        }
        if(tmpRootCounter == 1) { //If the tree has only one root
            return true;
        } else { //If the tree has no root
            return false;
        }
    }

    /**
     * Adds another ScaffoldTree to the existing one if possible.
     * The new tree is inserted at the node that both trees have in common.
     * All children of the new tree at this node are taken over if they do not already exist.
     *
     * The new tree is simply taken over if the existing tree is empty.
     * If there is no match between the two trees, false is returned and the old tree is not changed.
     *
     * If a molecule does not generate a SchuffenhauerScaffold, it is stored as a node with empty SMILES and is treated normally.
     * All other empty nodes are then added to this tree accordingly.
     * @param aScaffoldTree tree to be inserted into the existing ScaffoldTree.
     * @return false if the new tree has nodes and cannot be inserted because the two trees have no common nodes.
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public boolean mergeTree(ScaffoldTree aScaffoldTree) throws CDKException {
        /*If the old ScaffoldTree is empty, transfer the new ScaffoldTree to be added.*/
        if(this.hasOneSingleRootNode() == false) {
            for(ScaffoldNodeBase tmpBaseNode : aScaffoldTree.getAllNodes()) {
                TreeNode tmpNode = (TreeNode) tmpBaseNode;
                this.addNode(tmpNode);
            }
            /*The new tree was inserted*/
            return true;
        }
        /*If the old Scaffold tree is not empty*/
        else {
            boolean tmpAreTreesOverlapping = true;
            /*Go through each level of the tree starting at the root*/
            for(int i = 0; i <= aScaffoldTree.getMaxLevel(); i++) {
                tmpAreTreesOverlapping = false;
                /*Compare all nodes of the old tree on this level with all nodes of the new tree on this level*/
                for(ScaffoldNodeBase tmpOldTreeNodeBase : this.getAllNodesOnLevel(i)) {
                    TreeNode tmpOldTreeNode = (TreeNode) tmpOldTreeNodeBase;
                    for(ScaffoldNodeBase tmpNewTreeNodeBase : aScaffoldTree.getAllNodesOnLevel(i)) {
                        TreeNode tmpNewTreeNode = (TreeNode) tmpNewTreeNodeBase;
                        /*Generate the corresponding SMILES of the molecules*/
                        String tmpOldSmiles = this.smilesGenerator.create((IAtomContainer) tmpOldTreeNode.getMolecule());
                        String tmpNewSmiles = this.smilesGenerator.create((IAtomContainer) tmpNewTreeNode.getMolecule());
                        /*Check whether a fragment occurs in both trees*/
                        if(tmpOldSmiles.equals(tmpNewSmiles)) {
                            /*Add the origin smiles to the OldSmilesTree fragment*/
                            for(Object tmpOriginSmiles : tmpNewTreeNode.getOriginSmilesList()) {
                                tmpOldTreeNode.addOriginSmiles((String) tmpOriginSmiles);
                            }
                            /*Trees are overlapping if a fragment occurs in both trees*/
                            tmpAreTreesOverlapping = true;
                            /*Get the children of the overlapping node*/
                            for(Object tmpNewChild : tmpNewTreeNode.getChildren()) {
                                TreeNode tmpNewChildNode = (TreeNode) tmpNewChild;
                                IAtomContainer tmpNewChildMolecule = (IAtomContainer) tmpNewChildNode.getMolecule();
                                /*Add the child if it is not already in the tree*/
                                if(!this.isMoleculeContained(tmpNewChildMolecule)) {
                                    int tmpChildrenNumber = tmpOldTreeNode.getChildren().size();
                                    tmpOldTreeNode.addChild(tmpNewChildMolecule);
                                    TreeNode tmpAddedNode = (TreeNode) tmpOldTreeNode.getChildren().get(tmpChildrenNumber);
                                    this.addNode(tmpAddedNode);
                                }
                            }
                        }
                    }
                }
                /*If there were no overlaps, there is no need to continue the search, as there will be no more in the future*/
                if(!tmpAreTreesOverlapping) {
                    /*If there was already no overlap at the root, the new tree could not be inserted*/
                    if(i == 0) {
                        /*The new tree was not inserted*/
                        return false;
                    }
                    /*End the search*/
                    break;
                }
            }
        }
        /*The new tree was inserted*/
        return true;
    }

    //<editor-fold desc="get/set">
    /**
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1.
     * All others are marked with 0. Each row and column number in the matrix is assigned to a node.
     * The assignment can be requested with getMatrixNodes/getMatrixNode.
     * only works with connected trees. Can be checked with isTreeConnected.
     * @return the adjacency matrix
     * @throws IllegalStateException if the tree is not connected
     */
    public Integer[][] getTreeAsMatrix() throws IllegalStateException {
        int tmpSize = this.nodeMap.size();
        Integer[][] tmpMatrix = new Integer[tmpSize][tmpSize];
        if(!this.isTreeConnected()) { //Only a connected matrix is calculated
            throw new IllegalStateException("Tree is not connected");
        }
        /*Set all values of the matrix to 0*/
        for (int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            for (int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) {
                tmpMatrix[tmpRow][tmpCol] = 0;
            }
        }
        /*Insert a 1 for each parent node*/
        int tmpCounter = 0;
        for (ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            TreeNode tmpNode = (TreeNode) tmpNodeBase;
            if (tmpNode.getParent() != null) {
                //Set a 1 at the level of the parent and at the level of the node
                tmpMatrix[tmpCounter][this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpNode.getParent()))] = 1;
                //Set a 1 at the level of the node and at the level of the parent
                tmpMatrix[this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpNode.getParent()))][tmpCounter] = 1;
            }
            tmpCounter++;
        }
        return tmpMatrix;
    }

    /**
     * Outputs root node of the tree.
     * @return root node of the tree
     * @throws IllegalStateException if the tree has no clear root
     */
    public TreeNode getRoot()  throws IllegalStateException {
        if(!this.hasOneSingleRootNode()) { //Checks whether the tree has a single root
            throw new IllegalStateException("Tree has no clear root");
        }
        for(ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            TreeNode tmpNode = (TreeNode) tmpNodeBase;
            /*Is the parent of the node in the tree*/
            if(!this.reverseNodeMap.containsKey(tmpNode.getParent())) {
                return tmpNode; //Return the first node without parent in tree
            }
        }
        throw new IllegalStateException("Tree has no clear root");
    }
    //</editor-fold>
}
