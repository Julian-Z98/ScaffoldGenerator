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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.*;


/**
 * Top-level class to organise the TreeNodes
 */
public class ScaffoldTree {

    /**
     * Saves all TreeNodes and numbers them in ascending order. Starts at 0. Key:Number, Value:TreeNode
     */
    private HashMap<Integer, TreeNode> nodeMap;
    /**
     * Saves all TreeNodes and numbers them in ascending order. Starts at 0.
     * nodeMap with key and value swapped. Key:TreeNode, Value:Number
     */
    private HashMap<TreeNode, Integer> reverseNodeMap;
    /**
     * Saves all TreeNodes according to their Unique SMILES. Key:SMILES, Value:TreeNode
     */
    private HashMap<String, TreeNode> smilesMap;
    /**
     * Saves all TreeNodes according to their level. Key:Level, Value:TreeNode
     */
    private ListMultimap<Integer, TreeNode> levelMap;
    /**
     * Generator for the creation of Unique SMILES
     * Unique SMILES: canonical SMILES string, different atom ordering produces the same (apart from rare exceptions) SMILES.
     *                No isotope or stereochemistry encoded.
     */
    private SmilesGenerator smilesGenerator;
    /**
     * Shows how many nodes have been added so far. Removing nodes has no effect on it.
     */
    private int nodeCounter;

    /**
     * Constructor
     */
    public ScaffoldTree() {
        this.nodeMap = new HashMap<Integer, TreeNode>();
        this.reverseNodeMap = new HashMap<TreeNode, Integer>();
        this.smilesMap = new HashMap<String, TreeNode>();
        this.levelMap = ArrayListMultimap.create();
        this.smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
        this.nodeCounter = 0;
    }

    /**
     * Add TreeNode to the ScaffoldTree
     * @param aNode TreeNode to be added
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public void addNode(TreeNode aNode) throws CDKException {
        Objects.requireNonNull(aNode, "Given TreeNode is 'null'");
        Objects.requireNonNull(aNode.getData(), "Given Data is 'null'");
        //Add to nodeMap
        this.nodeMap.put(this.nodeCounter, aNode);
        //Add to reverseNodeMap
        this.reverseNodeMap.put(aNode, this.nodeCounter);
        /*Add to smilesMap*/
        IAtomContainer tmpMolecule = (IAtomContainer) aNode.getData();
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
        this.smilesMap.put(tmpSmiles, aNode);
        //Add to levelMap
        this.levelMap.put(aNode.getLevel(), aNode);
        //Increase nodeCounter
        this.nodeCounter++;
    }

    public void removeNode(TreeNode aNode) throws CDKException {
        Objects.requireNonNull(aNode, "Given TreeNode is 'null'");
        if(!this.reverseNodeMap.containsKey(aNode)) { //Check if the node exists in the tree
            throw new IllegalArgumentException("Node is not in tree");
        }
        /*Remove from nodeMap and reverseNodeMap*/
        int tmpNumberInNodeMap = this.reverseNodeMap.get(aNode); //get number in nodeMap
        this.nodeMap.remove(tmpNumberInNodeMap);
        this.reverseNodeMap.remove(aNode);
        /*Remove from smilesMap*/
        this.smilesMap.clear();
        for(TreeNode tmpTreeNode : this.nodeMap.values()) {
            IAtomContainer tmpMolecule = (IAtomContainer) tmpTreeNode.getData();
            String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
            this.smilesMap.put(tmpSmiles , tmpTreeNode);
        }
        /*Remove from levelMap*/
        this.levelMap.clear();
        for(TreeNode tmpTreeNode : this.nodeMap.values()) {
            this.levelMap.put(tmpTreeNode.getLevel(), tmpTreeNode);
        }
    }

    /**
     * Checks whether the molecule is already present in the ScaffoldTree
     * Check whether it is the same molecule using the Unique SMILES.
     * Unique SMILES: canonical SMILES string, different atom ordering produces the same (apart from rare exceptions) SMILES.
     *                No isotope or stereochemistry encoded.
     * @param aMolecule Molecule to check
     * @return Whether the molecule is located in the ScaffoldTree
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public boolean isMoleculeInTree(IAtomContainer aMolecule) throws CDKException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'");
        //Generate SMILES
        String tmpSmiles = this.smilesGenerator.create(aMolecule);
        /*Check in smilesMap*/
        if(this.smilesMap.containsKey(tmpSmiles)) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Return the TreeNode that belongs to a specific molecule.
     * Check whether it is the same molecule using the Unique SMILES.
     * Unique SMILES: canonical SMILES string, different atom ordering produces the same (apart from rare exceptions) SMILES.
     *                No isotope or stereochemistry encoded.
     * @param aMolecule molecule that is being searched for
     * @return TreeNode of the searched molecule
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public TreeNode getTreeNode(IAtomContainer aMolecule) throws CDKException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'");
        if(!this.isMoleculeInTree(aMolecule)) { //Check if the molecule exists in the tree
            throw new IllegalArgumentException("Molecule is not in tree");
        }
        return this.smilesMap.get(this.smilesGenerator.create(aMolecule));
    }

    /**
     * Returns all TreeNodes of unique molecules.
     * If a molecule occurs more than once in the tree, only one corresponding TreeNode is returned.
     * @return all TreeNodes of unique molecules
     */
    public List<TreeNode> getUniqueTreeNodes() {
        List<TreeNode> tmpNodeList = new ArrayList<>();
        tmpNodeList.addAll(this.smilesMap.values());
        return tmpNodeList;
    }

    /**
     * Outputs the maximum level of the tree
     * @return the maximum level of the tree
     */
    public int getMaxLevel() {
        List<Integer> tmpLevelList = new ArrayList<>();
        tmpLevelList.addAll(this.levelMap.keys());
        return Collections.max(tmpLevelList);
    }

    /**
     * Return all TreeNodes that are at a certain level in the tree.
     * @param aLevel Level whose TreeNodes are to be returned
     * @return TreeNodes that are at a certain level
     */
    public List<TreeNode> getAllNodesOnLevel(int aLevel) {
        Objects.requireNonNull(aLevel, "Given number is 'null'");
        if(this.getMaxLevel() >= aLevel) { //Level must be less than or equal to the maximum level
            return this.levelMap.get(aLevel);
        }
        throw new IllegalArgumentException("Level does not exist: " + aLevel);
    }

    /**
     * Get all TreeNodes of the ScaffoldTree
     * @return all TreeNodes of the ScaffoldTree
     */
    public List<TreeNode> getAllNodes() {
        List<TreeNode> tmpList = new ArrayList<>();
        tmpList.addAll(this.levelMap.values());
        return tmpList;
    }

    /**
     * Checks if the parent of each node in the tree is present (root excluded).
     * If this is not the case, the tree is disjointed. Then it is actually more than one tree.
     * @return whether the tree is connected
     */
    public boolean isTreeConnected() {
        int tmpCounter = 0;
        boolean tmpIsTreeConnected = true;
        for(TreeNode tmpNode : this.nodeMap.values()) {
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
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1.
     * All others are marked with 0. Each row and column number in the matrix is assigned to a node.
     * The assignment can be requested with getMatrixNodes/getMatrixNode.
     * only works with connected trees. Can be checked with isTreeConnected.
     * @return the adjacency matrix
     */
    public Integer[][] getTreeAsMatrix() {
        int tmpSize = this.nodeMap.size();
        Integer[][] tmpMatrix = new Integer[tmpSize][tmpSize];
        if(this.isTreeConnected()) { //Only a connected matrix is calculated
            /*Set all values of the matrix to 0*/
            for (int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
                for (int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) {
                    tmpMatrix[tmpRow][tmpCol] = 0;
                }
            }
            /*Insert a 1 for each parent node*/
            int tmpCounter = 0;
            for (TreeNode tmpNode : this.nodeMap.values()) {
                if (tmpNode.getParent() != null) {
                    //Set a 1 at the level of the parent and at the level of the node
                    tmpMatrix[tmpCounter][this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpNode.getParent()))] = 1;
                    //Set a 1 at the level of the node and at the level of the parent
                    tmpMatrix[this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpNode.getParent()))][tmpCounter] = 1;
                }
                tmpCounter++;
            }
        }
        return tmpMatrix;
    }

    /**
     * Returns the number of the nodes and the nodes of the matrix in ascending order.
     * @return HashMap with the number of the nodes and the nodes in the order they appear in the matrix
     */
    public HashMap<Integer, TreeNode> getMatrixNodes() {
        return this.nodeMap;
    }

    /**
     * Gives the number of nodes as they occur in the matrix. Missing numbers have been removed.
     * @return List with all node numbers of the matrix
     */
    public List<Integer> getMatrixNodesNumbers() {
        List<Integer> tmpList = new ArrayList<>();
        tmpList.addAll(this.nodeMap.keySet());
        return tmpList;
    }

    /**
     * Indicates whether there is a node with this number in the tree.
     * @param aNumber Number of the node to be checked
     * @return Is the node in the tree
     */
    public boolean isNodeInTree(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.containsKey(aNumber);
    }

    /**
     * Returns the node that belongs to a certain row and column number of the matrix.
     * Returns zero if the node is not in the tree. Throwing an exception makes it difficult to build a tree
     * @param aNumber Row and column number whose node is requested
     * @return Node that belongs to the requested column and row number
     */
    public TreeNode getMatrixNode(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.get(aNumber);
    }

    /**
     * Checks whether the tree has a single root. Returns false in any other case.
     * @return Whether the tree has only one root
     */
    public boolean hasOneSingleRootNode() {
        int tmpRootCounter = 0;
        for(TreeNode tmpNode : this.nodeMap.values()) {
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
     * Outputs root node of the tree.
     * @return root node of the tree
     */
    public TreeNode getRoot() {
        if(!this.hasOneSingleRootNode()) { //Checks whether the tree has a single root
            throw new IllegalArgumentException("Tree has no clear root");
        }
        for(TreeNode tmpNode : this.nodeMap.values()) {
            /*Is the parent of the node in the tree*/
            if(!this.reverseNodeMap.containsKey(tmpNode.getParent())) {
                return tmpNode; //Return the first node without parent in tree
            }
        }
        throw new IllegalArgumentException("Tree has no clear root");
    }
}
