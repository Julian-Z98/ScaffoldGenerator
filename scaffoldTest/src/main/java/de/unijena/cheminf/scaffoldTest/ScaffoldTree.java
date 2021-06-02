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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


/**
 * Top-level class to organise the TreeNodes
 */
public class ScaffoldTree {

    /**
     * Saves all TreeNodes and numbers them in ascending order. Starts at 0. Key:Number, Value:TreeNode
     */
    private HashMap<Integer, TreeNode> nodeMap;
    /**
     * Saves all TreeNodes and numbers them in ascending order. Starts at 0. nodeMap with key and value swapped. Key:TreeNode, Value:Number
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
     * Generator for the creation of unique SMILES
     */
    private SmilesGenerator smilesGenerator;

    /**
     * Constructor
     */
    public ScaffoldTree() {
        this.nodeMap = new HashMap<Integer, TreeNode>();
        this.reverseNodeMap = new HashMap<TreeNode, Integer>();
        this.smilesMap = new HashMap<String, TreeNode>();
        this.levelMap = ArrayListMultimap.create();
        this.smilesGenerator = new SmilesGenerator(SmiFlavor.Unique);
    }

    /**
     * Add TreeNode to the ScaffoldTree
     * @param tmpNode TreeNode to be added
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public void addNode(TreeNode tmpNode) throws CDKException {
        //Add to nodeMap
        this.nodeMap.put(this.nodeMap.size(), tmpNode);
        //Add to reverseNodeMap
        this.reverseNodeMap.put(tmpNode, this.reverseNodeMap.size());
        //Add to smilesMap
        IAtomContainer tmpMolecule = (IAtomContainer) tmpNode.getData();
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
        this.smilesMap.put(tmpSmiles, tmpNode);
        //Add to levelMap
        this.levelMap.put(tmpNode.getLevel(), tmpNode);
    }

    /**
     * Checks whether the molecule is already present in the ScaffoldTree
     * @param tmpMolecule Molecule to check
     * @return Whether the molecule is located in the ScaffoldTree
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public boolean isMoleculeInTree(IAtomContainer tmpMolecule) throws CDKException {
        //Generate SMILES
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule);
        //Check in smilesMap
        if(this.smilesMap.containsKey(tmpSmiles)) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Return the TreeNode that belongs to a specific molecule
     * @param tmpMolecule molecule that is being searched for
     * @return TreeNode of the searched molecule
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public TreeNode getTreeNode(IAtomContainer tmpMolecule) throws CDKException {
        return this.smilesMap.get(this.smilesGenerator.create(tmpMolecule));
    }

    /**
     * Returns all TreeNodes of unique molecules. If a molecule occurs more than once in the tree, only one corresponding TreeNode is returned.
     * @return all TreeNodes of unique molecules
     */
    public List<TreeNode> getUniqueTreeNodes() {
        List<TreeNode> tmpNodeList = new ArrayList<>();
        tmpNodeList.addAll(this.smilesMap.values());
        return tmpNodeList;
    }

    /**
     * Return all TreeNodes that are at a certain level in the tree.
     * @param tmpLevel Level whose TreeNodes are to be returned
     * @return TreeNodes that are at a certain level
     */
    public List<TreeNode> getLevel(int tmpLevel) {
        return this.levelMap.get(tmpLevel);
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
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1. All others are marked with 0.
     * Each row and column number in the matrix is assigned to a node. The assignment can be requested with getMatrixNodes/getMatrixNode.
     * @return the adjascence matrix
     */
    public Integer[][] getTreeAsMatrix() {
        int tmpSize = this.nodeMap.size();
        Integer[][] tmpMatrix = new Integer[tmpSize][tmpSize];
        //Set all values of the matrix to 0
        for(int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            for(int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) {
                tmpMatrix[tmpRow][tmpCol] = 0;
            }
        }
        //Insert a 1 for each parent node
        int tmpCounter = 0;
        for(TreeNode tmpNode : this.nodeMap.values()) {
            //The root has no parent and is therefore skipped
            if(tmpCounter == 0) {
                tmpCounter++;
                continue;
            }
            //Set a 1 at the level of the parent and at the level of the node
            tmpMatrix[tmpCounter][this.reverseNodeMap.get(tmpNode.getParent())] = 1;
            tmpCounter++;
        }
        return tmpMatrix;
    }

    /**
     * Returns the nodes of the matrix in ascending order.
     * @return List of nodes in the order they appear in the matrix
     */
    public List<TreeNode> getMatrixNodes() {
        List<TreeNode> tmpList = new ArrayList<>();
        tmpList.addAll(this.nodeMap.values());
        return tmpList;
    }

    /**
     * Returns the node that belongs to a certain row and column number of the matrix.
     * @param tmpNumber Row and column number whose node is requested
     * @return Node that belongs to the requested column and row number
     */
    public TreeNode getMatrixNode(int tmpNumber) {
        return this.nodeMap.get(tmpNumber);
    }
}
