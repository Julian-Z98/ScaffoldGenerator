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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.*;

/**
 * Base class of node collection objects.
 * Top-level class to organise the ScaffoldNodeBase objects.
 *
 */
public abstract class ScaffoldNodeCollectionBase {

    /**
     * Saves all ScaffoldNodes and numbers them in ascending order. Starts at 0. Key:Number, Value:ScaffoldNode
     */
    protected HashMap<Integer, ScaffoldNodeBase> nodeMap;

    /**
     * Saves all ScaffoldNodes and numbers them in ascending order. Starts at 0.
     * nodeMap with key and value swapped. Key:ScaffoldNode, Value:Number
     */
    protected HashMap<ScaffoldNodeBase, Integer> reverseNodeMap;

    /**
     * Saves all ScaffoldNodes according to their SMILES. Key:SMILES, Value:ScaffoldNode
     */
    protected ListMultimap<String, ScaffoldNodeBase> smilesMap;

    /**
     * Saves all ScaffoldNodes according to their level. Key:Level, Value:ScaffoldNode
     */
    protected ListMultimap<Integer, ScaffoldNodeBase> levelMap;

    /**
     * Generator for the creation of SMILES
     */
    protected SmilesGenerator smilesGenerator;

    /**
     * Shows how many nodes have been added so far. Removing nodes has no effect on it.
     */
    protected int nodeCounter;

    /**
     * Constructor
     * @param aSmilesGenerator Used SMILES Generator
     */
    public ScaffoldNodeCollectionBase(SmilesGenerator aSmilesGenerator) {
        this.nodeMap = new HashMap<Integer, ScaffoldNodeBase>();
        this.reverseNodeMap = new HashMap<ScaffoldNodeBase, Integer>();
        this.smilesMap = ArrayListMultimap.create();
        this.levelMap = ArrayListMultimap.create();
        this.smilesGenerator = aSmilesGenerator;
        this.nodeCounter = 0;
    }

    /**
     * Default Constructor
     */
    public ScaffoldNodeCollectionBase() {
        this(new SmilesGenerator(SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols));
    }

    /**
     * Add ScaffoldNodeBase to the ScaffoldNodeCollectionBase
     * @param aNode Node to be added. Must match the tree object. For example, a ScaffoldTree requires a TreeNode.
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public abstract void addNode(ScaffoldNodeBase aNode) throws CDKException;

    /**
     * Removes a Node. This does not change the order. The numbering does not move up.
     * @param aNode Node to remove. Must match the tree object. For example, a ScaffoldTree requires a TreeNode.
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the Scaffold
     */
    public abstract void removeNode(ScaffoldNodeBase aNode) throws CDKException, IllegalArgumentException;

    /**
     * Checks whether the molecule is already present in the Scaffold Collection
     * Check whether it is the same molecule using the SMILES string.
     * @param aMolecule Molecule to check
     * @return Whether the molecule is located in the Scaffold Collection
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public boolean containsMolecule(IAtomContainer aMolecule) throws CDKException {
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
     * Indicates whether there is a node with this number in the collection.
     * @param aNumber Number of the node to be checked
     * @return Is the node in the collection
     */
    public boolean isNodeContained(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.containsKey(aNumber);
    }

    //<editor-fold desc="get/set">
    /**
     * Return the ScaffoldNode that belongs to a specific molecule.
     * Check whether it is the same molecule using the SMILES.
     * If a molecule occurs more than once in the Scaffold, only one corresponding ScaffoldNode is returned.
     * @param aMolecule molecule that is being searched for
     * @return ScaffoldNode of the searched molecule
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the ScaffoldCollection
     */
    public ScaffoldNodeBase getNode(IAtomContainer aMolecule) throws CDKException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'");
        if(!this.containsMolecule(aMolecule)) { //Check if the molecule exists in the ScaffoldCollection
            throw new IllegalArgumentException("Molecule is not in ScaffoldCollection");
        }
        return this.smilesMap.get(this.smilesGenerator.create(aMolecule)).get(0);
    }

    /**
     * Outputs the maximum level of the collection
     * @return the maximum level of the collection
     */
    public int getMaxLevel() {
        List<Integer> tmpLevelList = new ArrayList<>();
        tmpLevelList.addAll(this.levelMap.keys());
        return Collections.max(tmpLevelList);
    }

    /**
     * Return all ScaffoldNodes that are at a certain level in the ScaffoldCollection.
     * @param aLevel Level whose ScaffoldNodes are to be returned
     * @return ScaffoldNodes that are at a certain level
     * @throws IllegalArgumentException if the level does not exist
     */
    public List<ScaffoldNodeBase> getAllNodesOnLevel(int aLevel) throws IllegalArgumentException {
        Objects.requireNonNull(aLevel, "Given number is 'null'");
        if(this.getMaxLevel() >= aLevel) { //Level must be less than or equal to the maximum level
            return this.levelMap.get(aLevel);
        }
        throw new IllegalArgumentException("Level does not exist: " + aLevel);
    }

    /**
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1.
     * All others are marked with 0. Each row and column number in the matrix is assigned to a node.
     * The assignment can be requested with getMatrixNodes/getMatrixNode.
     * @return the adjacency matrix
     * @throws IllegalStateException if the tree is not connected
     */
    public abstract Integer[][] getMatrix() throws IllegalStateException;

    /**
     * Returns the number of the nodes and the nodes of the matrix in ascending order.
     * @return HashMap with the number of the nodes and the nodes in the order they appear in the matrix
     */
    public HashMap<Integer, ScaffoldNodeBase> getMatrixNodes() {
        return this.nodeMap;
    }

    /**
     * Returns the node that belongs to a certain row and column number of the matrix.
     * Returns zero if the node is not in the ScaffoldCollection. Throwing an exception makes it difficult to build a ScaffoldCollection
     * @param aNumber Row and column number whose node is requested
     * @return Node that belongs to the requested column and row number
     */
    public ScaffoldNodeBase getMatrixNode(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.get(aNumber);
    }

    /**
     * Get all ScaffoldNodes of the Scaffold
     * @return all ScaffoldNodes of the Scaffold
     */
    public List<ScaffoldNodeBase> getAllNodes() {
        List<ScaffoldNodeBase> tmpList = new ArrayList<>();
        tmpList.addAll(this.levelMap.values());
        return tmpList;
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
    //</editor-fold>
}
