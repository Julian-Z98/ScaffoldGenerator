/*
 * Copyright (c) 2021 Julian Zander, Jonas Schaub, Achim Zielesny
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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.*;

/**
 * Top-level class to organise the NetworkNodes
 */
public class ScaffoldNetwork {
    /**
     * Saves all NetworkNodes and numbers them in ascending order. Starts at 0. Key:Number, Value:NetworkNode
     */
    private HashMap<Integer, NetworkNode> nodeMap;
    /**
     * Saves all NetworkNodes and numbers them in ascending order. Starts at 0.
     * nodeMap with key and value swapped. Key:NetworkNode, Value:Number
     */
    private HashMap<NetworkNode, Integer> reverseNodeMap;
    /**
     * Saves all NetworkNodes according to their Unique SMILES. Key:SMILES, Value:NetworkNode
     */
    private ListMultimap<String, NetworkNode> smilesMap;
    /**
     * Saves all NetworkNodes according to their level. Key:Level, Value:NetworkNode
     */
    private ListMultimap<Integer, NetworkNode> levelMap;
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
    public ScaffoldNetwork(SmilesGenerator aSmilesGenerator) {
        this.nodeMap = new HashMap<Integer, NetworkNode>();
        this.reverseNodeMap = new HashMap<NetworkNode, Integer>();
        this.smilesMap = ArrayListMultimap.create();
        this.levelMap = ArrayListMultimap.create();
        this.smilesGenerator = aSmilesGenerator;
        this.nodeCounter = 0;
    }

    /**
     * Default Constructor
     */
    public ScaffoldNetwork() {
        this(new SmilesGenerator(SmiFlavor.Unique));
    }

    /**
     * Add a new NetworkNode to the ScaffoldNetwork
     * @param aNode NetworkNode to be added
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public NetworkNode addNode(NetworkNode aNode) throws CDKException {
        Objects.requireNonNull(aNode, "Given NetworkNode is 'null'");
        Objects.requireNonNull(aNode.getMolecule(), "Given Data is 'null'");
        //Add to nodeMap
        this.nodeMap.put(this.nodeCounter, aNode);
        //Add to reverseNodeMap
        this.reverseNodeMap.put(aNode, this.nodeCounter);
        /*Add to smilesMap*/
        IAtomContainer tmpMolecule = (IAtomContainer) aNode.getMolecule();
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
        this.smilesMap.put(tmpSmiles, aNode);
        /*Since the network is built from the leaves to the root,
         the levels of all nodes in the network must be re-determined for every node added.*/
        this.levelMap.put(aNode.getLevel(), aNode);
        ListMultimap<Integer, NetworkNode> tmpLevelMap = ArrayListMultimap.create();
        for(NetworkNode tmpNode : this.getAllNodes()) {
            tmpLevelMap.put(tmpNode.getLevel(), tmpNode);
        }
        this.levelMap = tmpLevelMap;
        //Increase nodeCounter
        this.nodeCounter++;
        return this.getNetworkNode((IAtomContainer) aNode.getMolecule());
    }

    /**
     * Removes a Node. This does not change the order. The numbering does not move up.
     * @param aNode Node to remove
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the network
     */
    public void removeNode(NetworkNode aNode) throws CDKException, IllegalArgumentException {
        Objects.requireNonNull(aNode, "Given NetworkNode is 'null'");
        if(!this.reverseNodeMap.containsKey(aNode)) { //Check if the node exists in the network
            throw new IllegalArgumentException("Node is not in network");
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
     * Checks whether the molecule is already present in the Scaffold Network
     * Check whether it is the same molecule using the Unique SMILES.
     * Unique SMILES: canonical SMILES string, different atom ordering produces the same (apart from rare exceptions) SMILES.
     *                No isotope or stereochemistry encoded.
     * @param aMolecule Molecule to check
     * @return Whether the molecule is located in the Scaffold Network
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public boolean isMoleculeInNetwork(IAtomContainer aMolecule) throws CDKException {
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
     * Return the NetworkNode that belongs to a specific molecule.
     * Check whether it is the same molecule using the Unique SMILES.
     * If a molecule occurs more than once in the network, only one corresponding NetworkNode is returned.
     * Unique SMILES: canonical SMILES string, different atom ordering produces the same (apart from rare exceptions) SMILES.
     *                No isotope or stereochemistry encoded.
     * @param aMolecule molecule that is being searched for
     * @return NetworkNode of the searched molecule
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the network
     */
    public NetworkNode getNetworkNode(IAtomContainer aMolecule) throws CDKException, IllegalArgumentException {
        Objects.requireNonNull(aMolecule, "Given atom container is 'null'");
        if(!this.isMoleculeInNetwork(aMolecule)) { //Check if the molecule exists in the network
            throw new IllegalArgumentException("Molecule is not in network");
        }
        return this.smilesMap.get(this.smilesGenerator.create(aMolecule)).get(0);
    }

    /**
     * Outputs the maximum level of the network
     * @return the maximum level of the network
     */
    public int getMaxLevel() {
        List<Integer> tmpLevelList = new ArrayList<>();
        tmpLevelList.addAll(this.levelMap.keys());
        return Collections.max(tmpLevelList);
    }

    /**
     * Return all NetworkNodes that are at a certain level in the network.
     * @param aLevel Level whose NetworkNodes are to be returned
     * @return NetworkNodes that are at a certain level
     * @throws IllegalArgumentException if the level does not exist
     */
    public List<NetworkNode> getAllNodesOnLevel(int aLevel) throws IllegalArgumentException {
        Objects.requireNonNull(aLevel, "Given number is 'null'");
        if(this.getMaxLevel() >= aLevel) { //Level must be less than or equal to the maximum level
            return this.levelMap.get(aLevel);
        }
        throw new IllegalArgumentException("Level does not exist: " + aLevel);
    }

    /**
     * Get all NetworkNodes of the ScaffoldNetwork
     * @return all NetworkNodes of the ScaffoldNetwork
     */
    public List<NetworkNode> getAllNodes() {
        List<NetworkNode> tmpList = new ArrayList<>();
        tmpList.addAll(this.levelMap.values());
        return tmpList;
    }

    /**
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1.
     * All others are marked with 0. Each row and column number in the matrix is assigned to a node.
     * The assignment can be requested with getMatrixNodes/getMatrixNode.
     * only works with connected networks. Can be checked with isNetworkConnected.
     * @return the adjacency matrix
     * @throws IllegalStateException if the network is not connected
     */
    public Integer[][] getNetworkAsMatrix() throws IllegalStateException {
        int tmpSize = this.nodeMap.size();
        Integer[][] tmpMatrix = new Integer[tmpSize][tmpSize];
        /*Set all values of the matrix to 0*/
        for (int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            for (int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) {
                tmpMatrix[tmpRow][tmpCol] = 0;
            }
        }
        /*Insert a 1 for each parent node*/
        int tmpCounter = 0;
        for (NetworkNode tmpNode : this.nodeMap.values()) {
            if (tmpNode.getParents() != null) {
                for(Object tmpParentNode : tmpNode.getParents()) {
                    //Set a 1 at the level of the parent and at the level of the node
                    tmpMatrix[tmpCounter][this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpParentNode))] = 1;
                    //Set a 1 at the level of the node and at the level of the parent
                    tmpMatrix[this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpParentNode))][tmpCounter] = 1;
                }
            }
            tmpCounter++;
        }
        return tmpMatrix;
    }

    /**
     * Returns the number of the nodes and the nodes of the matrix in ascending order.
     * @return HashMap with the number of the nodes and the nodes in the order they appear in the matrix
     */
    public HashMap<Integer, NetworkNode> getMatrixNodes() {
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
     * Indicates whether there is a node with this number in the network.
     * @param aNumber Number of the node to be checked
     * @return Is the node in the network
     */
    public boolean isNodeInNetwork(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.containsKey(aNumber);
    }

    /**
     * Returns the node that belongs to a certain row and column number of the matrix.
     * Returns zero if the node is not in the network. Throwing an exception makes it difficult to build a network
     * @param aNumber Row and column number whose node is requested
     * @return Node that belongs to the requested column and row number
     */
    public NetworkNode getMatrixNode(int aNumber) {
        Objects.requireNonNull(aNumber, "Given int is 'null'");
        return this.nodeMap.get(aNumber);
    }

    /**
     * Outputs root nodes of the network.
     * @return root nodes of the network
     */
    public List<NetworkNode> getRoots() {
        List<NetworkNode> tmpNodeList = new ArrayList<>();
        for(NetworkNode tmpNode : this.nodeMap.values()) {
            /*Is the parent of the node in the network*/
            if(tmpNode.getLevel() == 0) {
                //If the node has no parent, it is a root
                tmpNodeList.add(tmpNode);
            }
        }
        return tmpNodeList;
    }

    /**
     * Adds another ScaffoldNetwork to the existing one if possible.
     * The new network is inserted at the node that both networks have in common.
     * All children of the new network at this node and there linkages are taken over if they do not already exist.
     *
     * The new network is simply taken over if the existing network is empty.
     * If there is no match between the two networks the new network is inserted without linkages.
     *
     * If a molecule does not generate a SchuffenhauerScaffold, it is stored as a node with empty SMILES and is treated normally.
     * All other empty nodes are then added to this network accordingly.
     * @param aScaffoldNetwork network to be inserted into the existing ScaffoldNetwork.
     * @return false if the new network has nodes and cannot be inserted because the two networks have no common nodes.
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public void mergeNetwork(ScaffoldNetwork aScaffoldNetwork) throws CDKException {
        /*If the old ScaffoldNetwork is empty, transfer the new ScaffoldNetwork to be added.*/
        if(this.getRoots().isEmpty()) {
            for(Object tmpNode : aScaffoldNetwork.getAllNodes()) {
                this.addNode((NetworkNode) tmpNode);
            }
        }
        /*If the old Scaffold network is not empty*/
        else {
            int tmpMaxNodeNumber = aScaffoldNetwork.getAllNodes().size();
            ArrayList<NetworkNode> tmpNodesToAdd = new ArrayList(tmpMaxNodeNumber);
            ArrayList<IAtomContainer> tmpMoleculeList = new ArrayList(tmpMaxNodeNumber);
            /*Go through each level of the network starting at the root*/
            for(int i = 0; i <= aScaffoldNetwork.getMaxLevel(); i++) {
                /*Compare all nodes of the old network on this level with all nodes of the new network on this level*/
                for(NetworkNode tmpOldNetworkNode : this.getAllNodesOnLevel(i)) {
                    for(Object tmpNewNetworkObject : aScaffoldNetwork.getAllNodesOnLevel(i)) {
                        NetworkNode tmpNewNetworkNode = (NetworkNode) tmpNewNetworkObject;
                        /*Generate the corresponding SMILES of the molecules*/
                        String tmpOldSmiles = this.smilesGenerator.create((IAtomContainer) tmpOldNetworkNode.getMolecule());
                        String tmpNewSmiles = this.smilesGenerator.create((IAtomContainer) tmpNewNetworkNode.getMolecule());
                        /*Check whether a fragment occurs in both networks*/
                        if(tmpOldSmiles.equals(tmpNewSmiles)) {
                            /*Add the origin smiles to the OldSmilesNetwork fragment*/
                            for(Object tmpOriginSmiles : tmpNewNetworkNode.getOriginSmilesList()) {
                                tmpOldNetworkNode.addOriginSmiles((String) tmpOriginSmiles);
                            }
                            /*Networks are overlapping if a fragment occurs in both networks
                            Get the children of the overlapping node*/
                            for(Object tmpNewChild : tmpNewNetworkNode.getChildren()) {
                                NetworkNode tmpNewChildNode = (NetworkNode) tmpNewChild;
                                IAtomContainer tmpNewChildMolecule = (IAtomContainer) tmpNewChildNode.getMolecule();
                                /*Add the child if it is not already in the network*/
                                if(!this.isMoleculeInNetwork(tmpNewChildMolecule)) {
                                    NetworkNode tmpNewNode = new NetworkNode<>(tmpNewChildMolecule);
                                    tmpNewNode.addParent(tmpOldNetworkNode);
                                    tmpMoleculeList.add(tmpNewChildMolecule);
                                    this.addNode(tmpNewNode);
                                }
                            }
                        } else { /*Node is not in network*/
                            /*Add node to lists so that it is added to the network later*/
                            tmpNodesToAdd.add(tmpNewNetworkNode);
                            tmpMoleculeList.add((IAtomContainer) tmpNewNetworkNode.getMolecule());
                        }

                    }
                }
            }
            /*Add the nodes of the new network that are not yet in the whole network*/
            for(NetworkNode tmpNode : tmpNodesToAdd) {
                if(!this.isMoleculeInNetwork((IAtomContainer) tmpNode.getMolecule())) {
                    NetworkNode tmpNewNode = new NetworkNode<>((IAtomContainer) tmpNode.getMolecule());
                    this.addNode(tmpNewNode);
                }
            }
            /*Add the matching parents to the newly added nodes. Childs are automatically set when the parents are set.*/
            for(IAtomContainer tmpMolecule : tmpMoleculeList) {
                ArrayList<NetworkNode> tmpParentList = (ArrayList<NetworkNode>) aScaffoldNetwork.getNetworkNode(tmpMolecule).getParents();
                for(NetworkNode tmpParentNode : tmpParentList) {
                    /*Only molecules that are in the network*/
                    if(this.isMoleculeInNetwork((IAtomContainer) tmpParentNode.getMolecule())) {
                        NetworkNode tmpOldParentNode = this.getNetworkNode((IAtomContainer) tmpParentNode.getMolecule());
                        this.getNetworkNode(tmpMolecule).addParent(tmpOldParentNode);
                    }
                }
            }
        }
    }
}
