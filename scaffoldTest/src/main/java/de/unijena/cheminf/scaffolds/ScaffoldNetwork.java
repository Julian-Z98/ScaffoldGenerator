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

package de.unijena.cheminf.scaffolds;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

/**
 * Top-level class to organise the NetworkNodes
 *
 * @version 1.0
 */
public class ScaffoldNetwork extends ScaffoldNodeCollectionBase {

    //<editor-fold desc="Constructors">
    /**
     * Constructor
     * @param aSmilesGenerator Used SMILES Generator
     */
    public ScaffoldNetwork(SmilesGenerator aSmilesGenerator) {
        super(aSmilesGenerator);
    }

    /**
     * Default Constructor
     */
    public ScaffoldNetwork() {
        this(new SmilesGenerator(SmiFlavor.Unique));
    }
    //</editor-fold>

    //<editor-fold desc="Public methods">
    /**
     * Add a new NetworkNode to the ScaffoldNetwork
     * @param aNode NetworkNode to be added
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    @Override
    public void addNode(ScaffoldNodeBase aNode) throws CDKException {
        /*Parameter checks*/
        Objects.requireNonNull(aNode, "Given NetworkNode is 'null'");
        NetworkNode tmpNode = null;
        try {
            tmpNode = (NetworkNode) aNode;
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println("Node cant be added to NetworkTree");
            System.out.println("Parameter must be a NetworkNode");
        }
        //Add to nodeMap
        this.nodeMap.put(this.nodeCounter, tmpNode);
        //Add to reverseNodeMap
        this.reverseNodeMap.put(tmpNode, this.nodeCounter);
        /*Add to smilesMap*/
        IAtomContainer tmpMolecule = (IAtomContainer) tmpNode.getMolecule();
        String tmpSmiles = this.smilesGenerator.create(tmpMolecule); //Convert molecule to SMILES
        this.smilesMap.put(tmpSmiles, tmpNode);
        /*Since the network is built from the leaves to the root,
         the levels of all nodes in the network must be re-determined for every node added.*/
        this.levelMap.put(tmpNode.getLevel(), tmpNode);
        ListMultimap<Integer, ScaffoldNodeBase> tmpLevelMap = ArrayListMultimap.create();
        for(ScaffoldNodeBase tmpNodeBase : this.getAllNodes()) {
            NetworkNode tmpNetworkNode = (NetworkNode) tmpNodeBase;
            tmpLevelMap.put(tmpNode.getLevel(), tmpNetworkNode);
        }
        this.levelMap = tmpLevelMap;
        //Increase nodeCounter
        this.nodeCounter++;
    }

    /**
     * Removes a Node. This does not change the order. The numbering does not change.
     * @param aNode Node to remove
     * @throws CDKException In case of a problem with the SmilesGenerator
     * @throws IllegalArgumentException if the node is not in the Scaffold
     */
    @Override
    public void removeNode(ScaffoldNodeBase aNode) throws CDKException, IllegalArgumentException {
        /*Parameter checks*/
        Objects.requireNonNull(aNode, "Given ScaffoldNode is 'null'");
        NetworkNode tmpNode = null;
        try {
            tmpNode = (NetworkNode) aNode;
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println("Node cant be removed to NetworkTree");
            System.out.println("Parameter must be a NetworkNode");
        }
        if(!this.reverseNodeMap.containsKey(tmpNode)) { //Check if the node exists in the Scaffold
            throw new IllegalArgumentException("Node is not in Scaffold");
        }
        /*Remove from nodeMap and reverseNodeMap*/
        int tmpNumberInNodeMap = this.reverseNodeMap.get(tmpNode); //get number in nodeMap
        this.nodeMap.remove(tmpNumberInNodeMap);
        this.reverseNodeMap.remove(tmpNode);
        /*Remove from smilesMap*/
        String tmpSmiles = this.smilesGenerator.create((IAtomContainer) tmpNode.getMolecule()); //Convert molecule to SMILES
        this.smilesMap.remove(tmpSmiles, tmpNode);
        /*Remove from levelMap*/
        levelMap.remove(tmpNode.getLevel(), tmpNode);
    }

    /**
     * Adds another ScaffoldNetwork to the existing one if possible.
     * The new network is inserted at the node that both networks have in common.
     * All children of the new network at this node and there linkages are taken over if they do not already exist.
     *
     * The new network is simply taken over if the existing network is empty.
     * If there is no match between the two networks the new network is inserted without linkages.
     *
     * If a molecule generates an empty scaffold, it is stored as a node with empty SMILES and is treated normally.
     * All other empty nodes are then added to this network accordingly.
     * By querying the origins of this node, all molecules that do not produce a scaffold can be returned.
     * @param aScaffoldNetwork network to be inserted into the existing ScaffoldNetwork.
     * @throws CDKException In case of a problem with the SmilesGenerator
     */
    public void mergeNetwork(ScaffoldNetwork aScaffoldNetwork) throws CDKException {
        /*If the old ScaffoldNetwork is empty, transfer the new ScaffoldNetwork to be added.*/
        if(this.getAllNodes().size() == 0  || this.getAllNodes().isEmpty()) {
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
                for(ScaffoldNodeBase tmpOldNetworkNodeBase : this.getAllNodesOnLevel(i)) {
                    NetworkNode tmpOldNetworkNode = (NetworkNode) tmpOldNetworkNodeBase;
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
                                if(!this.containsMolecule(tmpNewChildMolecule)) {
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
                if(!this.containsMolecule((IAtomContainer) tmpNode.getMolecule())) {
                    NetworkNode tmpNewNode = new NetworkNode<>((IAtomContainer) tmpNode.getMolecule());
                    this.addNode(tmpNewNode);
                }
            }
            /*Add the matching parents to the newly added nodes. Children are automatically set when the parents are set.*/
            for(IAtomContainer tmpMolecule : tmpMoleculeList) {
                ScaffoldNodeBase tmpChildBase = aScaffoldNetwork.getNode(tmpMolecule);
                NetworkNode tmpChild = (NetworkNode) tmpChildBase;
                ArrayList<NetworkNode> tmpParentList = (ArrayList<NetworkNode>) tmpChild.getParents();
                for(NetworkNode tmpParentNode : tmpParentList) {
                    /*Only molecules that are in the network*/
                    if(this.containsMolecule((IAtomContainer) tmpParentNode.getMolecule())) {
                        NetworkNode tmpOldParentNode = (NetworkNode) this.getNode((IAtomContainer) tmpParentNode.getMolecule());
                        NetworkNode tmpOldChild = (NetworkNode) this.getNode(tmpMolecule);
                        tmpOldChild.addParent(tmpOldParentNode);
                    }
                }
            }
        }
    }

    //<editor-fold desc="get/set">
    /**
     * Outputs an adjacency matrix in which the parent node of each node is marked with a 1.
     * All others are marked with 0. Each row and column number in the matrix is assigned to a node.
     * The assignment can be requested with getMatrixNodes/getMatrixNode.
     * only works with connected networks. Can be checked with isConnected.
     * @return the adjacency matrix
     * @throws IllegalStateException if the network is not connected
     */
    @Override
    public Integer[][] getMatrix() throws IllegalStateException {
        int tmpSize = this.nodeMap.size();
        Integer[][] tmpMatrix = new Integer[tmpSize][tmpSize];
        /*Set all values of the matrix to 0*/
        for (int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) {
            Arrays.fill(tmpMatrix[tmpRow], 0);
        }
        /*Insert a 1 for each parent node*/
        int tmpCounter = 0;
        for (ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            NetworkNode tmpNode = (NetworkNode) tmpNodeBase;
            if (tmpNode.getParents() != null) {
                for(Object tmpParentNode : tmpNode.getParents()) {
                    /*Check if a node has been removed*/
                    if(this.reverseNodeMap.get(tmpParentNode) != null){
                        //Set a 1 at the level of the parent and at the level of the node
                        tmpMatrix[tmpCounter][this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpParentNode))] = 1;
                        //Set a 1 at the level of the node and at the level of the parent
                        tmpMatrix[this.getMatrixNodesNumbers().indexOf(this.reverseNodeMap.get(tmpParentNode))][tmpCounter] = 1;
                    }
                }
            }
            tmpCounter++;
        }
        return tmpMatrix;
    }

    /**
     * Outputs root nodes of the network.
     * @return root nodes of the network
     */
    public List<NetworkNode> getRoots() {
        List<NetworkNode> tmpNodeList = new ArrayList<>();
        for(ScaffoldNodeBase tmpNodeBase : this.nodeMap.values()) {
            NetworkNode tmpNode = (NetworkNode) tmpNodeBase;
            /*Is the parent of the node in the network*/
            if(tmpNode.getLevel() == 0) {
                //If the node has no parent, it is a root
                tmpNodeList.add(tmpNode);
            }
        }
        return tmpNodeList;
    }
    //</editor-fold>
    //</editor-fold>
}
