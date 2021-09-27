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
import java.util.List;
import java.util.Objects;

/**
 * The NetworkNodes are nodes from which a Network can be built.
 * It is used to organise the IAtomContainers and enables a relationship between the different objects.
 * A NetworkNode can have multiple children and parents
 * @param <MoleculeType> As MoleculeType, any data type can be defined.
 *                     In our scenario, the nodes contain molecules.
 */
public class NetworkNode <MoleculeType> extends ScaffoldNodeBase<MoleculeType> {

    /**
     * parents of the node
     */
    public List<NetworkNode<MoleculeType>> parents;

    /**
     * Shows if the node has parents
     * @return Whether the node has parents
     */
    public boolean isOrphan() {
        return parents.isEmpty();
    }

    /**
     * Creates a NetworkNode
     * @param aMolecule molecule of the NetworkNode
     */
    public NetworkNode(MoleculeType aMolecule) {
        super(aMolecule);
        this.parents =  new ArrayList<NetworkNode<MoleculeType>>();
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
     * Add the parents node and add this node as child to the parent node if not already done.
     * @param aParent parent that are added
     */
    public void addParent(NetworkNode<MoleculeType> aParent) {
        Objects.requireNonNull(aParent, "Given NetworkNode is 'null'");
        /*Add child if not already added*/
        boolean tmpIsAlreadyChild = false;
        for(ScaffoldNodeBase<MoleculeType> tmpBaseNode : aParent.getChildren()) {
            NetworkNode<MoleculeType> tmpNode = (NetworkNode<MoleculeType>) tmpBaseNode;
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

    //<editor-fold desc="get/set">
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
    //</editor-fold>
}
