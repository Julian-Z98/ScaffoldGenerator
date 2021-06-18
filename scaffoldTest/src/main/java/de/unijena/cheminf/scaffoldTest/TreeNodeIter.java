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
import java.util.Objects;

/**
 * Iterator class for the TreeNode class. Allows the tree to be walked through completely.
 *
 * Inspired by: https://github.com/gt4dev/yet-another-tree-structure
 */
public class TreeNodeIter<IAtomContainer> implements Iterator<TreeNode<IAtomContainer>> {

    /**
     * Indicates which part of the tree is currently being processed.
     */
    enum PROCESS_STAGES {
        ProcessParent, ProcessChildCurNode, ProcessChildSubNode
    }

    /**
     * Just handled node of the tree
     */
    private TreeNode<IAtomContainer> treeNode;

    /**
     * Next PROCESS_STAGE
     */
    private PROCESS_STAGES doNext;

    /**
     * Next IAtomContainer being processed
     */
    private TreeNode<IAtomContainer> next;

    /**
     * Shows whether the nextNode() has already been calculated.
     * nextNode() must have been executed before next() can be executed.
     */
    private boolean nextNodeHasBeenCalculated;

    /**
     * Currently treated TreeNode
     */
    private Iterator<TreeNode<IAtomContainer>> childrenCurNodeIter;
    /**
     * Next treated TreeNode
     */
    private Iterator<TreeNode<IAtomContainer>> childrenSubNodeIter;

    /**
     * Constructor
     * Iteration through the tree of a specific node
     * @param aTreeNode Node through whose tree it is to be iterated
     */
    public TreeNodeIter(TreeNode<IAtomContainer> aTreeNode) {
        Objects.requireNonNull(aTreeNode, "Given TreeNode is 'null'");
        this.treeNode = aTreeNode;
        this.doNext = PROCESS_STAGES.ProcessParent;
        this.childrenCurNodeIter = treeNode.getChildren().iterator();
    }

    /**
     * Indicates whether there is a next TreeNode. Returns false if the tree has already been traversed completely.
     * @return Is there a next TreeNode in the tree
     */
    @Override
    public boolean hasNext() {
        this.nextNodeHasBeenCalculated = true;
        if (this.doNext == PROCESS_STAGES.ProcessParent) {
            this.next = this.treeNode;
            this.doNext = PROCESS_STAGES.ProcessChildCurNode;
            return true;
        }

        if (this.doNext == PROCESS_STAGES.ProcessChildCurNode) {
            if (childrenCurNodeIter.hasNext()) {
                TreeNode<IAtomContainer> childDirect = childrenCurNodeIter.next();
                childrenSubNodeIter = childDirect.iterator();
                this.doNext = PROCESS_STAGES.ProcessChildSubNode;
                return hasNext();
            }

            else {
                this.doNext = null;
                return false;
            }
        }

        if (this.doNext == PROCESS_STAGES.ProcessChildSubNode) {
            if (childrenSubNodeIter.hasNext()) {
                this.next = childrenSubNodeIter.next();
                return true;
            }
            else {
                this.next = null;
                this.doNext = PROCESS_STAGES.ProcessChildCurNode;
                return hasNext();
            }
        }
        return false;
    }

    /**
     * Outputs the next TreeNode in the tree
     * @return next TreeNode in the tree
     */
    @Override
    public TreeNode<IAtomContainer> next() {
        if(this.nextNodeHasBeenCalculated == false) { //execute hasnext() if it has not been executed yet
            this.hasNext();
        }
        this.nextNodeHasBeenCalculated = false;
        return this.next;
    }
}