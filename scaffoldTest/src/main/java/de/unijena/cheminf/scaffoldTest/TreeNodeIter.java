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
 * Source: https://github.com/gt4dev/yet-another-tree-structure
 */
public class TreeNodeIter<T> implements Iterator<TreeNode<T>> {

    enum ProcessStages {
        ProcessParent, ProcessChildCurNode, ProcessChildSubNode
    }

    private TreeNode<T> treeNode;

    public TreeNodeIter(TreeNode<T> aTreeNode) {
        Objects.requireNonNull(aTreeNode, "Given TreeNode is 'null'");
        this.treeNode = aTreeNode;
        this.doNext = ProcessStages.ProcessParent;
        this.childrenCurNodeIter = treeNode.getChildren().iterator();
    }

    private ProcessStages doNext;
    private TreeNode<T> next;
    private boolean nextNodeHasBeenCalculated; //nextNode() must have been executed before next() can be executed.
    private Iterator<TreeNode<T>> childrenCurNodeIter;
    private Iterator<TreeNode<T>> childrenSubNodeIter;

    @Override
    public boolean hasNext() {
        this.nextNodeHasBeenCalculated = true;
        if (this.doNext == ProcessStages.ProcessParent) {
            this.next = this.treeNode;
            this.doNext = ProcessStages.ProcessChildCurNode;
            return true;
        }

        if (this.doNext == ProcessStages.ProcessChildCurNode) {
            if (childrenCurNodeIter.hasNext()) {
                TreeNode<T> childDirect = childrenCurNodeIter.next();
                childrenSubNodeIter = childDirect.iterator();
                this.doNext = ProcessStages.ProcessChildSubNode;
                return hasNext();
            }

            else {
                this.doNext = null;
                return false;
            }
        }

        if (this.doNext == ProcessStages.ProcessChildSubNode) {
            if (childrenSubNodeIter.hasNext()) {
                this.next = childrenSubNodeIter.next();
                return true;
            }
            else {
                this.next = null;
                this.doNext = ProcessStages.ProcessChildCurNode;
                return hasNext();
            }
        }
        return false;
    }

    @Override
    public TreeNode<T> next() {
        if(this.nextNodeHasBeenCalculated == false) { //execute hasnext() if it has not been executed yet
            this.hasNext();
        }
        this.nextNodeHasBeenCalculated = false;
        return this.next;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}