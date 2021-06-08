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
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;

/**
 * Source: https://github.com/gt4dev/yet-another-tree-structure
 */
public class TreeNode<T> implements Iterable<TreeNode<T>> {

    private T data;
    private TreeNode<T> parent;
    private List<TreeNode<T>> children;

    public boolean isRoot() {
        return parent == null;
    }

    public boolean isLeaf() {
        return children.size() == 0;
    }

    private List<TreeNode<T>> elementsIndex;

    public TreeNode(T aData) {
        Objects.requireNonNull(aData, "Given data is 'null'");
        this.data = aData;
        this.children = new LinkedList<TreeNode<T>>();
        this.elementsIndex = new LinkedList<TreeNode<T>>();
        this.elementsIndex.add(this);
    }

    public TreeNode<T> addChild(T aChild) {
        Objects.requireNonNull(aChild, "Given child is 'null'");
        TreeNode<T> childNode = new TreeNode<T>(aChild);
        childNode.parent = this;
        this.children.add(childNode);
        this.registerChildForSearch(childNode);
        return childNode;
    }

    public int getLevel() {
        if (this.isRoot())
            return 0;
        else
            return parent.getLevel() + 1;
    }

    private void registerChildForSearch(TreeNode<T> aNode) {
        Objects.requireNonNull(aNode, "Given Tree Node is 'null'");
        elementsIndex.add(aNode);
        if (parent != null)
            parent.registerChildForSearch(aNode);
    }

    public TreeNode<T> findTreeNode(Comparable<T> aCompData) {
        Objects.requireNonNull(aCompData, "Given comparable data is 'null'");
        for (TreeNode<T> element : this.elementsIndex) {
            T elData = element.data;
            if (aCompData.compareTo(elData) == 0)
                return element;
        }

        return null;
    }

    @Override
    public String toString() {
        return data != null ? data.toString() : "[data null]";
    }

    @Override
    public Iterator<TreeNode<T>> iterator() {
        TreeNodeIter<T> iter = new TreeNodeIter<T>(this);
        return iter;
    }

    //<editor-fold desc="get/set">
    /**
     * Get the node data.
     * @return node data
     */
    public T getData() {
        return this.data;
    }

    /**
     * Set the node data.
     * @param aData data that are set
     */
    public void setData(T aData) {
        Objects.requireNonNull(aData, "Given Data is 'null'");
        this.data = aData;
    }

    /**
     * Get the parent node.
     * @return parent node
     */
    public TreeNode<T> getParent() {
        return this.parent;
    }

    /**
     * Set the parent node.
     * @param aParent parent that are set
     */
    public void setParent(TreeNode<T> aParent) {
        Objects.requireNonNull(aParent, "Given TreeNode is 'null'");
        this.parent = aParent;
    }

    /**
     * Get the children nodes.
     * @return children nodes
     */
    public List<TreeNode<T>> getChildren() {
        return this.children;
    }

    /**
     * Set the children node.
     * @param aChildren children that are set
     */
    public void setChildren(List<TreeNode<T>> aChildren) {
        Objects.requireNonNull(aChildren, "Given TreeNode List is 'null'");
        this.children = aChildren;
    }
    //</editor-fold>




}