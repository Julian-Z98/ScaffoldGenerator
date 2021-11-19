/*
 * Copyright (c) 2021 Julian Zander, Jonas Schaub, Achim Zielesny, Christoph Steinbeck
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

import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.SingleGraph;
import org.graphstream.stream.file.FileSinkImages;
import org.graphstream.stream.file.images.Resolutions;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.concurrent.TimeUnit;

/**
 * Contains useful functions such as visualisation with Graphstream.
 * @author Julian Zander, Jonas Schaub (zanderjulian@gmx.de, jonas.schaub@uni-jena.de)
 * @version 1.0.1.0
 */
public final class GraphStreamUtility {
    /**
     * The ScaffoldNodeCollectionBase is displayed in a window with GraphStream.
     * The optional numbering of the nodes reflects the number of the node and the respective level.
     * @param aScaffoldCollection displayed Collection
     * @param aShowLabel Adds a label with node level and node number if true.
     * @throws Exception if anything goes wrong
     */
    public static void displayWithGraphStream(ScaffoldNodeCollectionBase aScaffoldCollection, boolean aShowLabel) throws Exception {
        /*Create a graph from the ScaffoldCollection*/
        Graph tmpGraph = new SingleGraph("TestGraph");
        System.setProperty("org.graphstream.ui", "swing");
        FileSinkImages tmpFileSinkImages = FileSinkImages.createDefault();
        tmpFileSinkImages.setOutputPolicy(FileSinkImages.OutputPolicy.BY_GRAPH_EVENT);
        tmpFileSinkImages.setOutputType(FileSinkImages.OutputType.jpg);
        tmpFileSinkImages.setQuality(FileSinkImages.Quality.HIGH);
        tmpFileSinkImages.setViewCenter(10000, 10000);
        //tmpFileSinkImages.setResolution(2048, 2048);
        tmpFileSinkImages.setResolution(Resolutions.UHD_4K);
        tmpFileSinkImages.setAutofit(false);
        tmpFileSinkImages.setLayoutPolicy(FileSinkImages.LayoutPolicy.COMPUTED_FULLY_AT_NEW_IMAGE);
        tmpGraph.setAttribute("ui.stylesheet", "node { shape: rounded-box; size-mode: fit; padding: 60px; } graph { shape: box; size-mode: fit; padding: 70px; }");
        tmpGraph.setAttribute("ui.quality");
        tmpGraph.setAttribute("ui.antialias");
        /*Add edges and nodes*/
        int tmpEdgeCount = 0;
        DepictionGenerator tmpGenerator = new DepictionGenerator().withSize(2048,2048).withFillToFit();
        Integer[][] tmpMatrix = aScaffoldCollection.getMatrix(); //Create the adjacency matrix
        for(int tmpRow = 0; tmpRow < tmpMatrix.length; tmpRow++) { //Create a node for each row
            /*Add the ScaffoldCollection nodes to the graph*/
            tmpGraph.addNode(String.valueOf(tmpRow));
            Node tmpNode = tmpGraph.getNode(String.valueOf(tmpRow));
            tmpNode.setAttribute("Node", aScaffoldCollection.getMatrixNode(tmpRow));
            /*Add a label to each node that corresponds to the level in the collection. 0 is the root.*/
            ScaffoldNodeBase tmpCollectionLevelNode =  aScaffoldCollection.getMatrixNode(tmpRow);
            /*Add a label if true*/
            if(aShowLabel) {
                String tmpLabel = tmpCollectionLevelNode.getLevel() + " " + tmpRow;
                //tmpNode.setAttribute("ui.label", tmpLabel);
            }
            /*Add the images*/
            ScaffoldNodeBase tmpCollectionNode = aScaffoldCollection.getMatrixNode(aScaffoldCollection.getMatrixNodesNumbers().get(tmpRow));
            IAtomContainer tmpCollectionNodeMolecule = (IAtomContainer) tmpCollectionNode.getMolecule();
            BufferedImage tmpNodeImg = tmpGenerator.depict(tmpCollectionNodeMolecule).toImg();
            /*The images are stored temporarily*/
            new File(System.getProperty("user.dir") + "//target/test-classes").mkdirs();
            File tmpSecOutputRemove = new File(System.getProperty("user.dir") + "//target/test-classes/GraphStream" + tmpRow + ".png");
            ImageIO.write(tmpNodeImg, "png", tmpSecOutputRemove);
            //set the images
            tmpNode.setAttribute("ui.style", "fill-mode: image-scaled;" + "fill-image: url('GraphStream" + tmpRow + ".png');");

            //tmpNode.setAttribute("ui.style", "fill-mode: image-scaled-ratio-max;" + "fill-image: url('GraphStream" + tmpRow + ".png');");
            /*Add edges*/
            for(int tmpCol = 0; tmpCol < tmpMatrix[tmpRow].length; tmpCol++) { //Go through each column of the row
                if(tmpRow < tmpCol) { //Skip a diagonal half to get edges in one direction only.
                    continue;
                }
                if(tmpMatrix[tmpRow][tmpCol].equals(1)) { //Insert an edge if there is a 1 in it
                    tmpGraph.addEdge("Edge" + tmpEdgeCount, tmpRow, tmpCol);
                    tmpEdgeCount++;
                }
            }
        }
        /*Display graph*/
        System.setProperty("org.graphstream.ui", "swing");
        tmpGraph.display();
        //BufferedImage bufferedImage = new BufferedImage(2048, 2048, BufferedImage.TYPE_INT_RGB);
        //Graphics2D tmpTestImage = bufferedImage.createGraphics();
        File tmpImageFile = new File("outputFile.jpg");
        //ImageIO.write(bufferedImage, "png", tmpImageFile);
        tmpFileSinkImages.writeAll(tmpGraph, tmpImageFile.getAbsolutePath());
        //tmpFileSinkImages.writeAll(tmpGraph, tmpImageFile.getAbsolutePath());
        //tmpGraph.setAttribute("ui.screenshot", "screenshot.jpg"); // Saved at: C:\Users\zande\IdeaProjects\ScaffoldGenerator\ScaffoldGenerator
        TimeUnit.SECONDS.sleep(300);
    }
}
