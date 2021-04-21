import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;
import org.openscience.cdk.renderer.AtomContainerRenderer;
import org.openscience.cdk.renderer.RendererModel;
import org.openscience.cdk.renderer.SymbolVisibility;
import org.openscience.cdk.renderer.color.UniColor;
import org.openscience.cdk.renderer.elements.IRenderingVisitor;
import org.openscience.cdk.renderer.font.AWTFontManager;
import org.openscience.cdk.renderer.generators.BasicSceneGenerator;
import org.openscience.cdk.renderer.generators.standard.StandardGenerator;
import org.openscience.cdk.renderer.visitor.IDrawVisitor;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.Arrays;

public class murckoTest {
    public static void main(String[] args) throws IOException, CDKException {
        //Load molecule
        InputStream tmpInputStream = new FileInputStream(
                new File("D:/Studium/Scaffold/MOL_Files/Test3.mol"));
        MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
        IAtomContainer tmpTestMol =(IAtomContainer) tmpReader.read(new AtomContainer());
        //Generate picture
        DepictionGenerator dptgen = new DepictionGenerator();
        //Genereate window
        dptgen.withSize(300, 350).withMolTitle().withTitleColor(Color.DARK_GRAY);
        BufferedImage tmpImg = dptgen.depict(tmpTestMol).toImg();
        ImageIcon tmpIcon=new ImageIcon(tmpImg);
        JFrame frame=new JFrame();
        frame.setLayout(new FlowLayout());
        frame.setSize(320,370);
        JLabel lbl=new JLabel();
        lbl.setIcon(tmpIcon);
        frame.add(lbl);
        frame.setTitle("Original");
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter();
    }
}
