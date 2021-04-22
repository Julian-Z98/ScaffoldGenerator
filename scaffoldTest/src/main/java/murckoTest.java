import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV3000Reader;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

public class murckoTest {
    public static void main(String[] args) throws IOException, CDKException {
        //Load molecule
        InputStream inputStream = new FileInputStream("MOL_Files/TestO.mol"); // Try Test1-Test7
        MDLV2000Reader reader = new MDLV2000Reader(inputStream);
        IAtomContainer testMol = reader.read(new AtomContainer());
        //Generate picture of the original molecule
        DepictionGenerator generator = new DepictionGenerator();
        generator.withSize(300, 350).withMolTitle().withTitleColor(Color.BLACK);
        BufferedImage imgOri = generator.depict(testMol).toImg();
        ImageIcon iconOri = new ImageIcon(imgOri);
        JFrame frame= new JFrame();
        frame.setLayout(new FlowLayout());
        frame.setSize(1000,800);
        JLabel lblOri = new JLabel();
        lblOri.setIcon(iconOri);
        lblOri.setText("Original");
        frame.add(lblOri);
        //Generate fragments, rings and frameworks
        MurckoFragmenter murckoFragmenter = new MurckoFragmenter(false,1);
        murckoFragmenter.setComputeRingFragments(true);
        murckoFragmenter.generateFragments(testMol);
        IAtomContainer[] fragments = murckoFragmenter.getFragmentsAsContainers();
        IAtomContainer[] rings = murckoFragmenter.getRingSystemsAsContainers();
        IAtomContainer[] frameworks = murckoFragmenter.getFrameworksAsContainers();
        //Generate pictures of the fragments
        int tmpCountFra = 0;
        for(IAtomContainer tmpFragment : fragments) {
            tmpCountFra++;
            BufferedImage tmpImgFra = generator.withBackgroundColor(Color.LIGHT_GRAY).depict(tmpFragment).toImg();
            ImageIcon tmpIconFra = new ImageIcon(tmpImgFra);
            JFrame tmpFrameFra = new JFrame();
            tmpFrameFra.setLayout(new FlowLayout());
            tmpFrameFra.setSize(320,370);
            JLabel tmpLblFra = new JLabel();
            tmpLblFra.setIcon(tmpIconFra);
            tmpLblFra.setText("Fragment "+tmpCountFra);
            frame.add(tmpLblFra);
        }
        //Generate pictures of the rings
        int tmpCountRgs = 0;
        for(IAtomContainer tmpRing : rings) {
            tmpCountRgs++;
            BufferedImage tmpImgRgs = generator.withBackgroundColor(Color.GRAY).depict(tmpRing).toImg();
            ImageIcon tmpIconRgs = new ImageIcon(tmpImgRgs);
            JFrame tmpFrameRgs = new JFrame();
            tmpFrameRgs.setLayout(new FlowLayout());
            tmpFrameRgs.setSize(320,370);
            JLabel tmpLblRgs = new JLabel();
            tmpLblRgs.setIcon(tmpIconRgs);
            tmpLblRgs.setText("Ring "+tmpCountRgs);
            frame.add(tmpLblRgs);
        }
        //Generate pictures of the frameworks
        int tmpCountFrw = 0;
        for(IAtomContainer tmpFramework : frameworks) {
            tmpCountFrw++;
            BufferedImage tmpImgFrw = generator.withBackgroundColor(Color.PINK).depict(tmpFramework).toImg();
            ImageIcon tmpIconFrw = new ImageIcon(tmpImgFrw);
            JFrame tmpFrameFrw = new JFrame();
            tmpFrameFrw.setLayout(new FlowLayout());
            tmpFrameFrw.setSize(320,370);
            JLabel tmpLblFrw = new JLabel();
            tmpLblFrw.setIcon(tmpIconFrw);
            tmpLblFrw.setText("Framework "+tmpCountFrw);
            frame.add(tmpLblFrw);
        }
        //Show window
        frame.setVisible(true);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
}
