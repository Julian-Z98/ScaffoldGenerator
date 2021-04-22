import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fragment.MurckoFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV3000Reader;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;

public class murckoTest {
    public static void main(String[] args) throws IOException, CDKException {
        //Load molecule
        InputStream tmpInputStream = new FileInputStream("MOL_Files/Test3.mol");
        MDLV3000Reader tmpReader = new MDLV3000Reader(tmpInputStream);
        IAtomContainer tmpTestMol = tmpReader.read(new AtomContainer());
        //Generate picture of the original molecule
        DepictionGenerator tmpGenerator = new DepictionGenerator();
        tmpGenerator.withSize(300, 350).withMolTitle().withTitleColor(Color.BLACK);
        BufferedImage tmpImg = tmpGenerator.depict(tmpTestMol).toImg();
        ImageIcon tmpIcon = new ImageIcon(tmpImg);
        JFrame frame= new JFrame();
        frame.setLayout(new FlowLayout());
        frame.setSize(1000,800);
        JLabel tmpOriLbl = new JLabel();
        tmpOriLbl.setIcon(tmpIcon);
        tmpOriLbl.setText("Original");
        frame.add(tmpOriLbl);
        //Generate fragments, rings and frameworks
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter(false,1);
        tmpMurckoFragmenter.setComputeRingFragments(true);
        tmpMurckoFragmenter.generateFragments(tmpTestMol);
        IAtomContainer[] tmpFragments = tmpMurckoFragmenter.getFragmentsAsContainers();
        IAtomContainer[] tmpRings = tmpMurckoFragmenter.getRingSystemsAsContainers();
        IAtomContainer[] tmpFrameworks = tmpMurckoFragmenter.getFrameworksAsContainers();
        //Generate pictures of the fragments
        int tmpCountFra = 0;
        for(IAtomContainer tmpFragment : tmpFragments) {
            tmpCountFra++;
            BufferedImage tmpImgFra = tmpGenerator.withBackgroundColor(Color.LIGHT_GRAY).depict(tmpFragment).toImg();
            ImageIcon tmpIconFra = new ImageIcon(tmpImgFra);
            JFrame frameFra = new JFrame();
            frameFra.setLayout(new FlowLayout());
            frameFra.setSize(320,370);
            JLabel tmpLblFra = new JLabel();
            tmpLblFra.setIcon(tmpIconFra);
            tmpLblFra.setText("Fragment "+tmpCountFra);
            frame.add(tmpLblFra);
        }
        //Generate pictures of the rings
        int tmpCountRgs = 0;
        for(IAtomContainer tmpRing : tmpRings) {
            tmpCountRgs++;
            BufferedImage tmpImgRgs = tmpGenerator.withBackgroundColor(Color.GRAY).depict(tmpRing).toImg();
            ImageIcon tmpIconRgs = new ImageIcon(tmpImgRgs);
            JFrame frameRgs = new JFrame();
            frameRgs.setLayout(new FlowLayout());
            frameRgs.setSize(320,370);
            JLabel tmpLblRgs = new JLabel();
            tmpLblRgs.setIcon(tmpIconRgs);
            tmpLblRgs.setText("Ring "+tmpCountRgs);
            frame.add(tmpLblRgs);
        }
        //Generate pictures of the frameworks
        int tmpCountFrw = 0;
        for(IAtomContainer tmpFramework : tmpFrameworks) {
            tmpCountFrw++;
            BufferedImage tmpImgFrw = tmpGenerator.withBackgroundColor(Color.PINK).depict(tmpFramework).toImg();
            ImageIcon tmpIconFrw = new ImageIcon(tmpImgFrw);
            JFrame frameFrw = new JFrame();
            frameFrw.setLayout(new FlowLayout());
            frameFrw.setSize(320,370);
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
