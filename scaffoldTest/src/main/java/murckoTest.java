import org.openscience.cdk.Atom;
import org.openscience.cdk.fragment.MurckoFragmenter;

public class murckoTest {
    public static void main(String[] args) {
        System.out.println("Hello World"); //Erster Test
        Atom atom = new Atom(6);
        atom.setCharge(21.0);
        MurckoFragmenter tmpMurckoFragmenter = new MurckoFragmenter();
        System.out.println(atom.getCharge()); //Erster Test
    }
}
