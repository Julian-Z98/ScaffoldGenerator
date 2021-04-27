/*
 *
 *  * MIT License
 *  *
 *  * Copyright (c) 2021 Julian Zander, Jonas Schaub,  Achim Zielesny
 *  *
 *  * Permission is hereby granted, free of charge, to any person obtaining a copy
 *  * of this software and associated documentation files (the "Software"), to deal
 *  * in the Software without restriction, including without limitation the rights
 *  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  * copies of the Software, and to permit persons to whom the Software is
 *  * furnished to do so, subject to the following conditions:
 *  *
 *  * The above copyright notice and this permission notice shall be included in all
 *  * copies or substantial portions of the Software.
 *  *
 *  * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  * SOFTWARE.
 *
 */

package de.unijena.cheminf.scaffoldTest;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;
import java.io.IOException;

/**
 * The method murckoFragmenter.getMurckoFraments is tested with the 9 files(Test1-Test9) from the resources folder.
 */
public class murckoFragmenterTest {

    private murckoFragmenter fragmenter;
    @Before
    public void initFragmenter() {
        fragmenter = new murckoFragmenter();
    }
    @Test
    public void testFragmenter() throws IOException, CDKException {
        for (int tmpCount = 1; tmpCount < 10; tmpCount++) {
            fragmenter.getMurckoFragments("Test" + tmpCount);
        }
    }
}
