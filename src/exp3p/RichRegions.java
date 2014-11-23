/*
 * Copyright 2013 Tomasz Konopka.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package exp3p;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.BitSet;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.file.OutputStreamMaker;
import jsequtils.sequence.FastaReader;

/**
 *
 * @author tkonopka
 */
public class RichRegions implements Runnable {

    private boolean isok = false;
    private final int[] mincounts = new int[4];
    private int windowlen = 10;
    private File genomefile = null;
    private String output = "stdout";
    private static int indexA = 0, indexT = 1, indexC = 2, indexG = 3;

    public static void printRichRegionsHelp() {
        System.out.println("Exp3p richregions: find regions rich in some nucleotides (e.g. A-rich regions)");
        System.out.println();
        System.out.println("Usage: java -jar Exp3p.jar richregions [options]");
        System.out.println();
        System.out.println("  --genome <File>     - fasta file for genome ");
        System.out.println("  --window <int>      - length of window");
        System.out.println("  --minA <int>        - minimum count of A");
        System.out.println("  --minT <int>        - minimum count of T");
        System.out.println("  --minC <int>        - minimum count of C");
        System.out.println("  --minG <int>        - minimum count of G");
        System.out.println("  --output <String>   - file where output will be stored (default is stdout)");
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@cemm.oeaw.ac.be)\n");
    }

    private boolean parseRichRegionsParameters(String[] args) {

        OptionParser prs = new OptionParser();

        prs.accepts("genome").withRequiredArg().ofType(File.class);
        prs.accepts("output").withRequiredArg().ofType(String.class);

        prs.accepts("window").withRequiredArg().ofType(Integer.class);
        prs.accepts("minA").withRequiredArg().ofType(Integer.class);
        prs.accepts("minT").withRequiredArg().ofType(Integer.class);
        prs.accepts("minC").withRequiredArg().ofType(Integer.class);
        prs.accepts("minG").withRequiredArg().ofType(Integer.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // determine the richness criteria
        if (options.has("minA")) {
            mincounts[indexA] = (int) (Integer) options.valueOf("minA");
        }
        if (options.has("minT")) {
            mincounts[indexT] = (int) (Integer) options.valueOf("minT");
        }
        if (options.has("minC")) {
            mincounts[indexC] = (int) (Integer) options.valueOf("minC");
        }
        if (options.has("minG")) {
            mincounts[indexG] = (int) (Integer) options.valueOf("minG");
        }

        // determine the sliding window length
        if (options.has("window")) {
            windowlen = (int) (Integer) options.valueOf("window");
        }

        // get the genome to search
        if (options.has("genome")) {
            genomefile = (File) options.valueOf("genome");
            if (!genomefile.exists() || !genomefile.canRead()) {
                System.out.println("Cannot read genome file: " + genomefile.getAbsolutePath());
                return false;
            }
        }

        // output the 
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        return true;
    }

    RichRegions(String[] args) {
        if (args == null || args.length == 0) {
            printRichRegionsHelp();
            return;
        }
        if (parseRichRegionsParameters(args)) {
            isok = true;
        }
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        FastaReader fr;
        OutputStream os;
        try {
            fr = new FastaReader(BufferedReaderMaker.makeBufferedReader(genomefile));
            os = OutputStreamMaker.makeOutputStream(output);
        } catch (IOException ex) {
            System.out.println("Trouble setup: " + ex.getMessage());
            return;
        }
        
        while (fr.hasNext()) {
            try {
                fr.readNext(true);                
                BitSet bs = findRichRegions(fr, os);
                outputBed(bs, fr.getChromosomeName(), os);
            } catch (IOException ex) {
                System.out.println("Trouble reading genome: " + ex.getMessage());
                return;
            }
        }

        try {
            fr.close();
            os.close();
        } catch (Exception ex) {
            System.out.println("Trouble closing genome: " + ex.getMessage());
        }

    }

    private void updateCounts(int[] counts, byte base, int change) {
        switch (base) {
            case 'A':
                counts[indexA] += change;
                return;
            case 'T':
                counts[indexT] += change;
                return;
            case 'C':
                counts[indexC] += change;
                return;
            case 'G':
                counts[indexG] += change;
        }
    }

    /**
     * Compare base counts to thresholds.
     *
     * @param counts
     *
     * must be array of length 4 or more.
     *
     * @param thresholds
     *
     * must be array of length 4 or more.
     *
     * @return
     *
     * true if all the thresholds are satisfied. false otherwise.
     *
     */
    private boolean decideOnRichRegion(int[] counts, int[] thresholds) {        
        for (int i = 0; i < 4; i++) {
            if (counts[i] < thresholds[i]) {
                return false;
            }
        }
        return true;
    }

    private BitSet findRichRegions(FastaReader fr, OutputStream os) {

        int chrlen = fr.getChromosomeLength();
        BitSet bs = new BitSet(chrlen);

        int[] counts = new int[4];

        for (int i = 0; i < chrlen; i++) {
            // shift the window one to the right
            updateCounts(counts, fr.getBaseAtPositionBase0(i), 1);
            int imin = i - windowlen;
            if (imin >= 0) {
                updateCounts(counts, fr.getBaseAtPositionBase0(imin), -1);
            }
            // check if the interval satisfies richness criteria
            if (decideOnRichRegion(counts, mincounts)) {                
                bs.set(Math.max(0, imin+1), i+1);
            }
        }

        return bs;
    }

    private void outputBed(BitSet bs, String chr, OutputStream os) throws IOException {
        
        int nextset = bs.nextSetBit(0);        
        while (nextset >= 0) {
            int nextclear = bs.nextClearBit(nextset);  
            
            if (nextclear < 0) {
                nextclear = bs.size();
            }
            String out = chr + "\t" + nextset + "\t" + nextclear + "\n";
            os.write(out.getBytes());

            if (nextclear < 0) {
                nextset = -1;
            } else {
                nextset = bs.nextSetBit(nextclear);
            }
        }

    }
}
