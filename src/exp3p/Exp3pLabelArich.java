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
import java.util.BitSet;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.BufferedReaderMaker;
import jsequtils.sequence.FastaReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * Exp3pLabelArich is one of the main utilities of Exp3p. It processes a bam file 
 * and labels reads that are possibly due to capture by A-rich DNA/RNA regions,
 * i.e. reads whose ends are very close to poly-A or poly-T stretches in the reference
 * genome.
 *
 * @author tkonopka
 */
public class Exp3pLabelArich implements Runnable {

    // variables with values determined from the command line
    private boolean isok = false;
    private File bamfile;
    private File genomefile = null;
    private String output = "out.bam";
    private int windowlen = 18;
    private int mincount = 12;        
    // private internal constants for book-keeping
    private static int indexA = 0, indexT = 1, indexC = 2, indexG = 3;
    // remove determines whether AR hits are removed from the output bam or not
    private boolean remove = false;

    public void printHelp() {
        System.out.println("Exp3p labelArich: labels reads close to A rich regions");
        System.out.println();
        System.out.println("Usage: java -jar Exp3p.jar labelArich [options]");
        System.out.println();
        System.out.println("  --bam <File>          - alignment file to process ");
        System.out.println("  --output <String>     - output file");
        System.out.println("  --genome <File>       - genome fasta file");
        System.out.println("  --window <int>        - window length [default "+windowlen+"]");
        System.out.println("  --mincount <int>      - minimum number of A bases in window [default "+mincount+"]");        
        System.out.println("  --remove <boolean>    - remove A rich reads [default "+remove+"]");        
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@ludwig.ox.ac.uk)\n");
    }

    private boolean parseParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // bam - input bam file
        prs.accepts("bam").withRequiredArg().ofType(File.class);
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // output - output file
        prs.accepts("output").withRequiredArg().ofType(String.class);

        // window and minA are for detecting A-rich regions
        prs.accepts("window").withRequiredArg().ofType(Integer.class);
        prs.accepts("mincount").withRequiredArg().ofType(Integer.class);
        prs.accepts("remove").withRequiredArg().ofType(Boolean.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // determine preferences for how to catch A rich regions
        if (options.has("window")) {
            windowlen = (int) (Integer) options.valueOf("window");
            if (windowlen <= 0) {
                System.out.println("Parameter window must greater than 0");
                return false;
            }
        }

        if (options.has("mincount")) {
            mincount = (int) (Integer) options.valueOf("mincount");
            if (mincount < 0) {
                System.out.println("Parameter mincount must be non-negative");
                return false;
            }
        }                       

        // check for inputs for alignment/labels
        if (options.has("bam")) {
            bamfile = (File) options.valueOf("bam");
            if (!bamfile.exists() || !bamfile.canRead()) {
                System.out.println("Cannot read bam file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument: bam");
            return false;
        }

        // check for inputs for alignment/labels
        if (options.has("genome")) {
            genomefile = (File) options.valueOf("genome");
            if (!genomefile.exists() || !genomefile.canRead()) {
                System.out.println("Cannot read genome file, or file does not exist");
                return false;
            }
        } else {
            System.out.println("Missing required argument: genome");
            return false;
        }

        // figure out where to write (bam) output
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        if (options.has("remove")) {
            remove = (Boolean) options.valueOf("remove");
        }

        return true;
    }

    public Exp3pLabelArich(String[] args) {
        if (args == null || args.length == 0) {
            printHelp();
            return;
        }
        if (parseParameters(args)) {
            isok = true;
        }
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        FastaReader fr;
        try {
            fr = new FastaReader(BufferedReaderMaker.makeBufferedReader(genomefile));
        } catch (IOException ex) {
            System.out.println("Trouble setup: " + ex.getMessage());
            return;
        }

        try {
            processBam(bamfile, output, fr);
        } catch (Exception ex) {
            System.out.println("Trouble processing bam file: " + ex.getMessage());
        }

        try {
            fr.close();
        } catch (Exception ex) {
            System.out.println("Trouble closing genome: " + ex.getMessage());
        }
    }

    private void processBam(File bamfile, String output, FastaReader fr) throws IOException {

        // open the input for reading
        SAMFileReader bamreader = new SAMFileReader(bamfile);
        bamreader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMFileWriter outSam;
        SAMFileHeader outHeader = bamreader.getFileHeader();
        outHeader.addComment("Exp3 - Marked reads in  AT rich regions (" + mincount + " A bases in window of size " + windowlen + ")");
        outSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, true, new File(output));

        // prep - determine starting chromosome, etc.
        int curchrindex = -1;
        SAMRecord record;

        BitSet bsA = null;
        BitSet bsT = null;

        int readcount = 0, readARpositive = 0;

        // start to iterate the 
        SAMRecordIterator bamiterator = bamreader.iterator();
        while (bamiterator.hasNext()) {
            record = bamiterator.next();
            record.setAttribute("AR", 0);
            readcount++;

            int recordRefIndex = record.getReferenceIndex();
            if (recordRefIndex >= 0 && !record.getReadUnmappedFlag()) {
                if (recordRefIndex != curchrindex) {
                    String recordRefName = record.getReferenceName();
                    while (!recordRefName.equals(fr.getChromosomeName()) && fr.hasNext()) {
                        fr.readNext(true);
                    }
                    if (recordRefName.equals(fr.getChromosomeName())) {
                        bsA = findRichRegions(fr, windowlen, indexA, mincount);
                        bsT = findRichRegions(fr, windowlen, indexT, mincount);
                    } else {
                        System.out.println("Cannot find chromosome " + recordRefName + " in genome.");
                        break;
                    }
                    curchrindex = recordRefIndex;
                }

                // check the read coordinates against the AT-rich region bitsets
                if (bsA != null && bsT != null) {
                    // perform a different check if a read is on minus or plus strand
                    if (record.getReadNegativeStrandFlag()) {
                        if (bsA.get(record.getAlignmentEnd() - 1)) {
                            record.setAttribute("AR", 1);
                            readARpositive++;
                        }
                    } else {
                        if (bsT.get(record.getAlignmentStart() - 1)) {
                            record.setAttribute("AR", 1);
                            readARpositive++;
                        }
                    }

                } else {
                    System.out.println("BitSets are null! This should not happen - something is wrong.");
                    break;
                }

            }

            // output the record (but only if it is good or user wants to see all reads)
            if (!remove || (Integer) record.getAttribute("AR") == 0) {
                outSam.addAlignment(record);
            }
        }

        outSam.close();
        bamreader.close();

        System.out.println("TotalReads\t" + readcount);
        System.out.println("ARpositive\t" + readARpositive);
        System.out.println("ARnegative\t" + (readcount - readARpositive));
    }

    /**
     *
     * book-keeping function. Counts the number of bases of a particular type
     *
     * @param counts
     * @param base
     * @param change
     */
    private static void updateCounts(int[] counts, byte base, int change) {
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
     *
     * make a bitset labeling regions with high content in one nucleotide
     *
     * @param fr
     *
     * a reader for the genome.
     *
     * @param indexA
     *
     * use one of indexA, indexT, indexC, indexG.
     *
     * @param indexG
     *
     * use one of indexA, indexT, indexC, indexG.
     *
     *
     * @return
     */
    private static BitSet findRichRegions(FastaReader fr, int windowlen, int indexA, int minA) {

        int chrlen = fr.getChromosomeLength();
        BitSet bs = new BitSet(chrlen);

        int[] counts = new int[4];

        // deal with the start of the chromosome
        // (this requires extra checking for negative indexes)
        for (int i = 0; i < chrlen && i < windowlen; i++) {
            // shift the window one to the right
            updateCounts(counts, fr.getBaseAtPositionBase0(i), 1);
            int imin = i - windowlen;
            if (imin >= 0) {
                updateCounts(counts, fr.getBaseAtPositionBase0(imin), -1);
            }
            if (counts[indexA] >= minA) {
                bs.set(Math.max(0, imin + 1), i + 1);
            }
        }

        // loop over the inside of the chromosome, less checking involved here
        for (int i = windowlen; i < chrlen; i++) {
            // shift the window one to the right
            updateCounts(counts, fr.getBaseAtPositionBase0(i), 1);
            int imin = i - windowlen;
            updateCounts(counts, fr.getBaseAtPositionBase0(imin), -1);
            if (counts[indexA] >= minA) {
                bs.set(imin, i + 1);
            }
        }

        return bs;
    }
}
