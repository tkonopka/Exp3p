/*
 * Copyright 2014 Tomasz Konopka.
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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;
import joptsimple.OptionParser;
import joptsimple.OptionSet;
import jsequtils.file.OutputStreamMaker;
import jsequtils.genome.GenomeInfo;
import jsequtils.genome.GenomicPositionComparator;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import static net.sf.samtools.CigarOperator.D;
import static net.sf.samtools.CigarOperator.I;
import static net.sf.samtools.CigarOperator.M;
import static net.sf.samtools.CigarOperator.N;
import static net.sf.samtools.CigarOperator.S;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

/**
 * Utility scans a bam file and calls pre-poly-A sites, i.e. either
 * transcription start sites for RNA-seq or sites before DNA poly-A stretches in
 * DNA-seq.
 *
 *
 * @author tkonopka
 */
public class Exp3pCallPPA implements Runnable {

    // variables with values determined from the command line   
    private String command;
    private File annofile = null, bedfile = null;
    private File genomefile = null;
    private ArrayList<File> bamfiles;
    private ArrayList<String> labels;
    private String output = "stdout";
    // whether to filter out labeled reads
    private boolean filterAR = false;
    // isok determines if class is initialized correctly
    private boolean isok = false;
    private int depththreshold = 5;
    private int readlen = 50;

    public void printCallPPAHelp() {
        System.out.println("Exp3p callPPA: call pre-poly-A (PPA) sites from exp3p sequencing data");
        System.out.println();
        System.out.println("Usage: java -jar Exp3p.jar callPPA [options]");
        System.out.println();
        System.out.println("  --anno <File>         - refGene-style annotation of transcripts");
        System.out.println("  --bam <File>          - alignment file to process ");
        System.out.println("                          (can specify multiple times)");
        System.out.println("  --bed <File>          - regions file (alternative to --anno and --ppa)");
        System.out.println("  --depth <int>         - depth threshld for calling PPA [default " + depththreshold + "]");
        System.out.println("  --genome <File>       - genome file");
        System.out.println("  --label <String>      - sample identifier");
        System.out.println("                          (must specify one label for each bam file)");
        System.out.println("  --filterAR <boolean>  - filter out reads labeld by AR tag [default " + filterAR + "]");
        //System.out.println("  --stranded <String>   - stranded protocol, use NONE (default), SAME, or REVERSE");
        System.out.println("  --output <String>     - file where output will be stored");
        System.out.println("  --readlen <int>       - read length [default " + readlen + "]");
        System.out.println();
        System.out.println("\nAuthor: Tomasz Konopka (tkonopka@cemm.oeaw.ac.be)\n");
    }

    private boolean parseCallPPAParameters(String[] args) {
        OptionParser prs = new OptionParser();

        // bam - input bam file
        prs.accepts("bam").withRequiredArg().ofType(File.class);
        // samplelabel - samplelabel that appear in vcf file
        prs.accepts("label").withRequiredArg().ofType(String.class);

        // for labeling ppa sites as belonging to regions or genes
        prs.accepts("anno").withRequiredArg().ofType(File.class);
        prs.accepts("bed").withRequiredArg().ofType(File.class);

        // threshold for calling ppa sites
        prs.accepts("depth").withRequiredArg().ofType(Integer.class);

        // stranded - stranded protocol
        //prs.accepts("stranded").withRequiredArg().ofType(String.class);

        // genome - genome annotation file
        prs.accepts("genome").withRequiredArg().ofType(File.class);

        // output - output file
        prs.accepts("output").withRequiredArg().ofType(String.class);

        // if reads labeled by the AR (A-rich) tag should be ignored or not
        prs.accepts("filterAR").withRequiredArg().ofType(Boolean.class);

        // now use OptionSet to parse the command line
        OptionSet options;
        try {
            options = prs.parse(args);
        } catch (Exception ex) {
            System.out.println("Error parsing command line parameters\n" + ex.getMessage());
            return false;
        }

        // determine the regions that the user wants to look at
        if (options.has("anno") || options.has("bed")) {
            if (options.has("anno") && options.has("bed")) {
                System.out.println("Cannot use both --bed and --anno");
                return false;
            }
            if (options.has("anno")) {
                annofile = (File) options.valueOf("anno");
                if (!annofile.exists() || !annofile.canRead()) {
                    System.out.println("Cannot read annotation file: " + annofile.getAbsolutePath());
                    return false;
                }
            }
            if (options.has("bed")) {
                bedfile = (File) options.valueOf("bed");
                if (!bedfile.exists() || !bedfile.canRead()) {
                    System.out.println("Cannot read bed file: " + bedfile.getAbsolutePath());
                    return false;
                }
            }
        }

        // choose strandedness protocol        

        // check for inputs for alignment/labels
        if (options.has("bam")) {
            bamfiles = new ArrayList<File>((List<File>) options.valuesOf("bam"));
        } else {
            System.out.println("Missing required argument: bam");
            return false;
        }
        if (options.has("label")) {
            labels = new ArrayList<String>((List<String>) options.valuesOf("label"));
        } else {
            System.out.println("Missing required argument: label");
            return false;
        }

        // the number of input files and labels must match to uniquely associated         
        if (labels.size() != bamfiles.size()) {
            System.out.println("number of bam and label arguments do not match");
            return false;
        }

        // figure out where to write the output
        if (options.has("output")) {
            output = (String) options.valueOf("output");
        }

        // figure out where to write the output
        if (options.has("genome")) {
            genomefile = (File) options.valueOf("genome");
            if (!genomefile.exists() || !genomefile.canRead()) {
                System.out.println("genome file is not readable");
                return false;
            }
        } else {
            System.out.println("Missing required argument: genome");
            return false;
        }

        if (options.has("filterAR")) {
            filterAR = (Boolean) options.valueOf("filterAR");
        }

        if (options.has("depth")) {
            depththreshold = (Integer) options.valueOf("depth");
            if (depththreshold < 1) {
                System.out.println("Parameter depth must be positive");
                return false;
            }
        }

        if (options.has("readlen")) {
            readlen = (Integer) options.valueOf("readlen");
            if (readlen < 1) {
                System.out.println("Parameter readlen must be positive");
                return false;
            }
        }

        return true;
    }

    public Exp3pCallPPA(String[] args) {
        if (args == null || args.length == 0) {
            printCallPPAHelp();
            return;
        }

        if (parseCallPPAParameters(args)) {
            isok = true;
        }

        // record the command given into a variable
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < args.length; i++) {
            sb.append(" ").append(args[i]);
        }
        command = sb.toString();
    }

    @Override
    public void run() {
        if (!isok) {
            return;
        }

        GenomeInfo ginfo = null;
        try {
            ginfo = new GenomeInfo(genomefile);
        } catch (IOException ex) {
            System.out.println("Error reading genome info: " + ex.getMessage());
            return;
        }

        Exp3pAnnotationMap emap = null;
        try {
            if (annofile != null) {
                emap = new Exp3pAnnotationMap(annofile, Exp3pAnnotationMap.AnnoFileTypeEnum.REFGENE, Integer.MAX_VALUE, 0);
            } else if (bedfile != null) {
                emap = new Exp3pAnnotationMap(bedfile, Exp3pAnnotationMap.AnnoFileTypeEnum.BED, Integer.MAX_VALUE, 0);
            }
        } catch (Exception ex) {
            System.out.println("Exception while loading annotation: " + ex.getMessage());
            return;
        }

        // now loop through all the bam files, collect relevant expression info,
        // and report expression sites into lists
        ArrayList<GenomicInterval> sitesplus = new ArrayList<GenomicInterval>();
        ArrayList<GenomicInterval> sitesminus = new ArrayList<GenomicInterval>();
        for (int i = 0; i < bamfiles.size(); i++) {
            callInOneBam(emap, bamfiles.get(i), sitesplus, sitesminus, readlen);
        }

        // convert called positions into PPA sites        
        ArrayList<Exp3pPPAInterval> ppaintervals = makePPAintervals(emap, sitesplus, sitesminus, ginfo);

        // print out all the recorded expression information                 
        try {
            outputPPAIntervalCalls(emap, ppaintervals, ginfo);
        } catch (Exception ex) {
            System.out.println("Exception while outputing expression values: " + ex.getMessage());
        }
    }

    /**
     * Scans one bam file, records coverage in memory and then,
     * chromosome-by-chromosome, find peaks that mark pre-poly-A intervals.
     *
     * @param emap
     * @param infile
     * @param fileindex
     * @param sitesplus
     * @param sitesminus
     * @param readlen
     */
    private void callInOneBam(Exp3pAnnotationMap emap, File infile,
            ArrayList<GenomicInterval> sitesplus, ArrayList<GenomicInterval> sitesminus, int readlen) {

        // set up input/output objects using SAM library
        SAMFileReader inputSam = new SAMFileReader(infile);
        inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        SAMSequenceDictionary seqdict = inputSam.getFileHeader().getSequenceDictionary();

        // keep track of current chromosome 
        int lastchr = -1;
        String lastchrname = "NA";

        // make objects for coverage, they will be redefined for each chromosome
        KilobaseIntArray covplus = null;
        KilobaseIntArray covminus = null;

        for (final SAMRecord record : inputSam) {
            // avoid unmapped reads, but otherwise process each record at a time
            if (!record.getReadUnmappedFlag()) {
                removeSplice(record);
                int nowchr = record.getReferenceIndex();
                String nowchrname = record.getReferenceName();

                // figure out whether to use this read or not, based on its AR tag            
                boolean useThisRead = true;
                if (filterAR) {
                    Object ARtag = record.getAttribute("AR");
                    if (ARtag != null) {
                        if ((Integer) ARtag == 1) {
                            useThisRead = false;
                        }
                    }
                }

                // check if list of transcript at this locus should be reloaded            
                if (nowchr != lastchr || covplus == null) {

                    if (lastchr >= 0) {
                        callFromTrack(sitesplus, lastchrname, covplus, true, readlen);
                        callFromTrack(sitesminus, lastchrname, covminus, false, readlen);
                    }
                    int chromlength = seqdict.getSequence(nowchrname).getSequenceLength();
                    covplus = new KilobaseIntArray(chromlength);
                    covminus = new KilobaseIntArray(chromlength);
                }

                if (useThisRead) {
                    updateCovTracks(record, covplus, covminus);
                }

                lastchr = nowchr;
                lastchrname = nowchrname;
            }
        }

        // call from the last chromosome
        if (lastchr >= 0) {
            callFromTrack(sitesplus, lastchrname, covplus, true, readlen);
            callFromTrack(sitesminus, lastchrname, covminus, false, readlen);
        }

        inputSam.close();
    }

    private void callFromTrack(ArrayList<GenomicInterval> sitelist, String chr,
            KilobaseIntArray track, boolean plusstrand, int readlen) {

        // bookkeeping devices
        int lastcov = 0;
        int tocall = 0;
        int tracklen = track.size();

        // on plus strand coverage, call PPA sites from the left hand side of coverage peaks
        // on minus strand coverage, called PPA sites from right side of coverage peaks        
        if (plusstrand) {
            for (int i = 0; i < tracklen; i++) {
                int nowcov = track.get(i);
                if (nowcov == 0) {
                    tocall = 0;
                }
                // check if site has sufficient coverage to be a candidate
                if (nowcov - tocall >= depththreshold) { // && lastcov < depththreshold) {
                    // look ahead to see if there is a more covered position ahead
                    int maxcov = nowcov;
                    int maxdiff = nowcov;
                    if (i > 0) {
                        maxdiff = nowcov - track.get(i - 1);
                    }
                    int maxpos = i;
                    int maxdiffpos = i;
                    lastcov = nowcov;
                    for (int j = 1; j < readlen && i + j < tracklen; j++) {
                        int temp = track.get(i + j);
                        int nowdiff = temp - lastcov;
                        if (temp >= nowcov && nowdiff >= 0) {
                            if (temp >= maxcov) {
                                maxcov = temp;
                                maxpos = i + j;
                            }
                            if (nowdiff > maxdiff) {
                                maxdiff = nowdiff;
                                maxdiffpos = i + j;
                            }
                        } else {
                            // break the inner for loop
                            j = readlen;
                        }
                        lastcov = temp;
                    }
                    // define the maximal position as the peak/PPA site
                    int sitestart = Math.min(maxdiffpos, maxpos);
                    int siteend = Math.min(tracklen, Math.max(Math.max(maxdiffpos, maxpos) + 1, sitestart + readlen));
                    sitelist.add(new GenomicInterval(chr, sitestart, siteend));
                    // skip the i counter ahead
                    i = Math.max(maxdiffpos, maxpos);
                    // make sure that after a successful PPA find, set the threshold for next one high
                    tocall = track.get(i);
                }
                lastcov = track.get(i);
            }
        } else {
            for (int i = tracklen - 1; i >= 0; i--) {
                int nowcov = track.get(i);
                if (nowcov == 0) {
                    tocall = 0;
                }
                // check if site has sufficient coverage to be a candidate
                if (nowcov - tocall >= depththreshold) { // && lastcov < depththreshold) {
                    // look ahead to see if there is more covered position, i.e. a stronger peak position                     
                    int maxcov = nowcov;
                    int maxdiff = nowcov;
                    if (i < tracklen - 1) {
                        maxdiff = nowcov - track.get(i + 1);
                    }
                    int maxpos = i;
                    int maxdiffpos = i;
                    lastcov = nowcov;
                    for (int j = 1; j < readlen && i - j >= 0; j++) {
                        int temp = track.get(i - j);
                        int nowdiff = temp - lastcov;
                        if (temp >= nowcov && nowdiff >= 0) {
                            if (temp >= maxcov) {
                                maxcov = temp;
                                maxpos = i - j;
                            }
                            if (nowdiff > maxdiff) {
                                maxdiff = nowdiff;
                                maxdiffpos = i - j;
                            }
                        } else {
                            // break out from the for loop
                            j = readlen;
                        }
                        lastcov = temp;
                    }
                    int siteend = Math.max(maxpos, maxdiffpos) + 1;
                    int sitestart = Math.max(1, Math.min(siteend - readlen, Math.min(maxpos, maxdiffpos)));
                    sitelist.add(new GenomicInterval(chr, sitestart, siteend));
                    //skip the i counter ahead
                    i = Math.min(maxdiffpos, maxpos);
                    tocall = track.get(i);
                }
                lastcov = track.get(i);
            }
        }
    }

    /**
     *
     * assumes the intervals have been sorted.
     *
     * @param intervals
     *
     * a set of intervals
     *
     * @return
     *
     * a new list in which intervals do not overlap with one another
     *
     */
    private ArrayList<GenomicInterval> mergeOverlappingIntervals(ArrayList<GenomicInterval> intervals) {
        int isize = intervals.size();
        ArrayList<GenomicInterval> ans = new ArrayList<GenomicInterval>(isize);

        if (isize == 0) {
            return ans;
        }

        GenomicInterval lastinterval = new GenomicInterval(intervals.get(0));
        ans.add(lastinterval);
        for (int i = 1; i < isize; i++) {
            if (!lastinterval.mergeWith(intervals.get(i))) {
                lastinterval = new GenomicInterval(intervals.get(i));
                ans.add(lastinterval);
            }
        }

        return ans;
    }

    private ArrayList<Exp3pPPAInterval> makePPAintervals(Exp3pAnnotationMap emap,
            ArrayList<GenomicInterval> sitesplus, ArrayList<GenomicInterval> sitesminus,
            GenomeInfo ginfo) {

        ArrayList<Exp3pPPAInterval> ans = new ArrayList<Exp3pPPAInterval>();
        GenomicPositionComparator gpc = new GenomicPositionComparator(ginfo);

        // prep by sorting the sites
        Collections.sort(sitesplus, gpc);
        Collections.sort(sitesminus, gpc);

        // used to merge overlapping intervals - is that really necessary?
        //sitesplus = mergeOverlappingIntervals(sitesplus);
        //sitesminus = mergeOverlappingIntervals(sitesminus);

        // now go through the plus sites and create the annotations
        if (sitesplus.size() > 0) {
            for (int i = 0; i < sitesplus.size(); i++) {
                GenomicInterval nowinterval = sitesplus.get(i);
                ans.addAll(makePPAIntervalsFromOneInterval(emap, nowinterval, true));
            }
        }

        // repeat for the sites called on negative strand
        if (sitesminus.size() > 0) {
            for (int i = 0; i < sitesminus.size(); i++) {
                GenomicInterval nowinterval = sitesminus.get(i);
                ans.addAll(makePPAIntervalsFromOneInterval(emap, nowinterval, false));
            }
        }

        // for good measure, sort the output
        // (it should already be sorted, but just to make sure...)
        Collections.sort(ans, gpc);

        return ans;
    }

    /**
     * Use a called position together with annotations to create an object
     * describing a PPA site
     *
     *
     * @param emap
     * @param ginterval
     * @param plusstrand
     *
     * true if reads contributing to the interval are aligned on plus strand
     *
     * @return
     */
    private ArrayList<Exp3pPPAInterval> makePPAIntervalsFromOneInterval(Exp3pAnnotationMap emap,
            GenomicInterval ginterval,
            boolean plusstrand) {

        ArrayList<Exp3pPPAInterval> ans = new ArrayList<Exp3pPPAInterval>();
        String gposchr = ginterval.getChr();

        // get a list of transcripts that match the genomic interval
        ArrayList<String> rawhits = new ArrayList<String>(4);
        rawhits.addAll(emap.getValuesFull(gposchr, ginterval.getPosition()));
        rawhits.addAll(emap.getValuesFull(gposchr, ginterval.getEnd()));
        Collections.sort(rawhits);

        // deal with the simple case first, if the peak does not match with any declared genic region
        if (rawhits.isEmpty()) {
            ans.add(makeOneBlankPPAInterval(ginterval, plusstrand));
            return ans;
        }

        // replace transcript names with gene names
        for (int i = 0; i < rawhits.size(); i++) {
            String nowtx = rawhits.get(i);
            String nowgene = emap.getGeneFromTx(nowtx);
            rawhits.set(i, nowtx.charAt(0) + nowgene);
        }
        Collections.sort(rawhits);

        // get unique strings
        ArrayList<String> hits = new ArrayList<String>();
        hits.add(rawhits.get(0));
        for (int i = 1; i < rawhits.size(); i++) {
            String nowh = rawhits.get(i);
            if (!nowh.equals(rawhits.get(i - 1))) {
                hits.add(nowh);
            }
        }

        // process each hit, without worrying about uniqueness now        
        int numadded = 0;
        if (plusstrand) {
            for (int i = 0; i < hits.size(); i++) {
                String nowgene = hits.get(i);
                char nowstrand = nowgene.charAt(0);
                nowgene = nowgene.substring(1);
                if (nowstrand == '-') {
                    Exp3pPPAInterval newinterval = new Exp3pPPAInterval(ginterval);
                    newinterval.setStrand('-');
                    newinterval.setLabel(nowgene + ":" + gposchr + ":" + ginterval.getPosition());
                    ans.add(newinterval);
                    numadded++;
                }
            }
        } else {
            for (int i = 0; i < hits.size(); i++) {
                String nowgene = hits.get(i);
                char nowstrand = nowgene.charAt(0);
                nowgene = nowgene.substring(1);
                if (nowstrand == '+') {
                    Exp3pPPAInterval newinterval = new Exp3pPPAInterval(ginterval);
                    newinterval.setStrand('+');
                    newinterval.setLabel(nowgene + ":" + gposchr + ":" + ginterval.getEnd());
                    ans.add(newinterval);
                    numadded++;
                }
            }
        }

        if (numadded == 0) {
            ans.add(makeOneBlankPPAInterval(ginterval, plusstrand));
        }

        return ans;
    }

    /**
     * make a single Exp3pPPAInterval without gene name.
     *
     * @param ginterval
     * @param plusstrand
     * @return
     */
    private Exp3pPPAInterval makeOneBlankPPAInterval(GenomicInterval ginterval, boolean plusstrand) {
        String gposchr = ginterval.getChr();
        Exp3pPPAInterval ans = new Exp3pPPAInterval(ginterval);
        if (plusstrand) {
            ans.setStrand('-');
            ans.setLabel(".:" + gposchr + ":" + (ginterval.getStart() + 1));
        } else {
            ans.setStrand('+');
            ans.setLabel(".:" + gposchr + ":" + ginterval.getEnd());
        }

        return ans;
    }

    /**
     * Adjust coverage using coordinates from one SAM record.
     *
     * The mappability multiplicity is ignored in this function.
     *
     *
     * @param record
     * @param covplus
     * @param covminus
     */
    static void updateCovTracks(SAMRecord record, KilobaseIntArray covplus, KilobaseIntArray covminus) {

        int[] recordpos = getGenomicPositions(record);
        KilobaseIntArray usecov;
        if (record.getReadNegativeStrandFlag()) {
            usecov = covminus;
        } else {
            usecov = covplus;
        }

        // increment coverage counts
        int rpsize = recordpos.length;
        for (int i = 0; i < rpsize; i++) {
            int nowpos = recordpos[i] - 1;
            if (nowpos >= 0) {
                usecov.increment(nowpos);
            }
        }
    }

    /**
     * Removes splice elements from the cigar and move the alignment start
     * position accordingly. For reads on the plus strand, the alignment start
     * position is unchanged. For reads on the negative strand, the alignment
     * end position is unchanged.
     *
     * Thus, the 5prime end of the reads are always in the same place as before
     * the function call
     *
     * The cigar string will not be properly formed after this function, e.g. may contain multiple M elements
     * in a row, eg. 20M100N30M -> 20M30M. However, this should not be a problem for applications within Exp3p.
     * 
     * @param record
     */
    static void removeSplice(SAMRecord record) {

        Cigar cigar = record.getCigar();

        ArrayList<CigarElement> oldElements = new ArrayList<CigarElement>();
        oldElements.addAll(cigar.getCigarElements());

        // check if one of the cigar elements is a splice
        boolean hassplice = false;
        for (int i = 0; i < oldElements.size(); i++) {
            if (oldElements.get(i).getOperator() == CigarOperator.N) {
                hassplice = true;
            }
        }
        if (!hassplice) {
            return;
        }

        // create a new cigar array and keep track of shift in start position
        ArrayList<CigarElement> newElements = new ArrayList<CigarElement>();
        int shiftstart = 0;

        for (int i = 0; i < oldElements.size(); i++) {
            if (oldElements.get(i).getOperator() != CigarOperator.N) {
                newElements.add(oldElements.get(i));
            } else {
                shiftstart += oldElements.get(i).getLength();
            }
        }

        // set the new cigar
        record.setCigar(new Cigar(newElements));
        // for negative strand , also set new start position
        if (record.getReadNegativeStrandFlag()) {
            record.setAlignmentStart(record.getAlignmentStart() + shiftstart);
        }

        // that's it, the changes will be reflected in the record object        
    }

    /**
     *
     * @param record
     *
     * @return
     *
     * an array the length of the read. Elements hold the genomic coordinates
     * (1-based) of the alignment
     *
     */
    static int[] getGenomicPositions(SAMRecord record) {

        Cigar cigar = record.getCigar();
        int cigarsize = cigar.numCigarElements();
        cigar.getCigarElements();

        int offset = 0;
        int nowpos = record.getAlignmentStart();
        int[] ans = new int[record.getBaseQualities().length];

        for (int i = 0; i < cigarsize; i++) {
            CigarElement ce = cigar.getCigarElement(i);
            int celen = ce.getLength();

            switch (ce.getOperator()) {
                case M:
                    // for matches transfer position in series into the answer array
                    for (int j = 0; j < celen; j++) {
                        ans[offset] = nowpos;
                        offset++;
                        nowpos++;
                    }
                    break;
                case I:
                    // for insertion, label bases as aligning onto negative positions
                    // this indicates an insertion and also where the insertion starts
                    for (int j = 0; j < celen; j++) {
                        ans[offset] = -nowpos;
                        offset++;
                    }
                    break;
                case D:
                    // for deletions, skip ahead 
                    nowpos += celen;
                    break;
                case S:
                    // for soft-clips, that is like insertions
                    for (int j = 0; j < celen; j++) {
                        ans[offset] = -nowpos;
                        offset++;
                    }
                    break;
                case N:
                    // for splices, like deletions, skip ahead in the genome
                    nowpos += celen;
                    break;
                case H:
                    // for hard clipping, don't do anything
                    break;
                default:
                    // Anything here? I think I have all the options covered
                    System.out.println("record: " + record.getReadName() + " (" + record.getReferenceName() + ":" + record.getAlignmentStart() + ")");
                    System.out.println("problematic cigar: " + record.getCigarString());
                    break;
            }
        }
        return ans;
    }

    /**
     *
     *
     * @return
     *
     * A string that contains multiple lines, all beginning with # (as a comment
     * on a table) The comment contains information about the samples
     *
     *
     *
     */
    private String getOutputHeader() {
        StringBuilder sb = new StringBuilder();
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd hh:mm:ss");


        sb.append("# exp3p callppa ").append("\n");
        sb.append("# ").append(sdf.format(new Date())).append("\n");
        sb.append("#\n");
        sb.append("# command:\tcallppa ").append(command).append("\n");
        sb.append("#\n");
        if (bedfile == null) {
            sb.append("# annotation\t").append(annofile.getAbsolutePath()).append("\n");
        } else {
            sb.append("# bed\t").append(bedfile.getAbsolutePath()).append("\n");
        }
        sb.append("# genome\t").append(genomefile).append("\n");
        int numsamples = bamfiles.size();
        for (int i = 0; i < numsamples; i++) {
            sb.append("# Sample\t").append(labels.get(i)).append("\t").append(bamfiles.get(i).getAbsolutePath()).append("\n");
        }
        sb.append("# depth\t").append(depththreshold).append("\n");
        sb.append("# readlen\t").append(readlen).append("\n");
        sb.append("#\n");

        return sb.toString();
    }

    /**
     * Output the calls for pre-poly-A sites
     *
     * @param emap
     * @throws FileNotFoundException
     * @throws IOException
     */
    private void outputPPAIntervalCalls(Exp3pAnnotationMap emap,
            ArrayList<Exp3pPPAInterval> ppaintervals,
            GenomeInfo ginfo) throws FileNotFoundException, IOException {

        // startyp
        StringBuilder sb;
        OutputStream os = OutputStreamMaker.makeOutputStream(output);

        // output the header
        sb = new StringBuilder();
        sb.append(getOutputHeader());
        sb.append("# Chr\tStart\tEnd\tStrand.Region.id\n");
        os.write(sb.toString().getBytes());

        GenomicPositionComparator gpc = new GenomicPositionComparator(ginfo);
        Collections.sort(ppaintervals, gpc);

        for (int i = 0; i < ppaintervals.size(); i++) {
            Exp3pPPAInterval nn = ppaintervals.get(i);
            os.write((nn.toString() + "\n").getBytes());
        }

        os.close();
    }
}
