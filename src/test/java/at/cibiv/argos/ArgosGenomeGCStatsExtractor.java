package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import at.cibiv.ngs.tools.bed.SimpleBEDFile;
import at.cibiv.ngs.tools.fasta.FastaSequence;
import at.cibiv.ngs.tools.fasta.FastaSequenceIterator;
import at.cibiv.ngs.tools.fasta.MultiFastaSequence;
import at.cibiv.ngs.tools.lds.GenomicITree;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.GenomicPosition.COORD_TYPE;
import at.cibiv.ngs.tools.util.Statistics;

public class ArgosGenomeGCStatsExtractor {

    /**
     * Calculate base composition statistics within/without the given ROI for
     * the passed genome mfasta.
     * 
     * @param mfastaFile
     * @param out
     * @throws IOException
     */
    public static void calcBaseProfile(File mfastaFile, File bedRoiFile, PrintStream out) throws IOException {

	SimpleBEDFile roi = new SimpleBEDFile(bedRoiFile);
	GenomicITree git = roi.getGenomicITree();
	MultiFastaSequence mfasta = new MultiFastaSequence(mfastaFile, true);

	Statistics stats = new Statistics();

	out.println("# GC profile for " + mfastaFile + " within " + bedRoiFile);

	for (String chr : mfasta.getChromosomes()) {
	    stats.inc("Handled Chromosomes");
	    FastaSequence ff = mfasta.getSequence(chr);
	    ff.validate(null);
	    FastaSequenceIterator it = ff.iterator();
	    chr = mfasta.toOrig(chr);

	    long pos = 0;
	    while (it.hasNext()) {
		stats.inc("Overall positions");

		Character c = it.next();
		pos++;
		GenomicPosition gpos = new GenomicPosition(chr, pos, COORD_TYPE.ONEBASED);

		boolean inRoi = git.contains(gpos);

		if (inRoi) {
		    stats.inc("ROI positions");
		}

		if (c != null) {
		    if (c == 'a' || c == 'A')
			if (inRoi)
			    stats.inc(chr + "inROI-A");
			else
			    stats.inc(chr + "outROI-A");
		    else if (c == 'g' || c == 'G')
			if (inRoi)
			    stats.inc(chr + "inROI-G");
			else
			    stats.inc(chr + "outROI-G");
		    else if (c == 't' || c == 'T')
			if (inRoi)
			    stats.inc(chr + "inROI-T");
			else
			    stats.inc(chr + "outROI-T");
		    else if (c == 'c' || c == 'C')
			if (inRoi)
			    stats.inc(chr + "inROI-C");
			else
			    stats.inc(chr + "outROI-C");
		    else if (inRoi)
			stats.inc(chr + "inROI-N");
		    else
			stats.inc(chr + "outROI-N");
		}

	    }
	    // inc overall stats
	    System.out.println("Chr: " + chr);
	    stats.inc("inROI-A", stats.get(chr + "inROI-A"));
	    stats.inc("outROI-A", stats.get(chr + "outROI-A"));
	    stats.inc("inROI-G", stats.get(chr + "inROI-G"));
	    stats.inc("outROI-G", stats.get(chr + "outROI-G"));
	    stats.inc("inROI-C", stats.get(chr + "inROI-C"));
	    stats.inc("outROI-C", stats.get(chr + "outROI-C"));
	    stats.inc("inROI-T", stats.get(chr + "inROI-T"));
	    stats.inc("outROI-T", stats.get(chr + "outROI-T"));
	    stats.inc("inROI-N", stats.get(chr + "inROI-N", 0));
	    stats.inc("outROI-N", stats.get(chr + "outROI-N", 0));
	}

	out.println("----------------------------------------");
	out.println(stats);
	out.println("----------------------------------------");

	double inROIACTG = stats.get("inROI-A") + stats.get("inROI-C") + stats.get("inROI-T") + stats.get("inROI-G");
	double outROIACTG = stats.get("outROI-A") + stats.get("outROI-C") + stats.get("outROI-T") + stats.get("outROI-G");
	double overallACTG = inROIACTG + outROIACTG;
	double inROIGC = stats.get("inROI-C") + stats.get("inROI-G");
	double outROIGC = stats.get("outROI-C") + stats.get("outROI-G");
	double overallGC = inROIGC + outROIGC;
	double inROIGCPerc = inROIGC / inROIACTG;
	double outROIGCPerc = outROIGC / outROIACTG;
	double overallGCPerc = overallGC / overallACTG;

	out.format("GC %%\toverall (%.2f)\tin ROI (%.2f)\toutside ROI (%.2f)\n", overallGCPerc, inROIGCPerc, outROIGCPerc);
	System.out.println("Finished.");
    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

	if (args.length < 2) {
	    System.out.println("Usage: " + ArgosGenomeGCStatsExtractor.class + " <mfasta file> <bed file> ");
	    System.exit(1);
	}

	File mfastaFile = new File(args[0]);
	File bedRoiFile = new File(args[1]);

	// File mfastaFile = new
	// File("C:/data/genomicAmbiguity/droso/ref/droso_r5_partial.fasta");
	// File bedRoiFile = new
	// File("C:/data/genomicAmbiguity/droso/droso_r5_partial-GENOME.ISS.wig.codoc.above50.bed");

	calcBaseProfile(mfastaFile, bedRoiFile, new PrintStream(new File(bedRoiFile.getAbsolutePath() + ".GC-hist.txt")));

    }

}
