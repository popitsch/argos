package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import at.cibiv.codoc.CompressedCoverageIterator;
import at.cibiv.codoc.CoverageDecompressor;
import at.cibiv.ngs.tools.exon.Gene;
import at.cibiv.ngs.tools.exon.RefSeqDb;
import at.cibiv.ngs.tools.lds.GenomicInterval;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.Statistics;
import at.cibiv.ngs.tools.util.StringUtils;
import at.cibiv.ngs.tools.util.TabIterator;

public class ArgosGeneAnalysis {
	public static boolean debug = true;

	// private static String toStr(Object o) {
	// if (o == null)
	// return "n/a";
	// return o + "";
	// }
	//
	// private static String toFloat(Float f, int comma) {
	// if (f == null)
	// return "n/a";
	// return String.format("%." + comma + "f", f);
	// }
	//
	/**
	 * Calculates the base-coverage per feature in a passed UCSC refseq data
	 * file. The output will be sorted by genomic coordinates. FIXME:
	 * experimental method - use with care.
	 * 
	 * @param covFile
	 * @param covVcfFile
	 * @param bedFile
	 * @param covOut
	 * @throws Throwable
	 */
	public static void doGeneAnalysis2(File covFile, File covVcfFile, File ucscFlatFile, PrintStream out, PrintStream statsOut) throws Throwable {
		CoverageDecompressor cov = null;
		Statistics stats = new Statistics();

		try {
			if (debug)
				System.out.println("Load " + covFile + " / " + ucscFlatFile);

			RefSeqDb db = new RefSeqDb(ucscFlatFile);

			// get a sorted and annotated list of all genomic intervals.
			List<GenomicInterval> allIntervalsSorted = new ArrayList<>();
			for (Gene g : db.getGenesSorted()) {
				Map<String, List<GenomicInterval>> data = g.getUtrExIntronLists();
				allIntervalsSorted.addAll(data.get("5pUTR"));
				allIntervalsSorted.addAll(data.get("3pUTR"));
				allIntervalsSorted.addAll(data.get("exon"));
				allIntervalsSorted.addAll(data.get("intron"));
			}

			Collections.sort(allIntervalsSorted);
			// create lists [per chrom.
			Map<String, List<GenomicInterval>> allIntervalsSortedPerChrom = new HashMap<String, List<GenomicInterval>>();
			for (GenomicInterval g : allIntervalsSorted) {
				List<GenomicInterval> l = allIntervalsSortedPerChrom.get(g.getOriginalChrom());
				if (l == null)
					l = new ArrayList<>();
				l.add(g);
				allIntervalsSortedPerChrom.put(g.getOriginalChrom(), l);
			}
			Double igSum = 0d;
			Double igWidth = 0d;
			List<String> consideredChromosomes = new ArrayList<>();

			// annotate intervals with score sum.
			try {
				cov = CoverageDecompressor.loadFromFile(covFile, covVcfFile);
				for (String cc : cov.getChromosomes())
					consideredChromosomes.add(StringUtils.prefixedChr(cc));
				if (debug)
					System.out.println("Handled chroms: " + Arrays.toString(cov.getChromosomes().toArray()));

				CompressedCoverageIterator it = cov.getCoverageIterator();

				String currentChr = null;
				Iterator<GenomicInterval> bfit = null;
				GenomicInterval nextInterval = null;
				List<GenomicInterval> overlapping = new ArrayList<GenomicInterval>();

				while (it.hasNext()) {

					Float coverage = it.next();
					GenomicPosition pos = it.getGenomicPosition();
					int category = 0;

					if (currentChr == null || !currentChr.equals(pos.getChromosomeOriginal())) {
						overlapping.clear();
						currentChr = pos.getChromosomeOriginal();
						if (allIntervalsSortedPerChrom.get(currentChr) == null) {
							System.err.println("Cannot load any intervals for " + currentChr);
							bfit = null;
							nextInterval = null;
						} else {
							bfit = allIntervalsSortedPerChrom.get(currentChr).iterator();
							nextInterval = null;
							System.err.println("switch to " + currentChr);
						}
					}

					// del intervals that are finished
					Iterator<GenomicInterval> delit = overlapping.iterator();
					while (delit.hasNext()) {
						GenomicInterval iv = delit.next();
						if (!iv.contains(pos)) {
							delit.remove();
							// System.out.println("fin" + iv);
						}
					}

					// calc overlapping intervals
					if (bfit != null) {
						if (nextInterval == null && bfit.hasNext())
							nextInterval = bfit.next();
						if (nextInterval != null)
							while (nextInterval.contains(pos)) {
								overlapping.add(nextInterval);
								if (!bfit.hasNext()) {
									nextInterval = null;
									break;
								} else {
									nextInterval = bfit.next();
								}
							}
					}
					// increment score sums.
					for (GenomicInterval gi : overlapping) {
						String type = (String) gi.getAnnotation("type");
						if (type.equals("exon"))
							category = Math.max(category, 4);
						else if (type.equals("5pUTR"))
							category = Math.max(category, 3);
						else if (type.equals("3pUTR"))
							category = Math.max(category, 2);
						else if (type.equals("intron"))
							category = Math.max(category, 1);

						Double sum = (Double) gi.getAnnotation("scoreSum");
						if (sum == null)
							sum = 0d;
						sum += coverage;
						gi.setAnnotation("scoreSum", sum);
					}
					// score intergenic regions
					if (overlapping.size() == 0) {
						igSum += coverage;
						igWidth++;
					}

					// stats
					switch (category) {
					case 4:
						stats.inc("Count-Exon");
						break;
					case 3:
						stats.inc("Count-5pUTR");
						break;
					case 2:
						stats.inc("Count-3pUTR");
						break;
					case 1:
						stats.inc("Count-Intron");
						break;
					default:
						stats.inc("Count-IG");
						break;
					}
					stats.inc("Count-ALL");

				}
			} finally {
				if (cov != null)
					cov.close();
			}

			// annotate genes with total score sum/width per category
			if (debug)
				System.out.println("Annotate genes");
			SortedSet<Gene> annotatedGenes = new TreeSet<Gene>();
			for (GenomicInterval g : allIntervalsSorted) {
				Double isum = (Double) g.getAnnotation("scoreSum");
				if (isum == null)
					isum = 0d;
				Gene gene = (Gene) g.getAnnotation("gene");
				String type = (String) g.getAnnotation("type");

				Double genesum = (Double) gene.getAnnotation("scoreSum" + type);
				if (genesum == null) {
					// outside of considered range? => skip!
					if (!consideredChromosomes.contains(StringUtils.prefixedChr(g.getChr()))) {
						if (debug)
							System.err.println("Skipped gene " + g + " as not considered chromosome");
						continue;
					}

					genesum = 0d;
				}
				genesum += isum;
				gene.setAnnotation("scoreSum" + type, genesum);

				Double genewidth = (Double) gene.getAnnotation("scoreWidth" + type);
				if (genewidth == null)
					genewidth = 0d;
				genewidth += (g.getWidth() + 1);
				gene.setAnnotation("scoreWidth" + type, genewidth);

				annotatedGenes.add(gene);
			}

			// output all genes
			if (debug)
				System.out.println("Output results: " + annotatedGenes.size() + " genes and " + allIntervalsSorted.size() + " intervals");
			out.print("Gene\tPos\t");
			out.print("Ex-AvgScore\tEx-Width\tEx-Sum\t");
			out.print("In-AvgScore\tIn-Width\tIn-Sum\t");
			out.print("3pUTR-AvgScore\t3pUTR-Width\t3pUTR-Sum\t");
			out.print("5pUTR-AvgScore\t5pUTR-Width\t5pUTR-Sum\t");
			out.print("IG-AvgScore\tIG-Width\tIG-Sum");
			out.println();
			Double igAvg = igSum / igWidth;

			String[] types = new String[] { "exon", "intron", "3pUTR", "5pUTR" };
			for (Gene g : annotatedGenes) {
				out.print(g.getName() + "\t" + g.toPosition1());
				for (String type : types) {
					Double genesum = (Double) g.getAnnotation("scoreSum" + type);
					Double genewidth = (Double) g.getAnnotation("scoreWidth" + type);
					if ((genesum == null) || (genewidth == null))
						out.print("\t-\t-\t-");
					else {
						// System.out.println(genesum + "?"+genewidth);
						out.print("\t" + (genesum / genewidth) + "\t" + genewidth + "\t" + genesum);
					}
				}
				out.println("\t" + igAvg + "\t" + igWidth + "\t" + igSum);
			}

			// print stats
			stats.toFile(statsOut);

		} finally {
			if (cov != null)
				cov.close();
		}
		if (debug)
			System.out.println("Finished.");

	}

	/**
	 * For debugging. Counts the number of distinct gene names in two files...
	 * 
	 * @throws IOException
	 */
	public static void smallCheck() throws IOException {

		String[] fns = new String[] { "/project/ngs-work/meta/annotations/exons/hg19/refseq/niko/UCSC-RefSeq-genes-exons-only.bed",
				"/project/ngs-work/meta/annotations/exons/hg19/refseq/niko/UCSC-RefSeq-genes-introns-only.bed" };

		for (String fn : fns) {

			TabIterator ti = new TabIterator(new File(fn), "#");
			String last = null;
			Set<String> seen = new HashSet<>();
			while (ti.hasNext()) {
				String[] t = ti.next();
				String geneName = t[3].substring(0, t[3].lastIndexOf("_"));
				if (last == null)
					last = geneName;
				if (!geneName.equals(last)) {
					if (seen.contains(geneName))
						System.err.println("Dupl name: " + geneName);
					seen.add(geneName);
					last = geneName;
				}
			}
			System.out.println("seen " + seen.size() + " gene names in " + fn);
		}
	}

	/**
	 * Creates a simple merged table
	 * 
	 * @param exons
	 * @param introns
	 * @param merged
	 * @throws IOException
	 */
	public static void mergeExonIntronTable(File exons, File introns, File merged) throws IOException {
		Map<String, Double> exonScore = new HashMap<>();
		TabIterator ti = new TabIterator(exons, "#");
		while (ti.hasNext()) {
			String[] t = ti.next();
			String geneName = t[0];
			Double avcov = Double.parseDouble(t[5]);
			exonScore.put(geneName, avcov);
		}
		Map<String, Double> intronScore = new HashMap<>();
		ti = new TabIterator(introns, "#");
		while (ti.hasNext()) {
			String[] t = ti.next();
			String geneName = t[0];
			if (!exonScore.keySet().contains(geneName))
				throw new IOException("No exon entry for gene " + geneName);
			Double avcov = Double.parseDouble(t[5]);
			intronScore.put(geneName, avcov);
		}

		PrintStream out = new PrintStream(merged);
		for (String geneName : exonScore.keySet()) {
			Double exScore = exonScore.get(geneName);
			Double inScore = intronScore.get(geneName);
			if (inScore == null)
				inScore = -100d;
			out.println(geneName + "\t" + exScore + "\t" + inScore);
		}
		out.close();
		System.exit(0);
	}

	/**
	 * @param args
	 * @throws Throwable
	 */
	public static void main(String[] args) throws Throwable {

//		args = new String[] {
//				"/project2/oesi/genAmb/output/hg19/hg19-GENOME.ISS.wig.codoc",
//				"/project2/oesi/genAmb/smalltest.delme.txt",
//				"/project2/oesi/genAmb/smalltest.delme.txt.csv"
//		};
//		
		
//		 args = new String[] {
//		 "/project2/oesi/genAmb/output/ecK12/eck12_MG1655_ecoli-chr-GENOME.AMB.wig.codoc",
//		 "/project2/oesi/genAmb/ref/ecK12/e.coli-K12-ucsc-refseq.txt",
//		 "/project2/oesi/genAmb/output/ecK12/eck12_MG1655_ecoli-chr-GENOME.AMB.FEATURES.csv"
//		 };
		 
//		 args = new String[] {
//				 "/project2/oesi/genAmb/output/dros_partial-rl100/droso_r5_partial-GENOME.AMB.wig.codoc",
//				 "/project2/oesi/genAmb/ref/dmel/droso_R5_dm3_ucsc.properchrom.txt",
//				 "/project2/oesi/genAmb/output/dros_partial-rl100/droso_r5_partial-GENOME.AMB.FEATURES.csv"
//				 };

		// File covFile = new
		// File("/project2/oesi/genAmb/output/dros_partial-rl50/droso_r5_partial-GENOME.AMB.wig.codoc");
		// File ucscFlatFile = new
		// File("/project2/oesi/genAmb/ref/dmel/droso_R5_dm3_ucsc.properchrom.txt");
		// PrintStream out = new
		// PrintStream("/project2/oesi/genAmb/output/dros_partial-rl50/droso_r5_partial-GENOME.AMB.wig.codoc.TEST.csv");

		if (args.length < 3) {
			System.err.println("Creates a FEATURE statistics file for the passed ISS/AMB/MSD signal using the passed gene annotations.");
			System.err.println("Usage: " + ArgosGeneAnalysis.class + " <codoc-score-file> <ucsc-flat-file> <output.csv> [<stats.txt>]");
			System.exit(1);
		}

		File covFile = new File(args[0]);
		File ucscFlatFile = new File(args[1]);
		PrintStream out = new PrintStream(args[2]);
		PrintStream stats = System.out;
		if (args.length > 3)
			stats = new PrintStream(args[3]);

		doGeneAnalysis2(covFile, null, ucscFlatFile, out, stats);

		out.close();

	}

}
