package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import at.cibiv.codoc.CoverageTools;
import at.cibiv.ngs.tools.lds.GenomicInterval;
import at.cibiv.ngs.tools.util.BinnedHistogram;
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
	// /**
	// * Calculates the base-coverage per feature in a passed UCSC refseq data
	// * file. The output will be sorted by genomic coordinates. FIXME:
	// * experimental method - use with care.
	// *
	// * @param covFile
	// * @param covVcfFile
	// * @param bedFile
	// * @param covOut
	// * @throws Throwable
	// */
	// public static void doGeneAnalysis2(File covFile, File covVcfFile, File
	// ucscFlatFile, PrintStream covOut) throws Throwable {
	// CoverageDecompressor cov = null;
	//
	// try {
	// if (debug)
	// System.out.println("Load " + covFile + " / " + ucscFlatFile);
	//
	// RefSeqDb db = new RefSeqDb(ucscFlatFile);
	// Map<Gene, Float> geneCoverageMap = new HashMap<Gene, Float>();
	// Map<Gene, Float> genePosUnmappedMap = new HashMap<Gene, Float>();
	// for (Gene gene : db.getGenesSorted()) {
	// for (ExonInterval ex : gene.getExonsSorted())
	// ex.setAnnotation("gene", gene);
	// // System.out.println( gene );
	// }
	//
	// cov = CoverageDecompressor.loadFromFile(covFile, covVcfFile);
	// CoverageTools.debug = false;
	// if (debug)
	// System.out.println(Arrays.toString(cov.getChromosomes().toArray()));
	//
	// covOut.println("min/max are gene positions; cov, width, avg.cov is calc only over exons.");
	// covOut.println("#chr\tmin\tmax\tname\tstrand\twidth\tcov\tavg.cov\tperc.uncov");
	//
	// Map<String, ITree<Long>> exonTrees = db.buildItrees();
	// Map<String, ITree<Long>> intronTrees = db.buildIntronItrees();
	//
	// CompressedCoverageIterator it = cov.getCoverageIterator();
	// // list of chromosome names that could not be matched.
	// Set<String> ignoredChromosomes = new HashSet<String>();
	//
	// while (it.hasNext()) {
	// Float h = it.next();
	// float coverage = (h == null) ? 0f : h;
	// GenomicPosition pos = it.getGenomicPosition();
	// if (pos == null) {
	// break;
	// }
	// ExonChromosomeTree tree = (ExonChromosomeTree)
	// exonTrees.get(pos.getChromosome());
	// if (tree == null)
	// tree = (ExonChromosomeTree) exonTrees.get(pos.getChromosomeOriginal());
	// if (tree == null) {
	// ignoredChromosomes.add(pos.getChromosome());
	// continue;
	// }
	// List<Interval<Long>> res = tree.query(pos.get0Position());
	// for (Interval<Long> ex : res) {
	// Gene g = (Gene) ex.getAnnotation("gene");
	// Float geneCov = geneCoverageMap.get(g);
	// if (geneCov == null)
	// geneCov = 0f;
	// geneCov += coverage;
	// geneCoverageMap.put(g, geneCov);
	//
	// if (geneCov == 0f) {
	// Float posUnmapped = genePosUnmappedMap.get(g);
	// if (posUnmapped == null)
	// posUnmapped = 0f;
	// posUnmapped++;
	// geneCoverageMap.put(g, geneCov);
	// }
	// }
	// }
	//
	// if (ignoredChromosomes.size() > 0)
	// System.err.println("Ignored chromosomes: " + ignoredChromosomes);
	//
	// for (Gene gene : db.getGenesSorted()) {
	// Float geneCov = geneCoverageMap.get(gene);
	// Float posUnmapped = genePosUnmappedMap.get(gene);
	// covOut.print(toStr(gene.getChromosome()) + "\t");
	// covOut.print(toStr(gene.getMin()) + "\t");
	// covOut.print(toStr(gene.getMax()) + "\t");
	// covOut.print(toStr(gene.getName()) + "\t");
	// covOut.print(toStr(gene.getStrand()) + "\t");
	// covOut.print(toFloat(gene.getWidth(), 0) + "\t");
	// covOut.print(toFloat(geneCov, 0) + "\t");
	// if (geneCov == null)
	// covOut.print("n/a\t");
	// else
	// covOut.print(toFloat((gene.getWidth() == null ? 0f : (geneCov /
	// gene.getWidth())), 1) + "\t");
	// covOut.print(toFloat((posUnmapped == null ? 0f : (posUnmapped /
	// gene.getWidth())), 1) + "\t");
	// covOut.println();
	// }
	//
	// } finally {
	// if (cov != null)
	// cov.close();
	// }
	// if (debug)
	// System.out.println("Finished.");
	//
	// }

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
	public static void doGeneAnalysis(File covFile, File covVcfFile, File bedFile, File outFile, File histFile, int binSize) throws Throwable {

		BinnedHistogram bh = new BinnedHistogram(binSize);

		List<File> covFiles = new ArrayList<File>();
		covFiles.add(covFile);

		System.out.println("calculate interval coverage for " + bedFile);
		Map<GenomicInterval, Map<File, Double>> intervalsAnnotated = CoverageTools.calculateCoveragePerBedFeature(covFiles, null, bedFile, null);

		Map<String, double[]> geneCoverageMap = new HashMap<String, double[]>();
		Map<String, GenomicInterval> geneMap = new HashMap<String, GenomicInterval>();

		for (GenomicInterval iv : intervalsAnnotated.keySet()) {

			double avgCov = intervalsAnnotated.get(iv).get(covFile);
			String geneName = iv.getUri().substring(0, iv.getUri().lastIndexOf("_"));
			iv.setAnnotation("gene", geneName);

			double exonWidth = iv.getWidth() + 1;
			double exonCov = avgCov * exonWidth;

			double[] cov = geneCoverageMap.get(geneName);

			if (cov == null)
				cov = new double[3];
			cov[0] += exonCov; // total cov
			cov[1] += exonWidth; // total width
			cov[2] += 1; // num of found exons
			geneCoverageMap.put(geneName, cov);

			GenomicInterval gene = geneMap.get(geneName);
			if (gene == null)
				gene = iv;
			String chr = gene.getOriginalChrom();
			gene = new GenomicInterval(gene.getOriginalChrom(), Math.min(gene.getMin(), iv.getMin()), Math.max(gene.getMax(), iv.getMax()), geneName);
			gene.setOriginalChrom(chr);
			geneMap.put(geneName, gene);
			// System.out.println(iv + "\t" + iv.getAnnotation("gene") + "\t" +
			// avgCov);
		}

		// write out.
		PrintStream out = new PrintStream(outFile);
		out.println("#gene\tpos\ttotal-cov\ttotal-width\tfound-exons\tav.cov");
		for (String geneName : geneCoverageMap.keySet()) {
			double[] cov = geneCoverageMap.get(geneName);
			GenomicInterval gene = geneMap.get(geneName);
			double avgCov = 0f;
			out.print(geneName + "\t" + gene.toCoordString() + "\t");
			if (cov == null) {
				out.println("no exon\tno exon\tno exon\tno exon");
			} else {
				avgCov = (cov[0] / cov[1]);
				out.println(cov[0] + "\t" + cov[1] + "\t" + cov[2] + "\t" + avgCov);
			}

			bh.push((int) Math.round(avgCov));
		}
		out.close();

		// write hist.
		PrintStream histOut = new PrintStream(histFile);
		histOut.println(bh.toCSV());
		histOut.close();

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
	 * @param exons
	 * @param introns
	 * @param merged
	 * @throws IOException
	 */
	public static void mergeExonIntronTable(File exons, File introns, File merged ) throws IOException {
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
			if ( ! exonScore.keySet().contains(geneName))
				throw new IOException("No exon entry for gene " + geneName);
			Double avcov = Double.parseDouble(t[5]);
			intronScore.put(geneName, avcov);
		}
			
		PrintStream out = new PrintStream(merged);
		for ( String geneName : exonScore.keySet() ) {
			Double exScore = exonScore.get(geneName);
			Double inScore = intronScore.get(geneName);
			if ( inScore == null ) inScore = -100d;
			out.println(geneName + "\t" + exScore + "\t" + inScore );
		}
		out.close();
		System.exit(0);
	}
	

	/**
	 * @param args
	 * @throws Throwable
	 */
	public static void main(String[] args) throws Throwable {

		
		mergeExonIntronTable(
				new File( "/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-ISS.csv"),
				new File( "/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-INTRONS-ISS.csv"),
				new File( "/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-MERGED-ISS.csv")
				);	
//		mergeExonIntronTable(
//				new File( "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT-EXONS-ISS.csv"),
//				new File( "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT-INTRONS-ISS.csv"),
//				new File( "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT-MERGED-ISS.csv")
//				);
		
		// smallCheck();
		// if ( 1==1) return;
		// args = new String[] {
		// "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-GENOME.ISS.wig.bw.wig.codoc",
		// "/project/ngs-work/meta/annotations/exons/hg19/refseq/niko/UCSC-RefSeq-genes-exons-only-sorted.bed",
		// "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/delme.csv",
		// "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/delme.stats.csv",
		// };

//		args = new String[] { "/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-GENOME.ISS.wig.codoc",
//				"/project/oesi/genomicAmbiguity/droso/anno/droso_R5_dm3_ucsc.EXONS.sorted.bed",
//				"/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-ISS.csv",
//				"/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-EXONS-ISS.hist.csv", };

//		 args = new String[] {
//		 "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-GENOME.ISS.wig.bw.wig.codoc",
//		 "/project/ngs-work/meta/annotations/exons/hg19/refseq/niko/UCSC-RefSeq-genes-introns-only-sorted.bed",
//		 "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT-INTRONS-ISS.csv",
//		 "/project2/oesi/genAmb/output/hg19-NO_SCORE_LIMIT/hg19-NO_SCORE_LIMIT-INTRONS-ISS.hist.csv",
//		 "100"
//		 };


		File codocFile = new File(args[0]);
		File codocVcfFile = null;
		File exons = new File(args[1]);
		File outFile = new File(args[2]);
		File csvFile = new File(args[3]);
		int binSize = Integer.parseInt(args[4]);
		doGeneAnalysis(codocFile, codocVcfFile, exons, outFile, csvFile, binSize);

	}

}
