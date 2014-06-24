package at.cibiv.argos.mapq;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import at.cibiv.codoc.CompressedCoverageIterator;
import at.cibiv.codoc.CoverageDecompressor;
import at.cibiv.codoc.utils.CodocException;
import at.cibiv.ngs.tools.util.StringUtils;

/**
 * compares mapping qualities to ISS/AMB values.
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class MapqAnalyzer {

	final static int GRID = 10;
	private static boolean debug = false;

	
	
	/**
	 * Extracts a table that lists average mapping quality per position and the CODOC score for this
	 * position. Only reads from the passed chromosome list will be included. Only GRID positions will
	 * be examined.
	 * @param sortedSamFile
	 * @param issCodoc
	 * @param out
	 * @param chrs
	 * @param maxReads
	 * @throws IOException
	 * @throws CodocException
	 */
	public static void analyzeMapq(File sortedSamFile, File issCodoc, PrintStream out, List<String> chrs, Integer maxReads) throws IOException, CodocException {
		SAMFileReader inputSam = new SAMFileReader(sortedSamFile);
		inputSam.setValidationStringency(ValidationStringency.LENIENT);
		SAMRecordIterator it = inputSam.iterator();

		Map<String, Map<Integer, List<Integer>>> scores = new HashMap<String, Map<Integer, List<Integer>>>();
		ArrayList<String> prefixedChrs = null;
		if (chrs != null) {
			prefixedChrs = new ArrayList<String>();
			for (String c : chrs)
				prefixedChrs.add(StringUtils.prefixedChr(c));
		}

		String oldChr = null;
		Map<Integer, List<Integer>> chrMap = null;
		int c = 0;
		while (it.hasNext()) {
			SAMRecord rec = it.next();
			int mapq = rec.getMappingQuality();
			// skip secondary alignments (occurs, e.g., with bwa mem skipped-read alignment) 
			if ( rec.getNotPrimaryAlignmentFlag() )
				continue;
			// skip reads not on the supported chroms
			String chr = StringUtils.prefixedChr(rec.getReferenceName());
			if (prefixedChrs != null && !prefixedChrs.contains(chr))
				continue;

			if (oldChr == null || !chr.equals(oldChr)) {
				if (chrMap != null) {
					scores.put(oldChr, chrMap);
				}
				chrMap = scores.get(chr);
				if (chrMap == null) {
					chrMap = new HashMap<Integer, List<Integer>>();
				}
				oldChr = chr;

			}
			int start = (rec.getAlignmentStart() / GRID) * GRID;
			int end = (rec.getAlignmentEnd() / GRID) * GRID;
			for (int x = start; x < end; x += GRID) {
				List<Integer> pos = chrMap.get(x);
				if (pos == null)
					pos = new ArrayList<Integer>();
				pos.add(mapq);
				chrMap.put(x, pos);
				// out.println(chr + ":" + x + " -> " + mapq);
			}
			if ((debug) && (++c % 100000 == 0))
				System.out.println("scores: " + c);
			if (maxReads != null)
				if (c > maxReads)
					break;
		}
		scores.put(oldChr, chrMap);

		System.out.println("Chromosomes: " + scores.keySet());

		CoverageDecompressor decomp = null;
		StringBuffer buf = new StringBuffer();
		try {
			decomp = CoverageDecompressor.decompress(issCodoc, null);
			CompressedCoverageIterator issIt = decomp.getCoverageIterator();
			c = 0;
			while (issIt.hasNext()) {
				Double cov = issIt.nextCoveragePrecise();
				String chr = StringUtils.prefixedChr(issIt.getCurrentReference());

				int off = (int) issIt.getOffset();
				if (off % GRID != 0)
					continue;
				chrMap = scores.get(chr);
				if (chrMap != null) {
					List<Integer> pos = chrMap.get(off);
					if (pos != null) {
						double avgMapq = 0d;
						for (Integer mapq : pos)
							avgMapq += mapq;
						avgMapq = avgMapq / (double) pos.size();

						if (avgMapq > 60) {
							float sum = 0f;
							for (Integer mapq : pos)
								sum += mapq;
							//System.out.println(sum + "/" + (double) pos.size() + "/" + avgMapq);
						}

						buf.append(chr + ":" + off + "\t" + avgMapq + "\t" + cov + "\n");
					}
				}

				if (buf.length() > 1000000) {
					out.print(buf.toString());
					buf = new StringBuffer();
				}
				if ((debug) && (++c % 1000000 == 0))
					System.out.println("output: " + c);

			}
			out.print(buf.toString());

			System.out.println("Finished.");

		} finally {
			if (decomp != null)
				decomp.close();
		}

	}

	/**
	 * Print usage information.
	 * 
	 * @param options
	 */
	private static void usage(Options options, String e) {

		HelpFormatter hf = new HelpFormatter();
		hf.setLeftPadding(10);
		hf.setDescPadding(2);
		hf.setWidth(160);
		hf.setSyntaxPrefix("Usage:    ");
		hf.printHelp("java -jar x.jar " + MapqAnalyzer.class + ", Params:", options, true);

		if (e != null)
			System.out.println("\nError: " + e);

		System.exit(1);
	}

	/**
	 * @param args
	 * @throws IOException
	 * @throws CodocException
	 */
	public static void main(String[] args) throws IOException, CodocException {
		// args = new String[] {
		// "/project/oesi/Project_Oesi2/GIST/GIST_1_N/work/GIST_1_N.ngm-FINAL.bam",
		// "/project/oesi/genomicAmbiguity/hg19-genomereads-100bp-step5.bam.data.wig.p0.codoc",
		// "/project/oesi/genomicAmbiguity/mapq-analysis/GIST_1_N.ngm-FINAL.mapq.ISS.dat"
		// };
		// args = new String[] {
		// "/project/ngs-work/philipp/hfq-pd/reads-ngm-chr/3xF-PD-chr.bam",
		// "/project2/oesi/genAmb/output/ecK12/eck12_MG1655_ecoli-chr-GENOME.AMB.wig.codoc",
		// "/project/ngs-work/philipp/hfq-pd/reads-ngm-chr/mapq-analysis/3xF-PD-chr.MAPQ_vs_AMB.dat"
		// };

		// args = new String[] {
		// "-bam",
		// "/project2/oesi/genAmb-OLD/test_data/human_wgs/mapping/bwamem-ERR125588.sort.bam",
		// "-codoc",
		// "/project2/oesi/genAmb/output/hg19/hg19-GENOME.AMB.wig.codoc",
		// "-o",
		// "/project2/oesi/genAmb-OLD/test_data/human_wgs/mapping/bwamem-ERR125588.sort.bam.MAPQ_vs_AMB.dat",
		// "-maxreads", "1000000", "-v" };

		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption("h", "help", false, "Print this usage information.");

		try {
			// parse the command line arguments
			CommandLine line = parser.parse(options, args, true);

			Option o = new Option("bam", true, "Input BAM file");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("codoc", true, "Input CODOC file");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("o", true, "Output file");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("chr", true, "List of chromosomes that should be included");
			o.setRequired(false);
			options.addOption(o);

			o = new Option("maxreads", true, "Maximum reads (optional)");
			o.setRequired(false);
			options.addOption(o);

			options.addOption("v", "verbose", false, "be verbose.");

			// validate that block-size has been set
			if (line.getArgs().length == 0 || line.hasOption("h")) {
				usage(options, null);
			}

			line = parser.parse(options, args);

			if (line.hasOption("v"))
				debug = true;
			else
				debug = false;

			ArrayList<String> chr = null;
			for (String c : line.getOptionValues("chr")) {
				if (chr == null)
					chr = new ArrayList<String>();
				chr.add(c);
			}
			Integer maxReads = (line.hasOption("maxreads") ? Integer.parseInt(line.getOptionValue("maxreads")) : null);

			analyzeMapq(new File(line.getOptionValue("bam")), new File(line.getOptionValue("codoc")), new PrintStream(line.getOptionValue("o")), chr, maxReads);

		} catch (Exception e) {
			e.printStackTrace();
			usage(options, e.toString());
		}

	}

}
