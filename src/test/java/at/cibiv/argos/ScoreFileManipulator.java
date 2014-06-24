package at.cibiv.argos;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.compress.compressors.gzip.GzipCompressorOutputStream;
import org.apache.commons.lang.ArrayUtils;

import at.cibiv.argos.mapq.MapqAnalyzer;
import at.cibiv.ngs.tools.util.StringUtils;
import at.cibiv.ngs.tools.util.TabIterator;

public class ScoreFileManipulator {

	private static boolean debug;

	public static void reduceScoreFile(File scoreFile, int rl, int step, List<String> chrs, File outFile) throws IOException {

		PrintStream out = new PrintStream(new GzipCompressorOutputStream(new FileOutputStream(outFile)));
		
		ArrayList<String> prefixedChrs = new ArrayList<String>();
		for (String c : chrs)
			prefixedChrs.add(StringUtils.prefixedChr(c));

		TabIterator ti = new TabIterator(scoreFile, null);
		Map<String, String> chromMap = new HashMap<String, String>();
		int c = 0;
		while (ti.hasNext()) {
			String[] t = ti.next();
			// parse header
			if (t[0].startsWith("#")) {
				out.println(StringUtils.concat(t, "\t"));
				t = ti.next();
				out.println(StringUtils.concat(t, "\t"));
				t[0] = t[0].substring(1);
				for (String tmp : t) {
					String[] x = tmp.split(":");
					chromMap.put(x[0], StringUtils.prefixedChr(x[1]));
				}
				continue;
			}
			// print contents
			String chr = StringUtils.prefixedChr(t[0].substring(t[0].indexOf("_")+1, t[0].indexOf(":")));
			if (!prefixedChrs.contains(chr)) {
				//System.out.println("no " + chr + " in " + ArrayUtils.toString(prefixedChrs));
				continue;
			}
			//System.out.println(t[0]);
			out.print(t[0]);
			for (int i = 1; i < t.length; i++) {
				String chrIdx = t[i].substring(0, t[i].indexOf(":"));
				String translatedChr = chromMap.get(chrIdx);
				//System.out.println(translatedChr + prefixedChrs.contains(translatedChr));
				if (prefixedChrs.contains(translatedChr))
					out.print("\t" + t[i]);
			}
			out.println();
			if (debug) {
				if (++c % 1000000 == 0)
					System.out.print(".");
				if (c % 100000000 == 0)
					System.out.println();
			}
		}
		out.close();
		System.out.println("Finished.");
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
	 */
	public static void main(String[] args) {

		args = new String[] { "-sc", "/project2/oesi/genAmb/output/hg19/hg19.scores.gz", "-rl", "200", "-step", "20", "-chr", "1", "-chr", "2", "-chr", "3",
				"-chr", "4", "-chr", "5", "-chr", "6", "-chr", "7", "-chr", "8", "-chr", "9", "-chr", "10", "-chr", "11", "-chr", "12", "-chr", "13", "-chr",
				"14", "-chr", "15", "-chr", "16", "-chr", "17", "-chr", "18", "-chr", "19", "-chr", "20", "-chr", "21", "-chr", "22", "-chr", "X", "-chr", "Y",
				"-o", "/project2/oesi/genAmb/output/hg19/hg19.scores.REDUCED-CHRS.gz", "-v" };

		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption("h", "help", false, "Print this usage information.");

		try {
			// parse the command line arguments
			CommandLine line = parser.parse(options, args, true);

			Option o = new Option("sc", true, "Input score file");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("rl", true, "Read length");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("step", true, "Step size");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("o", true, "Output file");
			o.setRequired(true);
			options.addOption(o);

			o = new Option("chr", true, "List of chromosomes that should be included");
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

			ArrayList<String> chrs = null;
			for (String c : line.getOptionValues("chr")) {
				if (chrs == null)
					chrs = new ArrayList<String>();
				chrs.add(c);
			}

			reduceScoreFile(new File(line.getOptionValue("sc")), Integer.parseInt(line.getOptionValue("rl")), Integer.parseInt(line.getOptionValue("step")),
					chrs, new File(line.getOptionValue("o")));

		} catch (Exception e) {
			e.printStackTrace();
			usage(options, e.toString());
		}
	}

}
