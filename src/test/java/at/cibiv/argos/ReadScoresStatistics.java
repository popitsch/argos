package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import at.cibiv.ngs.tools.util.BinnedHistogram;

/**
 * This simple helper class creates a histogram of the number of score entries
 * per read.
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class ReadScoresStatistics {

	public static void calcStats(File scoreFileGzip, int rl, int step, int binSize, PrintStream out) throws IOException {

		out.println("# ReadScoresStatistics.calcStats( scoreFileGzip=" + scoreFileGzip + ", rl=" + rl + ", step=" + step + ", binSize=" + binSize + ", out)");

		NGMScoreFileIterator it = new NGMScoreFileIterator(scoreFileGzip, rl, step);
		int c = 0;
		BinnedHistogram h = new BinnedHistogram(binSize);
		String uniqueExample = null;
		int unique = 0;
		while (it.hasNext()) {
			// get next entry in score file
			ReadScores sc = it.next();
			c++;
			h.push(sc.getEntities());
			if (sc.getEntities() == 1) {
				unique++;
				if (uniqueExample == null)
					uniqueExample = sc.toString();
			}

			if (c % 100000 == 0) {
				System.out.print(".");
			}
			if (c % 10000000 == 0)
				System.out.println();
		}
		out.println("# unique entries (i.e., only one score entry): " + unique);
		out.println("# unique example: " + uniqueExample);
		out.println(h.toCSV());
		out.close();
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		if (args.length < 5) {
			System.out.println("Usage: " + ReadScoresStatistics.class + " <ARGOS-score.gz file> <rl> <step> <bin-size> <Histogram file> ");
			System.exit(1);
		}

		// args = new String[] {
		// "/project2/oesi/genAmb/output/hg19/hg19.scores.gz" };
		// args = new String[] {
		// "c:/data/genomicAmbiguity/ecK12/eck12_MG1655_ecoli-chr.scores.gz" };
		int rl = Integer.parseInt(args[1]);
		int step = Integer.parseInt(args[2]);
		int binSize = Integer.parseInt(args[3]);
		calcStats(new File(args[0]), rl, step, binSize, new PrintStream(new File(args[4])));
	}

}
