package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import at.cibiv.ngs.tools.util.BinnedHistogram;

public class ReadScoresStatistics {

	public static void calcStats(File scoreFileGzip, int rl, int step, PrintStream out) throws IOException {

		NGMScoreFileIterator it = new NGMScoreFileIterator(scoreFileGzip, rl, step);
		int c = 0;
		BinnedHistogram h = new BinnedHistogram(100);

		while (it.hasNext()) {
			// get next entry in score file
			ReadScores sc = it.next();
			c++;
			h.push(sc.getEntities());

			if ( c % 100000 == 0 ) {
				System.out.print(".");
				break;
			}
			if ( c % 10000000 == 0 )
				System.out.println();
		}
		System.out.println(h);
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		args = new String[] { "/project2/oesi/genAmb/output/hg19/hg19.scores.gz" };
		int rl = 100;
		int step = 10;

		calcStats(new File(args[0]), rl, step, System.out);
	}

}
