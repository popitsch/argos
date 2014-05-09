package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import at.cibiv.ngs.tools.util.BinnedHistogram;

public class ReadScoresStatistics {

    public static void calcStats(File scoreFileGzip, int rl, int step, PrintStream out) throws IOException {

	NGMScoreFileIterator it = new NGMScoreFileIterator(scoreFileGzip, rl, step);
	BinnedHistogram bh = new BinnedHistogram(10);
	int c = 0;
	while (it.hasNext()) {
	    // get next entry in score file
	    ReadScores sc = it.next();
	    c++;
	    bh.push( sc.getEntities() );
	}
	System.out.println(bh);
    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

	//args = new String[] { "c:/data/genomicAmbiguity/ecK12/eck12_MG1655_ecoli-chr.scores.gz" };
	int rl = 100;
	int step = 10;

	calcStats(new File(args[0]), rl, step, System.out);
    }

}
