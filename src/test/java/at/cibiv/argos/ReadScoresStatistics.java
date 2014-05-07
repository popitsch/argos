package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

public class ReadScoresStatistics {

    public static void calcStats(File scoreFileGzip, int rl, int step, PrintStream out) throws IOException {

	NGMScoreFileIterator it = new NGMScoreFileIterator(scoreFileGzip, rl, step);
	int c = 0;
	while (it.hasNext()) {
	    // get next entry in score file
	    ReadScores sc = it.next();
	    c++;
	    if (c > 100)
		break;
	    System.out.println(sc);
	}
    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

	args = new String[] { "" };
	int rl = 100;
	int step = 10;

	calcStats(new File(args[0]), rl, step, System.out);
    }

}
