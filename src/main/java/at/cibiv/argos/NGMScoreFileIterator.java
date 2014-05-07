package at.cibiv.argos;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Reads the NextGenMap score file and iterates over its values.
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class NGMScoreFileIterator implements Iterator<ReadScores> {

	public static boolean debug = false;
	private int totalReads;
	private Map<Integer, String> chrMap = new HashMap<Integer, String>();
	private Map<String, Integer> chrMapReverse = new HashMap<String, Integer>();
	BufferedReader r = null;
	String cached = null;
	// the used read length (used to normalize the scores)
	private int rl;
	private int step;

	/**
	 * Constructor.
	 * 
	 * @param fn
	 * @throws IOException
	 */
	public NGMScoreFileIterator(File scores, int rl, int step) throws IOException {
		if (!scores.exists())
			throw new FileNotFoundException(scores.getAbsolutePath());
		this.rl = rl;
		this.step = step;

		// read header
		if (scores.getAbsolutePath().endsWith(".gz"))
			this.r = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(scores)), "UTF-8"));
		else
			this.r = new BufferedReader(new FileReader(scores));

		cached = r.readLine().substring(1).trim();
		this.totalReads = Integer.parseInt(cached);
		cached = r.readLine().substring(1).trim();
		String tmp2[] = cached.split("\t");
		for (String t : tmp2) {
			chrMap.put(Integer.parseInt(t.split(":")[0]), t.split(":")[1]);
			chrMapReverse.put(t.split(":")[1], Integer.parseInt(t.split(":")[0]));
		}
		if (debug) {
			System.out.println("total reads:\t" + totalReads);
			System.out.println("chrMap:\t" + chrMap);
		}
		cached = r.readLine();
	}

	/**
	 * Convert chromosome indices to chr strings.
	 * 
	 * @param idx
	 * @return
	 */
	public String chrIdx2Str(int idx) {
		return chrMap.get(idx);
	}

	/**
	 * Convert chr strings to chromosome indices.
	 * 
	 * @param idx
	 * @return
	 */
	public Integer str2ChrIdx(String c) {
		return chrMapReverse.get(c);
	}

	@Override
	public boolean hasNext() {
		return (cached != null);
	}

	@Override
	public ReadScores next() {
		ReadScores rc = null;
		try {
			rc = new ReadScores(cached, rl, step);
			rc.realChrIdx = str2ChrIdx(rc.realChr);
			cached = r.readLine();
		} catch (Exception e) {
			// System.err.println(e);
			e.printStackTrace();
			cached = null;
			return null;
		}
		return rc;
	}

	@Override
	public void remove() {
	}

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		NGMScoreFileIterator it = new NGMScoreFileIterator(new File("/project/oesi/genomicAmbiguity/pipeline-test/temp-test/tmpuBWPQI.bam-ngm.scores"), 100, 5);
		while (it.hasNext()) {
			ReadScores sc = it.next();
			System.out.println(sc);
		}
		System.out.println("Finished.");
	}

}
