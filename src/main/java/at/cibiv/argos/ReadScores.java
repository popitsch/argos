package at.cibiv.argos;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.GenomicPosition.COORD_TYPE;

/**
 * Represents a read in the NextGenMap read-scores file.
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class ReadScores {

	protected String realChr;
	protected int realChrIdx = -1; // has to be set manually!
	protected long realPosition0;
	protected float[] scoresSorted;
	protected boolean[] isReverseComplementSorted;
	protected int[] chrIdxSorted;
	protected long[] posSorted;
	protected String readName;

	Map<Integer, Float> issCache = new HashMap<Integer, Float>();
	Map<Integer, Integer> ambCache = new HashMap<Integer, Integer>();
	Map<Integer, Float> msdCache = new HashMap<Integer, Float>();

	/**
	 * The maximum achievable alignment score
	 */
	private Float maxScore = null;

	/**
	 * Cutoff used to determine whether a regions is considered "similar"
	 */
	private Float simCutoff;
	
	/**
	 * Used step size.
	 */
	private int step;

	/**
	 * Constructor.
	 * 
	 * @param t
	 * @throws ParseException
	 */
	public ReadScores(String t, int rl, int step) throws ParseException {
		this.step = step;
		String[] tags = t.split("\t");
		// parse the read name and the originating genome position
		this.readName = tags[0].substring(1, tags[0].indexOf("_"));
		String tmp = tags[0].substring(tags[0].indexOf("_") + 1);
		GenomicPosition pos = GenomicPosition.fromString(tmp, COORD_TYPE.ONEBASED);
		this.realChr = pos.getChromosomeOriginal();
		this.realPosition0 = pos.get0Position();

		this.scoresSorted = new float[tags.length - 1];
		this.isReverseComplementSorted = new boolean[tags.length - 1];
		this.chrIdxSorted = new int[tags.length - 1];
		this.posSorted = new long[tags.length - 1];

		// scores will be divided by this
		this.maxScore = rl * 10f;
		this.simCutoff = 0.5f;

		if (tags.length < 2)
			throw new ParseException("Error parsing line " + Arrays.toString(tags), 0);

		String[] e = tags[1].split(":");
		this.chrIdxSorted[0] = Integer.parseInt(e[0]);
		this.posSorted[0] = Long.parseLong(e[1]);
		this.isReverseComplementSorted[0] = (e[2].equals("1"));
		this.scoresSorted[0] = Float.parseFloat(e[3]) / maxScore;

		float lastscore = scoresSorted[0];
		for (int i = 2; i < tags.length; i++) {
			e = tags[i].split(":");
			if (e.length != 4)
				throw new ParseException("Error parsing tag " + Arrays.toString(e), i);
			this.chrIdxSorted[i - 1] = Integer.parseInt(e[0]);
			this.posSorted[i - 1] = Long.parseLong(e[1]);
			this.isReverseComplementSorted[i - 1] = (e[2].equals("1"));
			this.scoresSorted[i - 1] = lastscore - (Float.parseFloat(e[3]) / maxScore);
			lastscore = this.scoresSorted[i - 1];
		}
	}

	/**
	 * @return the best alignment score
	 */
	public float getBestScore() {
		return scoresSorted[0];
	}

	/**
	 * Check whether the score with the given idx is within the context window.
	 * 
	 * @param ctxSize
	 * @param idx
	 * @return
	 */
	private boolean isInContext(int idx, Integer ctxSize) {

		if (ctxSize == null)
			return true;

		if (realChrIdx != chrIdxSorted[idx])
			return false; // skip as on other chrom
		int diff = (int) Math.abs(realPosition0 - posSorted[idx]);
		if (diff > ctxSize)
			return false; // skip as too far away!
		return true;
	}

	/**
	 * @return the second best alignment score (if any). Skips the current position.
	 */
	public Float getSecondBestScore(Integer ctxSize) {
		for (int i = 0; i < scoresSorted.length; i++) {
		
			// skip out-of context scores
			if (!isInContext(i, ctxSize))
				continue;

			// skip original score FIXME: check this!
			if ((realChrIdx == chrIdxSorted[i]) && (realPosition0-2*step<=posSorted[i] ) && (realPosition0+2*step>=posSorted[i])) {
				continue;
			}
			
//			if ( ctxSize == null )
//				//System.out.println("["+ctxSize+"] " + realChrIdx+":"+(realPosition0-realPosition0%step) + "  vs  " + chrIdxSorted[i]+":"+(posSorted[i]-posSorted[i]%step));
//				System.out.println(scoresSorted[i] + " -> " + this);
			
			return scoresSorted[i];
		}

//		if ( ctxSize == null )
//			System.out.println("NULL -> " + this);

		return null;
	}

	// /**
	// * The average over all scores. If similarOnly is set, then only scores
	// * above the simCutoff are taken into account.
	// *
	// * @return the average score
	// */
	// public Float getAvgOtherScores(Integer ctxSize, boolean similarOnly) {
	// if (scoresSorted.length < 2)
	// return null;
	// float count = 0f, sum = 0f;
	// for (int i = 1; i < scoresSorted.length; i++) {
	// // stop if scores are below threshold.
	// if (similarOnly)
	// if (scoresSorted[i] < simCutoff)
	// break;
	// // skip out-of context scores
	// if (!isInContext(i, ctxSize))
	// continue;
	// count++;
	// sum += scoresSorted[i];
	// }
	// if (count == 0f)
	// return null;
	// return (sum / count);
	// }

	/**
	 * This score informs you about the (maximum) similarity of the current
	 * region to any region in the genome/the provided genomic context, i.e. it
	 * is the subtraction of the maximum score and the second-best score (in the
	 * whole genome or within a genomic window +/- ctxSize).
	 * 
	 * @param ctxSize
	 * @return
	 */
	public float getISS(Integer ctxSize) {
		Float iss = issCache.get(ctxSize);
		if (iss != null) {
			return iss;
		}
		Float sb = getSecondBestScore(ctxSize);
		if (sb == null)
			sb = 0f; // no 2nd hit -> ISS should become max.
		iss = getBestScore() - sb;
		issCache.put(ctxSize, iss);
		
//		if ( ctxSize == null )System.out.println("{ } iss " + iss);	
//		try {
//			System.in.read();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		return iss;
	}

	/**
	 * This score informs you about how the similarities of similar regions are
	 * distributed by providing the median of the scores of all regions that are
	 * similar (50% of the maximum alignment score)
	 * 
	 * @param ctxSize
	 * @return
	 */
	public float getMSD(Integer ctxSize) {
		Float msd = msdCache.get(ctxSize);
		if (msd != null)
			return msd;

		List<Float> relevantScores = new ArrayList<Float>();

		for (int i = 0; i < scoresSorted.length; i++) {
			// stop if scores are below threshold.
			if (scoresSorted[i] < simCutoff)
				break;
			// skip out-of context scores
			if (!isInContext(i, ctxSize))
				continue;
			relevantScores.add(scoresSorted[i]);
		}

		if (relevantScores.size() == 0)
			msd = 0f;
		else {
			// calc median (use mean between middle scores if required)
			int mid = relevantScores.size() / 2;
			if (relevantScores.size() % 2 == 1) {
				// System.out.println("MED " + mid + " ? " +
				// relevantScores.get(mid)
				// );
				msd = getBestScore() - relevantScores.get(mid);
			} else {
				msd = getBestScore() - (relevantScores.get(mid - 1) + relevantScores.get(mid)) / 2.0f;
			}
		}
		msdCache.put(ctxSize, msd);
		return msd;
	}

	/**
	 * This score informs you about the number of "similar" regions in the
	 * genome, i.e. it is the count of regions that have a minimum similarity
	 * (usually some percentage of the maximum alignment score)
	 * 
	 * @param ctxSize
	 * @return
	 */
	public float getAMB(Integer ctxSize) {
		Integer amb = ambCache.get(ctxSize);
		if (amb != null)
			return amb;

		amb = 0;
		for (int i = 0; i < scoresSorted.length; i++) {
			// stop if scores are below threshold.
			if (scoresSorted[i] < simCutoff)
				break;
			// skip out-of context scores (also skips the current score!)
			if (!isInContext(i, ctxSize))
				continue;
			amb++;
		}

		ambCache.put(ctxSize, amb);
		return amb;
	}

	/**
	 * @return the genomic position where the read stems from that was used to
	 *         calculate the alignment scores..
	 */
	public GenomicPosition getRealPosition() {
		return new GenomicPosition(realChr, realPosition0, COORD_TYPE.ZEROBASED);
	}

	/**
	 * @return the number of stored entities (position/score pairs)
	 */
	public int getEntities() {
		return scoresSorted.length;
	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(getRealPosition().toString1basedOrig() + "\t");
		for (int i = 0; i < scoresSorted.length; i++)
			sb.append(chrIdxSorted[i] + ":" + posSorted[i] + " " + scoresSorted[i] + (isReverseComplementSorted[i] ? " (T)" : " (F)") + "\t");
		return sb.toString();
	}

}
