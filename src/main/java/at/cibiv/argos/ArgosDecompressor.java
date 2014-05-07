package at.cibiv.argos;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import at.cibiv.argos.io.MixedInputStream;
import at.cibiv.ngs.tools.lds.GenomicITree;
import at.cibiv.ngs.tools.lds.GenomicInterval;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.PropertyConfiguration;
import at.cibiv.ngs.tools.util.StringUtils;

public class ArgosDecompressor {

	public static boolean debug = true;
	public static final String CMD = "decompress";
	public static final String CMD_INFO = "Decompresses ARGOS data.";

	/**
	 * Configuration
	 */
	private PropertyConfiguration config;

	/**
	 * data file
	 */
	private File dataFile;

	/**
	 * Convert from chr indices to strings
	 */
	private Map<Integer, String> chrMap = new HashMap<Integer, String>();
	private Map<String, Integer> chrMapReverse = new HashMap<String, Integer>();

	/**
	 * block index
	 */
	GenomicITree blockIndex = new GenomicITree(null);
	private long blockIndexLength;

	/**
	 * parameters read from data file
	 */
	int ADDR_GRID_W;

	/**
	 * Constructor.
	 * 
	 * @param dataFileName
	 * @throws IOException
	 */
	public ArgosDecompressor(String dataFileName) throws IOException {

		this.dataFile = new File(dataFileName);
		if (!dataFile.exists())
			throw new FileNotFoundException();
		if (debug)
			System.out.println("Read data from " + dataFileName);
		InputStream in = null;
		try {
			in = new FileInputStream(dataFileName);
			MixedInputStream mis = new MixedInputStream(in, false);

			// read byte-length of header
			this.blockIndexLength = mis.popLong();
			// read configuration
			this.config = PropertyConfiguration.fromString(mis.popString());
			if (debug) {
				System.out.println("------------------------\nConfiguration:");
				System.out.println(config);
				System.out.println("------------------------");
			}
			ADDR_GRID_W = Integer.parseInt(config.getProperty(ArgosCompressor.OPT_ADDR_GRID_W, ArgosCompressor.DEF_ADDR_GRID_W));
			// read block index
			do {
				String chr = mis.popString();
				int chrIdx = mis.popInteger();
				long min = mis.popLong();
				long max = mis.popLong();
				long blockByteStart = mis.popLong();

				chrMap.put(chrIdx, chr);
				chrMapReverse.put(chr, chrIdx);
				GenomicInterval gi = new GenomicInterval(StringUtils.prefixedChr(chr), min, max, "block");
				gi.setOriginalChrom(chr);
				gi.setAnnotation("chrIdx", chrIdx);
				gi.setAnnotation("blockByteStart", blockByteStart);
				blockIndex.insert(gi);
			} while (mis.getReadBytes() < blockIndexLength);

			blockIndex.buildTree();
			if (debug) {
				System.out.println("------------------------\nBlock index:");
				blockIndex.dump();
				System.out.println("------------------------");
			}

		} finally {
			if (in != null)
				in.close();
		}

	}

	/**
	 * Query a given genomic position.
	 * 
	 * @param q
	 * @throws IOException
	 */
	public ArgosHit query(GenomicPosition q) throws IOException {
		System.out.println("query " + q);
		List<? extends GenomicInterval> l = blockIndex.queryList(q);
		if (l == null || l.size() == 0) {
			return null;
		}
		if (l.size() != 1)
			throw new RuntimeException("Format error - overlapping blocks found");
		GenomicInterval gi = l.get(0);
		long offsetDataPortion = ((Long) gi.getAnnotation("blockByteStart")) + blockIndexLength;

		FileInputStream fin = null;
		try {
			fin = new FileInputStream(dataFile);
			long skipped = fin.skip(offsetDataPortion);
			if (skipped != offsetDataPortion)
				System.out.println("ERR");

			DataBlock db1 = DataBlock.fromDisc(fin);
			long gpos = q.get0Position() / ADDR_GRID_W;
			return new ArgosHit(db1.getTileAtGridPos(gpos));
		} finally {
			if (fin != null)
				fin.close();
		}

	}

	/**
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws IOException {
		// ArgosDecompressor gd = new
		// ArgosDecompressor("c:/data/genomicAmbiguity/tmpuBWPQI.bam-ngm.scores.gz.genAmb.data");
		ArgosDecompressor gd = new ArgosDecompressor("/project/oesi/genomicAmbiguity/pipeline-test/temp-test/tmpuBWPQI.bam-ngm.scores.gz.genAmb.data");

		for (int i = 0; i < 400; i += 20)
			System.out.println(gd.query(new GenomicPosition("2L", i)));
		// System.out.println( gd.query(new GenomicPosition("2L", 350)) );
	}

}
