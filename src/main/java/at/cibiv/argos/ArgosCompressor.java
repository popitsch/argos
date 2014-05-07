package at.cibiv.argos;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;

import at.cibiv.argos.io.MixedOutputStream;
import at.cibiv.ngs.tools.lds.GenomicITree;
import at.cibiv.ngs.tools.lds.GenomicInterval;
import at.cibiv.ngs.tools.util.FileUtils;
import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.MapUtils;
import at.cibiv.ngs.tools.util.PropertyConfiguration;
import at.cibiv.ngs.tools.wig.WigOutputStream;

/**
 * Compresses ARGOS data.
 * 
 * TODO: add GZIP block compression?
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class ArgosCompressor {

    public static boolean debug = true;
    public static final String CMD = "compress";
    public static final String CMD_INFO = "Compresses ARGOS data.";
    public static final SimpleDateFormat sdf = new SimpleDateFormat("yyyy.MM.dd G 'at' HH:mm:ss z");

    /**
     * Current configuration.
     */
    private PropertyConfiguration config;

    /**
     * configuration parameters
     */
    public static final String OPT_SCORE_FILE = "scores";
    public static final String OPT_READ_LENGTH = "readLength";
    public static final String OPT_STEP_SIZE = "stepSize";
    public static final String OPT_ADDR_GRID_W = "addressGridWidth";
    public static final String OPT_TMP_DIR = "tmpDir";
    public static final String OPT_KEEP_WORKDIR = "keepWorkDir";
    public static final String OPT_VERBOSE = "v";
    public static final String OPT_OUT_DATA = "outData";
    public static final String OPT_OUT_WIG_AMB = "outWigAMB";
    public static final String OPT_OUT_WIG_ISS = "outWigISS";
    public static final String OPT_OUT_WIG_MSD = "outWigMSD";
    public static final String OPT_CTX_SIZE = "ctxSize";

    /**
     * Default values
     */
    public static final String DEF_READ_LENGTH = "100";
    public static final String DEF_STEP_SIZE = "5";
    public static final String DEF_ADDR_GRID_W = "20";
    public static final String DEF_AMBCTX_SIZE = "1000";

    /**
     * Default configuration
     * 
     * @return
     * @throws IOException
     */
    public static PropertyConfiguration getDefaultConfiguration() {
	StringBuffer sb = new StringBuffer();

	// set default compression algorithm
	sb.append("argos-version=" + Main.VERSION + "\n");

	// set creation time
	sb.append("created=" + sdf.format(new Date()) + "\n");

	return PropertyConfiguration.fromString(sb.toString());
    }

    public PropertyConfiguration getConfig() {
	return config;
    }

    /**
     * Constructor.
     * 
     * @throws IOException
     */
    public ArgosCompressor(PropertyConfiguration config) throws IOException {
	this.config = config;

	if (!config.hasProperty(OPT_SCORE_FILE))
	    throw new FileNotFoundException("No valid score file passed");
	// parse paramteres
	File scoreFile = new File(config.getProperty(OPT_SCORE_FILE));
	int rl = Integer.parseInt(config.getProperty(OPT_READ_LENGTH, DEF_READ_LENGTH));
	int step = Integer.parseInt(config.getProperty(OPT_STEP_SIZE, DEF_STEP_SIZE));
	int addressWith = Integer.parseInt(config.getProperty(OPT_ADDR_GRID_W, DEF_ADDR_GRID_W));
	File workDir = config.hasProperty(OPT_TMP_DIR) ? new File(config.getProperty(OPT_TMP_DIR)) : new File(scoreFile.getParentFile(), "" + Math.random());
	boolean keepWorkDir = config.getBooleanProperty(OPT_KEEP_WORKDIR, false);
	ArgosCompressor.debug = config.getBooleanProperty(OPT_VERBOSE, false);

	File outData = new File(config.getProperty(OPT_OUT_DATA));
	File outWigISS = new File(config.getProperty(OPT_OUT_WIG_ISS));
	File outWigAMB = new File(config.getProperty(OPT_OUT_WIG_AMB));
	File outWigMSD = new File(config.getProperty(OPT_OUT_WIG_MSD));
	int ctxSize = Integer.parseInt(config.getProperty(OPT_CTX_SIZE, "" + DEF_AMBCTX_SIZE));

	compress(scoreFile, rl, step, addressWith, outData, outWigISS, outWigAMB, outWigMSD, ctxSize, workDir, keepWorkDir);
    }

    /**
     * Write block to disk
     * 
     * @param currentBlock
     * @param argosDataOut
     * @param it
     * @param ADDR_GRID_W
     * @param blockIndex
     * @param blockByteStart
     * @return
     * @throws IOException
     */
    private long writeBlock(DataBlock currentBlock, FileOutputStream argosDataOut, NGMScoreFileIterator it, int ADDR_GRID_W, GenomicITree blockIndex,
	    long blockByteStart) throws IOException {

	// write block to disc
	long writtenBytes = currentBlock.toDisc(argosDataOut);

	// add to block index
	GenomicInterval gi = new GenomicInterval(it.chrIdx2Str(currentBlock.chrIdx), currentBlock.gridStart * ADDR_GRID_W, (currentBlock.gridEnd + 1)
		* ADDR_GRID_W - 1, "block");
	gi.setOriginalChrom(it.chrIdx2Str(currentBlock.chrIdx));
	gi.setAnnotation("chrIdx", currentBlock.chrIdx);
	gi.setAnnotation("blockByteStart", blockByteStart);
	blockIndex.insert(gi);
	return writtenBytes;
    }

    /**
     * Constructor.
     * 
     * @param scoreFile
     * @param rl
     * @param step
     * @param tempDir
     * @throws IOException
     */
    public void compress(File scoreFile, int rl, int step, int ADDR_GRID_W, File outData, File outWigISS, File outWigAMB, File outWigMSD, int ctxSize,
	    File tempDir, boolean keepWorkDir) throws IOException {

	boolean writeTiles = false; // FIXME later

	if (!tempDir.exists()) {
	    if (!tempDir.mkdirs())
		throw new IOException("Could not create temporary directory " + tempDir);
	}

	File outWigISS_loc = new File(outWigISS.getAbsolutePath() + ".ctx" + ctxSize + ".wig");
	File outWigAMB_loc = new File(outWigAMB.getAbsolutePath() + ".ctx" + ctxSize + ".wig");
	File outWigMSD_loc = new File(outWigMSD.getAbsolutePath() + ".ctx" + ctxSize + ".wig");

	// iterates over score entries
	NGMScoreFileIterator it = new NGMScoreFileIterator(scoreFile, rl, step);
	// Stores the scores per genomic position
	Map<GenomicPosition, float[]> scoreMap = new TreeMap<GenomicPosition, float[]>();
	// output streams
	WigOutputStream outISS = null;
	WigOutputStream outAMB = null;
	WigOutputStream outMSD = null;
	WigOutputStream outISS_loc = null;
	WigOutputStream outAMB_loc = null;
	WigOutputStream outMSD_loc = null;

	File argosData = null;
	FileOutputStream argosDataOut = null;
	File argosIdx = null;
	FileOutputStream argosIdxOut = null;

	// Tree-index
	GenomicITree blockIndex = new GenomicITree(null);

	// will be set to true if there were errors.
	boolean hadErrors = false;
	boolean readWasSkipped = false;
	int skippedReads = 0;
	int count = 0;

	try {
	    outISS = new WigOutputStream(outWigISS);
	    outAMB = new WigOutputStream(outWigAMB);
	    outMSD = new WigOutputStream(outWigMSD);
	    outISS_loc = new WigOutputStream(outWigISS_loc);
	    outAMB_loc = new WigOutputStream(outWigAMB_loc);
	    outMSD_loc = new WigOutputStream(outWigMSD_loc);

	    argosData = FileUtils.createTempFile(tempDir, "data");
	    argosIdx = FileUtils.createTempFile(tempDir, "idx");

	    argosDataOut = new FileOutputStream(argosData);
	    argosIdxOut = new FileOutputStream(argosIdx);

	    // current position.
	    int chrIdx = -1;
	    long pos0 = 0;

	    // map of target_addresses -> score
	    Map<GenomicPointer, Float> targetAddrMap = new TreeMap<GenomicPointer, Float>();

	    // index of last block
	    DataBlock currentBlock = null;
	    Long lastPos0Grid = null;

	    // increasing counter of written bytes to data portion
	    long blockByteStart = 0l;

	    while (it.hasNext()) {
		// get next entry in score file
		ReadScores sc = it.next();

		if (++count % 1000000 == 0)
		    System.out.print(".");
		if (++count % 100000000 == 0)
		    System.out.println();

		// check for chrom change.
		boolean chromWasChanged = (chrIdx != sc.realChrIdx);
		if (chromWasChanged) {
		    if (debug)
			System.out.println("CHROM CHANGE");
		    // FIXME: check border-effects!
		    scoreMap = new TreeMap<GenomicPosition, float[]>();
		}

		// check whether a read was skipped (because it could not be
		// mapped)
		readWasSkipped = (!chromWasChanged && (pos0 + step != sc.realPosition0));
		if (readWasSkipped) {
		    // FIXME readlength and stepsize should be read from the
		    // scores archive.
		    if (debug)
			System.out.println("Read was skipped - please check whether you provided the proper read-length/step parameters!! " + sc.realChr + ":"
				+ (pos0 + step));
		    skippedReads++;
		}

		// update current positions
		pos0 = sc.realPosition0;
		chrIdx = sc.realChrIdx;
		long pos0Grid = (long) (pos0 / ADDR_GRID_W);
		if (lastPos0Grid == null)
		    lastPos0Grid = pos0Grid;

		// update score maps: add score of read to region spanned by
		// (original) read
		for (long pos1 = pos0 + 1; pos1 < pos0 + rl + 1; pos1 += step) {
		    GenomicPosition gp = new GenomicPosition(sc.realChr, pos1);
		    float[] scores = scoreMap.get(gp);
		    if (scores == null)
			scores = new float[7];
		    scores[0]++;
		    scores[1] += sc.getISS(null);
		    scores[2] += sc.getAMB(null);
		    scores[3] += sc.getMSD(null);
		    scores[4] += sc.getISS(ctxSize);
		    scores[5] += sc.getAMB(ctxSize);
		    scores[6] += sc.getMSD(ctxSize);
		    scoreMap.put(gp, scores);
		}

		// finish scores
		SortedSet<GenomicPosition> keys = new TreeSet<GenomicPosition>(scoreMap.keySet());
		List<GenomicPosition> toFinish = new ArrayList<GenomicPosition>();
		for (GenomicPosition p : keys)
		    if (p.compareTo(sc.getRealPosition()) <= 0) {
			toFinish.add(p);
		    } else
			break;
		for (GenomicPosition p : toFinish) {
		    float[] scores = scoreMap.get(p);

		    outISS.push(p, scores[1] / scores[0], step);
		    outAMB.push(p, scores[2] / scores[0], step);
		    outMSD.push(p, scores[3] / scores[0], step);
		    outISS_loc.push(p, scores[4] / scores[0], step);
		    outAMB_loc.push(p, scores[5] / scores[0], step);
		    outMSD_loc.push(p, scores[6] / scores[0], step);
		    scoreMap.remove(p);
		}

		// finish tiles
		if (writeTiles) {
		    if (pos0Grid != lastPos0Grid) {
			// sort by score. FIXME: test this and
			// measureperformance!
			targetAddrMap = MapUtils.sortMapByValues(targetAddrMap, true);
			if (debug) {
			    System.out.println("------------ TILE " + lastPos0Grid + "--------------------");
			    for (GenomicPointer ta : targetAddrMap.keySet()) {
				System.out.println(ta + "=> cum.score: " + targetAddrMap.get(ta) + (ta.isSelfReference() ? "*" : ""));
			    }
			    System.out.println("--------------------------------------");
			}
			// push to block
			if (currentBlock == null) {
			    // first block
			    currentBlock = new DataBlock(chrIdx, lastPos0Grid);
			}

			DataTile dt = new DataTile(chrIdx, lastPos0Grid);
			for (GenomicPointer ta : targetAddrMap.keySet()) {
			    dt.addPointerSorted(ta, targetAddrMap.get(ta));
			}
			currentBlock.addTileSorted(dt);

			targetAddrMap = new TreeMap<GenomicPointer, Float>();
			lastPos0Grid = pos0Grid;

			if (currentBlock.chrIdx != chrIdx || currentBlock.isFull()) {
			    // finish block
			    if (debug)
				System.out.println("FINISH BLOCK");
			    if (debug)
				System.out.println(currentBlock);

			    blockByteStart += writeBlock(currentBlock, argosDataOut, it, ADDR_GRID_W, blockIndex, blockByteStart);

			    currentBlock = new DataBlock(chrIdx, pos0Grid);
			}

		    }

		    // update target address map
		    for (int i = 0; i < sc.chrIdxSorted.length; i++) {
			GenomicPointer ptr = new GenomicPointer(chrIdx, pos0Grid, sc.chrIdxSorted[i], (long) (sc.posSorted[i] / ADDR_GRID_W),
				sc.scoresSorted[i]);
			if (targetAddrMap.containsKey(ptr)) {
			    // cumulate score!
			    ptr.score += targetAddrMap.get(ptr);
			}
			// System.out.println(ta + "/" + targetScore);
			targetAddrMap.put(ptr, ptr.score);
		    }
		}

	    } // while

	    // finish final positions
	    SortedSet<GenomicPosition> keys = new TreeSet<GenomicPosition>(scoreMap.keySet());
	    for (GenomicPosition p : keys) {
		float[] scores = scoreMap.get(p);
		outISS.push(p, scores[1] / scores[0], step);
		outAMB.push(p, scores[2] / scores[0], step);
		outMSD.push(p, scores[3] / scores[0], step);
		outISS_loc.push(p, scores[4] / scores[0], step);
		outAMB_loc.push(p, scores[5] / scores[0], step);
		outMSD_loc.push(p, scores[6] / scores[0], step);
	    }

	    // finish tiles & blocks
	    targetAddrMap = MapUtils.sortMapByValues(targetAddrMap, true);
	    if (chrIdx < 0)
		throw new IOException("Chridx unset - maybe scores file is empty?");
	    DataTile dt = new DataTile(chrIdx, lastPos0Grid);
	    for (GenomicPointer ta : targetAddrMap.keySet()) {
		dt.addPointerSorted(ta, targetAddrMap.get(ta));
	    }
	    if (currentBlock != null) {
		currentBlock.addTileSorted(dt);
		blockByteStart += writeBlock(currentBlock, argosDataOut, it, ADDR_GRID_W, blockIndex, blockByteStart);
	    }

	    // serialize block Index.
	    blockIndex.buildTree();
	    if (debug)
		blockIndex.dump();
	    ByteArrayOutputStream blockIndexBA = new ByteArrayOutputStream();
	    MixedOutputStream mout = new MixedOutputStream(blockIndexBA, false);
	    // write configuration to buffer
	    PropertyConfiguration decompressionConfig = new PropertyConfiguration();
	    decompressionConfig.setProperty(OPT_READ_LENGTH, rl);
	    decompressionConfig.setProperty(OPT_STEP_SIZE, step);
	    decompressionConfig.setProperty(OPT_ADDR_GRID_W, ADDR_GRID_W);
	    mout.push(decompressionConfig.toString());
	    // write block index to buffer
	    for (GenomicInterval gi : blockIndex.getIntervalsSorted()) {
		mout.push(gi.getOriginalChrom());
		mout.push((Integer) gi.getAnnotation("chrIdx"));
		mout.push(gi.getMin());
		mout.push(gi.getMax());
		mout.push((Long) gi.getAnnotation("blockByteStart"));
	    }
	    mout.flush();
	    // write size of header (add 8 bytes for the length itself!)
	    new MixedOutputStream(argosIdxOut, false).push(mout.getWrittenBytes() + 8);
	    // write buffer to header
	    argosIdxOut.write(blockIndexBA.toByteArray());

	    // close temp output streams
	    argosDataOut.flush();
	    argosDataOut.close();
	    argosDataOut = null;
	    argosIdxOut.flush();
	    argosIdxOut.close();
	    argosIdxOut = null;

	    // concatenate header and data blocks
	    FileOutputStream fout = new FileOutputStream(outData);
	    FileUtils.concatenateFile(fout, false, argosIdx, argosData);
	    fout.close();

	    System.out.println("Finished. Skipped reads: " + skippedReads);

	} catch (Exception e) {
	    e.printStackTrace();
	    hadErrors = true;

	} finally {

	    if (outISS != null)
		outISS.close();
	    if (outAMB != null)
		outAMB.close();
	    if (outMSD != null)
		outMSD.close();
	    if (outISS_loc != null)
		outISS_loc.close();
	    if (outAMB_loc != null)
		outAMB_loc.close();
	    if (outMSD_loc != null)
		outMSD_loc.close();
	    if (argosDataOut != null)
		argosDataOut.close();
	    if (argosData != null)
		argosData.delete();
	    if (argosIdxOut != null)
		argosIdxOut.close();
	    if (argosIdx != null)
		argosIdx.delete();

	    if (hadErrors) {
		if (debug)
		    System.err.println("ERRORS - removing files.");

		// remove (partially created) result files!
		if ((outData != null) && (outData.exists()))
		    outData.delete();

		if ((outWigISS != null) && (outWigISS.exists()))
		    outWigISS.delete();
		if ((outWigAMB != null) && (outWigAMB.exists()))
		    outWigAMB.delete();
		if ((outWigMSD != null) && (outWigMSD.exists()))
		    outWigMSD.delete();
		if ((outWigISS_loc != null) && (outWigISS_loc.exists()))
		    outWigISS_loc.delete();
		if ((outWigAMB_loc != null) && (outWigAMB_loc.exists()))
		    outWigAMB_loc.delete();
		if ((outWigMSD_loc != null) && (outWigMSD_loc.exists()))
		    outWigMSD_loc.delete();
	    }

	    // delete temporary files.
	    if (!keepWorkDir) {
		// delete temporary directory
		if (!tempDir.delete())
		    System.err.println("Could not remove temporary directory " + tempDir);
	    }

	}
	if (debug)
	    System.out.println("Finished.");

    }

    /**
     * Print usage information.
     * 
     * @param options
     */
    @SuppressWarnings({ "unchecked" })
    private static void usage(Options options, String subcommand, String e) {
	int w = 120;

	Options mandatory = new Options();

	for (Object o : options.getRequiredOptions())
	    if (o instanceof Option) {
		if (((Option) o).isRequired())
		    mandatory.addOption((Option) o);
	    } else if (o instanceof OptionGroup) {
		if (((OptionGroup) o).isRequired())
		    mandatory.addOptionGroup((OptionGroup) o);
	    }
	for (Option o : (Collection<Option>) options.getOptions())
	    if (o.isRequired())
		mandatory.addOption(o);
	PrintWriter pw = new PrintWriter(System.out);
	HelpFormatter hf = new HelpFormatter();
	hf.setSyntaxPrefix("");

	hf.printUsage(pw, w, "Usage:  java -jar x.jar " + CMD, options);
	pw.println();
	pw.flush();

	hf.printHelp(pw, w, "Mandatory Params", "", mandatory, 8, 2, "");
	pw.flush();
	pw.println();

	hf.printHelp(pw, w, "Params", "", options, 8, 2, "");

	pw.flush();

	System.out.println();
	System.out.println("\tArgosCompressor " + (Main.VERSION == null ? "" : Main.VERSION) + " (c) 2014");

	if (e != null)
	    System.out.println("\nError: " + e);
	System.exit(1);
    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

	// args = new String[] { "compress",
	// "-s",
	// "/project/oesi/genAmb/output/B31/B31_Schutzer_reference-goodnames.scores.gz",
	// "-rl", "50", "-step", "2"
	// //
	// "/project/oesi/genomicAmbiguity/pipeline-test/temp-test/tmpuBWPQI.bam-ngm.scores.gz"
	//
	// };

	// args = new String[] { "compress",
	// "-s",
	// "/project/oesi/genomicAmbiguity/pipeline-test/temp-test/tmpr3eOVT.bam-ngm.scores"
	//
	// };

	// new ArgosCompressor(new
	// File("/project/oesi/genomicAmbiguity/pipeline-test/temp-test/tmpuBWPQI.bam-ngm.scores"),
	// 100, 5, new File(
	// "/project/oesi/genomicAmbiguity/pipeline-test/temp-test/"));
	// new ArgosCompressor(new
	// File("c:/data/genomicAmbiguity/tmpuBWPQI.bam-ngm.scores.gz"), 100, 5,
	// new File(
	// "c:/data/genomicAmbiguity/"));

	CommandLineParser parser = new PosixParser();

	// create the Options
	Options options = new Options();
	options.addOption("h", "help", false, "Print this usage information.");

	try {
	    // parse the command line arguments
	    CommandLine line = parser.parse(options, args, true);

	    options = new Options();

	    Option opt = new Option("s", true, "Input score file.");
	    opt.setLongOpt(OPT_SCORE_FILE);
	    opt.setRequired(true);
	    options.addOption(opt);

	    opt = new Option("rl", true, "Used read length (default:" + DEF_READ_LENGTH + ")");
	    opt.setLongOpt(OPT_READ_LENGTH);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("step", true, "Used step size (default:" + DEF_STEP_SIZE + ")");
	    opt.setLongOpt(OPT_STEP_SIZE);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("gw", true, "Used address grid width (default:" + DEF_ADDR_GRID_W + ")");
	    opt.setLongOpt(OPT_ADDR_GRID_W);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("od", true, "Output file for the compressed data (default: <scorefile>+\".argos.data\")");
	    opt.setLongOpt(OPT_OUT_DATA);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("oi", true, "Output file for the ISS signal (default: <scorefile>+\".argos.ISS.wig\")");
	    opt.setLongOpt(OPT_OUT_WIG_ISS);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("oa", true, "Output file for the AMB signal (default: <scorefile>+\".argos.AMB.wig\")");
	    opt.setLongOpt(OPT_OUT_WIG_AMB);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("om", true, "Output file for the MSD signal (default: <scorefile>+\".argos.MSD.wig\")");
	    opt.setLongOpt(OPT_OUT_WIG_MSD);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("cx", true, "Size of the genomic window used for calculating the context-signals (default:" + DEF_AMBCTX_SIZE + ")");
	    opt.setLongOpt(OPT_CTX_SIZE);
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option("t", true, "Temporary working directory (default is current dir)");
	    opt.setLongOpt("temp");
	    opt.setRequired(false);
	    options.addOption(opt);

	    opt = new Option(OPT_KEEP_WORKDIR, false, "Do not delete the work directory (e.g., for debugging purposes)");
	    opt.setRequired(false);
	    options.addOption(opt);

	    options.addOption(OPT_VERBOSE, "verbose", false, "be verbose");

	    if (line.hasOption("h")) {
		usage(options, null, null);
	    }

	    line = parser.parse(options, args);

	    // parse configuration
	    PropertyConfiguration conf = new PropertyConfiguration();
	    conf.setProperty(OPT_SCORE_FILE, line.getOptionValue(OPT_SCORE_FILE));
	    if (line.hasOption(OPT_READ_LENGTH))
		conf.setProperty(OPT_READ_LENGTH, line.getOptionValue(OPT_READ_LENGTH));
	    if (line.hasOption(OPT_STEP_SIZE))
		conf.setProperty(OPT_STEP_SIZE, line.getOptionValue(OPT_STEP_SIZE));
	    if (line.hasOption(OPT_ADDR_GRID_W))
		conf.setProperty(OPT_ADDR_GRID_W, line.getOptionValue(OPT_ADDR_GRID_W));
	    if (line.hasOption(OPT_TMP_DIR))
		conf.setProperty(OPT_TMP_DIR, line.getOptionValue(OPT_TMP_DIR));
	    if (line.hasOption(OPT_CTX_SIZE))
		conf.setProperty(OPT_CTX_SIZE, line.getOptionValue(OPT_CTX_SIZE));
	    conf.setProperty(OPT_OUT_DATA, line.getOptionValue(OPT_OUT_DATA, line.getOptionValue(OPT_SCORE_FILE) + ".argos.data"));
	    conf.setProperty(OPT_OUT_WIG_ISS, line.getOptionValue(OPT_OUT_WIG_ISS, line.getOptionValue(OPT_SCORE_FILE) + ".argos.ISS.wig"));
	    conf.setProperty(OPT_OUT_WIG_AMB, line.getOptionValue(OPT_OUT_WIG_AMB, line.getOptionValue(OPT_SCORE_FILE) + ".argos.AMB.wig"));
	    conf.setProperty(OPT_OUT_WIG_MSD, line.getOptionValue(OPT_OUT_WIG_MSD, line.getOptionValue(OPT_SCORE_FILE) + ".argos.MSD.wig"));
	    conf.setProperty(OPT_KEEP_WORKDIR, line.hasOption(OPT_KEEP_WORKDIR) + "");

	    conf.setProperty(OPT_VERBOSE, line.hasOption(OPT_VERBOSE) + "");

	    new ArgosCompressor(conf);

	    System.out.println("Finished.");
	    System.exit(0);

	} catch (MissingOptionException e) {
	    System.err.println("Error: " + e.getMessage());
	    System.err.println();
	    usage(options, null, e.toString());
	} catch (Throwable e) {
	    e.printStackTrace();
	    System.err.println();
	    usage(options, null, e.toString());
	}

    }

}
