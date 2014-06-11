package at.cibiv.argos;

import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 * Main class
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class Main {

	static String VERSION = "0.1";
	public static final String CMD = "argos";
	public static final String CMD_INFO = "Process argos data.";
	
	/**
	 * Usage
	 * 
	 * @param options
	 */
	private static void usage(Options options) {
		System.out.println("Welcome to argos " + VERSION);
		System.out.println();
		System.out.println("Usage:\t\tjava -jar x.jar <command> [options]:\t");
		System.out.println("\t\t" + ArgosCompressor.CMD + "\t" + ArgosCompressor.CMD_INFO);
		System.out.println("\t\t" + ArgosDecompressor.CMD + "\t" + ArgosDecompressor.CMD_INFO);
		System.exit(1);
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws Throwable
	 */
	public static void main(String[] args) throws IOException {

		// create the command line parser
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption("h", "help", false, "Print this usage information.");

		try {
			// parse the command line arguments
			CommandLine line = parser.parse(options, args, true);

			// validate that block-size has been set
			if (line.getArgs().length == 0 || line.hasOption("h")) {
				usage(options);
			}
			String[] args2 = Arrays.copyOfRange(args, 1, args.length);

			if (line.getArgs()[0].equalsIgnoreCase(ArgosCompressor.CMD)) {
				ArgosCompressor.main(args2);
			} else if (line.getArgs()[0].equalsIgnoreCase(ArgosDecompressor.CMD)) {
				ArgosDecompressor.main(args2);
			} else
				usage(options);

		} catch (ParseException exp) {
			System.out.println("Unexpected exception:" + exp.getMessage());
		}
	}

}
