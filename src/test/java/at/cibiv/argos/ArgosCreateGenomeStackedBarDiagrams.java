package at.cibiv.argos;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import at.cibiv.codoc.CompressedCoverageIterator;
import at.cibiv.codoc.CoverageDecompressor;
import at.cibiv.codoc.utils.CodocException;
import at.cibiv.ngs.tools.util.BinnedHistogram;

public class ArgosCreateGenomeStackedBarDiagrams {

    public static void createGenomeStackedBarDiagram(List<File> scoreCodocFiles, List<String> scoreLabels, int binSize, PrintStream outAbs, PrintStream outPerc)
	    throws CodocException, IOException {
	List<BinnedHistogram> hist = new ArrayList<BinnedHistogram>();
	List<Long> histScales = new ArrayList<Long>();
	int maxbin = 100/binSize+1;
	for (File scoreCodocFile : scoreCodocFiles) {
	    CoverageDecompressor decomp = null;
	    try {
		decomp = CoverageDecompressor.loadFromFile(scoreCodocFile, null);
		CompressedCoverageIterator it = decomp.getCoverageIterator();
		BinnedHistogram bh = new BinnedHistogram(binSize);
		long countAll = 0;
		while (it.hasNext()) {
		    int cov = it.nextCoverage();
		    bh.push(cov);
		    countAll++;
		}
		hist.add(bh);
		histScales.add(countAll);
		
		System.out.println(bh);
	    } finally {
		if (decomp != null)
		    decomp.close();
	    }
	}

	// print hist
	for (String l : scoreLabels) {
	    outAbs.print(" \t" + l);
	    outPerc.print(" \t" + l);
	}
	outAbs.println();
	outPerc.println();
	for (int i = 0; i < maxbin; i++) {
	    outAbs.print(i);
	    outPerc.print(i);
	    for (int x = 0; x < hist.size(); x++) {
		BinnedHistogram bh = hist.get(x);
		long countAll = histScales.get(x);
		long countBin = bh.getCount(i)!=null?bh.getCount(i):0;
		outAbs.print("\t" + countBin);
		outPerc.print("\t" + ((double) countBin / (double) countAll));
	    }
	    outAbs.println();
	    outPerc.println();
	}

	System.out.println("Finished.");
    }

    /**
     * @param args
     * @throws IOException
     * @throws CodocException
     */
    public static void main(String[] args) throws CodocException, IOException {
	List<File> scoreCodocFiles = new ArrayList<File>();
	List<String> scoreLabels = new ArrayList<String>();
	scoreCodocFiles.add(new File("/project2/oesi/genAmb/output/ecK12/eck12_MG1655_ecoli-chr-GENOME.ISS.wig.codoc"));
	scoreLabels.add("e.coli");
//	scoreCodocFiles.add(new File("/project2/oesi/genAmb/output/dros_partial/droso_r5_partial-GENOME.ISS.wig.codoc"));
//	scoreLabels.add("d.mel");
//	scoreCodocFiles.add(new File("/project2/oesi/genAmb/output/hg19/hg19-GENOME.ISS.wig.codoc"));
//	scoreLabels.add("h.sapiens");

	File outAbs = new File("/project2/oesi/genAmb/output/POSTER_DATA/stackedHistAbs.csv");
	File outPrec = new File("/project2/oesi/genAmb/output/POSTER_DATA/stackedHistPrec.csv");
	createGenomeStackedBarDiagram(scoreCodocFiles, scoreLabels, 1, new PrintStream(outAbs), new PrintStream(outPrec));
    }

}
