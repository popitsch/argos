package at.cibiv.argos;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import at.cibiv.argos.io.MixedInputStream;
import at.cibiv.argos.io.MixedOutputStream;

public class DataTile {

    protected List<GenomicPointer> pointersSorted = new ArrayList<GenomicPointer>();
    protected List<Float> weightsSorted = new ArrayList<Float>();
    int chrIdx;
    long gridPos;

    public DataTile(int chrIdx, long gridPos) {
	this.chrIdx = chrIdx;
	this.gridPos = gridPos;
    }

    public void addPointerSorted(GenomicPointer ta, Float weight) {
	pointersSorted.add(ta);
	weightsSorted.add(weight);
    }

    @Override
    public String toString() {
	return "[T " + chrIdx + ":" + gridPos + " " + Arrays.toString(pointersSorted.toArray()) + "/" + Arrays.toString(weightsSorted.toArray()) + "]";
    }

    public long toDisc(OutputStream dataStream) throws IOException {
	MixedOutputStream mout = new MixedOutputStream(dataStream, false);
	for (int i = 0; i < pointersSorted.size(); i++) {
	    mout.push(pointersSorted.get(i).gridTarget);
	    mout.push(weightsSorted.get(i));
	}
	mout.flush();
	return mout.getWrittenBytes();
    }

    public static List<DataTile> fromDisc(InputStream in, List<Long> blockLengths, int chrIdx, long gridPos) throws IOException {
	List<DataTile> ret = new ArrayList<DataTile>();
	MixedInputStream min = new MixedInputStream(in, false);
	for (Long blockLen : blockLengths) {
	    DataTile dt = new DataTile(chrIdx, gridPos);
	    long startBytes = min.getReadBytes();
	    do {
		long gridTarget = min.popLong();
		Float weight = min.popFloat();
		GenomicPointer gp = new GenomicPointer(chrIdx, gridPos, chrIdx, gridTarget, weight);
		dt.addPointerSorted(gp, weight);
	    } while (min.getReadBytes() < startBytes + blockLen);
	    ret.add(dt);
	    gridPos++;
	}
	return ret;
    }

}
