package at.cibiv.argos;

import java.io.ByteArrayOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import at.cibiv.argos.io.MixedInputStream;
import at.cibiv.argos.io.MixedOutputStream;

public class DataBlock {

    /**
     * Amount of entities stored in this block.
     */
    private static final int MAX_PTR_PER_BLOCK = 10;
    protected int chrIdx;
    protected long gridStart;
    protected long gridEnd;
    protected List<DataTile> tiles = new ArrayList<DataTile>();
    protected Map<Long, DataTile> tileIndex = new HashMap<Long, DataTile>();
    int containedPointers = 0;

    public DataBlock(int chrIdx, long gridStart) {
	this.chrIdx = chrIdx;
	this.gridStart = gridStart;
    }

    public boolean isFull() {
	return containedPointers >= MAX_PTR_PER_BLOCK;
    }

    public void addTileSorted(DataTile dt) {
	this.tiles.add(dt);
	this.tileIndex.put(dt.gridPos, dt);
	gridEnd = dt.gridPos;
	containedPointers += dt.pointersSorted.size();
    }

    public DataTile getTileAtGridPos(long gpos) {
	return tileIndex.get(gpos);
    }

    @Override
    public String toString() {
	StringBuilder sb = new StringBuilder();
	sb.append("[B " + chrIdx + ":" + gridStart + "-" + chrIdx + ":" + gridEnd + " ");
	for (DataTile dt : tiles)
	    sb.append(dt);
	sb.append("]");
	return sb.toString();
    }

    public long toDisc(FileOutputStream argosDataOut) throws IOException {

	ByteArrayOutputStream idxStream = new ByteArrayOutputStream();
	ByteArrayOutputStream dataStream = new ByteArrayOutputStream();

	MixedOutputStream iout = new MixedOutputStream(idxStream, false);
	iout.push(chrIdx);
	iout.push(gridStart);
	iout.push(gridEnd);

	long writtenTotal = 0l;
	for (DataTile dt : tiles) {
	    long written = dt.toDisc(dataStream);
	    writtenTotal += written;
	    iout.push(written);
	}

	argosDataOut.write(idxStream.toByteArray());
	argosDataOut.write(dataStream.toByteArray());

	return writtenTotal + iout.getWrittenBytes();
    }

    public static DataBlock fromDisc(InputStream in) throws IOException {
	MixedInputStream min = new MixedInputStream(in, false);
	int chrIdx = min.popInteger();
	long gridStart = min.popLong();
	long gridEnd = min.popLong();

	List<Long> blockLengths = new ArrayList<Long>();
	for (int i = 0; i <= gridEnd - gridStart; i++) {
	    long blockSize = min.popLong();
	    blockLengths.add(blockSize);
	}
	DataBlock db = new DataBlock(chrIdx, gridStart);

	List<DataTile> dts = DataTile.fromDisc(in, blockLengths, chrIdx, gridStart);
	for (DataTile dt : dts) {
	    db.addTileSorted(dt);
	}

	return db;
    }

}
