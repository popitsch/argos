package at.cibiv.argos;

import java.util.ArrayList;
import java.util.List;

/**
 * Describes a result from querying a ARGOS data file.
 * 
 * @author niko.popitsch@univie.ac.at
 * 
 */
public class ArgosHit {

	private List<GenomicPointer> pointersSorted = new ArrayList<GenomicPointer>();

	/**
	 * Constructor.
	 * 
	 * @param dt
	 */
	public ArgosHit(DataTile dt) {
		if (dt != null)
			for (int i = 0; i < dt.pointersSorted.size(); i++) {
				GenomicPointer p = dt.pointersSorted.get(i);
				p.score = dt.weightsSorted.get(i);
				pointersSorted.add(p);
			}
	}

	public List<GenomicPointer> getPointersSorted() {
		return pointersSorted;
	}

	public void setPointersSorted(List<GenomicPointer> pointersSorted) {
		this.pointersSorted = pointersSorted;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (GenomicPointer p : pointersSorted)
			sb.append(p.toString() + "\n");
		return sb.toString();
	}

}
