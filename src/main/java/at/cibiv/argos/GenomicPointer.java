package at.cibiv.argos;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.lang.builder.HashCodeBuilder;

import at.cibiv.ngs.tools.util.GenomicPosition;
import at.cibiv.ngs.tools.util.GenomicPosition.COORD_TYPE;
import at.cibiv.ngs.tools.util.MapUtils;

public class GenomicPointer implements Comparable<GenomicPointer> {

	protected int chrIdxSource;
	protected long gridSource;
	protected int chrIdxTarget;
	protected long gridTarget;
	protected float score;

	public GenomicPointer(int chrIdxSource, long gridSource, int chrIdxTarget, long gridTarget, float score) {
		this.chrIdxSource = chrIdxSource;
		this.gridSource = gridSource;
		this.chrIdxTarget = chrIdxTarget;
		this.gridTarget = gridTarget;
		this.score = score;
	}

	public GenomicPosition getGridSsource() {
		return new GenomicPosition(chrIdxSource + "", gridSource, COORD_TYPE.ZEROBASED);
	}

	public boolean isSelfReference() {
		return (chrIdxSource == chrIdxTarget) && (gridSource == gridTarget);
	}

	/**
	 * {@inheritDoc} Two pointers are considered equal if their start and end
	 * positions match.
	 */
	@Override
	public int compareTo(GenomicPointer o) {
		if (o == null)
			return 1;

		if (gridSource < o.gridSource)
			return -1;
		else if (gridSource > o.gridSource)
			return 1;

		if (chrIdxSource < o.chrIdxSource)
			return -1;
		else if (chrIdxSource > o.chrIdxSource)
			return 1;

		if (gridTarget < o.gridTarget)
			return -1;
		else if (gridTarget > o.gridTarget)
			return 1;

		if (chrIdxTarget < o.chrIdxTarget)
			return -1;
		else if (chrIdxTarget > o.chrIdxTarget)
			return 1;

		return 0;
	}

	/**
	 * Two pointers are considered equal if their start and end positions match.
	 */
	@Override
	public boolean equals(Object o) {
		if (o == null)
			return false;
		if (!(o instanceof GenomicPointer))
			return false;
		GenomicPointer ptr = (GenomicPointer) o;
		if (chrIdxSource != ptr.chrIdxSource)
			return false;
		if (gridSource != ptr.gridSource)
			return false;
		if (chrIdxTarget != ptr.chrIdxTarget)
			return false;
		if (gridTarget != ptr.gridTarget)
			return false;
		return true;
	}

	@Override
	public int hashCode() {
		return new HashCodeBuilder(17, 31).append(chrIdxSource).append(gridSource).append(chrIdxTarget).append(gridTarget).toHashCode();
	}

	@Override
	public String toString() {
		return "[ptr " + chrIdxSource + ":" + gridSource + " -> " + chrIdxTarget + ":" + gridTarget + " (" + score + ") " + (isSelfReference() ? "*" : "")
				+ "]";
	}

	public static void main(String[] args) {
		HashMap<GenomicPointer, Integer> m = new HashMap<GenomicPointer, Integer>();
		GenomicPointer ga = new GenomicPointer(0, 0, 1, 1, 100);
		GenomicPointer ga2 = new GenomicPointer(0, 0, 1, 1, 200);
		m.put(ga, 900);
		System.out.println(m.get(ga));
		System.out.println(m.get(ga2));

		GenomicPointer t1 = new GenomicPointer(0, 0, 1, 1, 100);
		GenomicPointer t2 = new GenomicPointer(0, 0, 1, 8, 200);
		GenomicPointer t3 = new GenomicPointer(0, 0, 1, 2, 400);
		GenomicPointer t4 = new GenomicPointer(0, 0, 2, 8, 210);
		System.out.println(t1.equals(t1));
		System.out.println(t1.equals(t2));
		System.out.println(t2.equals(t1));
		System.out.println(t1.compareTo(t1));
		System.out.println(t1.compareTo(t2));
		System.out.println(t2.compareTo(t1));

		System.out.println(t1.hashCode());
		System.out.println(t2.hashCode());

		Set<GenomicPointer> s = new HashSet<GenomicPointer>();
		s.add(t1);
		s.add(t2);
		System.out.println(s);

		Map<GenomicPointer, Float> targetAddrMap = new TreeMap<GenomicPointer, Float>();
		System.out.println(targetAddrMap.put(t1, t1.score));
		System.out.println(targetAddrMap.put(t2, t2.score));
		System.out.println(targetAddrMap.put(t3, t3.score));
		System.out.println(targetAddrMap.put(t4, t4.score));

		System.out.println(targetAddrMap);
		targetAddrMap = MapUtils.sortMapByValues(targetAddrMap, true);
		System.out.println(targetAddrMap);

	}

}
