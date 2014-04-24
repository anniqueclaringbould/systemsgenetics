package org.molgenis.vcf;

import java.util.ArrayList;
import java.util.List;

public class VcfSample
{
	private static final String FIELD_GT = "GT";
	private static final char GENOTYPE_UNPHASED = '/';
	private static final char GENOTYPE_PHASED = '|';
	
	private final VcfRecord vcfRecord;
	private String[] tokens;

	private transient List<String> cachedAlleles;
	
	public VcfSample(VcfRecord vcfRecord) {
		this(vcfRecord, null);
	}
	
	public VcfSample(VcfRecord vcfRecord, String[] tokens) {
		if(vcfRecord == null) throw new IllegalArgumentException("vcfRecord is null");
		this.vcfRecord = vcfRecord;
		this.tokens = tokens;
	}
	
	public String getData(int idx) {
		return tokens[idx];
	}
	
	public List<Boolean> getPhasings() {
		// the first sub-field must always be the genotype if it is present
		String[] dataTypes = vcfRecord.getFormat();
		if(dataTypes.length == 0 || !dataTypes[0].equals(FIELD_GT)) return null;
		String genotype = tokens[0];
		
		// parse phasing
		List<Boolean> phasings = new ArrayList<Boolean>(1);
		final int nrChars = genotype.length();
		for (int i = 0; i < nrChars; ++i){
			switch(genotype.charAt(i)) {
				case GENOTYPE_PHASED:
					phasings.add(Boolean.TRUE);
					break;
				case GENOTYPE_UNPHASED:
					phasings.add(Boolean.FALSE);
					break;
				default:
					break;
			}
		}
		return phasings;
	}
	
	public List<String> getAlleles() {
		if(cachedAlleles == null) {
			// the first sub-field must always be the genotype if it is present
			String[] dataTypes = vcfRecord.getFormat();
			if(dataTypes.length == 0 || !dataTypes[0].equals(FIELD_GT)) return null;
			String genotype = tokens[0];
			
			String referenceAllele = vcfRecord.getReferenceAllele();
			List<String> alternateAlleles = vcfRecord.getAlternateAlleles();
			
			// performance optimization for the common case that a sample consists of two alleles
			cachedAlleles = new ArrayList<String>(2);
			final int nrGenotypeChars = genotype.length();
			for (int j = 0, start = 0; j < nrGenotypeChars; ++j)
			{
				char c = genotype.charAt(j);
				if (c == GENOTYPE_PHASED || c == GENOTYPE_UNPHASED || j == nrGenotypeChars - 1)
				{
					if(j - start == 1) {
						// performance optimization for the common case that an allele is described by one char
						char alleleChar = j == nrGenotypeChars - 1 ?  c : genotype.charAt(j - 1);
						if(alleleChar != '.') {
							int alleleIndex = Character.digit(alleleChar, 10);
							if(alleleIndex == 0)
								cachedAlleles.add(referenceAllele);
							else
								cachedAlleles.add(alternateAlleles.get(alleleIndex - 1));
						} else {
							cachedAlleles.add(null);
						}
					} else {
						String alleleIndexStr = j == nrGenotypeChars - 1 ? genotype.substring(start) : genotype.substring(start, j);
						if (!alleleIndexStr.equals(".")) {
							int alleleIndex = Integer.parseInt(alleleIndexStr);
							if(alleleIndex == 0)
								cachedAlleles.add(referenceAllele);
							else
								cachedAlleles.add(alternateAlleles.get(alleleIndex - 1));
						} else {
							cachedAlleles.add(null);
						}
					}
					start = j + 1;
				}
			}
		}
		return cachedAlleles;
	}
	
	public VcfSample createClone() {
		return new VcfSample(vcfRecord, tokens);
	}
	
	public void reset(String[] tokens) {
		this.tokens = tokens;
		this.cachedAlleles = null;
	}
}
