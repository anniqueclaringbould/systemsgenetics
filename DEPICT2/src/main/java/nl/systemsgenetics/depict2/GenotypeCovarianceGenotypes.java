/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package nl.systemsgenetics.depict2;

import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

/**
 *
 * @author patri
 */
public class GenotypeCovarianceGenotypes implements GenotypeCovarianceSource {

	private final RandomAccessGenotypeData referenceGenotypes;

	public GenotypeCovarianceGenotypes(RandomAccessGenotypeData referenceGenotypes) {
		this.referenceGenotypes = referenceGenotypes;
	}

	@Override
	public GenotypieCovarianceResult getCovarianceMatrixForRange(String chr, int start, int stop, double doubleMaxR) {

		final SimpleRegression regression = new SimpleRegression();

		final ArrayList<GeneticVariant> includedVariantsList = new ArrayList<>();
		

		newVariants:
		for (GeneticVariant newVariant : referenceGenotypes.getVariantsByRange(chr, start, stop)) {
			final float[] newVariantDosages = newVariant.getSampleDosages();
			for (GeneticVariant selectedVariant : includedVariantsList) {
				final float[] selectedVariantDosages = selectedVariant.getSampleDosages();

				//This loop is used to calculate correlation between variant dosages
				for (int i = 0; i < newVariantDosages.length; ++i) {
					regression.addData(newVariantDosages[i], selectedVariantDosages[i]);
				}

				//If correlation is tool large stop with current newVariant and move to next variant
				if (regression.getR() >= doubleMaxR) {
					continue newVariants;
				}
			}
		}

		if (includedVariantsList.isEmpty()) {
			return new GenotypieCovarianceResult(new double[0][0], new TIntHashSet(0));
		} else {
			
			final TIntHashSet includedVariantsPositions = new TIntHashSet(includedVariantsList.size());

			final double[][] cov = new double[includedVariantsList.size()][includedVariantsList.size()];
			for (int p = 0; p < includedVariantsList.size(); p++) {
				
				GeneticVariant pVariant = includedVariantsList.get(p);
				
				if(includedVariantsPositions.add(pVariant.getStartPos())){
					throw new RuntimeException("Cannot handle reference genotype data with multiple variants at same position: " + pVariant.getSequenceName() + ":" + pVariant.getStartPos());
				}
				
				final float[] pVariantDosages = pVariant.getSampleDosages();
				cov[p][p] = 1;
				for (int q = p + 1; q < includedVariantsList.size(); q++) {
					final float[] qVariantDosages = includedVariantsList.get(q).getSampleDosages();

					//This loop is used to calculate correlation between variant dosages
					for (int i = 0; i < pVariantDosages.length; ++i) {
						regression.addData(pVariantDosages[i], qVariantDosages[i]);
					}

					cov[p][q] = regression.getR();
					cov[q][p] = cov[p][q];
				}
			}

			return new GenotypieCovarianceResult(cov, includedVariantsPositions);

		}

	}

}
