package org.rdkit;

import java.io.IOException;
import java.util.Map;

import org.RDKit.ExplicitBitVect;
import org.RDKit.RDKFuncs;
import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.queries.function.FunctionValues;
import org.apache.lucene.queries.function.ValueSource;
import org.apache.lucene.search.IndexSearcher;
import org.apache.solr.search.FunctionQParser;
import org.apache.solr.search.SyntaxError;
import org.apache.solr.search.ValueSourceParser;

/**
 * Function plugin providing 'tanimoto' function.
 *
 * Use this as tanimoto(field,"FingerPrintStringLiteral"). The result is the Tanimoto similarity
 * metric between the value(s) in the field, and the BitVector represented by the string literal.
 */
public class TanimotoSimilarityParser extends ValueSourceParser {

	static {
		// Explicitly load the RDKit JNI library
		System.loadLibrary("GraphMolWrap");
	}

	private class TanimotoSimilarity extends ValueSource {
		final String NAME = "tanimoto";
		private final ValueSource source;
		private final String fps;

		TanimotoSimilarity(ValueSource source, String fps) {
			this.source = source;
			this.fps = fps;
		}

		@Override
		public FunctionValues getValues(Map context, LeafReaderContext readerContext) throws IOException {
			final FunctionValues values = source.getValues(context, readerContext);
			final ExplicitBitVect fp1 = new ExplicitBitVect(fps);
			return new FunctionValues() {
				@Override
				public float floatVal(int doc) {
					return (float)doubleVal(doc);
				}
				@Override
				public int intVal(int doc) {
					return (int)doubleVal(doc);
				}
				@Override
				public long longVal(int doc) {
					return (long)doubleVal(doc);
				}
				@Override
				public double doubleVal(int doc) {
					final String fps2 = values.strVal(doc);
					final ExplicitBitVect fp2 = new ExplicitBitVect(fps2);
					return RDKFuncs.TanimotoSimilarityEBV(fp1, fp2);
				}
				@Override
				public String strVal(int doc) {
					return Double.toString(doubleVal(doc));
				}
				@Override
				public String toString(int doc) {
					return NAME + '(' + values.toString(doc) + ',' + fps + ')';
				}
			};
		}

		@Override
		public void createWeight(Map context, IndexSearcher searcher) throws IOException {
			source.createWeight(context, searcher);
		}

		@Override
		public boolean equals(Object o) {
			if (null == o) return false;
			if (!(o instanceof TanimotoSimilarity)) {
				return false;
			}
			TanimotoSimilarity other = (TanimotoSimilarity)o;
			return (other.fps == null ? this.fps == null : other.fps.equals(this.fps)) &&
				(other.source == null? this.source == null : other.source.equals(this.source));
		}

		@Override
		public int hashCode() {
			int result = 12;
			if (this.source != null) {
				result = 37 * result + this.source.hashCode();
			}
			if (this.fps != null) {
				result = 37 * result + this.fps.hashCode();
			}
			return result;
		}

		@Override
		public String description() {
			return NAME + '(' + source.description() + ',' + "string)";
		}
	}
	@Override
	public ValueSource parse(FunctionQParser functionQParser) throws SyntaxError {
		ValueSource source = functionQParser.parseValueSource();
		String fps = functionQParser.parseArg();
		return new TanimotoSimilarity(source, fps);
	}
}
