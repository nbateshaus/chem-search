package org.rdkit;

import java.io.IOException;
import java.util.Map;

import org.apache.lucene.index.LeafReaderContext;
import org.apache.lucene.queries.function.FunctionValues;
import org.apache.lucene.queries.function.ValueSource;
import org.apache.lucene.search.IndexSearcher;
import org.apache.solr.search.FunctionQParser;
import org.apache.solr.search.SyntaxError;
import org.apache.solr.search.ValueSourceParser;

/**
 * Created by nik on 4/6/16.
 */
public class TanimotoSimilarityParser extends ValueSourceParser {
	public class TanimotoSimilarity extends ValueSource {
		public final String NAME = "tanimoto";
		private final ValueSource source;
		private final String fps;

		public TanimotoSimilarity(ValueSource source, String fps) {
			this.source = source;
			this.fps = fps;
		}

		@Override
		public FunctionValues getValues(Map context, LeafReaderContext readerContext) throws IOException {
			final FunctionValues values = source.getValues(context, readerContext);
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
					String fps2 = values.strVal(doc);
					return fps.length() - fps2.length(); // TODO nik: RDKit Tanimoto Similarity here
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
