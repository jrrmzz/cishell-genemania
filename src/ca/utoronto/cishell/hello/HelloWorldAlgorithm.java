package ca.utoronto.cishell.hello;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Dictionary;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.cishell.framework.CIShellContext;
import org.cishell.framework.algorithm.Algorithm;
import org.cishell.framework.algorithm.AlgorithmExecutionException;
import org.cishell.framework.data.BasicData;
import org.cishell.framework.data.Data;
import org.cishell.framework.data.DataProperty;
import org.genemania.domain.AttributeGroup;
import org.genemania.domain.Interaction;
import org.genemania.domain.InteractionNetwork;
import org.genemania.domain.InteractionNetworkGroup;
import org.genemania.domain.Node;
import org.genemania.domain.Organism;
import org.genemania.dto.RelatedGenesEngineRequestDto;
import org.genemania.dto.RelatedGenesEngineResponseDto;
import org.genemania.engine.Mania2;
import org.genemania.engine.cache.DataCache;
import org.genemania.engine.cache.MemObjectCache;
import org.genemania.engine.cache.SynchronizedObjectCache;
import org.genemania.exception.ApplicationException;
import org.genemania.exception.DataStoreException;
import org.genemania.plugin.FileUtils;
import org.genemania.plugin.GeneMania;
import org.genemania.plugin.NetworkUtils;
import org.genemania.plugin.apps.IQueryErrorHandler;
import org.genemania.plugin.controllers.DefaultGeneProvider;
import org.genemania.plugin.controllers.IGeneProvider;
import org.genemania.plugin.cytoscape.NullCytoscapeUtils;
import org.genemania.plugin.data.DataSet;
import org.genemania.plugin.data.DataSetManager;
import org.genemania.plugin.data.lucene.LuceneDataSetFactory;
import org.genemania.plugin.model.Group;
import org.genemania.plugin.model.Network;
import org.genemania.plugin.model.SearchResult;
import org.genemania.plugin.parsers.Query;
import org.genemania.plugin.parsers.TabDelimitedQueryParser;
import org.genemania.type.ScoringMethod;
import org.genemania.util.NullProgressReporter;
import org.osgi.service.log.LogService;
import org.xml.sax.SAXException;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.impl.UndirectedSparseEdge;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.impl.UndirectedSparseVertex;
import edu.uci.ics.jung.utils.UserData;

public class HelloWorldAlgorithm implements Algorithm {
    private Dictionary<?, ?> parameters;
	private LogService logger;
    
    public HelloWorldAlgorithm(Data[] data,
    				  Dictionary<?, ?> parameters,
    				  CIShellContext ciShellContext) {
        this.parameters = parameters;
        
        logger = (LogService) ciShellContext.getService(LogService.class.getName());
    }

    public Data[] execute() throws AlgorithmExecutionException {
    	try {
	    	String geneList = (String) parameters.get("geneList");
	    	
	    	DataSetManager dataSetManager = createDataSetManager();
	    	DataSet data = dataSetManager.open(new File("/Users/jay/genemania_plugin/gmdata-2012-08-02-core"));
	    	
	    	Query query = createQuery(geneList, data, parameters);
			NetworkUtils networkUtils = new NetworkUtils();
			SearchResult result = runAlgorithm(data, networkUtils, query);
	    	IGeneProvider geneProvider = new DefaultGeneProvider(networkUtils);
			Graph graph = buildGraph(result, geneProvider);
	        return new Data[] { prepareOutputData(graph, "GeneMANIA result") };
    	} catch (DataStoreException e) {
    		throw new AlgorithmExecutionException(e);
		} catch (SAXException e) {
    		throw new AlgorithmExecutionException(e);
    	} catch (ApplicationException e) {
    		throw new AlgorithmExecutionException(e);
		} catch (IOException e) {
    		throw new AlgorithmExecutionException(e);
		}
    }
    
    private Query createQuery(String geneList, DataSet data, Dictionary parameters) throws DataStoreException, IOException {
    	boolean first = true;
    	StringBuilder builder = new StringBuilder();
    	for (String name : geneList.split("\\s")) {
    		if (first) {
    			first = false;
    		} else {
    			builder.append("\t");
    		}
    		builder.append(name.trim());
    	}
    	TabDelimitedQueryParser parser = new TabDelimitedQueryParser();
    	String queryText = "9606\n" + builder.toString() + "\ndefault\n20\nbp";
		Reader reader = new StringReader(queryText);
		return parser.parse(data, reader, new IQueryErrorHandler() {
			@Override
			public void warn(String arg0) {
			}
			
			@Override
			public void handleUnrecognizedNetwork(String arg0) {
			}
			
			@Override
			public void handleUnrecognizedGene(String arg0) {
			}
			
			@Override
			public void handleSynonym(String arg0) {
			}
			
			@Override
			public void handleNetwork(InteractionNetwork arg0) {
			}
		});
	}

	private Graph buildGraph(SearchResult result, IGeneProvider geneProvider) {
		UndirectedSparseGraph graph = new UndirectedSparseGraph();
		Map<Long, Vertex> vertexCache = new HashMap<Long, Vertex>();
		Set<Pair> edges = new HashSet<Pair>();
		
		for (InteractionNetworkGroup group : result.getInteractionNetworkGroups().values()) {
			for (InteractionNetwork network : group.getInteractionNetworks()) {
				for (Interaction interaction : network.getInteractions()) {
					Node fromNode = interaction.getFromNode();
					Node toNode = interaction.getToNode();
		
					Pair pair = new Pair(fromNode.getId(), toNode.getId());
					if (edges.contains(pair)) {
						continue;
					}
					edges.add(pair);
					
					Vertex from = getVertex(vertexCache, fromNode, graph, geneProvider);
					Vertex to = getVertex(vertexCache, toNode, graph, geneProvider);
					
			    	UndirectedSparseEdge edge = new UndirectedSparseEdge(from, to);
			    	edge.removeUserDatum("labelvisible");
			    	edge.removeUserDatum("label");
			    	edge.setUserDatum("labelvisible", "false", UserData.CLONE);
			    	edge.setUserDatum("weight", interaction.getWeight(), UserData.CLONE);
			    	graph.addEdge(edge);
				}
			}
		}
    	return graph;
	}

	private Vertex getVertex(Map<Long, Vertex> vertexCache, Node node, Graph graph, IGeneProvider geneProvider) {
		Vertex vertex = vertexCache.get(node.getId());
		if (vertex != null) {
			return vertex;
		}
		vertex = new UndirectedSparseVertex();
		vertex.setUserDatum("label", geneProvider.getGene(node).getSymbol(), UserData.CLONE);
		graph.addVertex(vertex);
		vertexCache.put(node.getId(), vertex);
		return vertex;
	}

	private Data prepareOutputData(Graph outputGraph, String label) {
        Data outputData = new BasicData(outputGraph, Graph.class.getName());
 
        Dictionary metadata = outputData.getMetadata();
        metadata.put(DataProperty.TYPE, DataProperty.NETWORK_TYPE);
        metadata.put(DataProperty.LABEL, label);
 
        return outputData;
    }

	private SearchResult runAlgorithm(DataSet data, NetworkUtils networkUtils, Query query) throws DataStoreException, ApplicationException {
		RelatedGenesEngineRequestDto request = createRequest(query);
		RelatedGenesEngineResponseDto response = runQuery(data, networkUtils, request);
		
		List<String> queryGenes = query.getGenes();
		Organism organism = query.getOrganism();
		SearchResult options = networkUtils.createSearchOptions(organism, request, response, null, data, queryGenes);
		return options;
	}

	RelatedGenesEngineResponseDto runQuery(DataSet data, NetworkUtils networkUtils, RelatedGenesEngineRequestDto request) throws DataStoreException, ApplicationException {
		request.setProgressReporter(NullProgressReporter.instance());
		
		DataCache cache = new DataCache(new SynchronizedObjectCache(new MemObjectCache(data.getObjectCache(NullProgressReporter.instance(), false))));
		Mania2 mania = new Mania2(cache);
		RelatedGenesEngineResponseDto result = mania.findRelated(request);
		request.setCombiningMethod(result.getCombiningMethodApplied());
		networkUtils.normalizeNetworkWeights(result);
		return result;
	}

	RelatedGenesEngineRequestDto createRequest(Query query) throws ApplicationException {
		RelatedGenesEngineRequestDto request = new RelatedGenesEngineRequestDto();
		request.setNamespace(GeneMania.DEFAULT_NAMESPACE);
		request.setOrganismId(query.getOrganism().getId());
		
		
		request.setInteractionNetworks(collapseNetworks(query.getNetworks()));
		request.setPositiveNodes(query.getNodes());
		request.setLimitResults(query.getGeneLimit());
		request.setAttributesLimit(query.getAttributeLimit());
		request.setCombiningMethod(query.getCombiningMethod());
		
		if (query.getScoringMethod() != null) {
			request.setScoringMethod(query.getScoringMethod());
		} else {
			request.setScoringMethod(ScoringMethod.DISCRIMINANT);
		}
		return request;
	}

	protected Collection<Collection<Long>> collapseNetworks(Collection<? extends Group<?, ?>> networks) {
		Collection<Collection<Long>> result = new ArrayList<Collection<Long>>();
		
		List<Group<InteractionNetworkGroup, InteractionNetwork>> groups = new ArrayList<Group<InteractionNetworkGroup, InteractionNetwork>>();
		for (Group<?, ?> group : networks) {
			Group<InteractionNetworkGroup, InteractionNetwork> adapted = group.adapt(InteractionNetworkGroup.class, InteractionNetwork.class);
			if (adapted == null) {
				continue;
			}
			groups.add(adapted);
		}
		
		Collections.sort(groups, new Comparator<Group<InteractionNetworkGroup, InteractionNetwork>>() {
			@Override
			public int compare(Group<InteractionNetworkGroup, InteractionNetwork> g1, Group<InteractionNetworkGroup, InteractionNetwork> g2) {
				return (int) Math.signum(g1.getModel().getId() - g2.getModel().getId());
			}
		});
		
		for (Group<InteractionNetworkGroup, InteractionNetwork> group : groups) {
			Collection<Long> groupMembers = new HashSet<Long>();
			for (Network<InteractionNetwork> network : group.getNetworks()) {
				groupMembers.add(network.getModel().getId());
			}
			if (!groupMembers.isEmpty()) {
				List<Long> sorted = new ArrayList<Long>(groupMembers);
				Collections.sort(sorted);
				result.add(sorted);
			}
		}
		return result;
	}
	
	protected Collection<Long> collapseAttributeGroups(Collection<Group<?, ?>> groups) {
		List<Long> result = new ArrayList<Long>();
		for (Group<?, ?> group : groups) {
			Group<Object, AttributeGroup> adapted = group.adapt(Object.class, AttributeGroup.class);
			if (adapted == null) {
				continue;
			}
			for (Network<AttributeGroup> network : adapted.getNetworks()) {
				result.add(network.getModel().getId());
			}
		}
		return result;
	}

	static DataSetManager createDataSetManager() {
		DataSetManager dataSetManager = new DataSetManager();
		dataSetManager.addDataSetFactory(new LuceneDataSetFactory<Object, Object, Object>(dataSetManager, null, new FileUtils(), new NullCytoscapeUtils<Object, Object, Object>(), null), Collections.emptyMap());
		return dataSetManager;
	}
	
	static class Pair {
		private long from;
		private long to;

		public Pair(long from, long to) {
			if (to < from) {
				this.from = to;
				this.to = from;
			} else {
				this.from = from;
				this.to = to;
			}
		}
		
		@Override
		public boolean equals(Object obj) {
			Pair other = (Pair) obj;
			return other.from == from && other.to == to;
		}
		
		@Override
		public int hashCode() {
			return (int) (from * 31 + to);
		}
	}
}