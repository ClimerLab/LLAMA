#---------------------------------------------------------------------------------------------------
# Compiler selection
#---------------------------------------------------------------------------------------------------

CXX = g++

#---------------------------------------------------------------------------------------------------
# Directories
#---------------------------------------------------------------------------------------------------

OBJDIR = build
SRCDIR = src

#---------------------------------------------------------------------------------------------------
# Executables
#---------------------------------------------------------------------------------------------------

EXE = CombClustOut SepComp FindBestComp LLAMA CombComp ConGraphToCA

#---------------------------------------------------------------------------------------------------
# Object files
#---------------------------------------------------------------------------------------------------
CCOOBJ = CombClustOut.o
AJIOBJ = AverageJI.o
SEPOBJ = SeperateGraphIntoComponents.o Graph.o Vertex.o Edge.o
FBCOBJ = FindBestComp.o
COMBOBJ = CombComp.o Chromosome.o Cluster.o Node.o Graph.o Vertex.o Edge.o UtilityFunctions.o ConfigParser.o SplitMem.o
LLAMAOBJ = LLAMA.o PopulationStat.o ClusteringResults.o Chromosome.o Cluster.o IslandStat.o Graph.o Vertex.o Edge.o Node.o UtilityFunctions.o SplitMem.o ConfigParser.o
CGCAOBJ = ConGraphToCA.o Graph.o Vertex.o Edge.o

#---------------------------------------------------------------------------------------------------
# Compiler options
#---------------------------------------------------------------------------------------------------

CXXFLAGS = -O3 -fPIC -fexceptions -DIL_STD -std=c++11#-fno-strict-aliasing

#---------------------------------------------------------------------------------------------------
all: CXXFLAGS += -DNDEBUG
all: $(EXE)

debug: CXXFLAGS += -g
debug: $(EXE)


LLAMA: $(addprefix $(OBJDIR)/, LLAMA.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(LLAMAOBJ))

SepComp: $(addprefix $(OBJDIR)/, SeperateGraphIntoComponents.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(SEPOBJ))
	
ConGraphToCA: $(addprefix $(OBJDIR)/, ConGraphToCA.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(CGCAOBJ))
	
FindBestComp: $(addprefix $(OBJDIR)/, FindBestComp.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(FBCOBJ))
	
CombComp: $(addprefix $(OBJDIR)/, CombComp.o)
	$(CXX) -o $@ $(addprefix $(OBJDIR)/, $(COMBOBJ))

CombClustOut: $(OBJDIR)/CombClustOut.o
	$(CXX) $(CXXLNDIRS) -o $@ $(addprefix $(OBJDIR)/, $(CCOOBJ)) $(CXXLNFLAGS)
	
$(OBJDIR)/LLAMA.o: $(addprefix $(SRCDIR)/, LLAMA.cpp LLAMA.h) \
						   $(addprefix $(OBJDIR)/, Graph.o ClusteringResults.o IslandStat.o UtilityFunctions.o Chromosome.o ConfigParser.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/ConGraphToCA.o: $(addprefix $(SRCDIR)/, ConvertGraphToClustAssign.cpp) \
										 $(addprefix $(OBJDIR)/, Graph.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/SeperateGraphIntoComponents.o: $(addprefix $(SRCDIR)/, SeperateGraphIntoComponents.cpp) \
										 $(addprefix $(OBJDIR)/, Graph.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/FindBestComp.o: $(SRCDIR)/FindBestComponent.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/CombComp.o: $(addprefix $(SRCDIR)/, CombineComponents.cpp) \
					  $(addprefix $(OBJDIR)/, Graph.o Chromosome.o ConfigParser.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/CombClustOut.o: $(addprefix $(SRCDIR)/, CombClustOut.cpp)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/ClusteringResults.o: $(addprefix $(SRCDIR)/, ClusteringResults.cpp ClusteringResults.h) \
							   $(addprefix $(OBJDIR)/, Chromosome.o IslandStat.o PopulationStat.o Cluster.o UtilityFunctions.o SplitMem.o Graph.o ConfigParser.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Chromosome.o: $(addprefix $(SRCDIR)/, Chromosome.cpp Chromosome.h) \
						$(addprefix $(OBJDIR)/, Cluster.o Node.o Vertex.o Edge.o Graph.o UtilityFunctions.o SplitMem.o ConfigParser.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Cluster.o: $(addprefix $(SRCDIR)/, Cluster.cpp Cluster.h) \
					 $(addprefix $(OBJDIR)/, Vertex.o Graph.o Node.o UtilityFunctions.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/IslandStat.o: $(addprefix $(SRCDIR)/, IslandStat.cpp IslandStat.h) \
						$(addprefix $(OBJDIR)/, Chromosome.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/Node.o: $(addprefix $(SRCDIR)/, Node.cpp Node.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/UtilityFunctions.o: $(addprefix $(SRCDIR)/, UtilityFunctions.cpp UtilityFunction.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/PopulationStat.o: $(addprefix $(SRCDIR)/, PopulationStat.cpp PopulationStat.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/Graph.o: $(addprefix $(SRCDIR)/, Graph.cpp Graph.h) \
				   $(addprefix $(OBJDIR)/, Vertex.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Vertex.o: $(addprefix $(SRCDIR)/, Vertex.cpp Vertex.h ) \
					$(addprefix $(OBJDIR)/, Edge.o)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Edge.o: $(addprefix $(SRCDIR)/, Edge.cpp Edge.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJDIR)/Timer.o: $(addprefix $(SRCDIR)/, Timer.cpp Timer.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/SplitMem.o: $(addprefix $(SRCDIR)/, SplitMem.cpp SplitMem.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(OBJDIR)/ConfigParser.o: $(addprefix $(SRCDIR)/, ConfigParser.cpp ConfigParser.h)
	$(CXX) $(CXXFLAGS) -c -o $@ $<


#---------------------------------------------------------------------------------------------------
.PHONY: clean
clean:
	/bin/rm -f $(OBJDIR)/*.o
#---------------------------------------------------------------------------------------------------
