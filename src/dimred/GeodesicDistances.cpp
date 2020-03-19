/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 Copyright (c) 2015-2019 The plumed team
 (see the PEOPLE file at the root of the distribution for a list of names)
 
 See http://www.plumed.org for more information.
 
 This file is part of plumed, version 2.
 
 plumed is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 plumed is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with plumed.  If not, see <http://www.gnu.org/licenses/>.
 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "analysis/AnalysisBase.h"
#include "tools/PDB.h"
#include "core/ActionRegister.h"
#include "reference/MetricRegister.h"
#include "reference/ReferenceConfiguration.h"
#include <stdlib.h>
#include <limits>
#include "tools/Matrix.h"

//+PLUMEDOC ANALYSIS GEODESIC_DISTANCES
/*
 Calculate the geodesic distances between data points (nodes) using Djkstra's algorithm, which finds the shortest path between configurations. In few words, the Dijkstra's algorithm takes the first node as starting point and search for the neighbors node the one with the lowest distance. When found, the first node (start) and the new node (second) are marked as processed and the distance between the nodes is stored in a variable representing the total length of the oath. Then, the distances between the second node and its neighbors are calculated and the node with the minumum distance is marked as processed (third node). The distance between the second and the third node is added to the total path lenght. Again from the new node (third), algorithm search for the smaller distance between the neighbors until all nodes have been processed, meaning that a path connecting the first and the last configuration has been obtained. For more information, please read the relative Wikipedia page.
 
 GEODESIC_DISTANCES calculates the matrix of geodesic distances between a trajectory of atomic configurations.
 Three CVs describing the system are required.
 
 The geodesic distance matrix gives an indication of the distances between two configurations in order to characterize the
 pathways of conformational changes of the system.
 
 Geodesic distances preserve the dynamics proximity in the low dimension.
 
 The geodesic distances could be calculated using two methods:
 - the neighbors method (NEIGHBORS): the distances between neighbors are computed and sorted in ascending order to introduce an order parameter. The first NEIGHBORS neighbors are connected.
 - the epsilon method (EPSILON): the epsilon value to be used as a cut-off distance between pairs should be passed and all pairs with a distance less than this value are computed.
 
 Compulsory keywords:
 CUTOFF = this is the cutoff used in the Dijkstra's algorithm. If greater than 0 the keyword MAXCONNECT must be used, otherwise NEIGHBORS.
 Optional:
 MAXCONNECT = maximum number of connections that can be formed by any given node in the graph if you are using CUTOFF.
 NEIGHBORS = the number of neighbors to use when constructing the graph that will be used when we run Djkstra's algorithm. This number should be smaller than the half of the total number of nodes. GEODESIC_DISTANCES will tell if the choice of the value is appropriate or not.
 
 
 The values obtained from the methods above are passed to the Dijkstra's algorithm that computes the the shortest
 paths between configurations.
 
 BEWARE: some configurations should be not connected each other, so the choice of NEIGHBORS and EPSILON must be done with accuracy.
 
 ERROR CODES:
 1: Your MAXCONNECT value is equal or lower than 0 and this has no physical meaning. Change the value to a positive value
 2: The value entered for NEIGHBORS is not an integer
 
 
 \par Examples
 
 The following input tells plumed to perform the GEODESIC_DISTANCES on the EUCLIDEAN_DISSIMILARITIES and
 to use CLASSICAL_MDS. This routine must be used in conjungtion with EUCLIDEAN_DISSIMILARITIES and with CLASSICAL_MDS routines.
 
 EPSILON method
 
 \plumedfile
 cv1: READ FILE=data.dat VALUES=cv1
 cv2: READ FILE=data.dat VALUES=cv2
 cv3: READ FILE=data.dat VALUES=cv3
 
 ff: COLLECT_FRAMES ARG=cv1,cv2,cv3
 oo: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=ff
 gd: GEODESIC_DISTANCES USE_OUTPUT_DATA_FROM=oo CUTOFF=2 MAXCONNECT=1
 
 mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=gd NLOW_DIM=2
 
 OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=mds ARG=mds.* FILE=list_embed FMT=%8.4f
 \endplumedfile
 
 NEIGHBORS method
 
 \plumedfile
 cv1: READ FILE=data.dat VALUES=cv1
 cv2: READ FILE=data.dat VALUES=cv2
 cv3: READ FILE=data.dat VALUES=cv3
 
 ff: COLLECT_FRAMES ARG=cv1,cv2,cv3
 oo: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=ff
 gd: GEODESIC_DISTANCES USE_OUTPUT_DATA_FROM=oo CUTOFF=0 NEIGHBORS=5
 
 mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=gd NLOW_DIM=2
 
 OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=mds ARG=mds.* FILE=list_embed FMT=%8.4f
 \endplumedfile
 
 */
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class GeodesicDistances : public AnalysisBase {
private:
    double cutoff, cutoff2;
//    double distance_neighbors;
    unsigned nneighbors;
    std::vector<unsigned> nneigh;
    std::vector<unsigned> new_parent;
    std::vector<unsigned> old_parent;
    std::vector<unsigned> dist;
    std::vector<unsigned> sptSet;
    Matrix<double> dissimilarities;
    Matrix<std::pair<double,double> > adj_list;
    Matrix<std::pair<double,double> > h;
    Matrix<std::pair<double,double> > indices;
public:
    static void registerKeywords( Keywords& keys );
    GeodesicDistances( const ActionOptions& ao );
    int n = getNumberOfDataPoints();
    
    /// Do the analysis
    void performAnalysis();
    
    /// This ensures that classes that use this data know that dissimilarities were set
    bool dissimilaritiesWereSet() const { return true; }
    
    /// Get the squared dissimilarity between two reference configurations
    double getDissimilarity( const unsigned& i, const unsigned& j );
    void performTask( const unsigned&, const unsigned&, MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(GeodesicDistances,"GEODESIC_DISTANCES")

void GeodesicDistances::registerKeywords( Keywords& keys ) {
    AnalysisBase::registerKeywords( keys );
    
    keys.add("compulsory","NEIGHBORS","the number of neighbors to use when constructing the graph that will be used when we run Djkstra's algorithm");
    keys.add("compulsory","CUTOFF","this is the distance range over which to assume that the geodesic distance can be assumed to be the same as the euclidean distance. "
             "If this keyword is specified, you do not need to include the NEIGHBORS keyword");
    keys.add("compulsory","MAXCONNECT","0","maximum number of connections that can be formed by any given node in the graph if you are using CUTOFF mode "
             "By default this is set equal to zero and the number of connections is set equal to the number "
             "of nodes. You only really need to set this if you are working with a very large system and "
             "memory is at a premium");
    
}

GeodesicDistances::GeodesicDistances( const ActionOptions& ao ):
Action(ao),
AnalysisBase(ao),
cutoff(-1),
cutoff2(-1),
nneighbors(0)
{
    if( !my_input_data ) error("no input data was found");
    if( !dissimilaritiesWereSet() ) error("dissimilarities have not been calculated in input actions");
    int n=getNumberOfDataPoints();
    parse("CUTOFF",cutoff);
    if( cutoff>0 ) {
        log.printf("  connecting nodes that are within %f of each other \n", cutoff );
        parse("MAXCONNECT",nneighbors); cutoff2 = cutoff*cutoff; cutoff=cutoff2;
        if(nneighbors<=0)
        {
            log.printf( "Your MAXCONNECT value is equal or lower than 0 and this has no physical meaning. Change the value to a positive value : ERROR 1");
            error("Error 1");
            exit(1);
        }
    } else {
        parse("NEIGHBORS",nneighbors);
        log.printf("  connecting each node to its %d nearest neighbors \n", nneighbors );
        if(nneighbors>(n/2))
        {
            log.printf("  the number of nearest neighbors chosen (%d) is greater than the half of number of nodes (%d) \n", nneighbors,n );
            log.printf( " this will probably cause a wrong conection between nodes, so check it carefully ");
        }
        if(floor(nneighbors)!=nneighbors)
        {
            log.printf( "The value entered for NEIGHBORS is not an integer : ERROR 2");
            error("Error 2");
            exit(2);
        }
    }
}


void GeodesicDistances::performAnalysis() {
    // Resize all bookeeping arrays
    dissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); dissimilarities=0;
    
    for(unsigned i=0; i<getNumberOfDataPoints(); ++i) {
        for(unsigned j=0; j<=i; ++j) my_input_data->getDissimilarity(i,j);
    }
    
    nneigh.resize( getNumberOfDataPoints() );
    
    std::vector<std::vector<std::pair<double,int> > > neighborhood;
    neighborhood.resize(getNumberOfDataPoints());
    int k=0;
    double distance_neighbors;
    double dist;
    Matrix<double> h( getNumberOfDataPoints(), getNumberOfDataPoints() );
    double a = std::numeric_limits<double>::infinity();
    Matrix<double> indices(getNumberOfDataPoints(),getNumberOfDataPoints());
    int m=getNumberOfDataPoints();
    
    //        ##################################
    //        This performs the epsilon method
    //        ##################################
    if(cutoff>0) {
        for (int i=0;i<m;i++){
            for (int j=0;j<m;j++){
                dist=getDissimilarity(i,j);
                
                if(dist<=cutoff){
                    //                    matrix h will contain the distance between neighbors if the
                    //                    distance is lower than the cutoff, otherwise if INF
                    h[i][j]=getDissimilarity(i,j);
                }
                else{
                    h[i][j]=a;
                }
            }
        }
        // end for
    }
    // end if
    
    
    //        ##################################
    //        This performs the Neighbors method
    //        ##################################
    else {
        if( nneighbors==0 ) adj_list.resize( getNumberOfDataPoints(), getNumberOfDataPoints() );
        else adj_list.resize( getNumberOfDataPoints(), nneighbors );
        
        for (int i=0; i<getNumberOfDataPoints(); i++) {
            for (int j=0; j<getNumberOfDataPoints(); j++) {
                //            set all distances in the graph to infinity
                h[i][j]=1*a;
                //            create the neighbors list
                neighborhood[i].push_back(std::make_pair(getDissimilarity(i,j),i*n+j));
                
            }
            //        sort the neighbors list in distance increasing order
            sort(neighborhood[i].begin(),neighborhood[i].begin()+n);
        }
        for (int i=0; i<getNumberOfDataPoints(); i++) {
            for (int j=0; j<getNumberOfDataPoints(); j++) {
                distance_neighbors=neighborhood[i][j].first;
                k=(neighborhood[i][j].second)-n*i;
                //          states contains for each index the nearest neighbors
                //          ordered with increaing distance, i.e.
                //          each row contains indices of points sorted by their distance in ascending order
                //          so in the first row the first number says that the first neighbor is itself, the second number
                //          is the nearest negighbor, the second number the second neighbor etc...
                
                if(j<=nneighbors){
                    h[i][k]=getDissimilarity(i,j);
                    indices[i][k]=false;
                }
                else{
                    h[i][k]=a;
                    indices[i][k]=true;
                }
            }
        }
    }
    // end if-else
    
    // Resize dissimilarities matrix and set all elements to zero
    if( !usingLowMem() ) {
        dissimilarities.resize( getNumberOfDataPoints(), getNumberOfDataPoints() ); dissimilarities=0;
    }
}
// end function


double GeodesicDistances::getDissimilarity( const unsigned& iframe, const unsigned& jframe ) {
    
    if( dissimilarities(iframe,jframe)>0. ) { return dissimilarities(iframe,jframe); }
    double tmp_distances(getNumberOfDataPoints());
    
    Matrix<double> h( getNumberOfDataPoints(), getNumberOfDataPoints() );
    new_parent.resize( getNumberOfDataPoints() );
    old_parent.resize( getNumberOfDataPoints() );
    
    plumed_dbg_assert( iframe<getNumberOfDataPoints() && jframe<getNumberOfDataPoints() );
    if( !usingLowMem() ) {
        if( dissimilarities(iframe,jframe)>0. ) { return dissimilarities(iframe,jframe); }
    }
    if( iframe!=jframe ) {
        
        //        ##################################
        // Function that implements Dijkstra's single source shortest path algorithm
        // for a graph represented using adjacency matrix representation
        
        int j=0;
        
        dist.resize( getNumberOfDataPoints() );
        double a = std::numeric_limits<double>::infinity();
        // The output array.  dist[i] will hold the shortest
        // distance from j to i
        
        sptSet.resize( getNumberOfDataPoints() );
        // sptSet[i] will be true if vertex i is included in shortest
        // path tree or shortest distance from j to i is finalized
        // Initialize all distances as INFINITE and stpSet[] as false
        for (int i = 0 ; i < getNumberOfDataPoints() ; i++){
            dist[i] = a;
            sptSet[i] = false;
        }
        
        // Distance of source vertex from itself is always 0
        dist[j] = 0;
        
        // Find shortest path for all vertices
        for (int count = 0; count < getNumberOfDataPoints()-1; count++){
            // Pick the minimum distance vertex from the set of vertices not
            // yet processed. u is always equal to j in the first iteration.
            int u,min,min_index;
            for (int v = 0; v < getNumberOfDataPoints(); v++){
                if (sptSet[v] == false && dist[v] <= INT_MAX){
                    min = dist[v], min_index = v;
                }
            }
            u=min_index;
            
            // Mark the picked vertex as processed
            sptSet[u] = true;
            
            // Update dist value of the adjacent vertices of the picked vertex.
            for (int v = 0; v < n; v++){
                
                // Update dist[v] only if is not in sptSet, there is an edge from
                // u to v, and total weight of path from src to  v through u is
                // smaller than current value of dist[v]
                if (!sptSet[v] && h[u][v] && dist[u] != INT_MAX
                    && dist[u]+h[u][v] < dist[v]){
                    
                    dist[v] = dist[u] + h[u][v];
                    tmp_distances=dist[v];
                    old_parent[v]=u;
                    new_parent[v]=v;
                    if(old_parent[v]==v)
                    {
                        break;
                    }
                    break;
                }
            }
        }
        
        
        if( !usingLowMem() ) {
            for(unsigned j=0; j<getNumberOfDataPoints(); ++j ) dissimilarities(iframe,j) = tmp_distances*tmp_distances;
        }
        // Notice that the squares of the distances are required for the other programs that might use these distances
        return tmp_distances*tmp_distances;
    }
    return 0;
}

}
// end analysis

}
//  end PLMD
