    //!-------------------------------------------------------------------------
    //! Copyright (c) 2016 MOSAIC Group (MPI-CBG Dresden)
    //!
    //!
    //! This file is part of the PPM_RC program.
    //!
    //! PPM_RC is free software: you can redistribute it and/or modify
    //! it under the terms of the GNU Lesser General Public License
    //! as published by the Free Software Foundation, either
    //! version 3 of the License, or (at your option) any later
    //! version.
    //!
    //! PPM_RC is distributed in the hope that it will be useful,
    //! but WITHOUT ANY WARRANTY; without even the implied warranty of
    //! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    //! GNU General Public License for more details.
    //!
    //! You should have received a copy of the GNU General Public License
    //! and the GNU Lesser General Public License along with PPM_RC. If not,
    //! see <http://www.gnu.org/licenses/>.
    //!-------------------------------------------------------------------------
    //!  MOSAIC Group
    //!  Max Planck Institute of Molecular Cell Biology and Genetics
    //!  Pfotenhauerstr. 108, 01307 Dresden, Germany
    //!
    //!  Author           - y.afshar           June   2014
    //!-------------------------------------------------------------------------

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/foreach.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/pending/disjoint_sets.hpp>

using namespace boost;

// create a typedef for the Graph type
typedef boost::adjacency_list <
  boost::vecS            // edge list
, boost::vecS            // vertex list
, boost::undirectedS     // directedness
> Graph;

typedef graph_traits<Graph>::vertices_size_type VertexIndex;
typedef VertexIndex* Rank;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef Vertex* Parent;

int *createdgraph=NULL;

extern "C" {
  void * creategraph3d(int *npart, int *vertices_neigh, int *nvertices_neigh, int *childsize)
  {
    int nvertices = *npart;
    Graph g(nvertices);

    std::vector<VertexIndex> rank(num_vertices(g));
    std::vector<Vertex> parent(num_vertices(g));

    disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);

    initialize_incremental_components(g, ds);
    incremental_components(g, ds);

    // create edge between connected components based on mother list
    int i,j,n;
    for (i=0;i<nvertices;++i) {
      j=i*6;
      for (n=0;n<nvertices_neigh[i];++n) {
        if (vertices_neigh[j+n] == 0) break;
        boost::add_edge(i,vertices_neigh[j+n]-1,g);
        ds.union_set(i,vertices_neigh[j+n]-1);
      }
    }

    std::vector<int> child;

    typedef component_index<VertexIndex> Components;
    Components components(parent.begin(), parent.end());
    BOOST_FOREACH(VertexIndex current_index, components) {
      child.push_back(-2);
      BOOST_FOREACH(VertexIndex child_index, components[current_index]) {
        child.push_back(int(child_index));
      }
    }

    /*      for (j=0;j<child.size();++j) {std::cout << child[j] << std::endl;}*/

    *childsize=child.size();
    createdgraph=new int [*childsize];
    for (j=0;j<*childsize;++j) {
      createdgraph[j]=child[j]+1;
    }
    return createdgraph;
  }

  void * creategraph2d(int *npart, int *vertices_neigh, int *nvertices_neigh, int *childsize)
  {
    int nvertices = *npart;
    Graph g(nvertices);

    std::vector<VertexIndex> rank(num_vertices(g));
    std::vector<Vertex> parent(num_vertices(g));

    disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);

    initialize_incremental_components(g, ds);
    incremental_components(g, ds);

    // create edge between connected components based on mother list
    int i,j,n;
    for (i=0;i<nvertices;++i) {
      j=i*4;
      for (n=0;n<nvertices_neigh[i];++n) {
        if (vertices_neigh[j+n] == 0) break;
        boost::add_edge(i,vertices_neigh[j+n]-1,g);
        ds.union_set(i,vertices_neigh[j+n]-1);
      }
    }

    std::vector<int> child;

    typedef component_index<VertexIndex> Components;
    Components components(parent.begin(), parent.end());
    BOOST_FOREACH(VertexIndex current_index, components) {
      child.push_back(-2);
      BOOST_FOREACH(VertexIndex child_index, components[current_index]) {
        child.push_back(int(child_index));
      }
    }
    /*for (j=0;j<child.size();++j) {std::cout << child[j] << std::endl;}*/
    *childsize=child.size();
    createdgraph=new int [*childsize];
    for (j=0;j<*childsize;++j) {
      createdgraph[j]=child[j]+1;
    }
    return createdgraph;
  }

  int destroygraph()
  {
    if (createdgraph != NULL) {
      delete[] createdgraph;
      createdgraph = NULL;
    }
    return 0;
  }
}


