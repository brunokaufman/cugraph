#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <string>

#include <thrust/device_vector.h>

#include "sp_matrix.h"

template<typename N, typename T>
struct graph{
	SPM<N,T> edges_matrix;//La matriz de adyacencia.
	std::vector<std::string> node_labels;//Los labels de los nodos de cada lado.

	N n_labels, n_edges;

	graph(){//Constructor vacio.
		resize(0, 0);
	}

	graph(const std::string& filename, const char &delim = ',', bool weighted = false, bool directed = false){//Constructor desde archivo.
		read_csv<N,T>(filename, node_labels, edges_matrix, delim, weighted, directed);
		n_labels = node_labels.size();
		n_edges = edges_matrix.n_entries;
	}

	graph(const graph<N,T>& _graph){//Constructor de copia.
		int i=0;
		resize(_graph.n_labels, _graph.n_edges);
		edges_matrix.copy(_graph.edges_matrix);
		std::copy(_graph.node_labels.begin(), _graph.node_labels.end(), node_labels.begin());
	}

	void resize(const N& _n_labels, const N& _n_edges){
		n_labels = _n_labels;
		n_edges = _n_edges;

		edges_matrix.resize(n_labels, n_labels, n_edges);//Es cuadrada por ser adyacencia, y no incidencia. Va a tener muchos ceros.
		node_labels.resize(n_labels);
	}
};

#endif
