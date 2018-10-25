#ifndef IMPORT_H
#define IMPORT_H

#include<iostream>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<vector>
#include<string>

#include"sp_matrix.h"

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
	str.erase(0, str.find_first_not_of(chars));
	return str;
}

std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
	str.erase(str.find_last_not_of(chars) + 1);
	return str;
}

std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r "){
	return ltrim(rtrim(str, chars), chars);
}

template <typename N, typename T>
void read_csv(const std::string &filename_edgelist, std::vector<std::string>& node_labels, SPM<N,T>& edges_matrix, const char &delim = ',', bool weighted = false, bool directed = false){//Version donde toma en cuenta a los pesos de las edges como la tercera columna del CSV.
	//Almaceno los valores de los edges en el mismo indice de estos dos arrays para luego escribir a la matriz. Esto evita hacer resize() dinamico sobre la matriz, que en ublas/matrix_sparse no esta implementado preservando los valores anteriormente escritos.
	std::vector<N> edges_1(0);
	std::vector<N> edges_2(0);
	std::vector<T> edges_values(0);

	node_labels.resize(0);
	edges_matrix.resize(0, 0, 0);

	std::ifstream edgelist_file(filename_edgelist.c_str(), std::ifstream::in);
	if(edgelist_file.is_open()){//Se fija si pudo leer el archivo.
		std::cout << "Archivo " << filename_edgelist << " abierto." << std::endl;

		std::string line;
		while(std::getline(edgelist_file, line)){//Lee lineas hasta que no pueda (EOF).
			std::istringstream line_stream(line);//Arma un stream con la linea.
			std::string dummy;

			std::getline(line_stream, dummy, delim);//Lee el stream linea hasta el delim.
			dummy = trim(dummy);//Quito el whitespace al pedo.
			N i = 0;
			//Ahora voy a buscar si existe dummy ya en mi vector de labels de nodos. Todos los nodos (aunque sea bipartita) van a este conjunto de labels. La matriz de adyacencia (o incidencia) es lo que se encarga de separarlos.
			if(node_labels.size() > 0){
				i = std::find(node_labels.begin(), node_labels.end(), dummy) - node_labels.begin();
			}

			if(i == node_labels.size()){//Si no existe, lo agrego y ya tengo el valor nuevo.
				node_labels.push_back(dummy);//Agrego a la lista de nodos.
			}

			//Hago lo mismo para el lado opuesto del edge.
			std::getline(line_stream, dummy, delim);
			dummy = trim(dummy);
			N j = 0;
			if(node_labels.size() > 0){
				j = std::find(node_labels.begin(), node_labels.end(), dummy) - node_labels.begin();
			}

                        if(j == node_labels.size()){//Si no existe, lo agrego y ya tengo el valor nuevo.
				node_labels.push_back(dummy);
                        }

			edges_1.push_back(i);
			edges_2.push_back(j);

			if(weighted){
				//En este caso (weighted) leo tambien el valor del peso de los edges.
				std::getline(line_stream, dummy, delim);
				edges_values.push_back(static_cast<T>(strtof(dummy.c_str(), 0)));
			}
			else{
				edges_values.push_back(static_cast<T>(1.0));
			}
		}
	}
	else{//Avisa si no pudo leer el archivo.
		std::cout << "Error abriendo archivo " << filename_edgelist << "." << std::endl;
	}

	edges_matrix.copy(node_labels.size(), node_labels.size(), edges_1.size(), edges_1, edges_2, edges_values, !directed);
}

#endif
