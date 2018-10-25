#ifndef SPM_H
#define SPM_H

#include<iostream>
#include<iterator>

#include<thrust/execution_policy.h>
#include<thrust/host_vector.h>
#include<thrust/device_vector.h>
#include<thrust/tuple.h>
#include<thrust/pair.h>

#include<thrust/fill.h>
#include<thrust/copy.h>
#include<thrust/sort.h>
#include<thrust/reduce.h>

#include<thrust/iterator/zip_iterator.h>
#include<thrust/iterator/constant_iterator.h>

typedef typename thrust::device_vector<N> NDVec;
typedef typename thrust::device_vector<T> TDVec;
typedef typename thrust::host_vector<N> NHVec;

template<typename N, typename T>
struct triplet_compare{
	bool row_first;

	triplet_compare(bool _row_first){
		row_first = _row_first;
	}

	__device__ bool operator()(const thrust::tuple<N,N,T>& x, const thrust::tuple<N,N,T>& y){
		if(row_first){
			return (thrust::get<0>(x) < thrust::get<0>(y)) || ((thrust::get<0>(x) == thrust::get<0>(y)) && (thrust::get<1>(x) < thrust::get<1>(y)));
		}
		else{
			return (thrust::get<1>(x) < thrust::get<1>(y)) || ((thrust::get<1>(x) == thrust::get<1>(y)) && (thrust::get<0>(x) < thrust::get<0>(y)));
		}
	}
};

template<typename N, typename T>
struct SPM{
	private:
	N n_rows, n_cols, n_entries;//Informacion basica de la matriz.

	NDVec row_indices, col_indices, col_offsets;//Informacion indicial relevante a formato COO y CSC (voy a tener ambos al mismo tiempo por facilidad).
	TDVec values;//Los valores a utilizar.

	NHVec h_col_offsets;//Utilizo esto para pasar de offsets a indices (de CSC a COO).
	NDVec col_degrees;//Utilizo esto para guardar los grados en fila (para pasar de COO a CSC).

	//Ordeno las filas y columnas (default, fila primero).
	void sort(bool row_first = true){
		typedef typename NDVec::iterator IndexIterator;
		typedef typename TDVec::iterator ValueIterator;
		typedef thrust::tuple<IndexIterator, IndexIterator, ValueIterator> TripletTuple;
		typedef thrust::zip_iterator<TripletTuple> ZipIndex;

		ZipIndex zip_triplet(thrust::make_tuple(row_indices.begin(), column_indices.begin(), values.begin()));
		thrust::sort(thrust::device, zip_triplet, zip_triplet + n_entries, triplet_compare<N,T>(row_first));
	}

	//Funcion resize estandar, incluyendo el resize adecuado para COO y para CSC (indices columna en el primero, offsets en el segundo).
	void resize(N _n_rows, N _n_cols, N _n_entries){
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_entries = _n_entries;

		row_indices.resize(n_entries);
		col_indices.resize(n_entries);
		values.resize(n_entries);
		col_offsets.resize(n_cols);

		h_col_offsets.resize(n_cols);
		row_degrees.resize(n_cols);
	}

	//Funcion para calcular los grados en columna. En una red dirigida, seria el in-degree.
	void calculate_col_degrees(){
		typedef NDVec::iterator IndexIterator;
		thrust::pair<IndexIterator, IndexIterator> new_ends;
		sort_column_first();//Ordeno las columnas para que el reduce_by_key me las procese bien.
		node_keys.resize(n_cols);//Esto me va a almacenar los indices (importante si llego a tener vacia una columna).
		new_ends = thrust::reduce_by_key(thrust::device, column_indices.begin(), column_indices.end(), thrust::make_constant_iterator<N>(1), node_keys.begin(), col_degrees.begin());

	}

	//Transformacion de offsets (CSC) a indices (COO).
	void offsets_to_indices(){
		for(int i=0; i<n_cols-1; i++){
			thrust::fill(col_indices.begin() + h_col_offsets[i], col_indices.begin() + h_col_offsets[i+1], i);
		}
		thrust::fill(col_indices.begin() + h_col_offsets[n_cols-1], col_indices.end(), n_cols-1);
	}

	//Transformacion de indices (COO) a offsets (CSC).
	void indices_to_offsets(){
		thrust::copy(row_degrees.begin(), row_degrees.end() - 1, col_offsets.begin() + 1);
		thrust::inclusive_scan(thrust::device, col_offsets.begin(), col_offsets.end(), col_offsets.begin());
		thrust::copy(col_offsets.begin(), col_offsets.end(), h_col_offsets.begin());
	}

	public:
	//Null constructor si sos valiente.
	SPM(){}

	//Constructor desde tamanos.
	SPM(N _n_rows, N _n_cols, N _n_entries){
		this.SPM();
		resize(_n_rows, _n_cols, _n_entries);
	}

	//Constructor desde otra SPM.
	SPM(SPM<N,T> _mat){
		copy(_mat);
	}

	//Constructor basado en la copia basica.
	SPM(N _n_rows, N _n_cols, NDVec _row_indices, NDVec _col_ind_off, TDVec _values, bool offsets = false){
		copy(N _n_rows, N _n_cols, NDVec _row_indices, NDVec _col_ind_off, TDVec _values, bool offsets = false);
	}

	//Funcion de copia desde otra matriz.
	void copy(SPM<N,T> _mat){
		copy(_mat.n_rows, _mat.n_cols, _mat.row_indices, _mat.col_indices, _mat.values);//Uso la default, donde uso los indices de columna (COO) en vez de los offsets (CSC).
	}

	//Funcion copia basica desde indices/offsets y valores relevantes. Por default el input es en indices (COO).
	void copy(N _n_rows, N _n_cols, NDVec _row_indices, NDVec _col_ind_off, TDVec _values, bool offsets = false){
		resize(_n_rows, _n_cols, _values.size());

		thrust::copy(_row_indices.begin(), _row_indices.end(), row_indices.begin());
		thrust::copy(_values.begin(), _values.end(), values.begin());

		if(offsets){
			thrust::copy(_col_ind_off.begin(), _col_ind_off.end(), col_offsets.begin());
			thrust::copy(col_offsets.begin(), col_offsets.end(), h_col_offsets.begin());
			offsets_to_indices();
			calculate_row_degrees();
		}
		else{
			thrust::copy(_col_ind_off.begin(), _col_ind_off.end(), col_indices.begin());
			calculate_col_degrees();
			indices_to_offsets();
			thrust::copy(col_offsets.begin(), col_offsets.end(), h_col_offsets.begin());
		}
	}

	//Trasponer cambiando los elementos de fila y de columna.
	void transpose(){
		thrust::swap_ranges(row_indices.begin(), row_indices.end(), column_indices.begin());
		resize(n_cols, n_rows, n_entries);//No me olvido de cambiar los tamanos de fila y columna.
		update_col_offsets();
	}

	//Funcion impresion para diagnosticos.
	void print(int _n){
		int n = values.size();
		if(_n < n){
			n = _n;
		}

		std::cout << "Matriz de (" << n_rows << "x" << n_cols << ") con " << n_entries << " entradas:" << std::endl;

		std::cout << "   Indices fila: ";
		std::copy_n(row_indices.begin(), n, std::ostream_iterator<N>(std::cout, " "));
		std::cout << std::endl << "Indices columna: ";
		std::copy_n(column_indices.begin(), n, std::ostream_iterator<N>(std::cout, " "));
		std::cout << std::endl << "        Valores: ";
		std::copy_n(values.begin(), n, std::ostream_iterator<T>(std::cout, " "));
		std::cout << std::endl;
	}
};

#endif
