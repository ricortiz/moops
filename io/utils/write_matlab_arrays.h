#ifndef WRITE_MATLAB_ARRAYS_H
#define WRITE_MATLAB_ARRAYS_H
/***************************************************************************
 *   Copyright (C) 2007 by Ricardo Ortiz                                   *
 *   ricardo.ortiz@tulane.edu                                              *
 *                                                                         *
 ***************************************************************************/
#include <fstream>
#include <sstream>
#include <cmath>

/// @param A sparse matrix in storage format
/// @return string stream
template<typename SpMatrix>
bool write_matlab_matrix( std::string const &file, const SpMatrix & A )
{

    std::ofstream output( file.c_str() );

//     output.open( file.c_str() );

    if ( !output.is_open() || output.bad() )
        throw std::invalid_argument( "write_matrix(): Invalid file argument" );

    typedef typename SpMatrix::size_type size_type;
    typedef typename SpMatrix::value_type value_type;
    //--- get dimensions
    size_type m = A.size1();
    size_type n = A.size2();
    std::stringstream array;
    array.precision(30);
    for (size_type row = 0; row < m; ++row)
    {
        array << "[";
        for (size_type col = 0; col < n-1; ++col)
            array << A(row,col) << ",";
        array << A(row,n-1);
        array << "];";
    }

    output << "A = [" << array.str() << "];" << std::endl;
    output.flush();
    output.close();
    return true;
}

/**
 * @brief Output Default Vector.
 *  Writes the vector out in Matlab notation.
 *
 * @param v    The vector to be dumped.
 * @return   A string.
 */
template<typename Vector>
bool write_matlab_vector( std::string const &file, const Vector & v )
{
    std::stringstream stream;
    stream.precision( 30 );

    typedef typename Vector::const_iterator iterator;
//     stream << v.size() << std::endl;
    iterator vi = v.begin();
    for ( ; vi != v.end(); ++vi )
        stream << *vi << std::endl;

    std::ofstream output( file.c_str() );
    if ( !output.is_open() || output.bad() )
        throw std::invalid_argument( "write_vector(): Invalid file argument" );

    output.precision( 30 );
    output << stream.str() << std::endl;

    output.flush();
    output.close();
    return true;
};

#endif // WRITE_DATA
