#ifndef BITS_HPP
#define BITS_HPP

/// @brief Extracts the n least significant bits of x and puts them in
///         the array "bits" replacing 0 by -1
/// bits[0] = most significant bit, bits[1] = next bit, etc.
void getbits ( int x, int n, int bits[] )
{
    int pattern;

    pattern = 1 << n-1;
    for ( int i = 0; i < n; ++i, pattern >>= 1 )
        bits[i] = ( x & pattern ) ? 1 : -1;

}


/// @brief Convert array of bits to an integer
/// bits[] has same format as in getbits()
/// (i.e., bits[0] is m.s.b. and bits[n-1] is l.s.b.)
int bits2int ( int n, int bits[] )
{
    int x;
    int pattern = 1 << n-1;
    x = 0;
    for ( int i = 0; i < n; i++, pattern >>= 1 )
        if ( bits[i] == 1 )
            x |= pattern;

    return x;
}

#endif
