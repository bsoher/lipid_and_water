function [ floatVax ] = uint32le_to_VAXfloat( uint32le )
%UINT32LE_TO_VAXFLOAT Converts from IEEE-LE (UINT32) to VAX (float)
%  This function takes a raw 32bit unsigned integer (little endian)
%  and converts it into the equivalent floating point number it represents
%  in the VAX file format for floating point numbers. This is the VAXF
%  format and is the 32bit floating point format used in both VAXD and VAXG
%  files. (VAXD and VAXG are 64bit floating point formats).

%% Define floating value properties for VAX architecture
% The generic equation for a floating point number is:
% (-1)^double(S) * (F+C) * A^(double(E)-B);
% Different operating systems and file formats utilize different values
% for A, B, and C. F, E, and S are computed from the appropriate bits in
% the number as stored on disk.

    A = 2   ;%VAX specific
    B = 128 ;%VAX specific
    C = 0.5 ;%VAX specific

%% Convert raw unsigned number into right answer
% Flip the upper and lower bits (based on how Vax data storage format)
% VAX      <-----WORD1-----><-----WORD2----->
% IEEE-LE  <-----WORD2-----><-----WORD1----->

    word2  = bitshift(bitshift(uint32le,  0), -16);%mask FFFF0000
    word1  = bitshift(bitshift(uint32le, 16), -16);%mask 0000FFFF
    vaxInt = bitor(bitshift(word1,16), bitshift(word2, 0));
    
% Pull out the sign, exponent, and fractional component
% VAX FLOAT BYTES  <-----WORD1----><-----WORD2---->
% VAX FLOAT BITS   0123456789ABCDEF0123456789ABCDEF
% Sign Exp Fract   SEEEEEEEEFFFFFFFFFFFFFFFFFFFFFFF

    S = bitshift(bitshift(vaxInt , 0), -31);%2147483648=hex2dec('80000000')
    E = bitshift(bitshift(vaxInt , 1), -24);%2139095040=hex2dec('7F800000')
    F = bitshift(bitshift(vaxInt , 9),  -9);%   8388607=hex2dec('007FFFFF')

% Construct the floating point number from SEF (Sign, Exp, Fract)
% http://www.codeproject.com/KB/applications/libnumber.aspx
% http://www.opengroup.org/onlinepubs/9629399/chap14.htm#tagfcjh_20

    M = C+double(F)./2^24;%VAX Specific 16777216=2^24
    floatVax = (-1).^double(S) .* M .* A.^(double(E)-B);%Generic

end
%#eml