
/*============================================================================

This C source file is part of the SoftFloat IEEE Floating-Point Arithmetic
Package, Release 3b, by John R. Hauser.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 The Regents of the University of
California.  All Rights Reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions, and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions, and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the University nor the names of its contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS "AS IS", AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, ARE
DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=============================================================================*/

#include <xen/types.h>
#include <xen/float/platform.h>
#include <xen/float/internals.h>
#include <xen/float/softfloat.h>

#ifdef SOFTFLOAT_FAST_INT64

void ui32_to_extF80M( uint32_t a, extFloat80_t *zPtr )
{

    *zPtr = ui32_to_extF80( a );

}

#else

void ui32_to_extF80M( uint32_t a, extFloat80_t *zPtr )
{
    struct extFloat80M *zSPtr;
    uint_fast16_t uiZ64;
    uint64_t sigZ;
    int_fast8_t shiftDist;

    zSPtr = (struct extFloat80M *) zPtr;
    uiZ64 = 0;
    sigZ = 0;
    if ( a ) {
        shiftDist = softfloat_countLeadingZeros32( a );
        uiZ64 = packToExtF80UI64( 0, 0x401E - shiftDist );
        sigZ = (uint64_t) (a<<shiftDist)<<32;
    }
    zSPtr->signExp = uiZ64;
    zSPtr->signif = sigZ;

}

#endif
