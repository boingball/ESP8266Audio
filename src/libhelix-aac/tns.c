/* ***** BEGIN LICENSE BLOCK *****
    Source last modified: $Id: tns.c,v 1.2 2005/05/24 16:01:55 albertofloyd Exp $

    Portions Copyright (c) 1995-2005 RealNetworks, Inc. All Rights Reserved.

    The contents of this file, and the files included with this file,
    are subject to the current version of the RealNetworks Public
    Source License (the "RPSL") available at
    http://www.helixcommunity.org/content/rpsl unless you have licensed
    the file under the current version of the RealNetworks Community
    Source License (the "RCSL") available at
    http://www.helixcommunity.org/content/rcsl, in which case the RCSL
    will apply. You may also obtain the license terms directly from
    RealNetworks.  You may not use this file except in compliance with
    the RPSL or, if you have a valid RCSL with RealNetworks applicable
    to this file, the RCSL.  Please see the applicable RPSL or RCSL for
    the rights, obligations and limitations governing use of the
    contents of the file.

    This file is part of the Helix DNA Technology. RealNetworks is the
    developer of the Original Code and owns the copyrights in the
    portions it created.

    This file, and the files included with this file, is distributed
    and made available on an 'AS IS' basis, WITHOUT WARRANTY OF ANY
    KIND, EITHER EXPRESS OR IMPLIED, AND REALNETWORKS HEREBY DISCLAIMS
    ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION, ANY WARRANTIES
    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, QUIET
    ENJOYMENT OR NON-INFRINGEMENT.

    Technology Compatibility Kit Test Suite(s) Location:
      http://www.helixcommunity.org/content/tck

    Contributor(s):

 * ***** END LICENSE BLOCK ***** */

/**************************************************************************************
    Fixed-point HE-AAC decoder
    Jon Recker (jrecker@real.com)
    February 2005

    tns.c - apply TNS to spectrum
 **************************************************************************************/

#include "coder.h"
#include "assembly.h"

#define FBITS_LPC_COEFS	20

//fb
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"

/*  inverse quantization tables for TNS filter coefficients, format = Q31
    see bottom of file for table generation
    negative (vs. spec) since we use MADD for filter kernel
*/
static const int invQuant3[16] PROGMEM = {
    0x00000000, 0xc8767f65, 0x9becf22c, 0x83358feb, 0x83358feb, 0x9becf22c, 0xc8767f65, 0x00000000,
    0x2bc750e9, 0x5246dd49, 0x6ed9eba1, 0x7e0e2e32, 0x7e0e2e32, 0x6ed9eba1, 0x5246dd49, 0x2bc750e9,
};

static const int invQuant4[16] PROGMEM = {
    0x00000000, 0xe5632654, 0xcbf00dbe, 0xb4c373ee, 0xa0e0a15f, 0x9126145f, 0x8643c7b3, 0x80b381ac,
    0x7f7437ad, 0x7b1d1a49, 0x7294b5f2, 0x66256db2, 0x563ba8aa, 0x4362210e, 0x2e3d2abb, 0x17851aad,
};

#pragma GCC diagnostic pop

/**************************************************************************************
    Function:    DecodeLPCCoefs

    Description: decode LPC coefficients for TNS

    Inputs:      order of TNS filter
                resolution of coefficients (3 or 4 bits)
                coefficients unpacked from bitstream
                scratch buffer (b) of size >= order

    Outputs:     LPC coefficients in Q(FBITS_LPC_COEFS), in 'a'

    Return:      none

    Notes:       assumes no guard bits in input transform coefficients
                a[i] = Q(FBITS_LPC_COEFS), don't store a0 = 1.0
                  (so a[0] = first delay tap, etc.)
                max abs(a[i]) < log2(order), so for max order = 20 a[i] < 4.4
                  (up to 3 bits of gain) so a[i] has at least 31 - FBITS_LPC_COEFS - 3
                  guard bits
                to ensure no intermediate overflow in all-pole filter, set
                  FBITS_LPC_COEFS such that number of guard bits >= log2(max order)
 **************************************************************************************/
static int TNSMulShift32(int x, int y);
static int TNSQ20ToInt(long long v);
static int TNSClipQ20ToInt(long long v);

/* === m68k / Amiga optimised TNS helpers ===
 *
 * The 68030 MULS.L instruction produces a signed 32x32→64-bit product in
 * two data registers (Dh:Dl) in a single opcode, making it well-suited to
 * the fixed-point MAC loops in DecodeLPCCoefs and FilterRegion.
 *
 * These helpers compile only when AMIGA_M68K is defined.  The portable
 * long-long C path (TNSMulShift32 / TNSClipQ20ToInt below) remains the
 * fallback for all other targets.  Output values are bit-exact with the
 * portable path.
 *
 * Fixed-point conventions:
 *   a[] LPC coefficients : Q(FBITS_LPC_COEFS) = Q20, |a[i]| < 4.4 * 2^20
 *   hist[] samples       : integer (no fractional bits), fits in 32 bits
 *   sum accumulator      : Q20; right-shift by FBITS_LPC_COEFS to recover sample
 *
 * 68030 register roles used in the asm blocks below:
 *   Dl  "+d"  — read-write data register: multiplicand on entry, low 32 bits
 *               of the 64-bit product on exit  (same register, ISA requirement)
 *   Dh  "=d"  — write-only data register: high 32 bits of the 64-bit product
 *   Dm  "d"   — read-only  data register: the other multiplicand (unchanged)
 *
 * Carry chain: ADD.L sets the m68k X (extend) flag; ADDX.L consumes it for
 * correct 64-bit addition.  Both instructions live in one asm block so the
 * compiler cannot insert code between them that would corrupt X.
 *
 * No FPU, no 64-bit division, 68030-compatible.
 */
#ifdef AMIGA_M68K

/* tns_m68k_mulshift31 - signed 32x32 multiply, return bits [62:31]
 *
 * Computes ((int64)x * y) >> 31 in one MULS.L opcode plus two shifts.
 * This replaces the pattern TNSMulShift32(x, y) << 1 used in the
 * DecodeLPCCoefs Levinson-Durbin update loop, saving one instruction.
 *
 * Input  x : Q31  (invQuantTab entry, full 32-bit signed)
 * Input  y : Q(FBITS_LPC_COEFS)  (current LPC coefficient)
 * Output   : Q(FBITS_LPC_COEFS)  (updated coefficient, same format)
 */
static inline int tns_m68k_mulshift31(int x, int y)
{
    int hi;
    /* MULS.L Dm, Dh:Dl
     *   x  lives in Dl ("+d"): input multiplicand, receives low 32-bit result
     *   hi lives in Dh ("=d"): receives high 32-bit result
     *   y  lives in Dm ("d") : multiplier, left unchanged
     * Result: hi:x = signed 64-bit product of x * y
     */
    /* "=&d" early-clobber ensures GCC allocates hi to a register distinct
     * from x and y, satisfying the MULS.L encoding requirement Dh != Dl.
     */
    __asm__ (
        "muls.l %[y],%[hi]:%[x]"
        : [x] "+d"(x), [hi] "=&d"(hi)
        : [y] "d"(y)
    );
    /* (hi:x) >> 31 = bits [62:31] of the 64-bit product:
     *   hi << 1        extracts bits [62:32]  (upper half into position)
     *   (uint)x >> 31  extracts bit  [31]     (MSB of lower half)
     * Casting to unsigned avoids left-shifting a potentially negative int.
     */
    return (int)(((unsigned int)hi << 1) | ((unsigned int)x >> 31));
}

/* tns_m68k_mac - 64-bit signed accumulate: *sum_hi:*sum_lo += x * y
 *
 * Uses MULS.L for the 32x32→64 multiply, then ADD.L + ADDX.L to fold the
 * 64-bit product into the running accumulator with correct carry propagation.
 *
 * All three instructions are in one asm block so the compiler cannot insert
 * code between ADD.L (which sets X) and ADDX.L (which reads X).
 *
 * Inputs  x, y    : arbitrary signed 32-bit values
 * In/out *sum_hi  : high 32 bits of the 64-bit accumulator (signed)
 * In/out *sum_lo  : low  32 bits of the 64-bit accumulator (unsigned)
 */
static inline void tns_m68k_mac(int *sum_hi, unsigned int *sum_lo, int x, int y)
{
    int ph;
    /* "=&d" early-clobber on ph: MULS.L requires Dh != Dl, so ph must be
     * allocated to a register different from x (Dl), y (Dm), and the
     * accumulator halves.
     */
    __asm__ volatile (
        "muls.l %[y],%[ph]:%[x]\n\t"   /* ph:x = x*y  (Dl=x→lo, Dh=ph, Dm=y) */
        "add.l  %[x],%[slo]\n\t"        /* sum_lo += low product; sets X flag   */
        "addx.l %[ph],%[shi]"            /* sum_hi += high product + X           */
        : [x]   "+d"(x),    [ph]  "=&d"(ph),
          [slo] "+d"(*sum_lo), [shi] "+d"(*sum_hi)
        : [y] "d"(y)
        : "cc"
    );
}

#endif /* AMIGA_M68K */

static void DecodeLPCCoefs(int order, int res, signed char *filtCoef, int *a, int *b) {
    int i, m, t;
    const int *invQuantTab;

    if (res == 3) {
        invQuantTab = invQuant3;
    } else if (res == 4) {
        invQuantTab = invQuant4;
    } else {
        return;
    }

    for (m = 0; m < order; m++) {
        t = invQuantTab[filtCoef[m] & 0x0f];	/* t = Q31 */
        for (i = 0; i < m; i++) {
            /* Levinson-Durbin update: b[i] = a[i] - t * a[m-i-1] / 2^31
             * The Q31 coefficient t is multiplied with an a[] entry (Q20),
             * and the product is right-shifted by 31 to stay in Q20.
             * On m68k: one MULS.L + two shifts instead of two separate shifts.
             */
#ifdef AMIGA_M68K
            b[i] = a[i] - tns_m68k_mulshift31(t, a[m - i - 1]);
#else
            b[i] = a[i] - (TNSMulShift32(t, a[m - i - 1]) << 1);
#endif
        }
        for (i = 0; i < m; i++) {
            a[i] = b[i];
        }
        a[m] = t >> (31 - FBITS_LPC_COEFS);
    }
}

/**************************************************************************************
    Function:    FilterRegion

    Description: apply LPC filter to one region of coefficients

    Inputs:      number of transform coefficients in this region
                direction flag (forward = 1, backward = -1)
                order of filter
                'size' transform coefficients
                'order' LPC coefficients in Q(FBITS_LPC_COEFS)
                scratch buffer for history (must be >= order samples long)

    Outputs:     filtered transform coefficients

    Return:      guard bit mask (OR of abs value of all filtered transform coefs)

    Notes:       assumes no guard bits in input transform coefficients
                gains 0 int bits
                history buffer does not need to be preserved between regions
 **************************************************************************************/

/*
 * Amiga/m68k-safe TNS fixed-point helpers.
 *
 * Avoid endian-sensitive U64 hi32/lo32 field access and avoid left-shifting
 * negative signed values in the Q20 accumulator path.
 */
static int TNSMulShift32(int x, int y)
{
    return (int)(((long long)x * (long long)y) >> 32);
}

static int TNSQ20ToInt(long long v)
{
    unsigned long long u = (unsigned long long)v;
    return (int)(u >> FBITS_LPC_COEFS);
}

static int TNSClipQ20ToInt(long long v)
{
    long long hi = v >> 32;
    int y = TNSQ20ToInt(v);

    if ((hi >> 31) != (hi >> (FBITS_LPC_COEFS - 1))) {
        y = (int)((hi >> 31) ^ 0x7fffffff);
    }

    return y;
}

static int FilterRegion(int size, int dir, int order, int *audioCoef, int *a, int *hist) {
    int i, j, y, inc, gbMask;
#ifndef AMIGA_M68K
    long long sum64;
#endif

    /* init history to 0 every time */
    for (i = 0; i < order; i++) {
        hist[i] = 0;
    }

    gbMask = 0;
    inc = (dir ? -1 : 1);
    do {
        y = *audioCoef;

#ifdef AMIGA_M68K
        /* m68k path: explicit 64-bit accumulator as two 32-bit halves.
         *
         * Initialise sum_hi:sum_lo = (int64)y << FBITS_LPC_COEFS without a
         * 64-bit multiply.  Two arithmetic/logical shifts suffice:
         *   sum_hi = y >> (32 - FBITS_LPC_COEFS)   sign-extends the value
         *   sum_lo = (uint32)y << FBITS_LPC_COEFS   fills lower bits with 0
         * This is bit-exact with (long long)y * (1LL << FBITS_LPC_COEFS).
         */
        {
            int sum_hi = y >> (32 - FBITS_LPC_COEFS);
            unsigned int sum_lo = (unsigned int)y << FBITS_LPC_COEFS;
            int ph;

            /* Accumulate: sum_hi:sum_lo += hist[j] * a[j]  for j = order-1..0
             * History is shifted down (hist[j] = hist[j-1]) during the same pass.
             * Each MULS.L + ADD.L + ADDX.L triple is self-contained: ADD.L sets
             * the X flag and ADDX.L immediately consumes it in the same asm block.
             */
            for (j = order - 1; j > 0; j--) {
                int hj = hist[j];           /* read before the shift overwrites */
                hist[j] = hist[j - 1];     /* shift history toward older taps  */
                /* "=&d" on ph: MULS.L requires Dh != Dl, so ph must not share
                 * a register with hj (Dl), aj (Dm), or the accumulator halves.
                 */
                __asm__ volatile (
                    "muls.l %[aj],%[ph]:%[hj]\n\t" /* ph:hj = hist[j] * a[j]   */
                    "add.l  %[hj],%[slo]\n\t"       /* sum_lo += low;  sets X   */
                    "addx.l %[ph],%[shi]"            /* sum_hi += high + X       */
                    : [hj] "+d"(hj), [ph] "=&d"(ph),
                      [slo] "+d"(sum_lo), [shi] "+d"(sum_hi)
                    : [aj] "d"(a[j])
                    : "cc"
                );
            }
            /* j == 0: oldest tap — accumulate hist[0] without a history shift
             * (hist[0] is overwritten with y after the loop, below).
             */
            {
                int hj = hist[0];
                __asm__ volatile (
                    "muls.l %[aj],%[ph]:%[hj]\n\t"
                    "add.l  %[hj],%[slo]\n\t"
                    "addx.l %[ph],%[shi]"
                    : [hj] "+d"(hj), [ph] "=&d"(ph),
                      [slo] "+d"(sum_lo), [shi] "+d"(sum_hi)
                    : [aj] "d"(a[0])
                    : "cc"
                );
            }

            /* Extract y = (sum_hi:sum_lo) >> FBITS_LPC_COEFS with saturation.
             *
             * Overflow check: if the top (32 - FBITS_LPC_COEFS + 1) = 13 bits
             * of sum_hi are not all equal to the sign bit, the result does not
             * fit in a signed 32-bit integer and must be saturated.
             *
             *   sum_hi >> 31           : sign bit of the accumulator
             *   sum_hi >> (FBITS_LPC_COEFS - 1) : bits [31:19] uniform-sign test
             *
             * No-overflow extraction:
             *   y = (uint)sum_hi << (32 - FBITS_LPC_COEFS) | sum_lo >> FBITS_LPC_COEFS
             * Casting sum_hi to unsigned avoids UB from left-shifting a negative int.
             */
            if ((sum_hi >> 31) != (sum_hi >> (FBITS_LPC_COEFS - 1))) {
                y = (sum_hi >> 31) ^ 0x7fffffff;
            } else {
                y = (int)(((unsigned int)sum_hi << (32 - FBITS_LPC_COEFS)) |
                          (sum_lo >> FBITS_LPC_COEFS));
            }
        }
#else
        /* Portable path: 64-bit long long accumulator.
         * Use multiply instead of left-shifting a possibly negative int.
         */
        sum64 = (long long)y * (1LL << FBITS_LPC_COEFS);

        /* sum64 += (a1*y[n-1] + a2*y[n-2] + ... + a[order-1]*y[n-(order-1)]) */
        for (j = order - 1; j > 0; j--) {
            sum64 += (long long)hist[j] * (long long)a[j];
            hist[j] = hist[j - 1];
        }
        sum64 += (long long)hist[0] * (long long)a[0];

        y = TNSClipQ20ToInt(sum64);
#endif /* AMIGA_M68K */

        hist[0] = y;
        *audioCoef = y;
        audioCoef += inc;
        gbMask |= FASTABS(y);
    } while (--size);

    return gbMask;
}

/**************************************************************************************
    Function:    TNSFilter

    Description: apply temporal noise shaping, if enabled

    Inputs:      valid AACDecInfo struct
                index of current channel

    Outputs:     updated transform coefficients
                updated minimum guard bit count for this channel

    Return:      0 if successful, -1 if error
 **************************************************************************************/
int TNSFilter(AACDecInfo *aacDecInfo, int ch) {
#ifdef AMIGA_AAC_DISABLE_TNS_TEST
    (void)aacDecInfo;
    (void)ch;
    return 0;
#endif
    int win, winLen, nWindows, nSFB, filt, bottom, top, order, maxOrder, dir;
    int start, end, size, tnsMaxBand, numFilt, gbMask;
    int *audioCoef;
    unsigned char *filtLength, *filtOrder, *filtRes, *filtDir;
    signed char *filtCoef;
    const unsigned /*char*/ int *tnsMaxBandTab;
    const /*short*/ int *sfbTab;
    ICSInfo *icsInfo;
    TNSInfo *ti;
    PSInfoBase *psi;

    /* validate pointers */
    if (!aacDecInfo || !aacDecInfo->psInfoBase) {
        return -1;
    }
    psi = (PSInfoBase *)(aacDecInfo->psInfoBase);
    icsInfo = (ch == 1 && psi->commonWin == 1) ? &(psi->icsInfo[0]) : &(psi->icsInfo[ch]);
    ti = &psi->tnsInfo[ch];

    if (!ti->tnsDataPresent) {
        return 0;
    }

    if (icsInfo->winSequence == 2) {
        nWindows = NWINDOWS_SHORT;
        winLen = NSAMPS_SHORT;
        nSFB = sfBandTotalShort[psi->sampRateIdx];
        maxOrder = tnsMaxOrderShort[aacDecInfo->profile];
        sfbTab = sfBandTabShort + sfBandTabShortOffset[psi->sampRateIdx];
        tnsMaxBandTab = tnsMaxBandsShort + tnsMaxBandsShortOffset[aacDecInfo->profile];
        tnsMaxBand = tnsMaxBandTab[psi->sampRateIdx];
    } else {
        nWindows = NWINDOWS_LONG;
        winLen = NSAMPS_LONG;
        nSFB = sfBandTotalLong[psi->sampRateIdx];
        maxOrder = tnsMaxOrderLong[aacDecInfo->profile];
        sfbTab = sfBandTabLong + sfBandTabLongOffset[psi->sampRateIdx];
        tnsMaxBandTab = tnsMaxBandsLong + tnsMaxBandsLongOffset[aacDecInfo->profile];
        tnsMaxBand = tnsMaxBandTab[psi->sampRateIdx];
    }

    if (tnsMaxBand > icsInfo->maxSFB) {
        tnsMaxBand = icsInfo->maxSFB;
    }

    filtRes =    ti->coefRes;
    filtLength = ti->length;
    filtOrder =  ti->order;
    filtDir =    ti->dir;
    filtCoef =   ti->coef;

    gbMask = 0;
    audioCoef =  psi->coef[ch];
    for (win = 0; win < nWindows; win++) {
        bottom = nSFB;
        numFilt = ti->numFilt[win];
        for (filt = 0; filt < numFilt; filt++) {
            top = bottom;
            bottom = top - *filtLength++;
            bottom = MAX(bottom, 0);
            order = *filtOrder++;
            order = MIN(order, maxOrder);

            if (order) {
                start = sfbTab[MIN(bottom, tnsMaxBand)];
                end   = sfbTab[MIN(top, tnsMaxBand)];
                size = end - start;
                if (size > 0) {
                    dir = *filtDir++;
                    if (dir) {
                        start = end - 1;
                    }

                    DecodeLPCCoefs(order, filtRes[win], filtCoef, psi->tnsLPCBuf, psi->tnsWorkBuf);
                    gbMask |= FilterRegion(size, dir, order, audioCoef + start, psi->tnsLPCBuf, psi->tnsWorkBuf);
                }
                filtCoef += order;
            }
        }
        audioCoef += winLen;
    }

    /* update guard bit count if necessary */
    size = CLZ(gbMask) - 1;
    if (psi->gbCurrent[ch] > size) {
        psi->gbCurrent[ch] = size;
    }

    return 0;
}

/*	 Code to generate invQuantXXX[] tables
    {
      int res, i, t;
      double powScale, iqfac, iqfac_m, d;

      powScale = pow(2.0, 31) * -1.0;	/ ** make coefficients negative for using MADD in kernel ** /
      for (res = 3; res <= 4; res++) {
        iqfac =   ( ((1 << (res-1)) - 0.5) * (2.0 / M_PI) );
        iqfac_m = ( ((1 << (res-1)) + 0.5) * (2.0 / M_PI) );
        printf("static const int invQuant%d[16] = {\n", res);
        for (i = 0; i < 16; i++) {
          / ** extend bottom 4 bits into signed, 2's complement number ** /
          t = (i << 28) >> 28;

          if (t >= 0)   d = sin(t / iqfac);
          else          d = sin(t / iqfac_m);

          d *= powScale;
          printf("0x%08x, ", (int)(d > 0 ? d + 0.5 : d - 0.5));
          if ((i & 0x07) == 0x07)
             printf("\n");
        }
        printf("};\n\n");
 	  }
    }
*/

