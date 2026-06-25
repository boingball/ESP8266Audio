/**************************************************************************************
    Amiga / m68k optimisation helpers for the Helix AAC decoder.

    This header provides small, bit-exact inline-assembly helpers for the
    Motorola 68k family (68020 and above, which have the 32x32->64 MULS.L /
    MULU.L instructions and tolerate misaligned long accesses).

    Design rules (see project task description):
      * Every asm optimisation has a C fallback.
      * Every asm optimisation is behind an explicit per-stage flag
        (AMIGA_M68K_ASM_AAC_IMDCT / _DEQUANT / _STEREO / _HUFFMAN). The flags
        only *enable* use of these helpers; the helpers themselves are only
        compiled when we are actually building for a 68020+ target, so it is
        always safe to define the flags on any platform.
      * The helpers are bit-exact replacements for the portable C code, so
        decoder output is unchanged (no TNS / format / ABI impact).

    Nothing in this header changes any struct layout, symbol, or calling
    convention, so it does not affect the decoder module ABI or the position
    of DecoderModuleEntry.
 **************************************************************************************/

#ifndef _AMIGA_M68K_AAC_H
#define _AMIGA_M68K_AAC_H

/*  Detect a 68020+ GNU C target. Only these cores have the long (32x32->64)
    MULS.L form and handle misaligned 32-bit memory accesses in hardware, so
    the helpers below are restricted to them. Plain 68000/68010 fall back to
    the portable C paths in assembly.h / bitstream.c.
*/
#if defined(__GNUC__) && ( \
        defined(__mc68020__) || defined(__mc68030__) || \
        defined(__mc68040__) || defined(__mc68060__) || \
        defined(mc68020) || defined(mc68030) || \
        defined(mc68040) || defined(mc68060) )

#define AAC_M68K_HAVE_ASM 1

/*  AAC_M68K_MULSHIFT32(x, y)

    Signed 32x32 multiply returning the top 32 bits of the 64-bit product.
    Bit-exact equivalent of the portable C:

        (int)(((long long)x * (long long)y) >> 32)

    MULS.L <ea>,Dh:Dl computes Dh:Dl = Dl * <ea> (signed, full 64-bit result,
    Dh = high word). We seed Dl with x, multiply by y, and return Dh.
*/
static __inline int AAC_M68K_MULSHIFT32(int x, int y)
{
    int hi, lo;

    __asm__ ("muls.l %3,%0:%1"
             : "=d" (hi), "=d" (lo)
             : "1"  (x),  "d"  (y));
    (void)lo;					/* low 32 bits intentionally discarded */
    return hi;
}

/*  AAC_M68K_MADD64(sum, x, y)

    64-bit multiply-accumulate: sum += (long long)x * (long long)y, computed
    with a single MULS.L plus ADD/ADDX. Big-endian safe (the union is only used
    to name the high/low halves for the carry-propagating add). Bit-exact with
    the portable C 'sum += (Word64)x * (Word64)y'.
*/
static __inline long long AAC_M68K_MADD64(long long sum, int x, int y)
{
    union { long long w; struct { int hi; unsigned int lo; } r; } u;
    int phi;
    unsigned int plo;

    u.w = sum;
    __asm__ ("muls.l %3,%0:%1"
             : "=d" (phi), "=d" (plo)
             : "1"  (x),   "d"  (y));
    __asm__ ("add.l  %3,%1\n\t"
             "addx.l %2,%0"
             : "+d" (u.r.hi), "+d" (u.r.lo)
             : "d"  (phi),    "d"  (plo));
    return u.w;
}

/*  AAC_M68K_LOAD_BE32(p)

    Load 4 consecutive bytes as a big-endian 32-bit word, i.e.

        p[0]<<24 | p[1]<<16 | p[2]<<8 | p[3]

    On big-endian 68020+ this is a single (possibly misaligned) MOVE.L. Using
    inline asm avoids any strict-aliasing concern from casting the byte
    pointer. The caller guarantees at least 4 readable bytes.
*/
static __inline unsigned int AAC_M68K_LOAD_BE32(const unsigned char *p)
{
    unsigned int v;

    __asm__ ("move.l (%1),%0"
             : "=d" (v)
             : "a"  (p));
    return v;
}

/*  Misaligned long loads are only safe (and the load helper only meaningful)
    on the same 68020+ big-endian targets. */
#define AAC_M68K_BIG_ENDIAN_LOAD 1

#endif	/* 68020+ GNU C */

#endif	/* _AMIGA_M68K_AAC_H */
