/*
 * bn256-20080525/parameters.h 
 * Peter Schwabe
 * Public domain
*/

#include <gmp.h>

// Parameters used in fp.c
#define BN_P "21888242871839275222246405745257275088696311157297823662689037894645226208583"
#define BN_PINV32 3834012553ULL
#define BN_PINV64 9786893198990664585ULL

// Parameters used in fpe.c

#if(GMP_LIMB_BITS == 32)
#define N_LIMBS 8
#elif(GMP_LIMB_BITS == 64)
#define N_LIMBS 4
#else
#error "Only 32 and 64 bit architectures are supported"
#endif

// Parameters used in fp2.c
#define ALPHA (-1) // constant coefficient in the irreducible polynomial x^2 - alpha, used to construct F_{p^2}

// Parameters used in fp6.c
#define BN_XI "1", "9"
#define BN_YPMINUS1 "10307601595873709700152284273816112264069230130616436755625194854815875713954", "21575463638280843010398324269430826099269044274347216827212613867836435027261"
#define BN_ZETA "2203960485148121921418603742825762020974279258880205651966"

#define BN_XI2 "2", "82434016654300679721217353503190038836571781811386228921167322412819029493182"
#define BN_1O27XI3 "33584229007307684330866329205003349155640355552786982153068168390407752756482", "79380904926363617509320414484553370731513567670223775998161125286418324697139"
#define BN_1O3XI3 "54956011102867119814144902335460025891047854540924152614111548275212686328789", "54956011102867119814144902335460025891047854540924152614111548275212686328787"
#define BN_1O3XI "54956011102867119814144902335460025891047854540924152614111548275212686328789", "54956011102867119814144902335460025891047854540924152614111548275212686328789"
#define BN_1O3MODP "54956011102867119814144902335460025891047854540924152614111548275212686328789"
#define BN_COMETA_C0_CONST "82434016654300679717245125061265641166427769693017678351449950981238451124296"
#define BN_COMETA_C1_CONST "82434016654300679719231239282227840001499775752201953636308636697028740308739"

// Parameters used in fp12.c
#define BN_TAU "0", "0", "0", "1", "0", "0" // constant tau used to construct F_p^12 as F_p^6[Z]/ (Z^2 - tau)
#define BN_ZPMINUS1 "16469823323077808223889137241176536799009286646108169935659301613961712198316", "8376118865763821496583973867626364092589906065868298776909617916018768340080" // Z^(p-1)
#define BN_ZPMINUS1INV "5722266937896532885780051958958348231143373700109372999374820235121374419868", "18566938241244942414004596690298913868373833782006617400804628704885040364344" // Z^(1-p)

// Parameters used in curve.c
#define BN_X "4965661367192848881" // parameter x used to generate the curve (see "Pairing-Friendly Elliptic Curves of Prime Order")
#define BN_N "21888242871839275222246405745257275088548364400416034343698204186575808495617" // prime order of E(F_p)
#define BN_TRACE "147946756881789318990833708069417712967" // trace of Frobenius of the curve
#define BN_CHI "552484233613224096312617126783173147097382103762957654188882734314196910839907541213974502761540629817009608548654680343627701153829446747810907373256841551006201639677726139946029199968412598804882391702273019083653272047566316584365559776493027495458238373902875937659943504873220554161550525926302303331747463515644711876653177129578303191095900909191624817826566688241804408081892785725967931714097716709526092261278071952560171111444072049229123565057483750161460024353346284167282452756217662335528813519139808291170539072125381230815729071544861602750936964829313608137325426383735122175229541155376346436093930287402089517426973178917569713384748081827255472576937471496195752727188261435633271238710131736096299798168852925540549342330775279877006784354801422249722573783561685179618816480037695005515426162362431072245638324744480" // exponent of final exponentiation
#define BN_CHICOMP "4965661367192848881" // exponent of final exponentiation for compressed pairing
#define BN_LOOPLENGTH_ETA "11916685325773193009570696613837024235910666865622157670373"
#define BN_B "3" // parameter b in the curve equation x^2 = y^3 + b

// Parameters used in points.c
#define BN_CURVEGEN "1", "2", "1"
#define BN_TWISTGEN_X "11559732032986387107991004021392285783925812861821192530917403151452391805634", "10857046999023057135944570762232829481370756359578518086990519993285655852781"
#define BN_TWISTGEN_Y "4082367875863433681332203403145435568316851327593401208105741076214120093531", "8495653923123431417604973247489272438418190587263600148770280649306958101930"

// Parameters used for OptAte computation in ate_optate.c
#define BN_ZETA2 "21888242871839275220042445260109153167277707414472061641714758635765020556616" // zeta^2
#define BN_Z2P "10307601595873709700152284273816112264069230130616436755625194854815875713954", "21575463638280843010398324269430826099269044274347216827212613867836435027261" // Z^(2p)
#define BN_Z3P "3505843767911556378687030309984248845540243509899259641013678093033130930403", "2821565182194536844548159561693502659359617185244120367078079554186484126554" // Z^(3p)
