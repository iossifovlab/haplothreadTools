#!/usr/bin/env python
import sys

if len(sys.argv) < 4:
    print("Usage: createBim.py <map file> <chip summary file>")
    exit(1)

mapFn = sys.argv[1]
allelesFn = sys.argv[2]
outFn = sys.argv[3]
print(mapFn, file=sys.stderr)
print(allelesFn, file=sys.stderr)
print(outFn, file=sys.stderr)

n = 0
out=open(outFn, 'w')
MF = open(mapFn)
AF = open(allelesFn)

#chrom	pos	refA	sscAlleles	majorA	minorA	MM_N	Mm_N	mm_N00_N	majorA_N	minorA_N	0_N	MM_parents	Mm_parents	mm_parents	OO_parents	majorA_parents	minorA_parents	O_parents	percentCalled	MAF	HW_p	inSSC	flip

### new
#chrom	varId	pos	refA	sscAlleles	agreAlleles	majorA	minorA	MM_N	Mm_N	mm_N	00_N	majorA_N	minorA_N	0_N	MM_dad	Mm_dad	mm_dad	OO_dad	MM_mom	Mm_mom	mm_mom	OO_mom	majorA_parents	minorA_parents	O_parents	percentCalled	MAF	HW_p	inSSC	inAGRE	flip	reject	reason


HD = {v:i for i,v in 
      enumerate(AF.readline().strip('\n\rq').split('\t'))}

compM = {x:y for x,y in zip(list('ACGT0ID'),list('TGCA0ID'))}
chroms = ['chr' + str(i) for i in range(1,24)]

while True:
    ml = MF.readline()
    al = AF.readline()
    assert (ml and al) or ((not ml) and (not al))

    if not ml: break

    al = al.strip("\n\r").split("\t")
    als = [al[HD['majorA']], al[HD['minorA']]]
    assert len(als) < 3
    chipA = set(''.join(als))
    chipAcomp = set([compM[x] for x in ''.join(als)])
    if als[1]=="":
        als[1] = '0'
    als.reverse()

    cs = ml.strip("\n\r").split('\t',-1)
    assert len(cs) == 5
    ch,posId,posCMS,posS,refA = cs

    ch = 'chr' + ch
    pos = int(posS)

    flip = al[HD['flip']]
    if al[HD['reject']] == 'yes':
        pos = -1
    
    out.write('\t'.join(map(str,[ch, posId, posCMS, pos] 
                            + als + [flip,refA]))+'\n')
out.close()
MF.close()
AF.close()

