import sys
import bisect
import multiprocessing
from multiprocessing import Pool

class Record(object):
    def __init__(self, chrom, start, end, mC):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.mC = mC

bins = int(sys.argv[3])

def split_integer(m, n):
    assert n > 0
    quotient = int(m / n)
    remainder = m % n
    if remainder > 0:
        return [quotient] * (n - remainder) + [quotient + 1] * remainder
    if remainder < 0:
        return [quotient - 1] * -remainder + [quotient] * (n + remainder)
    return [quotient] * n

all_bdgs = {}
with open(sys.argv[1]) as infilebdg:
    for line in infilebdg:
        line = line.strip()
        line_s = line.split('\t')
        locus = str(int(int(line_s[1])/1000)) + '-' + str(int(int(line_s[2])/1000))
        record1 = Record(line_s[0],int(line_s[1]),int(line_s[2]),float(line_s[3]))
        if not line_s[0] in all_bdgs:
            all_bdgs[line_s[0]] = {}
        if not locus in all_bdgs[line_s[0]]:
            all_bdgs[line_s[0]][locus] = []
        all_bdgs[line_s[0]][locus].append(record1)

#print(all_bdgs)

def by_position(t):
    return t.start


def multi(line):
    line = line.strip()
    line_s = line.split('\t')
    if line_s[0] in all_bdgs:
        #print(111)
        start = int(line_s[1])
        end = int(line_s[2])
        strand = line_s[5]
        locus_start = int(start/1000) - 1
        locus_end = int(end/1000) + 1
        Records = []
        for x in range(locus_start,locus_end):
            for y in range(x,locus_end):
                if str(x)+'-'+str(y) in all_bdgs[line_s[0]]:
                    Records.extend(all_bdgs[line_s[0]][str(x)+'-'+str(y)])
        #print(Records)
        #strand = "+"
        if end-start > bins:
            size = split_integer(end-start, bins)
            start1 = start
            bed_split_bins = []
            for i in size:
                bed_split_bins.append((start1,start1+i))
                start1 = start1+i
            #print(bed_split_bins)
            Record_sort = sorted(Records,key=by_position)
            starts_sort = [record_temp.start for record_temp in Record_sort]
            strs_out = ""
            allC_total = 0
            for start_bed,end_bed in bed_split_bins:
                insert_left = bisect.bisect_left(starts_sort,start_bed)
                insert_right = bisect.bisect_right(starts_sort,end_bed)
                if insert_right == insert_left:
                    score = 'NA'
                elif insert_right > insert_left:
                    Record_used = Record_sort[insert_left:insert_right]
                    mC_total = 0
                    allC_total = 0
                    for record_temp in Record_used:
                        mC_total += record_temp.mC
                        allC_total += 1
                    score = float(mC_total)/allC_total
                strs_out = strs_out + str(score) + "\t"
            strs_out = strs_out.strip()
            if strand == '-':
                strs_out_s = strs_out.split('\t')
                strs_out_s.reverse()
                strs_out = '\t'.join(strs_out_s)
            results[line] = (strs_out+'\t'+str(allC_total))
        else:
            results[line] = ("shorter than bin number")
    else:
        results[line] = ("chrom not in bdg file")

reads = []
with open(sys.argv[2]) as infilebed:
    for line in infilebed:
        reads.append(line)

with multiprocessing.Manager() as MG:
    results = multiprocessing.Manager().dict()
    pool = Pool(20)
    for line in reads:
        pool.apply_async(multi, args=(line,))
    pool.close()
    pool.join()

for line, result in results.items():
    print(line+'\t'+result)


