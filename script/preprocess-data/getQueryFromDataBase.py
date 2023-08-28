import linecache
import random
filename = "/home/wangziyuan/bio/benchmark/mmseqs2-benchmark-pub/db/targetdb.fasta"
file = open(filename)
cnt = 0
res = ""
sum = 0
while sum < 1000000:
    title = file.readline().strip()
    text = file.readline().strip()
    if not title:
        break
    le = len(text)
    if le == 1000:
        cnt += 1
        sum += le
        print(cnt, sum)
        res += title + "\n" + text + "\n"

open("/home/wangziyuan/fxw/work/TASK03/query3/A.fasta", "w").write(res)

print("total", sum)
# >tr|C6E221|C6E221_GEOSM Transcriptional regulator, IclR family OX=443144 OS=Geobacter sp. (strain M21). GN= PE=4 SV=1
# YGLAASIEAAGQKLLPVLEQEIRERHFRQSPGLMSVAGIVRSTYDRIPAGVSTVGLNLEESEVAYGQRAIERLQDKWSEADTITHPTHKEFKNNSTYKAIAEDAFHAMFSKGPATCYAPLQVGLAPMVRLNHDCEVSDLNIVNFDSMISVCTTENCERVLREMVTRAGGLGRQGLAKQALEVTMFGLQYRETVPDKSVYNRLQLTALLRFINNKPLKLRRSLETLGIEHDDRLFQELVDLARSVSRNYDTTKKVTKM
