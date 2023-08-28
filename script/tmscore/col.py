import subprocess

res = subprocess.run(
    "ls /home/wangziyuan/bio/af2/models/colab".split(" "), capture_output=True, text=True)

ls = list(filter(lambda x: x != "", res.stdout.split("\n")))
output_colab = dict()
for i in ls:
    print(i)
    ans = subprocess.run(
        f"/home/wangziyuan/bio/af2/TMscore_cpp  /home/wangziyuan/bio/af2/models/colab/{i} /home/wangziyuan/bio/af2/models/real/{i}".split(" "), capture_output=True, text=True)
    # print("error", ans.stderr)
    ans = ans.stdout.split("\n")
    for j in ans:
        if "TM-score" in j and "TM-score" == j[:8]:
            k = j.split("=")[1].split(" ")[1]
            output_colab[i.split(".")[0]] = k


output_alphafold = dict()
res = subprocess.run(
    "ls /home/wangziyuan/bio/af2/alphafold2_res/caps14".split(" "), capture_output=True, text=True)
ls = list(filter(lambda x: x != "", res.stdout.split("\n")))
for i in ls:
    print(i)
    ans = subprocess.run(
        f"/home/wangziyuan/bio/af2/TMscore_cpp /home/wangziyuan/bio/af2/alphafold2_res/caps14/{i} /home/wangziyuan/bio/af2/models/real/{i}".split(" "), capture_output=True, text=True)
    ans = ans.stdout.split("\n")
    for j in ans:
        if "TM-score" in j and "TM-score" == j[:8]:
            k = j.split("=")[1].split(" ")[1]
            output_alphafold[i.split(".")[0]] = k


F = open("/home/wangziyuan/fxw/MyTest/af/tmscore/res/go.log", "w")
F2 = open("/home/wangziyuan/fxw/MyTest/af/tmscore/res/compare.log", "r").readlines()

F.write("name\ttantan\tno_tantan\t\tcolab\talphafold\n")
for i in F2:
    j = [l.strip() for l in i.split("\t")]
    if j[0] in output_colab and j[0] in output_alphafold:
        j.append("\t" + output_colab[j[0]] + "\t" + output_alphafold[j[0]])
        F.write("\t".join(j) + "\n")
