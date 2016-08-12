import sys
from collections import namedtuple, defaultdict
Sequence = namedtuple('Sequence', 'label sequence')

def iterate_fasta(input_file):
    sequences = []
    input_file.seek(0)
    ignore_lines = False
    for line in input_file:
        line = line.strip()
        if line[0] == ">": # FASTA Label
            sequences.append(Sequence(label=line, sequence=[]))
        elif line[0] == "+": # FASTQ beginning of quality info
            ignore_lines = True
        elif line[0] == "@": # FASTQ Label
            sequences.append(Sequence(label=line, sequence=[]))
            ignore_lines = False # We are at a new sequence
        elif not ignore_lines:
            sequences[-1].sequence.extend(list(line))
    return sequences

fasta = iterate_fasta(open("CAP256_4180_082wpi_v1v3_5t.fasta"))
pids = iterate_fasta(open("CAP256_4180_082wpi_v1v3_pid_found.fasta"))
nopids = iterate_fasta(open("CAP256_4180_082wpi_v1v3_pid_not_found.fasta"))

print(len(fasta), file=sys.stderr)
print("%d + %d = %d" % (len(pids), len(nopids), len(pids) + len(nopids)), file=sys.stderr)

pid_s = []
for p in pids:
    pid_s.append("".join(p.sequence))
nid_s = []
for n in nopids:
    nid_s.append("".join(n.sequence))

pid_count = [0] * len(pids)

p_i = 0
n_i = 0
for ref in fasta:
    #print(ref.label)
    r = "".join(ref.sequence)
    for n in nid_s:
        if n == r:
            #print("no")
            break
    else:
        for p_i, p in enumerate(pid_s):
            if p in r:
                pid_count[p_i] += 1
                if pid_count[p_i] > 1:
                    print(pid_count[p_i])
                #print("pid")
                break
        else:
            print("---")
print(pid_count)
