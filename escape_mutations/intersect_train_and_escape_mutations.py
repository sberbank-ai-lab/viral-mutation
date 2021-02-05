from Bio import SeqIO

RBD_LENGTH = 1273
TRAIN_FILE = 'data/cov/cov_all.fa'
MUT_ESCAPE_FILE = 'escape_mutations/data/aggregated_mut_escapes.fasta'

def main():
    train_rbd = set()
    train_other = set()
    for record in SeqIO.parse(TRAIN_FILE, 'fasta'):
        if len(record.seq) < 1000:
            continue
        if str(record.seq).count('X') > 0:
            continue

        if len(record.seq) == RBD_LENGTH:
            train_rbd.add(record.seq)
        else:
            train_other.add(record.seq)

    count = 0
    for record in SeqIO.parse(MUT_ESCAPE_FILE, 'fasta'):
        if len(record.seq) < 1000:
            continue
        if str(record.seq).count('X') > 0:
            continue

        if record.seq in train_rbd:
            count += 1

        for seq in train_other:
            if seq.find(record.seq) >= 0:
                count += 1
                break

    print(count, len(train_rbd), len(train_rbd) + len(train_other))


if __name__ == '__main__':
    main()
