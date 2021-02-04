"""Generates mutation sequences and writes them to a .fasta file along with the
mut_escape values for each antibody, plus other provided metadata (site, reference,
mutation).
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REFERENCE_SEQUENCE_FILE = 'escape_mutations/data/cov2_spike_wt.fasta'
AGGREGATED_MUT_ESCAPES_FILE = 'escape_mutations/data/aggregated_mut_escapes.csv'

OUTPUT_FASTA_FILE = 'escape_mutations/data/aggregated_mut_escapes.fasta'

def read_reference_sequence():
    for record in SeqIO.parse(REFERENCE_SEQUENCE_FILE, 'fasta'):
        return record.seq


def main():
    reference_sequence = list(read_reference_sequence())
    mut_escapes_df = pd.read_csv(AGGREGATED_MUT_ESCAPES_FILE)

    # Add header to the sequence records.
    records = [
        SeqRecord(
            Seq('XXXXXXXX'),
            id='id',
            name='',
            description='|'.join(mut_escapes_df.columns)
        )
    ]

    for index, row in mut_escapes_df.iterrows():
        mutation_index = row['site'] - 1

        assert reference_sequence[mutation_index] == row['wildtype']

        sequence = reference_sequence.copy()
        sequence[mutation_index] = row['mutation']

        description = '|'.join(row.astype(str).values.flatten().tolist())

        records.append(
            SeqRecord(
                Seq(''.join(sequence)),
                id=str(index),
                name='',
                description=description)
        )

    SeqIO.write(records, OUTPUT_FASTA_FILE, 'fasta')


if __name__ == '__main__':
    main()
