import tensorflow as tf
import pandas as pd

from lib.mutation import *
from lib.sequence_parser_utils import (parse_aggregated_mut_escapes,
                                       parse_gisaid, parse_nih, parse_viprbrc)

np.random.seed(1)
random.seed(1)
AAs = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H',
    'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W',
    'Y', 'V', 'X', 'Z', 'J', 'U', 'B',
]

TRAIN_TEST_FNAMES = ['data/cov/sars_cov2_seqs.fa',
                     'data/cov/viprbrc_db.fasta',
                     'data/cov/gisaid.fasta']

AGGREGATED_MUT_ESCAPES_FNAMES = ['escape_mutations/data/aggregated_mut_escapes.fasta']
AGGREGATED_MUT_ESCAPES_EMBEDDINGS_OUTPUT_FNAME = 'escape_mutations/data/mut_escape_embeddings.csv'


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Coronavirus sequence analysis')
    parser.add_argument('model_name', type=str,
                        help='Type of language model (e.g., hmm, lstm)')
    parser.add_argument('--namespace', type=str, default='cov',
                        help='Model namespace')
    parser.add_argument('--dim', type=int, default=512,
                        help='Embedding dimension')
    parser.add_argument('--batch-size', type=int, default=500,
                        help='Training and inference minibatch size')
    parser.add_argument('--n-epochs', type=int, default=11,
                        help='Number of training epochs')
    parser.add_argument('--seed', type=int, default=1,
                        help='Random seed')
    parser.add_argument('--checkpoint', type=str, default=None,
                        help='Model checkpoint')
    parser.add_argument('--train', action='store_true',
                        help='Train model')
    parser.add_argument('--train-split', action='store_true',
                        help='Train model on portion of data')
    parser.add_argument('--test', action='store_true',
                        help='Test model')
    parser.add_argument('--embed', action='store_true',
                        help='Analyze embeddings')
    parser.add_argument('--semantics', action='store_true',
                        help='Analyze mutational semantic change')
    parser.add_argument('--combfit', action='store_true',
                        help='Analyze combinatorial fitness')
    parser.add_argument('--reinfection', action='store_true',
                        help='Analyze reinfection cases')
    parser.add_argument('--mut_escapes', action='store_true',
                        help='Analyze mutation escapes')
    args = parser.parse_args()
    return args

def process(fnames):
    parse_fn_by_fname = {
        'data/cov/sars_cov2_seqs.fa': parse_nih,
        'data/cov/viprbrc_db.fasta': parse_viprbrc,
        'data/cov/gisaid.fasta': parse_gisaid,
        'escape_mutations/data/aggregated_mut_escapes.fasta': parse_aggregated_mut_escapes,
        'escape_mutations/data/aggregated_mut_escapes_test.fasta': parse_aggregated_mut_escapes,
    }

    seqs = {}
    for fname in fnames:
        for record in SeqIO.parse(fname, 'fasta'):
            if len(record.seq) < 1000:
                continue
            if str(record.seq).count('X') > 0:
                continue
            if record.seq not in seqs:
                seqs[record.seq] = []

            meta = parse_fn_by_fname[fname](record)

            meta['accession'] = record.description
            seqs[record.seq].append(meta)

    with open('data/cov/cov_all.fa', 'w') as of:
        for seq in seqs:
            metas = seqs[seq]
            for meta in metas:
                of.write('>{}\n'.format(meta['accession']))
                of.write('{}\n'.format(str(seq)))

    tprint('Total different seqs read: {}'.format(len(seqs)))

    return seqs

def split_seqs(seqs, split_method='random'):
    train_seqs, test_seqs = {}, {}

    tprint('Splitting seqs...')
    for idx, seq in enumerate(seqs):
        if idx % 10 < 2:
            test_seqs[seq] = seqs[seq]
        else:
            train_seqs[seq] = seqs[seq]
    tprint('{} train seqs, {} test seqs.'
           .format(len(train_seqs), len(test_seqs)))

    return train_seqs, test_seqs

def setup(args):
    fnames = AGGREGATED_MUT_ESCAPES_FNAMES if args.mut_escapes else TRAIN_TEST_FNAMES
    seqs = process(fnames)

    seq_len = max([ len(seq) for seq in seqs ]) + 2
    vocab_size = len(AAs) + 2

    model = get_model(args, seq_len, vocab_size,
                      inference_batch_size=args.batch_size)

    return model, seqs

def interpret_clusters(adata):
    clusters = sorted(set(adata.obs['louvain']))
    for cluster in clusters:
        tprint('Cluster {}'.format(cluster))
        adata_cluster = adata[adata.obs['louvain'] == cluster]
        for var in [ 'host', 'country', 'strain' ]:
            tprint('\t{}:'.format(var))
            counts = Counter(adata_cluster.obs[var])
            for val, count in counts.most_common():
                tprint('\t\t{}: {}'.format(val, count))
        tprint('')

def plot_umap(adata, categories, namespace='cov'):
    for category in categories:
        sc.pl.umap(adata, color=category,
                   save='_{}_{}.png'.format(namespace, category))

def analyze_embedding(args, model, seqs, vocabulary):
    seqs = embed_seqs(args, model, seqs, vocabulary, use_cache=False)

    X, obs = [], {}
    obs['n_seq'] = []
    obs['seq'] = []
    for seq in seqs:
        meta = seqs[seq][0]
        X.append(meta['embedding'].mean(0))
        for key in meta:
            if key == 'embedding':
                continue
            if key not in obs:
                obs[key] = []
            obs[key].append(Counter([
                meta[key] for meta in seqs[seq]
            ]).most_common(1)[0][0])
        obs['n_seq'].append(len(seqs[seq]))
        obs['seq'].append(str(seq))
    X = np.array(X)

    adata = AnnData(X)
    for key in obs:
        adata.obs[key] = obs[key]

    sc.pp.neighbors(adata, n_neighbors=20, use_rep='X')
    sc.tl.louvain(adata, resolution=1.)
    sc.tl.umap(adata, min_dist=1.)

    sc.set_figure_params(dpi_save=500)
    plot_umap(adata, [ 'host', 'group', 'continent', 'louvain' ])

    interpret_clusters(adata)

    adata_cov2 = adata[(adata.obs['louvain'] == '0') |
                       (adata.obs['louvain'] == '2')]
    plot_umap(adata_cov2, [ 'host', 'group', 'country' ],
              namespace='cov7')


def create_mut_escape_embeddings(args, model, seqs, vocabulary):
    seqs = embed_seqs(args, model, seqs, vocabulary, use_cache=True)

    data = []
    for seq in seqs:
        meta = seqs[seq][0]
        row = [meta['id']]
        row.extend(meta['embedding'].mean(0))
        data.append(row)

    header = ['id']
    for i in range(len(data[0]) - 1):
        header.append('x_{}'.format(i))

    df = pd.DataFrame(data, columns=header)
    df.to_csv(
        AGGREGATED_MUT_ESCAPES_EMBEDDINGS_OUTPUT_FNAME,
        index=False,
        header=True
    )


def main():
    args = parse_args()

    vocabulary = { aa: idx + 1 for idx, aa in enumerate(sorted(AAs)) }

    model, seqs = setup(args)

    if 'esm' in args.model_name:
        args.checkpoint = args.model_name
    elif args.checkpoint is not None:
        model.model_.load_weights(args.checkpoint)
        tprint('Model summary:')
        tprint(model.model_.summary())

    if args.train:
        batch_train(args, model, seqs, vocabulary, batch_size=1000)

    if args.train_split or args.test:
        train_test(args, model, seqs, vocabulary, split_seqs)

    if args.embed:
        if args.checkpoint is None and not args.train:
            raise ValueError('Model must be trained or loaded '
                             'from checkpoint.')
        no_embed = { 'hmm' }
        if args.model_name in no_embed:
            raise ValueError('Embeddings not available for models: {}'
                             .format(', '.join(no_embed)))
        analyze_embedding(args, model, seqs, vocabulary)

    if args.mut_escapes:
        if args.checkpoint is None and not args.train:
            raise ValueError('Model must be trained or loaded '
                             'from checkpoint.')
        no_embed = { 'hmm' }
        if args.model_name in no_embed:
            raise ValueError('Embeddings not available for models: {}'
                             .format(', '.join(no_embed)))
        create_mut_escape_embeddings(args, model, seqs, vocabulary)


    if args.semantics:
        if args.checkpoint is None and not args.train:
            raise ValueError('Model must be trained or loaded '
                             'from checkpoint.')

        from escape import load_baum2020, load_greaney2020
        tprint('Baum et al. 2020...')
        seq_to_mutate, seqs_escape = load_baum2020()
        analyze_semantics(args, model, vocabulary,
                          seq_to_mutate, seqs_escape, comb_batch=10000,
                          prob_cutoff=0, beta=1., plot_acquisition=True,)
        tprint('Greaney et al. 2020...')
        seq_to_mutate, seqs_escape = load_greaney2020()
        analyze_semantics(args, model, vocabulary,
                          seq_to_mutate, seqs_escape, comb_batch=10000,
                          min_pos=318, max_pos=540, # Restrict to RBD.
                          prob_cutoff=0, beta=1., plot_acquisition=True,
                          plot_namespace='cov2rbd')

    if args.combfit:
        from combinatorial_fitness import load_starr2020
        tprint('Starr et al. 2020...')
        wt_seqs, seqs_fitness = load_starr2020()
        strains = sorted(wt_seqs.keys())
        for strain in strains:
            analyze_comb_fitness(args, model, vocabulary,
                                 strain, wt_seqs[strain], seqs_fitness,
                                 comb_batch=10000, prob_cutoff=0., beta=1.)

    if args.reinfection:
        from plot_reinfection import plot_reinfection
        from reinfection import load_ratg13, load_sarscov1, load_to2020

        tprint('To et al. 2020...')
        wt_seq, mutants = load_to2020()
        analyze_reinfection(args, model, seqs, vocabulary, wt_seq, mutants,
                            namespace='to2020')
        plot_reinfection(namespace='to2020')
        null_combinatorial_fitness(args, model, seqs, vocabulary,
                                   wt_seq, mutants, n_permutations=100000000,
                                   namespace='to2020')

        tprint('Positive controls...')
        wt_seq, mutants = load_ratg13()
        analyze_reinfection(args, model, seqs, vocabulary, wt_seq, mutants,
                            namespace='ratg13')
        plot_reinfection(namespace='ratg13')
        wt_seq, mutants = load_sarscov1()
        analyze_reinfection(args, model, seqs, vocabulary, wt_seq, mutants,
                            namespace='sarscov1')
        plot_reinfection(namespace='sarscov1')


if __name__ == '__main__':
    gpus = tf.config.experimental.list_physical_devices('GPU')
    if gpus:
        try:
            for gpu in gpus:
                tf.config.experimental.set_memory_growth(gpu, True)
        except RuntimeError as e:
            print(e)

    strategy = tf.distribute.MirroredStrategy()

    options = tf.data.Options()
    options.experimental_distribute.auto_shard_policy = tf.data.experimental.AutoShardPolicy.DATA

    with strategy.scope():
        main()
