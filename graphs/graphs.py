from collections import namedtuple, defaultdict
try:
    from matplotlib import pyplot as plt
except ImportError:
    print("Matplotlib does not seem to be available. Try 'pip install matplotlib'\nError:")
    raise
import argparse

def tag_dist(input_file):
    tags = defaultdict(lambda: 0)
    for line in input_file:
        l = line.strip()
        if len(l) == 0:
            continue
        tags[l] = tags[l] + 1
    max_count = 0
    count_dist = defaultdict(lambda: 0)
    for count in tags.values():
        max_count = max(max_count, count)
        count_dist[count] = count_dist[count] + 1
    count_counts = [count_dist[i] for i in range(1, max_count + 1)]
    bars = plt.bar(range(max_count), count_counts, width=1.0, linewidth=0, log=True)
    plt.xlabel('Number of copies in bin')
    plt.ylabel('Unique ID\'s with bin size')
    plt.xticks([0.5 + x for x in range(max_count)], range(1, max_count+1))
    plt.savefig('bin_sizes')
    #plt.show()
    plt.cla()

def likelihood_cutoffs(input_file):
    num_bins = 20

    templates = defaultdict(list)
    for line in input_file:
        parts = line.strip().split()
        if len(parts) <= 1: continue
        template = parts[0]
        score = float(parts[1])
        templates[template].append(score)
    for template, scores in sorted(templates.items()):
        n, bins, patches = plt.hist(scores, num_bins)
        plt.xlabel('Likelihood')
        plt.ylabel('Frequency')
        plt.title(template)
        plt.savefig(template)
        #plt.show()
        plt.cla()

parser = argparse.ArgumentParser(description="Get info on PrimerID results")
parser.add_argument('command', type=str, choices=["tag_dist", "likelihoods"], default="tag_dist")
parser.add_argument('input', type=argparse.FileType('r'), help="the location of primer id results file to visualise")
args = parser.parse_args()

if args.command == "tag_dist":
    tag_dist(args.input)
elif args.command == "likelihoods":
    likelihood_cutoffs(args.input)
