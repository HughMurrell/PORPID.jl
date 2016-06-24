from collections import namedtuple, defaultdict
from matplotlib import pyplot as plt

def tag_dist(file_name):
    tags = defaultdict(lambda: 0)
    file = open(file_name)
    for line in file:
        l = line.strip()
        if len(l) == 0:
            continue
        tags[l] = tags[l] + 1
    max_count = 0
    count_dist = defaultdict(lambda: 0)
    for count in tags.values():
        max_count = max(max_count, count)
        count_dist[count] = count_dist[count] + 1
    for i in range(max_count + 1):
        print(i, count_dist[i])
    #for tag, count in sorted(tags.items(), key=lambda x:x[1]):
    #    print(tag, count)

def likelihood_cutoffs(file_name):
    num_bins = 20

    templates = defaultdict(list)
    input_file = open(file_name)
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

#tag_dist("public_data_tags.out")
likelihood_cutoffs("likelihood_data.out")
