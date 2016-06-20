from collections import namedtuple, defaultdict

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

tag_dist("public_data_tags.out")
