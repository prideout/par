import json

flare = json.load(open('flare.json'))
print flare

things = []

def traverse(node, parent):

    me = len(things)
    print '{:3} {}'.format(me, node['name'])
    things.append(parent)

    children = node.get('children', [])
    for child in children:
        traverse(child, me)

traverse(flare, 0)

for i in xrange(len(things)):
    print '{:3},'.format(things[i]),
    if (i + 1) % 12 == 0: print;

print '---\n', len(things)