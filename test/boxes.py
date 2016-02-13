import random

r = random.Random(1)

def generate_box():
    a = r.random()
    b = r.random()
    c = r.random()
    d = r.random()
    minx = min(a, b)
    maxx = max(a, b)
    miny = min(c, d)
    maxy = max(c, d)
    cx = 0.5 * (minx + maxx)
    cy = 0.5 * (miny + maxy)
    w = maxx - minx
    h = maxy - miny
    sz = r.random() * r.random()
    w = w * sz
    h = h * sz
    minx = cx - w * 0.5
    maxx = cx + w * 0.5
    miny = cy - h * 0.5
    maxy = cy + h * 0.5
    return minx, miny, maxx, maxy

for i in xrange(20):
    box = generate_box()
    print '{:.2}, {:.2}, {:.2}, {:.2},'.format(*box)
print
for i in xrange(10):
    box = generate_box()
    print '{:.2}, {:.2}, {:.2}, {:.2},'.format(*box)
