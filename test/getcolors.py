import os.path
import requests
from PIL import Image
from collections import defaultdict

filename = 'msquares_color.png'

if not os.path.exists(filename):
    baseurl = 'https://prideout.net/assets/'
    url = baseurl + filename
    r = requests.get(url, stream=True)
    with open(filename, 'wb') as fd:
        for chunk in r.iter_content():
            fd.write(chunk)

im = Image.open(filename)
im.split()[3].save('alpha.png')
Image.merge('RGB', im.split()[0:3]).save('rgb.png')
cols = defaultdict(set)
for pixel in im.getdata():
    r, g, b, a = pixel
    argb = '%0.2x%0.2x%0.2x%0.2x' % (a, r, g, b)
    cols[a].add(argb)
alphas = cols.keys()
alphas.sort()
final = []
for alpha in alphas:
    for col in cols[alpha]:
        final.append(col)
x = 0
for col in final:
    print '0x' + col + ',',
    x = x + 1
    if x == 5:
        print
        x = 0
