from PIL import Image
from PIL import ImageDraw

im = Image.new('RGB', (256,256))
draw = ImageDraw.Draw(im)
draw.ellipse([128,128,512-128,512-128], (255,0,0))
draw.ellipse([150,150,512-150,512-150], (0,255,0))
draw.ellipse([128-64,128-64,128+64,128+64], (255,0,0))
draw.ellipse([128-32,128-32,128+32,128+32], (0,0,0))
del draw
im.save('rgb.png', 'PNG')

im = Image.new('RGB', (256,256))
draw = ImageDraw.Draw(im)
draw.rectangle([0, 0, 256, 256], (0,0,255,20))
for x in xrange(32, 256-32, 32):
    draw.rectangle([x, 0, x + 16, 50], (0,255,255,20))
    draw.rectangle([x, 256-50, x + 16, 256], (0,255,255,20))
    if x > 32:
        draw.rectangle([0, x, 50, x + 16], (0,255,255,20))
        draw.rectangle([256-50, x, 256, x + 16], (0,255,255,20))
draw.rectangle([50, 50, 256-50, 256-50], (255,0,255,20))
del draw
im = im.resize((128,128))
im.save('tverts.png', 'PNG')

im = Image.new('RGBA', (256,256))
draw = ImageDraw.Draw(im)
draw.rectangle([0, 0, 256, 256], (0,0,255,20))
draw.pieslice([40,40,256-40,256-40], 45, -45, (255,0,0,100))
draw.ellipse([80,80,256-80,256-80], (0,255,0,200))
draw.ellipse([35+80,80,35+256-80,256-80], (0,200,200,200))
del draw
im.save('rgba.png', 'PNG')
