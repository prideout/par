from PIL import Image
from PIL import ImageDraw
im = Image.new('RGB', (256,256))
draw = ImageDraw.Draw(im)
draw.ellipse([128,128,512-128,512-128], (255,0,0))
draw.ellipse([150,150,512-150,512-150], (0,255,0))
del draw
im.save('rgb.png', 'PNG')
