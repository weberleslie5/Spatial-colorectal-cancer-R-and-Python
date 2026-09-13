import os
from PIL import Image, ImageChops
home = os.path.expanduser("~")
fig = f"{home}/projects/crc_visiumhd/figures"

def trim(im):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    return im.crop(ImageChops.difference(im, bg).getbbox())

imgs = [trim(Image.open(f"{fig}/03bC_{s}.png").convert("RGB")) for s in ["P1CRC", "P2CRC", "P5CRC"]]
imgs.append(trim(Image.open(f"{fig}/03bC_legend.png").convert("RGB")))
h = max(im.height for im in imgs[:3])
imgs = [im.resize((round(im.width * h / im.height), h)) for im in imgs]
pad, margin = 40, 60    # pad between panels, margin all around
row = Image.new("RGB", (sum(im.width for im in imgs) + pad * 3 + margin * 2, h + margin * 2), (255, 255, 255))
x = margin
for im in imgs:
    row.paste(im, (x, margin + (h - im.height) // 2)); x += im.width + pad
row.save(f"{fig}/fig3bC_row.png")
print("saved fig3bC_row.png")
